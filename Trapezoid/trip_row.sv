// PE row for Dense-IP (D x D) and TrIP (MS x D, MS x MS) -- Trapezoid Sec. III-A/B.
// One job = a group of up to RA consecutive A rows, processed against all N columns of B
// and all K tiles (KT wide). Per cycle the row:
//   1. builds A/B bitmasks, 2. picks #B columns (dynamic packing, <= CB, effectual MACs <= NMUL),
//   3. MFIU intersects + routes, 4. NMUL multipliers, 5. segmented reduction per (A row, B col),
//   6. accumulates into a local banked output buffer (cbuf).
// A/B "distribution networks" are modeled as index-based gathers (functionally a Benes network).
// Dense mode: ra=1, cb=1, bitmasks forced to all-ones (sparsity is not exploited).
module trip_row #(
    parameter MAXM = 32, MAXK = 32, MAXN = 32,
    parameter DW = 16, ACC = 40,
    parameter KT = 16, RA = 4, CB = 4, NMUL = 16
)(
    input  wire clk, rst, start,
    input  wire dense_mode,
    input  wire [15:0] M, K, N,
    input  wire [15:0] m_base,
    input  wire [3:0]  m_rows,                       // valid rows in group (<= RA)
    input  wire [MAXM*MAXK*DW-1:0] A_flat,
    input  wire [MAXK*MAXN*DW-1:0] B_flat,
    output reg  busy, done,
    output reg  [RA*MAXN*ACC-1:0] c_out,
    output reg  [31:0] cnt_stream_cycles,            // cycles spent streaming B
    output reg  [31:0] cnt_eff_macs                  // effectual multiplications issued
);
    localparam RW = $clog2(RA), CW = $clog2(CB), KW = $clog2(KT);
    localparam S_IDLE = 0, S_PICK = 1, S_STREAM = 2, S_DONE = 3;

    reg [1:0]  state;
    reg [15:0] kt, nkt, n;
    reg [RA-1:0] rows_rem, sub;
    reg signed [ACC-1:0] cbuf [0:RA-1][0:MAXN-1];

    `define GETA(m,k) $signed(A_flat[((m)*MAXK + (k))*DW +: DW])
    `define GETB(k,nn) $signed(B_flat[((k)*MAXN + (nn))*DW +: DW])
    function automatic integer pop(input [KT-1:0] v);
        integer i; begin pop = 0; for (i = 0; i < KT; i = i + 1) pop = pop + v[i]; end
    endfunction

    // ---- bitmasks of the current K tile ----
    reg [RA*KT-1:0] a_mask_all;                      // all valid rows (before subset selection)
    reg [RA*KT-1:0] a_mask_sub;
    reg [CB*KT-1:0] b_mask;
    reg [$clog2(CB+1)-1:0] cb;
    integer r, c, k, cum, cnt, tmp;
    reg stop;
    always @* begin
        a_mask_all = '0; a_mask_sub = '0; b_mask = '0;
        for (r = 0; r < RA; r = r + 1)
            for (k = 0; k < KT; k = k + 1)
                if (r < m_rows && (kt*KT + k) < K)
                    a_mask_all[r*KT + k] = dense_mode ? 1'b1 : (`GETA(m_base + r, kt*KT + k) != 0);
        for (r = 0; r < RA; r = r + 1)
            if (sub[r]) a_mask_sub[r*KT +: KT] = a_mask_all[r*KT +: KT];
        for (c = 0; c < CB; c = c + 1)
            for (k = 0; k < KT; k = k + 1)
                if ((n + c) < N && (kt*KT + k) < K)
                    b_mask[c*KT + k] = dense_mode ? 1'b1 : (`GETB(kt*KT + k, n + c) != 0);
        // dynamic packing of B columns: largest cb with effectual MACs <= NMUL
        cb = 1; cum = 0; stop = 0;
        for (c = 0; c < CB; c = c + 1) begin
            tmp = 0;
            for (r = 0; r < RA; r = r + 1)
                tmp = tmp + pop(a_mask_sub[r*KT +: KT] & b_mask[c*KT +: KT]);
            if (!stop && !dense_mode && (n + c) < N && (cum + tmp) <= NMUL) begin
                cum = cum + tmp; cb = c + 1;
            end else stop = 1;   // columns must be contiguous: stop at the first that does not fit
        end
        if (dense_mode) cb = 1;
    end

    // ---- MFIU ----
    wire [$clog2(RA*CB*KT+1)-1:0] n_eff;
    wire [NMUL-1:0] slot_v;
    wire [NMUL*RW-1:0] slot_r;
    wire [NMUL*CW-1:0] slot_c;
    wire [NMUL*KW-1:0] slot_k;
    mfiu #(.K(KT), .RA(RA), .CB(CB), .NMUL(NMUL)) u_mfiu (
        .a_mask(a_mask_sub), .b_mask(b_mask), .b_active(cb),
        .n_eff(n_eff), .slot_v(slot_v), .slot_r(slot_r), .slot_c(slot_c), .slot_k(slot_k));

    // ---- multipliers + segmented reduction ----
    reg signed [ACC-1:0] pair_sum [0:RA*CB-1];
    integer s, p;
    reg [RW-1:0] sr; reg [CW-1:0] sc; reg [KW-1:0] sk;
    always @* begin
        for (p = 0; p < RA*CB; p = p + 1) pair_sum[p] = 0;
        for (s = 0; s < NMUL; s = s + 1)
            if (slot_v[s]) begin
                sr = slot_r[s*RW +: RW]; sc = slot_c[s*CW +: CW]; sk = slot_k[s*KW +: KW];
                pair_sum[sr*CB + sc] = pair_sum[sr*CB + sc] +
                    `GETA(m_base + sr, kt*KT + sk) * `GETB(kt*KT + sk, n + sc);
            end
    end

    // ---- picking the next row subset that fits in the multipliers ----
    reg [RA-1:0] pick; integer pcum, pcnt, rp;
    always @* begin
        pick = '0; pcum = 0; pcnt = 0;
        for (rp = 0; rp < RA; rp = rp + 1)
            if (rows_rem[rp] && !(dense_mode && pcnt > 0) &&
                (pcum + pop(a_mask_all[rp*KT +: KT])) <= NMUL) begin
                pick[rp] = 1'b1; pcum = pcum + pop(a_mask_all[rp*KT +: KT]); pcnt = pcnt + 1;
            end
    end

    integer i, j;
    always @* begin
        for (i = 0; i < RA; i = i + 1)
            for (j = 0; j < MAXN; j = j + 1)
                c_out[(i*MAXN + j)*ACC +: ACC] = cbuf[i][j];
    end

    integer ii, jj, rr, cc;
    wire last_n = (n + cb) >= N;
    always @(posedge clk) begin
        if (rst) begin
            state <= S_IDLE; busy <= 0; done <= 0; kt <= 0; n <= 0; sub <= 0; rows_rem <= 0;
            cnt_stream_cycles <= 0; cnt_eff_macs <= 0;
        end else begin
            done <= 0;
            case (state)
            S_IDLE: if (start) begin
                for (ii = 0; ii < RA; ii = ii + 1) for (jj = 0; jj < MAXN; jj = jj + 1) cbuf[ii][jj] <= 0;
                nkt <= (K + KT - 1) / KT; kt <= 0; n <= 0;
                rows_rem <= (1 << m_rows) - 1;
                busy <= 1; state <= S_PICK;
            end
            S_PICK: begin
                sub <= pick; n <= 0;
                // skip subsets with no effectual A elements (all-zero tile rows)
                if (pcum == 0 && !dense_mode) begin
                    rows_rem <= rows_rem & ~pick;
                    if ((rows_rem & ~pick) == 0) begin
                        if (kt + 1 >= nkt) state <= S_DONE;
                        else begin kt <= kt + 1; rows_rem <= (1 << m_rows) - 1; end
                    end
                end else state <= S_STREAM;
            end
            S_STREAM: begin
                cnt_stream_cycles <= cnt_stream_cycles + 1;
                cnt_eff_macs <= cnt_eff_macs + n_eff;
                for (rr = 0; rr < RA; rr = rr + 1)
                    for (cc = 0; cc < CB; cc = cc + 1)
                        if (sub[rr] && cc < cb && (n + cc) < N)
                            cbuf[rr][n + cc] <= cbuf[rr][n + cc] + pair_sum[rr*CB + cc];
                if (last_n) begin
                    n <= 0; rows_rem <= rows_rem & ~sub;
                    if ((rows_rem & ~sub) == 0) begin
                        if (kt + 1 >= nkt) state <= S_DONE;
                        else begin kt <= kt + 1; rows_rem <= (1 << m_rows) - 1; state <= S_PICK; end
                    end else state <= S_PICK;
                end else n <= n + cb;
            end
            S_DONE: begin done <= 1; busy <= 0; state <= S_IDLE; end
            endcase
        end
    end

endmodule

`undef GETA
`undef GETB
