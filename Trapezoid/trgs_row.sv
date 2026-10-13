// TrGS PE row (HS x MS, HS x D) -- Trapezoid Sec. III-D, Fig. 16/17.
// Gustavson dataflow with B rows processed *spatially*: one PE row computes one C row.
// For each nonzero A[m,k] (broadcast through the A distribution network) the row streams up to
// W nonzeros of B row k per cycle (W = global-cache word width), the B distribution network routes
// each to the multiplier owning its n coordinate, and the products accumulate into per-n registers.
module trgs_row #(
    parameter MAXM = 32, MAXK = 32, MAXN = 32,
    parameter DW = 16, ACC = 40,
    parameter W = 16,                          // B nonzeros delivered per cycle (cache word)
    parameter NT = 16                          // N tile = multipliers in the PE row (coordinate tiling on N)
)(
    input  wire clk, rst, start,
    input  wire [15:0] K, N,
    input  wire [15:0] m,
    input  wire [MAXM*MAXK*DW-1:0] A_flat,
    input  wire [MAXK*MAXN*DW-1:0] B_flat,
    output reg  busy, done,
    output reg  [MAXN*ACC-1:0] c_out,
    output reg  [31:0] cnt_cycles,
    output reg  [31:0] cnt_macs,
    output reg  [31:0] cnt_brows,
    output reg  [31:0] cnt_bnnz
);
    `define GETA(mm,kk) $signed(A_flat[((mm)*MAXK + (kk))*DW +: DW])
    `define GETB(kk,nn) $signed(B_flat[((kk)*MAXN + (nn))*DW +: DW])

    localparam S_IDLE = 0, S_LOAD = 1, S_RUN = 2, S_DONE = 3;
    reg [1:0] state;
    reg [MAXK-1:0] arem;                        // A nonzeros not yet scheduled
    reg [15:0]     cur_k, nt;
    reg [MAXK-1:0] anz;                         // A nonzero mask (reloaded per N tile)
    reg [MAXN-1:0] brem;                        // remaining nonzeros of the current B row
    reg signed [ACC-1:0] cbuf [0:MAXN-1];

    // ---- next A nonzero and its B-row mask ----
    integer q, nn_, kq;
    reg nk_valid; reg [15:0] nk; reg [MAXN-1:0] nmask;
    always @* begin
        nk_valid = 0; nk = 0; nmask = '0;
        for (q = MAXK-1; q >= 0; q = q - 1) if (arem[q]) begin nk_valid = 1; nk = q; end
        if (nk_valid)
            for (nn_ = 0; nn_ < MAXN; nn_ = nn_ + 1)
                nmask[nn_] = (nn_ < N) && (nn_ >= nt*NT) && (nn_ < (nt+1)*NT) && (`GETB(nk, nn_) != 0);
    end

    // ---- this cycle's chunk: lowest W set bits of brem ----
    reg [MAXN-1:0] chunk; integer cc, cntw;
    reg signed [DW-1:0] a_cur;
    always @* begin
        chunk = '0; cntw = 0;
        for (cc = 0; cc < MAXN; cc = cc + 1)
            if (brem[cc] && cntw < W) begin chunk[cc] = 1'b1; cntw = cntw + 1; end
        a_cur = `GETA(m, cur_k);
    end

    integer i, jj, kk;
    always @* for (i = 0; i < MAXN; i = i + 1) c_out[i*ACC +: ACC] = cbuf[i];

    wire [MAXN-1:0] brem_after = brem & ~chunk;
    always @(posedge clk) begin
        if (rst) begin
            state <= S_IDLE; busy <= 0; done <= 0;
            cnt_cycles <= 0; cnt_macs <= 0; cnt_brows <= 0; cnt_bnnz <= 0;
        end else begin
            done <= 0;
            case (state)
            S_IDLE: if (start) begin
                for (jj = 0; jj < MAXN; jj = jj + 1) cbuf[jj] <= 0;
                for (kk = 0; kk < MAXK; kk = kk + 1) begin
                    arem[kk] <= (kk < K) && (`GETA(m, kk) != 0);
                    anz[kk]  <= (kk < K) && (`GETA(m, kk) != 0);
                end
                nt <= 0; busy <= 1; state <= S_LOAD;
            end
            S_LOAD: begin
                if (!nk_valid) state <= S_DONE;   // (A row empty)
                else begin
                    cur_k <= nk; brem <= nmask; arem[nk] <= 1'b0;
                    cnt_brows <= cnt_brows + 1; state <= S_RUN;
                end
            end
            S_RUN: begin
                cnt_cycles <= cnt_cycles + 1;
                cnt_macs   <= cnt_macs + cntw;
                cnt_bnnz   <= cnt_bnnz + cntw;
                for (cc = 0; cc < MAXN; cc = cc + 1)
                    if (chunk[cc]) cbuf[cc] <= cbuf[cc] + a_cur * `GETB(cur_k, cc);
                if (brem_after == 0) begin
                    if (nk_valid) begin     // overlap next B-row fetch with the last chunk
                        cur_k <= nk; brem <= nmask; arem[nk] <= 1'b0;
                        cnt_brows <= cnt_brows + 1;
                    end else if ((nt + 1) * NT < N) begin   // next N tile: re-stream A's nonzeros
                        nt <= nt + 1; arem <= anz; state <= S_LOAD;
                    end else state <= S_DONE;
                end else brem <= brem_after;
            end
            S_DONE: begin done <= 1; busy <= 0; state <= S_IDLE; end
            endcase
        end
    end
endmodule
`undef GETA
`undef GETB
