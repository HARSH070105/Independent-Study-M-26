// TrGT subrow (HS x HS) -- Trapezoid Sec. III-C, Fig. 14/15.
// Gustavson dataflow: one A row -> one C row, built as a linear combination of the B rows
// selected by A's nonzeros. Each selected B row is a coordinate-ordered fiber (the FIFO fed by the
// local buffer). The merge tree compares the heads of all fibers, emits the smallest n coordinate and
// accumulates every head that matches it. A subrow owns only MS multipliers (= NMUL/4 in the paper:
// K registers, K multipliers, radix-K merge tree), so an A row with more than MS nonzeros is processed
// in passes of MS nonzeros; each pass's merged output is added into the C-row buffer.
module trgt_sub #(
    parameter MAXM = 32, MAXK = 32, MAXN = 32,
    parameter DW = 16, ACC = 40,
    parameter MS = 4                            // multipliers (= fibers merged at once) per subrow
)(
    input  wire clk, rst, start,
    input  wire [15:0] K, N,
    input  wire [15:0] m,                       // A row handled by this subrow
    input  wire [MAXM*MAXK*DW-1:0] A_flat,
    input  wire [MAXK*MAXN*DW-1:0] B_flat,
    output reg  busy, done,
    output reg  [MAXN*ACC-1:0] c_out,
    output reg  [31:0] cnt_cycles,              // pass-setup + merge cycles
    output reg  [31:0] cnt_macs,                // multiplications issued
    output reg  [31:0] cnt_brows,               // B-row fiber fetches (one per A nonzero)
    output reg  [31:0] cnt_bnnz                 // B nonzeros streamed from the memory hierarchy
);
    `define GETA(mm,kk) $signed(A_flat[((mm)*MAXK + (kk))*DW +: DW])
    `define GETB(kk,nn) $signed(B_flat[((kk)*MAXN + (nn))*DW +: DW])

    localparam S_IDLE = 0, S_PASS = 1, S_RUN = 2, S_DONE = 3;
    reg [1:0] state;
    reg [MAXK-1:0]       arem;                  // A nonzeros not yet assigned to a pass
    reg [MAXK-1:0]       amask;                 // fibers merged in the current pass (<= MS)
    reg [MAXK*MAXN-1:0]  rem;                   // remaining (unconsumed) B nonzeros per fiber
    reg signed [ACC-1:0] cbuf [0:MAXN-1];

    // ---- next pass: first MS nonzeros of the remaining A row ----
    integer q, pc;
    reg [MAXK-1:0] pmask;
    always @* begin
        pmask = '0; pc = 0;
        for (q = 0; q < MAXK; q = q + 1)
            if (arem[q] && pc < MS) begin pmask[q] = 1'b1; pc = pc + 1; end
    end

    // ---- merge step: find min head coordinate, accumulate all matching heads ----
    integer k, h, nmin, nmatch;
    reg any;
    reg signed [ACC-1:0] acc;
    reg [MAXK*MAXN-1:0] clr;
    always @* begin
        any = 0; nmin = MAXN; acc = 0; nmatch = 0; clr = '0;
        for (k = 0; k < MAXK; k = k + 1)
            if (amask[k] && rem[k*MAXN +: MAXN] != 0) begin
                h = 0;
                while (!rem[k*MAXN + h]) h = h + 1;         // head = lowest set bit
                if (h < nmin) nmin = h;
                any = 1;
            end
        if (any)
            for (k = 0; k < MAXK; k = k + 1)
                if (amask[k] && rem[k*MAXN +: MAXN] != 0 && rem[k*MAXN + nmin]) begin
                    acc = acc + `GETA(m, k) * `GETB(k, nmin);
                    clr[k*MAXN + nmin] = 1'b1;
                    nmatch = nmatch + 1;
                end
    end

    integer i, kk, nn;
    always @* for (i = 0; i < MAXN; i = i + 1) c_out[i*ACC +: ACC] = cbuf[i];

    always @(posedge clk) begin
        if (rst) begin
            state <= S_IDLE; busy <= 0; done <= 0;
            cnt_cycles <= 0; cnt_macs <= 0; cnt_brows <= 0; cnt_bnnz <= 0;
        end else begin
            done <= 0;
            case (state)
            S_IDLE: if (start) begin
                for (nn = 0; nn < MAXN; nn = nn + 1) cbuf[nn] <= 0;
                for (kk = 0; kk < MAXK; kk = kk + 1)
                    arem[kk] <= (kk < K) && (`GETA(m, kk) != 0);
                busy <= 1; state <= S_PASS;
            end
            S_PASS: begin                                   // fetch the heads of the next MS fibers
                if (pmask == 0) state <= S_DONE;
                else begin
                    cnt_cycles <= cnt_cycles + 1;
                    cnt_brows  <= cnt_brows + pc;
                    amask <= pmask; arem <= arem & ~pmask;
                    for (kk = 0; kk < MAXK; kk = kk + 1)
                        for (nn = 0; nn < MAXN; nn = nn + 1)
                            rem[kk*MAXN + nn] <= pmask[kk] && (nn < N) && (`GETB(kk, nn) != 0);
                    state <= S_RUN;
                end
            end
            S_RUN: begin
                if (!any) state <= S_PASS;
                else begin
                    cnt_cycles <= cnt_cycles + 1;
                    cnt_macs   <= cnt_macs + nmatch;
                    cnt_bnnz   <= cnt_bnnz + nmatch;
                    cbuf[nmin] <= cbuf[nmin] + acc;         // adds across passes
                    rem <= rem & ~clr;
                end
            end
            S_DONE: begin done <= 1; busy <= 0; state <= S_IDLE; end
            endcase
        end
    end
endmodule
`undef GETA
`undef GETB
