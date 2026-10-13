// Trapezoid top level (scaled-down, parameterized).
//   mode 0: Dense IP     (D x D)            -- trip_row, one A row per PE row
//   mode 1: TrIP         (MS x D, MS x MS)  -- trip_row, up to RA A rows x CB B cols per PE row
//   mode 2: TrGT         (HS x HS)          -- trgt_sub, 4 subrows per PE row (NSUB = 4*NROWS)
//   mode 3: TrGS         (HS x MS, HS x D)  -- trgs_row, one PE row per A row
// A and B are presented as flat dense arrays; every dataflow derives its own sparse view
// (bitmask for D/MS, fiber lists for HS) and only touches/streams nonzeros as the paper describes.
// Rows are dispatched dynamically to the first idle unit (greedy load balancing).
// Reset (rst) clears counters; pulse rst before each run.
module trapezoid_top #(
    parameter MAXM = 32, MAXK = 32, MAXN = 32,
    parameter DW = 16, ACC = 40,
    parameter NROWS = 4,                 // PE rows (paper: 128)
    parameter KT = 16, RA = 4, CB = 4,
    parameter NMUL = 16,                 // multipliers per PE row (paper: 128)
    parameter W = 16                     // cache word (nonzeros/cycle) for TrGS
)(
    input  wire clk, rst, start,
    input  wire [1:0]  mode,
    input  wire [15:0] M, K, N,
    input  wire [MAXM*MAXK*DW-1:0] A_flat,
    input  wire [MAXK*MAXN*DW-1:0] B_flat,
    output reg  busy, done,
    output reg  [MAXM*MAXN*ACC-1:0] C_flat,
    output reg  [31:0] cyc_total,        // end-to-end cycles
    output wire [31:0] cnt_macs,         // multiplications issued by all engines
    output wire [31:0] cnt_engine_cyc,   // sum over engines of active cycles
    output wire [31:0] cnt_brows,        // B-row fiber fetches (Gustavson modes)
    output wire [31:0] cnt_bnnz          // B nonzeros streamed (Gustavson modes)
);
    localparam NSUB = 4*NROWS;
    localparam NU = NSUB;                // max units over all modes
    localparam ROWB = MAXN*ACC;

    // ---------------- per-unit dispatch registers ----------------
    reg [NU-1:0] ustart, ubusy;
    reg [15:0] ujob  [0:NU-1];
    reg [3:0]  urows [0:NU-1];

    // ---------------- engines ----------------
    wire [NROWS-1:0] td_trip, td_gs;
    wire [NSUB-1:0]  td_gt;
    wire [NROWS*RA*ROWB-1:0] c_trip;
    wire [NSUB*ROWB-1:0]     c_gt;
    wire [NROWS*ROWB-1:0]    c_gs;
    wire [NROWS*32-1:0] trip_cyc, trip_mac;
    wire [NSUB*32-1:0]  gt_cyc, gt_mac, gt_br, gt_bn;
    wire [NROWS*32-1:0] gs_cyc, gs_mac, gs_br, gs_bn;

    genvar g;
    generate
        for (g = 0; g < NROWS; g = g + 1) begin : g_trip
            wire bsy;
            trip_row #(MAXM,MAXK,MAXN,DW,ACC,KT,RA,CB,NMUL) u (
                .clk(clk), .rst(rst), .start(ustart[g] && mode <= 2'd1),
                .dense_mode(mode == 2'd0), .M(M), .K(K), .N(N),
                .m_base(ujob[g]), .m_rows(urows[g]), .A_flat(A_flat), .B_flat(B_flat),
                .busy(bsy), .done(td_trip[g]), .c_out(c_trip[g*RA*ROWB +: RA*ROWB]),
                .cnt_stream_cycles(trip_cyc[g*32 +: 32]), .cnt_eff_macs(trip_mac[g*32 +: 32]));
        end
        for (g = 0; g < NSUB; g = g + 1) begin : g_gt
            wire bsy;
            trgt_sub #(MAXM,MAXK,MAXN,DW,ACC,NMUL/4) u (
                .clk(clk), .rst(rst), .start(ustart[g] && mode == 2'd2),
                .K(K), .N(N), .m(ujob[g]), .A_flat(A_flat), .B_flat(B_flat),
                .busy(bsy), .done(td_gt[g]), .c_out(c_gt[g*ROWB +: ROWB]),
                .cnt_cycles(gt_cyc[g*32 +: 32]), .cnt_macs(gt_mac[g*32 +: 32]),
                .cnt_brows(gt_br[g*32 +: 32]), .cnt_bnnz(gt_bn[g*32 +: 32]));
        end
        for (g = 0; g < NROWS; g = g + 1) begin : g_gs
            wire bsy;
            trgs_row #(MAXM,MAXK,MAXN,DW,ACC,W,NMUL) u (
                .clk(clk), .rst(rst), .start(ustart[g] && mode == 2'd3),
                .K(K), .N(N), .m(ujob[g]), .A_flat(A_flat), .B_flat(B_flat),
                .busy(bsy), .done(td_gs[g]), .c_out(c_gs[g*ROWB +: ROWB]),
                .cnt_cycles(gs_cyc[g*32 +: 32]), .cnt_macs(gs_mac[g*32 +: 32]),
                .cnt_brows(gs_br[g*32 +: 32]), .cnt_bnnz(gs_bn[g*32 +: 32]));
        end
    endgenerate

    // ---------------- aggregate performance counters ----------------
    integer a;
    reg [31:0] s_mac, s_cyc, s_br, s_bn;
    always @* begin
        s_mac = 0; s_cyc = 0; s_br = 0; s_bn = 0;
        if (mode <= 2'd1) for (a = 0; a < NROWS; a = a + 1) begin
            s_mac = s_mac + trip_mac[a*32 +: 32]; s_cyc = s_cyc + trip_cyc[a*32 +: 32]; end
        else if (mode == 2'd2) for (a = 0; a < NSUB; a = a + 1) begin
            s_mac = s_mac + gt_mac[a*32 +: 32]; s_cyc = s_cyc + gt_cyc[a*32 +: 32];
            s_br = s_br + gt_br[a*32 +: 32]; s_bn = s_bn + gt_bn[a*32 +: 32]; end
        else for (a = 0; a < NROWS; a = a + 1) begin
            s_mac = s_mac + gs_mac[a*32 +: 32]; s_cyc = s_cyc + gs_cyc[a*32 +: 32];
            s_br = s_br + gs_br[a*32 +: 32]; s_bn = s_bn + gs_bn[a*32 +: 32]; end
    end
    assign cnt_macs = s_mac; assign cnt_engine_cyc = s_cyc;
    assign cnt_brows = s_br; assign cnt_bnnz = s_bn;

    // ---------------- scheduler ----------------
    localparam T_IDLE = 0, T_RUN = 1;
    reg tstate;
    integer nj, njobs, nu, u, i, mb, mr;
    reg anyb, avail, ud;
    always @(posedge clk) begin
        if (rst) begin
            tstate <= T_IDLE; busy <= 0; done <= 0; ustart <= 0; ubusy <= 0;
            cyc_total <= 0; nj = 0; C_flat <= '0;
        end else begin
            done <= 0;
            case (tstate)
            T_IDLE: if (start) begin
                C_flat <= '0; nj = 0; cyc_total <= 0; busy <= 1; tstate <= T_RUN;
            end
            T_RUN: begin
                cyc_total <= cyc_total + 1;
                nu    = (mode == 2'd2) ? NSUB : NROWS;
                njobs = (mode == 2'd1) ? (M + RA - 1) / RA : M;
                anyb = 0;
                for (u = 0; u < NU; u = u + 1) begin
                    if (u < nu) begin
                        ustart[u] <= 1'b0;
                        ud = (mode <= 2'd1) ? td_trip[u % NROWS] : (mode == 2'd2) ? td_gt[u] : td_gs[u % NROWS];
                        if (ubusy[u] && ud) begin               // job finished: write back result rows
                            ubusy[u] <= 1'b0;
                            if (mode <= 2'd1) begin
                                for (i = 0; i < urows[u]; i = i + 1)
                                    C_flat[(ujob[u] + i)*ROWB +: ROWB] <= c_trip[(u*RA + i)*ROWB +: ROWB];
                            end else if (mode == 2'd2)
                                C_flat[ujob[u]*ROWB +: ROWB] <= c_gt[u*ROWB +: ROWB];
                            else
                                C_flat[ujob[u]*ROWB +: ROWB] <= c_gs[u*ROWB +: ROWB];
                        end
                        avail = !ubusy[u] || ud;
                        if (avail && nj < njobs) begin          // dynamic dispatch of the next job
                            mb = (mode == 2'd1) ? nj*RA : nj;
                            mr = (mode == 2'd1) ? ((M - mb) < RA ? (M - mb) : RA) : 1;
                            ujob[u] <= mb; urows[u] <= mr;
                            ubusy[u] <= 1'b1; ustart[u] <= 1'b1;
                            nj = nj + 1; anyb = 1;
                        end else if (ubusy[u] && !ud) anyb = 1;
                    end
                end
                if (nj >= njobs && !anyb) begin
                    done <= 1; busy <= 0; tstate <= T_IDLE;
                end
            end
            endcase
        end
    end
endmodule
