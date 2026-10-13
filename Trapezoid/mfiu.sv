// Multi-fiber intersection unit (MFIU) -- Trapezoid Sec. III-B, Fig. 11.
// Combinational. Intersects RA A-row bitmasks with CB B-column bitmasks
// (RA*CB pairwise ANDs), computes a cumulative (prefix) sum of matching bits
// across all pairs, and emits per-slot routing metadata:
//   slot s (= multiplier index) -> {pair r,c ; coordinate k}
// Slots are assigned in pair order (r-major, then c), increasing k, so all
// partial products of one output element C[r][c] are contiguous (needed by
// the segmented merge-reduction tree).
module mfiu #(
    parameter K    = 8,    // coiteration-tile length (bitmask width)
    parameter RA   = 4,    // A rows intersected at once
    parameter CB   = 4,    // B columns intersected at once
    parameter NMUL = 16    // multipliers per PE row
)(
    input  wire [RA*K-1:0]            a_mask,   // row r at [r*K +: K]
    input  wire [CB*K-1:0]            b_mask,   // col c at [c*K +: K]
    input  wire [$clog2(CB+1)-1:0]    b_active, // number of B columns used (<=CB)
    output reg  [$clog2(RA*CB*K+1)-1:0] n_eff,  // total effectual MACs (may exceed NMUL)
    output reg  [NMUL-1:0]            slot_v,
    output reg  [NMUL*$clog2(RA)-1:0] slot_r,   // (RA==1 -> width 0; keep RA>=2)
    output reg  [NMUL*$clog2(CB)-1:0] slot_c,
    output reg  [NMUL*$clog2(K)-1:0]  slot_k
);
    localparam RW = $clog2(RA), CW = $clog2(CB), KW = $clog2(K);

    integer r, c, k, s;
    reg [K-1:0] and_m;
    integer base;

    always @* begin
        slot_v = '0; slot_r = '0; slot_c = '0; slot_k = '0;
        base = 0;
        for (r = 0; r < RA; r = r + 1)
            for (c = 0; c < CB; c = c + 1) begin
                if (c < b_active) begin
                    and_m = a_mask[r*K +: K] & b_mask[c*K +: K];   // step 1: AND
                    s = base;
                    for (k = 0; k < K; k = k + 1)                  // step 2/3: prefix-sum index
                        if (and_m[k]) begin
                            if (s < NMUL) begin
                                slot_v[s]            = 1'b1;
                                slot_r[s*RW +: RW]   = r[RW-1:0];
                                slot_c[s*CW +: CW]   = c[CW-1:0];
                                slot_k[s*KW +: KW]   = k[KW-1:0];
                            end
                            s = s + 1;
                        end
                    base = s;
                end
            end
        n_eff = base;
    end
endmodule
