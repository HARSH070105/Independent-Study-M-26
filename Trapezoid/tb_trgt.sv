// Unit test: trgt_sub (TrGT Temporal Gustavson subrow) against a reference matmul.
`timescale 1ns/1ps
module tb_trgt;
    localparam MAXM=32, MAXK=32, MAXN=32, DW=16, ACC=40, MS=4;
    reg clk = 0, rst = 1, start = 0;
    reg [15:0] K = 32, N = 32, m = 0;
    reg [MAXM*MAXK*DW-1:0] A_flat;
    reg [MAXK*MAXN*DW-1:0] B_flat;
    wire busy, done;
    wire [MAXN*ACC-1:0] c_out;
    wire [31:0] cyc, macs, brows, bnnz;
    always #5 clk = ~clk;

    trgt_sub #(MAXM,MAXK,MAXN,DW,ACC,MS) dut (
        .clk(clk), .rst(rst), .start(start),
        .K(K), .N(N), .m(m),
        .A_flat(A_flat), .B_flat(B_flat),
        .busy(busy), .done(done), .c_out(c_out),
        .cnt_cycles(cyc), .cnt_macs(macs), .cnt_brows(brows), .cnt_bnnz(bnnz)
    );

    integer Am [0:MAXM-1][0:MAXK-1], Bm [0:MAXK-1][0:MAXN-1];
    integer ref_c, got, errs, tot_errs = 0;
    integer i, j, kk;

    task gen(input real dA, input real dB);
        begin
            for (i=0;i<MAXM;i=i+1) for (j=0;j<MAXK;j=j+1) begin
                Am[i][j] = ($urandom % 1000) < dA*1000 ? ($urandom % 15) - 7 : 0;
                if (dA >= 1.0 && Am[i][j] == 0) Am[i][j] = 1;
                A_flat[(i*MAXK+j)*DW +: DW] = Am[i][j];
            end
            for (i=0;i<MAXK;i=i+1) for (j=0;j<MAXN;j=j+1) begin
                Bm[i][j] = ($urandom % 1000) < dB*1000 ? ($urandom % 15) - 7 : 0;
                if (dB >= 1.0 && Bm[i][j] == 0) Bm[i][j] = 1;
                B_flat[(i*MAXN+j)*DW +: DW] = Bm[i][j];
            end
        end
    endtask

    task run_row(input integer row_m);
        begin
            m = row_m;
            rst = 1; @(negedge clk); rst = 0;
            @(negedge clk); start = 1; @(negedge clk); start = 0;
            wait (done); @(negedge clk);
            errs = 0;
            for (j=0;j<N;j=j+1) begin
                ref_c = 0;
                for (kk=0;kk<K;kk=kk+1) ref_c = ref_c + Am[row_m][kk]*Bm[kk][j];
                got = $signed(c_out[j*ACC +: ACC]);
                if (got !== ref_c) begin
                    errs = errs + 1;
                    if (errs < 4) $display("  MISMATCH row %0d col %0d got %0d exp %0d", row_m, j, got, ref_c);
                end
            end
            tot_errs = tot_errs + errs;
            $display("  row=%0d K=%0d N=%0d cycles=%0d macs=%0d brows=%0d bnnz=%0d errs=%0d",
                     row_m, K, N, cyc, macs, brows, bnnz, errs);
        end
    endtask

    initial begin
        repeat (3) @(negedge clk); rst = 0;
        $display("TRGT");
        // HS x HS (single-pass & few-pass rows)
        K=32; N=32; gen(0.05, 0.05); run_row(0); run_row(1);
        gen(0.1, 0.1); run_row(0); run_row(2);
        // MS x MS (multi-pass rows: > MS=4 nonzeros in A row)
        gen(0.3, 0.3); run_row(0);
        // Dense x Dense (maximum 32/4 = 8 passes)
        gen(1.0, 1.0); run_row(0);
        // Empty A row
        gen(0.0, 0.5); run_row(0);
        // Ragged K, N
        K=21; N=19; gen(0.25, 0.25); run_row(3);
        if (tot_errs == 0) $display("TB_TRGT PASS"); else $display("TB_TRGT FAIL errs=%0d", tot_errs);
        $finish;
    end
    initial begin #50000000; $display("TIMEOUT"); $finish; end
endmodule

