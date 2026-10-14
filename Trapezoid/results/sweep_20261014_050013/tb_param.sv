// AUTO-GENERATED TESTBENCH FROM parameters.json
`timescale 1ns/1ps
module tb_param;
    localparam MAXM=64, MAXK=64, MAXN=64, DW=16, ACC=40;
    localparam NROWS=4, KT=16, RA=4, CB=4, NMUL=16, W=16;
    localparam ROWB = MAXN*ACC;

    reg clk = 0, rst = 1, start = 0; reg [1:0] mode = 0;
    reg [15:0] M = 16, K = 32, N = 32;
    reg [MAXM*MAXK*DW-1:0] A_flat; reg [MAXK*MAXN*DW-1:0] B_flat;
    wire busy, done; wire [MAXM*MAXN*ACC-1:0] C_flat;
    wire [31:0] cyc_total, cnt_macs, cnt_engine_cyc, cnt_brows, cnt_bnnz;
    always #5 clk = ~clk;

    trapezoid_top #(MAXM,MAXK,MAXN,DW,ACC,NROWS,KT,RA,CB,NMUL,W) dut (
        .clk(clk), .rst(rst), .start(start), .mode(mode), .M(M), .K(K), .N(N),
        .A_flat(A_flat), .B_flat(B_flat), .busy(busy), .done(done), .C_flat(C_flat),
        .cyc_total(cyc_total), .cnt_macs(cnt_macs), .cnt_engine_cyc(cnt_engine_cyc),
        .cnt_brows(cnt_brows), .cnt_bnnz(cnt_bnnz));

    integer Am [0:MAXM-1][0:MAXK-1], Bm [0:MAXK-1][0:MAXN-1];
    integer i, j, kk, errs, tot_errs = 0, nruns = 0, eff, ref_c, got, fd;
    integer a_i, b_i, m_i, run_idx;
    integer nzA, nzB, nzC;
    real densA [0:2];
    real densB [0:2];
    string mname [0:3];

    task gen(input real dA, input real dB);
        begin
            nzA = 0; nzB = 0;
            for (i=0;i<MAXM;i=i+1) for (j=0;j<MAXK;j=j+1) begin
                Am[i][j] = ($urandom % 10000) < dA*10000 ? (($urandom % 15) - 7) : 0;
                if (dA >= 1.0 && Am[i][j] == 0) Am[i][j] = 1;
                A_flat[(i*MAXK+j)*DW +: DW] = Am[i][j];
            end
            for (i=0;i<MAXK;i=i+1) for (j=0;j<MAXN;j=j+1) begin
                Bm[i][j] = ($urandom % 10000) < dB*10000 ? (($urandom % 15) - 7) : 0;
                if (dB >= 1.0 && Bm[i][j] == 0) Bm[i][j] = 1;
                B_flat[(i*MAXN+j)*DW +: DW] = Bm[i][j];
            end
        end
    endtask

    task run(input [1:0] md, input real dA, input real dB);
        begin
            mode = md;
            rst = 1; repeat (2) @(negedge clk); rst = 0; @(negedge clk);
            start = 1; @(negedge clk); start = 0;
            wait (done); @(negedge clk);
            // Verify against reference matmul
            errs = 0; eff = 0; nzA = 0; nzB = 0; nzC = 0;
            for (i=0;i<M;i=i+1) for (kk=0;kk<K;kk=kk+1) if (Am[i][kk] != 0) nzA = nzA + 1;
            for (kk=0;kk<K;kk=kk+1) for (j=0;j<N;j=j+1) if (Bm[kk][j] != 0) nzB = nzB + 1;
            for (i=0;i<M;i=i+1) for (j=0;j<N;j=j+1) begin
                ref_c = 0;
                for (kk=0;kk<K;kk=kk+1) begin
                    ref_c = ref_c + Am[i][kk]*Bm[kk][j];
                    if (Am[i][kk] != 0 && Bm[kk][j] != 0) eff = eff + 1;
                end
                if (ref_c != 0) nzC = nzC + 1;
                got = $signed(C_flat[(i*MAXN+j)*ACC +: ACC]);
                if (got !== ref_c) begin
                    errs = errs + 1;
                    if (errs < 4) $display("  MISMATCH mode %0d C[%0d][%0d] got %0d exp %0d", md, i, j, got, ref_c);
                end
            end
            tot_errs = tot_errs + errs; nruns = nruns + 1;
            $display("%-6s dA=%4.2f dB=%4.2f M=%0d K=%0d N=%0d cyc=%6d issMAC=%6d effMAC=%6d %s",
                     mname[md], dA, dB, M, K, N, cyc_total, cnt_macs, eff, errs ? "FAIL" : "ok");
            $fdisplay(fd, "%0s,%f,%f,%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d", mname[md], dA, dB, M, K, N,
                      cyc_total, cnt_macs, eff, cnt_bnnz, nzA, nzB, nzC);
        end
    endtask

    initial begin
        mname[0]="DENSE"; mname[1]="TRIP"; mname[2]="TRGT"; mname[3]="TRGS";
        densA[0] = 1.0;
        densA[1] = 0.5;
        densA[2] = 0.1;
        densB[0] = 1.0;
        densB[1] = 0.5;
        densB[2] = 0.1;

        fd = $fopen("sweep.csv", "w");
        $fdisplay(fd, "mode,dA,dB,M,K,N,cycles,issued_macs,eff_macs,b_nnz_streamed,nnzA,nnzB,nnzC");
        repeat (3) @(negedge clk);

        for (a_i = 0; a_i < 3; a_i = a_i + 1) begin
            for (b_i = 0; b_i < 3; b_i = b_i + 1) begin
                for (run_idx = 0; run_idx < 1; run_idx = run_idx + 1) begin
                    gen(densA[a_i], densB[b_i]);
                    run(2'd0, densA[a_i], densB[b_i]);
                    run(2'd1, densA[a_i], densB[b_i]);
                    run(2'd2, densA[a_i], densB[b_i]);
                    run(2'd3, densA[a_i], densB[b_i]);
                end
            end
        end

        // Ragged test case: ragged_prime
        M = 13; K = 21; N = 19;
        gen(0.3, 0.3);
        run(2'd0, 0.3, 0.3);
        run(2'd1, 0.3, 0.3);
        run(2'd2, 0.3, 0.3);
        run(2'd3, 0.3, 0.3);

        $fclose(fd);
        if (tot_errs == 0) $display("TB_PARAM PASS (%0d runs)", nruns);
        else $display("TB_PARAM FAIL (%0d errors in %0d runs)", tot_errs, nruns);
        $finish;
    end
    initial begin #5000000000; $display("TIMEOUT"); $finish; end
endmodule
