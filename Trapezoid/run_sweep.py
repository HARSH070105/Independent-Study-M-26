#!/usr/bin/env python3
"""Configurable sweep runner for Trapezoid accelerator.

Reads a parameter configuration file (JSON) to set matrix dimensions, densities,
methods (dataflows), and number of runs per point. Dynamically generates
the SystemVerilog testbench, compiles with iverilog, executes with vvp,
verifies outputs against reference matmul, logs results to sweep.csv,
and generates summary tables and publication figures.

Usage:
    python run_sweep.py [--config parameters.json]
    python run_sweep.py --M 16 --K 32 --N 32 --densities-A 0.5 0.2 --densities-B 0.5 0.2
"""

import argparse
import csv
import datetime
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Check for iverilog / vvp in standard PATH or miniconda env
def find_tool(name):
    # Check if tool directory in conda needs to be added to PATH
    user_home = Path.home()
    conda_bin = user_home / "miniconda3" / "envs" / "test" / "Library" / "bin"
    conda_ivl = user_home / "miniconda3" / "envs" / "test" / "Library" / "lib" / "ivl"
    if conda_bin.exists():
        path_str = os.environ.get("PATH", "")
        if str(conda_bin) not in path_str:
            os.environ["PATH"] = f"{conda_bin};{conda_ivl};{path_str}"

    p = shutil.which(name)
    if p:
        return p
    # Search common conda paths
    conda_candidates = [
        user_home / "miniconda3" / "envs" / "test" / "Library" / "bin" / f"{name}.exe",
        user_home / "miniconda3" / "Library" / "bin" / f"{name}.exe",
        user_home / "miniconda" / "Library" / "bin" / f"{name}.exe",
        user_home / "anaconda3" / "Library" / "bin" / f"{name}.exe",
    ]
    for c in conda_candidates:
        if c.exists():
            bin_dir = str(c.parent)
            ivl_dir = str(c.parent.parent / "lib" / "ivl")
            path_str = os.environ.get("PATH", "")
            if bin_dir not in path_str:
                os.environ["PATH"] = f"{bin_dir};{ivl_dir};{path_str}"
            return str(c)
    return None

def find_python_with_plotting():
    # Test current python
    try:
        import numpy, matplotlib
        return sys.executable
    except ImportError:
        pass
    user_home = Path.home()
    py_candidates = [
        user_home / "miniconda3" / "envs" / "test" / "python.exe",
        user_home / "miniconda3" / "python.exe",
        user_home / "anaconda3" / "python.exe",
    ]
    for c in py_candidates:
        if c.exists():
            return str(c)
    return sys.executable


def generate_testbench(params, tb_path):
    hw = params["hardware"]
    mat = params["matrix"]
    sweep = params["sweep"]
    ragged = params.get("ragged_cases", [])

    dens_a = sweep["densities_A"]
    dens_b = sweep["densities_B"]
    methods = sweep["methods"]
    num_runs = sweep.get("num_runs_per_point", 1)

    method_to_mode = {"DENSE": 0, "TRIP": 1, "TRGT": 2, "TRGS": 3}
    active_modes = [method_to_mode[m] for m in methods if m in method_to_mode]

    # Ensure MAXM >= M, MAXK >= K, MAXN >= N
    max_m = max(hw["MAXM"], mat["M"], *(r.get("M", 0) for r in ragged))
    max_k = max(hw["MAXK"], mat["K"], *(r.get("K", 0) for r in ragged))
    max_n = max(hw["MAXN"], mat["N"], *(r.get("N", 0) for r in ragged))

    tb_sv = f"""// AUTO-GENERATED TESTBENCH FROM parameters.json
`timescale 1ns/1ps
module tb_param;
    localparam MAXM={max_m}, MAXK={max_k}, MAXN={max_n}, DW=16, ACC=40;
    localparam NROWS={hw['NROWS']}, KT={hw['KT']}, RA={hw['RA']}, CB={hw['CB']}, NMUL={hw['NMUL']}, W={hw['W']};
    localparam ROWB = MAXN*ACC;

    reg clk = 0, rst = 1, start = 0; reg [1:0] mode = 0;
    reg [15:0] M = {mat['M']}, K = {mat['K']}, N = {mat['N']};
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
    real densA [0:{len(dens_a)-1}];
    real densB [0:{len(dens_b)-1}];
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
"""

    for idx, d in enumerate(dens_a):
        tb_sv += f"        densA[{idx}] = {d};\n"
    for idx, d in enumerate(dens_b):
        tb_sv += f"        densB[{idx}] = {d};\n"

    tb_sv += """
        fd = $fopen("sweep.csv", "w");
        $fdisplay(fd, "mode,dA,dB,M,K,N,cycles,issued_macs,eff_macs,b_nnz_streamed,nnzA,nnzB,nnzC");
        repeat (3) @(negedge clk);

        for (a_i = 0; a_i < """ + str(len(dens_a)) + """; a_i = a_i + 1) begin
            for (b_i = 0; b_i < """ + str(len(dens_b)) + """; b_i = b_i + 1) begin
                for (run_idx = 0; run_idx < """ + str(num_runs) + """; run_idx = run_idx + 1) begin
                    gen(densA[a_i], densB[b_i]);
"""
    for mode in active_modes:
        tb_sv += f"                    run(2'd{mode}, densA[a_i], densB[b_i]);\n"

    tb_sv += """                end
            end
        end
"""

    # Add ragged test cases if any
    for r in ragged:
        m_val = r.get("M", 13)
        k_val = r.get("K", 21)
        n_val = r.get("N", 19)
        da_val = r.get("dA", 0.3)
        db_val = r.get("dB", 0.3)
        tb_sv += f"""
        // Ragged test case: {r.get('name', 'ragged')}
        M = {m_val}; K = {k_val}; N = {n_val};
        gen({da_val}, {db_val});
"""
        for mode in active_modes:
            tb_sv += f"        run(2'd{mode}, {da_val}, {db_val});\n"

    tb_sv += """
        $fclose(fd);
        if (tot_errs == 0) $display("TB_PARAM PASS (%0d runs)", nruns);
        else $display("TB_PARAM FAIL (%0d errors in %0d runs)", tot_errs, nruns);
        $finish;
    end
    initial begin #5000000000; $display("TIMEOUT"); $finish; end
endmodule
"""
    tb_path.write_text(tb_sv)


def write_summary(csv_path, out_dir, hw, mat):
    rows = list(csv.DictReader(open(csv_path)))
    if not rows:
        return
    for r in rows:
        r["dA"] = round(float(r["dA"]), 2)
        r["dB"] = round(float(r["dB"]), 2)
        for k in ("M", "K", "N", "cycles", "issued_macs", "eff_macs", "b_nnz_streamed"):
            r[k] = int(r[k])

    target_shape = (mat["M"], mat["K"], mat["N"])
    grid = [r for r in rows if (r["M"], r["K"], r["N"]) == target_shape]
    ragged = [r for r in rows if (r["M"], r["K"], r["N"]) != target_shape]

    modes_present = []
    for m in ["DENSE", "TRIP", "TRGT", "TRGS"]:
        if any(r["mode"] == m for r in grid):
            modes_present.append(m)

    densA_set = sorted({r["dA"] for r in grid}, reverse=True)
    densB_set = sorted({r["dB"] for r in grid}, reverse=True)

    # Average if multiple runs
    from collections import defaultdict
    cyc_dict = defaultdict(list)
    eff_dict = defaultdict(list)
    iss_dict = defaultdict(list)
    for r in grid:
        key = (r["mode"], r["dA"], r["dB"])
        cyc_dict[key].append(r["cycles"])
        eff_dict[(r["dA"], r["dB"])].append(r["eff_macs"])
        iss_dict[key].append(r["issued_macs"])

    avg_cyc = {k: sum(v)/len(v) for k, v in cyc_dict.items()}
    avg_eff = {k: sum(v)/len(v) for k, v in eff_dict.items()}

    nrows = hw["NROWS"]
    nmul = hw["NMUL"]
    peak_macs_per_cyc = nrows * nmul

    lines = [
        "# Parameterized Trapezoid Sweep Results",
        "",
        f"- **Matrix Shape**: M={mat['M']}, K={mat['K']}, N={mat['N']}",
        f"- **Hardware Array**: {nrows} PE rows × {nmul} multipliers per row ({peak_macs_per_cyc} peak MACs/cycle)",
        f"- **Total Benchmark Runs**: {len(rows)} runs (checked vs reference matmul)",
        "",
    ]

    for m in modes_present:
        lines.append(f"## {m}: Cycles")
        lines.append("")
        lines.append("| dA \\ dB | " + " | ".join(f"{b:g}" for b in densB_set) + " |")
        lines.append("|---" * (len(densB_set) + 1) + "|")
        for a in densA_set:
            row_str = f"| **{a:g}** |"
            for b in densB_set:
                c = avg_cyc.get((m, a, b), 0)
                row_str += f" {int(round(c))} |"
            lines.append(row_str)
        lines.append("")

    if "DENSE" in modes_present:
        for m in modes_present:
            if m == "DENSE":
                continue
            lines.append(f"## {m}: Speedup vs DENSE")
            lines.append("")
            lines.append("| dA \\ dB | " + " | ".join(f"{b:g}" for b in densB_set) + " |")
            lines.append("|---" * (len(densB_set) + 1) + "|")
            for a in densA_set:
                row_str = f"| **{a:g}** |"
                for b in densB_set:
                    c_dense = avg_cyc.get(("DENSE", a, b), 0)
                    c_m = avg_cyc.get((m, a, b), 0)
                    sp = (c_dense / c_m) if c_m > 0 else 0
                    row_str += f" {sp:.1f}x |"
                lines.append(row_str)
            lines.append("")

    lines.append("## Multiplier Utilization (Effectual MACs / Peak MACs)")
    lines.append("")
    for m in modes_present:
        lines.append(f"### {m} Utilization")
        lines.append("")
        lines.append("| dA \\ dB | " + " | ".join(f"{b:g}" for b in densB_set) + " |")
        lines.append("|---" * (len(densB_set) + 1) + "|")
        for a in densA_set:
            row_str = f"| **{a:g}** |"
            for b in densB_set:
                eff = avg_eff.get((a, b), 0)
                cyc = avg_cyc.get((m, a, b), 0)
                util = (eff / (cyc * peak_macs_per_cyc)) * 100 if cyc > 0 else 0
                row_str += f" {util:.1f}% |"
            lines.append(row_str)
        lines.append("")

    lines.append("## Best Dataflow per Density Pair (Fewest Cycles)")
    lines.append("")
    lines.append("| dA \\ dB | " + " | ".join(f"{b:g}" for b in densB_set) + " |")
    lines.append("|---" * (len(densB_set) + 1) + "|")
    for a in densA_set:
        row_str = f"| **{a:g}** |"
        for b in densB_set:
            best_m = None
            best_c = float("inf")
            for m in modes_present:
                c = avg_cyc.get((m, a, b), float("inf"))
                if c < best_c:
                    best_c = c
                    best_m = m
            row_str += f" {best_m} ({int(round(best_c))}) |"
        lines.append(row_str)
    lines.append("")

    if ragged:
        lines.append("## Ragged Test Cases")
        lines.append("")
        lines.append("| Case | Shape (M,K,N) | dA | dB | Mode | Cycles | Issued MACs | Effectual MACs |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for r in ragged:
            lines.append(f"| ragged | ({r['M']},{r['K']},{r['N']}) | {r['dA']} | {r['dB']} | {r['mode']} | {r['cycles']} | {r['issued_macs']} | {r['eff_macs']} |")
        lines.append("")

    figs = sorted(p.name for p in out_dir.glob("fig_*.png"))
    if figs:
        lines.append("## Generated Plots")
        lines.append("")
        for f in figs:
            lines.append(f"![{f}]({f})")
        lines.append("")

    (out_dir / "summary.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description="Trapezoid Configurable Sweep Runner")
    parser.add_argument("--config", default=str(HERE / "parameters.json"), help="Path to parameters.json")
    parser.add_argument("--M", type=int, help="Override matrix dimension M")
    parser.add_argument("--K", type=int, help="Override matrix dimension K")
    parser.add_argument("--N", type=int, help="Override matrix dimension N")
    parser.add_argument("--densities-A", nargs="+", type=float, help="Override A densities list")
    parser.add_argument("--densities-B", nargs="+", type=float, help="Override B densities list")
    parser.add_argument("--methods", nargs="+", choices=["DENSE", "TRIP", "TRGT", "TRGS"], help="Override methods")
    parser.add_argument("--runs", type=int, help="Override number of runs per point")
    args = parser.parse_args()

    cfg_file = Path(args.config)
    if not cfg_file.exists():
        sys.exit(f"Config file not found: {cfg_file}")

    with open(cfg_file) as f:
        params = json.load(f)

    # Apply command-line overrides
    if args.M: params["matrix"]["M"] = args.M
    if args.K: params["matrix"]["K"] = args.K
    if args.N: params["matrix"]["N"] = args.N
    if args.densities_A: params["sweep"]["densities_A"] = args.densities_A
    if args.densities_B: params["sweep"]["densities_B"] = args.densities_B
    if args.methods: params["sweep"]["methods"] = args.methods
    if args.runs: params["sweep"]["num_runs_per_point"] = args.runs

    iverilog_bin = find_tool("iverilog")
    vvp_bin = find_tool("vvp")

    if not (iverilog_bin and vvp_bin):
        sys.exit("Error: iverilog or vvp not found on PATH or in conda environments.")

    # Create timestamped results directory
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = HERE / "results" / f"sweep_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Generate SV testbench
    gen_tb_path = out_dir / "tb_param.sv"
    print(f"Generating testbench: {gen_tb_path}")
    generate_testbench(params, gen_tb_path)

    # Compile with iverilog
    vvp_out = out_dir / "tb_param.vvp"
    compile_log = out_dir / "compile.log"
    sources = [
        str(HERE / "mfiu.sv"),
        str(HERE / "trip_row.sv"),
        str(HERE / "trgt_sub.sv"),
        str(HERE / "trgs_row.sv"),
        str(HERE / "trapezoid_top.sv"),
        str(gen_tb_path)
    ]
    print(f"Compiling with iverilog ({iverilog_bin})...")
    res = subprocess.run([iverilog_bin, "-g2012", "-o", str(vvp_out)] + sources,
                         capture_output=True, text=True)
    compile_log.write_text(res.stdout + res.stderr)
    if res.returncode != 0:
        print(f"Compilation FAILED:\n{res.stdout}\n{res.stderr}")
        sys.exit(1)
    print("Compilation successful.")

    # Run simulation
    print(f"Running simulation with vvp ({vvp_bin})...")
    sim_log = out_dir / "sim.log"
    res = subprocess.run([vvp_bin, str(vvp_out)], cwd=str(out_dir),
                         capture_output=True, text=True)
    sim_log.write_text(res.stdout + res.stderr)
    print(res.stdout)
    if "TB_PARAM PASS" not in res.stdout:
        print("Simulation encountered mismatches or failed to complete!")
    else:
        print("Simulation verification PASSED: All outputs match golden reference matrix multiplication.")

    csv_path = out_dir / "sweep.csv"
    if not csv_path.exists():
        sys.exit(f"Expected output file not generated: {csv_path}")

    # Generate plots using python with matplotlib/numpy
    py_plot = find_python_with_plotting()
    plot_script = HERE / "plot_results.py"
    if plot_script.exists():
        print(f"Generating publication plots using {py_plot}...")
        hw = params["hardware"]
        p_res = subprocess.run([
            py_plot, str(plot_script),
            "--csv", str(csv_path),
            "--out", str(out_dir),
            "--nrows", str(hw["NROWS"]),
            "--nmul", str(hw["NMUL"])
        ], capture_output=True, text=True)
        (out_dir / "plots.log").write_text(p_res.stdout + p_res.stderr)
        if p_res.returncode == 0:
            print("Plots generated successfully:")
            for p in sorted(out_dir.glob("fig_*.png")):
                print(f"  - {p.name}")
        else:
            print(f"Plot generation notice: {p_res.stderr.strip() or p_res.stdout.strip()}")

    # Write summary
    print(f"Writing summary report: {out_dir / 'summary.md'}")
    write_summary(csv_path, out_dir, params["hardware"], params["matrix"])

    # Update results/latest symlink / copy
    latest_dir = HERE / "results" / "latest"
    if latest_dir.exists():
        shutil.rmtree(latest_dir)
    shutil.copytree(out_dir, latest_dir)

    print(f"\nAll operations complete! Results saved in:\n  {out_dir}\n  (also mirrored to {latest_dir})")


if __name__ == "__main__":
    main()
