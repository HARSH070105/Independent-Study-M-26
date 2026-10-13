#!/usr/bin/env python3
"""Compile + run every Trapezoid testbench and save all results.

Usage:  python3 run_all.py            (from anywhere; needs iverilog + vvp on PATH)
Output: results/<YYYYmmdd_HHMMSS>/ and a copy at results/latest/
          unit_test.log      tb_trip output (dense + TrIP unit test)
          sweep.log          tb_top output (every run, PASS/FAIL)
          sweep.csv          raw numbers: mode,dA,dB,M,K,N,cycles,issued_macs,eff_macs,b_nnz_streamed
          summary.md         cycle / speedup-vs-dense / utilization tables per mode
          fig_*.png          3-D surfaces, heat-maps and paper-style bar charts (needs matplotlib + numpy;
                             produced by plot_results.py, skipped with a notice if they are missing)
Exit code is 0 only if every check passed.
"""
import csv, datetime, shutil, subprocess, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = ["mfiu.sv", "trip_row.sv", "trgt_sub.sv", "trgs_row.sv", "trapezoid_top.sv"]
NROWS, NMUL = 4, 16            # keep in sync with tb_top.sv (used only for utilization)
MODES = ["DENSE", "TRIP", "TRGT", "TRGS"]


def sh(cmd, cwd, log):
    """Run cmd, tee to log file, return (returncode, text)."""
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    text = p.stdout + p.stderr
    log.write_text(text)
    return p.returncode, text


def main():
    if not (shutil.which("iverilog") and shutil.which("vvp")):
        sys.exit("iverilog/vvp not found. Install Icarus Verilog (e.g. `sudo apt install iverilog`).")

    out = HERE / "results" / datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out.mkdir(parents=True)
    src = [str(HERE / f) for f in SRC]
    ok = True

    # ---- 1. unit tests (individual engines: Dense, TrIP, TrGT, TrGS) ----
    print("[1/3] unit tests (tb_dense, tb_trip, tb_trgt, tb_trgs) ...", flush=True)
    unit_specs = [
        ("DENSE", "tb_dense", ["mfiu.sv", "trip_row.sv", "tb_dense.sv"], "TB_DENSE PASS"),
        ("TRIP",  "tb_trip",  ["mfiu.sv", "trip_row.sv", "tb_trip.sv"],  "TB_TRIP PASS"),
        ("TRGT",  "tb_trgt",  ["trgt_sub.sv", "tb_trgt.sv"],             "TB_TRGT PASS"),
        ("TRGS",  "tb_trgs",  ["trgs_row.sv", "tb_trgs.sv"],             "TB_TRGS PASS"),
    ]
    unit_ok = True
    unit_logs = []
    compile_logs = []
    for label, name, files, token in unit_specs:
        vvp = out / f"{name}.vvp"
        clog = out / f"compile_{name}.log"
        ulog = out / f"{name}.log"
        rc, ctxt = sh(["iverilog", "-g2012", "-o", str(vvp), *[str(HERE / f) for f in files]], out, clog)
        compile_logs.append(ctxt)
        if rc == 0:
            _, utxt = sh(["vvp", str(vvp)], out, ulog)
            passed = token in utxt
            unit_logs.append(f"=== {label} ({name}) ===\n" + utxt)
        else:
            passed = False
            unit_logs.append(f"=== {label} ({name}) COMPILE FAILED ===\n" + ctxt)
        print(f"       {label:5s} ({name}): {'PASS' if passed else 'FAIL'}", flush=True)
        unit_ok &= passed
    (out / "compile_unit.log").write_text("".join(compile_logs))
    (out / "unit_test.log").write_text("\n".join(unit_logs))
    ok &= unit_ok


    # ---- 2. full sweep (all four dataflows); tb_top writes sweep.csv into cwd ----
    print("[2/3] density sweep, 4 dataflows (tb_top, ~2 min) ...", flush=True)
    rc, _ = sh(["iverilog", "-g2012", "-o", str(out / "tb_top.vvp"), *src, str(HERE / "tb_top.sv")],
               out, out / "compile_sweep.log")
    if rc == 0:
        _, txt = sh(["vvp", str(out / "tb_top.vvp")], out, out / "sweep.log")
        sweep_ok = "TB_TOP PASS" in txt
    else:
        sweep_ok = False
    print("      ", "PASS" if sweep_ok else "FAIL"); ok &= sweep_ok

    # ---- 3. summary ----
    csv_path = out / "sweep.csv"
    if csv_path.exists():
        print("[3/3] plots + summary ...", flush=True)
        rc, txt = sh([sys.executable, str(HERE / "plot_results.py"), "--csv", str(csv_path), "--out", str(out),
                      "--nrows", str(NROWS), "--nmul", str(NMUL)], out, out / "plots.log")
        print("      ", "figures written" if rc == 0 else "plotting skipped/failed (see plots.log): " + txt.strip().splitlines()[-1])
        write_summary(csv_path, out, unit_ok, sweep_ok)

    for f in out.glob("*.vvp"):
        f.unlink()
    latest = HERE / "results" / "latest"
    if latest.exists():
        shutil.rmtree(latest)
    shutil.copytree(out, latest)
    print(f"\nResults saved in {out}\n(also copied to {latest})")
    print("ALL PASSED" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)


def write_summary(csv_path, out, unit_ok, sweep_ok):
    rows = list(csv.DictReader(csv_path.open()))
    for r in rows:
        for k in ("dA", "dB"):
            r[k] = round(float(r[k]), 2)
        for k in ("M", "K", "N", "cycles", "issued_macs", "eff_macs", "b_nnz_streamed"):
            r[k] = int(r[k])
    grid = [r for r in rows if (r["M"], r["K"], r["N"]) == (16, 32, 32)]
    ragged = [r for r in rows if (r["M"], r["K"], r["N"]) != (16, 32, 32)]
    dens = sorted({r["dA"] for r in grid}, reverse=True)
    cyc = {(r["mode"], r["dA"], r["dB"]): r["cycles"] for r in grid}
    eff = {(r["dA"], r["dB"]): r["eff_macs"] for r in grid}
    iss = {(r["mode"], r["dA"], r["dB"]): r["issued_macs"] for r in grid}
    bnz = {(r["mode"], r["dA"], r["dB"]): r["b_nnz_streamed"] for r in grid}

    L = ["# Trapezoid RTL sweep results", "",
         f"- unit test: **{'PASS' if unit_ok else 'FAIL'}**, sweep ({len(rows)} runs, all checked vs reference matmul): "
         f"**{'PASS' if sweep_ok else 'FAIL'}**",
         f"- config: {NROWS} PE rows x {NMUL} multipliers, M=16 K=32 N=32 (ragged case listed at the end)",
         "- rows of every table = density of A, columns = density of B", ""]

    def table(title, fn, fmt):
        L.extend([f"## {title}", "", "| dA \\ dB | " + " | ".join(f"{d:g}" for d in dens) + " |",
                  "|---" * (len(dens) + 1) + "|"])
        for a in dens:
            L.append(f"| **{a:g}** | " + " | ".join(fmt(fn(a, b)) for b in dens) + " |")
        L.append("")

    for m in MODES:
        table(f"{m}: cycles", lambda a, b, m=m: cyc[(m, a, b)], lambda v: f"{v}")
    for m in MODES[1:]:
        table(f"{m}: speedup vs Dense", lambda a, b, m=m: cyc[("DENSE", a, b)] / cyc[(m, a, b)],
              lambda v: f"{v:.1f}x")
    for m in MODES:
        table(f"{m}: multiplier utilization (effectual MACs / (cycles x {NROWS * NMUL}))",
              lambda a, b, m=m: eff[(a, b)] / (cyc[(m, a, b)] * NROWS * NMUL), lambda v: f"{100 * v:.1f}%")
    L.append("## Best dataflow per (dA, dB) point (fewest cycles)\n")
    L += ["| dA \\ dB | " + " | ".join(f"{d:g}" for d in dens) + " |", "|---" * (len(dens) + 1) + "|"]
    for a in dens:
        cells = []
        for b in dens:
            best = min(MODES, key=lambda m: cyc[(m, a, b)])
            cells.append(f"{best} ({cyc[(best, a, b)]})")
        L.append(f"| **{a:g}** | " + " | ".join(cells) + " |")
    L.append("")
    L.append("## Gustavson traffic proxy: B nonzeros streamed (TrGT / TrGS), e.g. at dB=0.2\n")
    L += ["| dA | TrGT | TrGS |", "|---|---|---|"]
    for a in dens:
        L.append(f"| {a:g} | {bnz[('TRGT', a, 0.2)]} | {bnz[('TRGS', a, 0.2)]} |")
    L.append("")
    if ragged:
        r0 = ragged[0]
        L += [f"## Ragged shape M={r0['M']} K={r0['K']} N={r0['N']}, dA=dB=0.3", "",
              "| mode | cycles | issued MACs | effectual MACs |", "|---|---|---|---|"]
        L += [f"| {r['mode']} | {r['cycles']} | {r['issued_macs']} | {r['eff_macs']} |" for r in ragged]
        L.append("")
    figs = sorted(p.name for p in out.glob("fig_*.png"))
    if figs:
        L += ["## Figures", ""] + [f"![{f}]({f})" for f in figs] + [""]
    (out / "summary.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
