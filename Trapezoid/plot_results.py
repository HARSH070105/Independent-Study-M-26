#!/usr/bin/env python3
"""Plots for the Trapezoid RTL sweep (reads sweep.csv written by tb_top.sv).

  fig_3d_cycles_surface.png   3-D surface: cycles over (A density, B density), one panel per dataflow,
                              with the same surface projected on the floor as a heat-map
  fig_3d_best_dataflow.png    3-D bars: cycles of the fastest dataflow per (dA, dB), coloured by winner
  fig_heatmaps.png            annotated heat-maps: cycles per dataflow, best-dataflow map, speedup vs dense
  fig_perf_by_dataflow.png    paper Fig. 19/22/25 style: grouped bars, normalised to the best dataflow,
                              grouped by D x D / MS x D / MS x MS / HS x D / HS x MS / HS x HS, + gmean
  fig_utilization.png         paper Fig. 23 style: compute utilisation (effectual MACs / peak) per workload

Usage: python3 plot_results.py --csv sweep.csv --out DIR [--nrows 4 --nmul 16]
"""
import argparse, csv, math
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize, to_rgb
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)

MODES = ["DENSE", "TRIP", "TRGT", "TRGS"]
LABEL = {"DENSE": "Dense IP (TPU-like)", "TRIP": "TrIP", "TRGT": "TrGT", "TRGS": "TrGS"}
SHORT = {"DENSE": "Dense", "TRIP": "TrIP", "TRGT": "TrGT", "TRGS": "TrGS"}
# validated categorical slots 1,2,3,7 (validate_palette.js: all checks pass; aqua contrast WARN -> labels + table view shipped)
COLOR = {"DENSE": "#2a78d6", "TRIP": "#eb6834", "TRGT": "#1baf7a", "TRGS": "#6250d6"}
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
       "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]       # blue 100 -> 700
CMAP = LinearSegmentedColormap.from_list("seq_blue", SEQ)
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "text.color": INK, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.linewidth": 0.8,
})

# Representative workload classes (paper Sec. IV): D = dense, MS = moderately sparse, HS = hyper-sparse.
# The simulated matrices are tiny (16x32x32), so "HS" here means density <= 5% (a few dozen nonzeros).
WORKLOADS = [
    ("D×D",   [(1.0, 1.0)]),
    ("MS×D",  [(0.5, 1.0), (0.3, 1.0), (0.2, 1.0)]),
    ("MS×MS", [(0.5, 0.5), (0.3, 0.3), (0.2, 0.2)]),
    ("HS×D",  [(0.05, 1.0), (0.03, 1.0), (0.03, 0.7)]),
    ("HS×MS", [(0.05, 0.5), (0.03, 0.3), (0.03, 0.2)]),
    ("HS×HS", [(0.05, 0.05), (0.03, 0.05), (0.03, 0.03)]),
]


def luminance(rgb):
    def lin(c):
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(c) for c in rgb[:3])
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def ink_on(rgb):
    L = luminance(rgb)
    return "white" if 1.05 / (L + 0.05) > (L + 0.05) / (luminance(to_rgb(INK)) + 0.05) else INK


def load(csv_path):
    rows = []
    for r in csv.DictReader(open(csv_path)):
        rows.append(dict(mode=r["mode"], dA=round(float(r["dA"]), 2), dB=round(float(r["dB"]), 2),
                         shape=(int(r["M"]), int(r["K"]), int(r["N"])), cycles=int(r["cycles"]),
                         eff=int(r["eff_macs"]), issued=int(r["issued_macs"])))
    shape = Counter(r["shape"] for r in rows).most_common(1)[0][0]     # main grid (ragged runs excluded)
    rows = [r for r in rows if r["shape"] == shape]
    cyc = {m: {} for m in MODES}
    eff = {}
    for r in rows:
        cyc[r["mode"]][(r["dA"], r["dB"])] = r["cycles"]
        eff[(r["dA"], r["dB"])] = r["eff"]
    dens = sorted({r["dA"] for r in rows}, reverse=True)
    return cyc, eff, dens, shape


# ----------------------------------------------------------------------------- 3-D
def style3d(ax, n, dens):
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.988, 0.988, 0.984, 1.0))
        try:
            axis._axinfo["grid"].update(color=GRID, linewidth=0.6)
        except Exception:
            pass
    ticks = [f"{d:g}" for d in dens]
    ax.set_xticks(range(n)); ax.set_xticklabels(ticks, fontsize=7, color=MUTED)
    ax.set_yticks(range(n)); ax.set_yticklabels(ticks, fontsize=7, color=MUTED)
    ax.tick_params(axis="z", labelsize=7, colors=MUTED)
    ax.set_xlabel("A density", labelpad=6, fontsize=9, color=INK2)
    ax.set_ylabel("B density", labelpad=6, fontsize=9, color=INK2)
    ax.set_zlabel("cycles", labelpad=2, fontsize=9, color=INK2)


def fig_3d_surface(cyc, dens, cfg, out):
    n = len(dens)
    X, Y = np.meshgrid(np.arange(n), np.arange(n))
    Z = {m: np.array([[cyc[m][(dens[ia], dens[ib])] for ia in range(n)] for ib in range(n)], float) for m in MODES}
    zmax = max(z.max() for z in Z.values())
    fig = plt.figure(figsize=(14, 11))
    axes = []
    for i, m in enumerate(MODES):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        axes.append(ax)
        ax.plot_surface(X, Y, Z[m], cmap=CMAP, vmin=0, vmax=zmax, rstride=1, cstride=1,
                        linewidth=0.5, edgecolor=(1, 1, 1, 0.55), antialiased=True)
        ax.contourf(X, Y, Z[m], zdir="z", offset=0, levels=np.linspace(0, zmax, 14), cmap=CMAP, vmin=0, vmax=zmax)
        ax.set_zlim(0, zmax * 1.04)
        ax.set_box_aspect((1, 1, 0.75))
        ax.view_init(elev=26, azim=124)
        style3d(ax, n, dens)
        ax.set_title(LABEL[m], fontsize=12, pad=-6)
    fig.suptitle("Cycles vs. operand density — one surface per dataflow (floor = same data as a heat-map)",
                 fontsize=13, y=0.985)
    fig.text(0.5, 0.955, cfg + " · shared colour/height scale · axes = density grid of the sweep (1.0 = dense)",
             ha="center", fontsize=9, color=INK2)
    cax = fig.add_axes([0.93, 0.30, 0.014, 0.40])
    fig.colorbar(ScalarMappable(Normalize(0, zmax), CMAP), cax=cax, label="cycles")
    fig.subplots_adjust(left=0.0, right=0.92, top=0.89, bottom=0.02, wspace=0.0, hspace=0.04)
    fig.savefig(out, dpi=150); plt.close(fig)


def fig_3d_best(cyc, dens, cfg, out):
    n = len(dens)
    best = np.zeros((n, n), int); height = np.zeros((n, n))
    for ia in range(n):
        for ib in range(n):
            c = [cyc[m][(dens[ia], dens[ib])] for m in MODES]
            best[ia, ib] = int(np.argmin(c)); height[ia, ib] = min(c)
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_subplot(projection="3d")
    xs, ys = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    cols = [COLOR[MODES[k]] for k in best.ravel()]
    w = 0.72
    ax.bar3d(xs.ravel() - w / 2, ys.ravel() - w / 2, np.zeros(n * n), w, w, height.ravel(), color=cols, shade=True)
    ax.set_box_aspect((1, 1, 0.65)); ax.view_init(elev=30, azim=124)
    style3d(ax, n, dens)
    ax.set_zlim(0, height.max() * 1.05)
    ax.legend(handles=[Patch(color=COLOR[m], label=LABEL[m]) for m in MODES], loc="upper left",
              bbox_to_anchor=(0.0, 0.90), frameon=True, facecolor=SURF, edgecolor=AXIS, title="fastest dataflow")
    ax.set_title("Cycles of the fastest dataflow at each (A density, B density)", fontsize=13, pad=14)
    fig.text(0.5, 0.90, cfg, ha="center", fontsize=9, color=INK2)
    fig.subplots_adjust(left=0.0, right=1.0, top=0.95, bottom=0.02)
    fig.savefig(out, dpi=150); plt.close(fig)


# ----------------------------------------------------------------------------- heat-maps
def fig_heatmaps(cyc, dens, cfg, out):
    n = len(dens)
    Z = {m: np.array([[cyc[m][(dens[ia], dens[ib])] for ib in range(n)] for ia in range(n)], float) for m in MODES}
    allc = np.stack([Z[m] for m in MODES])
    best = allc.argmin(0); bestc = allc.min(0)
    speed = Z["DENSE"] / bestc
    norm = Normalize(0, allc.max())
    fig, axs = plt.subplots(2, 3, figsize=(17, 11))
    ticks = [f"{d:g}" for d in dens]

    def frame(ax, title):
        ax.set_xticks(range(n)); ax.set_xticklabels(ticks, fontsize=8)
        ax.set_yticks(range(n)); ax.set_yticklabels(ticks, fontsize=8)
        ax.set_xticks(np.arange(-.5, n, 1), minor=True); ax.set_yticks(np.arange(-.5, n, 1), minor=True)
        ax.grid(which="minor", color=SURF, linewidth=2); ax.grid(which="major", visible=False)
        ax.tick_params(which="minor", length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlabel("B density"); ax.set_ylabel("A density")
        ax.set_title(title, fontsize=12, loc="left")

    def annotate(ax, M, cmap, nrm, fmt):
        for i in range(n):
            for j in range(n):
                ax.text(j, i, fmt(i, j), ha="center", va="center", fontsize=8, color=ink_on(cmap(nrm(M[i, j]))))

    for k, m in enumerate(MODES):
        ax = axs[k // 3, k % 3]
        ax.imshow(Z[m], cmap=CMAP, norm=norm)
        annotate(ax, Z[m], CMAP, norm, lambda i, j, m=m: f"{int(Z[m][i, j])}")
        frame(ax, f"{LABEL[m]} — cycles")
    ax = axs[1, 1]                                                      # best-dataflow map
    ax.imshow(best, cmap=ListedColormap([COLOR[m] for m in MODES]), vmin=-.5, vmax=3.5)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{SHORT[MODES[best[i, j]]]}\n{int(bestc[i, j])}", ha="center", va="center", fontsize=7.5,
                    color=ink_on(to_rgb(COLOR[MODES[best[i, j]]])))
    frame(ax, "Fastest dataflow (cycles)")
    ax = axs[1, 2]                                                      # speedup over dense
    snorm = Normalize(1, speed.max())
    im = ax.imshow(speed, cmap=CMAP, norm=snorm)
    annotate(ax, speed, CMAP, snorm, lambda i, j: f"{speed[i, j]:.1f}×")
    frame(ax, "Best dataflow: speedup vs Dense IP")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle("Trapezoid RTL — performance across the sparsity plane (rows: A density, columns: B density)",
                 fontsize=14, y=0.985)
    fig.text(0.5, 0.945, cfg + " · cycle panels share one scale (darker = more cycles)",
             ha="center", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=150); plt.close(fig)


# ----------------------------------------------------------------------------- paper-style grouped bars
def layout(workloads, gap=0.8):
    xs, seps, centers, x = [], [], [], 0.0
    for ci, (cat, items) in enumerate(workloads):
        if ci > 0:
            seps.append(xs[-1] + (gap + 1) / 2); x = xs[-1] + 1 + gap
        first = len(xs)
        for _ in items:
            xs.append(x); x += 1
        centers.append((cat, (xs[first] + xs[-1]) / 2))
        x = xs[-1]
    return xs, seps, centers


def grouped_bars(vals, workloads, ylabel, ylim, out, title, cfg, gmean=None):
    flat = [w for _, items in workloads for w in items]
    xs, seps, centers = layout(workloads)
    fig, ax = plt.subplots(figsize=(14.5, 4.9))
    bw, k = 0.185, len(MODES)
    for j, m in enumerate(MODES):
        ax.bar([x + (j - (k - 1) / 2) * 0.2 for x in xs], vals[m], width=bw, color=COLOR[m], zorder=3)
    labels = [f"{a:g}/{b:g}" for a, b in flat]
    xt = list(xs)
    if gmean is not None:
        xg = xs[-1] + 1.8
        ax.axvspan(xg - 0.55, xg + 0.55, color="#ebebe8", zorder=0)
        for j, m in enumerate(MODES):
            ax.bar(xg + (j - (k - 1) / 2) * 0.2, gmean[m], width=bw, color=COLOR[m], zorder=3)
        xt.append(xg); labels.append("gmean")
        seps.append(xs[-1] + 0.9)
        ax.set_xlim(xs[0] - 0.7, xg + 0.7)
    else:
        ax.set_xlim(xs[0] - 0.7, xs[-1] + 0.7)
    ax.set_xticks(xt); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    for s in seps:
        ax.axvline(s, color=MUTED, linestyle=(0, (4, 3)), linewidth=0.9, zorder=1)
    for cat, c in centers:
        ax.text(c, -0.30, cat, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=10, fontweight="bold",
                color=INK2)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontweight="bold", fontsize=9)
    ax.yaxis.grid(True, zorder=0); ax.set_axisbelow(True)
    ax.legend(handles=[Patch(color=COLOR[m], label=LABEL[m]) for m in MODES], loc="lower center",
              bbox_to_anchor=(0.5, 1.02), ncol=len(MODES), frameon=True, facecolor="#efefec", edgecolor=AXIS)
    ax.set_xlabel("workload = A density / B density", fontsize=8, color=MUTED, labelpad=34)
    fig.text(0.01, 0.005, cfg, fontsize=8, color=MUTED)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, y=1.0, color=INK)
    fig.subplots_adjust(left=0.065, right=0.995, top=0.86, bottom=0.30)
    fig.savefig(out, dpi=150); plt.close(fig)


def fig_perf(cyc, workloads, cfg, out):
    flat = [w for _, items in workloads for w in items]
    best = [min(cyc[m][w] for m in MODES) for w in flat]
    vals = {m: [b / cyc[m][w] for b, w in zip(best, flat)] for m in MODES}
    gm = {m: math.exp(sum(math.log(v) for v in vals[m]) / len(flat)) for m in MODES}
    grouped_bars(vals, workloads, "Normalized Performance\n(1/cycles, best dataflow = 1.0)", (0, 1.08), out,
                 "Performance of each Trapezoid dataflow, normalised to the best dataflow per workload "
                 "(cf. paper Figs. 19 / 22 / 25)", cfg, gmean=gm)
    return gm


def fig_util(cyc, eff, workloads, nrows, nmul, cfg, out):
    flat = [w for _, items in workloads for w in items]
    vals = {m: [eff[w] / (cyc[m][w] * nrows * nmul) for w in flat] for m in MODES}
    grouped_bars(vals, workloads, "Compute Utilization\n(effectual MACs / peak)", (0, 1.0), out,
                 "Compute utilisation per workload — only effectual (nonzero×nonzero) MACs count as useful work "
                 "(cf. paper Fig. 23)", cfg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nrows", type=int, default=4)
    ap.add_argument("--nmul", type=int, default=16)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cyc, eff, dens, (M, K, N) = load(a.csv)
    cfg = f"{a.nrows} PE rows × {a.nmul} multipliers · M={M}, K={K}, N={N} · INT16 data"
    wl = [(c, [w for w in items if w in eff]) for c, items in WORKLOADS]
    wl = [(c, items) for c, items in wl if items]
    fig_3d_surface(cyc, dens, cfg, out / "fig_3d_cycles_surface.png")
    fig_3d_best(cyc, dens, cfg, out / "fig_3d_best_dataflow.png")
    fig_heatmaps(cyc, dens, cfg, out / "fig_heatmaps.png")
    gm = fig_perf(cyc, wl, cfg, out / "fig_perf_by_dataflow.png")
    fig_util(cyc, eff, wl, a.nrows, a.nmul, cfg, out / "fig_utilization.png")
    print("gmean normalized performance:", {LABEL[m]: round(v, 3) for m, v in gm.items()})
    for f in sorted(out.glob("fig_*.png")):
        print("wrote", f)


if __name__ == "__main__":
    main()
