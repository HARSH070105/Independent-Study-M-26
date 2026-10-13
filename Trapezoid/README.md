# Trapezoid (ISCA'24) – scaled-down parameterized RTL model

Icarus Verilog 12, SystemVerilog-2012 subset. Run: `make` (unit test + 104-run sweep, writes `sweep.csv`).

## Files
| File | Paper section | Role |
|---|---|---|
| `mfiu.sv` | III-B, Fig. 11 | Multi-fiber intersection unit: pairwise AND, prefix-sum slot assignment, routing metadata |
| `trip_row.sv` | III-A/B | PE row for Dense IP (mode 0) and TrIP (mode 1): dynamic B-column packing, row-subset packing, segmented reduction, banked output buffer |
| `trgt_sub.sv` | III-C | TrGT subrow (HS×HS): fiber-head merge tree, MS multipliers per subrow, passes when A row has > MS nonzeros |
| `trgs_row.sv` | III-D | TrGS PE row (HS×MS, HS×D): spatial B rows, W nonzeros/cycle, N-tiling by #multipliers |
| `trapezoid_top.sv` | III | Dynamic row dispatch to idle units, C write-back, aggregated counters |
| `tb_trip.sv`, `tb_top.sv` | – | Reference-checked testbenches; `tb_top` sweeps density (A,B ∈ {1, .5, .2, .1, .05}) for all four modes + ragged shapes |

## Parameters (`trapezoid_top`)
`NROWS` PE rows (paper 128), `NMUL` multipliers/row (paper 128), `KT` K-tile (= NMUL), `RA`/`CB` rows/columns intersected at once (paper 4/4), `W` cache word (paper 16), `MAXM/K/N` memory dims, `DW=16` signed integer data, `ACC=40`.
Current testbench: NROWS=4, NMUL=16, M=16, K=N=32.

## Interface
Drive `A_flat`/`B_flat` (dense row-major, element (i,j) at `(i*MAXK+j)*DW`), set `mode`, `M,K,N`; pulse `rst` then `start`; wait `done`; read `C_flat`, `cyc_total`, `cnt_macs`, `cnt_engine_cyc`, `cnt_brows`, `cnt_bnnz`. **Pulse `rst` before every run** (it clears counters).

## Fidelity and deliberate simplifications
Faithful: dataflow loop nests, MFIU algorithm, dynamic packing limits (≤ NMUL effectual MACs, ≤ 4×4 intersections), 4-way subrow split in TrGT, TrGS wide-word streaming, per-mode multiplier/bandwidth limits.
Simplified (affects absolute cycle counts, not trends):
- Integer arithmetic instead of FP32 (exact checking).
- Benes distribution networks and the merge-reduction tree are behavioral (index gather / segmented sum); no pipeline latency modeled.
- No global cache / HBM model yet: `cnt_brows`/`cnt_bnnz` count fiber fetches and streamed nonzeros as a traffic proxy but there are no miss stalls or bandwidth limit.
- Matrices are dense arrays; each dataflow derives its sparse view (bitmask or fiber list) combinationally, instead of reading CSR/bitmask-compressed storage.
- Scheduler uses blocking variables (simulation-only, not synthesizable as written).

## Running
```
sudo apt install iverilog          # or: brew install icarus-verilog
python3 run_all.py                 # compile + run everything, save results (~1 min)
make unit | make sweep             # individual runs
```
`run_all.py` writes `results/<timestamp>/` (and `results/latest/`): `unit_test.log`, `sweep.log`, `sweep.csv`,
`summary.md` (cycles, speedup vs dense, utilization, best-dataflow tables), `cycles_vs_density.png`.
