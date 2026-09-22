#!/usr/bin/env python3
"""Paper B gate G2 (paperB/prl_gate.md "G2"): one ion alone on the P1 2 d
state; the reference pair of G1 rerun on the same seeds (so the reference
is G1's, seed for seed), and the three operator families on the shared
128-group exit tables: the local control at N_g = 8, the rank-k global
factorisation for k in K_GRID, and the exit-table truncations for f in
F_GRID. Everything as preregistered; this script only runs it.

    .venv/bin/python paperB/gate2/run_gate2.py --ion 58CeII
    .venv/bin/python paperB/gate2/run_gate2.py --ion 57LaII --n 2000 --seeds 1 --build-seeds 101 --k 1 --f 0.9   # smoke
"""
import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1", ROOT / "paper4/phase2_energy", ROOT / "paper3"):
    sys.path.insert(0, str(p))

import legs as L                                                           # noqa: E402
from run_gate1 import build, IONS, STATE, SHELL, TAU_MIN, SEEDS, BUILD_SEEDS, NG_FINE   # noqa: E402
import operators as OP                                                     # noqa: E402

K_GRID = (1, 2, 4, 8, 16, 32)
NG_LOCAL = (2, 4, 8, 16, 32)          # the local family, transported on the 128-group tables
F_GRID = (0.1, 0.2, 0.5, 0.9, 0.99, 0.999)
NG_CONTROL = 8
N_PACKETS = {"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000}   # each ion's final G1 record


def leg_specs(k_grid=K_GRID, f_grid=F_GRID, ng_control=NG_CONTROL, ng_fine=NG_FINE, build_seeds=BUILD_SEEDS,
              ng_local=NG_LOCAL):
    base = dict(packets="energy", scale="equilibrium")
    specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(build_seeds), collect_events=True, **base),
             "R2": dict(mode="sobolev_dmacro", collect_events=True, **base),
             f"A2_ng{ng_control}": dict(mode="sobolev_group", kernel="R2build", ng=int(ng_control), **base)}
    for n in ng_local:
        specs[f"L{ng_fine}_ng{n}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), transform=OP.local(n), **base)
    for k in k_grid:
        specs[f"G_k{k}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), transform=OP.nmf_rank(k), **base)
    for f in f_grid:
        specs[f"T_f{f:g}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), transform=OP.truncate(f), **base)
    specs[f"K{ng_fine}"] = dict(mode="sobolev_group", kernel="R2", ng=int(ng_fine), kernel_only=True, **base)
    specs[f"K{ng_fine}build"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), kernel_only=True, **base)
    return specs


def run(ion, n=None, seeds=SEEDS, build_seeds=BUILD_SEEDS, k_grid=K_GRID, f_grid=F_GRID, ng_control=NG_CONTROL, ng_fine=NG_FINE,
        budget=None, out=None, verbose=True, tau_min=TAU_MIN, ng_local=NG_LOCAL):
    n = N_PACKETS[ion] if n is None else int(n)
    t0 = time.time()
    st, zone, atom, n_ion = build(ion, tau_min)
    specs = leg_specs(k_grid, f_grid, ng_control, ng_fine, build_seeds, ng_local)
    if verbose:
        print(f"{ion}: n_ion {n_ion:.4e} cm^-3, {atom.n_lines_total} lines, {atom.n_opacity} opacity lines; "
              f"{len(specs)} legs x {len(seeds)} seeds x {n} packets", flush=True)
    row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=int(ng_fine), budget_s=budget, verbose=verbose)
    row.update(gate="G2", ion=ion, element=IONS[ion], n_ion=n_ion, state=str(STATE.relative_to(ROOT)), shell=SHELL,
               t_d=float(zone["t_exp"] / 86400.0), tau_min=tau_min, k_grid=list(k_grid), f_grid=list(f_grid),
               ng_control=int(ng_control), ng_fine=int(ng_fine), build_seeds=list(build_seeds), ng_local=list(ng_local),
               prereg="paperB/prl_gate.md#g2", t_wall_total=time.time() - t0)
    for tag, k in row["kernels"].items():
        k["table_kb"] = k["serialized_bytes"] / 1024.0
    out = Path(out) if out else HERE / f"gate2_{ion}.json"
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s (rss {row['rss_mb']:.0f} MB)", flush=True)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", required=True, choices=sorted(IONS))
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--seeds", default=",".join(str(s) for s in SEEDS))
    ap.add_argument("--build-seeds", default=",".join(str(s) for s in BUILD_SEEDS))
    ap.add_argument("--k", default=",".join(str(k) for k in K_GRID))
    ap.add_argument("--f", default=",".join(str(f) for f in F_GRID))
    ap.add_argument("--ng-local", default=",".join(str(n) for n in NG_LOCAL))
    ap.add_argument("--budget", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run(a.ion, a.n, tuple(int(s) for s in a.seeds.split(",")), tuple(int(s) for s in a.build_seeds.split(",")),
        tuple(int(k) for k in a.k.split(",")), tuple(float(f) for f in a.f.split(",")), budget=a.budget, out=a.out,
        ng_local=tuple(int(n) for n in a.ng_local.split(",")))


if __name__ == "__main__":
    main()
