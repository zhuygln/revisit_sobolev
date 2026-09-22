#!/usr/bin/env python3
"""Paper B gate G1: one ion, alone, on the P1 2 d photospheric state, through
the reference macroatom, the R_ij legs at N_g = 2..32 (plus the 128-group
kernel for the event-level metric) and the epsilon grid. Everything as
preregistered in paperB/prl_gate.md; this script only runs it.

    .venv/bin/python paperB/gate1/run_gate1.py --ion 57LaII
    .venv/bin/python paperB/gate1/run_gate1.py --ion 57LaII --n 2000 --seeds 1 --build-seeds 101 --eps 0,1 --ng 4   # smoke
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper4/phase2_energy",
          ROOT / "paper3/phase11_observables", ROOT / "paper3/phase12_grid"):
    sys.path.insert(0, str(p))

from sobolev.ejecta import EjectaState                  # noqa: E402
from sobolev.abundances import ATOMIC_MASS               # noqa: E402
from forest_mc import ForestAtom                          # noqa: E402
import legs as L                                          # noqa: E402

STATE = ROOT / "paper4/phase1_benchmarks/P1_t2.json"
SHELL = 28
IONS = {"57LaII": "La", "58CeII": "Ce", "60NdII": "Nd"}
NG_GRID = (2, 4, 8, 16, 32)
NG_FINE = 128
EPS_GRID = tuple(round(0.05 * k, 2) for k in range(21))
N_PACKETS = 300_000
SEEDS = (1, 2, 3)              # evaluation: the reference photometry, the independent fine matrix, every closure's transport
BUILD_SEEDS = (101, 102, 103)  # construction: the events every R_N kernel is built from (pooled), never evaluated against themselves
TAU_MIN = 1e-3


def leg_specs(ng_grid=NG_GRID, eps_grid=EPS_GRID, ng_fine=NG_FINE, build_seeds=BUILD_SEEDS):
    """The preregistered leg dict for run_legs (every leg energy packets).
    R2build (its own seeds) supplies the events the kernels are built from;
    R2 (the evaluation seeds) is the reference and supplies the independent
    128-group matrix K128 the event-level metric is scored against; K128build
    is the in-sample fine matrix, recorded as a diagnostic only."""
    specs = {"R2build": dict(mode="sobolev_dmacro", packets="energy", scale="equilibrium", seeds=list(build_seeds), collect_events=True),
             "R2": dict(mode="sobolev_dmacro", packets="energy", scale="equilibrium")}
    for e in eps_grid:
        specs[f"E{e:.2f}"] = dict(mode="sobolev_tla", packets="energy", scale="equilibrium", eps=float(e))
    for ng in ng_grid:
        specs[f"A2_ng{ng}"] = dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2build", ng=int(ng))
    specs[f"K{ng_fine}"] = dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2", ng=int(ng_fine), kernel_only=True)
    specs[f"K{ng_fine}build"] = dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2build", ng=int(ng_fine), kernel_only=True)
    return specs


def build(ion, tau_min=TAU_MIN):
    st = EjectaState.from_json(STATE)
    zone = st.local_zone(SHELL)
    el = IONS[ion]
    n_ion = float(st.n_ion(el, "II", SHELL, ATOMIC_MASS[el]))
    atom = ForestAtom.from_cached([(ion, n_ion)], zone["T_gas"], zone["t_exp"], tau_min=tau_min)
    return st, zone, atom, n_ion


def run(ion, n=N_PACKETS, seeds=SEEDS, ng_grid=NG_GRID, eps_grid=EPS_GRID, ng_fine=NG_FINE, budget=None, out=None, verbose=True, tau_min=TAU_MIN, build_seeds=BUILD_SEEDS):
    t0 = time.time()
    st, zone, atom, n_ion = build(ion, tau_min)
    specs = leg_specs(ng_grid, eps_grid, ng_fine, build_seeds)
    if verbose:
        print(f"{ion}: n_ion {n_ion:.4e} cm^-3, {atom.n_lines_total} lines, {atom.n_opacity} opacity lines, "
              f"tau_max {atom.op_tau.max():.1f}; {len(specs)} legs x {len(seeds)} seeds x {n} packets", flush=True)
    row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=max(ng_grid), budget_s=budget, verbose=verbose)
    row.update(gate="G1", ion=ion, element=IONS[ion], n_ion=n_ion, state=str(STATE.relative_to(ROOT)), shell=SHELL,
               t_d=float(zone["t_exp"] / 86400.0), tau_min=tau_min, ng_grid=list(ng_grid), ng_fine=int(ng_fine), build_seeds=list(build_seeds),
               eps_grid=list(eps_grid), prereg="paperB/prl_gate.md", t_wall_total=time.time() - t0)
    for tag, k in row["kernels"].items():
        k["table_kb"] = k["serialized_bytes"] / 1024.0          # the whole kernel as saved: matrix + discrete exit tables
    out = Path(out) if out else HERE / f"gate1_{ion}.json"
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s", flush=True)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", required=True, choices=sorted(IONS))
    ap.add_argument("--n", type=int, default=N_PACKETS)
    ap.add_argument("--seeds", default=",".join(str(s) for s in SEEDS))
    ap.add_argument("--build-seeds", default=",".join(str(s) for s in BUILD_SEEDS))
    ap.add_argument("--ng", default=",".join(str(g) for g in NG_GRID))
    ap.add_argument("--ng-fine", type=int, default=NG_FINE)
    ap.add_argument("--eps", default=None, help="comma list; default the preregistered 21-value grid")
    ap.add_argument("--budget", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    eps = EPS_GRID if a.eps is None else tuple(float(x) for x in a.eps.split(","))
    run(a.ion, a.n, tuple(int(s) for s in a.seeds.split(",")), tuple(int(g) for g in a.ng.split(",")), eps, a.ng_fine, a.budget, a.out,
        build_seeds=tuple(int(s) for s in a.build_seeds.split(",")))


if __name__ == "__main__":
    main()
