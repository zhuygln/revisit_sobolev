#!/usr/bin/env python3
"""Paper B, the R2M robustness check (post-pass, paperB/prl_gate.md "R2M
robustness"): one ion alone on the P1 2 d state through the
radiation-field-driven macroatom (Lucy's macroatom with internal upward
jumps under an imposed W = 1/2 B_nu(T_zone), `sobolev_macro`) as the
reference, the 8-group operator rebuilt from its events, and G1's
downward-trained 8-group operator scored against it. Not a gate.

    .venv/bin/python paperB/r2m/run_r2m.py --ion 57LaII
    .venv/bin/python paperB/r2m/run_r2m.py --ion 57LaII --n 2000 --seeds 1 --build-seeds 101   # smoke
"""
import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1", ROOT / "paper4/phase2_energy"):
    sys.path.insert(0, str(p))

import legs as L                                                   # noqa: E402
from run_gate1 import build, IONS, STATE, SHELL, TAU_MIN, N_PACKETS, SEEDS, BUILD_SEEDS   # noqa: E402

NG = 8
NG_FINE = 128
W_DIL = 0.5


def leg_specs(ng=NG, ng_fine=NG_FINE, build_seeds=BUILD_SEEDS, W=W_DIL):
    """The preregistered leg dict: G1's reference pair and its 8-group
    operator, the radiation-field-driven pair and the operator rebuilt from
    it, and the two fine matrices of the new reference."""
    return {
        "R2build": dict(mode="sobolev_dmacro", packets="energy", scale="equilibrium", seeds=list(build_seeds), collect_events=True),
        "R2": dict(mode="sobolev_dmacro", packets="energy", scale="equilibrium"),
        "R2Mbuild": dict(mode="sobolev_macro", packets="energy", scale="equilibrium", seeds=list(build_seeds), collect_events=True, macro_kw=dict(W=W)),
        "R2M": dict(mode="sobolev_macro", packets="energy", scale="equilibrium", collect_events=True, macro_kw=dict(W=W)),
        f"A2_ng{ng}": dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2build", ng=int(ng)),
        f"A2M_ng{ng}": dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2Mbuild", ng=int(ng)),
        f"K{ng_fine}M": dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2M", ng=int(ng_fine), kernel_only=True),
        f"K{ng_fine}Mbuild": dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2Mbuild", ng=int(ng_fine), kernel_only=True),
    }


def run(ion, n=N_PACKETS, seeds=SEEDS, build_seeds=BUILD_SEEDS, ng=NG, ng_fine=NG_FINE, W=W_DIL, budget=None, out=None,
        verbose=True, tau_min=TAU_MIN):
    t0 = time.time()
    st, zone, atom, n_ion = build(ion, tau_min)
    specs = leg_specs(ng, ng_fine, build_seeds, W)
    if verbose:
        print(f"{ion}: n_ion {n_ion:.4e} cm^-3, {atom.n_lines_total} lines, {atom.n_opacity} opacity lines; "
              f"{len(specs)} legs, {n} packets, W = {W}", flush=True)
    row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=int(ng), budget_s=budget, verbose=verbose)
    row.update(gate="R2M", ion=ion, element=IONS[ion], n_ion=n_ion, state=str(STATE.relative_to(ROOT)), shell=SHELL,
               t_d=float(zone["t_exp"] / 86400.0), tau_min=tau_min, ng_check=int(ng), ng_fine=int(ng_fine), macro_W=float(W),
               build_seeds=list(build_seeds), prereg="paperB/prl_gate.md#r2m-robustness", t_wall_total=time.time() - t0)
    for tag, k in row["kernels"].items():
        k["table_kb"] = k["serialized_bytes"] / 1024.0
    out = Path(out) if out else HERE / f"r2m_{ion}.json"
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s (rss {row['rss_mb']:.0f} MB)", flush=True)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", required=True, choices=sorted(IONS))
    ap.add_argument("--n", type=int, default=N_PACKETS)
    ap.add_argument("--seeds", default=",".join(str(s) for s in SEEDS))
    ap.add_argument("--build-seeds", default=",".join(str(s) for s in BUILD_SEEDS))
    ap.add_argument("--ng", type=int, default=NG)
    ap.add_argument("--ng-fine", type=int, default=NG_FINE)
    ap.add_argument("--budget", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run(a.ion, a.n, tuple(int(s) for s in a.seeds.split(",")), tuple(int(s) for s in a.build_seeds.split(",")),
        a.ng, a.ng_fine, budget=a.budget, out=a.out)


if __name__ == "__main__":
    main()
