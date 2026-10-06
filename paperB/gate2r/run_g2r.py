#!/usr/bin/env python3
"""G2R runner (paperB/prl_gate.md "G2R"): one ion at one block count N_g per
process -- G2's reference pair rerun seed for seed (R2build on the build
seeds, R2 on the evaluation seeds), the independent fine matrices, and the
32 orderings of the 128 fine groups block-coarsened at N_g and transported
on the shared 128-group tables. Ordering 0 is G2's own L128_ng{N_g} leg.

    .venv/bin/python paperB/gate2r/run_g2r.py --ion 58CeII --k 16
    .venv/bin/python paperB/gate2r/run_g2r.py --ion 57LaII --k 4 --n 2000 --seeds 1 --build-seeds 101 --n-scrambled 2   # smoke

After the transport the event-level loss of every ordering against K128 is
computed and stored, and the 128x128 matrix of every scrambled ordering is
dropped from the record (32 of them per record would be ~5 MB of text each
and carry nothing the stored loss and the permutation seed do not).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1", ROOT / "paperB/gate2r", ROOT / "paper4/phase2_energy", ROOT / "paper3"):
    sys.path.insert(0, str(p))

import legs as L                                                           # noqa: E402
from run_gate1 import build, IONS, STATE, SHELL, TAU_MIN, SEEDS, BUILD_SEEDS, NG_FINE   # noqa: E402
import operators_perm as OPP                                               # noqa: E402
import importlib.util as _ilu                                              # noqa: E402
_spec = _ilu.spec_from_file_location("paperB_gate1_analyse", ROOT / "paperB/gate1/analyse.py")
A = _ilu.module_from_spec(_spec); _spec.loader.exec_module(A)               # m_event

K_GRID = (2, 4, 8, 16, 32)
N_PACKETS = {"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000}   # G2's
ODIR = Path(os.environ.get("G2R_OUT", HERE))


def record_path(ion, k):
    return ODIR / f"g2r_{ion}_k{k}.json"


def leg_specs(k, orderings, build_seeds=BUILD_SEEDS, ng_fine=NG_FINE):
    base = dict(packets="energy", scale="equilibrium")
    specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(build_seeds), collect_events=True, **base),
             "R2": dict(mode="sobolev_dmacro", collect_events=True, **base)}
    for pid, perm, seed in orderings:
        specs[f"P{pid:03d}_ng{k}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine),
                                          transform=OPP.permuted_local(k, perm, pid, seed), **base)
    specs[f"K{ng_fine}"] = dict(mode="sobolev_group", kernel="R2", ng=int(ng_fine), kernel_only=True, **base)
    specs[f"K{ng_fine}build"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), kernel_only=True, **base)
    return specs


def run(ion, k, n=None, seeds=SEEDS, build_seeds=BUILD_SEEDS, n_scrambled=OPP.N_SCRAMBLED, ng_fine=NG_FINE, budget=None,
        out=None, verbose=True, tau_min=TAU_MIN):
    n = N_PACKETS[ion] if n is None else int(n)
    t0 = time.time()
    st, zone, atom, n_ion = build(ion, tau_min)
    orderings = OPP.orderings(n_scrambled, ng_fine)
    specs = leg_specs(k, orderings, build_seeds, ng_fine)
    if verbose:
        print(f"{ion} N_g = {k}: n_ion {n_ion:.4e} cm^-3, {atom.n_opacity} opacity lines; "
              f"{len(specs)} legs x {len(seeds)} seeds x {n} packets", flush=True)
    row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=int(ng_fine), budget_s=budget, verbose=verbose)
    fine = row["kernels"][f"K{ng_fine}"]
    for tag, kr in row["kernels"].items():
        kr["table_kb"] = kr["serialized_bytes"] / 1024.0
        if tag.startswith("P"):
            kr["event"] = A.m_event(kr, fine)                   # against the independent fine matrix, before R is dropped
            if kr["transform"]["perm_id"] != 0:
                del kr["R"]
    row.update(gate="G2R", ion=ion, element=IONS[ion], n_ion=n_ion, state=str(STATE.relative_to(ROOT)), shell=SHELL,
               t_d=float(zone["t_exp"] / 86400.0), tau_min=tau_min, k=int(k), k_grid=list(K_GRID), ng_fine=int(ng_fine),
               build_seeds=list(build_seeds), n_scrambled=int(n_scrambled), perm_seed0=OPP.PERM_SEED0,
               perm_seeds=[s for _, _, s in orderings], perms=[p.tolist() for _, p, _ in orderings],
               prereg="paperB/prl_gate.md#g2r", t_wall_total=time.time() - t0)
    out = Path(out) if out else record_path(ion, k)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    (out.parent / f"{out.stem}.done").write_text(f"{row['git']} {row['t_wall_total']:.0f}s {row['n']}\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s (rss {row['rss_mb']:.0f} MB)", flush=True)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", required=True, choices=sorted(IONS))
    ap.add_argument("--k", type=int, required=True, choices=K_GRID)
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--seeds", default=",".join(str(s) for s in SEEDS))
    ap.add_argument("--build-seeds", default=",".join(str(s) for s in BUILD_SEEDS))
    ap.add_argument("--n-scrambled", type=int, default=OPP.N_SCRAMBLED)
    ap.add_argument("--budget", type=float, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--skip-done", action="store_true")
    a = ap.parse_args()
    rp = Path(a.out) if a.out else record_path(a.ion, a.k)
    if a.skip_done and (rp.parent / f"{rp.stem}.done").exists():
        print(f"skip {rp.name}: done"); return
    run(a.ion, a.k, a.n, tuple(int(s) for s in a.seeds.split(",")), tuple(int(s) for s in a.build_seeds.split(",")),
        n_scrambled=a.n_scrambled, budget=a.budget, out=a.out)


if __name__ == "__main__":
    main()
