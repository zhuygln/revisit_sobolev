#!/usr/bin/env python3
"""G3U runner (paperB/prl_gate.md "G3U"): one frozen G3 record at a time,
its reference and its affected decision legs re-transported on the common
evaluation seed set 1..12 with the per-seed magnitudes kept
(legs.py `mags_per_seed`), the trained operators preserved -- rebuilt
from the build seeds or loaded from the saved kernels and interpolated
with the frozen lambda, then compared entry by entry with the matrix in the
frozen record (`operator_match`; a mismatch is gray).

    .venv/bin/python paperB/gate3u/run_g3u.py --record gate3_P_P3d_60NdII.json
    .venv/bin/python paperB/gate3u/run_g3u.py --record gate3_partc_p1blend.json
    ... --seeds 1,2 --n 2000 --legs Arec_ng16                                     # smoke
"""
import argparse
import importlib.util as _ilu
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1", ROOT / "paperB/gate3", ROOT / "paper4/phase2_energy", ROOT / "paper3", ROOT / "paper3/phase5_mixture",
          ROOT / "paper2/phase1"):
    sys.path.insert(0, str(p))

import legs as L                                                    # noqa: E402
from redistribution import RedistributionKernel                     # noqa: E402
import run_gate3 as R3                                              # noqa: E402  (states, kernel paths, interpolation; env G3_KDIR)
import states as S                                                  # noqa: E402
import operators3 as O3                                             # noqa: E402
_spec = _ilu.spec_from_file_location("paperB_gate3u_analyse", HERE / "analyse.py")
U = _ilu.module_from_spec(_spec); _spec.loader.exec_module(U)

SEEDS = tuple(U.PREREG["seeds"])
ODIR = Path(os.environ.get("G3U_OUT", HERE))


def out_path(name):
    return ODIR / name.replace("gate3_", "g3u_")


def _finish(row, frozen, name, affected, t0, verbose):
    match = {}
    for tag in affected:
        if tag in row["kernels"] and tag in frozen["kernels"]:
            match[tag] = float(np.abs(np.asarray(row["kernels"][tag]["R"]) - np.asarray(frozen["kernels"][tag]["R"])).max())
    for tag, k in row["kernels"].items():
        k["table_kb"] = k["serialized_bytes"] / 1024.0
    row.update(gate="G3U", frozen_record=name, affected_legs=list(affected), operator_match=match, prereg="paperB/prl_gate.md#g3u",
               t_wall_total=time.time() - t0)
    out = out_path(name); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    (out.parent / f"{out.stem}.done").write_text(f"{row['git']} {row['t_wall_total']:.0f}s {row['n']}\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s; operator match {match}", flush=True)
    return row


def run_state_record(name, frozen, affected, seeds=SEEDS, n=None, verbose=True, tau_min=R3.TAU_MIN):
    t0 = time.time()
    axis, value, ion = frozen["axis"], frozen["value"], frozen["ion"]
    n = int(frozen["n"]) if n is None else int(n)
    st = S.build_state(axis, value, ion, tau_min); sup = R3.support(ion); ng_t = frozen["ng_t"]
    base = dict(packets="energy", scale="equilibrium")
    specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(frozen["build_seeds"]), collect_events=True, **base),
             "R2": dict(mode="sobolev_dmacro", collect_events=True, **base)}
    for tag in affected:
        if tag.startswith("Arec_ng"):
            specs[tag] = dict(mode="sobolev_group", kernel="R2build", ng=int(tag.split("_ng")[1]), **base)
        elif tag == f"Afix_ng{ng_t}":
            anchor = RedistributionKernel.load(R3.kernel_path("ref", ion, ng_t))
            anchor.metadata["injected"] = dict(kind="anchor", built_at="ref", ng=int(ng_t))
            specs[tag] = dict(mode="sobolev_group", kernel_obj=anchor, source="anchor:ref", **base)
        elif tag in (f"Aint_ng{ng_t}", f"AintM_ng{ng_t}"):
            interp = R3.interpolation_for(axis, value, ion, ng_t)
            k = interp["whole" if tag.startswith("Aint_") else "matrix"]
            k.metadata["injected"] = dict(kind=k.metadata["transform"]["kind"], lam=interp["lam"], endpoints=interp["endpoints"])
            specs[tag] = dict(mode="sobolev_group", kernel_obj=k, source="interp_whole" if tag.startswith("Aint_") else "interp_matrix", **base)
        else:
            raise ValueError(f"unknown affected leg {tag}")
    specs["K128"] = dict(mode="sobolev_group", kernel="R2", ng=128, kernel_only=True, transform=O3.support_diag(), **base)
    specs["K128build"] = dict(mode="sobolev_group", kernel="R2build", ng=128, kernel_only=True, **base)
    if verbose:
        print(f"{name}: {ion} {axis}:{frozen['label']}  {len(specs)} legs; seeds {list(seeds)} x {n}; affected {affected}", flush=True)
    row = L.run_legs(st["zone"], st["atom"], n, legs=specs, seeds=tuple(seeds), ng=128, verbose=verbose, kernel_range=(sup["nu_lo"], sup["nu_hi"]))
    row.update(axis=axis, value=value, label=frozen["label"], ion=ion, ng_t=int(ng_t), build_seeds=list(frozen["build_seeds"]))
    return _finish(row, frozen, name, affected, t0, verbose)


def run_blend_record(name, frozen, affected, seeds=SEEDS, n=None, verbose=True, tau_min=R3.TAU_MIN):
    t0 = time.time()
    part = frozen["part"]; n = int(frozen["n"]) if n is None else int(n)
    base = dict(packets="energy", scale="equilibrium")
    specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(frozen["build_seeds"]), collect_events=True, **base),
             "R2": dict(mode="sobolev_dmacro", collect_events=True, **base)}
    kernel_range = None
    if part == "c":
        from sobolev.ejecta import EjectaState
        st = EjectaState.from_json(R3.STATE); zone = st.local_zone(R3.SHELL)
        atom, _ = L.atom_for_zone(st, R3.SHELL, tau_min=tau_min)
        for tag in affected:
            assert tag.startswith("Arec_ng"), tag
            specs[tag] = dict(mode="sobolev_group", kernel="R2build", ng=int(tag.split("_ng")[1]), **base)
    else:
        from mixture import composition_weights
        from forest_mc import ForestAtom
        sup = json.loads(R3.SUPPORT.read_text())["blend3"]; kernel_range = (sup["nu_lo"], sup["nu_hi"])
        zone, _, _ = S.ref_zone_and_density(next(iter(R3.IONS)))
        dens = {ion: S.ref_zone_and_density(ion)[1] for ion in R3.IONS}
        atom = ForestAtom.from_cached([(ion, dens[ion]) for ion in R3.IONS], zone["T_gas"], zone["t_exp"], tau_min=tau_min)
        for tag in affected:
            ng = int(tag.split("_ng")[1])
            if tag.startswith("Adirect_ng"):
                specs[tag] = dict(mode="sobolev_group", kernel="R2build", ng=ng, **base)
            elif tag.startswith("Amix_ng"):
                ks = [RedistributionKernel.load(R3.kernel_path("blend3", ion, ng, "single")) for ion in R3.IONS]
                km = RedistributionKernel.mix(ks, composition_weights(atom, ks[0].edges),
                                              metadata=dict(injected=dict(kind="mix", weights_rule="composition_weights", ng=ng)))
                specs[tag] = dict(mode="sobolev_group", kernel_obj=km, source="mix:composition_weights", **base)
            else:
                raise ValueError(f"unknown affected leg {tag}")
    specs["K128"] = dict(mode="sobolev_group", kernel="R2", ng=128, kernel_only=True, transform=O3.support_diag(), **base)
    specs["K128build"] = dict(mode="sobolev_group", kernel="R2build", ng=128, kernel_only=True, **base)
    if verbose:
        print(f"{name}: part {part}  {len(specs)} legs; seeds {list(seeds)} x {n}; affected {affected}", flush=True)
    row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=128, verbose=verbose, kernel_range=kernel_range)
    row.update(part=part, build_seeds=list(frozen["build_seeds"]))
    return _finish(row, frozen, name, affected, t0, verbose)


def run_record(name, seeds=SEEDS, n=None, legs=None, verbose=True):
    frozen = U.frozen_records()[name]
    affected = legs or sorted({a["leg"] for a in U.affected({name: frozen})})
    if not affected:
        raise ValueError(f"{name}: no affected leg")
    if "part" in frozen:
        return run_blend_record(name, frozen, affected, seeds, n, verbose)
    return run_state_record(name, frozen, affected, seeds, n, verbose)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", required=True, help="the frozen G3 record's file name, e.g. gate3_P_P3d_60NdII.json")
    ap.add_argument("--seeds", default=",".join(str(s) for s in SEEDS))
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--legs", default=None, help="override the affected legs (smoke only)")
    ap.add_argument("--skip-done", action="store_true")
    a = ap.parse_args()
    op = out_path(a.record)
    if a.skip_done and (op.parent / f"{op.stem}.done").exists():
        print(f"skip {op.name}: done"); return
    run_record(a.record, tuple(int(s) for s in a.seeds.split(",")), a.n, a.legs.split(",") if a.legs else None)


if __name__ == "__main__":
    main()
