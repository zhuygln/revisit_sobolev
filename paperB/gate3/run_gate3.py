#!/usr/bin/env python3
"""Paper B gate G3 runner (paperB/prl_gate.md "G3"): one (axis, state, ion)
at a time on the frozen per-ion support, with the two tests kept apart --
the anchor R_16 built at theta0 transported here (transfer) and the
operators rebuilt here (existence) -- and, at the interior point of an axis,
the whole-operator interpolation built from the two endpoints' kernels.
Parts (b) and (c) are the composability and realistic-mixture runs.

    .venv/bin/python paperB/gate3/run_gate3.py --ion 58CeII --axis ref        # theta0: the anchor
    .venv/bin/python paperB/gate3/run_gate3.py --ion 58CeII --axis T          # every state on the axis, interior last
    .venv/bin/python paperB/gate3/run_gate3.py --ion 58CeII --part b          # La+Ce+Nd composability (ion ignored)
    .venv/bin/python paperB/gate3/run_gate3.py --part c                       # the 13-ion P1 blend
    ... --n 2000 --seeds 1 --build-seeds 101 --ng 4,16                         # smoke
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1", ROOT / "paper4/phase2_energy", ROOT / "paper3", ROOT / "paper3/phase5_mixture",
          ROOT / "paper2/phase1"):
    sys.path.insert(0, str(p))

import legs as L                                                    # noqa: E402
from redistribution import RedistributionKernel                     # noqa: E402
from run_gate1 import IONS, STATE, SHELL, TAU_MIN, SEEDS, BUILD_SEEDS, NG_FINE, EPS_GRID   # noqa: E402
import states as S                                                  # noqa: E402
import operators3 as O3                                             # noqa: E402

NG_GRID = (2, 4, 8, 16, 32)
NG_T = 16                       # the one controlled representation size of the transfer and interpolation legs
N_PACKETS = {"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000}
N_BLEND = 1_000_000
import os
KDIR = Path(os.environ.get("G3_KDIR", HERE / "kernels"))        # kernels (gitignored); overridable for smoke runs and tests
ODIR = Path(os.environ.get("G3_OUT", HERE))                      # records
SUPPORT = HERE / "gate3_support.json"


def support(ion):
    return json.loads(SUPPORT.read_text())["per_ion"][ion]


def record_path(axis, label, ion, n=None):
    """A gray-remedy rerun at a raised packet count gets its own file; the
    original stays (G1's convention: gate1_60NdII_n1e6.json)."""
    ODIR.mkdir(parents=True, exist_ok=True)
    suffix = "" if n is None or n == N_PACKETS[ion] else f"_n{n:.0e}".replace("e+0", "e")
    return ODIR / (f"gate3_ref_{ion}{suffix}.json" if axis == "ref" else f"gate3_{axis}_{label}_{ion}{suffix}.json")


def kernel_path(label, ion, ng=NG_T, tag="rec"):
    return KDIR / f"{label}_{ion}_{tag}_ng{ng}.npz"


def leg_specs(anchor=None, interp=None, ng_grid=NG_GRID, ng_t=NG_T, ng_fine=NG_FINE, build_seeds=BUILD_SEEDS, save_to=None):
    base = dict(packets="energy", scale="equilibrium")
    specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(build_seeds), collect_events=True, **base),
             "R2": dict(mode="sobolev_dmacro", collect_events=True, **base)}
    for n in ng_grid:
        specs[f"Arec_ng{n}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(n), **base)
    if save_to is not None:
        specs[f"Arec_ng{ng_t}"]["save"] = str(save_to)
    if anchor is not None:
        specs[f"Afix_ng{ng_t}"] = dict(mode="sobolev_group", kernel_obj=anchor, source="anchor:ref", **base)
    if interp is not None:
        specs[f"Aint_ng{ng_t}"] = dict(mode="sobolev_group", kernel_obj=interp["whole"], source="interp_whole", **base)
        specs[f"AintM_ng{ng_t}"] = dict(mode="sobolev_group", kernel_obj=interp["matrix"], source="interp_matrix", **base)
    specs[f"K{ng_fine}"] = dict(mode="sobolev_group", kernel="R2", ng=int(ng_fine), kernel_only=True, transform=O3.support_diag(), **base)
    specs[f"K{ng_fine}build"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), kernel_only=True, **base)
    return specs


def interpolation_for(axis, value, ion, ng=NG_T):
    """At the interior point: the whole operator (and the matrix-only
    diagnostic) from the two bracketing endpoints' saved kernels."""
    a = S.AXES[axis]
    if value != a["interior"]:
        return None
    va, vb = a["bracket"]
    la = S.label(axis, va) if va != a["ref"] else "ref"; lb = S.label(axis, vb) if vb != a["ref"] else "ref"
    pa, pb = kernel_path(la, ion, ng), kernel_path(lb, ion, ng)
    if not (pa.exists() and pb.exists()):
        raise FileNotFoundError(f"interior {S.label(axis, value)} needs the endpoint kernels {pa.name} and {pb.name}: run the endpoints first")
    ka, kb = RedistributionKernel.load(pa), RedistributionKernel.load(pb)
    ca, cb, c = S.coord_of(axis, va, ion), S.coord_of(axis, vb, ion), S.coord_of(axis, value, ion)
    lam = (c - ca) / (cb - ca)
    ends = [la, lb]
    whole = O3.interpolate_whole(ka, kb, lam, coord=a["coord"], endpoints=ends)
    tables_from = "a" if la == "ref" else "b"
    matrix = O3.interpolate_matrix(ka, kb, lam, tables_from, coord=a["coord"], endpoints=ends)
    return dict(whole=whole, matrix=matrix, lam=lam, coord=a["coord"], endpoints=ends, coord_values=[ca, cb, c])


def run_state(axis, value, ion, n=None, seeds=SEEDS, build_seeds=BUILD_SEEDS, ng_grid=NG_GRID, ng_t=NG_T, ng_fine=NG_FINE,
              out=None, verbose=True, tau_min=TAU_MIN):
    n = N_PACKETS[ion] if n is None else int(n)
    t0 = time.time()
    st = S.build_state(axis, value, ion, tau_min)
    sup = support(ion)
    anchor = None
    if axis != "ref":
        pa = kernel_path("ref", ion, ng_t)
        if not pa.exists():
            raise FileNotFoundError(f"the anchor {pa.name} does not exist: run --axis ref first")
        anchor = RedistributionKernel.load(pa)
        anchor.metadata["injected"] = dict(kind="anchor", built_at="ref", ng=int(ng_t))
    interp = None if axis == "ref" else interpolation_for(axis, value, ion, ng_t)
    if interp:
        for k in ("whole", "matrix"):
            interp[k].metadata["injected"] = dict(kind=interp[k].metadata["transform"]["kind"], lam=interp["lam"], endpoints=interp["endpoints"])
    save_to = kernel_path(st["label"], ion, ng_t)
    if n != N_PACKETS[ion] and save_to.exists():
        save_to = None                                          # a gray-remedy rerun keeps the preregistered-count kernel on disk
    specs = leg_specs(anchor, interp, ng_grid, ng_t, ng_fine, build_seeds, save_to=save_to)
    if verbose:
        a = st["atom"]
        print(f"{ion} {axis}:{st['label']}  T_gas {st['T_gas']:.0f} t_core {st['t_core']:.0f} n {st['n_ion']:.3e}; "
              f"{a.n_opacity} opacity lines; support {sup['nu_lo']:.3e}..{sup['nu_hi']:.3e}; {len(specs)} legs x {len(seeds)} x {n}", flush=True)
    row = L.run_legs(st["zone"], st["atom"], n, legs=specs, seeds=tuple(seeds), ng=int(ng_fine), verbose=verbose,
                     kernel_range=(sup["nu_lo"], sup["nu_hi"]))
    row.update(gate="G3", axis=axis, value=value, label=st["label"], ion=ion, element=IONS[ion], coord=st["coord"],
               coord_value=st["coord_value"], support=sup, n_ion=st["n_ion"], T_gas=st["T_gas"], t_core=st["t_core"],
               day=st["day"], shell=st["shell"], ng_grid=list(ng_grid), ng_t=int(ng_t), ng_fine=int(ng_fine),
               build_seeds=list(build_seeds), tau_min=tau_min, prereg="paperB/prl_gate.md#g3",
               interpolation=None if not interp else dict(lam=interp["lam"], coord=interp["coord"], endpoints=interp["endpoints"],
                                                          coord_values=interp["coord_values"]),
               t_wall_total=time.time() - t0)
    for tag, k in row["kernels"].items():
        k["table_kb"] = k["serialized_bytes"] / 1024.0
    out = Path(out) if out else record_path(axis, st["label"], ion, n)
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    (out.parent / f"{out.stem}.done").write_text(f"{row['git']} {row['t_wall_total']:.0f}s {n}\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s (rss {row['rss_mb']:.0f} MB)", flush=True)
    return row


def run_axis(axis, ion, skip_done=False, **kw):
    """Every state on the axis, the interior last. `skip_done` resumes a
    chain after an infrastructure failure: a state whose completion marker
    exists is not rerun (the marker is written only after the record)."""
    if not kernel_path("ref", ion, kw.get("ng_t", NG_T)).exists():
        run_state("ref", None, ion, **kw)
    a = S.AXES[axis]
    order = [v for v in a["grid"] if v != a["interior"]] + [a["interior"]]
    out = []
    for v in order:
        rp = record_path(axis, S.label(axis, v), ion, kw.get("n"))
        if skip_done and (rp.parent / f"{rp.stem}.done").exists():
            print(f"skip {rp.name}: done", flush=True); continue
        out.append(run_state(axis, v, ion, **kw))
    return out


# ---- part (b): species composability; part (c): the realistic mixture ----
def run_part_b(n=N_BLEND, seeds=SEEDS, build_seeds=BUILD_SEEDS, ng_grid=NG_GRID, ng_fine=NG_FINE, out=None, verbose=True, tau_min=TAU_MIN):
    from mixture import composition_weights
    from forest_mc import ForestAtom
    t0 = time.time()
    sup = json.loads(SUPPORT.read_text())["blend3"]
    zone, _, sh = S.ref_zone_and_density(IONS.__iter__().__next__())
    dens = {ion: S.ref_zone_and_density(ion)[1] for ion in IONS}
    base = dict(packets="energy", scale="equilibrium")
    # the single-ion kernels on the three-ion support, one build run per ion
    per_ion = {}
    for ion in IONS:
        st = S.build_state("ref", None, ion, tau_min)
        specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(build_seeds), collect_events=True, **base)}
        for ng in ng_grid:
            specs[f"Kb_ng{ng}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng), kernel_only=True,
                                      save=str(kernel_path("blend3", ion, ng, "single")), **base)
        n_ion_pk = N_PACKETS[ion]
        L.run_legs(st["zone"], st["atom"], n_ion_pk, legs=specs, seeds=tuple(seeds), ng=int(ng_fine), verbose=verbose,
                   kernel_range=(sup["nu_lo"], sup["nu_hi"]))
        per_ion[ion] = {ng: RedistributionKernel.load(kernel_path("blend3", ion, ng, "single")) for ng in ng_grid}
    atom = ForestAtom.from_cached([(ion, dens[ion]) for ion in IONS], zone["T_gas"], zone["t_exp"], tau_min=tau_min)
    specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(build_seeds), collect_events=True, **base),
             "R2": dict(mode="sobolev_dmacro", collect_events=True, **base)}
    for ng in ng_grid:
        specs[f"Adirect_ng{ng}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng), **base)
        ks = [per_ion[ion][ng] for ion in IONS]
        w = composition_weights(atom, ks[0].edges)
        km = RedistributionKernel.mix(ks, w, metadata=dict(injected=dict(kind="mix", weights_rule="composition_weights", ng=int(ng))))
        specs[f"Amix_ng{ng}"] = dict(mode="sobolev_group", kernel_obj=km, source="mix:composition_weights", **base)
    specs[f"K{ng_fine}"] = dict(mode="sobolev_group", kernel="R2", ng=int(ng_fine), kernel_only=True, transform=O3.support_diag(), **base)
    specs[f"K{ng_fine}build"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), kernel_only=True, **base)
    row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=int(ng_fine), verbose=verbose, kernel_range=(sup["nu_lo"], sup["nu_hi"]))
    row.update(gate="G3", part="b", ions=list(IONS), n_ion=dens, support=sup, ng_grid=list(ng_grid), ng_fine=int(ng_fine),
               build_seeds=list(build_seeds), shell=sh, prereg="paperB/prl_gate.md#g3", t_wall_total=time.time() - t0)
    for tag, k in row["kernels"].items():
        k["table_kb"] = k["serialized_bytes"] / 1024.0
    out = Path(out) if out else ODIR / "gate3_partb_blend3.json"
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    (out.parent / f"{out.stem}.done").write_text(f"{row['git']} {row['t_wall_total']:.0f}s {n}\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s", flush=True)
    return row


def run_part_c(n=N_BLEND, seeds=SEEDS, build_seeds=BUILD_SEEDS, ng_grid=NG_GRID, ng_fine=NG_FINE, eps_grid=EPS_GRID, out=None,
               verbose=True, tau_min=TAU_MIN):
    from sobolev.ejecta import EjectaState
    t0 = time.time()
    st = EjectaState.from_json(STATE); zone = st.local_zone(SHELL)
    atom, n_ion = L.atom_for_zone(st, SHELL, tau_min=tau_min)
    base = dict(packets="energy", scale="equilibrium")
    specs = {"R2build": dict(mode="sobolev_dmacro", seeds=list(build_seeds), collect_events=True, **base),
             "R2": dict(mode="sobolev_dmacro", collect_events=True, **base)}
    for ng in ng_grid:
        specs[f"Arec_ng{ng}"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng), **base)
    for e in eps_grid:
        specs[f"E{e:.2f}"] = dict(mode="sobolev_tla", eps=float(e), **base)
    specs[f"K{ng_fine}"] = dict(mode="sobolev_group", kernel="R2", ng=int(ng_fine), kernel_only=True, transform=O3.support_diag(), **base)
    specs[f"K{ng_fine}build"] = dict(mode="sobolev_group", kernel="R2build", ng=int(ng_fine), kernel_only=True, **base)
    row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=int(ng_fine), verbose=verbose)
    row.update(gate="G3", part="c", ions=sorted(n_ion), n_ion=n_ion, ng_grid=list(ng_grid), ng_fine=int(ng_fine),
               eps_grid=list(eps_grid), build_seeds=list(build_seeds), shell=SHELL, state=str(STATE.relative_to(ROOT)),
               prereg="paperB/prl_gate.md#g3", t_wall_total=time.time() - t0)
    for tag, k in row["kernels"].items():
        k["table_kb"] = k["serialized_bytes"] / 1024.0
    out = Path(out) if out else ODIR / "gate3_partc_p1blend.json"
    out.write_text(json.dumps(row, indent=1, default=float) + "\n")
    (out.parent / f"{out.stem}.done").write_text(f"{row['git']} {row['t_wall_total']:.0f}s {n}\n")
    if verbose:
        print(f"wrote {out} in {row['t_wall_total']:.0f} s", flush=True)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", choices=sorted(IONS), default=None)
    ap.add_argument("--axis", choices=["ref"] + list(S.AXES), default=None)
    ap.add_argument("--state", default=None, help="one grid value on the axis (default: the whole axis)")
    ap.add_argument("--part", choices=["b", "c"], default=None)
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--seeds", default=",".join(str(s) for s in SEEDS))
    ap.add_argument("--build-seeds", default=",".join(str(s) for s in BUILD_SEEDS))
    ap.add_argument("--ng", default=",".join(str(g) for g in NG_GRID))
    ap.add_argument("--eps", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--skip-done", action="store_true", help="resume: skip states whose completion marker exists")
    a = ap.parse_args()
    seeds = tuple(int(s) for s in a.seeds.split(",")); bseeds = tuple(int(s) for s in a.build_seeds.split(","))
    ng = tuple(int(g) for g in a.ng.split(","))
    if a.part == "b":
        run_part_b(a.n or N_BLEND, seeds, bseeds, ng, out=a.out); return
    if a.part == "c":
        eps = EPS_GRID if a.eps is None else tuple(float(x) for x in a.eps.split(","))
        run_part_c(a.n or N_BLEND, seeds, bseeds, ng, eps_grid=eps, out=a.out); return
    if not a.ion or not a.axis:
        ap.error("--ion and --axis (or --part) are required")
    kw = dict(n=a.n, seeds=seeds, build_seeds=bseeds, ng_grid=ng)
    if a.axis == "ref":
        run_state("ref", None, a.ion, out=a.out, **kw)
    elif a.state is not None:
        v = float(a.state) if a.axis != "P" else int(a.state)
        run_state(a.axis, v, a.ion, out=a.out, **kw)
    else:
        run_axis(a.axis, a.ion, skip_done=a.skip_done, **kw)


if __name__ == "__main__":
    main()
