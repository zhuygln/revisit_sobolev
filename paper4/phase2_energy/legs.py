"""Paper IV: the energy ladder and the closure legs on a benchmark zone.

One zone (an `EjectaState` shell, the PI's local-shell reduction), one atom
built from the compact cache with the full lanthanide pattern, and the legs:

    R1    sobolev_branch     photon   Paper III's reference, grey "conserving" scale
    R1E   sobolev_branch     energy   the same histories, energy bookkeeping (exact scale)
    R2    sobolev_dmacro     energy   resolved opacity + downward macroatom   <- reference
    B1    expansion_branch   photon   grouped opacity + photon branching (Paper III's B)
    B1E   expansion_branch   energy
    B2    expansion_dmacro   energy   grouped opacity + the same macroatom   <- Gate 2 primary
    Bbin2 binned_dmacro      energy   exact-sum bins + the same macroatom
    A2    sobolev_group      energy   resolved opacity + kernel from R2 (control)
    C2    expansion_group    energy   grouped opacity + kernel from R2 (control)
    Cbin2 binned_group       energy
    C1    expansion_group    photon   kernel from R1 (Paper III's C_both)

Every energy leg is normalised by the exact geometric-series scale
(`sobolev/energy_balance.py`); the photon legs keep Paper III's grey
"conserving" scale so that R1 and B1 are the Paper III numbers. Per leg the
row records magnitudes through the DECam + 2MASS passbands, dm against R2
and against R1, the energy accounting (identity residual, comoving deposit,
work, dead ends, k-packets), and the run cost.

The adequacy ratio of the single zone -- band saturation in the zone
against the neighbouring shells -- is recorded when `--shells` lists them.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase11_observables",
          ROOT / "paper3/phase12_grid"):
    sys.path.insert(0, str(p))

from sobolev.constants import C                                   # noqa: E402
from sobolev import photometry as phot                            # noqa: E402
from sobolev import energy_balance as eb                          # noqa: E402
from sobolev.energy_packets import energy_accounting              # noqa: E402
from sobolev.abundances import ATOMIC_MASS, Z_OF                  # noqa: E402
from sobolev.ejecta import EjectaState                            # noqa: E402
from sobolev.forest_stats import band_saturation                  # noqa: E402
from forest_mc import ForestAtom, run_mc                          # noqa: E402
from redistribution import RedistributionKernel                   # noqa: E402
from observables import observe, LAM_WIN, N_SPEC, BAND3800        # noqa: E402
from grid import photometer, rss_mb                               # noqa: E402

DAY = 86400.0
SEEDS = (1, 2, 3)
NG = 32
MAX_STEPS = 1_000_000
CHAIN_MAX = 2000
LEGS = {
    "R1": dict(mode="sobolev_branch", packets="photon", scale="conserving"),
    "R1E": dict(mode="sobolev_branch", packets="energy", scale="equilibrium"),
    "R2": dict(mode="sobolev_dmacro", packets="energy", scale="equilibrium"),
    "B1": dict(mode="expansion_branch", packets="photon", scale="conserving"),
    "B1E": dict(mode="expansion_branch", packets="energy", scale="equilibrium"),
    "B2": dict(mode="expansion_dmacro", packets="energy", scale="equilibrium"),
    "Bbin2": dict(mode="binned_dmacro", packets="energy", scale="equilibrium"),
    "A2": dict(mode="sobolev_group", packets="energy", scale="equilibrium", kernel="R2"),
    "C2": dict(mode="expansion_group", packets="energy", scale="equilibrium", kernel="R2"),
    "Cbin2": dict(mode="binned_group", packets="energy", scale="equilibrium", kernel="R2"),
    "C1": dict(mode="expansion_group", packets="photon", scale="conserving", kernel="R1"),
}
LADDER = ("R1", "R1E", "R2", "B1", "B1E", "B2")
PHASE3 = ("R2", "A2", "B2", "Bbin2", "C2", "Cbin2")


def git_sha():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                              text=True, cwd=ROOT).stdout.strip()
    except Exception:
        return None


def ions_of(state, stages=("II",)):
    """[(ion name, element, stage)] for every lanthanide with mass in the state."""
    out = []
    for el in state.X:
        if el == "bulk" or el not in Z_OF:
            continue
        for st in stages:
            if f"{el} {st}" in state.f_ion or (st == "II" and not state.f_ion):
                out.append((f"{Z_OF[el]}{el}{st}", el, st))
    return out


def atom_for_zone(state, shell, stages=("II",), tau_min=1e-3, n_ion_min=1e-30):
    """The blend atom of one shell from the compact cache: (atom, {ion: n_ion})."""
    specs, n_ion = [], {}
    for ion, el, st in ions_of(state, stages):
        n = state.n_ion(el, st, shell, ATOMIC_MASS[el])
        if n > n_ion_min:
            specs.append((ion, n)); n_ion[ion] = n
    atom = ForestAtom.from_cached(specs, float(state.T_gas[shell]), float(state.t), tau_min=tau_min)
    return atom, n_ion


def run_legs(zone, atom, n, legs=LADDER, seeds=SEEDS, ng=NG, relativity="worldline",
             chain_max=CHAIN_MAX, budget_s=None, a_cut=None, verbose=True):
    """All `legs` on one zone dict (from `EjectaState.local_zone`) with one atom."""
    t0 = time.time()
    lo, hi = (float(x) for x in phot.nu_edges(*LAM_WIN, 1))
    l_core = phot.planck_luminosity(lo, hi, zone["r_core"], zone["t_core"])
    edges = phot.nu_edges(*LAM_WIN, N_SPEC)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    row = dict(zone={k: v for k, v in zone.items() if k not in ("X", "f_ion")},
               n=n, seeds=list(seeds), ng=ng, relativity=relativity, chain_max=chain_max,
               a_cut=a_cut, lam_window=list(LAM_WIN), n_spec=N_SPEC, L_core_window=l_core,
               n_opacity=int(atom.n_opacity), n_lines=int(atom.n_lines_total),
               tau_max=float(atom.op_tau.max()) if atom.n_opacity else 0.0, rss_mb_atom=rss_mb(),
               git=git_sha(), legs={}, timing={})
    nb = (C / (BAND3800[1] * 1e-8), C / (BAND3800[0] * 1e-8))
    if atom.n_opacity >= 2:
        row.update({f"band_{k}": v for k, v in band_saturation(atom, *nb).items()})

    def mc(spec, seed, **kw):
        return run_mc(atom, zone["r_core"], zone["r_out"], zone["t_exp"], lo, hi, n, spec["mode"],
                      seed=seed, t_core=zone["t_core"], relativity=relativity, max_steps=MAX_STEPS,
                      chain_max=chain_max, chain_overflow="absorb", packets=spec["packets"],
                      launch_weight="energy", wall_s=budget_s, a_cut=a_cut, **kw)

    kernels, results = {}, {}
    order = [l for l in legs if "kernel" not in LEGS[l]] + [l for l in legs if "kernel" in LEGS[l]]
    for tag in order:
        spec = LEGS[tag]
        tl = time.time()
        kw = {}
        if "kernel" in spec:
            src = spec["kernel"]
            if src not in kernels:
                if src not in results:
                    raise ValueError(f"leg {tag} needs the events of {src}; add it to the legs")
                ev = [r["events"] for r in results[src]]
                nu_in = np.concatenate([e[0] for e in ev]); nu_out = np.concatenate([e[1] for e in ev])
                w_in = np.concatenate([e[2] for e in ev])
                w_out = np.concatenate([e[3] for e in ev]) if len(ev[0]) == 4 else None
                k_lo, k_hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
                kernels[src] = RedistributionKernel.from_branching_mc(
                    nu_in, nu_out, w_in, ng, nu_lo=k_lo, nu_hi=k_hi, w_out=w_out)
            kw["kernel"] = kernels[src]
        collect = tag in ("R1", "R2")
        res = [mc(spec, s, collect_events=collect, **kw) for s in seeds]
        results[tag] = res
        o = photometer(observe(res, l_core, spec["scale"]), edges, nu_c, phot.D_40MPC)
        acc = [energy_accounting(r) for r in res]
        o["energy"] = {k: float(np.mean([a[k] for a in acc])) for k in acc[0] if k != "packets"}
        o["energy"]["packets"] = spec["packets"]
        o["events_per_packet"] = float(np.mean([r["n_events"].mean() for r in res]))
        o["reabs_per_packet"] = float(np.mean([r["n_reabs"].mean() for r in res]))
        o["n_trapped"] = int(sum(r["n_trapped"] for r in res))
        o["mode"] = spec["mode"]; o["scale"] = spec["scale"]
        o["t_wall"] = row["timing"][tag] = time.time() - tl
        row["legs"][tag] = o
        if verbose:
            e = o["energy"]
            print(f"  {tag:6s} {spec['mode']:17s} {o['t_wall']:7.1f}s ev/pkt={o['events_per_packet']:6.1f} "
                  f"esc={e['esc_frac']:.3f} core={e['core_frac']:.3f} dep_cm={e['dep_cm_frac']:+.3f} "
                  f"work={e['work_frac']:+.4f} dead={e['n_dead_end']} g={o['mags'].get('g', np.nan):.2f} "
                  f"K={o['mags'].get('K', np.nan):.2f}", flush=True)
    for ref_tag in ("R2", "R1"):
        if ref_tag not in row["legs"]:
            continue
        ref = row["legs"][ref_tag]
        for tag, o in row["legs"].items():
            o[f"dm_vs_{ref_tag}"] = phot.delta_mag(o["mags"], ref["mags"])
            o[f"dcolor_vs_{ref_tag}"] = {k: o["colors"][k] - ref["colors"][k] for k in ref["colors"]}
            o[f"dm_bol_vs_{ref_tag}"] = phot.bol_delta_mag(o["L_bol"], ref["L_bol"])
    row["t_wall"] = time.time() - t0
    row["rss_mb"] = rss_mb()
    return row


def adequacy(state, shell, neighbours, stages=("II",), tau_min=1e-3):
    """Band saturation S in the zone and in its neighbouring shells."""
    nb = (C / (BAND3800[1] * 1e-8), C / (BAND3800[0] * 1e-8))
    out = {}
    for s in [shell] + list(neighbours):
        atom, _ = atom_for_zone(state, s, stages, tau_min)
        sat = band_saturation(atom, *nb) if atom.n_opacity >= 2 else {}
        out[str(s)] = dict(shell=int(s), v_c=float(state.v_edges[s] / C), rho=float(state.rho[s]),
                           T=float(state.T_gas[s]), n_opacity=int(atom.n_opacity), **sat)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state", help="EjectaState JSON (paper4/phase1_benchmarks/P?_t?.json)")
    ap.add_argument("--shell", type=int, default=None, help="zone shell (default: the state's photospheric shell)")
    ap.add_argument("--legs", default="ladder", help="'ladder', 'phase3', 'all', or a comma list")
    ap.add_argument("--n", type=int, default=100000)
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--ng", type=int, default=NG)
    ap.add_argument("--relativity", default="worldline")
    ap.add_argument("--chain-max", type=int, default=CHAIN_MAX)
    ap.add_argument("--budget", type=float, default=None, help="wall seconds per run")
    ap.add_argument("--a-cut", type=float, default=None)
    ap.add_argument("--tau-min", type=float, default=1e-3)
    ap.add_argument("--neighbours", default="", help="comma list of shells for the adequacy ratio")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    state = EjectaState.from_json(a.state)
    shell = state.meta.get("photospheric_shell", 0) if a.shell is None else a.shell
    legs = {"ladder": LADDER, "phase3": PHASE3, "all": tuple(LEGS)}.get(a.legs, tuple(a.legs.split(",")))
    seeds = tuple(int(s) for s in a.seeds.split(","))
    zone = state.local_zone(shell)
    t0 = time.time()
    atom, n_ion = atom_for_zone(state, shell, tau_min=a.tau_min)
    print(f"{state.meta.get('name')} t={state.t / DAY:g} d shell {shell}: v={zone['v_core']:.4f}c rho={zone['rho']:.3e} "
          f"T={zone['T_gas']:.0f} K; {len(n_ion)} ions, {atom.n_lines_total} lines, {atom.n_opacity} opacity, "
          f"atom {time.time() - t0:.1f}s, rss {rss_mb():.0f} MB", flush=True)
    row = run_legs(zone, atom, a.n, legs, seeds, a.ng, a.relativity or None, a.chain_max, a.budget, a.a_cut)
    row.update(state=str(a.state), name=state.meta.get("name"), t_d=state.t / DAY, shell=shell, n_ion=n_ion,
               tau_min=a.tau_min, ions=sorted(n_ion))
    if a.neighbours:
        row["adequacy"] = adequacy(state, shell, [int(s) for s in a.neighbours.split(",")], tau_min=a.tau_min)
    out = a.out or HERE / f"legs_{state.meta.get('name')}_t{state.t / DAY:g}_s{shell}.json"
    Path(out).write_text(json.dumps(row, indent=1, default=float) + "\n")
    print(f"wrote {out} in {row['t_wall']:.0f}s")


if __name__ == "__main__":
    main()
