"""Paper III item 5: does a real kilonova cross the cancellation boundary?

The question is no longer "how inaccurate is grouped opacity?" but:

    can an approximate transport model appear accurate at one epoch
    only because two large errors happen to cancel there?

§4.32 found the grouped-closure error changes sign with band saturation, §4.34b
reproduced that boundary in a controlled forest, and §4.35 showed the two
constituent approximations are not additive and can flip the total's sign. A
kilonova sweeps its own saturation: homologous expansion gives rho ~ t^-3, so
tau ~ n t ~ t^-2, and the ejecta walk through the control parameter as they
evolve.

If the trajectory crosses the boundary, the closure is too opaque early, appears
excellent at some epoch, and is too bright late -- while the constituent errors
stay individually large throughout. That is a modelling systematic, not a
numerical curiosity: a closure validated at one epoch fails, with opposite sign,
at another.

Normalization is the ASTROPHYSICAL standard (`from_conditions`): fixed
composition and density history, each ion given whatever optical depths physics
hands it. This is deliberately NOT the controlled standard used for the
cross-ion audit -- different question, different normalization (§4.33).

Legs per epoch, all against that epoch's own `sobolev_branch`:
    A  sobolev_group     exact opacity, grouped redistribution
    B  expansion_branch  grouped opacity, exact A*beta
    C  expansion_group   both -- the practical closure
    C' binned_group      the exact-sum grouped variant

Usage: python trajectory.py [--ion 58CeII] [--n 400000]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase0_reference"):
    sys.path.insert(0, str(p))
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import R_CORE, R_OUT, T_CORE, nu_of
from redistribution import RedistributionKernel
from reference import BANDS, SEEDS

from sobolev.constants import C
from sobolev.forest_stats import band_saturation
from sobolev.normalization import from_conditions

DATA = ROOT / "data"
DAY = 86400.0
BAND = (3800.0, 3955.0)

# A lanthanide-rich ejecta trajectory. rho at 1 d chosen so the ion passes
# through the boundary within the observable window; homologous thereafter.
RHO_1D = 2.0e-17          # g cm^-3
X_LAN = 0.1               # lanthanide mass fraction
ION_FRAC = 1.0
T_GAS_1D = 5000.0         # K; cools as t^-1/2, a standard kilonova scaling
EPOCHS = (0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0)
A_OF = {"58CeII": 140.1, "57LaII": 138.9, "60NdII": 144.2, "59PrII": 140.9}


def state(t_d):
    """Homologous ejecta state at epoch t_d (days)."""
    return dict(t_exp=t_d * DAY,
                rho=RHO_1D * t_d ** -3,
                T_gas=T_GAS_1D * t_d ** -0.5,
                r_core=R_CORE * t_d, r_out=R_OUT * t_d)


def epoch(ion, t_d, n, ng=32):
    st = state(t_d)
    lev, tr = DATA / f"{ion}_levels_calib.txt", DATA / f"{ion}_transitions_calib.txt"
    n_ion, meta = from_conditions(lev, tr, st["T_gas"], st["t_exp"], st["rho"],
                                  X_LAN, ION_FRAC, A_OF[ion])
    atom = ForestAtom.from_gsi(lev, tr, st["T_gas"], n_ion, st["t_exp"], tau_min=1e-3)
    if atom.n_opacity < 10:
        return {"t_d": t_d, "skipped": f"{atom.n_opacity} opacity lines"}
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005

    ref, ei, eo = [], [], []
    for s in SEEDS:
        r = run_mc(atom, st["r_core"], st["r_out"], st["t_exp"], lo, hi, n,
                   "sobolev_branch", seed=s, t_core=T_CORE, collect_events=True)
        ref.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
        e = r["events"]; ei.append(e[0]); eo.append(e[1])
    refm = {b: float(np.mean([x[b] for x in ref])) for b in BANDS}
    nu_in, nu_out = np.concatenate(ei), np.concatenate(eo)
    kern = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size),
                                                  ng, nu_lo=lo, nu_hi=hi)
    live = [b for b in BANDS if refm[b] > 1e-3]

    legs = {}
    for tag, mode in (("A_redist", "sobolev_group"), ("B_opacity", "expansion_branch"),
                      ("C_both", "expansion_group"), ("C_binned", "binned_group")):
        rows = []
        for s in SEEDS:
            kw = {"kernel": kern} if mode.endswith("_group") else {}
            r = run_mc(atom, st["r_core"], st["r_out"], st["t_exp"], lo, hi, n,
                       mode, seed=s, t_core=T_CORE, **kw)
            rows.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0]
                         for b, w in BANDS.items()})
        dF = {b: float(np.mean([x[b] for x in rows]) / refm[b] - 1) for b in live}
        legs[tag] = {"dF": dF, "band3800": dF.get("band3800"),
                     "worst": max(abs(v) for v in dF.values())}
    nb = (C / (BAND[1] * 1e-8), C / (BAND[0] * 1e-8))
    out = {"t_d": t_d, "n_ion": float(n_ion), "T_gas": st["T_gas"],
           "rho": st["rho"], "n_opacity": int(atom.n_opacity),
           "tau_max": float(atom.op_tau.max()), "ref_bands": refm, "legs": legs}
    out.update({f"band_{k}": v for k, v in band_saturation(atom, *nb).items()})
    return out


def main(ion, n):
    print(f"{ion}: rho(1d) = {RHO_1D:.1e} g/cm3, X_lan = {X_LAN}, "
          f"T(1d) = {T_GAS_1D:.0f} K, homologous rho ~ t^-3, T ~ t^-1/2", flush=True)
    print(f"{'t/d':>5s}{'T_gas':>7s}{'n_ion':>11s}{'tau_max':>9s}{'S_band':>9s}"
          f"{'A':>9s}{'B':>9s}{'C':>9s}{'C_bin':>9s}", flush=True)
    rows = []
    for t_d in EPOCHS:
        r = epoch(ion, t_d, n)
        rows.append(r)
        if "skipped" in r:
            print(f"{t_d:5.2f}  skipped -- {r['skipped']}", flush=True); continue
        L = r["legs"]
        def f(t):
            v = L[t]["band3800"]
            return "    --  " if v is None else f"{100*v:+8.1f}%"
        print(f"{t_d:5.2f}{r['T_gas']:7.0f}{r['n_ion']:11.1f}{r['tau_max']:9.2f}"
              f"{r['band_S_band']:9.1f}{f('A_redist')}{f('B_opacity')}"
              f"{f('C_both')}{f('C_binned')}", flush=True)
    (HERE / f"trajectory_{ion}.json").write_text(json.dumps(
        {"ion": ion, "n": n, "rho_1d": RHO_1D, "X_lan": X_LAN,
         "T_gas_1d": T_GAS_1D, "rows": rows}, indent=1))
    print(f"wrote trajectory_{ion}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="58CeII")
    ap.add_argument("--n", type=float, default=4e5)
    a = ap.parse_args()
    main(a.ion, int(a.n))
