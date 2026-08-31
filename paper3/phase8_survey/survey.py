"""Paper III E3b: the collapse across many real ions.

F34 left the phase diagram partial -- band-local saturation clearly controls the
grouped-closure failure and the redistribution axes clearly do not, but only
three real atoms were on the curve and one sat 2.5x off it. Three points cannot
distinguish a law from a coincidence.

The GSI archive holds 27 ions and all of them are already on disk. This runs the
same measurement on as many as are affordable, so the collapse is tested against
a real population rather than against synthetic forests alone.

Per ion, exactly the recipe every previous ion used (setup.py's, verbatim):
n_ion such that the strongest classical Sobolev depth in 3850-3950 A at 3000 K
equals 5. Then the reference (`sobolev_branch`), a kernel trained on its own
events, and the grouped closure -- the same three legs as §4.30, so the numbers
are directly comparable to La/Ce/Nd.

Ions with no line in the normalization window are skipped and reported: the
recipe has nothing to anchor on, which is a property of the ion, not a failure.

Usage: python survey.py [--ions all|<name> ...] [--n 500000]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase0_reference"):
    sys.path.insert(0, str(p))
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of
from redistribution import RedistributionKernel
from reference import BANDS, SEEDS

from sobolev.atomic_data import load_gsi
from sobolev.constants import C
from sobolev.forest_stats import band_saturation, redistribution_range
from sobolev.optical_depth import SIGMA_CLASSICAL
from sobolev.populations import boltzmann_fractions_from_levels, statistical_weight

DATA = ROOT / "data"
WINDOW = (3850.0, 3950.0)     # the normalization window, as for La/Ce/Nd
BAND = (3800.0, 3955.0)       # the band every real-atom failure lives in
TAU_MAX_TARGET = 5.0
NG = 32


def n_ion_for(lev_path, tr_path):
    """n_ion such that the ion's STRONGEST line has tau = 5.

    setup.py's original recipe pins tau_max = 5 inside 3850-3950 A, which works
    for La/Ce/Nd because their strongest lines sit in or near that window. It
    fails for a survey: an ion whose window holds only weak lines needs an absurd
    density to reach tau = 5 there, and every other line then goes
    astronomically thick. Yb II wants n_ion = 1.7e12 cm^-3, giving tau_max = 1.7e8
    and beta_min = 1.5e-8 -- a packet entering that resonance never escapes, and
    the branch chain does not terminate.

    Normalizing on the ion's global maximum instead is scale-free, applies to
    every ion identically, and bounds beta_min at (1-e^-5)/5 = 0.199 by
    construction. The survey therefore re-runs La/Ce/Nd under this recipe too, so
    every point on the collapse is normalized the same way; those numbers differ
    slightly from §4.30's, which used the window recipe.

    Returns (n_ion, lines in the 3850-3950 A window, for reference).
    """
    levels = load_gsi(lev_path); lines = load_gsi(tr_path)
    lam_all = lines["WV_Transition"].to_numpy()
    frac = boltzmann_fractions_from_levels(levels, T_SHELL)
    pop = frac[lines["Lower"].to_numpy()]
    g_l = statistical_weight(lines["J_Lower"].to_numpy())
    f_lu = 10 ** lines["Log(gf)"].to_numpy() / g_l
    tau_per_n = SIGMA_CLASSICAL * f_lu * pop * lam_all * 1e-8 * T_EXP
    mx = np.nanmax(tau_per_n)
    n_win = int(((lam_all >= WINDOW[0]) & (lam_all < WINDOW[1])).sum())
    if not np.isfinite(mx) or mx <= 0:
        return None, n_win
    return TAU_MAX_TARGET / mx, n_win


def measure(ion, n):
    lev = DATA / f"{ion}_levels_calib.txt"
    tr = DATA / f"{ion}_transitions_calib.txt"
    if not (lev.exists() and tr.exists()):
        return {"ion": ion, "skipped": "files not extracted"}
    t0 = time.time()
    n_ion, n_win = n_ion_for(lev, tr)
    if n_ion is None:
        return {"ion": ion, "skipped": f"no usable line in {WINDOW[0]:.0f}-{WINDOW[1]:.0f} A",
                "n_window": n_win}
    atom = ForestAtom.from_gsi(lev, tr, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    if atom.n_opacity < 10:
        return {"ion": ion, "skipped": f"only {atom.n_opacity} opacity lines"}
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    build = time.time() - t0

    ref, ev_in, ev_out = [], [], []
    for s in SEEDS:
        r = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_branch",
                   seed=s, t_core=T_CORE, collect_events=True)
        ref.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
        e = r["events"]; ev_in.append(e[0]); ev_out.append(e[1])
    refm = {b: float(np.mean([x[b] for x in ref])) for b in BANDS}
    nu_in, nu_out = np.concatenate(ev_in), np.concatenate(ev_out)
    kern = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size),
                                                  NG, nu_lo=lo, nu_hi=hi)

    legs = {}
    for mode in ("sobolev_group", "binned_group", "expansion_group", "expansion_branch"):
        rows = []
        for s in SEEDS:
            kw = {"kernel": kern} if mode.endswith("_group") else {}
            r = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode,
                       seed=s, t_core=T_CORE, **kw)
            rows.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0]
                         for b, w in BANDS.items()})
        live = [b for b in BANDS if refm[b] > 1e-3]
        dF = {b: float(np.mean([x[b] for x in rows]) / refm[b] - 1) for b in live}
        legs[mode] = {"dF": dF, "worst": max(abs(v) for v in dF.values()),
                      "band3800": dF.get("band3800")}

    nu_b = (C / (BAND[1] * 1e-8), C / (BAND[0] * 1e-8))
    out = {"ion": ion, "n_ion": float(n_ion), "n_window": n_win,
           "n_lines_total": int(atom.n_lines_total), "n_opacity": int(atom.n_opacity),
           "tau_max": float(atom.op_tau.max()), "ref_bands": refm,
           "n_events": int(nu_in.size), "build_s": build,
           "wall_s": time.time() - t0, "legs": legs}
    out.update({f"band_{k}": v for k, v in
                band_saturation(atom, nu_b[0], nu_b[1]).items()})
    out.update(redistribution_range(nu_in, nu_out, edges=kern.edges))
    return out


ALL = ["57LaII", "58CeII", "60NdII", "59PrII", "59PrIII", "60NdIII",
       "66DyIII", "67HoIII", "68ErIII", "69TmII", "69TmIII", "70YbII",
       "58CeIII"]


def main(ions, n):
    rows = []
    for ion in ions:
        r = measure(ion, n)
        rows.append(r)
        if "skipped" in r:
            print(f"  {ion:9s} SKIPPED -- {r['skipped']}", flush=True)
            continue
        print(f"  {ion:9s} op {r['n_opacity']:6d}  band {r['band_n_band']:5d} "
              f"Nsat {r['band_n_sat_band']:4d}  S {r['band_S_band']:8.1f} | "
              f"Rij {100*r['legs']['sobolev_group']['worst']:6.2f}%  "
              f"binned {100*r['legs']['binned_group']['worst']:8.2f}%  "
              f"exp {100*r['legs']['expansion_group']['worst']:8.2f}%  "
              f"[{r['wall_s']:.0f}s]", flush=True)
    (HERE / "survey.json").write_text(json.dumps({"n": n, "ng": NG, "rows": rows}, indent=1))
    print("wrote survey.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ions", nargs="*", default=["all"])
    ap.add_argument("--n", type=float, default=5e5)
    a = ap.parse_args()
    main(ALL if a.ions == ["all"] else a.ions, int(a.n))
