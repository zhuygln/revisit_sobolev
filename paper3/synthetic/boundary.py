"""Paper III E4: locate the Delta F = 0 contour, and test whether it is 2-D.

F35 found the grouped-closure error changes SIGN -- too opaque at low band
saturation, too bright at high -- and located the crossing at S ~ 50 three ways.
F34 found the redistribution axes are null for the error's MAGNITUDE. Those are
different questions, and the one that matters now is:

    does the sign boundary MOVE when redistribution range changes?

If it does not, the boundary is one-dimensional in S and the criterion is a
number. If it does, the phase diagram is genuinely 2-D and the second axis is
required -- which would also explain La II vs Pr II, identical in S and 5x apart
in error, without invoking a sign flip they do not have (§4.32 correction).

Design: rather than a blind grid, scan tau at fixed (n_lines, dlnlam) to
BRACKET the crossing, then interpolate S at Delta F = 0 per row. That spends
the compute where the contour is instead of where it is not.

Usage: python boundary.py [--n 150000]
"""
import argparse, itertools, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", HERE):
    sys.path.insert(0, str(p))
import phase as P
from forest import synthetic_forest
from forest_mc import run_mc
from redistribution import RedistributionKernel

from sobolev.forest_stats import band_saturation, redistribution_range

SEEDS = P.SEEDS


def point(n_lines, tau, dlnlam, n, ng=32, seed=7, delocalize=0.0, n_exit=2,
          exit_tau=0.0):
    """One forest: signed core-band error for the binned and expansion closures.

    `delocalize` is the F37 fix: without it the model has no net inflow to the
    measured band and only the too-opaque branch exists.
    """
    atom, _ = synthetic_forest(n_lines=n_lines, tau=tau, tau_spread=1.8,
                               span=0.3, n_exit=n_exit, dlnlam=dlnlam,
                               f_return=0.5, jitter=0.5, seed=seed,
                               delocalize=delocalize, exit_tau=exit_tau)
    lo, hi = atom.op_nu.min() * 0.99, atom.op_nu.max() * 1.01
    ref, ev, _ = P.measure(atom, lo, hi, n, "sobolev_branch")
    if ev is None or ev[0].size < 1000 or ref["core"] <= 1e-3:
        return None
    kern = RedistributionKernel.from_branching_mc(ev[0], ev[1], np.ones(ev[0].size),
                                                  ng, nu_lo=lo, nu_hi=hi)
    core = P.bands_of(lo, hi)["core"]
    bs = band_saturation(atom, core[0], core[1])
    rr = redistribution_range(ev[0], ev[1], edges=kern.edges)
    out = {"n_lines": n_lines, "tau": tau, "dlnlam": dlnlam,
           "delocalize": delocalize, "n_exit": n_exit, "exit_tau": exit_tau,
           "S_band": bs["S_band"], "n_sat_band": bs["n_sat_band"],
           "n_band": bs["n_band"], "range": rr["mean_abs_dlnlam"],
           "same_group_frac": rr["same_group_frac"], "ref_core": ref["core"]}
    for mode, tag in (("binned_group", "binned"), ("expansion_group", "expansion")):
        b, _, ipp = P.measure(atom, lo, hi, n, mode, kernel=kern)
        out[tag] = b["core"] / ref["core"] - 1
        out[tag + "_ipp"] = ipp
    return out


def crossing(rows, key):
    """S at which the signed error crosses zero, by log-linear interpolation.

    Detects the crossing in EITHER direction and reports which. The first
    version tested only `a < 0 <= b` and so silently missed positive-to-negative
    crossings -- which is exactly what the delocalized model turned out to
    produce, the opposite direction from the real ions.
    """
    r = sorted([x for x in rows if x[key] is not None], key=lambda x: x["S_band"])
    for a, b in zip(r, r[1:]):
        if a["S_band"] <= 0 or b["S_band"] <= 0:
            continue
        if a[key] * b[key] < 0:
            f = -a[key] / (b[key] - a[key])
            S = float(np.exp(np.log(a["S_band"]) +
                             f * (np.log(b["S_band"]) - np.log(a["S_band"]))))
            return {"S": S, "direction": "neg_to_pos" if a[key] < 0 else "pos_to_neg"}
    return None


def main(n, out_name):
    N_LINES = (100, 300)
    DLNLAM = (0.005, 0.03, 0.15)
    TAU = (0.02, 0.06, 0.15, 0.4, 1.0, 2.5)
    combos = list(itertools.product(N_LINES, DLNLAM))
    print(f"{len(combos)} rows x {len(TAU)} tau, {n} packets x {len(SEEDS)} seeds", flush=True)
    all_rows, summary, t0 = [], [], time.time()
    for nl, dl in combos:
        rows = []
        for tau in TAU:
            r = point(nl, tau, dl, n)
            if r is None:
                continue
            rows.append(r); all_rows.append(r)
            print(f"  N={nl:3d} dln={dl:.3f} tau={tau:5.2f} | S {r['S_band']:8.1f} "
                  f"Nsat {r['n_sat_band']:4d} rng {r['range']:.3f} | "
                  f"binned {100*r['binned']:+8.1f}%  exp {100*r['expansion']:+8.1f}%",
                  flush=True)
        s = {"n_lines": nl, "dlnlam": dl,
             "S_cross_binned": crossing(rows, "binned"),
             "S_cross_expansion": crossing(rows, "expansion"),
             "range_median": float(np.median([r["range"] for r in rows])) if rows else None}
        summary.append(s)
        print(f"   -> zero crossing: binned S = {s['S_cross_binned']}, "
              f"expansion S = {s['S_cross_expansion']}", flush=True)
    (HERE / out_name).write_text(json.dumps(
        {"n": n, "rows": all_rows, "summary": summary}, indent=1))
    print(f"wrote {out_name}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=1.5e5)
    ap.add_argument("--out", default="boundary.json")
    a = ap.parse_args()
    main(int(a.n), a.out)
