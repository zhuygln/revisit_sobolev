"""Paper III audit: which cross-ion claims survive the normalization change?

F35 showed the project's `window_tau_max` recipe is ion-specific by accident.
Every claim that compares ions -- F24's density limit, F27's Gate 2, F30's
opacity decomposition, F31/F33's memory result -- was measured under it. This
re-measures the same quantities under the CONTROLLED standard (`global_tau_max`,
every ion at the same strongest-line depth) so each claim can be classified
invariant / conditional / superseded.

Measured per ion, all legs against that ion's own `sobolev_branch`:

  compression   N_g = 4, 8, 32          -> F25, F27
  opacity       binned / expansion      -> F30
  memory        m = 0, 1, 4 on binned   -> F31, F33

Usage: python audit.py [--ions ...] [--n 500000]
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

from sobolev.constants import C
from sobolev.forest_stats import band_saturation, redistribution_range
from sobolev.normalization import global_tau_max

DATA = ROOT / "data"
BAND = (3800.0, 3955.0)
IONS = ["57LaII", "58CeII", "60NdII", "59PrII", "70YbII"]


def legs(atom, lo, hi, n, kern, modes):
    out = {}
    for tag in modes:
        mode, kw = tag, {}
        if "+m" in tag:
            mode, m = tag.split("+m"); kw["line_memory"] = int(m)
        if mode.endswith("_group"):
            kw["kernel"] = kern
        rows = []
        for s in SEEDS:
            r = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode,
                       seed=s, t_core=T_CORE, **kw)
            rows.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0]
                         for b, w in BANDS.items()})
        out[tag] = {b: float(np.mean([x[b] for x in rows])) for b in BANDS}
    return out


def audit_ion(ion, n):
    lev, tr = DATA / f"{ion}_levels_calib.txt", DATA / f"{ion}_transitions_calib.txt"
    if not (lev.exists() and tr.exists()):
        return {"ion": ion, "skipped": "not extracted"}
    n_ion, meta = global_tau_max(lev, tr, T_SHELL, T_EXP)
    if n_ion is None:
        return {"ion": ion, "skipped": meta["reason"]}
    t0 = time.time()
    atom = ForestAtom.from_gsi(lev, tr, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005

    ref, ei, eo = [], [], []
    for s in SEEDS:
        r = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_branch",
                   seed=s, t_core=T_CORE, collect_events=True)
        ref.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
        e = r["events"]; ei.append(e[0]); eo.append(e[1])
    refm = {b: float(np.mean([x[b] for x in ref])) for b in BANDS}
    nu_in, nu_out = np.concatenate(ei), np.concatenate(eo)
    live = [b for b in BANDS if refm[b] > 1e-3]

    def score(bands):
        dF = {b: bands[b] / refm[b] - 1 for b in live}
        return {"dF": dF, "worst": max(abs(v) for v in dF.values()),
                "band3800": dF.get("band3800")}

    out = {"ion": ion, "n_ion": float(n_ion), "normalization": meta,
           "n_opacity": int(atom.n_opacity), "tau_max": float(atom.op_tau.max()),
           "ref_bands": refm, "n_events": int(nu_in.size)}
    nb = (C / (BAND[1] * 1e-8), C / (BAND[0] * 1e-8))
    out.update({f"band_{k}": v for k, v in band_saturation(atom, *nb).items()})

    # F25/F27: does a small matrix still reproduce branching?
    out["compression"] = {}
    for ng in (4, 8, 32):
        k = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size),
                                                   ng, nu_lo=lo, nu_hi=hi)
        out["compression"][str(ng)] = score(legs(atom, lo, hi, n, k, ["sobolev_group"])["sobolev_group"])
        if ng == 32:
            kern32 = k
    out.update(redistribution_range(nu_in, nu_out, edges=kern32.edges))

    # F30 opacity decomposition, F31/F33 memory -- all on the N_g = 32 kernel
    tags = ["binned_group", "expansion_group", "expansion_branch",
            "binned_group+m1", "binned_group+m4"]
    got = legs(atom, lo, hi, n, kern32, tags)
    out["legs"] = {t: score(got[t]) for t in tags}
    out["wall_s"] = time.time() - t0
    return out


def main(ions, n):
    rows = []
    for ion in ions:
        r = audit_ion(ion, n)
        rows.append(r)
        if "skipped" in r:
            print(f"  {ion:9s} SKIPPED -- {r['skipped']}", flush=True); continue
        c = r["compression"]; L = r["legs"]
        def pc(t):   # band3800 is None when the reference blacks that band out
            v = L[t]["band3800"]
            return "    --  " if v is None else f"{100*v:+7.1f}%"
        print(f"  {ion:9s} op {r['n_opacity']:6d} tau_max {r['tau_max']:5.2f} "
              f"Nsat {r['band_n_sat_band']:3d} S {r['band_S_band']:7.1f} | "
              f"Ng4 {100*c['4']['worst']:6.2f}% Ng32 {100*c['32']['worst']:5.2f}% | "
              f"bin {pc('binned_group')} exp {pc('expansion_group')} "
              f"m1 {pc('binned_group+m1')} m4 {pc('binned_group+m4')}  "
              f"[{r['wall_s']:.0f}s]", flush=True)
    (HERE / "audit.json").write_text(json.dumps({"n": n, "rows": rows}, indent=1))
    print("wrote audit.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ions", nargs="*", default=IONS)
    ap.add_argument("--n", type=float, default=5e5)
    a = ap.parse_args()
    main(a.ions, int(a.n))
