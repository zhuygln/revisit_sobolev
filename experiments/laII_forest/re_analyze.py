"""Reduce the radiative-equilibrium pair (re_run.py) to what the paper quotes.

Per variant (N = 1: redistribution at the input T; N = 15: self-consistent T)
and per seed: the emergent band ratio in 3800-3955 A (blue-margin
normalization), the paired resolved-vs-expansion differential, the fill-in
fraction F^RE / F^noRE, the band just redward (3955-3995 A, where
re-emission from the forest lands inside the transport window), the
temperature profile at the last iteration and its convergence over the last
three, and the null (tau_max = 0.05, rho x 0.01: both modes at the continuum).

The RE-off production pair gives the F^noRE denominators; their blue-margin
values are the control on the margin choice.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.spectra import band_ratio

R_CORE, T_CORE = 8.64e12, 6000.0
BAND, RED_BAND = (3800.0, 3955.0), (3955.0, 3978.0)  # red band stops well before the grid edge at 3997 A
BLUE_MARGIN = (3785.0, 3805.0)
res = json.loads((HERE / "re_results.json").read_text())


def t_profile(run, it):
    p = run / f"plt_{it:05d}.dat"
    if not p.exists():
        return None
    a = np.loadtxt(p, comments="#")
    return a[:, 0], a[:, 3]  # r, T_gas


def last_ratio(run, n_iter, band):
    spec = run / f"spectrum_{n_iter}.dat"
    return band_ratio(spec, band, BLUE_MARGIN, R_CORE, T_CORE) if spec.exists() else None


# RE-off reference (production), blue-margin normalized for like-for-like
noRE = {m: band_ratio(HERE / f"run_{m}" / "spectrum_1.dat", BAND, BLUE_MARGIN, R_CORE, T_CORE) for m in ("bb", "exp")}
noRE_red = {m: band_ratio(HERE / f"run_{m}" / "spectrum_1.dat", RED_BAND, BLUE_MARGIN, R_CORE, T_CORE) for m in ("bb", "exp")}
print(f"RE off (blue margin): bb {noRE['bb']:.4f}  exp {noRE['exp']:.4f}  -> Delta_exp {100*(noRE['exp']/noRE['bb']-1):+.1f}%   red band bb {noRE_red['bb']:.4f} exp {noRE_red['exp']:.4f}")

summary = {"noRE": dict(noRE), "noRE_red": dict(noRE_red), "variants": {}}
for tag, n_iter in (("prod", 1), ("prod", 15), ("null", 15)):
    pairs = []
    for seed in (1, 2, 3):
        rb = HERE / f"re_{tag}_N{n_iter}_bb_s{seed}"; re_ = HERE / f"re_{tag}_N{n_iter}_exp_s{seed}"
        fb, fe = last_ratio(rb, n_iter, BAND), last_ratio(re_, n_iter, BAND)
        if fb is None or fe is None:
            continue
        fbr, fer = last_ratio(rb, n_iter, RED_BAND), last_ratio(re_, n_iter, RED_BAND)
        entry = dict(seed=seed, f_bb=fb, f_exp=fe, d_exp=(fe - fb) / fb, f_bb_red=fbr, f_exp_red=fer)
        # T convergence from plt files
        for mode, run in (("bb", rb), ("exp", re_)):
            prof = [t_profile(run, it) for it in range(max(1, n_iter - 2), n_iter + 1)]
            prof = [x for x in prof if x is not None]
            if prof:
                T_last = prof[-1][1]
                entry[f"T_{mode}_inner"] = float(T_last[0]); entry[f"T_{mode}_outer"] = float(T_last[-1]); entry[f"T_{mode}_median"] = float(np.median(T_last))
                if len(prof) >= 2:
                    dT = max(float(np.max(np.abs(b[1] - a[1]) / a[1])) for a, b in zip(prof[:-1], prof[1:]))
                    entry[f"T_{mode}_max_frac_change_last3"] = dT
        # per-iteration band ratio trace (from re_results.json)
        key_b, key_e = f"{tag}_N{n_iter}_bb_s{seed}", f"{tag}_N{n_iter}_exp_s{seed}"
        entry["trace_bb"], entry["trace_exp"] = res["runs"].get(key_b), res["runs"].get(key_e)
        pairs.append(entry)
    if not pairs:
        continue
    d = np.array([p["d_exp"] for p in pairs]); fb = np.array([p["f_bb"] for p in pairs]); fe = np.array([p["f_exp"] for p in pairs])
    v = dict(n_pairs=len(pairs), f_bb_mean=float(fb.mean()), f_exp_mean=float(fe.mean()),
             d_exp_mean=float(d.mean()), d_exp_std=float(d.std(ddof=1)) if len(d) > 1 else None,
             fill_in_bb=float(fb.mean() / noRE["bb"]), fill_in_exp=float(fe.mean() / noRE["exp"]),
             pairs=pairs)
    summary["variants"][f"{tag}_N{n_iter}"] = v
    print(f"\n{tag} N={n_iter}: {len(pairs)} pairs")
    print(f"  F_bb {fb.mean():.4f}  F_exp {fe.mean():.4f}  ->  Delta_exp(RE) {100*d.mean():+.2f}%" + (f" +- {100*d.std(ddof=1):.2f}" if len(d) > 1 else ""))
    print(f"  fill-in: bb {fb.mean()/noRE['bb']:.2f}x  exp {fe.mean()/noRE['exp']:.2f}x   (RE off: {noRE['bb']:.4f}, {noRE['exp']:.4f})")
    p0 = pairs[0]
    if "T_bb_median" in p0:
        print(f"  T_gas (seed 1, last iter): bb inner {p0['T_bb_inner']:.0f} K median {p0['T_bb_median']:.0f} K | exp inner {p0.get('T_exp_inner', float('nan')):.0f} K median {p0.get('T_exp_median', float('nan')):.0f} K"
              + (f" | max |dT/T| over last 3 iters: bb {p0.get('T_bb_max_frac_change_last3', float('nan')):.3f} exp {p0.get('T_exp_max_frac_change_last3', float('nan')):.3f}" if n_iter > 1 else ""))
    if p0.get("f_bb_red") is not None:
        print(f"  red band 3955-3995: bb {np.mean([p['f_bb_red'] for p in pairs]):.4f}  exp {np.mean([p['f_exp_red'] for p in pairs]):.4f}  (RE off: {noRE_red['bb']:.4f}, {noRE_red['exp']:.4f})")
    if p0.get("trace_bb"):
        print(f"  band ratio by iteration (seed 1, bb): " + " ".join(f"{x:.3f}" for x in p0["trace_bb"][:15]))
(HERE / "re_summary.json").write_text(json.dumps(summary, indent=1))
print("\nwrote re_summary.json")
