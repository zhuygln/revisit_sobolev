"""Reduce the seed-matched pairs (seeds.py) to what the paper quotes.

Per case: seed mean and standard error of each leg, the PAIRED Delta_exp
(same seed in both legs) with its mean, std and sem, and the correlation
between the legs. The quadrature estimate is printed beside the paired one so
the gain from matching seeds is visible. Also the Poisson expectation for the
single-realization scatter, from the packet counts in the spectrum file
(column 3), as a null: the measured seed std should sit within ~1.3x of it.

Reads mc_noise_seeds.json (and folds in the five pairs of run.py's
mc_noise_results.json where they overlap); writes mc_noise_summary.json.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C

BAND = (3800.0, 3955.0)


def poisson_frac_err(spec_path):
    """sqrt(N)/N over the packets landing in the band, from column 3."""
    try:
        s = np.loadtxt(spec_path, comments="#")
    except Exception:
        return None
    if s.shape[1] < 3:
        return None
    lam = C / s[:, 0] * 1e8
    m = (lam > BAND[0]) & (lam < BAND[1])
    n = s[m, 2].sum()
    return float(1.0 / np.sqrt(n)) if n > 0 else None


seeds = json.loads((HERE / "mc_noise_seeds.json").read_text())
summary = {}
print(f"{'case':14s} {'n':>2s} {'F_bb':>8s} {'+-sem':>7s} {'F_exp':>8s} {'+-sem':>7s} | {'D_exp paired':>12s} {'std':>6s} {'sem':>6s} {'quad':>6s} {'corr':>6s} | {'poisson':>7s} {'seed-std':>8s}")
for case, legs in sorted(seeds.items()):
    bb = {int(k): v for k, v in legs.get("bb", {}).items() if v is not None}
    ex = {int(k): v for k, v in legs.get("exp", {}).items() if v is not None}
    # fold in the legacy five where this case is the headline
    if case == "tau5_vd100" or case == "tau5_vd10":
        legacy = json.loads((HERE / "mc_noise_results.json").read_text())
        for mode, d in (("bb", bb), ("exp", ex)):
            key = f"{case}_{mode}"
            if key in legacy:
                for i, v in enumerate(legacy[key]["values"], 1):
                    d.setdefault(i, v)
    common = sorted(set(bb) & set(ex))
    if len(bb) < 2:
        continue
    fb = np.array([bb[k] for k in sorted(bb)]); fe = np.array([ex[k] for k in sorted(ex)]) if ex else np.array([])
    out = dict(n_bb=len(fb), n_exp=len(fe),
               f_bb_mean=float(fb.mean()), f_bb_sem=float(fb.std(ddof=1) / np.sqrt(len(fb))),
               f_bb_std=float(fb.std(ddof=1)))
    if len(fe) >= 2:
        out.update(f_exp_mean=float(fe.mean()), f_exp_sem=float(fe.std(ddof=1) / np.sqrt(len(fe))),
                   f_exp_std=float(fe.std(ddof=1)))
    if len(common) >= 2:
        pb = np.array([bb[k] for k in common]); pe = np.array([ex[k] for k in common])
        d = (pe - pb) / pb
        out.update(n_pairs=len(common), d_exp_paired_mean=float(d.mean()),
                   d_exp_paired_std=float(d.std(ddof=1)), d_exp_paired_sem=float(d.std(ddof=1) / np.sqrt(len(d))),
                   corr=float(np.corrcoef(pb, pe)[0, 1]) if len(common) > 2 else None,
                   d_exp_quadrature_std=float(np.hypot(pb.std(ddof=1) / pb.mean(), pe.std(ddof=1) / pe.mean())))
    # poisson null from the first bb spectrum
    for k in sorted(bb):
        spec = HERE / f"run_{case}_bb_s{k}" / "spectrum_1.dat"
        if spec.exists():
            out["poisson_frac_err"] = poisson_frac_err(spec); break
    summary[case] = out
    pe_s = f"{out.get('poisson_frac_err', float('nan')):7.4f}" if out.get("poisson_frac_err") else "   n/a "
    corr_s = "  n/a " if out.get("corr") is None else f"{out['corr']:+.3f}"
    print(f"{case:14s} {len(fb):2d} {out['f_bb_mean']:8.5f} {out['f_bb_sem']:7.5f} "
          f"{out.get('f_exp_mean', float('nan')):8.5f} {out.get('f_exp_sem', float('nan')):7.5f} | "
          f"{100*out.get('d_exp_paired_mean', float('nan')):+11.2f}% {100*out.get('d_exp_paired_std', float('nan')):5.2f} "
          f"{100*out.get('d_exp_paired_sem', float('nan')):5.2f} {100*out.get('d_exp_quadrature_std', float('nan')):5.2f} "
          f"{corr_s:>6s} | {pe_s} {out['f_bb_std']/out['f_bb_mean']:8.4f}")
(HERE / "mc_noise_summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True))
print("wrote mc_noise_summary.json")
