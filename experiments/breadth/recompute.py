"""Recompute breadth band fluxes with a red margin clear of the grid edge.

sweep.py normalized against 3800-3955-equivalent red margins running to
hi*1.008, while the transport grid ended at hi*1.010. SEDONA's last spectrum
bin is a partial bin and collapses (L/L_cont 48.5 -> 15.4 in the line-free
control), so the margin straddled it, depressed the reference, and inflated
every band flux by ~7%.

Effect by quantity:
  Delta_expansion -- SAFE. Both fluxes carry the same bias and it cancels;
                     the line-free control windows correctly gave ~0%.
  Delta_Sobolev   -- CONTAMINATED. The analytic leg is correctly normalized
                     to 1, so it was compared against an inflated SEDONA
                     value and picked up a spurious ~-7%.

This is the F8 lesson in a new guise: same-code differentials are robust,
cross-code comparisons are only as good as the shared normalization.

The spectra are on disk, so no SEDONA re-runs are needed.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[1]))
from sobolev.constants import C

rows = json.loads((HERE / "breadth_results.json").read_text())
out = []
for r in rows:
    lo, hi = r["window"], r["window"] + 100.0
    tag = f"{r['mix']}_w{int(lo)}_t{r['t_day']:g}"
    band_lo, band_hi = lo * 0.9885, hi + 1.0
    # keep the margin inside the grid: grid ends at hi*1.010, the final bin
    # is unreliable, so stop at hi*1.006.
    red_lo, red_hi = hi * 1.002, hi * 1.006

    fluxes = {}
    for mode in ("bb", "exp"):
        p = HERE / f"run_{tag}_{mode}" / "spectrum_1.dat"
        if not p.exists():
            fluxes[mode] = None
            continue
        s = np.loadtxt(p, comments="#")
        nu, lum = s[:, 0], s[:, 1]
        lam = C / nu * 1e8
        red = (lam > red_lo) & (lam < red_hi) & (lum > 0)
        if red.sum() < 3:
            fluxes[mode] = None
            continue
        ratio = lum / np.mean(lum[red])
        m = (lam > band_lo) & (lam < band_hi)
        o = np.argsort(lam[m])
        fluxes[mode] = float(
            np.trapezoid(ratio[m][o], lam[m][o]) / (band_hi - band_lo)
        )

    n = dict(r)
    n["f_res"], n["f_exp"] = fluxes["bb"], fluxes["exp"]
    if n["f_res"] and n["f_exp"]:
        n["d_exp"] = (n["f_exp"] - n["f_res"]) / n["f_res"]
        n["d_sob"] = (n["f_sob"] - n["f_res"]) / n["f_res"]
    out.append(n)
    print(
        f"{tag:22s} tau_max={r['tau_max']:8.2f}  res={n['f_res']:.4f} "
        f"(was {r['f_res']:.4f})  D_sob={n['d_sob']:+7.1%} "
        f"(was {r['d_sob']:+.1%})  D_exp={n['d_exp']:+7.1%}"
    )

(HERE / "breadth_results_fixed.json").write_text(json.dumps(out, indent=1))

ctrl = [r for r in out if r["tau_max"] < 0.05]
if ctrl:
    v = np.array([r["f_res"] for r in ctrl])
    print(
        f"\nline-free controls (tau_max < 0.05): f_res = {v.mean():.4f} "
        f"+- {v.std():.4f}  (must be 1.000)"
    )
