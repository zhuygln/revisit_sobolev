"""Breadth sweep through the canonical normalization, with a deterministic
reference (referee Comments 3, 5, 8; repo finding: the stale Delta_Sob row).

recompute.py normalized SEDONA band fluxes by RAW luminosity in the red
margin (`ratio = lum / np.mean(lum[red])`), leaving the Planck slope across
the band in the answer -- the bug commit 1e2ba21 fixed everywhere else via
sobolev.spectra.band_ratio, and which had produced a spurious 5-8%
"v_D-independent Sobolev floor" of the same sign and size as the breadth
median the referee flagged. The breadth directory never migrated.

What moves and what must not:
  Delta_exp  -- SEDONA vs SEDONA, same-code: must change by < MC noise.
  Delta_sob  -- analytic vs SEDONA, cross-code: expected to fall.
  New: Delta_sob_det, analytic Sobolev vs the closed-form resolved leg on
  identical rays -- no Monte Carlo, no normalization, no emission convention.
  New: the stimulated-emission factor (SEDONA applies it; 5e-3 at 9100 A),
  reported with and without.
  New: the transmission-weighted saturation deficit per condition, the
  ensemble statistic that controls Delta_exp (Comment 8).

No SEDONA reruns: all 72 spectra are on disk.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C
from sobolev.rays import RaySet
from sobolev.sobolev_leg import crossing_depths, expansion_damp, resolved_attenuation, sobolev_attenuation
from sobolev.spectra import band_average, band_ratio
from lines import T_CORE, V_D, all_conditions

old = {f"{r['mix']}_w{int(r['window'])}_t{r['t_day']:g}": r
       for r in json.loads((HERE / "breadth_results_fixed.json").read_text())}

out = []
print(f"{'tag':22s} {'tau_max':>8s} {'f_res':>7s} {'D_exp':>7s} {'(old)':>7s} {'D_sob':>7s} {'(old)':>7s} {'D_sob_det':>9s} {'D_exp_det':>9s} {'lnEw[e^D]':>9s}")
for stim in (False, True):
    for c in all_conditions(stim=stim):
        tag = c["tag"]
        r_core, r_out, t_exp, rho = c["r_core"], c["r_out"], c["t_exp"], c["rho"]
        band, margin = c["band"], c["red_margin"]
        nu = np.geomspace(c["nu_lo"], c["nu_hi"], 2000)
        lam = C / nu * 1e8
        rays = RaySet.midpoint(r_core, r_out, 200)
        kw = dict(r_core=r_core, r_out=r_out, t_exp=t_exp, n_ref=rho, rays=rays)

        f_sob = band_average(lam, sobolev_attenuation(nu, c["lines"], relativity="first", **kw), band)
        f_exa = band_average(lam, sobolev_attenuation(nu, c["lines"], damp=expansion_damp, relativity="first", **kw), band)
        f_det = band_average(lam, resolved_attenuation(nu, c["lines"], v_doppler=V_D, sweep="first", **kw), band)
        # the exact ensemble statistic: F_exp/F_sob = E_w[e^{S-E}], w ~ p e^{-S}
        S, E, p, w = crossing_depths(nu, c["lines"], r_core, r_out, t_exp, rho, rays=rays)
        m = (lam > band[0]) & (lam < band[1])
        wt = (w[:, None] * np.exp(-S))[:, m]
        ln_ew = float(np.log(np.sum(wt * np.exp((S - E)[:, m])) / wt.sum()))
        d_w = float(np.sum(wt * (S - E)[:, m]) / wt.sum())          # <D>_w, the Jensen lower bound
        d_naive = float(np.mean((S - E)[:, m]))                       # the referee's unweighted deficit

        row = dict(old[tag]); row.update(stim=stim, f_sob=f_sob, f_exp_ana=f_exa, f_res_det=f_det,
                                        d_sob_det=(f_sob - f_det) / f_det, d_exp_det=(f_exa - f_det) / f_det,
                                        ln_ew_expD=ln_ew, d_w=d_w, d_naive=d_naive)
        if not stim:  # SEDONA fluxes are stim-independent (SEDONA always applies it); compute once
            fl = {}
            for mode in ("bb", "exp"):
                spec = HERE / f"run_{tag}_{mode}" / "spectrum_1.dat"
                fl[mode] = band_ratio(spec, band, margin, r_core, T_CORE) if spec.exists() else None
            row.update(f_res=fl["bb"], f_exp=fl["exp"])
        else:
            prev = next(x for x in out if x["tag"] == tag and not x["stim"])
            row.update(f_res=prev["f_res"], f_exp=prev["f_exp"])
        row["tag"] = tag
        if row["f_res"] and row["f_exp"]:
            row["d_exp"] = (row["f_exp"] - row["f_res"]) / row["f_res"]
            row["d_sob"] = (f_sob - row["f_res"]) / row["f_res"]
            row["d_det_vs_sedona"] = (f_det - row["f_res"]) / row["f_res"]
        out.append(row)
        if stim:
            o = old[tag]
            print(f"{tag:22s} {row['tau_max']:8.2f} {row['f_res']:7.4f} {100*row['d_exp']:+6.1f}% {100*o['d_exp']:+6.1f}% "
                  f"{100*row['d_sob']:+6.1f}% {100*o['d_sob']:+6.1f}% {100*row['d_sob_det']:+8.2f}% {100*row['d_exp_det']:+8.1f}% {ln_ew:9.4f}")

(HERE / "breadth_results_v2.json").write_text(json.dumps(out, indent=1))

def med(key, rows): return float(np.median([r[key] for r in rows]))
S_ = [r for r in out if r["stim"]]; N_ = [r for r in out if not r["stim"]]
print("\n=== summary (stim on) ===")
print(f"  D_exp     median {100*med('d_exp',S_):+.2f}%  max {100*max(r['d_exp'] for r in S_):+.1f}%   (old file: {100*med('d_exp',list(old.values())):+.2f}%, {100*max(r['d_exp'] for r in old.values()):+.1f}%)")
print(f"  D_sob     median {100*med('d_sob',S_):+.2f}%  max {100*max(r['d_sob'] for r in S_):+.1f}%   (old file: {100*med('d_sob',list(old.values())):+.2f}%, {100*max(r['d_sob'] for r in old.values()):+.1f}%)")
print(f"  D_sob_det median {100*med('d_sob_det',S_):+.2f}%  max {100*max(r['d_sob_det'] for r in S_):+.2f}%   <- deterministic reference, no MC")
print(f"  D_exp_det median {100*med('d_exp_det',S_):+.2f}%  max {100*max(r['d_exp_det'] for r in S_):+.1f}%")
hi = [r for r in S_ if r['tau_max'] > 3]; lo = [r for r in S_ if r['tau_max'] < 0.5]
print(f"  tau_max>3:  D_exp {100*med('d_exp',hi):+.1f}%  D_sob {100*med('d_sob',hi):+.2f}%  D_sob_det {100*med('d_sob_det',hi):+.2f}%")
print(f"  tau_max<0.5: D_exp {100*med('d_exp',lo):+.1f}%  D_sob {100*med('d_sob',lo):+.2f}%  D_sob_det {100*med('d_sob_det',lo):+.2f}%")
print(f"  det vs SEDONA resolved: median {100*med('d_det_vs_sedona',S_):+.2f}%, spread {100*np.std([r['d_det_vs_sedona'] for r in S_]):.2f}")
dd = [abs(r['d_exp'] - old[r['tag']]['d_exp']) for r in N_]
print(f"\n  NULL  D_exp change vs old file: max {100*max(dd):.3f} points (must be < MC noise ~0.4)")
ctrl = [r for r in S_ if r['tau_max'] < 0.05]
if ctrl:
    print(f"  NULL  line-free windows: f_res {np.mean([r['f_res'] for r in ctrl]):.4f}+-{np.std([r['f_res'] for r in ctrl]):.4f}, f_sob {np.mean([r['f_sob'] for r in ctrl]):.5f}, f_det {np.mean([r['f_res_det'] for r in ctrl]):.5f}")
st = [(r['tag'], r['f_sob']) for r in S_]; ns = {r['tag']: r['f_sob'] for r in N_}
shift = max(abs(fs - ns[t]) / ns[t] for t, fs in st)
print(f"  stim factor: max |shift| in f_sob {100*shift:.3f}% (9100 A window); at 4300 A {100*max(abs(r['f_sob']-ns[r['tag']])/ns[r['tag']] for r in S_ if r['window']==4300.0):.5f}%")
