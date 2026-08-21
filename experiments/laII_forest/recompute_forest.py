"""Recompute the La II forest three-leg table through the SHARED convention.

The numbers in the manuscript's forest table and in its thermal-emission (F8)
section were produced by compare.py, which predates sobolev.spectra and
carries two of the three normalization bugs this project has since found:

  * its red margin is 3952-3978 A, which straddles SEDONA's final spectrum
    bin -- that bin is partial, its flux collapses, and including it depresses
    the reference and inflates every band value (the bug of Section 4.9);
  * it band-averages by dividing by the NOMINAL band width rather than the
    span actually integrated.

The result was a forest table quoting SEDONA resolved = 0.3426 while the
separation section, computed through the shared helper, quoted 0.3435 for the
identical run -- the same measurement printed twice with different values two
pages apart.

This script recomputes all three legs, plus the F8 emission-off comparison,
through sobolev.spectra so that every forest number in the paper comes from
one pipeline. No new transport runs.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.optical_depth import stimulated_emission_factor
from sobolev.rays import RaySet
from sobolev.sobolev_leg import expansion_damp, resolved_attenuation, sobolev_attenuation
from sobolev.spectra import band_average, band_ratio

T_EXP, R_CORE, R_OUT, T_CORE, T_SHELL = 86400.0, 8.64e12, 2.592e13, 6000.0, 3000.0
V_D = 1.0e7
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)

d = np.load(HERE / "forest_lines.npz")
N_ION = float(d["n_ion"])
# stimulated-emission factor folded into pop, as SEDONA does (3e-6 here)
LINES = [(C / (lam * 1e-8), f, p * stimulated_emission_factor(C / (lam * 1e-8), T_SHELL))
         for lam, f, p in zip(d["lam"], d["f_lu"], d["pop"])]


def solver(t_shell):
    """Deterministic resolved leg, band-averaged the shared way."""
    nu = np.geomspace(7.50e14, 7.95e14, 1600)
    lum = emergent_luminosity(
        nu, LINES, lambda r: np.full_like(r, N_ION),
        lambda r: np.full_like(r, t_shell),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D, n_impact=150,
    )
    ratio = lum / (4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE))
    return band_average(C / nu * 1e8, ratio, BAND)


res = band_ratio(HERE / "run_bb" / "spectrum_1.dat", BAND, MARGIN,
                 R_CORE, T_CORE)
exp = band_ratio(HERE / "run_exp" / "spectrum_1.dat", BAND, MARGIN,
                 R_CORE, T_CORE)
sol_warm = solver(T_SHELL)
sol_cold = solver(0.0)

print("=== La II forest, T=3000 K, day 1, v_D=100 km/s, 3800-3955 A ===")
print(f"  deterministic solver (S = B_nu(T_shell))  {sol_warm:.4f}"
      f"   vs SEDONA {(sol_warm - res) / res:+.1%}")
print(f"  deterministic solver (T_shell -> 0)       {sol_cold:.4f}"
      f"   vs SEDONA {(sol_cold - res) / res:+.1%}")
print(f"  SEDONA resolved                           {res:.4f}")
print(f"  SEDONA expansion opacity                  {exp:.4f}")
print(f"  Delta_expansion                           {(exp - res) / res:+.1%}")

out = dict(res=res, exp=exp, solver_warm=sol_warm, solver_cold=sol_cold,
           d_exp=(exp - res) / res,
           d_solver_warm=(sol_warm - res) / res,
           d_solver_cold=(sol_cold - res) / res)

# --- the deterministic reference (referee Comment 5) ------------------------
# Sobolev and resolved legs on IDENTICAL rays, same source convention
# (T_shell -> 0), same populations, same transport treatment. The closed-form
# erf leg carries the first-order pair; the brute-force solver carries the
# frozen ("exact") and physical ("worldline" + dilution) pairs and also checks
# the erf leg. SEDONA's resolved run, at its seed mean, becomes the validation.
nu = np.geomspace(7.50e14, 7.95e14, 1600)
lam = C / nu * 1e8
rays = RaySet.midpoint(R_CORE, R_OUT, 200, n_env=0)
kw = dict(r_core=R_CORE, r_out=R_OUT, t_exp=T_EXP, n_ref=N_ION, rays=rays)
f_sob = band_average(lam, sobolev_attenuation(nu, LINES, relativity="first", **kw), BAND)
f_sob_cl = band_average(lam, sobolev_attenuation(nu, LINES, **kw), BAND)
f_exa = band_average(lam, sobolev_attenuation(nu, LINES, relativity="first", damp=expansion_damp, **kw), BAND)
f_det = band_average(lam, resolved_attenuation(nu, LINES, v_doppler=V_D, sweep="first", **kw), BAND)
f_det_cl = band_average(lam, resolved_attenuation(nu, LINES, v_doppler=V_D, sweep="classical", **kw), BAND)


def solver_on_rays(relativity, dilution=None):
    lum = emergent_luminosity(nu, LINES, lambda r: np.full_like(r, N_ION),
                              lambda r: np.full_like(r, 0.0), T_EXP, R_CORE, R_OUT,
                              T_CORE, V_D, relativity=relativity, dilution=dilution, rays=rays)
    return band_average(lam, lum / (4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)), BAND)


f_bf_first = solver_on_rays("first")
f_bf_exact = solver_on_rays("exact")
f_bf_world = solver_on_rays("worldline", dilution=True)
# Sobolev legs in the matching transport modes (exact/worldline change tau by
# <= 1% at beta_out = 0.01; the pair must share the mode)
f_sob_ex = band_average(lam, sobolev_attenuation(nu, LINES, relativity="exact", **kw), BAND)
f_sob_wl = band_average(lam, sobolev_attenuation(nu, LINES, relativity="worldline", **kw), BAND)

# SEDONA seed means: the 10 matched pairs of mc_noise/seeds.py when reduced
# (mc_noise_summary.json), else the original five of mc_noise/run.py.
summ = ROOT / "experiments/mc_noise/mc_noise_summary.json"
if summ.exists():
    sm = json.loads(summ.read_text())["tau5_vd100"]
    res_mean, res_sem, exp_mean = sm["f_bb_mean"], sm["f_bb_sem"], sm["f_exp_mean"]
    exp_sem, d_exp_paired, n_seeds = sm["f_exp_sem"], sm["d_exp_paired_mean"], sm["n_pairs"]
else:
    seeds = json.loads((ROOT / "experiments/mc_noise/mc_noise_results.json").read_text())
    res_mean, res_sem = seeds["tau5_vd100_bb"]["mean"], seeds["tau5_vd100_bb"]["std"] / np.sqrt(seeds["tau5_vd100_bb"]["n_seeds"])
    exp_mean, exp_sem, d_exp_paired, n_seeds = seeds["tau5_vd100_exp"]["mean"], None, None, 5

det = dict(
    rays="midpoint, 200 core", nu_points=1600, stim=True,
    f_sob_first=f_sob, f_sob_classical=f_sob_cl, f_exp_ana=f_exa,
    f_res_det_first=f_det, f_res_det_classical=f_det_cl,
    f_bf_first=f_bf_first, f_bf_exact=f_bf_exact, f_bf_worldline_dil=f_bf_world,
    f_sob_exact=f_sob_ex, f_sob_worldline=f_sob_wl,
    d_sob_det_first=(f_sob - f_det) / f_det,
    d_sob_det_classical=(f_sob_cl - f_det_cl) / f_det_cl,
    d_sob_det_exact=(f_sob_ex - f_bf_exact) / f_bf_exact,
    d_sob_det_worldline=(f_sob_wl - f_bf_world) / f_bf_world,
    d_exp_det_first=(f_exa - f_det) / f_det,
    erf_vs_bruteforce_first=(f_det - f_bf_first) / f_bf_first,
    sedona_res_seedmean=res_mean, sedona_res_sem=res_sem, sedona_exp_seedmean=exp_mean,
    sedona_exp_sem=exp_sem, sedona_n_seeds=n_seeds,
    d_sob_sedona_seedmean=(f_sob - res_mean) / res_mean,
    det_vs_sedona_seedmean=(f_det - res_mean) / res_mean,
    d_exp_sedona_seedmean=(d_exp_paired if d_exp_paired is not None else (exp_mean - res_mean) / res_mean),
)
out["deterministic"] = det
print("\n=== deterministic reference, identical rays, T_shell -> 0 ===")
print(f"  erf first {f_det:.5f}  vs brute-force first {f_bf_first:.5f}  ({100*det['erf_vs_bruteforce_first']:+.3f}%)")
print(f"  Delta_Sob^det  first {100*det['d_sob_det_first']:+.2f}%  classical {100*det['d_sob_det_classical']:+.2f}%  "
      f"exact {100*det['d_sob_det_exact']:+.2f}%  worldline+dilution {100*det['d_sob_det_worldline']:+.2f}%")
print(f"  Delta_exp^det (Poisson gap) {100*det['d_exp_det_first']:+.2f}%")
print(f"  SEDONA resolved seed mean {res_mean:.5f} +- {res_sem:.5f}: det leg vs it {100*det['det_vs_sedona_seedmean']:+.2f}%, "
      f"Sobolev vs it {100*det['d_sob_sedona_seedmean']:+.2f}%, Delta_exp(SEDONA) {100*det['d_exp_sedona_seedmean']:+.2f}%")
(HERE / "forest_table.json").write_text(json.dumps(out, indent=1))
print("\nwrote forest_table.json")
