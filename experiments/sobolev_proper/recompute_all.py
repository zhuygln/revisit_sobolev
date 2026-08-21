"""Recompute every Sobolev-vs-expansion number through the shared band-ratio
helper, and report what moved.

The previous values compared a correctly-normalized analytic leg against a
SEDONA band flux normalized by RAW LUMINOSITY, leaving the Planck slope in
the answer. That produced a spurious ~5-7% "v_D-independent Sobolev floor".
Delta_expansion, being a same-code differential, should barely move -- that
invariance is the check that this helper changed only what it should.

No new transport runs: every spectrum is already on disk.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.optical_depth import stimulated_emission_factor
from sobolev.rays import RaySet
from sobolev.sobolev_leg import sobolev_attenuation, expansion_damp, resolved_attenuation
from sobolev.spectra import band_average, band_ratio

FOREST = ROOT / "experiments/laII_forest"
T_EXP, R_CORE, R_OUT, T_CORE = 86400.0, 8.64e12, 2.592e13, 6000.0
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)

d = np.load(FOREST / "forest_lines.npz")
N_ION = float(d["n_ion"])
T_SHELL = 3000.0
# stimulated-emission factor folded into pop, as SEDONA does (3e-6 here)
LINES = [(C / (l * 1e-8), f, p * stimulated_emission_factor(C / (l * 1e-8), T_SHELL))
         for l, f, p in zip(d["lam"], d["f_lu"], d["pop"])]
nu = np.geomspace(7.50e14, 7.95e14, 4000)
lam_an = C / nu * 1e8
# Both analytic legs and the deterministic reference on IDENTICAL rays
# (sobolev/rays.py), in the first-order transport mode the closed-form
# resolved leg is exact in.
RAYS = RaySet.midpoint(R_CORE, R_OUT, 400)


def analytic(scale, damp=None):
    att = sobolev_attenuation(nu, LINES, R_CORE, R_OUT, T_EXP, N_ION * scale,
                              damp=damp, relativity="first", rays=RAYS)
    return band_average(lam_an, att, BAND)


def deterministic(scale, v_d_kms):
    """Closed-form resolved leg (Gaussian, uniform n_l, first-order Doppler):
    the noise-free reference for Delta_Sob, free at any v_D."""
    att = resolved_attenuation(nu, LINES, R_CORE, R_OUT, T_EXP, N_ION * scale,
                               v_d_kms * 1e5, rays=RAYS, sweep="first")
    return band_average(lam_an, att, BAND)


def sedona(run):
    return band_ratio(FOREST / run / "spectrum_1.dat", BAND, MARGIN,
                      R_CORE, T_CORE)


# Seed-matched pairs (mc_noise/seeds.py -> analyze.py): where a grid or
# frontier point has them, its SEDONA fluxes are the seed means and its
# Delta_exp the PAIRED mean, so the single production realization no longer
# carries the number.
SUMM = ROOT / "experiments/mc_noise/mc_noise_summary.json"
SEEDS = json.loads(SUMM.read_text()) if SUMM.exists() else {}


def sedona_pair(tau, vd, rb, re):
    """(f_res, f_exp, d_exp, n_pairs) -- seed means if available, else single."""
    key = f"tau{tau:g}_vd{vd:g}"
    sm = SEEDS.get(key)
    if sm and sm.get("n_pairs", 0) >= 2:
        return sm["f_bb_mean"], sm["f_exp_mean"], sm["d_exp_paired_mean"], sm["n_pairs"]
    fb, fe = sedona(rb), sedona(re)
    return fb, fe, (fe - fb) / fb, 1


rows = []
old = json.loads((FOREST / "sweep_results.json").read_text())

print("=== tau_max x v_D grid ===")
print("tau_max  v_D | D_Sob(SEDONA ref)  D_Sob^det  det-vs-SEDONA |  D_exp(SEDONA)  D_exp^det")
for tau in (0.5, 5.0, 50.0):
    f_sob = analytic(tau / 5.0)
    f_exp_ana = analytic(tau / 5.0, damp=expansion_damp)
    for vd in (10, 30, 100, 300):
        fb, fe, d_exp, npair = sedona_pair(tau, vd, f"sweep_tau{tau:g}_vd{vd:g}_bb", f"sweep_tau{tau:g}_vd{vd:g}_exp")
        f_det = deterministic(tau / 5.0, vd)
        d_sob = (f_sob - fb) / fb
        d_sob_det, d_exp_det = (f_sob - f_det) / f_det, (f_exp_ana - f_det) / f_det
        rows.append(dict(tau_max=tau, v_d=vd, f_res=fb, f_exp=fe, f_sob=f_sob,
                         f_exp_ana=f_exp_ana, f_res_det=f_det, d_sob=d_sob, d_exp=d_exp,
                         d_sob_det=d_sob_det, d_exp_det=d_exp_det, n_seed_pairs=npair,
                         det_vs_sedona=(f_det - fb) / fb))
        print(f"{tau:6g} {vd:4d} | {d_sob:+9.2%}        {d_sob_det:+8.2%}   {(f_det-fb)/fb:+8.2%}    |"
              f" {d_exp:+9.1%}      {d_exp_det:+8.1%}")

print("\n=== thermal-width frontier (tau_max = 5) ===")
print(" v_D | F_res(SEDONA) F_res^det | D_Sob(SEDONA) D_Sob^det | D_exp(SEDONA) D_exp^det")
frontier = []
f_sob5 = analytic(1.0)
f_exa5 = analytic(1.0, damp=expansion_damp)
for vd, rb, re in ((1, "tsweep_T3000_vd1_bb", "tsweep_T3000_vd1_exp"),
                   (3, "tsweep_T3000_vd3_bb", "tsweep_T3000_vd3_exp"),
                   (10, "sweep_tau5_vd10_bb", "sweep_tau5_vd10_exp"),
                   (30, "sweep_tau5_vd30_bb", "sweep_tau5_vd30_exp"),
                   (100, "sweep_tau5_vd100_bb", "sweep_tau5_vd100_exp"),
                   (300, "sweep_tau5_vd300_bb", "sweep_tau5_vd300_exp")):
    fb, fe, d_exp_f, npair = sedona_pair(5.0, vd, rb, re)
    f_det = deterministic(1.0, vd)
    frontier.append(dict(v_d=vd, f_res=fb, f_exp=fe, f_sob=f_sob5, f_res_det=f_det,
                         d_sob=(f_sob5 - fb) / fb, d_exp=d_exp_f, n_seed_pairs=npair,
                         d_sob_det=(f_sob5 - f_det) / f_det, d_exp_det=(f_exa5 - f_det) / f_det,
                         det_vs_sedona=(f_det - fb) / fb))
    print(f"{vd:4d} | {fb:10.4f}  {f_det:8.4f} | {(f_sob5-fb)/fb:+10.2%} {(f_sob5-f_det)/f_det:+8.2%} |"
          f" {d_exp_f:+10.1%} {(f_exa5-f_det)/f_det:+8.1%}")

print("\n=== temperature axis (v_D = 100, strength pinned) ===")
print("   T | F_resolved  D_expansion   old")
# The T axis fed the validity map from tsweep_results.json, which was written
# with the raw-luminosity normalization. D_exp is a same-code differential so
# the bias largely cancels, but "largely" is not "exactly", and having one
# table in the paper sourced from a different pipeline than the table beside
# it is how the two disagreed by up to 2.7 points. Recompute it here.
old_t = [x for x in old_all if x.get("axis") == "T"] if (
    old_all := json.loads((FOREST / "tsweep_results.json").read_text())) else []
t_axis = []
for t in (2500, 3000, 4000, 5000):
    fb_t, fe_t = sedona(f"tsweep_T{t}_bb"), sedona(f"tsweep_T{t}_exp")
    o = [x for x in old_t if x["T"] == float(t)]
    d = (fe_t - fb_t) / fb_t
    t_axis.append(dict(T=t, v_d=100, f_res=fb_t, f_exp=fe_t, d_exp=d))
    print(f"{t:5d} | {fb_t:10.4f} {d:+12.1%} {o[0]['delta_sob']:+8.1%}"
          if o else f"{t:5d} | {fb_t:10.4f} {d:+12.1%}")

# headline forest condition
fb = sedona("run_bb"); fe = sedona("run_exp")
print(f"\n=== La II forest headline (v_D = 100 km/s, tau_max = 5) ===")
print(f"  SEDONA resolved  {fb:.4f}")
print(f"  Sobolev proper   {f_sob5:.4f}   Delta = {(f_sob5-fb)/fb:+.1%}")
print(f"  expansion        {fe:.4f}   Delta = {(fe-fb)/fb:+.1%}")

(Path(__file__).parent / "separation_results.json").write_text(json.dumps(
    dict(grid=rows, frontier=frontier, t_axis=t_axis,
         headline=dict(res=fb, sob=f_sob5, exp=fe)), indent=1))
print("\nwrote separation_results.json")


# ---------------------------------------------------------------------------
# Figure, generated HERE from the numbers just computed.
#
# It used to live in compare.py, which recomputed everything from
# sweep_results.json -- the pre-fix file normalized by raw luminosity. So when
# this script corrected the numbers, the figure kept plotting the retracted
# ~5-8% "v_D-independent Sobolev floor", and sat in the manuscript directly
# opposite the paragraph retracting it. Nothing flagged the contradiction
# because the figure had its own private copy of the data.
#
# The figure is therefore built from `rows`/`frontier` in this process. There
# is no path by which it can disagree with separation_results.json again.
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

lam_grid = lam_an
sob_la = sobolev_attenuation(nu, LINES, R_CORE, R_OUT, T_EXP, N_ION)
exp_la = sobolev_attenuation(nu, LINES, R_CORE, R_OUT, T_EXP, N_ION,
                             damp=expansion_damp)
_, lam_bb, r_bb = band_ratio(FOREST / "run_bb" / "spectrum_1.dat", BAND,
                             MARGIN, R_CORE, T_CORE, return_spectrum=True)

fig = plt.figure(figsize=(11, 7.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1])

ax = fig.add_subplot(gs[0, :])
ax.plot(lam_bb, r_bb, "C1", lw=0.9, label="SEDONA resolved (truth)")
ax.plot(lam_grid, sob_la, "C2", lw=1.0,
        label=r"Sobolev proper: $e^{-\tau_S}$/line")
ax.plot(lam_grid, exp_la, "C0", lw=1.0,
        label=r"expansion: $e^{-(1-e^{-\tau_S})}$/crossing")
ax.set_xlim(3790, 3980)
ax.set_ylim(0, 1.3)
ax.set_xlabel(r"wavelength [$\AA$]")
ax.set_ylabel(r"$L_\lambda/L_\lambda^{\rm cont}$")
ax.set_title("La II forest: the Sobolev approximation vs its "
             "expansion-opacity implementation (T = 3000 K, day 1)")
ax.legend(fontsize=8, loc="lower left")

ax2 = fig.add_subplot(gs[1, 0])
for tau, color in [(0.5, "C2"), (5.0, "C1"), (50.0, "C3")]:
    sel = sorted([r for r in rows if r["tau_max"] == tau],
                 key=lambda r: r["v_d"])
    ax2.semilogx([r["v_d"] for r in sel], [100 * r["d_exp"] for r in sel],
                 "o-", color=color, label=rf"$\tau_{{\max}}$={tau:g}")
    ax2.semilogx([r["v_d"] for r in sel], [100 * r["d_sob_det"] for r in sel],
                 "s--", color=color, alpha=0.6)
# The frontier extends the tau_max = 5 curves down to the thermal end.
fr = sorted(frontier, key=lambda r: r["v_d"])
ax2.semilogx([r["v_d"] for r in fr], [100 * r["d_exp"] for r in fr],
             "o-", color="C1", alpha=0.5, lw=0.8)
ax2.semilogx([r["v_d"] for r in fr], [100 * r["d_sob_det"] for r in fr],
             "s--", color="C1", alpha=0.35, lw=0.8)
ax2.set_xlabel(r"$v_D$ [km/s]")
ax2.set_ylabel(r"$\Delta$ [%]")
ax2.set_title("solid: expansion (SEDONA)   dashed: Sobolev proper vs deterministic reference", fontsize=8)
ax2.axhline(0, color="k", lw=0.6)
ax2.legend(fontsize=7)
ax2.grid(alpha=0.3)

ax3 = fig.add_subplot(gs[1, 1])
width, taus = 0.35, [0.5, 5.0, 50.0]
at_100 = [next(r for r in rows if r["tau_max"] == t and r["v_d"] == 100)
          for t in taus]
x = np.arange(len(taus))
ax3.bar(x - width / 2, [100 * r["d_sob_det"] for r in at_100], width,
        color="C2", label="Sobolev proper (det. ref.)")
ax3.bar(x + width / 2, [100 * r["d_exp"] for r in at_100], width,
        color="C0", label="expansion opacity")
ax3.set_xticks(x)
ax3.set_xticklabels([f"{t:g}" for t in taus])
ax3.set_xlabel(r"$\tau_{\max}$")
ax3.set_ylabel(r"$\Delta$ [%]")
ax3.set_title(r"at $v_D$ = 100 km/s", fontsize=9)
ax3.axhline(0, color="k", lw=0.6)
ax3.legend(fontsize=7)
ax3.grid(alpha=0.3, axis="y")

fig.tight_layout()
for out in (ROOT / "outputs" / "fig10_sobolev_vs_expansion.png",
            ROOT / "docs/figures" / "fig10_sobolev_vs_expansion.png"):
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print("saved", out)
