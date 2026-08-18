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
from sobolev.sobolev_leg import sobolev_attenuation, expansion_damp
from sobolev.spectra import band_average, band_ratio

FOREST = ROOT / "experiments/laII_forest"
T_EXP, R_CORE, R_OUT, T_CORE = 86400.0, 8.64e12, 2.592e13, 6000.0
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)

d = np.load(FOREST / "forest_lines.npz")
N_ION = float(d["n_ion"])
LINES = [(C / (l * 1e-8), f, p) for l, f, p in zip(d["lam"], d["f_lu"], d["pop"])]
nu = np.geomspace(7.50e14, 7.95e14, 4000)
lam_an = C / nu * 1e8


def analytic(scale, damp=None):
    att = sobolev_attenuation(nu, LINES, R_CORE, R_OUT, T_EXP, N_ION * scale,
                              damp=damp)
    return band_average(lam_an, att, BAND)


def sedona(run):
    return band_ratio(FOREST / run / "spectrum_1.dat", BAND, MARGIN,
                      R_CORE, T_CORE)


rows = []
old = json.loads((FOREST / "sweep_results.json").read_text())

print("=== tau_max x v_D grid ===")
print("tau_max  v_D | D_Sob new    old |  D_exp new    old")
for tau in (0.5, 5.0, 50.0):
    f_sob = analytic(tau / 5.0)
    f_exp_ana = analytic(tau / 5.0, damp=expansion_damp)
    for vd in (10, 30, 100, 300):
        fb, fe = sedona(f"sweep_tau{tau:g}_vd{vd:g}_bb"), sedona(f"sweep_tau{tau:g}_vd{vd:g}_exp")
        o = [x for x in old if x["tau_max"] == tau and x["v_d_kms"] == vd][0]
        d_sob, d_exp = (f_sob - fb) / fb, (fe - fb) / fb
        d_sob_old = (f_sob - o["bb"]) / o["bb"]
        rows.append(dict(tau_max=tau, v_d=vd, f_res=fb, f_exp=fe, f_sob=f_sob,
                         f_exp_ana=f_exp_ana, d_sob=d_sob, d_exp=d_exp))
        print(f"{tau:6g} {vd:4d} | {d_sob:+8.1%} {d_sob_old:+8.1%} |"
              f" {d_exp:+8.1%} {o['delta_sob']:+8.1%}")

print("\n=== thermal-width frontier (tau_max = 5) ===")
print(" v_D | F_resolved  D_Sobolev   D_expansion")
frontier = []
f_sob5 = analytic(1.0)
for vd, rb, re in ((1, "tsweep_T3000_vd1_bb", "tsweep_T3000_vd1_exp"),
                   (3, "tsweep_T3000_vd3_bb", "tsweep_T3000_vd3_exp"),
                   (10, "sweep_tau5_vd10_bb", "sweep_tau5_vd10_exp"),
                   (30, "sweep_tau5_vd30_bb", "sweep_tau5_vd30_exp"),
                   (100, "sweep_tau5_vd100_bb", "sweep_tau5_vd100_exp"),
                   (300, "sweep_tau5_vd300_bb", "sweep_tau5_vd300_exp")):
    fb, fe = sedona(rb), sedona(re)
    frontier.append(dict(v_d=vd, f_res=fb, f_exp=fe, f_sob=f_sob5,
                         d_sob=(f_sob5 - fb) / fb, d_exp=(fe - fb) / fb))
    print(f"{vd:4d} | {fb:10.4f} {(f_sob5-fb)/fb:+10.1%} {(fe-fb)/fb:+12.1%}")

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
    ax2.semilogx([r["v_d"] for r in sel], [100 * r["d_sob"] for r in sel],
                 "s--", color=color, alpha=0.6)
# The frontier extends the tau_max = 5 curves down to the thermal end.
fr = sorted(frontier, key=lambda r: r["v_d"])
ax2.semilogx([r["v_d"] for r in fr], [100 * r["d_exp"] for r in fr],
             "o-", color="C1", alpha=0.5, lw=0.8)
ax2.semilogx([r["v_d"] for r in fr], [100 * r["d_sob"] for r in fr],
             "s--", color="C1", alpha=0.35, lw=0.8)
ax2.set_xlabel(r"$v_D$ [km/s]")
ax2.set_ylabel(r"$\Delta$ [%]")
ax2.set_title("solid: expansion   dashed: Sobolev proper", fontsize=9)
ax2.axhline(0, color="k", lw=0.6)
ax2.legend(fontsize=7)
ax2.grid(alpha=0.3)

ax3 = fig.add_subplot(gs[1, 1])
width, taus = 0.35, [0.5, 5.0, 50.0]
at_100 = [next(r for r in rows if r["tau_max"] == t and r["v_d"] == 100)
          for t in taus]
x = np.arange(len(taus))
ax3.bar(x - width / 2, [100 * r["d_sob"] for r in at_100], width,
        color="C2", label="Sobolev proper")
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
