"""Figure 10: separating Sobolev-approximation error from expansion-opacity
implementation error.

Three predictions of the same physical configuration:
  resolved      -- SEDONA bound-bound (truth)
  Sobolev       -- per-line exp(-tau_S), delta-function resonances
  expansion     -- per-crossing exp(-(1-e^-tau_S))

Convention (F8): the analytic legs carry no thermal emission, so they are
compared against SEDONA's resolved mode, which likewise deposits absorbed
energy without re-emitting at fixed temperature.

Note that the Sobolev leg has NO v_D dependence by construction -- a
delta-function resonance has no width. Any v_D dependence in the resolved
result is therefore un-modelled by the Sobolev approximation itself, which
is exactly what the sweep panel measures.
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.formal_transfer import planck_bnu
from sobolev.sobolev_leg import expansion_damp, sobolev_attenuation

T_EXP, R_CORE, R_OUT, T_CORE = 86400.0, 8.64e12, 2.592e13, 6000.0
FOREST = ROOT / "experiments/laII_forest"
BLEND = ROOT / "experiments/multiion"
BAND = (3800.0, 3955.0)


def band_avg(lam, ratio):
    m = (lam > BAND[0]) & (lam < BAND[1])
    o = np.argsort(lam[m])
    return np.trapezoid(ratio[m][o], lam[m][o]) / (BAND[1] - BAND[0])


def sedona_band(run_dir):
    s = np.loadtxt(Path(run_dir) / "spectrum_1.dat", comments="#")
    nu, lum = s[:, 0], s[:, 1]
    cont = 4 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    lam = C / nu * 1e8
    red = (lam > 3952) & (lam < 3978) & (lum > 0)
    return lam, lum / np.mean(lum[red] / cont[red]) / cont


def analytic_legs(lines, n_ref, nu_grid):
    sob = sobolev_attenuation(nu_grid, lines, R_CORE, R_OUT, T_EXP, n_ref)
    exp_ = sobolev_attenuation(
        nu_grid, lines, R_CORE, R_OUT, T_EXP, n_ref, damp=expansion_damp
    )
    return sob, exp_


nu_grid = np.geomspace(7.50e14, 7.95e14, 4000)
lam_grid = C / nu_grid * 1e8

# ---------------- single-ion La II forest ----------------
d = np.load(FOREST / "forest_lines.npz")
N_ION = float(d["n_ion"])
lines_la = [
    (C / (l * 1e-8), f, p) for l, f, p in zip(d["lam"], d["f_lu"], d["pop"])
]
sob_la, exp_la = analytic_legs(lines_la, N_ION, nu_grid)
lam_bb, r_bb = sedona_band(FOREST / "run_bb")
lam_ex, r_ex = sedona_band(FOREST / "run_exp")

F_res = band_avg(lam_bb, r_bb)
F_sob = band_avg(lam_grid, sob_la)
F_exp_sed = band_avg(lam_ex, r_ex)
F_exp_ana = band_avg(lam_grid, exp_la)

print("=== La II forest (T=3000 K, day 1, v_D=100 km/s) ===")
print(f"  SEDONA resolved (truth)     {F_res:.4f}")
print(f"  Sobolev proper              {F_sob:.4f}   Delta = {(F_sob-F_res)/F_res:+.1%}")
print(f"  expansion (analytic)        {F_exp_ana:.4f}   Delta = {(F_exp_ana-F_res)/F_res:+.1%}")
print(f"  expansion (SEDONA)          {F_exp_sed:.4f}   Delta = {(F_exp_sed-F_res)/F_res:+.1%}")

# ---------------- multi-ion blend ----------------
db = np.load(BLEND / "multiion_lines.npz")
RHO = float(db["rho"])
lines_mi = [
    (C / (l * 1e-8), f, p)
    for l, f, p in zip(db["lam"], db["f_lu"], db["popfrac_per_rho"])
]
sob_mi, exp_mi = analytic_legs(lines_mi, RHO, nu_grid)
lam_bb_m, r_bb_m = sedona_band(BLEND / "run_bb")
lam_ex_m, r_ex_m = sedona_band(BLEND / "run_exp")
Fm_res, Fm_sob = band_avg(lam_bb_m, r_bb_m), band_avg(lam_grid, sob_mi)
Fm_exp = band_avg(lam_ex_m, r_ex_m)
print("=== La II + Ce II blend ===")
print(f"  SEDONA resolved (truth)     {Fm_res:.4f}")
print(f"  Sobolev proper              {Fm_sob:.4f}   Delta = {(Fm_sob-Fm_res)/Fm_res:+.1%}")
print(f"  expansion (SEDONA)          {Fm_exp:.4f}   Delta = {(Fm_exp-Fm_res)/Fm_res:+.1%}")

# ---------------- sweep points: tau_max x v_D ----------------
sweep = json.loads((FOREST / "sweep_results.json").read_text())
rows = []
for r in sweep:
    if not (r.get("bb") and r.get("exp")):
        continue
    n_ref = N_ION * (r["tau_max"] / 5.0)
    sob, exp_ana = analytic_legs(lines_la, n_ref, nu_grid)
    f_sob, f_exp_ana = band_avg(lam_grid, sob), band_avg(lam_grid, exp_ana)
    rows.append(
        dict(
            tau_max=r["tau_max"], v_d=r["v_d_kms"], f_res=r["bb"],
            f_exp=r["exp"], f_sob=f_sob, f_exp_ana=f_exp_ana,
            d_sob=(f_sob - r["bb"]) / r["bb"],
            d_exp=(r["exp"] - r["bb"]) / r["bb"],
            d_exp_ana=(f_exp_ana - r["bb"]) / r["bb"],
        )
    )
    print(
        f"tau_max={r['tau_max']:5g} v_D={r['v_d_kms']:5g}  "
        f"res={r['bb']:.4f} sob={f_sob:.4f} exp_ana={f_exp_ana:.4f} "
        f"exp_sed={r['exp']:.4f} | D_sob={rows[-1]['d_sob']:+7.1%} "
        f"D_exp_ana={rows[-1]['d_exp_ana']:+7.1%} D_exp_sed={rows[-1]['d_exp']:+7.1%}"
    )

(Path(__file__).parent / "separation_results.json").write_text(
    json.dumps(
        dict(
            la_forest=dict(res=F_res, sob=F_sob, exp=F_exp_sed),
            blend=dict(res=Fm_res, sob=Fm_sob, exp=Fm_exp),
            sweep=rows,
        ),
        indent=1,
    )
)

# ---------------- figure ----------------
fig = plt.figure(figsize=(11, 7.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1])

ax = fig.add_subplot(gs[0, :])
ax.plot(lam_bb, r_bb, "C1", lw=0.9, label="SEDONA resolved (truth)")
ax.plot(lam_grid, sob_la, "C2", lw=1.0, label=r"Sobolev proper: $e^{-\tau_S}$/line")
ax.plot(lam_grid, exp_la, "C0", lw=1.0,
        label=r"expansion: $e^{-(1-e^{-\tau_S})}$/crossing")
ax.set_xlim(3790, 3980)
ax.set_ylim(0, 1.3)
ax.set_xlabel(r"wavelength [$\AA$]")
ax.set_ylabel(r"$L_\lambda/L_\lambda^{\rm cont}$")
ax.set_title(
    "La II forest: the Sobolev approximation vs its expansion-opacity "
    "implementation (T = 3000 K, day 1)"
)
ax.legend(fontsize=8, loc="lower left")

ax2 = fig.add_subplot(gs[1, 0])
for tau, color in [(0.5, "C2"), (5.0, "C1"), (50.0, "C3")]:
    sel = sorted([r for r in rows if r["tau_max"] == tau], key=lambda r: r["v_d"])
    ax2.semilogx([r["v_d"] for r in sel], [100 * r["d_exp"] for r in sel],
                 "o-", color=color, label=rf"$\tau_{{\max}}$={tau:g}")
    ax2.semilogx([r["v_d"] for r in sel], [100 * r["d_sob"] for r in sel],
                 "s--", color=color, alpha=0.6)
ax2.set_xlabel(r"$v_D$ [km/s]")
ax2.set_ylabel(r"$\Delta$ [%]")
ax2.set_title("solid: expansion   dashed: Sobolev proper", fontsize=9)
ax2.axhline(0, color="k", lw=0.6)
ax2.legend(fontsize=7)
ax2.grid(alpha=0.3)

ax3 = fig.add_subplot(gs[1, 1])
width = 0.35
taus = [0.5, 5.0, 50.0]
at_100 = [next(r for r in rows if r["tau_max"] == t and r["v_d"] == 100.0) for t in taus]
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
out = ROOT / "outputs" / "fig10_sobolev_vs_expansion.png"
fig.savefig(out, dpi=200)
print("saved", out)
