"""Figure 9: La II + Ce II blended forest, three-way comparison."""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu

T_EXP = 86400.0
R_CORE = 8.64e12
R_OUT = 2.592e13
T_CORE = 6000.0
T_SHELL = 3000.0
V_D = 1.0e7

d = np.load("multiion_lines.npz")
RHO = float(d["rho"])
# per-line pop_frac carries species density per gram; n_l_of_r supplies rho.
LINES = [
    (C / (lam * 1e-8), f, pf)
    for lam, f, pf in zip(d["lam"], d["f_lu"], d["popfrac_per_rho"])
]


def load_sedona(run):
    s = np.loadtxt(f"{run}/spectrum_1.dat", comments="#")
    nu, lum = s[:, 0], s[:, 1]
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    lam = C / nu * 1e8
    red = (lam > 3952) & (lam < 3978) & (lum > 0)
    return lam, lum / np.mean(lum[red] / cont[red]) / cont


lam_bb, r_bb = load_sedona("run_bb")
lam_exp, r_exp = load_sedona("run_exp")

cache = Path("solver_result.npz")
if cache.exists():
    c = np.load(cache)
    nu_py, l_py = c["nu"], c["lum"]
else:
    nu_py = np.geomspace(7.50e14, 7.95e14, 1600)
    l_py = emergent_luminosity(
        nu_py, LINES, lambda r: np.full_like(r, RHO),
        lambda r: np.full_like(r, T_SHELL),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D, n_impact=150,
    )
r_py = l_py / (4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu_py, T_CORE))
lam_py = C / nu_py * 1e8

fig, ax = plt.subplots(figsize=(9, 4.8))
ax.plot(lam_bb, r_bb, "C1", lw=0.8, label="SEDONA resolved bound-bound")
ax.plot(lam_exp, r_exp, "C0", lw=0.8, label="SEDONA expansion opacity")
ax.plot(lam_py, r_py, "k--", lw=1.1, label="Python formal solver")
for lam, tau, z in zip(d["lam"], d["tau"], d["z"]):
    if tau > 0.1:
        ax.axvline(lam, color=("C3" if z == 57 else "C2"), lw=0.5, alpha=0.45)
ax.set_xlim(3790, 3980)
ax.set_ylim(0, 1.35)
ax.set_xlabel(r"wavelength [$\AA$]")
ax.set_ylabel(r"$L_\lambda / L_\lambda^{\rm cont}$")
n_la = int(((d["z"] == 57) & (d["tau"] > 0.1)).sum())
n_ce = int(((d["z"] == 58) & (d["tau"] > 0.1)).sum())
ax.set_title(
    f"La II + Ce II blend, 3850-3950 $\\AA$, T = 3000 K, day 1 "
    f"(red lines: La II $\\tau>0.1$ [{n_la}]; green: Ce II [{n_ce}])"
)
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
out = Path(__file__).resolve().parents[2] / "outputs" / "fig9_multiion_forest.png"
fig.savefig(out, dpi=200)
print("saved", out)


def band_flux(lam_A, ratio, lo=3800.0, hi=3955.0):
    m = (lam_A > lo) & (lam_A < hi)
    order = np.argsort(lam_A[m])
    return np.trapezoid(ratio[m][order], lam_A[m][order]) / (hi - lo)


fb, fe, fp = band_flux(lam_bb, r_bb), band_flux(lam_exp, r_exp), band_flux(lam_py, r_py)
print(f"band-averaged L/L_cont, python solver   : {fp:.4f}")
print(f"band-averaged L/L_cont, sedona resolved : {fb:.4f}")
print(f"band-averaged L/L_cont, sedona expansion: {fe:.4f}")
print(f"Delta_Sob (expansion vs resolved)       : {(fe - fb) / fb:+.2%}")
print(f"solver vs sedona resolved               : {(fp - fb) / fb:+.2%}")
