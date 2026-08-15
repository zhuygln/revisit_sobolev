"""Figure 6: real GSI La II forest, 3850-3950 A -- three-way comparison.

Also computes the first validity-map datapoint: the frequency-integrated flux
error of the expansion-opacity treatment relative to the resolved one.
"""

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

d = np.load("forest_lines.npz")
N_ION = float(d["n_ion"])
LINES = [
    (C / (lam * 1e-8), f, p) for lam, f, p in zip(d["lam"], d["f_lu"], d["pop"])
]


def load_sedona(run):
    s = np.loadtxt(f"{run}/spectrum_1.dat", comments="#")
    nu, lum = s[:, 0], s[:, 1]
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    lam_A = C / nu * 1e8
    red = (lam_A > 3952) & (lam_A < 3978) & (lum > 0)  # line-free red margin
    return lam_A, lum / np.mean(lum[red] / cont[red]) / cont


lam_bb, r_bb = load_sedona("run_bb")
lam_exp, r_exp = load_sedona("run_exp")

nu_py = np.geomspace(7.50e14, 7.95e14, 1600)
l_py = emergent_luminosity(
    nu_py, LINES, lambda r: np.full_like(r, N_ION),
    lambda r: np.full_like(r, T_SHELL),
    T_EXP, R_CORE, R_OUT, T_CORE, V_D, n_impact=150,
)
r_py = l_py / (4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu_py, T_CORE))
lam_py = C / nu_py * 1e8

fig, ax = plt.subplots(figsize=(9, 4.8))
ax.plot(lam_bb, r_bb, "C1", lw=0.8, label="SEDONA resolved bound-bound")
ax.plot(lam_exp, r_exp, "C0", lw=0.8, label="SEDONA expansion opacity")
ax.plot(lam_py, r_py, "k--", lw=1.1, label="Python formal solver")
for lam, tau in zip(d["lam"], d["tau"]):
    if tau > 0.1:
        ax.axvline(lam, color="gray", lw=0.5, alpha=0.4)
ax.set_xlim(3790, 3980)
ax.set_ylim(0, 1.35)
ax.set_xlabel(r"wavelength [$\AA$]")
ax.set_ylabel(r"$L_\lambda / L_\lambda^{\rm cont}$")
ax.set_title(
    "GSI La II forest 3850-3950 $\\AA$: 153 lines, T = 3000 K, day 1, "
    "$v_D$ = 100 km/s (grey lines: $\\tau_S$ > 0.1)"
)
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
out = Path(__file__).resolve().parents[2] / "outputs" / "fig6_laII_forest.png"
fig.savefig(out, dpi=200)
print("saved", out)

# Validity-map datapoint: integrated flux error over the absorbed band.
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
