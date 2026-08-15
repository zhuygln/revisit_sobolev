"""Figure 4: minimal 1-line model -- Python formal solver vs SEDONA (both modes).

Run from experiments/minimal_1line/ after make_model.py and both SEDONA runs.
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

# --- identical physics to make_model.py / param_bb.lua ---
T_EXP = 20 * 86400.0
R_CORE = 1.728e14
R_OUT = 5.184e14
T_CORE = 2.0e4
T_SHELL = 2000.0
N_H = 5.3975
F_LU = 0.6647
LAMBDA0 = 12398.42 / 10.2 * 1e-8
NU0 = C / LAMBDA0
V_D = 1.0e7  # = SEDONA line_velocity_width
TAU_S = 2.0


def load_sedona(path):
    d = np.loadtxt(path, comments="#")
    return d[:, 0], d[:, 1]


nu_bb, l_bb = load_sedona("run_bb/spectrum_1.dat")
nu_exp, l_exp = load_sedona("run_exp/spectrum_1.dat")


def rescale_to_continuum(nu, lum):
    """SEDONA's lightbulb pours the WHOLE core luminosity into the transport
    frequency window, so its absolute L_nu exceeds the analytic continuum by
    1/f_window (~15 here). Rescale using the red side (v < -500 km/s), which
    is pure continuum: no resonance plane there intersects a core ray, and
    the cold shell adds no emission."""
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    v = (nu / NU0 - 1.0) * C / 1e5
    red = (v < -500) & (v > -25000) & (lum > 0)
    scale = np.mean(lum[red] / cont[red])
    return lum / scale, cont


l_bb, cont_bb = rescale_to_continuum(nu_bb, l_bb)
l_exp, cont_exp = rescale_to_continuum(nu_exp, l_exp)

nu_py = np.geomspace(2.30e15, 2.62e15, 400)
l_py = emergent_luminosity(
    nu_py,
    [(NU0, F_LU)],
    lambda r: np.full_like(r, N_H),
    lambda r: np.full_like(r, T_SHELL),
    T_EXP, R_CORE, R_OUT, T_CORE, V_D,
    n_impact=200,
)
l_cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu_py, T_CORE)

fig, (ax, ax2) = plt.subplots(
    2, 1, figsize=(7.5, 6.5), sharex=True, height_ratios=[2, 1]
)
v_bb = (nu_bb / NU0 - 1.0) * C / 1e5
v_exp = (nu_exp / NU0 - 1.0) * C / 1e5
v_py = (nu_py / NU0 - 1.0) * C / 1e5

ax.plot(v_bb, l_bb, "C1", lw=1.0, label="SEDONA resolved bound-bound")
ax.plot(v_exp, l_exp, "C0", lw=1.0, label="SEDONA expansion opacity (Sobolev)")
ax.plot(v_py, l_py, "k--", lw=1.4, label="Python formal solver")
ax.plot(v_py, l_cont, color="gray", lw=0.8, ls=":", label="continuum")
ax.set_ylabel(r"$L_\nu$ [erg s$^{-1}$ Hz$^{-1}$]")
ax.set_title(
    f"Minimal 1-line model: fake Ly$\\alpha$, $\\tau_S$ = {TAU_S:.0f}, day 20"
)
ax.legend(fontsize=8)

# Ratio panel against the analytic continuum, with the Sobolev trough marked.
ax2.plot(v_bb, l_bb / cont_bb, "C1", lw=1.0)
ax2.plot(v_exp, l_exp / cont_exp, "C0", lw=1.0)
ax2.plot(v_py, l_py / l_cont, "k--", lw=1.4)
ax2.axhline(np.exp(-TAU_S), color="red", ls=":", lw=1,
            label=r"$e^{-\tau_S}$ = " + f"{np.exp(-TAU_S):.3f}")
# The known single-line failure of the expansion-opacity treatment: a photon
# crossing one resonance is attenuated by exp(-(1 - e^-tau_S)), not e^-tau_S.
exp_pred = np.exp(-(1.0 - np.exp(-TAU_S)))
ax2.axhline(exp_pred, color="C0", ls=":", lw=1,
            label=r"$e^{-(1-e^{-\tau_S})}$ = " + f"{exp_pred:.3f}")
ax2.set_xlabel("velocity from line centre [km/s]")
ax2.set_ylabel(r"$L_\nu / L_\nu^{\rm cont}$")
ax2.set_ylim(0, 1.3)
ax2.legend(fontsize=8)
# The transport-grid edges (~+-20000 km/s) carry MC binning artifacts; the
# physics lives in the +-8000 km/s neighbourhood of the line.
ax.set_xlim(-8000, 8000)
fig.tight_layout()
out = Path(__file__).resolve().parents[2] / "outputs" / "fig4_minimal_1line_threeway.png"
fig.savefig(out, dpi=200)
print("saved", out)

# Quantitative trough comparison. The shell spans 1000-3000 km/s, and the
# trough saturates at exp(-tau_S) only where the resonance plane lies fully
# in front of the core (v > r_core/t = 1000 km/s), away from both edges:
sel = lambda v: (np.asarray(v) > 1400) & (np.asarray(v) < 2600)
print(f"analytic exp(-tau_S)      = {np.exp(-TAU_S):.4f}")
print(f"python solver trough      = {np.mean((l_py/l_cont)[sel(v_py)]):.4f}")
print(f"sedona bb trough          = {np.mean((l_bb/cont_bb)[sel(v_bb)]):.4f}")
print(f"sedona expansion trough   = {np.mean((l_exp/cont_exp)[sel(v_exp)]):.4f}")
