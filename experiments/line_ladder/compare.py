"""Figure 5: the 2- and 20-line ladder -- Python solver vs SEDONA both modes.

Also prints the analytic Sobolev staircase: at observer frequency nu, a core
photon is attenuated by exp(-sum tau_k) over lines whose resonance plane lies
in front of the core along its path; the expansion-opacity prediction replaces
each tau_k by (1 - e^-tau_k).
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.constants import C, SIGMA_CLASSICAL
from sobolev.formal_transfer import emergent_luminosity, planck_bnu

T_EXP = 20 * 86400.0
R_CORE = 1.728e14
R_OUT = 5.184e14
T_CORE = 2.0e4
T_SHELL = 2000.0
V_D = 1.0e7
NU_REF = C / (12398.42 / 10.2 * 1e-8)  # reddest line

d = np.load("ladder_lines.npz")
N_H = float(d["n_h"])
CASES = {
    "2 lines, dv = 1500 km/s": (d["lam2"], d["f2"], "2line"),
    "20 lines, dv = 750 km/s": (d["lam20"], d["f20"], "20line"),
}


def load_sedona(run):
    s = np.loadtxt(f"run_{run}/spectrum_1.dat", comments="#")
    nu, lum = s[:, 0], s[:, 1]
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    v = (nu / NU_REF - 1.0) * C / 1e5
    red = (v < -500) & (v > -12000) & (lum > 0)
    return nu, lum / np.mean(lum[red] / cont[red]) / cont


def sobolev_staircase(nu_grid, lams, fs, damp, n_p=200):
    """p-averaged Sobolev attenuation of the core light.

    A core ray at impact parameter p crosses line k's resonance plane iff
    sqrt(r_core^2 - p^2) < z_res,k < sqrt(r_out^2 - p^2); its attenuation is
    exp(-sum damp(tau_k)) over crossed lines. Averaging exp(...) over p with
    weight p is essential: planes at z_res < r_core are crossed by only the
    outer core rays, and dropping them (a plane-counting shortcut) predicts a
    visibly too-shallow trough.
    """
    p = np.linspace(0.0, R_CORE, n_p, endpoint=False) + R_CORE / (2 * n_p)
    z_lo = np.sqrt(R_CORE**2 - p**2)
    z_hi = np.sqrt(R_OUT**2 - p**2)
    att = np.zeros((p.size, nu_grid.size))
    for lam, f in zip(lams, fs):
        nu_k = C / (lam * 1e-8)
        tau_k = SIGMA_CLASSICAL * f * N_H * (lam * 1e-8) * T_EXP
        z_res = C * T_EXP * (1.0 - nu_k / nu_grid)
        crossed = (z_res[None, :] > z_lo[:, None]) & (z_res[None, :] < z_hi[:, None])
        att += np.where(crossed, damp(tau_k), 0.0)
    return np.sum(p[:, None] * np.exp(-att), axis=0) / np.sum(p)


fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
for ax, (title, (lams, fs, tag)) in zip(axes, CASES.items()):
    nu_bb, r_bb = load_sedona(f"{tag}_bb")
    nu_exp, r_exp = load_sedona(f"{tag}_exp")

    nu_py = np.geomspace(2.34e15, 2.70e15, 700)
    lines = [(C / (lam * 1e-8), f) for lam, f in zip(lams, fs)]
    l_py = emergent_luminosity(
        nu_py, lines, lambda r: np.full_like(r, N_H),
        lambda r: np.full_like(r, T_SHELL),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D, n_impact=150,
    )
    r_py = l_py / (4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu_py, T_CORE))

    v = lambda nu: (np.asarray(nu) / NU_REF - 1.0) * C / 1e5
    ax.plot(v(nu_bb), r_bb, "C1", lw=0.9, label="SEDONA resolved bound-bound")
    ax.plot(v(nu_exp), r_exp, "C0", lw=0.9, label="SEDONA expansion opacity")
    ax.plot(v(nu_py), r_py, "k--", lw=1.3, label="Python formal solver")
    ax.plot(v(nu_py), sobolev_staircase(nu_py, lams, fs, lambda t: t),
            "r:", lw=1.2, label=r"analytic Sobolev $e^{-\Sigma\tau_k}$")
    ax.plot(v(nu_py), sobolev_staircase(nu_py, lams, fs, lambda t: 1 - np.exp(-t)),
            "b:", lw=1.2, label=r"expansion pred. $e^{-\Sigma(1-e^{-\tau_k})}$")
    ax.set_ylabel(r"$L_\nu / L_\nu^{\rm cont}$")
    ax.set_ylim(0, 1.35)
    ax.set_title(title, fontsize=10)
axes[0].legend(fontsize=7, ncol=2)
axes[1].set_xlabel("velocity from reddest line [km/s]")
axes[1].set_xlim(-5000, 21000)
fig.suptitle(r"Line ladder, $\tau_S$ = 0.5 per line, day 20", y=0.995)
fig.tight_layout()
out = Path(__file__).resolve().parents[2] / "outputs" / "fig5_line_ladder.png"
fig.savefig(out, dpi=200)
print("saved", out)

# Quantitative: mean depth in the saturated centre of the 20-line forest.
sel20 = (v(nu_py) > 5000) & (v(nu_py) < 12000)
print(f"20-line forest centre, python solver : {np.mean(r_py[sel20]):.4f}")
nu_bb, r_bb = load_sedona("20line_bb")
nu_exp, r_exp = load_sedona("20line_exp")
sb = (v(nu_bb) > 5000) & (v(nu_bb) < 12000)
se = (v(nu_exp) > 5000) & (v(nu_exp) < 12000)
print(f"20-line forest centre, sedona bb     : {np.mean(r_bb[sb]):.4f}")
print(f"20-line forest centre, sedona exp    : {np.mean(r_exp[se]):.4f}")
stair = sobolev_staircase(nu_py, d["lam20"], d["f20"], lambda t: t)
stair_e = sobolev_staircase(nu_py, d["lam20"], d["f20"], lambda t: 1 - np.exp(-t))
print(f"analytic Sobolev staircase           : {np.mean(stair[sel20]):.4f}")
print(f"analytic expansion prediction        : {np.mean(stair_e[sel20]):.4f}")
