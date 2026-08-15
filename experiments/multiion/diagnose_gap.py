"""Remaining structural suspects for the solver-vs-SEDONA gap in the blend.

Both hypotheses so far are dead: Voigt wings (3e-5 effect) and z-resolution
(flat over 16x). What is left that differs BETWEEN the two codes and has the
right sign (solver brighter = absorbing less, or emitting more)?

  A) shell thermal emission. The solver adds S = B_nu(T_shell) along the
     ray. In a fixed-temperature SEDONA run with no radiative equilibrium,
     absorbed packets may simply be destroyed, with no thermal source from
     the gas -- which would make the solver brighter in saturated regions.
  B) ray count. Production used n_impact=150; the convergence sweep used 40
     and landed 0.7% lower, so the p-quadrature is not obviously converged.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu

T_EXP, R_CORE, R_OUT = 86400.0, 8.64e12, 2.592e13
T_CORE, V_D = 6000.0, 1.0e7

d = np.load("multiion_lines.npz")
RHO = float(d["rho"])
LINES = [
    (C / (lam * 1e-8), f, pf)
    for lam, f, pf in zip(d["lam"], d["f_lu"], d["popfrac_per_rho"])
]
nu = np.geomspace(7.50e14, 7.95e14, 400)
cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
lam = C / nu * 1e8
m = (lam > 3800) & (lam < 3955)
order = np.argsort(lam[m])


def band(t_shell, n_impact):
    l = emergent_luminosity(
        nu, LINES, lambda r: np.full_like(r, RHO),
        lambda r: np.full_like(r, t_shell),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D,
        n_impact=n_impact, n_z_per_doppler=8.0,
    )
    r = l / cont
    return np.trapezoid(r[m][order], lam[m][order]) / 155.0


print("A) shell emission (n_impact=40):", flush=True)
for t in [3000.0, 1000.0, 10.0]:
    print(f"   T_shell = {t:7g} K   band = {band(t, 40):.4f}", flush=True)

print("B) ray count (T_shell=3000):", flush=True)
for n in [20, 40, 80, 160, 320]:
    print(f"   n_impact = {n:4d}   band = {band(3000.0, n):.4f}", flush=True)
