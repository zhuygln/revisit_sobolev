"""Is the solver's 6.5% excess brightness in the dense blend a resolution
artifact? Phase 0A's lesson in reverse: with 2529 overlapping lines the
opacity along a ray is structured on the Doppler scale everywhere, not just
in isolated resonance regions, so the trapezoidal tau integral needs more
z-points than the isolated-line case.

Cheap settings (few rays, coarse nu grid) -- we want the TREND with
n_z_per_doppler, not an absolute number.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu

T_EXP, R_CORE, R_OUT = 86400.0, 8.64e12, 2.592e13
T_CORE, T_SHELL, V_D = 6000.0, 3000.0, 1.0e7

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

for n_z in [4.0, 8.0, 16.0, 32.0, 64.0]:
    l = emergent_luminosity(
        nu, LINES, lambda r: np.full_like(r, RHO),
        lambda r: np.full_like(r, T_SHELL),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D,
        n_impact=40, n_z_per_doppler=n_z,
    )
    r = l / cont
    band = np.trapezoid(r[m][order], lam[m][order]) / 155.0
    print(f"n_z_per_doppler = {n_z:5g}   band = {band:.4f}", flush=True)
