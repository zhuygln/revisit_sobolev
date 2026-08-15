"""Precompute the Python-solver leg for the multi-ion blend (slow: 2529 lines)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity

T_EXP = 86400.0
R_CORE = 8.64e12
R_OUT = 2.592e13
T_CORE = 6000.0
T_SHELL = 3000.0
V_D = 1.0e7

d = np.load("multiion_lines.npz")
RHO = float(d["rho"])
LINES = [
    (C / (lam * 1e-8), f, pf)
    for lam, f, pf in zip(d["lam"], d["f_lu"], d["popfrac_per_rho"])
]
print(f"{len(LINES)} lines")

nu_py = np.geomspace(7.50e14, 7.95e14, 1600)
l_py = emergent_luminosity(
    nu_py, LINES, lambda r: np.full_like(r, RHO),
    lambda r: np.full_like(r, T_SHELL),
    T_EXP, R_CORE, R_OUT, T_CORE, V_D, n_impact=150,
)
np.savez("solver_result.npz", nu=nu_py, lum=l_py)
print("wrote solver_result.npz")
