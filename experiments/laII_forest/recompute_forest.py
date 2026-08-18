"""Recompute the La II forest three-leg table through the SHARED convention.

The numbers in the manuscript's forest table and in its thermal-emission (F8)
section were produced by compare.py, which predates sobolev.spectra and
carries two of the three normalization bugs this project has since found:

  * its red margin is 3952-3978 A, which straddles SEDONA's final spectrum
    bin -- that bin is partial, its flux collapses, and including it depresses
    the reference and inflates every band value (the bug of Section 4.9);
  * it band-averages by dividing by the NOMINAL band width rather than the
    span actually integrated.

The result was a forest table quoting SEDONA resolved = 0.3426 while the
separation section, computed through the shared helper, quoted 0.3435 for the
identical run -- the same measurement printed twice with different values two
pages apart.

This script recomputes all three legs, plus the F8 emission-off comparison,
through sobolev.spectra so that every forest number in the paper comes from
one pipeline. No new transport runs.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.spectra import band_average, band_ratio

T_EXP, R_CORE, R_OUT, T_CORE, T_SHELL = 86400.0, 8.64e12, 2.592e13, 6000.0, 3000.0
V_D = 1.0e7
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)

d = np.load(HERE / "forest_lines.npz")
N_ION = float(d["n_ion"])
LINES = [(C / (lam * 1e-8), f, p)
         for lam, f, p in zip(d["lam"], d["f_lu"], d["pop"])]


def solver(t_shell):
    """Deterministic resolved leg, band-averaged the shared way."""
    nu = np.geomspace(7.50e14, 7.95e14, 1600)
    lum = emergent_luminosity(
        nu, LINES, lambda r: np.full_like(r, N_ION),
        lambda r: np.full_like(r, t_shell),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D, n_impact=150,
    )
    ratio = lum / (4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE))
    return band_average(C / nu * 1e8, ratio, BAND)


res = band_ratio(HERE / "run_bb" / "spectrum_1.dat", BAND, MARGIN,
                 R_CORE, T_CORE)
exp = band_ratio(HERE / "run_exp" / "spectrum_1.dat", BAND, MARGIN,
                 R_CORE, T_CORE)
sol_warm = solver(T_SHELL)
sol_cold = solver(0.0)

print("=== La II forest, T=3000 K, day 1, v_D=100 km/s, 3800-3955 A ===")
print(f"  deterministic solver (S = B_nu(T_shell))  {sol_warm:.4f}"
      f"   vs SEDONA {(sol_warm - res) / res:+.1%}")
print(f"  deterministic solver (T_shell -> 0)       {sol_cold:.4f}"
      f"   vs SEDONA {(sol_cold - res) / res:+.1%}")
print(f"  SEDONA resolved                           {res:.4f}")
print(f"  SEDONA expansion opacity                  {exp:.4f}")
print(f"  Delta_expansion                           {(exp - res) / res:+.1%}")

out = dict(res=res, exp=exp, solver_warm=sol_warm, solver_cold=sol_cold,
           d_exp=(exp - res) / res,
           d_solver_warm=(sol_warm - res) / res,
           d_solver_cold=(sol_cold - res) / res)
(HERE / "forest_table.json").write_text(json.dumps(out, indent=1))
print("\nwrote forest_table.json")
