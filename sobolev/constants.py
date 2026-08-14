"""Physical constants in CGS.

Kept local rather than pulled from astropy so that Phase 0 has no dependency
beyond numpy/scipy, and so every number entering a result is visible here.
"""

import numpy as np

C = 2.99792458e10  # speed of light, cm/s
E_ESU = 4.80320425e-10  # elementary charge, statcoulomb
M_E = 9.1093837015e-28  # electron mass, g
K_B = 1.380649e-16  # Boltzmann constant, erg/K
H = 6.62607015e-27  # Planck constant, erg s

# Classical line-absorption coefficient, pi e^2 / (m_e c), in cm^2 Hz.
# Derived, not hardcoded: this prefactor is shared by both the resolved and the
# Sobolev optical depth, so any error in it cancels in the comparison but corrupts
# every absolute tau. See docs/babystep_plan.md sections 3 and 21.
SIGMA_CLASSICAL = np.pi * E_ESU**2 / (M_E * C)  # ~0.02654 cm^2 Hz
