"""Line profile functions.

Phase 0 uses Gaussian only; Voigt arrives in Week 1 (babystep_plan.md section 18).
Every profile here must satisfy the normalization

    \\int phi(nu) d(nu) = 1

which is enforced by tests/test_profiles.py.
"""

import numpy as np

from .constants import C


def doppler_width_hz(nu0, v_doppler_cms):
    """Doppler width in Hz for rest frequency `nu0` and thermal velocity `v_doppler_cms`."""
    return nu0 * v_doppler_cms / C


def gaussian(dnu, dnu_doppler):
    """Normalized Gaussian profile evaluated at frequency offset `dnu`.

        phi(dnu) = exp(-(dnu/dnu_D)^2) / (sqrt(pi) dnu_D)
    """
    dnu = np.asarray(dnu, dtype=float)
    return np.exp(-((dnu / dnu_doppler) ** 2)) / (np.sqrt(np.pi) * dnu_doppler)


def voigt(dnu, dnu_doppler, gamma):
    """Normalized Voigt profile. Deferred to Week 1 -- Phase 0 is Gaussian-only."""
    raise NotImplementedError("Voigt profile: Week 1, not Phase 0")
