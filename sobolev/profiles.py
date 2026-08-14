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
    """Normalized Voigt profile evaluated at frequency offset `dnu`.

    Gaussian core of Doppler width `dnu_doppler` convolved with a Lorentzian of
    FWHM `gamma` (both in Hz), via the Faddeeva function:

        phi(dnu) = Re[w(z)] / (sqrt(pi) dnu_D),
        z = (dnu + i gamma/2) / dnu_D

    gamma is the full damping width (natural + collisional); for a purely
    radiative line gamma = A_ul / (2 pi). In the gamma -> 0 limit this reduces
    to `gaussian` exactly, which tests/test_profiles.py checks.
    """
    from scipy.special import wofz

    dnu = np.asarray(dnu, dtype=float)
    z = (dnu + 0.5j * gamma) / dnu_doppler
    return wofz(z).real / (np.sqrt(np.pi) * dnu_doppler)
