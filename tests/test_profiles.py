"""Profile normalization: the cheapest bug to catch and the most expensive to miss.

A profile that integrates to 0.98 instead of 1.0 shifts every optical depth by 2%,
which is the same size as the effect Phase 3 is trying to measure.
"""

import numpy as np

from sobolev.constants import C
from sobolev.profiles import doppler_width_hz, gaussian, voigt


def test_gaussian_normalizes_to_unity():
    nu0 = C / 4000e-8  # 4000 Angstrom
    dnu_d = doppler_width_hz(nu0, 3.0e5)  # 3 km/s
    dnu = np.linspace(-10 * dnu_d, 10 * dnu_d, 20001)
    assert np.isclose(np.trapezoid(gaussian(dnu, dnu_d), dnu), 1.0, rtol=1e-6)


def test_gaussian_peaks_at_line_centre():
    nu0 = C / 4000e-8
    dnu_d = doppler_width_hz(nu0, 3.0e5)
    assert gaussian(0.0, dnu_d) > gaussian(dnu_d, dnu_d)


def test_voigt_normalizes_to_unity():
    # The Lorentzian wings fall off as 1/dnu^2, so normalization converges much
    # more slowly than the Gaussian: integrate far out and allow the truncated
    # tail (~ gamma / (pi^2 * cutoff)) in the tolerance.
    nu0 = C / 4000e-8
    dnu_d = doppler_width_hz(nu0, 3.0e5)
    gamma = 2.0 * dnu_d
    cutoff = 2000 * dnu_d
    dnu = np.linspace(-cutoff, cutoff, 400001)
    norm = np.trapezoid(voigt(dnu, dnu_d, gamma), dnu)
    assert np.isclose(norm, 1.0, rtol=1e-3)


def test_voigt_reduces_to_gaussian():
    nu0 = C / 4000e-8
    dnu_d = doppler_width_hz(nu0, 3.0e5)
    dnu = np.linspace(-5 * dnu_d, 5 * dnu_d, 1001)
    assert np.allclose(voigt(dnu, dnu_d, 0.0), gaussian(dnu, dnu_d), rtol=1e-12)


def test_voigt_wings_are_lorentzian():
    # Far from line centre the profile must follow gamma / (4 pi^2 dnu^2).
    nu0 = C / 4000e-8
    dnu_d = doppler_width_hz(nu0, 3.0e5)
    gamma = 0.1 * dnu_d
    far = 50 * dnu_d
    expected = gamma / (4 * np.pi**2 * far**2)
    assert np.isclose(float(voigt(far, dnu_d, gamma)), expected, rtol=1e-3)
