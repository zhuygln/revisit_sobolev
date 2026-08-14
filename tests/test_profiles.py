"""Profile normalization: the cheapest bug to catch and the most expensive to miss.

A profile that integrates to 0.98 instead of 1.0 shifts every optical depth by 2%,
which is the same size as the effect Phase 3 is trying to measure.
"""

import numpy as np

from sobolev.constants import C
from sobolev.profiles import doppler_width_hz, gaussian


def test_gaussian_normalizes_to_unity():
    nu0 = C / 4000e-8  # 4000 Angstrom
    dnu_d = doppler_width_hz(nu0, 3.0e5)  # 3 km/s
    dnu = np.linspace(-10 * dnu_d, 10 * dnu_d, 20001)
    assert np.isclose(np.trapezoid(gaussian(dnu, dnu_d), dnu), 1.0, rtol=1e-6)


def test_gaussian_peaks_at_line_centre():
    nu0 = C / 4000e-8
    dnu_d = doppler_width_hz(nu0, 3.0e5)
    assert gaussian(0.0, dnu_d) > gaussian(dnu_d, dnu_d)
