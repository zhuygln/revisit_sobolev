"""The band-ratio convention, pinned against synthetic spectra.

These exist because three separate spurious results came from getting this
normalization wrong. A null spectrum must return exactly 1, and the helper
must be insensitive to the two things that broke before: an overall scale
factor, and the Planck slope across the band.
"""

import numpy as np
import pytest

from sobolev.constants import C
from sobolev.formal_transfer import planck_bnu
from sobolev.spectra import band_average, band_ratio

R_CORE, T_CORE = 8.64e12, 6000.0
BAND = (3800.0, 3955.0)
MARGIN = (3960.0, 3975.0)


def _write(tmp_path, ratio_fn, scale=1.0):
    lam = np.linspace(3700.0, 3990.0, 4000)
    nu = C / (lam * 1e-8)
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    lum = scale * cont * ratio_fn(lam)
    p = tmp_path / "spectrum_1.dat"
    np.savetxt(p, np.column_stack([nu, lum, np.ones_like(nu)]))
    return p


def test_null_spectrum_returns_unity(tmp_path):
    """No absorption anywhere -> exactly 1. This is the control that caught
    the breadth-sweep bug (it returned 1.075)."""
    p = _write(tmp_path, lambda lam: np.ones_like(lam))
    assert np.isclose(band_ratio(p, BAND, MARGIN, R_CORE, T_CORE), 1.0, rtol=1e-6)


def test_insensitive_to_overall_scale(tmp_path):
    """SEDONA's lightbulb inflates absolute L_nu by 1/f_window (~48x here);
    the helper must divide that out."""
    trough = lambda lam: np.where((lam > 3850) & (lam < 3900), 0.3, 1.0)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = band_ratio(_write(tmp_path / "a", trough, scale=1.0),
                   BAND, MARGIN, R_CORE, T_CORE)
    b = band_ratio(_write(tmp_path / "b", trough, scale=48.5),
                   BAND, MARGIN, R_CORE, T_CORE)
    assert np.isclose(a, b, rtol=1e-9)


def test_removes_the_planck_slope(tmp_path):
    """A flat-in-ratio spectrum must give 1 even though L_nu itself varies
    across the band -- normalizing by raw luminosity instead of the continuum
    ratio is exactly what produced the spurious Sobolev floor."""
    p = _write(tmp_path, lambda lam: np.ones_like(lam))
    got = band_ratio(p, BAND, MARGIN, R_CORE, T_CORE)
    # compare against the wrong convention to show it would NOT give 1
    nu, lum = np.loadtxt(p)[:, 0], np.loadtxt(p)[:, 1]
    lam = C / nu * 1e8
    red = (lam > MARGIN[0]) & (lam < MARGIN[1])
    m = (lam > BAND[0]) & (lam < BAND[1])
    o = np.argsort(lam[m])
    wrong = np.trapezoid((lum / np.mean(lum[red]))[m][o], lam[m][o]) / (
        BAND[1] - BAND[0]
    )
    assert np.isclose(got, 1.0, rtol=1e-6)
    assert abs(wrong - 1.0) > 0.01  # the old convention is off by >1%


def test_known_absorption_recovered(tmp_path):
    """A top-hat trough of known depth over a known fraction of the band."""
    depth, lo, hi = 0.25, 3850.0, 3900.0
    p = _write(tmp_path, lambda lam: np.where((lam > lo) & (lam < hi), depth, 1.0))
    frac = (hi - lo) / (BAND[1] - BAND[0])
    expected = 1.0 * (1 - frac) + depth * frac
    assert np.isclose(band_ratio(p, BAND, MARGIN, R_CORE, T_CORE), expected, rtol=2e-3)


def test_empty_margin_rejected(tmp_path):
    p = _write(tmp_path, lambda lam: np.ones_like(lam))
    with pytest.raises(ValueError):
        band_ratio(p, BAND, (5000.0, 5010.0), R_CORE, T_CORE)


def test_band_average_matches_band_ratio_convention(tmp_path):
    """The analytic-leg helper must use the same averaging as the SEDONA one,
    or cross-code comparisons reintroduce the bug."""
    depth, lo, hi = 0.4, 3820.0, 3910.0
    p = _write(tmp_path, lambda lam: np.where((lam > lo) & (lam < hi), depth, 1.0))
    sed = band_ratio(p, BAND, MARGIN, R_CORE, T_CORE)
    lam = np.linspace(3700.0, 3990.0, 6000)
    ana = band_average(lam, np.where((lam > lo) & (lam < hi), depth, 1.0), BAND)
    assert np.isclose(sed, ana, rtol=3e-3)
