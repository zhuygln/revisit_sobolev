"""Real transmission curves (`data/filters/`) and `photometry.Passband`."""
import hashlib
import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sobolev import photometry as ph
from sobolev.formal_transfer import planck_bnu

EDGES = ph.nu_edges(1000.0, 30000.0, 200)
NU_C = np.sqrt(EDGES[1:] * EDGES[:-1])


@pytest.fixture(scope="module")
def pb():
    return ph.load_passbands()


def test_curves_load_and_are_sane(pb):
    assert set(pb) == set(ph.BANDS_PHOT)
    for b, p in pb.items():
        assert p.T.min() >= 0.0 and 0.5 < p.T.max() <= 1.0
        assert 1000.0 < p.lam_lo < p.lam_hi < 30000.0


def test_lam_eff_matches_readme(pb):
    readme = (ph.FILTER_DIR / "README.md").read_text()
    for b, p in pb.items():
        m = re.search(rf"`{ph.FILTER_FILES[b]}\.dat`.*?\| (\d+) \| `([0-9a-f]{{16}})`", readme)
        assert m, b
        assert p.lam_eff == pytest.approx(float(m.group(1)), rel=0.01)
        sha = hashlib.sha256((ph.FILTER_DIR / f"{ph.FILTER_FILES[b]}.dat").read_bytes()).hexdigest()
        assert sha[:16] == m.group(2), f"{b}: curve changed since the README was written"


def test_flat_ab_spectrum_is_zero_mag(pb):
    lnu = np.full(NU_C.size, ph.F_AB * 4.0 * np.pi * ph.D_40MPC**2)
    for b, p in pb.items():
        assert ph.ab_magnitude(NU_C, lnu, p, edges=EDGES) == pytest.approx(0.0, abs=1e-6)


def test_blackbody_histogram_matches_fine_integration(pb):
    """A 6000 K blackbody through the 200-bin pipeline vs a 20001-point grid."""
    T, R = 6000.0, 1e14
    fine = np.geomspace(EDGES[0], EDGES[-1], 20001)
    lf = 4.0 * np.pi**2 * R**2 * planck_bnu(fine, T)
    lh = np.array([np.trapezoid(lf[(fine >= a) & (fine <= c)], fine[(fine >= a) & (fine <= c)]) / (c - a)
                   for a, c in zip(EDGES[:-1], EDGES[1:])])
    f = lf / (4.0 * np.pi * ph.D_40MPC**2)
    for b in ("g", "r", "i", "z", "J", "H", "K"):
        p = pb[b]
        Tf = np.interp(ph.C / fine * 1e8, p.lam, p.T, left=0.0, right=0.0)
        m_fine = -2.5 * np.log10(np.trapezoid(f * Tf / fine, fine) / np.trapezoid(Tf / fine, fine) / ph.F_AB)
        assert ph.ab_magnitude(NU_C, lh, p, edges=EDGES) == pytest.approx(m_fine, abs=0.01)


def test_delta_mag_invariant_under_scaling(pb):
    rng = np.random.default_rng(3)
    lnu = 1e38 * np.exp(rng.normal(0, 0.3, NU_C.size))
    m1 = ph.magnitudes(NU_C, lnu, pb, edges=EDGES)
    m2 = ph.magnitudes(NU_C, 10.0 * lnu, pb, edges=EDGES)
    for b in pb:
        assert m1[b] - m2[b] == pytest.approx(2.5, abs=1e-9)


def test_passband_outside_grid_is_nan(pb):
    edges = ph.nu_edges(4000.0, 9000.0, 50)
    lnu = np.ones(50)
    assert np.isnan(ph.ab_magnitude(np.sqrt(edges[1:] * edges[:-1]), lnu, pb["g"], edges=edges))
    assert np.isfinite(ph.ab_magnitude(np.sqrt(edges[1:] * edges[:-1]), lnu, pb["r"], edges=edges))
    with pytest.raises(ValueError):
        ph.band_flux_nu(NU_C, lnu, pb["r"])          # a Passband needs edges


def test_tophat_path_unchanged(pb):
    """The tuple path never sees `edges`; results are those of F41."""
    rng = np.random.default_rng(5)
    lnu = 1e38 * np.exp(rng.normal(0, 0.3, NU_C.size))
    a = ph.magnitudes(NU_C, lnu)
    b = ph.magnitudes(NU_C, lnu, edges=EDGES)
    assert a == b
    assert ph.band_flux_nu(NU_C, lnu, ph.BANDS_PHOT["g"]) == ph.band_flux_nu(NU_C, lnu, ph.BANDS_PHOT["g"], edges=EDGES)
