"""Synthetic photometry: the magnitudes the closure error gets quoted in.

Every §10 number is a difference of magnitudes, so the absolute scale has to be
right before any difference means anything. The limits pinned here are the ones
that would silently corrupt a result rather than crash it: an AB zero point off
by a constant, a Planck normalization off by 4pi, a spectrum that integrates to
the wrong luminosity, or a self-difference that is not exactly zero.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paper2/phase1"))
from forest_mc import ForestAtom, run_mc

from sobolev.constants import C, H
from sobolev.formal_transfer import planck_bnu
from sobolev import photometry as phot

SIGMA_SB = 5.670374419e-5     # Stefan-Boltzmann, erg cm^-2 s^-1 K^-4
T_EXP = 86400.0


def test_planck_luminosity_is_stefan_boltzmann():
    """Over the full spectrum, 4 pi^2 r^2 int B_nu dnu = 4 pi r^2 sigma T^4.

    This is the only check that the `4 pi^2` of sobolev/spectra.py:64 -- not
    4 pi, not pi -- is the convention actually implemented.
    """
    r, T = 1e13, 6000.0
    L = phot.planck_luminosity(1e10, 1e18, r, T)
    assert L == pytest.approx(4 * np.pi * r**2 * SIGMA_SB * T**4, rel=1e-4)


def test_ab_zero_point():
    """A flat f_nu = 3631 Jy source is AB 0.000 in every band, colours 0."""
    edges = phot.nu_edges(1000.0, 30000.0, 400)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    l_nu = np.full(nu_c.size, phot.F_AB * 4 * np.pi * phot.D_40MPC**2)
    mags = phot.magnitudes(nu_c, l_nu)
    assert set(mags) == set(phot.BANDS_PHOT)
    for b, m in mags.items():
        assert m == pytest.approx(0.0, abs=1e-9), b
    for c, v in phot.colors(mags).items():
        assert v == pytest.approx(0.0, abs=1e-9), c


def test_delta_mag_of_itself_is_exactly_zero():
    """The distance and zero point must cancel to the last bit, or every
    reported dm carries a floor set by the arbitrary distance."""
    edges = phot.nu_edges(1000.0, 30000.0, 200)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    l_nu = planck_bnu(nu_c, 5000.0) * 1e30
    m = phot.magnitudes(nu_c, l_nu)
    assert all(v == 0.0 for v in phot.delta_mag(m, m).values())


def test_magnitude_is_distance_dependent_but_colour_is_not():
    """Doubling the distance dims everything by the same 1.505 mag."""
    edges = phot.nu_edges(1000.0, 30000.0, 200)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    l_nu = planck_bnu(nu_c, 5000.0) * 1e30
    near = phot.magnitudes(nu_c, l_nu, distance_cm=phot.D_40MPC)
    far = phot.magnitudes(nu_c, l_nu, distance_cm=2 * phot.D_40MPC)
    for b in near:
        assert far[b] - near[b] == pytest.approx(2.5 * np.log10(4.0), abs=1e-9)
    for c, v in phot.colors(far).items():
        assert v == pytest.approx(phot.colors(near)[c], abs=1e-9)


def test_hotter_blackbody_is_bluer():
    """A sanity anchor with a known sign: g-r must decrease with temperature."""
    edges = phot.nu_edges(1000.0, 30000.0, 400)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    gr = [phot.colors(phot.magnitudes(nu_c, planck_bnu(nu_c, T) * 1e30))["g-r"]
          for T in (3000.0, 5000.0, 9000.0)]
    assert gr[0] > gr[1] > gr[2]


def _transparent_atom():
    """A shell that is transparent to the launch window.

    The single line sits at 100 A, blueward of every launched packet. Comoving
    frequency only ever decreases along a leg, so no packet can redshift up to
    it -- the shell is empty in practice while `op_nu` stays non-empty, which
    `run_mc`'s resonance search requires.
    """
    return ForestAtom(nu0=np.array([C / 100e-8]), f_osc=np.array([1e-3]),
                      n_lower=np.array([1e6]), n_upper=np.zeros(1),
                      A=np.array([1e8]), lower=np.zeros(1, int),
                      upper=np.ones(1, int), t_exp=T_EXP, tau_min=1e-3,
                      stim=False)


def test_transparent_shell_returns_the_core_blackbody():
    """With no opacity the emergent spectrum IS the core, so its colours must
    match the analytic Planck colours and L_bol must equal L_core."""
    T_core, r_core = 6000.0, 8.64e12
    lo, hi = phot.nu_edges(1000.0, 30000.0, 1)
    l_core = phot.planck_luminosity(lo, hi, r_core, T_core)
    res = run_mc(_transparent_atom(), r_core, 3 * r_core, T_EXP, lo, hi,
                 400000, "sobolev_branch", seed=1, t_core=T_core)
    assert res["n_interactions"] == 0

    edges = phot.nu_edges(1000.0, 30000.0, 200)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    l_nu = phot.emergent_lnu(res, edges, l_core)

    # integrating the MC spectrum recovers the escaping luminosity
    assert float(np.sum(l_nu * np.diff(edges))) == pytest.approx(
        phot.bolometric(res, l_core), rel=1e-9)
    # ... which, with nothing absorbing, is the whole core within the window
    assert phot.bolometric(res, l_core) == pytest.approx(l_core, rel=0.02)

    mc = phot.colors(phot.magnitudes(nu_c, l_nu))
    exact = phot.colors(phot.magnitudes(nu_c, planck_bnu(nu_c, T_core)))
    for c in exact:
        assert mc[c] == pytest.approx(exact[c], abs=0.02), c


def test_bol_delta_mag_sign():
    """Negative dm means too bright -- the convention every table uses."""
    assert phot.bol_delta_mag(2.0, 1.0) == pytest.approx(-2.5 * np.log10(2.0))
    assert phot.bol_delta_mag(1.0, 2.0) > 0
    assert np.isnan(phot.bol_delta_mag(0.0, 1.0))


def test_band_outside_the_grid_is_nan_not_wrong():
    """A band the launch window does not cover must refuse to answer rather
    than integrate over whatever happens to be at the edge."""
    edges = phot.nu_edges(4000.0, 5000.0, 50)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    l_nu = np.ones(nu_c.size)
    assert np.isnan(phot.ab_magnitude(nu_c, l_nu, phot.BANDS_PHOT["K"]))
    assert np.isfinite(phot.ab_magnitude(nu_c, l_nu, (4200.0, 4800.0)))


def test_absolute_spectrum_agrees_with_band_ratio():
    """The absolute L_nu divided by the analytic core continuum must reproduce
    `band_ratio`, which is computed from the packet list by a different route.

    This is the end-to-end check on the normalization: `emergent_lnu` scales by
    L_core/E_inj, `band_ratio` divides escaped by launched weights, and they can
    only agree if the `4 pi^2 r^2 B_nu` convention and the E_inj bookkeeping are
    both right. Resolution-limited, so it is asserted on bands many bins wide.
    """
    from forest_mc import band_ratio
    from sobolev.optical_depth import tau_sobolev

    T_core, r_core = 6000.0, 8.64e12
    lam = np.geomspace(3000.0, 9000.0, 60) * 1e-8          # a modest forest
    f = np.full(lam.size, 0.05)
    n_l = np.array([1.5 / tau_sobolev(fi, 1.0, l, T_EXP) for fi, l in zip(f, lam)])
    atom = ForestAtom(nu0=C / lam, f_osc=f, n_lower=n_l, n_upper=np.zeros(lam.size),
                      A=np.full(lam.size, 1e8), lower=np.zeros(lam.size, int),
                      upper=np.arange(1, lam.size + 1), t_exp=T_EXP,
                      tau_min=1e-3, stim=False)

    lo, hi = (float(x) for x in phot.nu_edges(1000.0, 30000.0, 1))
    l_core = phot.planck_luminosity(lo, hi, r_core, T_core)
    res = run_mc(atom, r_core, 3 * r_core, T_EXP, lo, hi, 300000,
                 "sobolev_branch", seed=3, t_core=T_core)

    edges = phot.nu_edges(1000.0, 30000.0, 400)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    l_nu = phot.emergent_lnu(res, edges, l_core)
    cont = 4.0 * np.pi**2 * r_core**2 * planck_bnu(nu_c, T_core)

    for lam_lo, lam_hi in ((4000.0, 5500.0), (5500.0, 6900.0), (11000.0, 14000.0)):
        nu_lo, nu_hi = C / (lam_hi * 1e-8), C / (lam_lo * 1e-8)
        m = (nu_c >= nu_lo) & (nu_c <= nu_hi)
        mine = np.trapezoid(l_nu[m], nu_c[m]) / np.trapezoid(cont[m], nu_c[m])
        ref = band_ratio(res, nu_lo, nu_hi, weight="energy")[0]
        assert mine == pytest.approx(ref, rel=0.02), (lam_lo, lam_hi)
