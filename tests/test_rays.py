"""Shared impact-parameter quadrature, the closed-form resolved leg, and the
Poisson/Bernoulli identity.

These pin the infrastructure the referee response rests on: (i) both legs
on identical rays, with the legacy paths bit-identical; (ii) a resolved
reference that is exact where the Sobolev leg is exact and reproduces the
boundary law where it is not; (iii) the identity that makes the expansion
differential a statement about ensemble statistics rather than a bug.
"""

import numpy as np
import pytest

from sobolev.constants import C, H, K_B
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.optical_depth import stimulated_emission_factor, tau_sobolev
from sobolev.rays import RaySet, core_rays_midpoint
from sobolev.sobolev_leg import (
    crossing_depths,
    expansion_damp,
    resolved_attenuation,
    sobolev_attenuation,
)

T_EXP = 864000.0
R_CORE = 1.0e14
R_OUT = 5.0e14
LAMBDA0 = 4000e-8
NU0 = C / LAMBDA0
F_OSC = 0.5
N0 = 2.0
CT = C * T_EXP


def nu_for_plane(z_res, nu0=NU0):
    return nu0 / (1.0 - z_res / CT)


# ------------------------------------------------------------ the rays


def test_core_rays_midpoint_is_the_legacy_expression_bit_for_bit():
    for n in (1, 7, 200, 400):
        legacy = np.linspace(0.0, R_CORE, n, endpoint=False) + R_CORE / (2 * n)
        assert np.array_equal(core_rays_midpoint(R_CORE, n), legacy)


def test_midpoint_weights_integrate_the_disk_exactly():
    rs = RaySet.midpoint(R_CORE, R_OUT, 64, n_env=32)
    assert np.isclose(rs.w[rs.is_core].sum(), 0.5 * R_CORE**2, rtol=1e-12)
    assert np.isclose(rs.w[~rs.is_core].sum(), 0.5 * (R_OUT**2 - R_CORE**2), rtol=1e-12)
    assert rs.n_core == 64 and rs.p.size == 96


def test_gauss_legendre_integrates_a_cubic_exactly():
    rs = RaySet.gauss_legendre(R_CORE, R_OUT, 4)
    # int_0^R p * p^2 dp = R^4 / 4
    assert np.isclose(np.sum(rs.w * rs.p**2), R_CORE**4 / 4.0, rtol=1e-12)


def test_sobolev_leg_default_and_explicit_midpoint_rays_agree_exactly():
    nu = np.array([nu_for_plane(z) for z in (0.5e14, 1.5e14, 3.0e14)])
    a = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, n_p=200)
    rs = RaySet.midpoint(R_CORE, R_OUT, 200)
    b = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, rays=rs)
    assert (a < 0.999).all()
    assert np.allclose(a, b, rtol=1e-14, atol=0)


# ------------------------------------------------- the resolved leg (erf)


def test_resolved_leg_reduces_to_sobolev_as_width_vanishes():
    """With v_D -> 0 the erf bracket becomes the crossing test."""
    nu = np.array([nu_for_plane(z) for z in (0.5e14, 1.5e14, 3.0e14, 4.5e14)])
    sob = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    res = resolved_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, 1.0e3)
    assert (sob < 0.999).all()
    assert np.allclose(res, sob, rtol=1e-6)


def test_resolved_leg_half_counts_a_plane_sitting_on_the_edge():
    """The F12 boundary law: a resonance exactly at the outer edge has half its
    profile outside the shell, so the resolved depth is tau_S/2 where the
    Sobolev step counts all of it. Radial ray (tiny core) isolates it."""
    r_core = 1.0e-4 * R_OUT
    nu = np.array([nu_for_plane(R_OUT)])  # plane at z = r_out for p = 0
    tau_s = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    res = resolved_attenuation(
        nu, [(NU0, F_OSC)], r_core, R_OUT, T_EXP, N0, 1.0e6, n_p=1
    )[0]
    assert np.isclose(res, np.exp(-0.5 * tau_s), rtol=2e-3)


def test_resolved_leg_first_mode_matches_the_brute_force_solver():
    """The closed form against `emergent_luminosity(relativity='first')` on
    identical rays: a 3-line toy, T_shell -> 0, must agree to 2e-4."""
    lines = [(NU0, F_OSC), (NU0 * 1.004, 0.3 * F_OSC), (NU0 * 1.009, 0.8 * F_OSC)]
    v_d = 3.0e7  # 300 km/s keeps the brute-force z grid cheap
    rs = RaySet.midpoint(R_CORE, R_OUT, 40, n_env=0)
    nu = np.linspace(NU0 * 0.999, NU0 * 1.014, 60)
    const = lambda v: (lambda r: np.full_like(np.asarray(r, dtype=float), v))
    lum = emergent_luminosity(
        nu, lines, const(N0), const(0.0), T_EXP, R_CORE, R_OUT, 2.0e4, v_d,
        relativity="first", rays=rs,
    )
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, 2.0e4)
    brute = lum / cont
    closed = resolved_attenuation(
        nu, lines, R_CORE, R_OUT, T_EXP, N0, v_d, rays=rs, sweep="first"
    )
    assert (closed < 0.99).any()
    assert np.allclose(closed, brute, rtol=2e-4, atol=2e-4)


def test_sobolev_first_mode_is_the_delta_limit_of_the_first_order_integral():
    nu = np.array([nu_for_plane(z) for z in (0.5e14, 1.5e14, 3.0e14)])
    sob = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, relativity="first")
    res = resolved_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, 1.0e3, sweep="first")
    assert np.allclose(res, sob, rtol=1e-6)
    # and it differs from the classical plane by the nu0/nu factor, O(beta)
    classical = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    assert not np.allclose(sob, classical, rtol=1e-4)
    assert np.allclose(sob, classical, rtol=2e-2)


def test_empty_line_list_is_transparent_in_every_leg():
    nu = np.linspace(NU0 * 0.99, NU0 * 1.02, 11)
    assert np.all(sobolev_attenuation(nu, [], R_CORE, R_OUT, T_EXP, N0) == 1.0)
    assert np.all(resolved_attenuation(nu, [], R_CORE, R_OUT, T_EXP, N0, 1e6) == 1.0)
    S, E, _, _ = crossing_depths(nu, [], R_CORE, R_OUT, T_EXP, N0)
    assert not S.any() and not E.any()


# ------------------------------------------------ the Poisson identity


def test_expansion_over_sobolev_is_the_transmission_weighted_mean_of_exp_deficit():
    """F_exp / F_Sob = E_w[e^{S-E}], w ~ p e^{-S}: exact, by construction, and
    the statement that turns the expansion differential into ensemble
    statistics. Pinned to 1e-13 on a 5-line forest."""
    rng = np.random.default_rng(3)
    lines = [(NU0 * (1.0 + 0.012 * rng.random()), F_OSC * rng.uniform(0.1, 3.0))
             for _ in range(5)]
    nu = np.linspace(NU0 * 0.999, NU0 * 1.016, 80)
    S, E, p, w = crossing_depths(nu, lines, R_CORE, R_OUT, T_EXP, N0)
    f_sob = np.sum(w[:, None] * np.exp(-S), axis=0) / w.sum()
    f_exp = np.sum(w[:, None] * np.exp(-E), axis=0) / w.sum()
    # the legs themselves agree with these sums
    assert np.allclose(f_sob, sobolev_attenuation(nu, lines, R_CORE, R_OUT, T_EXP, N0), rtol=1e-13)
    assert np.allclose(
        f_exp, sobolev_attenuation(nu, lines, R_CORE, R_OUT, T_EXP, N0, damp=expansion_damp),
        rtol=1e-13,
    )
    wt = w[:, None] * np.exp(-S)
    ew = np.sum(wt * np.exp(S - E), axis=0) / np.sum(wt, axis=0)
    assert np.allclose(f_exp / f_sob, ew, rtol=1e-13)
    # and E is what the expansion leg exponentiates: the expected count is
    # preserved, the transmission is not
    assert (E <= S + 1e-15).all()
    assert (f_exp >= f_sob - 1e-15).all()


def test_poisson_gap_closes_as_lines_weaken():
    lines = [(NU0 * (1.0 + 0.003 * k), F_OSC) for k in range(4)]
    nu = np.linspace(NU0 * 0.999, NU0 * 1.016, 40)
    gaps = []
    for scale in (1.0, 0.1, 0.01):
        f_sob = sobolev_attenuation(nu, lines, R_CORE, R_OUT, T_EXP, N0 * scale)
        f_exp = sobolev_attenuation(nu, lines, R_CORE, R_OUT, T_EXP, N0 * scale, damp=expansion_damp)
        gaps.append(np.max(f_exp / f_sob - 1.0))
    assert gaps[0] > 0.05 and gaps[2] < 1e-3
    # second order in tau: a 10x weaker forest has a ~100x smaller gap
    assert np.isclose(gaps[1] / gaps[2], 100.0, rtol=0.25)


# -------------------------------------------------- stimulated emission


def test_stimulated_emission_factor_values():
    assert np.isclose(stimulated_emission_factor(C / 3800e-8, 3000.0), 1.0 - 3.4e-6, atol=2e-7)
    assert np.isclose(stimulated_emission_factor(C / 9100e-8, 3000.0), 1.0 - 5.1e-3, atol=2e-4)
    nu = C / 7000e-8
    assert np.isclose(
        stimulated_emission_factor(nu, 3000.0), 1.0 - np.exp(-H * nu / (K_B * 3000.0)), rtol=1e-12
    )
