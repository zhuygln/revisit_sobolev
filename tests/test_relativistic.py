"""Frozen snapshot versus photon worldline.

The central fact pinned here: for a single line in homologous flow,

    frozen snapshot (t held fixed)   tau/tau_S = (1-beta)/gamma
    photon worldline (t advances)    tau/tau_S = 1/gamma

The first has an O(beta) term, the second does not. Finding F3 was read first
as a frame ambiguity and then as "the leading relativistic correction"; both
readings were wrong, and the offset is an artifact of integrating a frozen
snapshot. `test_ground_truth_*` establish this from first principles -- direct
integration of the resolved opacity along the path, with no Sobolev assumption
anywhere -- so the attribution rests on a calculation rather than on algebra.
"""

import numpy as np
import pytest

from sobolev.constants import C, SIGMA_CLASSICAL
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.optical_depth import tau_sobolev
from sobolev.sobolev_leg import tau_sobolev_frozen, tau_sobolev_relativistic

T_EXP = 864000.0
LAMBDA0 = 4000e-8
NU0 = C / LAMBDA0
F_OSC = 0.5
N0 = 2.0
V_D = 3.0e5


def const(v):
    return lambda r: np.full_like(np.asarray(r, dtype=float), v)


# ------------------------------------------------------- ground truth


def _ground_truth_tau(beta, worldline, v_d=1.0e6, half_widths=400, n=200001):
    """Integrate the resolved line opacity along the photon path directly.

    No Sobolev approximation is used: this is the quantity the candidate tau
    laws are trying to predict. `worldline=False` freezes the ejecta age.
    """
    dnu_d = NU0 * v_d / C
    g_res = 1.0 / np.sqrt(1.0 - beta**2)
    nu = NU0 / (g_res * (1.0 - beta))  # resonates exactly at beta
    r_res = beta * C * T_EXP
    half = half_widths * v_d * T_EXP
    s = np.linspace(-half, half, n)

    r = r_res + s
    t = T_EXP + s / C if worldline else T_EXP
    b = r / (C * t)
    lorentz = 1.0 / np.sqrt(1.0 - b**2)
    d = lorentz * (1.0 - b)
    phi = np.exp(-(((nu * d - NU0) / dnu_d) ** 2)) / (np.sqrt(np.pi) * dnu_d)
    tau = np.trapezoid(d * SIGMA_CLASSICAL * F_OSC * N0 * phi, s)
    return tau / tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)


@pytest.mark.parametrize("beta", [0.01, 0.05, 0.1, 0.2, 0.3])
def test_ground_truth_worldline_is_one_over_gamma(beta):
    lorentz = 1.0 / np.sqrt(1.0 - beta**2)
    assert np.isclose(_ground_truth_tau(beta, worldline=True), 1.0 / lorentz, rtol=1e-4)


@pytest.mark.parametrize("beta", [0.01, 0.05, 0.1, 0.2, 0.3])
def test_ground_truth_frozen_is_one_minus_beta_over_gamma(beta):
    lorentz = 1.0 / np.sqrt(1.0 - beta**2)
    assert np.isclose(
        _ground_truth_tau(beta, worldline=False), (1.0 - beta) / lorentz, rtol=1e-4
    )


def _ground_truth_tau_from_emission(
    beta, dilute, v_d=1.0e6, half_widths=400, n=200001
):
    """Same integration, but normalized against tau_S at the EMISSION epoch.

    `_ground_truth_tau` anchors its clock at the resonance, so it measures
    tau / tau_S(t_res) and is by construction blind to how the medium got
    there. A spectrum is labelled by the epoch the light left the source, not
    by the epoch each individual resonance was crossed, so this variant
    anchors at t0 instead.

    A photon leaving the origin at t0 is at z = c(t - t0), so the fluid it is
    passing has b = z/(ct) = 1 - t0/t, and it meets the material moving at
    beta at

        t_res = t0 / (1 - beta).

    Later than the frozen-geometry guess t0 (1 + beta), because the fluid
    element is itself moving outward while the photon chases it.

    dilute=False holds the density at N0 and so isolates the later crossing
    epoch on its own, which makes tau LARGER:

        tau / tau_S(t0) = 1 / [(1 - beta) gamma].

    dilute=True lets the homologous density fall as t^-3, which is the
    physical case and reverses the sign:

        tau / tau_S(t0) = (1 - beta)^2 / gamma.

    At beta = 0.3 those are 1.363 and 0.467. The dilution is therefore not a
    correction to the 1/gamma term, it is a factor of 2.9 that swamps it --
    see `test_dilution_dominates_the_lorentz_correction`.
    """
    t0 = T_EXP
    dnu_d = NU0 * v_d / C
    g_res = 1.0 / np.sqrt(1.0 - beta**2)
    nu = NU0 / (g_res * (1.0 - beta))  # resonates exactly at beta

    s_res = beta * C * t0 / (1.0 - beta)
    # db/ds = (1-b)^2/(c t0), so one Doppler width in velocity spans
    # v_d t0 / (1-beta)^2 of path length; keep the window in units of that.
    half = half_widths * v_d * t0 / (1.0 - beta) ** 2
    s = s_res + np.linspace(-half, half, n)

    t = t0 + s / C
    b = s / (C * t)  # radial ray from the origin, so r = s
    lorentz = 1.0 / np.sqrt(1.0 - b**2)
    d = lorentz * (1.0 - b)
    n_l = N0 * (t0 / t) ** 3 if dilute else N0
    phi = np.exp(-(((nu * d - NU0) / dnu_d) ** 2)) / (np.sqrt(np.pi) * dnu_d)
    tau = np.trapezoid(d * SIGMA_CLASSICAL * F_OSC * n_l * phi, s)
    return tau / tau_sobolev(F_OSC, N0, LAMBDA0, t0)


@pytest.mark.parametrize("beta", [0.01, 0.05, 0.1, 0.2, 0.3])
def test_from_emission_without_dilution_is_one_over_one_minus_beta_gamma(beta):
    """The crossing epoch alone. tau GROWS, because t_res > t0 and tau_S ~ t."""
    lorentz = 1.0 / np.sqrt(1.0 - beta**2)
    assert np.isclose(
        _ground_truth_tau_from_emission(beta, dilute=False),
        1.0 / ((1.0 - beta) * lorentz),
        rtol=1e-4,
    )


@pytest.mark.parametrize("beta", [0.01, 0.05, 0.1, 0.2, 0.3])
def test_from_emission_with_dilution_is_one_minus_beta_squared_over_gamma(beta):
    """The physical law against an observed spectrum's own epoch.

    tau_S ~ n t, n ~ t^-3, t_res = t0/(1-beta), so the medium contributes
    (1-beta)^3 * 1/(1-beta) = (1-beta)^2, and transport contributes 1/gamma.
    """
    lorentz = 1.0 / np.sqrt(1.0 - beta**2)
    assert np.isclose(
        _ground_truth_tau_from_emission(beta, dilute=True),
        (1.0 - beta) ** 2 / lorentz,
        rtol=1e-4,
    )


def test_dilution_dominates_the_lorentz_correction():
    """The effect this project calls "the relativistic correction" is the small one.

    At beta = 0.3 the 1/gamma deficit is 4.6% while the light-travel dilution
    deficit is 51%. Anything that models the first and not the second is
    wrong by an order of magnitude more than the term it kept.
    """
    beta = 0.3
    lorentz = 1.0 / np.sqrt(1.0 - beta**2)
    gamma_deficit = 1.0 - 1.0 / lorentz
    dilution_deficit = 1.0 - (1.0 - beta) ** 2
    assert np.isclose(gamma_deficit, 0.0461, atol=5e-4)
    assert np.isclose(dilution_deficit, 0.51, atol=5e-4)
    assert dilution_deficit > 10.0 * gamma_deficit


@pytest.mark.parametrize("beta", [0.05, 0.1, 0.2, 0.3])
def test_worldline_resonance_epoch_is_t0_over_one_minus_beta(beta):
    """Locate the crossing from the integrand itself, not from the algebra.

    Guards the sign error that motivated all of this: the frozen-geometry
    guess is t0 (1 + beta), which at beta = 0.3 is 1.30 t0 against the true
    1.43 t0 -- a difference far outside the tolerance here.
    """
    t0, v_d, n = T_EXP, 1.0e6, 200001
    dnu_d = NU0 * v_d / C
    g_res = 1.0 / np.sqrt(1.0 - beta**2)
    nu = NU0 / (g_res * (1.0 - beta))

    s_res = beta * C * t0 / (1.0 - beta)
    half = 400 * v_d * t0 / (1.0 - beta) ** 2
    s = s_res + np.linspace(-half, half, n)
    t = t0 + s / C
    b = s / (C * t)
    d = (1.0 - b) / np.sqrt(1.0 - b**2)
    phi = np.exp(-(((nu * d - NU0) / dnu_d) ** 2))

    # rtol is set by the grid spacing of the peak locator (dt/t0 ~ 1.5e-4),
    # not by the physics; it is still 250x tighter than the separation below.
    t_peak = t[np.argmax(phi)]
    assert np.isclose(t_peak / t0, 1.0 / (1.0 - beta), rtol=5e-4)

    # The frozen guess is excluded at every beta. The two laws separate by
    # 1/(1-b) - (1+b) = b^2/(1-b), which is only 0.26% at beta = 0.05, so a
    # fixed tolerance here would test the grid rather than the physics.
    separation = beta**2 / (1.0 - beta)
    assert abs(t_peak / t0 - (1.0 + beta)) > 0.5 * separation


@pytest.mark.parametrize("a", [0.0, 0.1, 0.3])
@pytest.mark.parametrize("y", [1.05, 1.2, 1.4])
def test_worldline_locus_is_linear_in_z(a, y):
    """z_res = (Z0/2)(y^2 - 1) + p^2/(2 Z0), one root, for every ray.

    Under worldline transport D = Z0/sqrt(Z0^2 + 2 Z0 z - p^2), so D = nu0/nu
    is LINEAR in z -- no discriminant and no branch to choose. The two-root
    quadratic that motivates Jeffery's CD/CP surfaces belongs to the frozen
    snapshot, not to the physical problem.
    """
    from scipy.optimize import brentq

    z0 = 0.0
    ct = C * T_EXP
    big_z = ct - z0
    p = a * ct

    def doppler_minus_x(z):
        return big_z / np.sqrt(big_z**2 + 2.0 * big_z * z - p**2) - 1.0 / y

    closed_form = 0.5 * big_z * (y**2 - 1.0) + p**2 / (2.0 * big_z)
    root = brentq(doppler_minus_x, 0.0, 10.0 * ct, xtol=1e-6, rtol=1e-14)
    assert np.isclose(root, closed_form, rtol=1e-12)


@pytest.mark.parametrize("a", [0.0, 0.1, 0.3])
def test_worldline_tau_is_impact_parameter_independent(a):
    """tau/tau_S(t_res) = 1/gamma exactly, for ANY impact parameter.

    D(z_res + Z0) = gamma Z0 and |dD/dz| = D^3/Z0, so the p-dependence cancels
    identically. This is why the plane-vs-CD/CP question does not arise in
    worldline transport.
    """
    y = 1.2
    ct = C * T_EXP
    big_z = ct
    p = a * ct

    z_res = 0.5 * big_z * (y**2 - 1.0) + p**2 / (2.0 * big_z)
    ct_res = z_res + big_z
    beta_res = np.sqrt(p**2 + z_res**2) / ct_res
    lorentz = 1.0 / np.sqrt(1.0 - beta_res**2)

    # tau = sigma f n_l Z0 / (nu0 D), tau_S(t_res) = sigma f n_l ct_res / nu0
    d_factor = big_z / np.sqrt(big_z**2 + 2.0 * big_z * z_res - p**2)
    ratio = (big_z / d_factor) / ct_res
    assert np.isclose(ratio, 1.0 / lorentz, rtol=1e-12)


# ------------------------------------------------------- the tau laws


def test_worldline_law_has_no_first_order_term():
    """1/gamma = 1 - beta^2/2: the deficit must scale as beta^2, not beta."""
    classical = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    b = np.array([0.01, 0.02, 0.04])
    deficit = 1.0 - tau_sobolev_relativistic(F_OSC, N0, LAMBDA0, T_EXP, b) / classical
    # quadrupling beta must roughly sixteen-fold the deficit
    assert np.isclose(deficit[2] / deficit[0], 16.0, rtol=0.05)


def test_frozen_law_has_a_first_order_term():
    classical = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    b = np.array([0.01, 0.04])
    deficit = 1.0 - tau_sobolev_frozen(F_OSC, N0, LAMBDA0, T_EXP, b) / classical
    assert np.isclose(deficit[1] / deficit[0], 4.0, rtol=0.05)


def test_frozen_law_equals_one_minus_beta_over_gamma_for_radial_rays():
    classical = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    b = np.array([0.01, 0.05, 0.1, 0.2, 0.3])
    lorentz = 1.0 / np.sqrt(1.0 - b**2)
    ratio = tau_sobolev_frozen(F_OSC, N0, LAMBDA0, T_EXP, b) / classical
    assert np.allclose(ratio, (1.0 - b) / lorentz, rtol=1e-9)


def test_both_laws_reduce_to_classical():
    classical = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    assert np.isclose(
        tau_sobolev_relativistic(F_OSC, N0, LAMBDA0, T_EXP, 0.0), classical, rtol=1e-12
    )
    assert np.isclose(
        tau_sobolev_frozen(F_OSC, N0, LAMBDA0, T_EXP, 0.0), classical, rtol=1e-12
    )


def test_worldline_absorbs_less_than_frozen():
    """1/gamma > (1-beta)/gamma, so the frozen treatment overstates how much
    the relativistic correction removes."""
    b = np.array([0.05, 0.1, 0.2, 0.3])
    rel = tau_sobolev_relativistic(F_OSC, N0, LAMBDA0, T_EXP, b)
    fro = tau_sobolev_frozen(F_OSC, N0, LAMBDA0, T_EXP, b)
    assert np.all(rel > fro)


# ------------------------------------------------------- solver modes

BETA_LO, BETA_HI = 5.0e-4, 0.35
R_CORE_W = BETA_LO * C * T_EXP
R_OUT_W = BETA_HI * C * T_EXP


def _trough(relativity, beta, n_impact=60):
    nu = np.array([NU0 / (1.0 - beta)])
    lum = emergent_luminosity(
        nu, [(NU0, F_OSC)], const(N0), const(10.0), T_EXP,
        R_CORE_W, R_OUT_W, 2.0e4, V_D, n_impact=n_impact, relativity=relativity,
    )
    cont = 4.0 * np.pi**2 * R_CORE_W**2 * planck_bnu(nu, 2.0e4)
    return float(lum[0] / cont[0])


def test_solver_modes_agree_at_low_beta():
    for mode in ("exact", "first"):
        assert np.isclose(
            _trough("worldline", 0.002), _trough(mode, 0.002), rtol=5e-3
        )


def test_solver_probes_actually_absorb():
    # Probes are labelled by the FROZEN resonance velocity nu = nu0/(1-beta).
    # Under worldline anchoring the ejecta age grows along the ray, so
    # b = z/(c t_exp + z - z0) saturates at ~0.26 for this shell and a
    # nominal beta = 0.3 has no resonance inside it at all. Stay below that:
    # the guard exists because a probe outside the shell absorbs nothing and
    # would make every comparison pass trivially.
    for beta in (0.002, 0.1, 0.2):
        assert _trough("worldline", beta) < 0.95


def test_solver_worldline_absorbs_more_than_frozen():
    """1/gamma > (1-beta)/gamma, so the worldline treatment gives the deeper
    trough; the frozen snapshot understates the absorption."""
    for beta in (0.05, 0.1, 0.2):
        assert _trough("worldline", beta) < _trough("exact", beta)


def test_invalid_mode_rejected():
    with pytest.raises(ValueError):
        _trough("newtonian", 0.01)
