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
