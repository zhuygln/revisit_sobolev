"""Frame-treatment tests: the beta -> 0 limit, and the size of the O(beta)
term that separates the solver's historical mode from exact relativity.

These exist because Finding F3 recorded an O(v_bulk/c) discrepancy as an
unexplained systematic. The derivation in sobolev_leg.tau_sobolev_relativistic
attributes it to a missing Doppler factor on the opacity; these tests pin both
the limit and the leading correction so the attribution cannot silently rot.
"""

import numpy as np

from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.optical_depth import tau_sobolev
from sobolev.sobolev_leg import tau_sobolev_relativistic

T_EXP = 864000.0
R_CORE = 1.0e14
R_OUT = 5.0e14
T_CORE = 2.0e4
LAMBDA0 = 4000e-8
NU0 = C / LAMBDA0
F_OSC = 0.5
N0 = 2.0
V_D = 3.0e5


def const(v):
    return lambda r: np.full_like(np.asarray(r, dtype=float), v)


# ---------------------------------------------------------------- tau formula


def test_relativistic_tau_reduces_to_classical():
    """beta -> 0 must recover the textbook Sobolev depth exactly."""
    classical = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    rel = tau_sobolev_relativistic(F_OSC, N0, LAMBDA0, T_EXP, 0.0)
    assert np.isclose(rel, classical, rtol=1e-12)


def test_relativistic_tau_leading_order_is_one_doppler_factor():
    """tau_rel ~ tau_S (1 - beta) at small beta. The opacity transformation
    contributes D^2 but the sweep rate cancels one factor, so the net leading
    behaviour is a SINGLE Doppler factor -- which is why the solver's
    first-order mode, despite omitting D entirely, agrees at O(beta)."""
    classical = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    for beta in (1e-4, 1e-3, 1e-2):
        rel = tau_sobolev_relativistic(F_OSC, N0, LAMBDA0, T_EXP, beta)
        assert np.isclose(rel / classical, 1.0 - beta, rtol=5 * beta**2 + 1e-9)


def test_relativistic_tau_departs_from_one_minus_beta_at_second_order():
    """At beta = 0.3 the exact result must differ measurably from (1-beta),
    otherwise the gamma terms are not actually being carried."""
    classical = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    rel = tau_sobolev_relativistic(F_OSC, N0, LAMBDA0, T_EXP, 0.3) / classical
    assert abs(rel - 0.7) > 0.01 * 0.7


def test_relativistic_tau_decreases_with_beta():
    betas = np.array([0.0, 0.01, 0.05, 0.1, 0.2, 0.3])
    taus = tau_sobolev_relativistic(F_OSC, N0, LAMBDA0, T_EXP, betas)
    assert np.all(np.diff(taus) < 0)


# ---------------------------------------------------------------- solver modes


# A shell wide enough in velocity to place resonances from beta ~ 1e-3 out to
# beta ~ 0.3. The earlier fixed geometry only reached beta = 0.019, so a
# "high beta" probe fell outside the shell and absorbed nothing.
BETA_LO, BETA_HI = 5.0e-4, 0.35
R_CORE_W = BETA_LO * C * T_EXP
R_OUT_W = BETA_HI * C * T_EXP


def _trough(relativity, beta, n_impact=60):
    """Transmitted fraction at the frequency whose resonance sits at beta."""
    z_res = beta * C * T_EXP
    nu = np.array([NU0 / (1.0 - beta)])
    lum = emergent_luminosity(
        nu, [(NU0, F_OSC)], const(N0), const(10.0), T_EXP,
        R_CORE_W, R_OUT_W, T_CORE, V_D, n_impact=n_impact,
        relativity=relativity,
    )
    cont = 4.0 * np.pi**2 * R_CORE_W**2 * planck_bnu(nu, T_CORE)
    return float(lum[0] / cont[0])


def test_modes_agree_at_low_beta():
    """At beta ~ 2e-3 the modes must agree to well under a percent, since
    they differ only at O(beta^2)."""
    assert np.isclose(_trough("first", 0.002), _trough("exact", 0.002), rtol=3e-3)


def test_modes_absorb_something_in_the_wide_shell():
    """Guard against the failure that fooled the first draft of this file:
    a resonance outside the shell absorbs nothing and every comparison
    trivially passes."""
    for beta in (0.002, 0.1, 0.3):
        assert _trough("exact", beta) < 0.95


def test_modes_separate_at_second_order():
    """The two modes agree at O(beta) but must visibly separate by beta = 0.3,
    where O(beta^2) is several percent."""
    lo = abs(_trough("exact", 0.002) - _trough("first", 0.002))
    hi = abs(_trough("exact", 0.3) - _trough("first", 0.3))
    assert hi > 10 * lo


def test_invalid_mode_rejected():
    import pytest

    with pytest.raises(ValueError):
        _trough("newtonian", 0.01)
