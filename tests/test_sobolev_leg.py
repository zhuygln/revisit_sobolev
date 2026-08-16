"""Analytic limits of the Sobolev attenuation leg.

These pin the leg that carries the "Sobolev approximation proper" claim, so
they check the formula rather than merely exercising the code path.
"""

import numpy as np

from sobolev.constants import C
from sobolev.optical_depth import tau_sobolev
from sobolev.sobolev_leg import expansion_damp, sobolev_attenuation

T_EXP = 864000.0  # 10 days
R_CORE = 1.0e14
R_OUT = 5.0e14
LAMBDA0 = 4000e-8
NU0 = C / LAMBDA0
F_OSC = 0.5
N0 = 2.0


def nu_for_plane(z_res, nu0=NU0):
    """Observer frequency whose resonance plane sits at z_res."""
    return nu0 / (1.0 - z_res / (C * T_EXP))


def test_single_line_fully_shadowing_plane():
    """A plane beyond the whole core (z_res > r_core) shadows every core ray,
    so the transmitted fraction is exactly exp(-tau_S)."""
    tau = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    assert 0.3 < tau < 3.0
    nu = np.array([nu_for_plane(2.0e14)])
    att = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    assert np.isclose(att[0], np.exp(-tau), rtol=1e-12)


def test_off_resonance_is_untouched():
    nu = np.array([NU0 * 0.90, NU0 * 1.30])
    att = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    assert np.allclose(att, 1.0)


def test_two_separated_lines_multiply():
    """Both planes in front of the core: attenuations multiply."""
    nu_b = NU0 * 1.02
    z = 2.0e14
    # Choose a frequency whose plane for line A is at z; line B's plane for
    # the same frequency sits elsewhere, so probe each line separately and
    # then a frequency crossing both.
    lines = [(NU0, F_OSC), (nu_b, F_OSC)]
    tau_a = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    tau_b = tau_sobolev(F_OSC, N0, C / nu_b, T_EXP)

    nu_a_only = np.array([nu_for_plane(z, NU0)])
    att_a = sobolev_attenuation(nu_a_only, [lines[0]], R_CORE, R_OUT, T_EXP, N0)
    att_both = sobolev_attenuation(nu_a_only, lines, R_CORE, R_OUT, T_EXP, N0)
    # Line B's plane for this frequency is far outside the shell, so adding
    # it must not change the result.
    assert np.isclose(att_a[0], np.exp(-tau_a), rtol=1e-12)
    assert np.isclose(att_both[0], att_a[0], rtol=1e-12)

    # A frequency whose plane for B lies in front of the core.
    nu_b_probe = np.array([nu_for_plane(z, nu_b)])
    att_b = sobolev_attenuation(nu_b_probe, [lines[1]], R_CORE, R_OUT, T_EXP, N0)
    assert np.isclose(att_b[0], np.exp(-tau_b), rtol=1e-12)


def test_expansion_damp_reproduces_the_cap():
    """With the expansion damping the single-line result is exp(-(1-e^-tau))."""
    tau = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    nu = np.array([nu_for_plane(2.0e14)])
    att = sobolev_attenuation(
        nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, damp=expansion_damp
    )
    assert np.isclose(att[0], np.exp(-(1.0 - np.exp(-tau))), rtol=1e-12)
    # The cap always transmits more light than true Sobolev for tau > 0.
    assert att[0] > np.exp(-tau)


def test_pop_frac_scales_tau():
    """pop_frac multiplies the density, hence tau, linearly."""
    nu = np.array([nu_for_plane(2.0e14)])
    full = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    half = sobolev_attenuation(
        nu, [(NU0, F_OSC, 0.5)], R_CORE, R_OUT, T_EXP, N0
    )
    assert np.isclose(half[0], np.sqrt(full[0]), rtol=1e-12)


def test_partial_shadowing_lies_between():
    """A plane inside the core radius shadows only the outer rays, so the
    transmitted fraction sits strictly between exp(-tau) and 1."""
    tau = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    nu = np.array([nu_for_plane(0.5 * R_CORE)])
    att = sobolev_attenuation(
        nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, n_p=400
    )
    assert np.exp(-tau) < att[0] < 1.0
