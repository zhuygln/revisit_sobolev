"""Formal transfer solver: analytic limits it must reproduce before use.

The validation sequence of babystep_plan.md section 13: bare continuum,
then a single line in the Sobolev limit, then two lines.
"""

import numpy as np

from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.optical_depth import tau_sobolev

# t_exp is chosen large so that the test resonance planes sit at
# v_bulk/c < 1%: the solver (observer frame) and the Sobolev formula
# (comoving frame) differ at O(v_bulk/c) -- see the module docstring.
T_EXP = 864000.0  # 10 days
R_CORE = 1.0e14
R_OUT = 5.0e14
T_CORE = 1.0e4
V_D = 3.0e5  # 3 km/s
LAMBDA0 = 4000e-8
NU0 = C / LAMBDA0
F_OSC = 0.5


def const_n(value):
    return lambda r: np.full_like(np.asarray(r, dtype=float), value)


def test_bare_core_is_blackbody_sphere():
    """No opacity: L_nu = 4 pi^2 r_core^2 B_nu(T_core) exactly."""
    nu = np.linspace(0.9 * NU0, 1.1 * NU0, 5)
    lum = emergent_luminosity(
        nu, [(NU0, F_OSC)], const_n(0.0), const_n(100.0), T_EXP,
        R_CORE, R_OUT, T_CORE, V_D,
    )
    expected = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    assert np.allclose(lum, expected, rtol=1e-3)


def test_single_line_trough_is_sobolev():
    """Cold shell (S ~ 0), constant n_l: the absorption trough must reach
    exp(-tau_S) for frequencies whose resonance plane lies in front of the
    whole core -- the Sobolev limit of the resolved calculation."""
    n0 = 2.0
    tau_s = tau_sobolev(F_OSC, n0, LAMBDA0, T_EXP)
    assert 0.5 < tau_s < 3.0  # keep the test in a regime where the trough is visible

    # Blueshifted frequencies resonate at z > 0; pick one whose plane sits
    # between the core front and the outer edge for every core ray.
    z_res = 2.0e14
    nu_probe = NU0 / (1.0 - z_res / (C * T_EXP))
    nu = np.array([0.97 * NU0, nu_probe, 1.05 * NU0])

    lum = emergent_luminosity(
        nu, [(NU0, F_OSC)], const_n(n0), const_n(10.0), T_EXP,
        R_CORE, R_OUT, T_CORE, V_D,
    )
    lum0 = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    depth = lum[1] / lum0[1]
    assert np.isclose(depth, np.exp(-tau_s), rtol=2e-2)
    # Far from resonance the continuum is untouched.
    assert np.isclose(lum[0] / lum0[0], 1.0, rtol=1e-3)
    assert np.isclose(lum[2] / lum0[2], 1.0, rtol=1e-3)


def test_two_lines_are_independent_when_separated():
    """Two well-separated lines: each trough matches its own exp(-tau_S)."""
    n0 = 2.0
    nu_b = NU0 * 1.02  # second line 2% blueward -- far beyond v_D/c ~ 1e-5
    z_res = 2.0e14
    shift = 1.0 / (1.0 - z_res / (C * T_EXP))
    nu = np.array([NU0 * shift, nu_b * shift])

    lum = emergent_luminosity(
        nu, [(NU0, F_OSC), (nu_b, 0.5 * F_OSC)], const_n(n0), const_n(10.0),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D,
    )
    lum0 = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    tau_a = tau_sobolev(F_OSC, n0, LAMBDA0, T_EXP)
    tau_b = tau_sobolev(0.5 * F_OSC, n0, C / nu_b, T_EXP)
    assert 0.3 < tau_a < 3.0 and 0.3 < tau_b < 3.0  # non-trivial troughs
    assert np.isclose(lum[0] / lum0[0], np.exp(-tau_a), rtol=2e-2)
    assert np.isclose(lum[1] / lum0[1], np.exp(-tau_b), rtol=2e-2)


def test_two_lines_blend_multiplicatively():
    """Overlapping troughs: a frequency whose path crosses BOTH resonance
    planes in front of the core is attenuated by exp(-(tau_1 + tau_2))."""
    n0 = 2.0
    dv = 1.5e8  # 1500 km/s separation -- well inside the shell's 4000 km/s span
    nu_b = NU0 / (1.0 - dv / C)
    # Probe so that line A resonates at z = 3e14 and line B at ~1.7e14: both
    # planes lie between the core front (1e14) and the outer edge.
    z_a = 3.0e14
    nu = np.array([NU0 / (1.0 - z_a / (C * T_EXP))])
    lum = emergent_luminosity(
        nu, [(NU0, F_OSC), (nu_b, F_OSC)], const_n(n0), const_n(10.0),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D,
    )
    lum0 = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    tau_a = tau_sobolev(F_OSC, n0, LAMBDA0, T_EXP)
    tau_b = tau_sobolev(F_OSC, n0, C / nu_b, T_EXP)
    assert np.isclose(lum[0] / lum0[0], np.exp(-(tau_a + tau_b)), rtol=3e-2)


def test_voigt_and_cutoff_options():
    """4-tuple lines select Voigt; a tiny gamma with a generous cutoff must
    reproduce the Gaussian trough, and a cutoff cannot deepen absorption."""
    n0 = 2.0
    z_res = 2.0e14
    nu = np.array([NU0 / (1.0 - z_res / (C * T_EXP))])
    base = emergent_luminosity(
        nu, [(NU0, F_OSC)], const_n(n0), const_n(10.0), T_EXP,
        R_CORE, R_OUT, T_CORE, V_D,
    )
    tiny_gamma = emergent_luminosity(
        nu, [(NU0, F_OSC, 1.0, 1e-6 * NU0 * V_D / C)], const_n(n0),
        const_n(10.0), T_EXP, R_CORE, R_OUT, T_CORE, V_D, cutoff_widths=30.0,
    )
    assert np.isclose(tiny_gamma[0], base[0], rtol=1e-3)
    truncated = emergent_luminosity(
        nu, [(NU0, F_OSC)], const_n(n0), const_n(10.0), T_EXP,
        R_CORE, R_OUT, T_CORE, V_D, cutoff_widths=5.0,
    )
    assert truncated[0] >= base[0] * (1 - 1e-9)  # less opacity, never more


def test_warm_shell_fills_in_the_trough():
    """With S = B_nu(T_shell) > 0 the trough must be shallower than pure
    absorption -- emission partially refills it."""
    n0 = 2.0
    z_res = 2.0e14
    nu = np.array([NU0 / (1.0 - z_res / (C * T_EXP))])
    cold = emergent_luminosity(
        nu, [(NU0, F_OSC)], const_n(n0), const_n(10.0), T_EXP,
        R_CORE, R_OUT, T_CORE, V_D,
    )
    warm = emergent_luminosity(
        nu, [(NU0, F_OSC)], const_n(n0), const_n(8000.0), T_EXP,
        R_CORE, R_OUT, T_CORE, V_D,
    )
    assert warm[0] > cold[0]
