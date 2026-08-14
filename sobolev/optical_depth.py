"""Resolved vs Sobolev optical depth for a 1-D homologous slab.

Both treatments live in one module on purpose: the entire Phase 0 result is the
comparison between them, and they must share the same prefactor, the same
oscillator strength, and the same population n_l(v). Splitting them across files
is how inconsistent inputs sneak into the comparison.

Geometry (babystep_plan.md section 3): v(x) = x / t, so a photon of lab frequency
nu sees the line at rest frequency nu0 come into resonance at the single point

    v_res = c (1 - nu0/nu).

The toy-model helpers (velocity grid, prescribed populations) also live here
because Phase 0 has exactly one model and does not warrant an ejecta package.
"""

import numpy as np

from .constants import C, SIGMA_CLASSICAL


# --------------------------------------------------------------------------
# Toy model: 1-D homologous slab with a prescribed lower-level population
# --------------------------------------------------------------------------


def homologous_velocity_grid(v_min, v_max, n_points):
    """Uniform velocity grid in cm/s. With v = x/t this is also the spatial grid."""
    return np.linspace(v_min, v_max, n_points)


def population_constant(v, n0):
    """Phase 0A: n_l(v) = n0. The Sobolev limit must be recovered exactly here."""
    return np.full_like(np.asarray(v, dtype=float), n0)


def population_tanh(v, n0, amplitude, v_res, v_scale):
    """Phase 0B: a controlled gradient that deliberately breaks Sobolev locality.

        n_l(v) = n0 [1 + A tanh((v - v_res) / v_scale)]

    The control parameter for the whole phase is epsilon = v_doppler / v_scale.
    """
    v = np.asarray(v, dtype=float)
    return n0 * (1.0 + amplitude * np.tanh((v - v_res) / v_scale))


# --------------------------------------------------------------------------
# The two optical depths under comparison
# --------------------------------------------------------------------------


def tau_sobolev(f_osc, n_l_at_resonance, lambda0_cm, t_seconds):
    """Sobolev optical depth (babystep_plan.md section 3).

        tau_S = (pi e^2 / m_e c) f n_l(v_res) lambda0 t

    Closed form, so it is implemented here. Note that it samples n_l at the
    resonance point only -- that locality is precisely what Phase 0B attacks.
    """
    return SIGMA_CLASSICAL * f_osc * n_l_at_resonance * lambda0_cm * t_seconds


def tau_exact(nu, nu0, f_osc, n_l_of_v, v_grid, t_seconds, v_doppler):
    """Resolved optical depth by direct integration along the slab.

        tau(nu) = \\int alpha_nu(x) dx,
        alpha_nu(x) = (pi e^2 / m_e c) f n_l(x) phi[nu'(x) - nu0],
        nu'(x) = nu (1 - v(x)/c).

    The grid must resolve v_doppler, not just v_scale: an unconverged grid
    looks exactly like a physical Sobolev deviation, which is the main way
    this experiment can fool you. tests/test_sobolev_limit.py checks
    convergence, not just agreement, for that reason.
    """
    from .profiles import doppler_width_hz, gaussian

    v_grid = np.asarray(v_grid, dtype=float)
    x = v_grid * t_seconds  # homologous: v = x/t, so dx = t dv
    nu_comoving = nu * (1.0 - v_grid / C)
    dnu_d = doppler_width_hz(nu0, v_doppler)
    alpha = SIGMA_CLASSICAL * f_osc * np.asarray(n_l_of_v, dtype=float) * gaussian(
        nu_comoving - nu0, dnu_d
    )
    return np.trapezoid(alpha, x)


def sobolev_error(tau_exact_value, tau_sobolev_value):
    """Relative Sobolev error E_Sob = |(tau_exact - tau_S) / tau_exact|."""
    return np.abs((tau_exact_value - tau_sobolev_value) / tau_exact_value)


def sweep_gradient_error(
    epsilons,
    nu0,
    f_osc,
    n0,
    amplitude,
    t_seconds,
    v_doppler,
    offset_in_scale=0.5,
    half_width_dopplers=15.0,
    points_per_doppler=200.0,
):
    """Phase 0B: E_Sob as a function of epsilon = v_D / v_scale.

    For each epsilon, a tanh population gradient of velocity scale
    v_scale = v_doppler / epsilon is placed so that the resonance point sits at
    u = offset_in_scale on the tanh curve, NOT at its centre. This offset is
    essential physics, not a tuning knob: with the resonance at the tanh centre
    the error cancels identically -- the Gaussian profile is even and the tanh
    is odd (first-order cancellation), and the centre is also the inflection
    point where n_l'' = 0 (second-order cancellation). A symmetric experiment
    would "show" that Sobolev never fails. Off-centre, the leading error is
    the curvature term, so E_Sob ~ epsilon^2 for small epsilon.

    The velocity grid spans +-half_width_dopplers around the resonance: the
    Gaussian confines all opacity contributions to a few Doppler widths, so the
    grid never needs to cover v_scale itself.

    Returns an array of E_Sob values, one per epsilon.
    """
    lambda0_cm = C / nu0
    errors = []
    for eps in np.asarray(epsilons, dtype=float):
        v_scale = v_doppler / eps
        n_points = int(2.0 * half_width_dopplers * points_per_doppler) + 1
        v_grid = np.linspace(
            -half_width_dopplers * v_doppler, half_width_dopplers * v_doppler, n_points
        )
        # Centre the tanh at -offset*v_scale so the resonance (v = 0 for
        # nu = nu0) sits at u = +offset_in_scale.
        v_centre = -offset_in_scale * v_scale
        n_l = population_tanh(v_grid, n0, amplitude, v_centre, v_scale)
        exact = tau_exact(nu0, nu0, f_osc, n_l, v_grid, t_seconds, v_doppler)
        n_l_res = float(population_tanh(0.0, n0, amplitude, v_centre, v_scale))
        analytic = tau_sobolev(f_osc, n_l_res, lambda0_cm, t_seconds)
        errors.append(sobolev_error(exact, analytic))
    return np.array(errors)
