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

    THIS IS THE SESSION 1 EXERCISE -- implement it yourself. The recipe:

      1. map v_grid -> x = v * t (so dx = t dv);
      2. comoving frequency nu' = nu (1 - v/c);
      3. profile offset dnu = nu' - nu0, width from profiles.doppler_width_hz;
      4. alpha = SIGMA_CLASSICAL * f_osc * n_l_of_v * profiles.gaussian(dnu, dnu_D);
      5. integrate over x with np.trapezoid.

    Watch the resolution: the grid must resolve v_doppler, not just v_scale.
    An unconverged grid looks exactly like a physical Sobolev deviation, which
    is the main way this experiment can fool you.
    """
    raise NotImplementedError("Phase 0A, Session 1 -- see docstring for the recipe")


def sobolev_error(tau_exact_value, tau_sobolev_value):
    """Relative Sobolev error E_Sob = |(tau_exact - tau_S) / tau_exact|."""
    return np.abs((tau_exact_value - tau_sobolev_value) / tau_exact_value)
