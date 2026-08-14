"""The Phase 0A gate.

babystep_plan.md section 3.1: "If this calculation cannot reproduce the analytic
Sobolev result, do not move on."

These tests fail until tau_exact is implemented. That is intentional -- they are
the specification for Session 1, not regression cover for finished code.
"""

import numpy as np
import pytest

from sobolev.constants import C
from sobolev.optical_depth import (
    homologous_velocity_grid,
    population_constant,
    sobolev_error,
    tau_exact,
    tau_sobolev,
)

LAMBDA0_CM = 4000e-8
NU0 = C / LAMBDA0_CM
F_OSC = 0.1
N0 = 1.0e6  # cm^-3
T_SECONDS = 86400.0  # 1 day
V_DOPPLER = 3.0e5  # 3 km/s


def test_constant_population_recovers_sobolev():
    """Phase 0A: with n_l constant, the resolved integral must converge to tau_S."""
    v_grid = homologous_velocity_grid(-3.0e8, 3.0e8, 200001)
    n_l = population_constant(v_grid, N0)

    exact = tau_exact(NU0, NU0, F_OSC, n_l, v_grid, T_SECONDS, V_DOPPLER)
    analytic = tau_sobolev(F_OSC, N0, LAMBDA0_CM, T_SECONDS)

    assert sobolev_error(exact, analytic) < 1e-3


def test_error_shrinks_under_grid_refinement():
    """Convergence, not just agreement: a coarse grid must not accidentally pass."""
    analytic = tau_sobolev(F_OSC, N0, LAMBDA0_CM, T_SECONDS)
    errors = []
    for n_points in (20001, 200001):
        v_grid = homologous_velocity_grid(-3.0e8, 3.0e8, n_points)
        n_l = population_constant(v_grid, N0)
        exact = tau_exact(NU0, NU0, F_OSC, n_l, v_grid, T_SECONDS, V_DOPPLER)
        errors.append(sobolev_error(exact, analytic))
    assert errors[1] < errors[0]


@pytest.mark.xfail(reason="Phase 0B result -- write this once the gradient sweep runs")
def test_gradient_breaks_sobolev():
    """Phase 0B: E_Sob must grow as epsilon = v_D / v_scale approaches unity."""
    raise NotImplementedError
