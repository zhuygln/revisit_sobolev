"""sobolev/xkn.py: the xkn-diff secular prescription (Ricigliano et al. 2024)
used as the P1 light curve's boundary condition (Paper IV Phase 10c)."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sobolev.constants import C                                     # noqa: E402
from sobolev.source import MSUN, DAY, thermalization                # noqa: E402
from sobolev import xkn                                             # noqa: E402


def test_photosphere_reaches_the_centre_at_eq28_and_sits_high_early():
    m = xkn.XknSecular()
    t2 = xkn.t_photosphere_centre(m.m, m.v_max_profile_c * C, m.kappa)
    assert abs(t2 / DAY - 119.0) < 2.0
    assert xkn.photosphere_x(t2 * 1.001, m.m, m.v_max_profile_c * C, m.kappa) == 0.0
    x1 = xkn.photosphere_x(t2 * 0.999, m.m, m.v_max_profile_c * C, m.kappa)
    assert 0.0 < x1 < 0.1
    x = [m.x_ph(t * DAY) for t in (0.5, 2.0, 8.0)]
    assert x[0] > x[1] > x[2] and 0.9 < x[1] < 0.92          # 0.908 at 2 d


def test_mass_outside_and_the_spherical_measure():
    assert xkn.mass_outside(0.0) == pytest.approx(1.0) and xkn.mass_outside(1.0) == pytest.approx(0.0)
    # against a numerical integral of x^2 (1-x^2)^3
    xs = np.linspace(0.5, 1.0, 200001)
    num = np.trapezoid(xs ** 2 * (1 - xs ** 2) ** 3, xs) / (16.0 / 315.0)
    assert abs(xkn.mass_outside(0.5) - num) < 1e-8


def test_thick_luminosity_modes_solve_the_ode_and_converge():
    m = xkn.XknSecular(n_modes=500)
    t = np.geomspace(0.5, 15, 60) * DAY
    L500 = m.L_diff(t); L2000 = xkn.XknSecular(n_modes=2000).L_diff(t)
    assert np.all(L2000 > 0) and np.max(np.abs(L500 / L2000 - 1.0)) < 0.03      # mode convergence (1.5 % at 0.5 d)
    # the n = 1 mode satisfies its ODE (finite differences on the internal grid)
    rho0 = m.m / (4.0 / 3.0 * np.pi * (m.v_max_diff_c * C * m.t0) ** 3)
    tau0 = 3.0 * m.kappa * rho0 * (m.v_max_diff_c * C * m.t0) ** 2 / C
    E0 = xkn.A_RAD * m.T0 ** 4
    tg = np.geomspace(m.t0, 15 * DAY, 4000)
    n = 1.0; sign = 1.0
    a = n ** 2 * np.pi ** 2 * tg ** 2 / (2 * m.t0 * tau0)
    src = (tg / m.t0) * m.heating(tg) * m.f_th0 * (tg / m.t0) ** (-m.beta)
    S = sign * rho0 * np.sqrt(2.0) / (n * np.pi * E0) * src
    phi = np.zeros(tg.size); phi[0] = 1.0
    for k in range(tg.size - 1):
        dt = tg[k + 1] - tg[k]; r = n ** 2 * np.pi ** 2 * 0.5 * (tg[k] + tg[k + 1]) / (m.t0 * tau0)
        phi[k + 1] = phi[k] * np.exp(-r * dt) + 0.5 * (S[k] + S[k + 1]) * (-np.expm1(-r * dt)) / r
    dphi = np.gradient(phi, tg)
    resid = dphi + (tg / (m.t0 * tau0)) * n ** 2 * np.pi ** 2 * phi - S
    mid = slice(200, -200)
    assert np.max(np.abs(resid[mid])) < 1e-3 * np.max(np.abs(S[mid]))


def test_thick_luminosity_is_bounded_by_the_heating_and_declines():
    m = xkn.XknSecular()
    t = np.geomspace(0.5, 15, 40) * DAY
    L = m.L_thick(t)
    assert np.all(np.diff(L) < 0)                                          # monotone decline after 0.5 d
    heat = m.m * m.heating(t) * m.f_th0 * (t / m.t0) ** (-m.beta)
    assert np.all(L < 3.0 * heat)                                          # the diffusion cannot exceed stored + current heating by much
    E_L = np.trapezoid(L, t); E_h = np.trapezoid(heat, t)
    assert 0.1 < E_L / E_h < 3.0


def test_temperatures_continue_at_the_photosphere_and_floor():
    m = xkn.XknSecular()
    t = 2.0 * DAY
    x_ph = m.x_ph(t); T_ph = float(m.T_ph(t)[0])
    assert 3500 < T_ph < 3900
    T = m.thin_temperatures(t, np.array([x_ph, 0.95, 0.999]))
    assert abs(T[0] - T_ph) < 1e-9 and T[1] < T_ph and T[2] == m.T_floor


def test_thin_thermalisation_reduces_to_barnes_at_x_zero():
    m = xkn.XknSecular()
    t = np.array([1.0, 3.0]) * DAY
    f0 = xkn.f_th_thin(t, 0.0, m.m_msun, m.v_max_profile_c)
    assert np.allclose(f0, thermalization(t, m.m_msun, m.v_max_profile_c))
    assert np.all(xkn.f_th_thin(t, 0.9, m.m_msun, m.v_max_profile_c) < f0)   # X = t/(1-x^2) larger -> less thermalised
    L_thin = m.L_thin(2 * DAY, np.linspace(0.75, 1.0, 25))
    assert 0 < L_thin < float(m.L_thick(2 * DAY)[0])
