"""Source Model Gate (plan of 2026-09-02, Step 3): the heating-powered one-zone
source must pass these before any grid point is run."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sobolev import source as src
from sobolev.source import DAY, MSUN


def test_gate1_constant_deposition_is_analytic():
    """Qdot = const -> L = Qdot (1 - exp(-(t/tau_d)^2))."""
    t = src.time_grid(1e-2, 30.0, 3000)
    tau_d = 2.0 * DAY
    q = np.full_like(t, 1e41)
    L = src.arnett_luminosity(t, q, tau_d)
    exact = -1e41 * np.expm1(-(t / tau_d) ** 2)
    assert np.max(np.abs(L / exact - 1.0)) < 1e-4


def test_gate1b_matches_direct_quadrature():
    """The ODE integrator vs the Arnett integral evaluated by trapezoid."""
    t = src.time_grid(1e-2, 30.0, 30000)      # the trapezoid needs the fine grid, the ODE does not
    tau_d = 2.0 * DAY
    q = src.deposited_power(t, 0.01, 0.1)
    L = src.arnett_luminosity(t, q, tau_d)
    z = (t / tau_d) ** 2
    ok = t < 5 * tau_d                        # e^{z} stays representable
    integ = cumulative_trapezoid(t[ok] * np.exp(z[ok]) * q[ok], t[ok], initial=0.0)
    direct = 2.0 / tau_d ** 2 * np.exp(-z[ok]) * integ
    late = t[ok] > 0.1 * DAY                  # the trapezoid start-up error is early
    assert np.max(np.abs(L[ok][late] / direct[late] - 1.0)) < 1e-3
    # and the coarse production grid agrees with the fine one
    s3 = src.SourceModel(0.01, 0.1, n=3000); s30 = src.SourceModel(0.01, 0.1, n=30000)
    assert s3.luminosity(DAY) == pytest.approx(s30.luminosity(DAY), rel=1e-3)


def test_gate2_powerlaw_deposition_arnett_rule():
    """Qdot ~ t^-1.3: L peaks near tau_d and joins Qdot after a few tau_d."""
    t = src.time_grid(1e-2, 30.0, 3000)
    tau_d = 1.5 * DAY
    q = 1e41 * (t / DAY) ** -1.3
    L = src.arnett_luminosity(t, q, tau_d)
    tp = t[np.argmax(L)]
    assert 0.7 * tau_d <= tp <= 1.5 * tau_d
    i5 = np.argmin(np.abs(t - 5 * tau_d))
    assert abs(L[i5] / q[i5] - 1.0) < 0.10


def test_gate3_peak_time_scales_as_sqrt_kappa_m_over_v():
    """Doubling kappa M at fixed v moves the peak by sqrt(2) (+-5 %), when the
    deposition's own shape is held fixed."""
    t = src.time_grid(1e-2, 30.0, 6000)
    q = 1e41 * (t / DAY) ** -1.3
    tp = []
    for kappa_m in (1.0, 2.0):
        tau_d = src.diffusion_time(0.01 * kappa_m, 0.1, kappa=1.0)
        tp.append(t[np.argmax(src.arnett_luminosity(t, q, tau_d))])
    assert tp[1] / tp[0] == pytest.approx(np.sqrt(2.0), rel=0.05)
    # and the model's own tau_d has the stated scaling
    assert src.diffusion_time(0.02, 0.1) / src.diffusion_time(0.01, 0.1) == pytest.approx(np.sqrt(2.0))
    assert src.diffusion_time(0.01, 0.2) / src.diffusion_time(0.01, 0.1) == pytest.approx(1 / np.sqrt(2.0))


def test_gate4_energy_conservation():
    """int L dt = int Qdot(t') eta(t'/tau_d) dt' (closed form), hence <= int Qdot dt,
    with equality as kappa -> 0."""
    s = src.SourceModel(0.01, 0.1)
    E_rad = np.trapezoid(s.L, s.t)
    E_dep = np.trapezoid(s.qdot, s.t)
    E_pred = np.trapezoid(s.qdot * src.radiated_fraction(s.t / s.tau_d), s.t)
    assert E_rad == pytest.approx(E_pred, rel=1e-3)
    assert E_rad < E_dep
    # kappa -> 0: tau_d = 2 ms, below the grid's first point, so nothing is
    # lost to expansion and L tracks Qdot (the early plateau included)
    s0 = src.SourceModel(0.01, 0.1, kappa=1e-16)
    assert np.trapezoid(s0.L, s0.t) == pytest.approx(np.trapezoid(s0.qdot, s0.t), rel=0.02)
    late = s0.t > 10.0                        # past the arctan drop at t0 = 1.3 s
    assert np.max(np.abs(s0.L[late] / s0.qdot[late] - 1.0)) < 0.02
    assert src.radiated_fraction(1.0) == pytest.approx(0.758, abs=2e-3)
    assert src.radiated_fraction(0.0) == 0.0 and src.radiated_fraction(50.0) > 0.999


def test_gate5_central_model_plausible():
    s = src.SourceModel(0.01, 0.1)
    st1, st3 = s.state(1.0), s.state(3.0)
    assert 3e40 <= st1["L"] <= 5e41
    assert 4000 <= st1["T_eff"] <= 12000
    assert 2000 <= st3["T_eff"] <= 6000
    assert s.tau_d / DAY == pytest.approx(2.0, abs=0.3)
    assert st1["t_core"] == st1["T_gas"] == st1["T_eff"]
    # grey photosphere (default): kappa rho (R_out - R_ph) = 2/3 sits at 0.99 v_ej at
    # 1 d for the central model and recedes to the v_ej/2 floor by ~7 d
    # (`v_ph` is the fraction v_ph / v_ej)
    assert st1["r_core"] == pytest.approx(st1["v_ph"] * 0.1 * src.C * DAY)
    assert 0.95 <= st1["v_ph"] < 1.0
    assert st1["tau_grey"] > 10 and not st1["v_ph_floored"]
    st7 = s.state(7.0)
    assert 0.5 <= st7["v_ph"] < 0.6
    assert st1["r_out"] == pytest.approx(0.1 * src.C * DAY)
    # the plan's fixed convention is still available and is what the floor uses
    s_half = src.SourceModel(0.01, 0.1, v_ph_frac=0.5)
    assert s_half.state(1.0)["r_core"] == pytest.approx(0.5 * 0.1 * src.C * DAY)
    assert s_half.state(1.0)["T_eff"] > st1["T_eff"]        # smaller R_ph, same L
    # rho at 1 d from a uniform sphere out to v_ej, as observables.rho_1d_from_mass
    assert st1["rho"] == pytest.approx(0.01 * MSUN / (4 / 3 * np.pi * (0.1 * src.C * DAY) ** 3))
    # heating at 1 d: eps ~ 2e10 erg/g/s (Korobkin), f_th ~ 0.5
    assert src.heating_rate(DAY) == pytest.approx(2.0e10, rel=0.15)
    assert st1["Qdot"] == pytest.approx(2e41, rel=0.3)


def test_gate6_barnes_parameters():
    p, clamped = src.barnes_params(0.01, 0.1)
    assert p == pytest.approx((0.56, 0.17, 0.74)) and not clamped
    p, clamped = src.barnes_params(0.05, 0.3)
    assert p == pytest.approx((0.95, 0.15, 1.13)) and not clamped
    _, clamped = src.barnes_params(0.01, 0.05)
    assert clamped
    assert src.SourceModel(0.01, 0.05).fth_clamped and not src.SourceModel(0.01, 0.1).fth_clamped
    t = np.geomspace(0.1, 20, 200) * DAY
    f = src.thermalization(t, 0.01, 0.1)
    assert 0.4 < src.thermalization(DAY, 0.01, 0.1) < 0.6
    assert np.all(np.diff(f) < 0) and np.all((f > 0) & (f <= 0.72))
    # interpolation is exact on the table nodes and monotone between them
    a_mid = src.barnes_params(np.sqrt(1e-3 * 5e-3), 0.1)[0][0]
    assert 0.81 < a_mid < 2.01


def test_t_scale_perturbs_the_launch_temperature_only():
    """§4.43: `t_scale` changes t_core/T_eff at fixed L, R_ph, v_ph; T_gas only with t_scale_gas."""
    base = src.SourceModel(0.01, 0.1).state(2.0)
    hot = src.SourceModel(0.01, 0.1, t_scale=1.25).state(2.0)
    gas = src.SourceModel(0.01, 0.1, t_scale=1.25, t_scale_gas=True).state(2.0)
    for k in ("L", "R_ph", "r_core", "v_ph", "rho", "tau_grey", "v_ph_floored"):
        assert hot[k] == base[k] and gas[k] == base[k]
    assert hot["t_core"] == pytest.approx(1.25 * base["t_core"])
    assert hot["T_eff"] == pytest.approx(1.25 * base["T_eff"])
    assert hot["T_eff_grey"] == base["T_eff"]
    assert hot["T_gas"] == base["T_gas"] and gas["T_gas"] == pytest.approx(1.25 * base["T_gas"])
    assert (hot["t_scale"], hot["t_scale_gas"], gas["t_scale_gas"]) == (1.25, False, True)
    assert base["t_scale"] == 1.0


def test_v_ph_floored_flag_matches_the_clip():
    """The floor flag is set exactly when the grey photosphere would sit below v_ej/2."""
    m = src.SourceModel(0.003, 0.2)
    flags = [m.state(t)["v_ph_floored"] for t in (0.5, 1.0, 2.0, 3.0, 5.0, 7.0)]
    vph = [m.state(t)["v_ph"] for t in (0.5, 1.0, 2.0, 3.0, 5.0, 7.0)]
    assert flags == [v == src.V_PH_MIN for v in vph]
    assert any(flags) and not all(flags)
