"""Paper IV WP1: EjectaState, the profiles, Gate 1's checks, and the
local-shell single-zone reduction."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sobolev.constants import C                     # noqa: E402
from sobolev.source import MSUN, SourceModel        # noqa: E402
from sobolev import ejecta as ej                    # noqa: E402

DAY = 86400.0


def _state(n_shell=16, profile="power", t=2.0 * DAY, m=0.01 * MSUN):
    if profile == "power":
        v_edges, rho = ej.power_law_profile(m, t, 0.15 * C, 0.35 * C, -3.0, n_shell)
    else:
        v_edges, rho = ej.xkn_profile(m, t, 0.115 * C, n_shell, v_min_frac=0.05)
    n = rho.size
    X = {"Ce": np.full(n, 1e-3), "Nd": np.full(n, 1.5e-3), "bulk": np.full(n, 1.0 - 2.5e-3)}
    f_ion = {"Ce II": np.full(n, 0.7), "Ce III": np.full(n, 0.3), "Nd II": np.ones(n)}
    return ej.EjectaState(t=t, v_edges=v_edges, rho=rho, T_gas=np.full(n, 3200.0),
                          T_rad=np.full(n, 3200.0), X=X, f_ion=f_ion, Y_e=0.29, meta={"name": "toy"})


@pytest.mark.parametrize("profile", ["power", "xkn"])
def test_profiles_integrate_to_the_target_mass(profile):
    m = 0.0264 * MSUN
    st = _state(profile=profile, m=m)
    assert abs(st.mass() - m) / m < 1e-12          # exact by construction
    out = st.check(m_target=m)
    assert out["mass_rel"] < 1e-12 and out["X_sum_max_dev"] < 1e-12 and out["f_ion_max_dev"] < 1e-12


def test_power_law_shells_follow_v_to_the_gamma():
    st = _state(n_shell=64)
    v = st.v_mid
    slope = np.polyfit(np.log(v), np.log(st.rho), 1)[0]
    assert abs(slope + 3.0) < 0.02


def test_gate_one_catches_a_broken_state():
    st = _state()
    st.X["bulk"] = st.X["bulk"] * 0.9
    with pytest.raises(ValueError):
        st.check()
    st = _state()
    st.f_ion["Ce III"] = st.f_ion["Ce III"] * 0.5
    with pytest.raises(ValueError):
        st.check()
    st = _state()
    with pytest.raises(ValueError):
        st.check(m_target=st.mass() * 1.01)


def test_xkn_vmax_reproduces_the_rms_velocity():
    v_max = ej.xkn_vmax_from_vrms(0.06 * C)
    assert abs(v_max / C - 0.1149) < 2e-4
    v_edges, rho = ej.xkn_profile(0.0264 * MSUN, 1.0 * DAY, v_max, 400)
    st = ej.EjectaState(t=DAY, v_edges=v_edges, rho=rho, T_gas=np.ones(400), T_rad=np.ones(400),
                        X={"bulk": np.ones(400)})
    m = st.shell_mass()
    v_rms = np.sqrt(np.sum(m * st.v_mid ** 2) / m.sum())
    assert abs(v_rms / C - 0.06) < 5e-4


def test_photospheric_shell_and_local_zone():
    st = _state(n_shell=32)
    kappa = 1.0
    tg = st.tau_grey(kappa)
    assert np.all(np.diff(tg) <= 0) and tg[-1] > 0
    s, ok = st.photospheric_shell(kappa)
    assert ok and tg[s] >= 2 / 3 and (s == st.n_shell - 1 or tg[s + 1] < 2 / 3)
    z = st.local_zone(s)
    assert z["rho"] == st.rho[s] and z["T_gas"] == st.T_gas[s] and z["r_core"] == st.r_edges[s]
    assert z["r_out"] == st.r_edges[-1] and z["core_law"] == "local_shell" and z["shell"] == s
    assert z["X"]["Ce"] == 1e-3 and z["f_ion"]["Ce III"] == 0.3
    # the local zone is the shell's own state, not an average
    assert z["rho"] != np.average(st.rho, weights=st.shell_mass())
    # n_ion of Ce II in that shell
    n = st.n_ion("Ce", "II", s, 140.1)
    assert np.isclose(n, st.rho[s] * 1e-3 * 0.7 / (140.1 * ej.M_U))
    # a thin ejecta never reaches 2/3
    thin = _state(m=1e-6 * MSUN)
    s2, ok2 = thin.photospheric_shell(kappa)
    assert not ok2 and s2 == 0


def test_grey_temperature_run_and_source_consistency():
    st = _state(n_shell=32)
    src = SourceModel(0.01, 0.35)
    L = src.luminosity(st.t)
    T, s, t_eff = ej.grey_temperature(st, 1.0, L)
    assert T.shape == (32,) and np.all(np.diff(T) <= 1e-9)        # cooler outward
    assert T[-1] < t_eff < T[s] * 1.001 * (1.5 * (st.tau_grey(1.0)[s] + 2 / 3)) ** 0.25


def test_json_round_trip(tmp_path):
    st = _state()
    st.n_e = np.full(st.n_shell, 1e8)
    st.to_json(tmp_path / "s.json")
    back = ej.EjectaState.from_json(tmp_path / "s.json")
    assert back.mass() == st.mass() and back.Y_e == 0.29 and back.meta["name"] == "toy"
    assert np.array_equal(back.rho, st.rho) and np.array_equal(back.f_ion["Ce II"], st.f_ion["Ce II"])
    assert np.array_equal(back.n_e, st.n_e)
    back.check(m_target=st.mass())


def test_plot_writes_a_figure(tmp_path):
    st = _state()
    out = st.plot(tmp_path / "g.png", kappa=1.0, atomic_mass={"Ce": 140.1, "Nd": 144.2})
    assert Path(out).stat().st_size > 1000
