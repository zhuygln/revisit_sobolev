"""Paper IV Phase 10b: the time-slab contract of run_mc (pause at t_stop,
resume, absolute injection energy, escape times, the event cap) and the
sobolev.timeslab helpers. Toy atoms from tests/test_zoned_run_mc.py, worldline
transport, energy packets."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1", ROOT / "paper2/phase0/three_level_atom"):
    sys.path.insert(0, str(p))

from sobolev.constants import C, H                          # noqa: E402
from sobolev import timeslab as ts                          # noqa: E402
from forest_mc import run_mc                                # noqa: E402
import test_zoned_run_mc as TZ                              # noqa: E402

R_CORE, R_OUT, T_EXP = TZ.R_CORE, TZ.R_OUT, TZ.T_EXP
CT = C * T_EXP
KW = dict(packets="energy", t_core=6000.0, launch_weight="energy", relativity="worldline")


def _zone3(scale=(1.0, 1.0, 1.0)):
    fa, a = TZ.forest()
    edges = np.array([R_CORE, 1.8 * R_CORE, 3.0 * R_CORE, R_OUT])
    return TZ.zoned(a, 3, edges, scale=scale), edges


def _chain(z3, lo, hi, first_res, dt, mode, max_slabs=40, seed=7):
    """Resume `first_res`'s paused packets through slabs of duration dt until
    none is paused; returns the per-slab results."""
    out = [first_res]
    t = T_EXP
    for k in range(max_slabs):
        pop = ts.carried_population(out[-1], z3.r_edges * ((t + dt) / T_EXP))
        if pop.n == 0:
            break
        # the atom is frozen (same tables); the geometry is the slab's own epoch
        r = run_mc(z3, z3.r_edges[0], z3.r_edges[-1], T_EXP, lo, hi, 0, mode, seed=seed + k, resume=pop.as_dict(),
                   t_stop=T_EXP + (k + 2) * dt, **KW)
        out.append(r); t += dt
    return out


@pytest.mark.parametrize("mode", ("sobolev_dmacro", "expansion_dmacro", "sobolev_thermal"))
def test_slab_split_reproduces_the_single_run_and_the_accounting_closes(mode):
    z3, edges = _zone3()
    lo, hi = TZ.tfm.pump_band()
    n = 40000
    single = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, n, mode, seed=5, **KW)
    dt = 400.0                                  # light crossing of the zone is ~1150 s
    first = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, n, mode, seed=5, t_stop=T_EXP + dt, **KW)
    assert first["n_paused"] > 0.2 * n          # the split is real
    slabs = _chain(z3, lo, hi, first, dt, mode)
    for r in slabs:
        assert abs(r["accounting"]["identity_residual"]) < 1e-12
    for a, b in zip(slabs[:-1], slabs[1:]):
        assert np.isclose(b["accounting"]["E_carried_in"], a["accounting"]["E_carried_out"], rtol=1e-12)
    e_in = single["accounting"]["E_inj"]
    e_esc_single = single["accounting"]["E_esc"]
    e_esc_slabs = sum(r["accounting"]["E_esc"] for r in slabs)
    p = e_esc_single / e_in
    assert abs(e_esc_slabs - e_esc_single) / e_in < 4 * np.sqrt(p * (1 - p) / n) + 0.01
    # the escape spectrum: chi^2 on the shared frequency histogram
    nu_s = single["nu_out"]; nu_c = np.concatenate([r["nu_out"] for r in slabs])
    b_edges = np.quantile(nu_s, np.linspace(0, 1, 21))
    h1, _ = np.histogram(nu_s, b_edges); h2, _ = np.histogram(nu_c, b_edges)
    chi2 = np.sum((h1 - h2) ** 2 / np.maximum(h1 + h2, 1))
    assert chi2 < 20 + 4 * np.sqrt(2 * 20)


def test_escape_time_along_radial_rays_through_an_opacity_free_zone():
    z3, edges = _zone3(scale=(1e-12, 1e-12, 1e-12))     # no line has tau > tau_min
    assert z3.n_opacity == 0
    lo, hi = TZ.tfm.pump_band()
    n = 1000
    r0 = np.linspace(edges[0] * 1.01, edges[-1] * 0.99, n)
    pop = ts.Population(r0, np.ones(n), np.full(n, 0.5 * (lo + hi)), np.ones(n), np.full(n, CT),
                        (np.searchsorted(edges, r0, side="right") - 1).astype(np.int32), np.zeros(n, np.int32))
    res = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, 0, "sobolev_absorb", seed=1, resume=pop.as_dict(), **KW)
    assert res["n_escaped"] == n
    b_out = res["b_out"]
    ct_expected = CT + (b_out * CT - r0) / (1.0 - b_out)
    assert np.allclose(res["ct_esc"], ct_expected, rtol=1e-9)
    obs = ts.escapes(res)
    assert np.allclose(obs["t_obs"], (CT - r0) / C, rtol=1e-9)
    assert np.allclose(obs["e"].sum(), res["accounting"]["E_esc"])


def test_heating_injection_energy_and_shell_shares():
    z3, edges = _zone3(scale=(1e-12, 1e-12, 1e-12))
    z3.rho = np.array([4.0, 1.0, 0.25])
    lo, hi = TZ.tfm.pump_band()
    E = 3.0e40
    n = 20000
    dt = 1e-3                                              # a slab so short nothing moves
    res = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, n, "sobolev_absorb", seed=2, launch="volume",
                 launch_energy=E, t_stop=T_EXP + dt, **KW)
    a = res["accounting"]
    assert abs(a["E_inj_cm"] - E) / E < 1e-12 and a["E_carried_in"] == 0.0
    assert res["n_paused"] == n and abs(a["identity_residual"]) < 1e-12
    v = edges ** 3; m = z3.rho * (v[1:] - v[:-1]); m /= m.sum()
    sh = np.searchsorted(edges * (1 + dt / T_EXP), res["r"], side="right") - 1
    frac = np.bincount(sh, minlength=3) / n
    assert np.abs(frac - m).max() < 5 * np.sqrt(m.max() * (1 - m.max()) / n)
    # closed-form heating energy
    m_shell = np.array([1e30, 2e30]); e_tot, per = ts.heating_energy(m_shell, 4 * 86400.0, 5 * 86400.0)
    num = np.trapezoid(ts.fontes_heating_rate(np.linspace(4 * 86400.0, 5 * 86400.0, 200001)), np.linspace(4 * 86400.0, 5 * 86400.0, 200001))
    assert abs(e_tot - m_shell.sum() * num) / e_tot < 1e-8 and np.allclose(per / per.sum(), m_shell / m_shell.sum())


def test_pause_is_exact_and_resume_is_bit_reproducible():
    z3, edges = _zone3()
    lo, hi = TZ.tfm.pump_band()
    n = 20000; dt = 300.0
    res = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, n, "sobolev_dmacro", seed=3, t_stop=T_EXP + dt, **KW)
    p = res["fate"] == 4
    assert p.sum() > 0
    assert np.all(res["ctime"][p] == C * (T_EXP + dt))
    assert np.all(res["s_acc"][p] == 0.0)
    e_t = edges * ((T_EXP + dt) / T_EXP)
    assert np.all(res["r"][p] >= e_t[0]) and np.all(res["r"][p] <= e_t[-1])
    sh = np.searchsorted(e_t, res["r"][p], side="right") - 1
    assert np.array_equal(np.clip(sh, 0, 2), res["shell_of"][p])
    pop = ts.carried_population(res, e_t)
    r1 = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, 0, "sobolev_dmacro", seed=9, resume=pop.as_dict(),
                t_stop=T_EXP + 2 * dt, **KW)
    r2 = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, 0, "sobolev_dmacro", seed=9, resume=pop.as_dict(),
                t_stop=T_EXP + 2 * dt, **KW)
    for k in ("nu_out_all", "fate", "w", "ctime", "ct_esc", "r", "mu"):
        assert np.array_equal(r1[k], r2[k], equal_nan=True), k


def test_event_cap_books_capped_energy():
    z3, edges = _zone3()
    lo, hi = TZ.tfm.pump_band()
    res = run_mc(z3, edges[0], edges[-1], T_EXP, lo, hi, 20000, "sobolev_dmacro", seed=4, max_events=1, **KW)
    a = res["accounting"]
    assert res["n_capped"] == res["n_interacted"] and res["n_capped"] > 0
    assert abs(a["identity_residual"]) < 1e-12 and a["E_capped"] > 0
    assert np.all(res["n_events"][res["fate"] == 5] == 1)


def test_initial_radiation_energy_and_the_radiation_temperature_round_trip():
    import test_ejecta as TE
    st = TE._state(n_shell=6, profile="power")
    st.T_rad = np.array([5000.0, 4500.0, 4000.0, 3500.0, 3000.0, 2500.0]); st.T_gas = st.T_rad.copy()
    rng = np.random.default_rng(0)
    pop, E_tot, E_s = ts.initial_radiation(st, 60000, rng, 1e14, 3e15, [1, 2, 3, 4])
    V = st.shell_volume()[1:5]
    assert np.allclose(E_s, ts.A_RAD * st.T_rad[1:5] ** 4 * V) and abs(pop.energy_cm() - E_tot) / E_tot < 1e-10
    T, n = ts.t_rad_from_population(pop, st, [1, 2, 3, 4], n_min=50)
    assert np.all(n >= 50) and np.allclose(T, st.T_rad[1:5], rtol=0.01)
    # the lab energy exceeds the comoving one on average by the Doppler boost
    assert pop.energy_lab() > pop.energy_cm()


def test_lightcurve_driver_smoke(tmp_path):
    """The Phase 10b driver end to end on a tiny configuration (8 shells, two
    0.15-d slabs, the resolved leg): checkpoints, atomic run.json, --analyse
    with the global closure, and --resume reproducing a deleted slab bit for
    bit. Needs the Nd cache."""
    import json, shutil, subprocess
    from sobolev.atomic_cache import CACHE_DIR
    if not (CACHE_DIR / "60NdII.npz").exists():
        pytest.skip("Nd cache not built")
    drv = ROOT / "paper4/phase10_fontes/lightcurve.py"
    py = sys.executable
    out = tmp_path / "lc"
    args = [py, str(drv), "--out", str(out), "--n-shell", "8", "--t1", "4.3", "--n-slabs", "2", "--n-init", "400",
            "--n-heat", "200", "--legs", "R2"]
    subprocess.run(args, check=True, capture_output=True, timeout=900)
    run = json.loads((out / "run.json").read_text())
    assert run["done"]["R2"] == [0, 1] and all(abs(t["identity"]) < 1e-12 for t in run["tallies"]["R2"])
    assert np.isclose(run["tallies"]["R2"][1]["E_carried_in"], run["tallies"]["R2"][0]["E_carried_out"], rtol=1e-12)
    subprocess.run([py, str(drv), "--out", str(out), "--analyse"], check=True, capture_output=True, timeout=300)
    s = json.loads((out / "summary.json").read_text())
    assert abs(s["legs"]["R2"]["closure"]) < 1e-10 and (out / "esc_R2_001.npz").exists()
    # resume: drop slab 1 and rerun it from the slab-0 checkpoint
    t1 = dict(run["tallies"]["R2"][1]); run["done"]["R2"] = [0]; run["tallies"]["R2"] = run["tallies"]["R2"][:1]
    shutil.copy(out / "pop_R2.prev.npz", out / "pop_R2.npz")
    (out / "run.json").write_text(json.dumps(run))
    subprocess.run(args + ["--resume"], check=True, capture_output=True, timeout=900)
    run2 = json.loads((out / "run.json").read_text())
    t1b = run2["tallies"]["R2"][1]
    assert all(t1[k] == t1b[k] for k in ("E_esc", "W", "E_carried_out", "seed", "n_paused", "events_max"))
