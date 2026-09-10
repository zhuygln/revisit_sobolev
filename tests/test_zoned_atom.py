"""Paper IV Phases 8-10: the zoned atom reproduces the single-zone atom at one
shell, bit for bit, and its per-shell tables behave."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1", ROOT / "paper3/synthetic"):
    sys.path.insert(0, str(p))

from sobolev.constants import C                                # noqa: E402
from sobolev.zoned_atom import ZonedAtom, ZonedMacroAtom        # noqa: E402
from sobolev.macroatom import DownwardMacroAtom                 # noqa: E402
from sobolev import atomic_cache as ac                          # noqa: E402
from forest_mc import ForestAtom                                # noqa: E402
from forest import synthetic_forest                             # noqa: E402
from test_macroatom import cascade_atom, LEVEL_E, NU_13, NU_32, NU_21, F_OSC, T_EXP   # noqa: E402
from sobolev.optical_depth import tau_sobolev                   # noqa: E402

R_EDGES = np.array([8.64e12, 4.32e13])


def _cascade_arrays(tau13, a31, a32, tau32=0.0):
    n1 = tau13 / tau_sobolev(F_OSC, 1.0, C / NU_13, T_EXP)
    n2 = tau32 / tau_sobolev(F_OSC, 1.0, C / NU_32, T_EXP) if tau32 > 0 else 0.0
    return dict(nu0=[NU_13, NU_32, NU_21], f_osc=[F_OSC, F_OSC, 0.0], n_lower=[n1, n2, 0.0],
                n_upper=[0.0, 0.0, 0.0], A=[a31, a32, 1.0], lower=[1, 2, 1], upper=[3, 3, 2])


def test_zoned_macroatom_matches_downward_macroatom_at_one_shell():
    fa = cascade_atom(3.0, 1.0, 1.0, tau32=6.0)
    dm = fa.dmacro()
    zm = ZonedMacroAtom(fa.nu0_all, fa.A_all, fa.lower_all, fa.upper_all, fa.beta_all[None, :], LEVEL_E)
    assert zm.n_entries == dm.n_entries and np.array_equal(zm.line, dm.line) and np.array_equal(zm.kind, dm.kind)
    assert np.array_equal(zm.off, dm.off) and np.array_equal(zm.has_exit, dm.has_exit)
    rng = np.random.default_rng(0)
    lev = rng.choice(np.flatnonzero(dm.has_exit), 100000); v = rng.uniform(size=100000)
    assert np.array_equal(zm.sample(lev, v, np.zeros(lev.size, int)), dm.sample(lev, v))
    ul, pd_, pi_ = zm.probabilities(3); ul2, pd2, pi2 = dm.probabilities(3)
    assert np.array_equal(ul, ul2) and np.allclose(pd_, pd2) and np.allclose(pi_, pi2)
    # the walk consumes the same draws
    ex, nj, dead = zm.walk(np.full(5000, 3), np.random.default_rng(3), np.zeros(5000, int))
    ex2, nj2, dead2 = dm.walk(np.full(5000, 3), np.random.default_rng(3))
    assert np.array_equal(ex, ex2) and np.array_equal(nj, nj2)


def test_identical_shells_share_every_segment_and_differ_when_beta_differs():
    fa = cascade_atom(3.0, 1.0, 1.0, tau32=6.0)
    b = fa.beta_all
    zm = ZonedMacroAtom(fa.nu0_all, fa.A_all, fa.lower_all, fa.upper_all, np.array([b, b]), LEVEL_E)
    assert not zm.dep.any() and zm.n_dep_entries == 0
    assert np.array_equal(zm.start[0], zm.start[1])
    b2 = b.copy(); b2[1] *= 0.5                      # line 3->2 more trapped in shell 1
    zm2 = ZonedMacroAtom(fa.nu0_all, fa.A_all, fa.lower_all, fa.upper_all, np.array([b, b2]), LEVEL_E)
    assert zm2.dep[3] and not zm2.dep[2]
    _, p0, _ = zm2.probabilities(3, 0); _, p1, _ = zm2.probabilities(3, 1)
    assert p1[0] > p0[0]                              # 3->1 wins more when 3->2 is trapped
    _, q0, _ = zm.probabilities(3, 0)
    assert np.array_equal(p0, q0)


def test_zoned_atom_one_shell_equals_forest_atom():
    a = _cascade_arrays(3.0, 1.0, 1.0, tau32=6.0)
    fa = ForestAtom(**a, t_exp=T_EXP, tau_min=1e-6, stim=False); fa.level_energy_cm = LEVEL_E
    za = ZonedAtom(a["nu0"], a["f_osc"], [a["n_lower"]], [a["n_upper"]], a["A"], a["lower"], a["upper"],
                   R_EDGES, T_EXP, tau_min=1e-6, stim=False, level_energy_cm=LEVEL_E)
    assert za.n_shell == 1 and za.n_opacity == fa.n_opacity
    assert np.array_equal(za.op_idx, fa.op_idx) and np.array_equal(za.op_nu, fa.op_nu)
    assert np.array_equal(za.op_tau[0], fa.op_tau) and np.array_equal(za.op_p[0], fa.op_p)
    assert np.array_equal(za.op_beta[0], fa.beta_all[fa.op_idx]) and np.array_equal(za.beta_all, fa.beta_all)
    assert np.array_equal(za.op_nxt[0], np.arange(za.n_opacity))
    assert np.array_equal(za.macro.cum, fa.dmacro().cum)


def test_two_shells_union_and_skip_table():
    a = _cascade_arrays(3.0, 1.0, 1.0, tau32=6.0)
    n_low = np.array([a["n_lower"], [0.0, a["n_lower"][1], 0.0]])      # shell 1: only 3->2 has opacity
    za = ZonedAtom(a["nu0"], a["f_osc"], n_low, [a["n_upper"]] * 2, a["A"], a["lower"], a["upper"],
                   np.array([8.64e12, 2e13, 4.32e13]), T_EXP, tau_min=1e-6, stim=False, level_energy_cm=LEVEL_E)
    assert za.n_shell == 2 and za.n_opacity == 2                        # union: both pumped lines
    assert list(za.n_opacity_shell) == [2, 1]
    k32 = int(np.flatnonzero(za.op_nu == NU_32)[0]); k13 = int(np.flatnonzero(za.op_nu == NU_13)[0])
    assert za.op_p[1, k13] == 0.0 and za.op_beta[1, k13] == 1.0 and za.op_p[1, k32] > 0
    assert za.op_nxt[1, k13] == (k32 if k32 < k13 else -1)
    assert za.macro.dep[3]                                              # beta of the pump differs


@pytest.mark.skipif(not (ac.CACHE_DIR / "57LaII.npz").exists(), reason="La II cache not built")
def test_from_state_one_shell_equals_from_cached(tmp_path):
    from sobolev.ejecta import EjectaState
    n = 3
    st = EjectaState(t=2 * 86400.0, v_edges=np.array([0.09, 0.10, 0.11, 0.115]) * C, rho=np.full(n, 3e-15),
                     T_gas=np.full(n, 3400.0), T_rad=np.full(n, 3400.0),
                     X={"La": np.full(n, 4e-3), "bulk": np.full(n, 0.996)}, f_ion={"La II": np.ones(n)})
    from sobolev.abundances import ATOMIC_MASS
    n_ion = st.n_ion("La", "II", 1, ATOMIC_MASS["La"])
    fa = ForestAtom.from_cached([("57LaII", n_ion)], 3400.0, st.t, tau_min=1e-3)
    za = ZonedAtom.from_state(st, [1], tau_min=1e-3, emis_cut=None)
    assert np.array_equal(za.op_idx, fa.op_idx) and np.array_equal(za.op_tau[0], fa.op_tau)
    assert np.array_equal(za.op_beta[0], fa.beta_all[fa.op_idx])
    assert np.array_equal(za.macro.cum, fa.dmacro().cum) and np.array_equal(za.macro.dead_end, fa.dmacro().dead_end)
    s = za.thermal_sampler(0, weight="energy_beta"); f = fa.thermal_sampler(weight="energy_beta")
    u = np.random.default_rng(1).uniform(size=100000)
    assert np.array_equal(s(u), f(u))
    # the cut sampler drops < emis_cut of the weight
    zc = ZonedAtom.from_state(st, [1], tau_min=1e-3, emis_cut=1e-4)
    sc = zc.thermal_sampler(0, weight="energy_beta")
    assert sc.n_lines < s.n_lines and np.mean(sc(u) != f(u)) < 1e-3
    # three shells: nested opacity sets, skip table non-trivial
    z3 = ZonedAtom.from_state(st, [0, 1, 2], tau_min=1e-3)
    assert z3.n_shell == 3 and z3.n_opacity == z3.n_opacity_shell.max()
