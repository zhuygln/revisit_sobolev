"""Paper IV Phase 5: the dual-role closure D (capped net reprocessing)."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1", ROOT / "paper3/synthetic"):
    sys.path.insert(0, str(p))

from sobolev.energy_balance import capped_reprocessing      # noqa: E402
from forest_mc import run_mc                                # noqa: E402
from forest import synthetic_forest                         # noqa: E402
import test_forest_mc as tfm                                # noqa: E402
from test_macroatom import cascade_atom                     # noqa: E402

R_CORE, R_OUT, T_EXP, T_CORE = 8.64e13, 2.592e14, 86400.0, 6000.0


def _forest():
    atom, _ = synthetic_forest(n_lines=60, tau=6.0, span=0.25, n_exit=2, dlnlam=0.05, seed=3)
    atom.level_energy_cm = np.zeros(atom.lower_all.max() + 1 if atom.upper_all.max() < atom.lower_all.max() else atom.upper_all.max() + 1)
    # level energies from the lines: upper = lower + h nu / hc (a tree; the
    # synthetic forest's exits share lower levels, so set from the pump lines)
    from sobolev.constants import C
    E = np.zeros(max(atom.upper_all.max(), atom.lower_all.max()) + 1)
    for l in np.argsort(atom.nu0_all):
        lo, up = atom.lower_all[l], atom.upper_all[l]
        E[up] = max(E[up], E[lo] + atom.nu0_all[l] / C)
    atom.level_energy_cm = E
    return atom


def test_capped_probability_is_the_stated_ratio():
    atom = _forest()
    edges, E = atom.expansion_bins(4.17e-5, weight="poisson")
    p = capped_reprocessing(atom, edges, E, 4.17e-5)
    assert p.shape == E.shape and np.all((p >= 0) & (p <= 1))
    assert np.all(p[E == 0] == 0)
    # every line here has tau ~ 6 >> dnu/nu, so cap = dnu/nu per line and p = n_lines dnu/nu / E
    b = np.clip(np.searchsorted(edges, atom.op_nu, side="right") - 1, 0, E.size - 1)
    for k in np.unique(b):
        n_l = np.sum(b == k)
        assert np.isclose(p[k], min(1.0, n_l * 4.17e-5 / E[k]))
    # with a cap above every tau the ratio is sum tau / sum (1 - e^-tau) >= 1 -> p = 1
    p1 = capped_reprocessing(atom, edges, E, 1e3)
    assert np.all(p1[E > 0] == 1.0)


def test_capped_closure_is_inert_when_off_and_scatters_when_on():
    atom = _forest()
    lo, hi = atom.op_nu.min() * 0.99, atom.op_nu.max() * 1.01
    kw = dict(seed=1, t_core=T_CORE, packets="energy", launch_weight="energy")
    a = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 20000, "expansion_dmacro", **kw)
    b = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 20000, "expansion_dmacro", reprocess=None, **kw)
    assert np.array_equal(a["nu_out_all"], b["nu_out_all"], equal_nan=True)
    d = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 20000, "expansion_dmacro", reprocess="capped", **kw)
    assert abs(d["accounting"]["identity_residual"]) < 1e-12 and d["accounting"]["E_dep_cm"] == 0.0
    # nearly every interaction is a coherent scattering: the exit-energy tally
    # (lined re-emissions only) is a small fraction of the interacting energy
    assert d["exit_energy"].sum() < 0.2 * a["exit_energy"].sum()
    assert d["n_interactions"] > 0
    # a cap above every tau reproduces B2 exactly (p = 1 everywhere, one extra draw per event)
    e = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 20000, "expansion_dmacro", reprocess="capped", tau_cap=1e3, **kw)
    assert np.isclose(e["exit_energy"].sum(), a["exit_energy"].sum(), rtol=0.05)


def test_capped_reprocessing_rejects_sobolev_legs():
    fa = cascade_atom(2.0, 1.0, 1.0)
    lo, hi = tfm.pump_band()
    with pytest.raises(ValueError):
        run_mc(fa, tfm.R_CORE, tfm.R_OUT, tfm.T_EXP, lo, hi, 100, "sobolev_dmacro", packets="energy", reprocess="capped")
