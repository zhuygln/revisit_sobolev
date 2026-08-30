"""Memory depth m for the grouped-opacity legs.

The m-sweep is only a controlled comparison if the depths share one RNG stream
and m = 0 is a literal no-op, so those are pinned exactly rather than
statistically. The credit draws no random numbers, which is what makes that
possible.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/synthetic"):
    sys.path.insert(0, str(p))
from forest import synthetic_forest
from forest_mc import band_ratio, run_mc
from redistribution import RedistributionKernel

R_CORE, R_OUT, T_EXP, T_CORE = 8.64e13, 2.592e14, 86400.0, 6000.0


@pytest.fixture(scope="module")
def setup():
    atom, _ = synthetic_forest(n_lines=60, tau=6.0, span=0.25, n_exit=2,
                               dlnlam=0.05, seed=3)
    lo, hi = atom.op_nu.min() * 0.99, atom.op_nu.max() * 1.01
    res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 40000, "sobolev_branch",
                 seed=1, t_core=T_CORE, collect_events=True)
    e = res["events"]
    kern = RedistributionKernel.from_branching_mc(e[0], e[1], np.ones(e[0].size),
                                                  16, nu_lo=lo, nu_hi=hi)
    return atom, lo, hi, kern


def _run(setup, mode, m, seed=1, n=30000):
    atom, lo, hi, kern = setup
    return run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode,
                  seed=seed, t_core=T_CORE, kernel=kern, line_memory=m)


def test_m0_is_bit_identical_to_no_memory(setup):
    """The credit consumes no RNG, so this is exact, not statistical."""
    a = _run(setup, "binned_group", 0)
    b = _run(setup, "binned_group", False)
    assert np.array_equal(a["nu_out_all"], b["nu_out_all"], equal_nan=True)
    assert np.array_equal(a["n_events"], b["n_events"])
    assert np.array_equal(a["fate"], b["fate"])


def test_m1_reproduces_the_boolean_form(setup):
    a = _run(setup, "binned_group", 1)
    b = _run(setup, "binned_group", True)
    assert np.array_equal(a["nu_out_all"], b["nu_out_all"], equal_nan=True)
    assert np.array_equal(a["n_events"], b["n_events"])


def test_memory_reduces_interactions(setup):
    """Crediting optical depth away can only make a packet travel further."""
    ev = [float(_run(setup, "binned_group", m)["n_events"].mean())
          for m in (0, 1, 2, 4)]
    assert ev[1] < ev[0]
    # deeper memory keeps crediting, so the count cannot rise back above m=0
    assert all(e <= ev[0] + 1e-12 for e in ev)


def test_deeper_memory_is_monotone_in_interactions(setup):
    ev = [float(_run(setup, "binned_group", m)["n_events"].mean())
          for m in (1, 2, 4, 8)]
    assert ev[0] >= ev[1] >= ev[2] >= ev[3] - 1e-12


def test_memory_is_inert_for_the_sobolev_legs(setup):
    """Sobolev transport has its own at-resonance skip and no tau_r at all."""
    a = _run(setup, "sobolev_group", 0)
    b = _run(setup, "sobolev_group", 8)
    assert np.array_equal(a["nu_out_all"], b["nu_out_all"], equal_nan=True)


def test_negative_depth_is_rejected(setup):
    with pytest.raises(ValueError, match="line_memory"):
        _run(setup, "binned_group", -1)


def test_all_depths_share_one_rng_stream(setup):
    """Different m must differ only through the credit, never through a
    different sequence of random draws -- otherwise the sweep confounds memory
    with Monte Carlo noise."""
    a = _run(setup, "binned_group", 0)
    b = _run(setup, "binned_group", 8)
    # launch is drawn before any interaction, so it must be untouched
    assert np.array_equal(a["nu_launch"], b["nu_launch"])
