"""Synthetic forests: the dials must be independent and faithful.

The phase diagram is only meaningful if crowding, saturation and redistribution
range can be varied one at a time. These tests pin that orthogonality, and pin
the limit that makes the range axis interpretable: n_exit = 1 must be exactly
coherent scattering, not approximately.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3/synthetic"):
    sys.path.insert(0, str(p))
from forest import synthetic_forest, synthetic_ladder
from forest_mc import run_mc

from sobolev.forest_stats import SE_sums, crowding, redistribution_range, saturation_stats

R_CORE, R_OUT, T_EXP, T_CORE = 8.64e13, 2.592e14, 86400.0, 6000.0


def _events(atom, n=60000, seed=1):
    lo, hi = atom.op_nu.min() * 0.99, atom.op_nu.max() * 1.01
    res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_branch",
                 seed=seed, t_core=T_CORE, collect_events=True)
    return res["events"], (lo, hi)


def test_exit_channels_carry_no_opacity():
    """n_exit > 1 adds branching channels but must not add opacity lines."""
    a1, _ = synthetic_forest(n_lines=50, n_exit=1)
    a4, _ = synthetic_forest(n_lines=50, n_exit=4)
    assert a1.n_opacity == 50 and a4.n_opacity == 50
    assert a4.n_lines_total > a1.n_lines_total


def test_target_tau_is_realized():
    for tau in (0.1, 1.0, 5.0, 25.0):
        a, _ = synthetic_forest(n_lines=30, tau=tau)
        assert saturation_stats(a)["tau_max"] == pytest.approx(tau, rel=1e-6)


def test_saturation_is_independent_of_the_redistribution_dial():
    """E/S must not move when only the exit range changes."""
    ref = None
    for d in (0.01, 0.1, 0.5):
        a, _ = synthetic_forest(n_lines=80, tau=5.0, n_exit=2, dlnlam=d)
        v = SE_sums(a)["E_over_S"]
        if ref is None:
            ref = v
        assert v == pytest.approx(ref, rel=1e-12)


def test_crowding_is_independent_of_the_redistribution_dial():
    ref = None
    for d in (0.01, 0.1, 0.5):
        a, _ = synthetic_forest(n_lines=80, tau=5.0, span=0.2, n_exit=2, dlnlam=d)
        v = crowding(a)["n_sat_per_lnlam"]
        if ref is None:
            ref = v
        assert v == pytest.approx(ref, rel=1e-12)


def test_crowding_scales_with_line_count_and_span():
    a, _ = synthetic_forest(n_lines=100, span=0.2, tau=5.0)
    b, _ = synthetic_forest(n_lines=200, span=0.2, tau=5.0)
    c, _ = synthetic_forest(n_lines=100, span=0.4, tau=5.0)
    ca, cb, cc = (crowding(x)["n_sat_per_lnlam"] for x in (a, b, c))
    assert cb == pytest.approx(2 * ca, rel=0.02)
    assert cc == pytest.approx(ca / 2, rel=0.02)


def test_star_topology_is_exactly_coherent():
    """n_exit = 1: the only way down is the absorbing line, so the photon must
    re-emerge at exactly the frequency it was absorbed at."""
    a, _ = synthetic_forest(n_lines=60, tau=5.0, n_exit=1)
    (nu_in, nu_out, w), (lo, hi) = _events(a)
    assert nu_in.size > 0
    assert np.array_equal(nu_in, nu_out)
    r = redistribution_range(nu_in, nu_out, w,
                             edges=np.geomspace(lo, hi, 33))
    assert r["mean_abs_dlnlam"] == pytest.approx(0.0, abs=1e-15)
    assert r["same_group_frac"] == pytest.approx(1.0)


def test_the_dial_controls_the_measured_range_monotonically():
    got = []
    for d in (0.02, 0.08, 0.3):
        a, _ = synthetic_forest(n_lines=80, tau=5.0, span=0.25, n_exit=2, dlnlam=d)
        (nu_in, nu_out, w), _ = _events(a)
        got.append(redistribution_range(nu_in, nu_out, w)["mean_abs_dlnlam"])
    assert got[0] < got[1] < got[2]
    # the A*beta weighting suppresses the return channel, so the measured range
    # is a fixed fraction of the dial rather than equal to it
    assert 0.5 < got[2] / 0.3 < 1.0


def test_f_return_controls_the_same_group_fraction():
    """Independently of how far exits reach: how often the photon comes back to
    its own resonance."""
    edges = None
    fracs = []
    for fr in (0.1, 0.9):
        a, _ = synthetic_forest(n_lines=80, tau=5.0, span=0.25, n_exit=2,
                                dlnlam=0.1, f_return=fr)
        (nu_in, nu_out, w), (lo, hi) = _events(a)
        edges = np.geomspace(lo, hi, 33)
        fracs.append(redistribution_range(nu_in, nu_out, w, edges=edges)["same_group_frac"])
    assert fracs[1] > fracs[0]


def test_ladder_builds_and_redistributes():
    a, info = synthetic_ladder(n_lines=60, tau=5.0, n_rungs=4)
    assert a.n_opacity == 60
    (nu_in, nu_out, w), _ = _events(a)
    assert nu_in.size > 0
    assert redistribution_range(nu_in, nu_out, w)["mean_abs_dlnlam"] > 0


def test_group_transport_runs_on_a_synthetic_forest():
    """The whole Paper III chain must work on a synthetic atom, or the phase
    diagram cannot be measured at all."""
    sys.path.insert(0, str(ROOT / "paper3"))
    from redistribution import RedistributionKernel
    a, _ = synthetic_forest(n_lines=80, tau=5.0, n_exit=2, dlnlam=0.1)
    (nu_in, nu_out, w), (lo, hi) = _events(a)
    k = RedistributionKernel.from_branching_mc(nu_in, nu_out, w, 16,
                                               nu_lo=lo, nu_hi=hi)
    assert k.validate_energy() < 1e-10
    res = run_mc(a, R_CORE, R_OUT, T_EXP, lo, hi, 20000, "sobolev_group",
                 seed=2, t_core=T_CORE, kernel=k)
    assert np.isfinite(res["nu_out_all"][res["fate"] == 1]).all()
