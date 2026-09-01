"""The ejecta model behind F40, and the crossing interpolator §10 added.

`state()` gained a `core_law` argument so §10 could cool the illuminating core.
The one thing that must not happen is F40 quietly moving underneath it, so the
first test pins the default against the committed trajectory JSON that produced
the published numbers.

`crossing_epoch` exists because F40's headline -- t = 1.17 d, S = 47.5 -- was
interpolated by hand and stored nowhere.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3",
          ROOT / "paper3/phase0_reference", ROOT / "paper3/phase10_kilonova",
          ROOT / "paper3/phase11_observables"):
    sys.path.insert(0, str(p))
from trajectory import EPOCHS, state
from observables import crossing_epoch

TRAJ = ROOT / "paper3/phase10_kilonova/trajectory_58CeII.json"


def test_default_core_law_reproduces_the_published_trajectory():
    """Every epoch of the committed F40 run must fall out of the default.

    This is the regression guard for the `core_law` change: §10 cools the core,
    §4.36 did not, and the published numbers belong to the frozen one.
    """
    rows = json.loads(TRAJ.read_text())["rows"]
    assert len(rows) == len(EPOCHS)
    for row in rows:
        st = state(row["t_d"])
        assert st["t_exp"] == pytest.approx(row["t_d"] * 86400.0)
        assert st["rho"] == pytest.approx(row["rho"], rel=1e-12)
        assert st["T_gas"] == pytest.approx(row["T_gas"], rel=1e-12)
        assert st["t_core"] == 6000.0        # frozen, as §4.36 ran it
        assert st["core_law"] == "fixed"


def test_cool_core_tracks_the_gas():
    for t_d in EPOCHS:
        st = state(t_d, "cool")
        assert st["t_core"] == pytest.approx(st["T_gas"])
        # everything else is untouched by the core law
        base = state(t_d)
        for k in ("t_exp", "rho", "T_gas", "r_core", "r_out"):
            assert st[k] == base[k]


def test_unknown_core_law_raises():
    with pytest.raises(ValueError):
        state(1.0, "photospheric")


def _rows(values, S=None):
    """Minimal rows in the shape crossing_epoch consumes."""
    S = S if S is not None else [10.0 * 2 ** -i for i in range(len(values))]
    return [{"t_d": 0.5 * (i + 1), "band_S_band": s,
             "legs": {"C": {"dF_b3800": v}}}
            for i, (v, s) in enumerate(zip(values, S))]


def test_crossing_detected_in_both_directions():
    """The first version of the sibling helper in synthetic/boundary.py tested
    `a < 0 <= b` only and silently missed every positive-to-negative crossing,
    in an experiment whose entire subject is the sign."""
    down = crossing_epoch(_rows([+0.2, -0.2]), "C")
    assert down["direction"] == "pos_to_neg"
    assert down["t_d"] == pytest.approx(0.75)

    up = crossing_epoch(_rows([-0.2, +0.2]), "C")
    assert up["direction"] == "neg_to_pos"
    assert up["t_d"] == pytest.approx(0.75)


def test_crossing_interpolates_linearly_in_t_and_logarithmically_in_S():
    r = crossing_epoch(_rows([+3.0, -1.0], S=[100.0, 1.0]), "C")
    assert r["t_d"] == pytest.approx(0.5 + 0.75 * 0.5)      # f = 3/4
    assert r["S_band"] == pytest.approx(100.0 * (0.01 ** 0.75))


def test_no_crossing_returns_none_and_skips_dead_epochs():
    assert crossing_epoch(_rows([-0.3, -0.1, -0.05]), "C") is None
    rows = _rows([-0.2, +0.2])
    rows.insert(1, {"t_d": 0.6, "skipped": "too few opacity lines"})
    rows.append({"t_d": 2.0, "band_S_band": 0.1,
                 "legs": {"C": {"dF_b3800": None}}})
    assert crossing_epoch(rows, "C")["direction"] == "neg_to_pos"


def test_only_the_first_crossing_is_reported():
    """Late epochs wander through zero on MC noise; the crossing that matters
    is the first, and reporting a noise crossing as physics would be worse than
    reporting none."""
    r = crossing_epoch(_rows([-0.4, +0.4, -0.01, +0.01]), "C")
    assert r["t_d"] == pytest.approx(0.75)
