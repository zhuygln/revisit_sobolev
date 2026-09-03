"""`grid_table.py` and `robustness.chain_summary` as data: the numbers the paper quotes."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
P12 = ROOT / "paper3" / "phase12_grid"
sys.path.insert(0, str(P12))
import grid_table as G
import robustness as RB

HAVE_GRID = len(list((P12 / "grid").glob("model_M*.json"))) == 27
HAVE_CHAIN = len(list((P12 / "robustness").glob("chain_model_*_t[0-9].json"))) == 4


def _cell(point, t, ran=True, n_used=200_000, redo=None, dm=None, dcolor=None, a_dm=None, floor=0.02):
    dm = dm or {"g": -1.0, "K": 2.0}
    dcolor = dcolor or {"g-r": -0.5, "i-J": -1.0, "J-K": -0.8}
    return {"point": list(point), "t_d": t, "status": "ok" if ran else "over_budget", "ran": ran,
            "floored": False, "redo_budget_s": redo, "n_used": n_used, "trapped_frac": 0.01,
            "live": list(dm), "floor": floor,
            "legs": {"C_both": {"dm": dm, "dcolor": dcolor, "worst_dcolor": max(abs(v) for v in dcolor.values()), "worst_key": "i-J"},
                     "C_binned": {"dm": dm, "dcolor": dcolor, "worst_dcolor": 1.0, "worst_key": "i-J"},
                     "A_redist": {"dm": a_dm or {"g": 0.01, "K": -0.02}, "dcolor": {}, "worst_dcolor": np.nan, "worst_key": ""}}}


def test_nir_sign_count_counts_live_nir_colours_only():
    cells = [_cell((0.01, 0.1, 0.01), 1.0, dcolor={"g-r": -0.5, "i-J": -1.0, "J-K": 0.3}),
             _cell((0.01, 0.1, 0.01), 2.0, dcolor={"g-r": -0.5, "J-K": -0.3}),
             _cell((0.01, 0.1, 0.01), 3.0, ran=False)]
    assert G.nir_sign_count(cells) == (2, 3)


def test_band_extremes_and_floor_stats_filter():
    cells = [_cell((0.01, 0.1, 0.01), 1.0, dm={"g": -2.0, "K": 3.0}, floor=0.02),
             _cell((0.01, 0.1, 0.01), 2.0, dm={"g": -0.3}, floor=0.03, n_used=20_000, redo=5400),
             _cell((0.01, 0.1, 0.01), 3.0, dm={"g": -0.5}, floor=np.nan),
             _cell((0.01, 0.1, 0.01), 5.0, ran=False)]
    be = G.band_extremes(cells)
    assert be["g"]["n"] == 3 and be["g"]["min"] == -2.0 and be["g"]["max"] == -0.3
    assert be["K"]["n"] == 1 and be["r"] is None
    fs = G.floor_stats(cells, n_min=100_000)
    assert fs["n_well"] == 1 and fs["median"] == 0.02 and fs["redone_range"] == [0.03, 0.03]


def test_trapped_fraction():
    assert G.trapped_fraction({"ref": {"n_trapped": 300}, "n_used": 1000}) == 0.1
    assert np.isnan(G.trapped_fraction({"ref": {"n_trapped": 0}}))


def test_per_point_aggregates_worst_cell_with_its_floor():
    cells = [_cell((0.01, 0.1, 0.01), 1.0, dcolor={"i-J": -0.4}, floor=0.02),
             _cell((0.01, 0.1, 0.01), 2.0, dcolor={"i-J": -1.4}, floor=0.05, n_used=20_000, redo=5400),
             _cell((0.01, 0.1, 0.01), 3.0, ran=False)]
    p = G.per_point(cells)[0]
    assert p["ran"] == 2 and p["over_budget"] == 1 and p["redone_epochs"] == [2.0]
    assert p["C_both"]["worst_dcolor"] == 1.4 and p["C_both"]["floor_at_cell"] == 0.05 and p["C_both"]["t_d"] == 2.0
    assert p["floor_max_well"] == 0.02 and p["floor_max_all"] == 0.05


@pytest.mark.skipif(not HAVE_GRID, reason="grid not present")
def test_committed_grid_numbers():
    models = G.load(P12 / "grid")
    s = G.summary(models)
    assert s["nir_negative"] == [195, 199]
    assert s["n_redone_cells"] == 9 and len(s["cells"]) == 162
    assert abs(s["floor"]["median"] - 0.021) < 0.001 and s["floor"]["n_well"] == 62
    assert s["band_extremes"]["C_both"]["g"]["max"] < 0 and s["band_extremes"]["C_both"]["K"]["max"] > 3.5
    # the Markdown renderers and the data agree on the per-point worst colour
    md = G.per_model(models)
    p = next(q for q in s["points"] if q["point"] == [0.01, 0.1, 0.01])
    assert f"{p['C_both']['worst_dcolor']:.2f} ± {p['C_both']['floor_at_cell']:.2f}" in md


@pytest.mark.skipif(not HAVE_CHAIN, reason="chain runs not present")
def test_chain_summary_reproduces_the_report():
    d = RB.chain_summary()
    s = d["summary"]
    assert s["n_cells"] == 4 and s["top_cap"] == "8000" and s["stored_reproduced"]
    assert s["C_both"]["signs_kept_top"] == [12, 12] and s["C_both"]["criterion_met_top"] == [4, 12]
    lo, hi = s["C_both"]["dm_change_range_top"]
    assert 0.13 < lo < 0.15 and 0.20 < hi < 0.22
    assert 0.10 < s["trapped_range_base"][0] < 0.11 and 0.138 < s["trapped_range_base"][1] < 0.139
    for cell in d["cells"]:
        for c, run in cell["runs"].items():
            if c != "2000":
                assert run["B_opacity_mags_identical"]
