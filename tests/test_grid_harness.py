"""Phase 12 harness flags added for §4.41 (chain-limit and budget reruns, the
T-direction check): `run_epoch(n_override, chain_max, atom)`, the row's
provenance keys, `model_name`'s t_scale suffix and `run_grid.merge_redo`."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper3/phase12_grid", ROOT / "paper3/synthetic"):
    sys.path.insert(0, str(p))
import grid                                  # noqa: E402
import run_grid                              # noqa: E402
from forest import synthetic_forest          # noqa: E402
from sobolev.source import SourceModel       # noqa: E402


@pytest.fixture(scope="module")
def cell():
    st = SourceModel(0.01, 0.1).state(2.0)
    atom, _ = synthetic_forest(n_lines=30, tau=3.0, t_exp=st["t_exp"], seed=3)
    return st, (atom, {"toy": 1.0})


def test_model_name_suffix():
    assert grid.model_name(0.01, 0.1, 0.01) == "model_M0.01_v0.1_X0.01"
    assert grid.model_name(0.01, 0.1, 0.01, 1.0) == "model_M0.01_v0.1_X0.01"
    assert grid.model_name(0.01, 0.1, 0.01, 1.25) == "model_M0.01_v0.1_X0.01_T1.25"


def test_run_epoch_override_skips_the_probe_and_records_provenance(cell):
    st, atom = cell
    r = grid.run_epoch(st, 0.01, 300_000, n_override=1500, chain_max=50, atom=atom)
    assert r["status"] == "reduced_n" and r["n_used"] == 1500 and r["n_override"] is True
    assert "probe_s_per_packet" not in r
    assert r["chain_max"] == 50 and r["budget_s"] == grid.BUDGET_S
    assert r["source"]["v_ph_floored"] is False and r["source"]["t_scale"] == 1.0
    assert set(r["legs"]) == {t for t, _ in grid.LEGS}
    assert np.isfinite(r["legs"]["C_both"]["dcolor"]["g-r"])


def test_run_epoch_same_seeds_reproduce_and_chain_max_reaches_run_mc(cell):
    st, atom = cell
    a = grid.run_epoch(st, 0.01, 300_000, n_override=1500, chain_max=50, atom=atom)
    b = grid.run_epoch(st, 0.01, 300_000, n_override=1500, chain_max=50, atom=atom)
    assert a["ref"]["mags"] == b["ref"]["mags"] and a["legs"]["B_opacity"]["mags"] == b["legs"]["B_opacity"]["mags"]
    # a chain cap of 1 thermalizes every re-absorption chain that did not leave
    # its line on the first draw: the reference's trapped count must rise
    c = grid.run_epoch(st, 0.01, 300_000, n_override=1500, chain_max=1, atom=atom)
    assert c["ref"]["n_trapped"] > a["ref"]["n_trapped"]


def test_merge_redo_replaces_one_row_and_keeps_the_header(tmp_path):
    rows = [{"t_d": 0.5, "status": "over_budget"}, {"t_d": 1.0, "status": "ok", "n_used": 100}]
    model = {"m_ej_msun": 0.01, "v_ej_c": 0.1, "x_lan": 0.1, "budget_s": 1500.0,
             "complete": True, "rows": rows}
    (tmp_path / "model_M0.01_v0.1_X0.1.json").write_text(json.dumps(model))
    redo = dict(model, budget_s=5400.0, git="abc",
                rows=[{"t_d": 0.5, "status": "reduced_n", "n_used": 20000}])
    (tmp_path / "redo").mkdir()
    (tmp_path / "redo" / "model_M0.01_v0.1_X0.1_t0.5.json").write_text(json.dumps(redo))
    assert run_grid.redo_cells(tmp_path) == [(0.01, 0.1, 0.1, 0.5, 0.0)]
    merged = run_grid.merge_redo(tmp_path, verbose=False)
    assert merged == [("model_M0.01_v0.1_X0.1.json", 0.5, "reduced_n", 20000)]
    d = json.loads((tmp_path / "model_M0.01_v0.1_X0.1.json").read_text())
    assert d["budget_s"] == 1500.0 and d["complete"] is True
    assert [r["status"] for r in d["rows"]] == ["reduced_n", "ok"]
    assert d["rows"][0]["redo"] == {"budget_s": 5400.0, "git": "abc",
                                    "file": "model_M0.01_v0.1_X0.1_t0.5.json",
                                    "previous_status": "over_budget"}
    assert run_grid.redo_cells(tmp_path) == []
