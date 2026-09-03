"""paper3/freeze.py: the numeric comparison, the manifests and the check tiers.
The suite never runs the full regeneration (that is `freeze.py --check`, ~75 s)."""
import json, math, sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paper3"))
import freeze  # noqa: E402

FROZEN = ROOT / "paper3" / "FROZEN.json"
needs_frozen = pytest.mark.skipif(not FROZEN.exists(), reason="FROZEN.json not present")


def test_compare_tolerance_and_nan():
    a = {"x": [1.0, float("nan"), 2.0], "s": "ok", "n": 3, "d": {"cond": 1e12}}
    b = {"x": [1.0 + 1e-12, float("nan"), 2.0 * (1 + 1e-10)], "s": "ok", "n": 3, "d": {"cond": 1e12 * (1 + 1e-7)}}
    assert freeze.compare(a, b) == []
    c = json.loads(json.dumps(b).replace("2.0000000002", "2.001"))
    c["x"][2] = 2.001
    assert any("/x[2]" in d for d in freeze.compare(a, c))
    assert any("cond" in d for d in freeze.compare(a, {**b, "d": {"cond": 1e12 * 1.01}}))
    assert any("keys differ" in d for d in freeze.compare(a, {**b, "extra": 1}))
    assert any("length" in d for d in freeze.compare(a, {**b, "x": [1.0]}))
    assert freeze.compare({"v": float("nan")}, {"v": 0.0}) != []
    assert freeze.compare({"v": 1}, {"v": True}) != []   # bool is not a number here


def test_input_manifest_complete():
    fs = [freeze.rel(f) for f in freeze.input_files()]
    assert sum(f.startswith("paper3/phase12_grid/grid/model_") for f in fs) == 27
    assert sum(f.startswith("paper3/phase12_grid/grid/tscale/model_") for f in fs) == 4
    assert sum(f.startswith("paper3/phase12_grid/robustness/chain_model_") for f in fs) == 5   # 4 cells + the probe
    assert sum(f.startswith("data/filters/") for f in fs) == 7
    assert "paper3/phase12_grid/sensitivity.py" in fs and "sobolev/photometry.py" in fs
    assert len(fs) == len(set(fs))


def test_check_outputs_detects_edit(tmp_path):
    (tmp_path / "a.json").write_text('{"v": 1}')
    (tmp_path / "f.pdf").write_bytes(b"%PDF")
    fz = {"outputs": {"a.json": freeze.sha(tmp_path / "a.json"), "f.pdf": freeze.sha(tmp_path / "f.pdf"), "gone.json": "0" * 64}}
    fails, warns = freeze.check_outputs(fz, root=tmp_path)
    assert fails == ["output missing: gone.json"] and warns == []
    (tmp_path / "a.json").write_text('{"v": 2}')
    (tmp_path / "f.pdf").write_bytes(b"%PDF-1.4")
    fails, warns = freeze.check_outputs(fz, root=tmp_path)
    assert "output differs from FROZEN: a.json" in fails and warns == ["output differs from FROZEN: f.pdf"]
    fails, warns = freeze.check_outputs(fz, root=tmp_path, strict=True)
    assert "output differs from FROZEN: f.pdf" in fails


def test_check_regenerated_detects_perturbation(tmp_path):
    """A perturbed derived JSON (one number nudged) is caught by tier (c); the
    headline comparison by tier (d)."""
    dest, sd = freeze.scratch_dest(tmp_path / "a"), freeze.scratch_dest(tmp_path / "b")
    for k in freeze.DERIVED:
        for d in (dest, sd):
            Path(d[k]).write_text(json.dumps({"x": [0.5, 1.0], "cond": 1e9}))
    for n in freeze.FIGS:
        for ext in ("png", "pdf"):
            for d in (dest, sd):
                (Path(d["figdir"]) / f"{n}.{ext}").write_bytes(b"fig")
    fz = {"headline": {"k": 1.0}}
    assert freeze.check_regenerated(fz, dest, sd, {"k": 1.0}) == ([], [])
    Path(sd["syserr"]).write_text(json.dumps({"x": [0.5, 1.0 + 1e-6], "cond": 1e9}))
    fails, _ = freeze.check_regenerated(fz, dest, sd, {"k": 1.0})
    assert fails and fails[0].startswith("syserr.json/x[1]")
    fails, _ = freeze.check_regenerated(fz, dest, dest, {"k": 1.1})
    assert fails == ["headline/k: 1.0 vs 1.1"]
    (Path(sd["figdir"]) / "fig2_bol_vs_colour.pdf").write_bytes(b"other")
    fails, warns = freeze.check_regenerated(fz, dest, sd, {"k": 1.0})
    assert any("fig2" in w for w in warns)
    fails, warns = freeze.check_regenerated(fz, dest, sd, {"k": 1.0}, strict=True)
    assert any("fig2" in f for f in fails)


@needs_frozen
def test_headline_matches_committed_json():
    """headline() recomputed from the committed derived files equals FROZEN.headline."""
    fz = json.loads(FROZEN.read_text())
    h = freeze.headline(freeze.canonical())
    assert freeze.compare(fz["headline"], h) == []
    for k in ("gate2.T0.cb", "gate2.T1.cb_by_X", "control.A_redist.cc", "colour.nir_negative", "floor.median",
              "gate3.dense.T0.survives", "gate3.sparse.T0.survives", "gate3.optical.T0.survives",
              "tscale.gas.cos", "chain.dm_change_8000_range", "chain.criterion_met",
              "chain.median_chi2_res_dof_override", "syserr.one_mode.C_both", "syserr.null_median.C_both"):
        assert k in h, k
    assert h["gate2.T0.cb"] == h["gate2.T0.n"] == 27
    assert h["gate2.T1.cb_by_X"] == {"0.001": [1, 9], "0.01": [6, 9], "0.1": [9, 9]}
    assert h["control.A_redist.cc"] == 21 and h["colour.nir_negative"] == [195, 199]
    assert h["chain.criterion_met"] == [4, 12] and h["chain.signs_kept"] == [12, 12]
    assert math.isclose(h["chain.median_chi2_res_dof_override"], 116.2, abs_tol=0.1)
    assert math.isclose(h["gate2.T0.median_chi2_res_dof"], 118.0, abs_tol=0.1)   # the two medians stay distinct
    assert 0.79 < h["syserr.one_mode.C_both"] < 0.81 and h["syserr.null_p95.C_both"] < 0.4


@needs_frozen
def test_frozen_manifests_and_git_fields():
    fz = json.loads(FROZEN.read_text())
    assert fz["tag"] == freeze.TAG
    assert set(fz["git"]["trees"]) == {"grid", "grid/tscale", "robustness", "data/filters"}
    assert set(fz["inputs"]) == {freeze.rel(f) for f in freeze.input_files()}
    assert freeze.check_inputs(fz) == []   # the working tree's inputs are the frozen ones
    outs = set(fz["outputs"])
    assert {"paper3/phase12_grid/sensitivity_chain8000.json", "paper3/phase13_observability/observability.json",
            "paper3/figures/fig5_observability.pdf"} <= outs
    assert all((ROOT / p).exists() for p in outs)
