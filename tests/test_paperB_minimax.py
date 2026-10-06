"""Paper B: the fair scalar comparator (joint minimax eps) from the existing
eps grids, pinned on the frozen records and checked on a synthetic curve."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


MM = _load("paperB_minimax", "paperB/scalar/minimax.py")


def test_minimax_on_the_frozen_records():
    name, row = MM.g1_record("60NdII")
    if row is None:
        pytest.skip("no G1 record")
    assert name == "gate1_60NdII_n1e6.json"
    r = MM.read(row)
    assert r["prereg"]["eps"] == 0.0 and abs(r["prereg"]["max_band"] - 1.0581) < 2e-4
    assert r["joint_minimax"]["eps"] == 0.10 and abs(r["joint_minimax"]["joint"] - 0.8219) < 2e-4
    assert r["joint_minimax"]["max_band"] >= r["joint_minimax"]["max_colour"]
    pc = ROOT / "paperB/gate3/gate3_partc_p1blend.json"
    if pc.exists():
        rb = MM.read(json.loads(pc.read_text()))
        assert rb["prereg"]["eps"] == 0.0 and abs(rb["prereg"]["max_band"] - 1.6906) < 2e-4
        assert rb["joint_minimax"]["eps"] == 0.05 and abs(rb["joint_minimax"]["joint"] - 1.5524) < 2e-4
    out = json.loads((ROOT / "paperB/scalar/minimax.json").read_text())
    assert out["per_record"]["60NdII"]["joint_minimax"]["eps"] == r["joint_minimax"]["eps"]


def test_minimax_picks_the_joint_optimum_not_the_mean_optimum(monkeypatch):
    # three eps values: the mean-optimal one has one bad band and a bad colour; the joint optimum is the middle one
    legs = {"E0.00": dict(eps=0.0, band=dict(mean=0.10, max=0.50), colour=dict(max=0.40)),
            "E0.50": dict(eps=0.5, band=dict(mean=0.20, max=0.25), colour=dict(max=0.20)),
            "E1.00": dict(eps=1.0, band=dict(mean=0.30, max=0.30), colour=dict(max=0.10))}
    monkeypatch.setattr(MM.A, "metrics", lambda row: dict(legs=legs, live_bands=["g", "r"]))
    r = MM.read({})
    assert r["prereg"]["eps"] == 0.0 and r["band_minimax"]["eps"] == 0.5 and r["joint_minimax"]["eps"] == 0.5
    assert r["joint_minimax"]["joint"] == 0.25 and r["joint_minimax"]["interior"] is True and r["prereg"]["interior"] is False
    assert abs(r["joint_gain"] - 0.25) < 1e-12
