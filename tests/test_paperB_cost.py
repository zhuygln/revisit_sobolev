"""Paper B: the cost layers are read from the existing records (no transport)."""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


C = _load("paperB_costs", "paperB/cost/costs.py")


def test_layers_read_the_record_fields():
    row = dict(n=1000, legs={"R2build": dict(t_wall=10.0, seeds=[101, 102], events_per_packet=5.0),
                             "R2": dict(t_wall=8.0, seeds=[1, 2, 3], events_per_packet=5.0, n_interactions=100),
                             "A2_ng4": dict(t_wall=9.0, seeds=[1, 2, 3], events_per_packet=4.0, n_interactions=90),
                             "E0.10": dict(t_wall=40.0, seeds=[1, 2, 3], events_per_packet=20.0, n_interactions=900)},
               kernels={"A2_ng4": dict(ng=4, n_events=5000, n_exit_samples=77, serialized_bytes=2048)})
    d = C.layers(row, "A2_ng4")
    assert d["offline"] == dict(build_wall_s=10.0, build_packets=2000, build_events_per_packet=5.0, kernel_events=5000)
    assert d["stored"]["matrix_bytes"] == 128 and d["stored"]["exit_lines"] == 77 and d["stored"]["serialized_kb"] == 2.0
    assert d["online"]["op_over_ref"] == 9.0 / 8.0 and d["online"]["scalar_over_ref"] == 5.0


def test_costs_json_is_the_regeneration_of_the_records():
    p = ROOT / "paperB/cost/costs.json"
    if not p.exists() or not (ROOT / "paperB/gate3/gate3_partc_p1blend.json").exists():
        pytest.skip("no records")
    out = json.loads(p.read_text())
    assert set(out["per_record"]) == {"57LaII", "58CeII", "60NdII", "blend13"}
    nd = out["per_record"]["60NdII"]
    assert nd["record"] == "paperB/gate1/gate1_60NdII_n1e6.json" and nd["stored"]["ng"] == 2 and nd["stored"]["exit_lines"] == 86187
    assert 0.9 < nd["online"]["op_over_ref"] < 1.2 and nd["online"]["scalar_over_ref"] > 3
    assert "does not remove the reference calculation" in out["note"]
