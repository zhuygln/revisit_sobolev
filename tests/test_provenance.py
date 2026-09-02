"""Provenance closures (2026-09-02): the three results that lived only in the
notebook now have drivers and JSON. These tests pin the committed JSON to the
numbers the report quotes and check the generators reproduce them."""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
P3 = ROOT / "paper3"
sys.path.insert(0, str(P3 / "phase9_audit"))
sys.path.insert(0, str(P3 / "phase8_survey"))
sys.path.insert(0, str(P3 / "synthetic"))


def _load(p):
    if not p.exists():
        pytest.skip(f"{p.name} not present")
    return json.load(open(p))


# --- F38: the counterfactual table is a function of survey.json --------------

def test_counterfactual_table_reproduces_section_4_35():
    from counterfactual_table import table
    src = _load(P3 / "phase8_survey" / "survey.json")
    t = {r["ion"]: r for r in table(src["rows"])}
    assert t["59PrII"]["interaction"] == pytest.approx(-0.064, abs=1e-3)
    assert t["58CeII"]["interaction"] == pytest.approx(-0.043, abs=1e-3)
    assert abs(t["60NdII"]["interaction"]) < 0.01 and abs(t["57LaII"]["interaction"]) < 0.01
    # the redistribution approximation contributes almost nothing
    assert all(abs(r["A"]) <= 0.015 for r in t.values())
    # Pr II: B alone is too bright, B with A is too opaque
    assert t["59PrII"]["B"] > 0 > t["59PrII"]["C"]


def test_counterfactual_json_matches_generator():
    from counterfactual_table import table
    src = _load(P3 / "phase8_survey" / "survey.json")
    saved = _load(P3 / "phase9_audit" / "counterfactual_table.json")
    for a, b in zip(table(src["rows"]), saved["rows"]):
        assert a["ion"] == b["ion"]
        assert a["interaction"] == pytest.approx(b["interaction"], abs=1e-9)


# --- F35: the Ce II density scan (§4.32) -------------------------------------

def test_density_scan_ce_crosses_between_45_and_67():
    d = _load(P3 / "phase8_survey" / "density_scan_58CeII.json")
    rows = sorted(d["rows"], key=lambda r: r["band_S_band"])
    assert len(rows) == 6
    signs = [r["binned"] < 0 for r in rows]
    assert signs == [True, True, False, False, False, False], signs
    assert rows[0]["binned"] == pytest.approx(-0.334, abs=0.05)   # notebook: -33.4 %, one seed
    assert rows[3]["binned"] == pytest.approx(1.246, abs=0.10)    # +124.6 %
    c = d["crossing"]["binned"]
    assert c["direction"] == "neg_to_pos" and 45 < c["S"] < 67
    # the third point is the survey's global normalization
    assert rows[2]["n_ion"] == pytest.approx(d["n_ion_global"], rel=1e-6)
    assert rows[2]["tau_max"] == pytest.approx(5.0, abs=0.01)


def test_density_scan_la_binned_negative_above_global():
    d = _load(P3 / "phase8_survey" / "density_scan_57LaII.json")
    rows = sorted(d["rows"], key=lambda r: r["band_S_band"])
    assert rows[0]["binned"] > 0 and all(r["binned"] < 0 for r in rows[1:])
    assert all(r["expansion"] > 0 for r in rows)          # the practical closure never crosses
    assert d["crossing"]["expansion"] is None


def test_density_scan_crossing_helper():
    from density_scan import crossing
    rows = [{"band_S_band": 10.0, "binned": -0.1}, {"band_S_band": 100.0, "binned": 0.1}]
    c = crossing(rows)
    assert c["direction"] == "neg_to_pos" and c["S"] == pytest.approx(31.62, rel=1e-3)
    assert crossing([{"band_S_band": 1.0, "binned": 0.1}, {"band_S_band": 2.0, "binned": 0.2}]) is None


# --- F39: the exit-tau scan (§4.34b) -----------------------------------------

def test_exit_tau_scan_crosses_neg_to_pos():
    from boundary import crossing
    d = _load(P3 / "synthetic" / "boundary_exit_tau.json")
    assert d["delocalize"] == 1.0 and d["n_exit"] == 6
    got = {}
    for xt in (0.5, 2.0):
        rows = [r for r in d["rows"] if r["exit_tau"] == xt and r["n_lines"] == 100]
        assert rows, f"no rows at exit_tau={xt}"
        c = crossing(rows, "binned")
        assert c is not None and c["direction"] == "neg_to_pos", (xt, c)
        got[xt] = c["S"]
    assert got[0.5] == pytest.approx(179, rel=0.1)
    assert got[2.0] == pytest.approx(688, rel=0.15)
    assert got[2.0] > got[0.5]
