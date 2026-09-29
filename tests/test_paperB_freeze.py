"""Paper B: the frozen record is a pure function of the committed gate
records, the macros are defined for everything the manuscript uses, and the
figures are generated from the frozen record alone."""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


F = _load("paperB_freeze", "paperB/freeze.py")
LT = _load("paperB_tables", "docs/paperB/latex_tables.py")


def test_frozen_json_is_the_regeneration_of_the_records():
    h = F.build()
    txt = json.dumps(h, indent=1, default=float, sort_keys=True) + "\n"
    assert (ROOT / "paperB/FROZEN.json").read_text() == txt, "run paperB/freeze.py"


def test_the_frozen_record_carries_the_headline_results():
    h = json.loads((ROOT / "paperB/FROZEN.json").read_text())
    assert h["g1"]["B1"] == "GREEN" and h["g1"]["decision"] == "CONTINUE"
    assert h["g2"]["H1"] == "GREEN" and h["g2"]["H2"] == "YELLOW" and h["g2"]["decision"] == "WRITE"
    g3 = h["g3"]
    assert (g3["C1"], g3["C2"], g3["C3"], g3["C4"], g3["C5"]) == ("RED", "RED", "YELLOW", "YELLOW", "GREEN")
    assert g3["per_ion"]["58CeII"]["n_fresh_fail"] == 7 and g3["per_ion"]["58CeII"]["k_rec_max"] == 32
    assert g3["per_ion"]["60NdII"]["n_fresh_fail"] == 0 and g3["per_ion"]["57LaII"]["n_fresh_fail"] == 0
    assert g3["per_ion"]["60NdII"]["transfer_holds_axes"] == ["J"] and "P" in g3["per_ion"]["60NdII"]["interpolation_fails"] or g3["per_ion"]["60NdII"]["interp_passes_axes"] == ["T", "D", "J"]
    assert not g3["nd_trajectory"]["matrix_only_passes"] and g3["nd_trajectory"]["interp_whole_band"] > h["g1"]["dm_max"]
    assert g3["part_b"]["k_mix"] is None and g3["part_b"]["k_direct"] == 4 and g3["part_c"]["k_rec"] == 4
    for ion in h["ions"]:
        assert h["r2m"]["per_ion"][ion]["survives"] == "YES"
        assert h["audit"]["per_ion"][ion]["exact"] == 1.0            # every stored exit is a real line
        # the local family never needs more archetypes than the global one
        kl, kg = h["g2"]["per_ion"][ion]["k_local"], h["g2"]["per_ion"][ion]["k_global"]
        assert kl is not None and (kg is None or kg >= kl)
    # the mechanism, on the decisive ions: better on events, worse on the light
    for ion in h["decisive"]:
        c = h["g2"]["contrast"][ion]
        assert c["event_global"] < c["event_local"] and c["band_global"] > c["band_local"]
        assert c["params_global"] > c["params_local"]


def test_every_macro_the_manuscript_uses_is_defined_and_no_literal_numbers():
    chk = _load("paperB_check", "docs/paperB/check_structure.py")
    chk.main()                                                        # exits non-zero on any problem


def test_macros_are_generated_and_ordered_lists_survive_the_json_round_trip():
    h = json.loads((ROOT / "paperB/FROZEN.json").read_text())
    names = [n for n, _, _ in LT.macros(h)]
    assert len(names) == len(set(names)) and "PBNdKlocal" in names
    for ion in h["ions"]:
        for fam in ("local", "glob"):
            ks = [e["k"] for e in h["g2"]["per_ion"][ion][fam]]
            assert ks == sorted(ks) and all(isinstance(k, int) for k in ks)
        rho = [e["rho"] for e in h["g2"]["per_ion"][ion]["trunc"]]
        assert len(rho) == 6


@pytest.mark.parametrize("which", [1, 2, 3])
def test_figures_build_from_the_frozen_record(tmp_path, which):
    fig = _load("paperB_figures", "docs/paperB/figures.py")
    written = fig.main(tmp_path, which=which)
    assert len(written) == 2 and all(p.stat().st_size > 5000 for p in written)
