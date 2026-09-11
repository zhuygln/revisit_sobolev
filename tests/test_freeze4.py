"""paper4/freeze.py and docs/paper4/latex_tables.py: the frozen headline is
the regeneration from the committed records, the macros cover every quoted
name once, and the fragments are byte-stable."""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "paper4" / "FROZEN.json"
needs_frozen = pytest.mark.skipif(not FROZEN.exists(), reason="FROZEN.json not present")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


@needs_frozen
def test_frozen_headline_regenerates_from_the_records():
    fz = _load("paper4_freeze", ROOT / "paper4/freeze.py")
    same = json.dumps(json.loads(FROZEN.read_text())["headline"], sort_keys=True, default=float) == json.dumps(fz.headline(), sort_keys=True, default=float)
    assert same, "FROZEN.json stale: run paper4/freeze.py"   # NaN-tolerant (json)


@needs_frozen
def test_macros_are_unique_and_the_fragments_are_committed():
    lt = _load("paper4_latex_tables", ROOT / "docs/paper4/latex_tables.py")
    h = json.loads(FROZEN.read_text())["headline"]
    m = lt.macros(h)
    names = [n for n, _ in m]
    assert len(names) == len(set(names)) and len(names) >= 40
    assert set(lt.quoted_names(h)) <= set(names)
    for name, fn in (("numbers.tex", lt.numbers_tex), ("tab_lightcurve.tex", lt.tab_lightcurve),
                     ("tab_grid.tex", lt.tab_grid), ("tab_matrix.tex", lt.tab_matrix)):
        assert (ROOT / "docs/paper4" / name).read_text() == fn(h), f"{name} stale: run make tables"
    tex = lt.numbers_tex(h)
    assert tex.count("\\newcommand") == len(names) and "do not edit" in tex


def test_formats_are_ranges_sorted_by_magnitude():
    lt = _load("paper4_latex_tables", ROOT / "docs/paper4/latex_tables.py")
    assert lt.rngabs2([-0.194, -0.146]) == "0.15--0.19"
    assert lt.pct0(0.153) == "15\\%" and lt.sf2(-0.1704) == "-0.17" and lt.sci(4.477e46) == "4.48 \\times 10^{46}"
