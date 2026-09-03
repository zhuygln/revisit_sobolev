"""docs/paper3/latex_tables.py: the number macros and table fragments are a
pure function of FROZEN.headline, every macro points at a real headline key,
and the committed numbers.tex is that function's output."""
import importlib.util, json, sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS3 = ROOT / "docs" / "paper3"
FROZEN = ROOT / "paper3" / "FROZEN.json"
needs_frozen = pytest.mark.skipif(not FROZEN.exists(), reason="FROZEN.json not present")


def _load():
    spec = importlib.util.spec_from_file_location("paper3_latex_tables", DOCS3 / "latex_tables.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _headline():
    return json.loads(FROZEN.read_text())["headline"]


def test_macro_names_unique_and_letters_only():
    lt = _load()
    names = [m[0] for m in lt.MACROS]
    # a "byx" entry shares its base name with the total's entry and expands to three suffixed macros
    byx = sum(m[2] == "byx" for m in lt.MACROS)
    assert len(set(names)) == len(names) - byx
    assert all(n.isalpha() for n in names), "LaTeX macro names cannot contain digits"
    assert all(fmt in lt.FMT or fmt == "byx" for _, _, fmt, _ in lt.MACROS)


@needs_frozen
def test_every_macro_key_exists_in_frozen():
    lt = _load()
    h = _headline()
    for name, key, fmt, _ in lt.MACROS:
        assert key in h, f"{name}: headline key {key} missing"
    out = lt.macros(h)
    byx = sum(m[2] == "byx" for m in lt.MACROS)
    assert len(out) == len(lt.MACROS) + 2 * byx
    assert len({n for n, _ in out}) == len(out)


@needs_frozen
def test_numbers_tex_is_the_regeneration():
    lt = _load()
    assert (DOCS3 / "numbers.tex").read_text() == lt.numbers_tex(_headline())


@needs_frozen
def test_regeneration_is_byte_identical(tmp_path):
    sys.path.insert(0, str(ROOT / "paper3"))
    import freeze
    lt = _load()
    h = _headline()
    a = lt.main(h, freeze.canonical(), tmp_path)
    first = {Path(p).name: Path(p).read_bytes() for p in a}
    b = lt.main(h, freeze.canonical(), tmp_path)
    assert {Path(p).name: Path(p).read_bytes() for p in b} == first
    for name in ("numbers.tex", "tab_verdict.tex", "tab_grid.tex", "tab_scenarios.tex"):
        assert first[name] == (DOCS3 / name).read_bytes(), f"committed {name} is stale"


@needs_frozen
def test_verdict_table_has_three_X_rows():
    sys.path.insert(0, str(ROOT / "paper3"))
    import freeze
    lt = _load()
    rows = lt.verdict_rows(_headline(), freeze.canonical())
    assert [r["X"] for r in rows] == ["$10^{-3}$", "$10^{-2}$", "$10^{-1}$"]
    tex = lt.tab_verdict(_headline(), freeze.canonical())
    assert tex.count(r"\\") == 2 + 3   # two header rows + three body rows


def test_formatters():
    lt = _load()
    assert lt.FMT["rngwhole"]([0.96, 2.84]) == "1--3"
    assert lt.FMT["of"]([27, 27]) == "27 of 27"
    assert lt.FMT["slash"]({"survives": 18, "eligible": 26}) == "18/26"
    assert lt.FMT["pow"](100000) == "$10^{5}$"
    assert lt.FMT["pct"]([292, 524]) == r"56\%"
    assert lt.FMT["rng2"]([-2.123, -0.239]).startswith("$-2.12$")


def test_column_check_rejects_ragged_rows():
    lt = _load()
    bad = "\\begin{tabular}{lcc}\n\\toprule\na & b & c & d \\\\\n\\bottomrule\n\\end{tabular}\n"
    with pytest.raises(AssertionError):
        lt._check_columns(bad)
