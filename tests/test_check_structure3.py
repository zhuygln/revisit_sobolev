"""docs/paper3/check_structure.py: the literal-number ban catches hand-typed
results and ignores macros; the budgets are counted on prose only."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS3 = ROOT / "docs" / "paper3"


def _load():
    spec = importlib.util.spec_from_file_location("paper3_check_structure", DOCS3 / "check_structure.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_literal_ban_catches_results_and_ignores_macros():
    cs = _load()
    text = ("\\begin{abstract}\n"
            "class C-B at 27 of 27 points, median 0.83 mag, 18/26 dense\n"
            "the range 0.96--2.84 mag and 1--3 mag\n"
            "class C-B at \\GateTwoCB\\ of \\nPoints\\ points  % literal-ok\n"
            "at $\\Xlan=10^{-3}$ with \\HeadlineColourRange\\ mag; 32 groups; a $1$~mag allowance\n"
            "\\end{document}\n")
    problems, exempt = cs.check_literals(text)
    hits = " ".join(problems)
    assert "27 of 27" in hits and "0.83 mag" in hits and "18/26" in hits
    assert "0.96--2" in hits and "2.84 mag" in hits and "1--3 mag" in hits
    assert not any("GateTwoCB" in p or "HeadlineColourRange" in p for p in problems)
    assert not any("10^{-3}" in p or "32 groups" in p for p in problems)
    assert exempt == []
    assert not any("line 4" in p for p in problems)


def test_literal_ok_is_printed_not_failed():
    cs = _load()
    text = "\\begin{abstract}\nthe depths 23.5 and 21.5 mag  % literal-ok\n\\end{document}\n"
    problems, exempt = cs.check_literals(text)
    assert problems == [] and len(exempt) == 1 and exempt[0][1] == "21.5 mag"


def test_retracted_values():
    cs = _load()
    assert cs.check_retracted("class C-B at 24 of 24 points") != []
    assert cs.check_retracted("the earlier draft quoted 24 of 24, retracted after the mask") == []
    assert cs.check_retracted("the 0.7 mag headline") != []
    assert cs.check_retracted("\\GateTwoCB\\ of \\nPoints") == []


def test_word_count_excludes_display_items_and_keys():
    cs = _load()
    text = ("\\begin{figure}\\caption{one two three four five six}\\end{figure}\n"
            "alpha beta \\cite{a,b,c} gamma \\ref{fig:x} \\GateTwoCB\\ delta")
    assert cs.words(text) == 5   # alpha beta gamma \GateTwoCB delta


def test_committed_manuscript_passes_everything_but_the_todos():
    cs = _load()
    import io, contextlib
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        rc = cs.main()
    problems = [l.strip("- ").strip() for l in err.getvalue().splitlines() if l.startswith("  - ")]
    assert all(p.startswith("unresolved TODO") for p in problems), problems
    assert rc == (1 if problems else 0)


def test_si_is_checked_with_the_same_rules(tmp_path, monkeypatch):
    cs = _load()
    si = tmp_path / "si.tex"
    si.write_text("\\input{si_tab_missing}\nC-B at 27 of 27 points.\n\\todo{x}\n")
    monkeypatch.setattr(cs, "SI", si)
    problems, exempt = cs.check_si()
    assert all(p.startswith("si.tex: ") for p in problems)
    assert any("si_tab_missing" in p for p in problems)
    assert any("27 of 27" in p for p in problems)
    assert any("TODO" in p for p in problems)
    assert exempt == []
    monkeypatch.setattr(cs, "SI", tmp_path / "absent.tex")
    assert cs.check_si() == ([], [])
