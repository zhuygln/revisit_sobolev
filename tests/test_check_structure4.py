"""docs/paper4/check_structure.py: the literal-number ban catches hand-typed
results and ignores macros; the section word count ignores display items."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("paper4_check_structure", ROOT / "docs/paper4/check_structure.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def test_literal_ban_catches_results_and_ignores_macros():
    cs = _load()
    text = ("\\begin{abstract}\n"
            "the line-binned closure within 0.13 mag and the expansion excess 15\\% at 2 of 3 epochs\n"
            "the range 0.13--0.45 mag and 2--4\\%\n"
            "within \\LcBbinthColour\\ mag and \\LcBthDL\\ at 2 d  % literal-ok\n"
            "at $X_{\\rm lan}=10^{-1}$ with \\BtwoRtwoZPone\\ mag; 24 fine shells; a $1$~mag scale\n"
            "\\end{document}\n")
    problems, exempt = cs.check_literals(text)
    hits = " ".join(problems)
    assert "0.13 mag" in hits and "15\\%" in hits and "2 of 3" in hits and "0.13--0" in hits and "2--4\\%" in hits
    assert not any("LcBbinthColour" in p or "BtwoRtwoZPone" in p or "10^{-1}" in p or "24 fine" in p for p in problems)
    assert exempt == [] and not any("line 4" in p for p in problems)


def test_words_ignore_display_items_and_keys():
    cs = _load()
    t = "\\section{Results}\\label{sec:results} one two three \\citep{fontes2020} \\begin{figure}\\includegraphics[width=1in]{x}\\caption{a b c d}\\end{figure} four"
    assert cs.words(t) == 5


def test_claim_check_substitutes_macros_and_flags_undefined_ones(tmp_path):
    """docs/paper4/check_claims.py: the read-through substitutes the frozen
    value into the sentence, and a macro that numbers.tex does not define is
    a problem (it would typeset as nothing)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("paper4_check_claims", ROOT / "docs/paper4/check_claims.py")
    cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
    vals = {"LcBthDL": "15\\%", "PxBtwoIK": "-0.33"}
    out, used = cc.substitute("The expansion closure is \\LcBthDL\\ brighter \\citep{x} and \\PxBtwoIK\\ mag in $i-K$.", vals)
    assert out == "The expansion closure is [15%] brighter and [-0.33] mag in $i-K$."
    assert used == ["LcBthDL", "PxBtwoIK"]
    body = "\\begin{abstract}\nOne \\LcBthDL\\ here.\n\\begin{figure}\\caption{\\LcBthDL}\\end{figure}\nTwo sentences.\n"
    ss = cc.sentences(body)
    assert any("\\LcBthDL" in x for x in ss) and not any("caption" in x for x in ss)
