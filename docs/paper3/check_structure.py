#!/usr/bin/env python3
"""Fail the Paper III build if the manuscript has lost content, quotes a
number that is not a macro, or exceeds the journal's budgets.

Inherits Paper I's checks (labels, minimum section length, figures on disk, no
\\todo outside its definition, retracted numbers) and adds the ones this paper
needs:

  * numbers.tex must be byte-identical to latex_tables.numbers_tex(FROZEN) --
    the file is committed, so a stale copy would compile silently;
  * every macro flagged `quoted` must be used at least once (manuscript or SI);
  * no literal result number in the prose: "k of n", "k/n", decimals with
    mag or %, and decimal ranges are banned outside numbers.tex unless the
    line ends with `% literal-ok`, in which case it is printed for review;
  * Nature Astronomy budgets: abstract <= 150 words, main text <= 4000,
    Methods <= 3500, display items <= 6, unique main-text citations <= 50.
"""

import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
TEX = HERE / "manuscript.tex"
SI = HERE / "si.tex"
FIG_DIR = HERE / "figures"

REQUIRED_LABELS = [
    "sec:intro", "sec:experiment", "sec:asymmetry", "sec:bias",
    "sec:notdegeneracy", "sec:nuisance", "sec:observable", "sec:allowance",
    "sec:discussion", "sec:methods", "sec:ed",
    "fig:one", "fig:grid", "fig:notparam", "fig:observe", "tab:verdict",
    "edfig:source", "edfig:tscale", "edfig:chain", "edfig:control",
    "edtab:grid", "edtab:scenarios",
]
MIN_WORDS = 60
BUDGET = {"abstract": 150, "main": 4000, "methods": 3500,
          "display_items": 6, "main_cites": 50}

# Literal numbers that look like results. Definitions (grid values, group
# counts, thresholds, depths) also match some of these; those lines carry
# `% literal-ok` and are printed so that the exemption stays visible.
LITERAL = [
    (r"(?<![\d.^{-])\d+\s*(?:of|/)\s*\d+(?![\d}])", "k of n / k/n"),
    (r"\d+\.\d+\$?\s*~?\\?(?:mag|\\%)", "decimal with mag or %"),
    (r"\d+\.\d+\s*(?:--|-)\s*\$?\d", "decimal range"),
    (r"\d+\s*--\s*\d+\s*~?mag", "integer range in mag"),
]
LITERAL_OK = "% literal-ok"

# Values the project has retracted (see docs/results_report.md 4.41-4.44 and
# the erratum trail in docs/lab_notebook.md). None may reappear.
RETRACTED = [
    (r"166 of 170|166/170", "the pre-mask 166 of 170 count"),
    (r"24 of 24|24/24", "the pre-mask 24/24 Gate-2 count"),
    (r"0\.7\s*~?mag", "the retracted 0.7 mag headline amplitude"),
    (r"40\s*\\%\s*hotter|\+0\.33", "the retracted 'T_eff 40% hotter' reading"),
    (r"18/18 dense|18 of 18 dense", "the pre-redo dense Gate-3 count as current"),
]
RETRACT_EXEMPT = ("retract", "earlier version", "earlier draft", "erratum")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def strip_comments(text):
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", l) for l in text.splitlines())


def words(text):
    """Word count of LaTeX prose: environments that are display items are
    removed, citation and reference keys are not words, macros are."""
    t = re.sub(r"\\begin\{(figure|table)\}.*?\\end\{\1\}", " ", text, flags=re.S)
    t = re.sub(r"\\(cite|ref|label|input|includegraphics)\*?(\[[^\]]*\])?\{[^}]*\}", " ", t)
    t = re.sub(r"\\(begin|end)\{[^}]*\}", " ", t)
    t = re.sub(r"\\(item|centering|small|footnotesize|scriptsize|clearpage|addcontentsline\{[^}]*\}\{[^}]*\})", " ", t)
    t = re.sub(r"\\section\*?\{([^}]*)\}", r" \1 ", t)
    return len([w for w in t.split() if re.search(r"[A-Za-z0-9]", w)])


def region(text, start_pat, end_pat):
    a = re.search(start_pat, text)
    b = re.search(end_pat, text)
    if not a or not b or b.start() < a.start():
        return None
    return text[a.start():b.start()]


def check_numbers_tex():
    freeze = _load("paper3_freeze", ROOT / "paper3" / "freeze.py")
    lt = _load("latex_tables", HERE / "latex_tables.py")
    if not freeze.FROZEN.exists():
        return [f"missing {freeze.FROZEN} -- run paper3/freeze.py"], lt, None
    h = json.loads(freeze.FROZEN.read_text())["headline"]
    want = lt.numbers_tex(h)
    have = (HERE / "numbers.tex").read_text()
    out = []
    if want != have:
        out.append("numbers.tex is not the regeneration from FROZEN.json -- "
                   "run `make tables`")
    return out, lt, h


def check_macro_usage(lt, used_text):
    out = []
    for name in lt.quoted_names():
        if not re.search(r"\\" + name + r"(?![A-Za-z])", used_text):
            out.append(f"macro \\{name} is generated but never used")
    return out


def check_literals(text):
    out, exempt = [], []
    body = region(text, r"\\begin\{abstract\}", r"\\end\{document\}") or text
    offset = text.find(body)
    for i, raw in enumerate(body.splitlines()):
        line_no = text[:offset].count("\n") + i + 1
        if raw.lstrip().startswith("%"):
            continue
        ok = raw.rstrip().endswith(LITERAL_OK)
        line = re.sub(r"(?<!\\)%.*$", "", raw)
        line = re.sub(r"\\(cite|ref|label)\{[^}]*\}", "", line)
        for pat, why in LITERAL:
            for m in re.finditer(pat, line):
                if ok:
                    exempt.append((line_no, m.group(0), why))
                else:
                    out.append(f"line {line_no}: literal number '{m.group(0)}' "
                               f"({why}) -- use a macro from numbers.tex or "
                               f"mark the line `{LITERAL_OK}`")
    return out, exempt


def check_retracted(text):
    out = []
    for pat, why in RETRACTED:
        for m in re.finditer(pat, text):
            near = text[max(0, m.start() - 400): m.end() + 400].lower()
            if any(k in near for k in RETRACT_EXEMPT):
                continue
            out.append(f"line {text[:m.start()].count(chr(10)) + 1}: "
                       f"retracted value reappears -- {why}")
    return out


def check_budgets(text):
    out = []
    counts = {}
    ab = region(text, r"\\begin\{abstract\}", r"\\end\{abstract\}")
    counts["abstract"] = words(ab) if ab else 0
    main = region(text, r"\\section\{Introduction\}", r"\\section\*\{Methods\}")
    counts["main"] = words(main) if main else 0
    meth = region(text, r"\\section\*\{Methods\}", r"\\section\*\{Extended Data\}")
    counts["methods"] = words(meth) if meth else 0
    counts["display_items"] = len(re.findall(r"\\begin\{(?:figure|table)\}", main or ""))
    keys = set()
    for m in re.finditer(r"\\cite\{([^}]*)\}", main or ""):
        keys.update(k.strip() for k in m.group(1).split(","))
    counts["main_cites"] = len(keys)
    for k, lim in BUDGET.items():
        if counts[k] > lim:
            out.append(f"{k}: {counts[k]} exceeds the budget of {lim}")
    return out, counts


def main():
    raw = TEX.read_text()
    text = strip_comments(raw)
    problems = []

    present = set(re.findall(r"\\label\{([^}]+)\}", text))
    for lab in REQUIRED_LABELS:
        if lab not in present:
            problems.append(f"missing \\label{{{lab}}}")

    heads = [(m.start(), 1 if m.group(1) is None else 2, m.group(2))
             for m in re.finditer(r"\\(sub)?section\*?\{([^}]*)\}", text)]
    for i, (pos, level, title) in enumerate(heads):
        end = len(text)
        for pos2, level2, _ in heads[i + 1:]:
            if level2 <= level:
                end = pos2
                break
        if title.lower().startswith(("extended data", "data availability",
                                     "code availability")):
            continue
        n = words(text[pos:end])
        if n < MIN_WORDS:
            problems.append(f"{'sub' if level == 2 else ''}section '{title}' "
                            f"has only {n} words (< {MIN_WORDS}) -- content loss?")

    for fig in re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", text):
        if not ((FIG_DIR / fig).exists() or any((FIG_DIR / (fig + e)).exists()
                                                 for e in (".pdf", ".png"))):
            problems.append(f"missing figure {fig} (no .pdf or .png)")
    for frag in re.findall(r"\\input\{([^}]+)\}", text):
        if not (HERE / (frag + ".tex")).exists():
            problems.append(f"missing \\input fragment {frag}.tex")

    for l in raw.splitlines():
        if "\\todo{" in l and "newcommand" not in l:
            problems.append(f"unresolved TODO: {l.strip()[:70]}")

    num_problems, lt, h = check_numbers_tex()
    problems.extend(num_problems)
    used = raw + (SI.read_text() if SI.exists() else "")
    problems.extend(check_macro_usage(lt, used))

    lit, exempt = check_literals(raw)
    problems.extend(lit)
    problems.extend(check_retracted(text))
    budget_problems, counts = check_budgets(text)
    problems.extend(budget_problems)

    if exempt:
        print(f"{len(exempt)} literal(s) exempted with `{LITERAL_OK}`:")
        for line_no, s, why in exempt:
            print(f"  line {line_no}: {s!r} ({why})")
    print("budgets: " + ", ".join(f"{k} {counts[k]}/{BUDGET[k]}" for k in BUDGET))
    if problems:
        print("MANUSCRIPT STRUCTURE CHECK FAILED:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    print(f"structure check OK: {len(REQUIRED_LABELS)} labels, all sections >= "
          f"{MIN_WORDS} words, figures present, numbers.tex frozen, "
          f"{len(lt.quoted_names())} macros used, no TODOs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
