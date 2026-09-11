#!/usr/bin/env python3
"""Structure check of the Paper IV manuscript (MNRAS): the required labels,
no thin sections, figures and fragments on disk, no TODOs, numbers.tex
byte-equal to the regeneration from paper4/FROZEN.json, every quoted macro
used, and no hand-typed result number in the prose. MNRAS has no hard word
budget; the counts are printed, not enforced."""
import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEX = HERE / "manuscript.tex"
FIG_DIR = HERE / "figures"

REQUIRED_LABELS = [
    "sec:intro", "sec:method", "sec:states", "sec:results", "sec:closures", "sec:redistribution",
    "sec:snapshot", "sec:lightcurve", "sec:p1xkn", "sec:macroatom", "sec:discussion", "sec:conclusions",
    "fig:matrix", "fig:snapshot", "fig:lightcurve", "fig:p1xkn", "fig:macroatom",
    "tab:matrix", "tab:lightcurve", "tab:grid", "tab:p1xkn",
]
MIN_WORDS = 60

# Literal numbers that look like results. Definitions (grid values, epochs,
# thresholds) also match some of these; those lines carry `% literal-ok`
# and are printed so that the exemption stays visible.
LITERAL = [
    (r"(?<![\d.^{-])\d+\s*(?:of|/)\s*\d+(?![\d}])", "k of n / k/n"),
    (r"\d+\.\d+\$?\s*~?\\?(?:mag|\\%)", "decimal with mag or %"),
    (r"\d+\.\d+\s*(?:--|-)\s*\$?\d", "decimal range"),
    (r"\d+\s*--\s*\d+\s*~?(?:mag|\\%)", "integer range in mag or %"),
    (r"(?<![\d.])\d+\s*\\%", "integer percentage"),
]
LITERAL_OK = "% literal-ok"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def strip_comments(text):
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", l) for l in text.splitlines())


def words(text):
    t = re.sub(r"\\begin\{(figure|table)\*?\}.*?\\end\{\1\*?\}", " ", text, flags=re.S)
    t = re.sub(r"\\(cite[pt]?|ref|label|input|includegraphics)\*?(\[[^\]]*\])?\{[^}]*\}", " ", t)
    t = re.sub(r"\\(begin|end)\{[^}]*\}", " ", t)
    t = re.sub(r"\\(item|centering|small|footnotesize|scriptsize|clearpage)", " ", t)
    t = re.sub(r"\\(sub)*section\*?\{([^}]*)\}", r" \2 ", t)
    return len([w for w in t.split() if re.search(r"[A-Za-z0-9]", w)])


def region(text, start_pat, end_pat):
    a = re.search(start_pat, text); b = re.search(end_pat, text)
    if not a or not b or b.start() < a.start():
        return None
    return text[a.start():b.start()]


def check_numbers_tex():
    freeze = _load("paper4_freeze", ROOT / "paper4" / "freeze.py")
    lt = _load("paper4_latex_tables", HERE / "latex_tables.py")
    if not freeze.OUT.exists():
        return [f"missing {freeze.OUT} -- run paper4/freeze.py"], lt, None
    h = json.loads(freeze.OUT.read_text())["headline"]
    out = []
    for name, fn in (("numbers.tex", lt.numbers_tex), ("tab_lightcurve.tex", lt.tab_lightcurve),
                     ("tab_grid.tex", lt.tab_grid), ("tab_matrix.tex", lt.tab_matrix), ("tab_p1xkn.tex", lt.tab_p1xkn)):
        p = HERE / name
        if not p.exists() or fn(h) != p.read_text():
            out.append(f"{name} is not the regeneration from FROZEN.json -- run `make tables`")
    return out, lt, h


def check_macro_usage(lt, h, used_text):
    return [f"macro \\{name} is generated but never used" for name in lt.quoted_names(h)
            if not re.search(r"\\" + name + r"(?![A-Za-z])", used_text)]


def check_literals(text):
    out, exempt = [], []
    body = region(text, r"\\begin\{abstract\}", r"\\end\{document\}") or text
    offset = text.find(body)
    for i, raw in enumerate(body.splitlines()):
        line_no = text[:offset].count("\n") + i + 1
        if raw.lstrip().startswith("%"):
            continue
        ok = raw.rstrip().endswith(LITERAL_OK)
        if LITERAL_OK in raw and not ok:
            out.append(f"line {line_no}: `{LITERAL_OK}` must end the line -- text after it is a comment")
        line = re.sub(r"(?<!\\)%.*$", "", raw)
        line = re.sub(r"\\(cite[pt]?|ref|label)\{[^}]*\}", "", line)
        for pat, why in LITERAL:
            for m in re.finditer(pat, line):
                if ok:
                    exempt.append((line_no, m.group(0), why))
                else:
                    out.append(f"line {line_no}: literal number '{m.group(0)}' ({why}) -- use a macro "
                               f"from numbers.tex or mark the line `{LITERAL_OK}`")
    return out, exempt


def counts_of(text):
    c = {}
    ab = region(text, r"\\begin\{abstract\}", r"\\end\{abstract\}")
    c["abstract"] = words(ab) if ab else 0
    main = region(text, r"\\section\{Introduction\}", r"\\section\*?\{(?:Data availability|Acknowledgements)\}")
    c["main"] = words(main) if main else 0
    c["display_items"] = len(re.findall(r"\\begin\{(?:figure|table)\*?\}", main or ""))
    keys = set()
    for m in re.finditer(r"\\cite[pt]?\{([^}]*)\}", main or ""):
        keys.update(k.strip() for k in m.group(1).split(","))
    c["cites"] = len(keys)
    return c


def main():
    if not TEX.exists():
        print(f"missing {TEX}", file=sys.stderr); return 1
    raw = TEX.read_text()
    text = strip_comments(raw)
    problems = []
    present = set(re.findall(r"\\label\{([^}]+)\}", text))
    problems += [f"missing \\label{{{lab}}}" for lab in REQUIRED_LABELS if lab not in present]
    heads = [(m.start(), 1 if m.group(1) is None else 2, m.group(2))
             for m in re.finditer(r"\\(sub)?section\*?\{([^}]*)\}", text)]
    for i, (pos, level, title) in enumerate(heads):
        end = len(text)
        for pos2, level2, _ in heads[i + 1:]:
            if level2 <= level:
                end = pos2; break
        if title.lower().startswith(("data availability", "acknowledgements")):
            continue
        n = words(text[pos:end])
        if n < MIN_WORDS:
            problems.append(f"{'sub' if level == 2 else ''}section '{title}' has only {n} words (< {MIN_WORDS})")
    for fig in re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", text):
        if not ((FIG_DIR / fig).exists() or any((FIG_DIR / (fig + e)).exists() for e in (".pdf", ".png"))):
            problems.append(f"missing figure {fig} (no .pdf or .png)")
    for frag in re.findall(r"\\input\{([^}]+)\}", text):
        if not (HERE / (frag + ".tex")).exists():
            problems.append(f"missing \\input fragment {frag}.tex")
    for l in raw.splitlines():
        if "\\todo{" in l and "newcommand" not in l:
            problems.append(f"unresolved TODO: {l.strip()[:70]}")
    num_problems, lt, h = check_numbers_tex()
    problems += num_problems
    if h is not None:
        problems += check_macro_usage(lt, h, raw)
    lit, exempt = check_literals(raw)
    problems += lit
    if exempt:
        print(f"{len(exempt)} literal(s) exempted with `{LITERAL_OK}`:")
        for line_no, s, why in exempt:
            print(f"  line {line_no}: {s!r} ({why})")
    c = counts_of(text)
    print("counts: " + ", ".join(f"{k} {v}" for k, v in c.items()))
    if problems:
        print("MANUSCRIPT STRUCTURE CHECK FAILED:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    print(f"structure check OK: {len(REQUIRED_LABELS)} labels, all sections >= {MIN_WORDS} words, figures present, "
          f"numbers.tex frozen, {len(lt.quoted_names(h))} macros used, no TODOs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
