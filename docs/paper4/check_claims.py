#!/usr/bin/env python3
"""Pre-submission claim check: print every sentence of the manuscript that
quotes a number, with the macro replaced by the value the freeze gives it,
so that each claim can be read against the record in the words a referee
will read. A macro naming the wrong quantity, or a sign convention slipped
between the prose and the table, is invisible in the LaTeX source and
obvious here.

    make claims            # the read-through
    make claims ARGS=-q    # assertions only (undefined or unused macros)

Exits non-zero when the prose uses a macro that numbers.tex does not define
(a silent empty value in the PDF). The PI's standing rule after the table
incident of 2026-09-11: every number in a table or a figure annotation is
machine-generated from paper4/FROZEN.json, none is transcribed by hand, and
this read-through runs before submission."""
import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEX = HERE / "manuscript.tex"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def sentences(text):
    """The abstract and body, split into sentences, display items dropped."""
    body = text[text.index("\\begin{abstract}"):]
    body = re.sub(r"%.*", "", body)
    body = re.sub(r"\\begin\{(figure|table)\*?\}.*?\\end\{\1\*?\}", " ", body, flags=re.S)
    return [re.sub(r"\s+", " ", s).strip() for s in re.split(r"(?<=[.!?])\s+", body)]


def substitute(s, vals):
    used = [m for m in re.findall(r"\\([A-Z][A-Za-z]+)(?![A-Za-z])", s) if m in vals]
    for m in used:
        s = re.sub(r"\\" + m + r"\\?(?![A-Za-z])", "[" + vals[m].replace("\\%", "%") + "]", s)
    s = re.sub(r"\\(cite[pt]?|ref|label)\{[^}]*\}", "", s)
    return re.sub(r"\s+", " ", s).strip(), used


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-q", "--quiet", action="store_true", help="assertions only")
    a = ap.parse_args()
    lt = _load("paper4_latex_tables", HERE / "latex_tables.py")
    h = json.loads((ROOT / "paper4" / "FROZEN.json").read_text())["headline"]
    vals = dict(lt.macros(h))
    raw = TEX.read_text()
    prose = re.sub(r"%.*", "", raw)
    prose = prose[prose.index("\\begin{abstract}"):]

    problems = []
    defined = set(vals)
    # every macro the prose uses that looks like one of ours but is not defined
    declared = set(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}", raw))
    for m in set(re.findall(r"\\([A-Z][A-Za-z]+)(?![A-Za-z])", prose)):
        if m in defined or m in declared or m in ("Pone", "Ptwo"):
            continue
        if m in ("Rres", "Bexp", "Bbin", "Msun", "Xlan", "Nd", "NdII", "NdIII"):
            continue
        problems.append(f"\\{m} is used in the prose but is neither a number macro nor declared in the preamble")
    n = 0
    for s in sentences(raw):
        out, used = substitute(s, vals)
        if not used:
            continue
        n += 1
        if not a.quiet:
            print(f"{n:3d}. {out}")
    if not a.quiet:
        print()
    print(f"claim check: {n} sentences quote {len(set(sum([substitute(s, vals)[1] for s in sentences(raw)], [])))} "
          f"of {len(vals)} macros; {len(problems)} problem(s)")
    for p in problems:
        print(f"  - {p}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
