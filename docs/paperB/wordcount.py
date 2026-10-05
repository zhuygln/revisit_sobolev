#!/usr/bin/env python3
"""Approximate PRL core word count: the abstract and body up to the End
Matter heading, comments and display items removed, macros counted as one
word each, with the APS figure/table word-equivalent rule applied
separately (150 words + 20 per cm of height for a single-column figure;
here every figure is counted at the two-column, full-width estimate the
Makefile layout produces). Reported, not enforced: the Makefile prints it.

    .venv/bin/python docs/paperB/wordcount.py
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE / "manuscript.tex"


def core_text(raw):
    s = re.sub(r"(?m)%.*$", "", raw)
    s = s[s.index("\\begin{abstract}"):]
    end = s.find("\\section*{End Matter}")
    if end >= 0:
        s = s[:end]
    s = re.sub(r"\\begin\{(figure|table)\*?\}.*?\\end\{\1\*?\}", " ", s, flags=re.S)
    s = re.sub(r"\\(cite[pt]?|ref|label|input|includegraphics)(\[[^\]]*\])?\{[^}]*\}", " ", s)
    s = re.sub(r"\\(begin|end)\{[^}]*\}", " ", s)
    s = re.sub(r"\\section\*?\{([^}]*)\}", r" \1 ", s)
    s = re.sub(r"\\paragraph\{([^}]*)\}", r" \1 ", s)
    s = re.sub(r"\$[^$]*\$", " EQN ", s)
    s = re.sub(r"\\[A-Za-z]+\*?", " M ", s)
    s = re.sub(r"[{}~\\]", " ", s)
    return s


def main():
    raw = TEX.read_text()
    words = len(re.findall(r"\S+", core_text(raw)))
    n_fig = len(re.findall(r"\\begin\{figure", raw[:raw.find("\\section*{End Matter}")] if "\\section*{End Matter}" in raw else raw))
    n_tab = len(re.findall(r"\\begin\{table", raw[:raw.find("\\section*{End Matter}")] if "\\section*{End Matter}" in raw else raw))
    n_disp = len(re.findall(r"\\begin\{equation|\\\[", raw))
    print(f"core text words (macros count one each): {words}")
    print(f"core figures: {n_fig}; core tables: {n_tab}; display equations: {n_disp}")
    print("APS: core limit 3750 words incl. figure/table equivalents (~150 + 20/cm each at single column; "
          "a full-width figure counts roughly double); End Matter up to two pages, not counted.")


if __name__ == "__main__":
    main()
