#!/usr/bin/env python3
"""Fail the build if the manuscript has silently lost content.

Written after an appendix was destroyed by a programmatic edit whose search
string was not unique: the entire worldline derivation vanished and the page
count did not move, because a newly added appendix replaced it almost
byte-for-byte. LaTeX compiled cleanly. Nothing flagged it.

Page count is not a content check. This is.

Checks, in order of how they would have caught that failure:
  1. every expected \\label is present;
  2. every section and subsection exceeds a minimum word count, so a section
     reduced to a stub fails;
  3. every \\includegraphics target exists on disk;
  4. no \\todo remains outside the macro definition.
"""

import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
TEX = HERE / "manuscript.tex"

REQUIRED_LABELS = [
    "sec:intro", "sec:primer", "sec:primer-expop", "sec:methods",
    "sec:conventions", "sec:popcontrol", "sec:conventions-f8", "sec:data",
    "sec:verification", "sec:results", "sec:crowding", "sec:forest",
    "sec:separation", "sec:map", "sec:binwidth", "sec:breadth",
    "sec:multiion", "sec:cost",
    "app:formal", "app:sobolev", "app:symmetry", "app:expop", "app:pavg",
    "app:boundary", "app:frame",
    "tab:findings",
]
# The destroyed appendix left ~30 words; the shortest legitimate subsection
# here is ~82. 60 separates them with margin at both ends.
MIN_WORDS = 60
FIG_DIR = HERE.parent / "figures"


def main():
    text = TEX.read_text()
    problems = []

    # 1. labels
    present = set(re.findall(r"\\label\{([^}]+)\}", text))
    for lab in REQUIRED_LABELS:
        if lab not in present:
            problems.append(f"missing \\label{{{lab}}}")

    # 2. section lengths -- the check that would have caught the loss.
    # A heading's body runs to the next heading of the SAME OR HIGHER level,
    # so a \section that merely introduces subsections is not counted as empty.
    heads = [
        (m.start(), 1 if m.group(1) is None else 2, m.group(2))
        for m in re.finditer(r"\\(sub)?section\*?\{([^}]*)\}", text)
    ]
    for i, (pos, level, title) in enumerate(heads):
        end = len(text)
        for pos2, level2, _ in heads[i + 1:]:
            if level2 <= level:
                end = pos2
                break
        if title.lower().startswith(("reproducib", "derivations")):
            continue
        words = len(text[pos:end].split())
        if words < MIN_WORDS:
            problems.append(
                f"{'sub' if level == 2 else ''}section '{title}' has only "
                f"{words} words (< {MIN_WORDS}) -- content loss?"
            )

    # 3. figures exist
    for fig in re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", text):
        if not (FIG_DIR / fig).exists():
            problems.append(f"missing figure {fig}")

    # 4. no leftover TODOs
    todos = [l for l in text.splitlines()
             if "\\todo{" in l and "newcommand" not in l]
    for t in todos:
        problems.append(f"unresolved TODO: {t.strip()[:70]}")

    if problems:
        print("MANUSCRIPT STRUCTURE CHECK FAILED:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    print(f"structure check OK: {len(REQUIRED_LABELS)} labels, "
          f"all sections >= {MIN_WORDS} words, figures present, no TODOs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
