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


SEPARATION_JSON = (
    HERE.parents[1] / "experiments/sobolev_proper/separation_results.json"
)


def check_headline_numbers(text):
    """The La II headline must agree with the file that produced it."""
    import json

    if not SEPARATION_JSON.exists():
        return [f"missing {SEPARATION_JSON.name} -- cannot verify numbers"]
    h = json.loads(SEPARATION_JSON.read_text())["headline"]
    d_sob = (h["sob"] - h["res"]) / h["res"]
    d_exp = (h["exp"] - h["res"]) / h["res"]

    out = []
    # Anchor on the headline sentence and search only inside it -- a loose
    # regex over the whole document happily matches a similar-looking number
    # from an unrelated section.
    anchor = re.search(
        r"For the La\\,\\textsc\{ii\} forest, against .{0,400}?"
        r"Across the sweep",
        text, re.S,
    )
    if not anchor:
        return ["headline sentence not found -- has it been reworded?"]
    span = anchor.group(0)

    nums = [float(x) for x in re.findall(r"\$=?([\d.]+)\$", span)]
    pcts = [float(x) for x in re.findall(r"\$([+-][\d.]+)\\%\$", span)]
    for label, value in (("resolved", h["res"]), ("Sobolev", h["sob"]),
                         ("expansion", h["exp"])):
        if not any(abs(n - value) < 5e-5 for n in nums):
            out.append(
                f"headline {label} = {value:.4f} (from "
                f"{SEPARATION_JSON.name}) not quoted; sentence has {nums}"
            )
    for label, value in (("Sobolev", 100 * d_sob), ("expansion", 100 * d_exp)):
        if not any(abs(p - value) < 0.05 for p in pcts):
            out.append(
                f"headline Delta {label} = {value:+.1f}% not quoted; "
                f"sentence has {pcts}"
            )
    return out


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

    # 5. headline numbers must match the file that generated them. Every
    # quoted result here has been wrong at least once from hand-copying, so
    # the numbers the abstract and separation section rest on are checked
    # against experiments/sobolev_proper/separation_results.json directly.
    problems.extend(check_headline_numbers(text))

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
