#!/usr/bin/env python3
"""The book's no-hand-typed-results rule, scoped as the PI asked (amendment 3):
pedagogical constants (10^5 packets, three groups, 0.1c, example optical
depths) are free everywhere; inside a `::: {.result}` ... `:::` region, or on
a line carrying `<!-- result -->`, every number in prose must come from an inline
`{python}` expression (formulas in `$...$` and code spans are not prose). A line ending in `<!-- literal-ok -->` is exempt and
printed. Exits non-zero on a violation."""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
CHAPTERS = sorted((HERE / "chapters").glob("*.qmd")) + [HERE / "index.qmd"]
NUMBER = re.compile(r"(?<![A-Za-z_\\{])[-+]?\d+(?:\.\d+)?(?:\s*(?:%|\\%|mag))?(?![\d.]*[A-Za-z_}])")


def strip_code(line):
    line = re.sub(r"`[^`]*`", "", line)            # inline code and {python} expressions
    return re.sub(r"\$[^$]*\$", "", line)          # inline math: a formula is not a result


def main():
    problems, exempt = [], []
    for path in CHAPTERS:
        in_result = False
        for n, raw in enumerate(path.read_text().splitlines(), 1):
            s = raw.strip()
            if re.match(r"^:::\s*\{\.result\}", s):
                in_result = True; continue
            if in_result and s == ":::":
                in_result = False; continue
            if not (in_result or "<!-- result -->" in raw):
                continue
            line = strip_code(raw)
            line = re.sub(r"<!--.*?-->", "", line)
            for m in NUMBER.finditer(line):
                if raw.rstrip().endswith("<!-- literal-ok -->"):
                    exempt.append((path.name, n, m.group(0)))
                else:
                    problems.append(f"{path.name}:{n}: literal '{m.group(0)}' in a result region -- quote it through `{{python}}`")
    for name, n, s in exempt:
        print(f"exempt: {name}:{n}: {s!r}")
    if problems:
        print("NO-LITERALS CHECK FAILED:", file=sys.stderr)
        for p in problems:
            print("  - " + p, file=sys.stderr)
        return 1
    print(f"no-literals check OK ({len(CHAPTERS)} files, {len(exempt)} exemptions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
