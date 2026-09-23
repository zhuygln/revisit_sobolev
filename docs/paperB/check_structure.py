#!/usr/bin/env python3
"""Paper B structural checks on the skeleton: the committed fragments are the
regeneration of FROZEN.json, every macro the prose uses is defined, the
figures the manuscript includes exist, and no literal result number has crept
into the prose (a line may carry `% literal-ok` for a genuine constant).

    .venv/bin/python docs/paperB/check_structure.py
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import latex_tables as LT                                        # noqa: E402

TEX = HERE / "manuscript.tex"
# a number with a decimal point, or >= 3 digits: the kind that can only be a result
LITERAL = re.compile(r"(?<![\w\\.])\d+\.\d+|(?<![\w\\.])\d{3,}")


def main():
    problems = []
    for name, text in (("numbers.tex", LT.numbers_tex(LT.json.loads(LT.FROZEN.read_text()))),
                       ("tab_gates.tex", LT.tab_gates(LT.json.loads(LT.FROZEN.read_text()))),
                       ("tab_contrast.tex", LT.tab_contrast(LT.json.loads(LT.FROZEN.read_text())))):
        p = HERE / name
        if not p.exists() or p.read_text() != text:
            problems.append(f"{name} is not the regeneration of paperB/FROZEN.json (run `make tables`)")
    tex = TEX.read_text()
    defined = set(re.findall(r"\\newcommand\{\\(\w+)\}", (HERE / "numbers.tex").read_text()))
    used = set(re.findall(r"\\(PB\w+)", tex)) | set(re.findall(r"\\(PB\w+)", (HERE / "tab_gates.tex").read_text())) \
        | set(re.findall(r"\\(PB\w+)", (HERE / "tab_contrast.tex").read_text()))
    for m in sorted(used - defined):
        problems.append(f"manuscript uses undefined macro \\{m}")
    for inc in re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", tex):
        if not (HERE / inc).exists():
            problems.append(f"missing figure {inc}")
    body = []
    for line in tex.splitlines():
        if line.lstrip().startswith("%") or "literal-ok" in line:
            continue
        body.append(re.sub(r"%.*$", "", line))
    for i, line in enumerate(body, 1):
        stripped = re.sub(r"\$[^$]*\$", "", line)               # inline math is not prose
        for m in LITERAL.findall(stripped):
            problems.append(f"literal number {m!r} in the prose: {line.strip()[:70]}")
    placeholders = tex.count("TO BE WRITTEN") + tex.count("TO BE COMPLETED")
    if problems:
        for p in problems:
            print("FAIL:", p)
        sys.exit(1)
    print(f"structure check OK: {len(defined)} macros defined, {len(used)} used, "
          f"{placeholders} section(s) still to write (skeleton)")


if __name__ == "__main__":
    main()
