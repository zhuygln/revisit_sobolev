"""Emit the manuscript's single-line ladder table (tab:ladder) from
ladder_summary.json, so the numbers in the paper are written by the file that
produced them. Prints LaTeX to stdout; --insert splices it into manuscript.tex
at the marker line `%%LADDER_TABLE` (or replaces an existing table with that
label)."""
import json, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = 0.1371
rows = {r["tag"]: r for r in json.loads((HERE / "ladder_summary.json").read_text())}


def cell(tag):
    r = rows.get(tag)
    if r is None:
        return "--"
    sem = f"\\pm{r['sem']:.4f}" if r["sem"] else ""
    return f"${r['mean']:.4f}{sem}$"


def pct(tag):
    r = rows.get(tag)
    return "--" if r is None or r["vs_frozen"] is None else f"${100*r['vs_frozen']:+.1f}\\%$"


lines = [
    r"\begin{table}",
    r"\centering",
    r"\small",
    r"\caption{Single-line benchmark: \textsc{sedona} resolved trough depth",
    r"($1400$--$2600\kms$ window) along one axis at a time, fixed seeds, three",
    r"seeds per rung (five at the anchor). The frozen-snapshot target is",
    r"$0.1371$; the production run sat at $0.1420$ with $2\times10^6$ packets",
    r"on a $2\times10^{-4}$ grid.}",
    r"\label{tab:ladder}",
    r"\begin{tabular}{@{}llcc@{}}",
    r"\toprule",
    r"axis & value & trough & vs frozen \\",
    r"\midrule",
]
axes = [
    ("packets ($\\mathrm{d}\\nu/\\nu=2\\times10^{-4}$)", [("A_emit5e+05", "$5\\times10^5$"), ("A_emit2e+06", "$2\\times10^6$"), ("A_emit8e+06", "$8\\times10^6$"), ("A_emit3e+07", "$3.2\\times10^7$")]),
    ("transport $\\mathrm{d}\\nu/\\nu$ ($8\\times10^6$)", [("B_dnut4.00e-04", "$4\\times10^{-4}$"), ("B_dnut2.00e-04", "$2\\times10^{-4}$"), ("B_dnut1.00e-04", "$1\\times10^{-4}$"), ("B_dnut4.17e-05", "$4.2\\times10^{-5}$"), ("B_dnut2.00e-05", "$2\\times10^{-5}$")]),
    ("spectrum $\\mathrm{d}\\nu/\\nu$", [("C_dnus5e-04", "$5\\times10^{-4}$"), ("C_dnus2e-04", "$2\\times10^{-4}$"), ("C_dnus1e-04", "$1\\times10^{-4}$")]),
    ("zones", [("D_nz25", "25"), ("D_nz100", "100"), ("D_nz101", "101 (half-zone shift)"), ("D_nz400", "400")]),
    ("anchor ($3.2\\times10^7$, $4.2\\times10^{-5}$, 5 seeds)", [("E_anchor_bb", "resolved"), ("E_anchor_exp", "expansion")]),
    ("zero opacity", [("Z_zero", "$\\rho\\times10^{-6}$")]),
]
for name, items in axes:
    first = True
    for tag, label in items:
        if tag not in rows:
            continue
        lead = name if first else ""
        first = False
        vs = pct(tag) if tag != "E_anchor_exp" else f"(Poisson value $0.4212$)"
        lines.append(f"{lead} & {label} & {cell(tag)} & {vs} \\\\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
table = "\n".join(lines)
print(table)

if "--insert" in sys.argv:
    p = ROOT / "docs/paper/manuscript.tex"
    s = p.read_text(encoding="utf-8")
    pat = re.compile(r"\\begin\{table\}\n\\centering\n\\small\n\\caption\{Single-line benchmark:.*?\\end\{table\}", re.S)
    if pat.search(s):
        s = pat.sub(lambda m: table, s)
    elif "%%LADDER_TABLE" in s:
        s = s.replace("%%LADDER_TABLE", table)
    else:
        sys.exit("no insertion point: add %%LADDER_TABLE to manuscript.tex")
    p.write_text(s, encoding="utf-8")
    print("inserted into manuscript.tex")
