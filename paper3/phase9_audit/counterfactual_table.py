"""Paper III: the F38 counterfactual table (§4.35) as a generator.

Reads the survey's `survey.json` (which carries all four legs; `audit.json`
carries no `sobolev_group` leg) and emits the A / B / C / A+B / C-(A+B) table
of signed band-3800 errors as markdown and JSON.

    A = sobolev_group     (exact opacity, grouped redistribution)
    B = expansion_branch  (grouped opacity, exact redistribution)
    C = expansion_group   (both)

Usage: python counterfactual_table.py [--json survey.json] [--out counterfactual_table.json]
"""
import argparse, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SURVEY = HERE.parent / "phase8_survey" / "survey.json"
IONS = ("60NdII", "57LaII", "59PrII", "58CeII")   # the four with a live band
LABEL = {"60NdII": "Nd II", "57LaII": "La II", "59PrII": "Pr II", "58CeII": "Ce II"}


def table(rows, ions=IONS):
    by = {r["ion"]: r for r in rows}
    out = []
    for ion in ions:
        r = by[ion]
        A = r["legs"]["sobolev_group"]["band3800"]
        B = r["legs"]["expansion_branch"]["band3800"]
        C = r["legs"]["expansion_group"]["band3800"]
        out.append({"ion": ion, "label": LABEL.get(ion, ion), "S": r["band_S_band"],
                    "A": A, "B": B, "C": C, "A_plus_B": A + B, "interaction": C - (A + B)})
    return sorted(out, key=lambda x: x["S"])


def markdown(t):
    lines = ["| ion | S | A | B | C | A+B | **C−(A+B)** |", "|---|---|---|---|---|---|---|"]
    for r in t:
        p = lambda k: f"{100*r[k]:+.1f}%".replace("-", "−")
        lines.append(f"| {r['label']} | {r['S']:.1f} | {p('A')} | {p('B')} | {p('C')} | "
                     f"{p('A_plus_B')} | {p('interaction')} |")
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=str(SURVEY))
    ap.add_argument("--out", default=str(HERE / "counterfactual_table.json"))
    a = ap.parse_args()
    src = json.load(open(a.json))
    t = table(src["rows"])
    print(markdown(t))
    Path(a.out).write_text(json.dumps({"source": a.json, "n": src.get("n"), "rows": t}, indent=1))
    print(f"wrote {a.out}")
