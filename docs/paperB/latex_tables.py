#!/usr/bin/env python3
"""Paper B: the LaTeX macros and table fragments, generated from
paperB/FROZEN.json. Every number the manuscript prints is a macro defined
here; none is typed into the prose (the standing rule after Paper IV).

    .venv/bin/python docs/paperB/latex_tables.py      # writes numbers.tex, tab_*.tex
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = ROOT / "paperB/FROZEN.json"
IONS = ("57LaII", "58CeII", "60NdII")
SHORT = {"57LaII": "La", "58CeII": "Ce", "60NdII": "Nd"}


def f2(v): return f"{v:.2f}"
def f3(v): return f"{v:.3f}"
def f1(v): return f"{v:.1f}"
def pct0(v): return f"{100 * v:.0f}\\%"
def intv(v): return f"{int(round(v))}"
def com(v): return f"{int(round(v)):,}".replace(",", "{,}")
def word(v): return str(v)
def wordlow(v): return str(v).lower()


FMT = dict(f1=f1, f2=f2, f3=f3, pct0=pct0, int=intv, com=com, word=word, wordlow=wordlow)


def _g1(h, ion, key):
    return h["g1"]["per_ion"][ion][key]


def macros(h):
    """[(name, value, fmt)] -- the name is the LaTeX macro without backslash."""
    out = []
    out.append(("PBdmMax", h["g1"]["dm_max"], "f2"))
    out.append(("PBdcolMax", h["g1"]["dcolour_max"], "f2"))
    out.append(("PBGoneB", h["g1"]["B1"], "word"))
    out.append(("PBGoneDecision", h["g1"]["decision"], "word"))
    out.append(("PBGtwoHone", h["g2"]["H1"], "word"))
    out.append(("PBGtwoHtwo", h["g2"]["H2"], "word"))
    for ion in IONS:
        s = SHORT[ion]
        out += [(f"PB{s}NgStar", _g1(h, ion, "ng_star"), "int"),
                (f"PB{s}EpsStar", _g1(h, ion, "eps_star"), "f2"),
                (f"PB{s}EpsErr", _g1(h, ion, "eps_max_dm"), "f2"),
                (f"PB{s}Rerr", _g1(h, ion, "e_R"), "f3"),
                (f"PB{s}EvPerPkt", _g1(h, ion, "events_per_packet"), "f1"),
                (f"PB{s}Klocal", h["g2"]["per_ion"][ion]["k_local"], "int"),
                (f"PB{s}Lstar", h["g2"]["per_ion"][ion]["L_star"], "f3"),
                (f"PB{s}Shift", abs(h["r2m"]["per_ion"][ion]["shift_max"]), "f2"),
                (f"PB{s}Rebuilt", h["r2m"]["per_ion"][ion]["rebuilt_band"], "f3"),
                (f"PB{s}Downward", h["r2m"]["per_ion"][ion]["downward_band"], "f2"),
                (f"PB{s}Survives", h["r2m"]["per_ion"][ion]["survives"], "wordlow"),
                (f"PB{s}NExit", h["audit"]["per_ion"][ion]["n_exit"], "com"),
                (f"PB{s}NLines", h["audit"]["per_ion"][ion]["n_lines"], "com"),
                (f"PB{s}KbTables", h["audit"]["per_ion"][ion]["kb_tables"], "int"),
                (f"PB{s}Conc", h["audit"]["per_ion"][ion]["conc90"], "com")]
        c = h["g2"]["contrast"][ion]
        out += [(f"PB{s}EventGlobal", c["event_global"], "f3"), (f"PB{s}EventLocal", c["event_local"], "f3"),
                (f"PB{s}BandGlobal", c["band_global"], "f3"), (f"PB{s}BandLocal", c["band_local"], "f3"),
                (f"PB{s}EventRatio", c["event_ratio"], "f1"), (f"PB{s}BandRatio", c["band_ratio"], "f1")]
        kg = h["g2"]["per_ion"][ion]["k_global"]
        out.append((f"PB{s}Kglobal", "undefined" if kg is None else str(int(kg)), "word"))
    out.append(("PBContrastK", h["g2"]["contrast"][IONS[0]]["k"], "int"))
    out.append(("PBParamsGlobal", h["g2"]["contrast"][IONS[0]]["params_global"], "com"))
    out.append(("PBParamsLocal", h["g2"]["contrast"][IONS[0]]["params_local"], "com"))
    # ---- G3 ----
    g3 = h["g3"]
    out += [("PBGthreeCone", g3["C1"], "word"), ("PBGthreeCtwo", g3["C2"], "word"), ("PBGthreeCthree", g3["C3"], "word"),
            ("PBGthreeCfour", g3["C4"], "word"), ("PBGthreeCfive", g3["C5"], "word"),
            ("PBNgT", g3["ng_t"], "int"), ("PBNAxes", g3["n_axes"], "int"), ("PBNStatesPerIon", g3["n_states_per_ion"], "int")]
    AXN = {"T": "T", "D": "D", "J": "J", "P": "P"}
    for ion in IONS:
        s = SHORT[ion]; r = g3["per_ion"][ion]
        out += [(f"PB{s}Cone", r["C1"], "word"), (f"PB{s}Ctwo", r["C2"], "word"), (f"PB{s}Cthree", r["C3"], "word"),
                (f"PB{s}TransferHolds", ", ".join(r["transfer_holds_axes"]) or "none", "word"),
                (f"PB{s}TransferFails", ", ".join(r["transfer_fails_axes"]) or "none", "word"),
                (f"PB{s}InterpPasses", ", ".join(r["interp_passes_axes"]) or "none", "word"),
                (f"PB{s}NStates", r["n_readable"], "int"), (f"PB{s}NGray", r["n_gray"], "int"),
                (f"PB{s}NFreshFail", r["n_fresh_fail"], "int"), (f"PB{s}NFreshFailBand", r["n_fresh_fail_band"], "int"),
                (f"PB{s}WorstFreshBand", r["worst_fresh_band"], "f3"), (f"PB{s}WorstFreshColour", r["worst_fresh_colour"], "f3"),
                (f"PB{s}WorstColourOver", r["worst_fail_colour_over"], "f3"), (f"PB{s}LeastColourOver", r["least_fail_colour_over"], "f3"),
                (f"PB{s}KrecMax", r["k_rec_max"], "int")]
        for st in r["states"]:
            if st["interp_whole_band"] is not None:
                ax = AXN[st["axis"]]
                out += [(f"PB{s}Fix{ax}", st["anchor_band"], "f3"), (f"PB{s}Fresh{ax}", st["fresh_band"], "f3"),
                        (f"PB{s}IntW{ax}", st["interp_whole_band"], "f3"), (f"PB{s}IntM{ax}", st["interp_matrix_band"], "f3"),
                        (f"PB{s}Lam{ax}", st["lam"], "f2"), (f"PB{s}Cls{ax}", st["cls"], "word")]
    nt = g3["nd_trajectory"]
    out += [("PBNdTrajAnchor", nt["anchor_band"], "f3"), ("PBNdTrajWhole", nt["interp_whole_band"], "f3"),
            ("PBNdTrajWholeColour", nt["interp_whole_colour"], "f3"), ("PBNdTrajMatrix", nt["interp_matrix_band"], "f3"),
            ("PBNdTrajMatrixColour", nt["interp_matrix_colour"], "f3"), ("PBNdTrajFresh", nt["fresh_band"], "f3"),
            ("PBNdTrajKrec", nt["k_rec"], "int"), ("PBNdTrajMatrixPasses", "passes" if nt["matrix_only_passes"] else "fails", "word")]
    pb, pc = g3["part_b"], g3["part_c"]
    out += [("PBKmix", "undefined" if pb["k_mix"] is None else str(int(pb["k_mix"])), "word"), ("PBKdirect", pb["k_direct"], "int"),
            ("PBAmixBestN", pb["amix_best_n"], "int"), ("PBAmixBestBand", pb["amix_best_band"], "f3"), ("PBAmixBestColour", pb["amix_best_colour"], "f3"),
            ("PBAdirectBand", pb["adirect_band"], "f3"), ("PBAdirectColour", pb["adirect_colour"], "f3"),
            ("PBBlendKrec", pc["k_rec"], "int"), ("PBBlendBand", pc["band"], "f3"), ("PBBlendColour", pc["colour"], "f3"),
            ("PBBlendEps", pc["eps_star"], "f2"), ("PBBlendEpsErr", pc["eps_err"], "f2"), ("PBBlendNIons", pc["n_ions"], "int")]
    return out


def numbers_tex(h):
    lines = ["% generated by docs/paperB/latex_tables.py from paperB/FROZEN.json -- do not edit"]
    for name, val, fmt in macros(h):
        lines.append(f"\\newcommand{{\\{name}}}{{{FMT[fmt](val)}}}")
    return "\n".join(lines) + "\n"


def tab_gates(h):
    """The result table: per ion, G1's N_g* and eps*, G2's counts, H2's L*."""
    L = ["% generated by docs/paperB/latex_tables.py -- do not edit",
         "\\begin{tabular}{lccccccc}", "\\hline",
         "ion & $N_g^\\star$ & $\\epsilon^\\star$ & error at $\\epsilon^\\star$ & $K^\\star_{\\rm local}$ & "
         "$K^\\star_{\\rm global}$ & $L^\\star$ & $R_8$ vs the new reference \\\\",
         " & & & [mag] & & & & [mag] \\\\", "\\hline"]
    for ion in IONS:
        s = SHORT[ion]
        L.append(f"{h['ion_names'][ion]} & \\PB{s}NgStar & \\PB{s}EpsStar & \\PB{s}EpsErr & \\PB{s}Klocal & "
                 f"\\PB{s}Kglobal & \\PB{s}Lstar & \\PB{s}Rebuilt \\\\")
    L += ["\\hline", "\\end{tabular}"]
    return "\n".join(L) + "\n"


def tab_contrast(h):
    """The mechanism table: events against observables at 32 archetypes."""
    L = ["% generated by docs/paperB/latex_tables.py -- do not edit",
         "\\begin{tabular}{lcccc}", "\\hline",
         "ion & $m_{\\rm event}$ global & $m_{\\rm event}$ local & max $|\\Delta m|$ global & max $|\\Delta m|$ local \\\\",
         " & & & [mag] & [mag] \\\\", "\\hline"]
    for ion in IONS:
        s = SHORT[ion]
        L.append(f"{h['ion_names'][ion]} & \\PB{s}EventGlobal & \\PB{s}EventLocal & \\PB{s}BandGlobal & \\PB{s}BandLocal \\\\")
    L += ["\\hline", "\\end{tabular}"]
    return "\n".join(L) + "\n"


def tab_generality(h):
    """G3 per ion: the readings, where transfer holds, how often the fresh
    operator fails, the largest resolution it needed, and the interpolation
    result at each axis' interior point (whole operator / matrix only)."""
    L = ["% generated by docs/paperB/latex_tables.py -- do not edit",
         "\\begin{tabular}{lccccccccc}", "\\hline",
         "ion & C1 & C2 & C3 & transfer holds & fresh $R_{16}$ fails & $\\max K^\\star_{\\rm rec}$ & "
         "interp.\\ $T$ & interp.\\ $D$ & interp.\\ $J$ \\\\",
         " & & & & (axes) & (of states) & & whole / matrix & whole / matrix & whole / matrix \\\\", "\\hline"]
    for ion in IONS:
        s = SHORT[ion]
        L.append(f"{h['ion_names'][ion]} & \\PB{s}Cone & \\PB{s}Ctwo & \\PB{s}Cthree & \\PB{s}TransferHolds & "
                 f"\\PB{s}NFreshFail / \\PB{s}NStates & \\PB{s}KrecMax & "
                 f"\\PB{s}IntWT / \\PB{s}IntMT & \\PB{s}IntWD / \\PB{s}IntMD & \\PB{s}IntWJ / \\PB{s}IntMJ \\\\")
    L += ["\\hline", "\\end{tabular}"]
    return "\n".join(L) + "\n"


def main():
    h = json.loads(FROZEN.read_text())
    written = []
    for name, text in (("numbers.tex", numbers_tex(h)), ("tab_gates.tex", tab_gates(h)), ("tab_contrast.tex", tab_contrast(h)),
                       ("tab_generality.tex", tab_generality(h))):
        p = HERE / name; p.write_text(text); written.append(p)
    return written


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
