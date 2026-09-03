"""Paper III: every number the manuscript quotes, and its three tables, from
the frozen analysis.

    numbers.tex        one \\newcommand per headline number (paper3/FROZEN.json
                       `headline`, through `MACROS`); rounding happens here
                       and nowhere else, so a number that is not in FROZEN
                       cannot appear in the paper
    tab_verdict.tex    Table 1, the verdict per lanthanide fraction
    tab_grid.tex       Extended Data Table 1, the 27 grid points
    tab_scenarios.tex  Extended Data Table 2, the observing scenarios and Gate 3

`main(h, dest, out_dir)` is what `paper3/freeze.py` calls (h = the headline
dict, dest = the derived-JSON paths); `check_structure.py` regenerates
`numbers.tex` in memory through `numbers_tex(h)` and requires byte equality
with the committed file.
"""
import ast, json, math, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
XS = ("0.001", "0.01", "0.1")
XSUF = ("Xlow", "Xmid", "Xhigh")


# ------------------------------------------------------------- formatting
def _num(x, nd):
    """A number for running text: minus signs in math mode, no digits lost."""
    s = f"{x:.{nd}f}"
    return f"$-{s[1:]}$" if s.startswith("-") else s


def _of(v):                 # [k, n] -> "27 of 27"
    return f"{v[0]} of {v[1]}"


def _slash(v):              # [k, n] or {survives, eligible} -> "27/27"
    if isinstance(v, dict):
        v = (v["survives"], v["eligible"])
    return f"{v[0]}/{v[1]}"


def _rng(v, nd, unit=""):   # [a, b] -> "0.96--2.84"
    a, b = (f"{x:.{nd}f}" for x in v)
    if a.startswith("-") or b.startswith("-"):
        return f"{_num(v[0], nd)} to {_num(v[1], nd)}{unit}"
    return f"{a}--{b}{unit}"


def _pct(v, nd=0):          # fraction or [k, n] -> "56\%"
    f = v[0] / v[1] if isinstance(v, (list, tuple)) else v
    return f"{100 * f:.{nd}f}\\%"


FMT = {
    "int": lambda v: f"{int(round(v))}",
    # exact powers of ten as $10^{k}$ (packet counts)
    "pow": lambda v: (f"$10^{{{int(round(math.log10(v)))}}}$" if abs(v - 10 ** round(math.log10(v))) < 1e-9 * v else f"{int(round(v))}"),
    "f1": lambda v: _num(v, 1), "f2": lambda v: _num(v, 2), "f3": lambda v: _num(v, 3),
    "of": _of, "slash": _slash,
    "rng1": lambda v: _rng(v, 1), "rng2": lambda v: _rng(v, 2),
    "pct": _pct, "pct1": lambda v: _pct(v, 1),
    "pctrng1": lambda v: _rng([100 * x for x in v], 1, "\\%"),
    # the headline amplitude, rounded to whole magnitudes -- the precision the
    # chain-cap test supports ("1--3 mag", never "1.0--2.8")
    "rngwhole": lambda v: f"{int(math.floor(v[0] + 0.5)):d}--{int(math.floor(v[1] + 0.5)):d}",
}

# (macro name, headline key, format, quoted-in-manuscript)
# "byx" expands to three macros <name>Xlow/Xmid/Xhigh, each "k of n".
MACROS = [
    # Gate 2 and the tangent spaces
    ("nPoints", "gate2.T0.n", "int", True),
    ("GateTwoCB", "gate2.T0.cb", "int", True),
    ("GateTwoMedianR", "gate2.T0.median_R", "f2", True),
    ("GateTwoMedianChi", "gate2.T0.median_chi2_res_dof", "int", True),
    ("GateTwoBinnedCB", "gate2.T0.C_binned.cb", "int", True),
    ("TOneCB", "gate2.T1.cb", "int", True),
    ("TOneCB", "gate2.T1.cb_by_X", "byx", True),
    ("TOneMedianChi", "gate2.T1.median_chi2_res_dof", "int", True),
    ("TOneLeftover", "gate2.T1.leftover_gt4", "int", True),
    ("TOneDetermined", "gate2.T1.n_determined", "int", True),
    ("TOneALMedian", "gate2.T1.abs_aL_median", "f1", True),
    ("TOneALMax", "gate2.T1.abs_aL_max", "f1", True),
    ("TTwoCB", "gate2.T2.cb", "int", True),
    ("TTwoCB", "gate2.T2.cb_by_X", "byx", True),
    ("TTwoMedianChi", "gate2.T2.median_chi2_res_dof", "int", True),
    ("TTwoLeftover", "gate2.T2.leftover_gt4", "int", True),
    ("TTwoDetermined", "gate2.T2.n_determined", "int", True),
    ("TTwoDlnT", "gate2.T2.aT_median", "f2", True),
    ("TThreeCB", "gate2.T3.cb", "int", True),
    ("TThreeCB", "gate2.T3.cb_by_X", "byx", True),
    ("TThreeMedianChi", "gate2.T3.median_chi2_res_dof", "int", True),
    ("TThreeLeftover", "gate2.T3.leftover_gt4", "int", True),
    ("TThreeDetermined", "gate2.T3.n_determined", "int", True),
    ("TThreeUnderdetermined", "gate2.T3.underdetermined", "int", True),
    ("TThreeBeyondLinear", "gate2.T3.lin2c_gt1", "of", True),
    ("ControlCC", "control.A_redist.cc", "int", True),
    ("ControlN", "control.A_redist.n", "int", True),
    ("BOpacityCB", "control.B_opacity.cb", "int", True),
    ("ThrChi", "thresholds.chi2_small", "int", True),
    ("ThrR", "thresholds.R_max", "f1", True),
    ("ThrSignif", "thresholds.signif_min", "int", True),
    # the grid
    ("nModels", "grid.n_models", "int", True),
    ("nCells", "grid.n_cells", "int", True),
    ("nCellsRan", "grid.n_cells_ran", "int", True),
    ("nRedone", "grid.n_redone_cells", "int", True),
    ("WorstColourRange", "colour.worst_per_model_range", "rng2", True),
    ("WorstColourRangeCoarse", "colour.worst_per_model_range", "rng1", True),
    ("HeadlineColourRange", "colour.worst_per_model_range", "rngwhole", True),
    ("gDmRange", "colour.g_dm_range", "rng1", True),
    ("gN", "colour.g_n", "int", True),
    ("KDmMax", "colour.K_dm_max", "f1", True),
    ("KN", "colour.K_n", "int", True),
    ("NIRNegative", "colour.nir_negative", "of", True),
    ("FloorNWell", "floor.n_well", "int", True),
    ("FloorNMin", "floor.n_min", "pow", True),
    ("FloorMedian", "floor.median", "f3", True),
    ("FloorNinety", "floor.p90", "f3", True),
    ("FloorMax", "floor.max", "f2", True),
    ("FloorRedoneRange", "floor.redone_range", "rng2", True),
    # Gate 3
    ("Distance", "gate3.distance_Mpc", "int", True),
    ("DenseSurvives", "gate3.dense.T0.survives", "int", True),
    ("DenseEligible", "gate3.dense.T0.eligible", "int", True),
    ("SparseSurvives", "gate3.sparse.T0.survives", "int", True),
    ("SparseEligible", "gate3.sparse.T0.eligible", "int", True),
    ("OpticalSurvives", "gate3.optical.T0.survives", "int", True),
    ("OpticalEligible", "gate3.optical.T0.eligible", "int", True),
    ("DenseTOne", "gate3.dense.T1.survives", "int", True),
    ("DenseTOne", "gate3.dense.T1.by_X", "byx", True),
    ("DenseTTwo", "gate3.dense.T2.survives", "int", True),
    ("DenseTThree", "gate3.dense.T3.survives", "int", True),
    ("SparseTOne", "gate3.sparse.T1.survives", "int", True),
    ("SparseTOneEligible", "gate3.sparse.T1.eligible", "int", True),
    ("OpticalTOne", "gate3.optical.T1.survives", "int", True),
    ("OpticalTOneEligible", "gate3.optical.T1.eligible", "int", True),
    ("DenseNIRShare", "gate3.dense.T0.median_nir_share", "pct", True),
    # T_eff validation
    ("IllumCos", "tscale.illum.cos", "f2", True),
    ("IllumNorm", "tscale.illum.norm_ratio", "f2", True),
    ("GasCos", "tscale.gas.cos", "f2", True),
    ("GasNorm", "tscale.gas.norm_ratio", "f2", True),
    # chain cap
    ("ChainCells", "chain.n_cells", "int", True),
    ("ChainBase", "chain.base_cap", "int", True),
    ("ChainTop", "chain.top_cap", "int", True),
    ("ChainTrapped", "chain.trapped_range", "pctrng1", True),
    ("ChainTrappedTop", "chain.trapped_range_top", "pctrng1", True),
    ("ChainDmChange", "chain.dm_change_8000_range", "rng2", True),
    ("ChainDmRefChange", "chain.dm_ref_change_8000_range", "rng2", True),
    ("ChainSigns", "chain.signs_kept", "of", True),
    ("ChainCriterion", "chain.criterion_met", "of", True),
    ("ChainRelMax", "chain.rel_max", "pct", True),
    ("ChainBinnedSigns", "chain.C_binned.signs_kept", "of", True),
    ("ChainBinnedCriterion", "chain.C_binned.criterion_met", "of", True),
    ("ChainCB", "chain.cb_with_override", "int", True),
    ("ChainN", "chain.n_with_override", "int", True),
    ("ChainMedianChi", "chain.median_chi2_res_dof_override", "int", True),
    ("ChainMedianR", "chain.median_R_override", "f2", True),
    # the allowance
    ("SysLive", "syserr.n_live.C_both", "int", True),
    ("SysKeys", "syserr.n_keys", "int", True),
    ("SysGtHalf", "syserr.frac_gt_0p5.C_both", "pct", True),
    ("SysGtHalfCount", "syserr.frac_gt_0p5.C_both", "of", True),
    ("SysGtOne", "syserr.frac_gt_1.C_both", "pct", True),
    ("SysGtOneCount", "syserr.frac_gt_1.C_both", "of", True),
    ("SysBinnedGtHalf", "syserr.frac_gt_0p5.C_binned", "pct", True),
    ("SysBinnedGtOne", "syserr.frac_gt_1.C_binned", "pct", True),
    ("SysControlGtHalf", "syserr.frac_gt_0p5.A_redist", "pct", True),
    ("SysSign", "syserr.sign_pattern.C_both", "of", True),
    ("SysBinnedSign", "syserr.sign_pattern.C_binned", "of", True),
    ("SysControlSign", "syserr.sign_pattern.A_redist", "of", True),
    ("SysOneMode", "syserr.one_mode.C_both", "f2", True),
    ("SysBinnedOneMode", "syserr.one_mode.C_binned", "f2", True),
    ("SysControlOneMode", "syserr.one_mode.A_redist", "f2", True),
    ("SysSVD", "syserr.svd_filled.C_both", "f2", True),
    ("SysSVDKeys", "syserr.svd_n_keys", "int", True),
    ("SysNullMedian", "syserr.null_median.C_both", "f2", True),
    ("SysNullNinetyFive", "syserr.null_p95.C_both", "f2", True),
    ("SysMPScale", "syserr.mp_scale", "f2", True),
    ("SysChiOne", "syserr.chi2_equiv_1mag_median.C_both", "f2", True),
    ("SysChiOneRange", "syserr.chi2_equiv_1mag_range.C_both", "rng2", True),
    ("SysChiHalf", "syserr.chi2_equiv_0p5mag_median.C_both", "f1", True),
    ("SysBinnedChiOne", "syserr.chi2_equiv_1mag_median.C_binned", "f2", True),
    ("SysControlChiOne", "syserr.chi2_equiv_1mag_median.A_redist", "f3", True),
    ("SysResidOneMode", "syserr.residual_T1.one_mode", "f2", True),
    ("SysResidGtHalf", "syserr.residual_T1.frac_gt_0p5", "pct", True),
    ("SysResidChiOne", "syserr.residual_T1.chi2_equiv_1mag_median", "f2", True),
]


def macros(h):
    """The (name, text) pairs `numbers.tex` defines, in MACROS order."""
    out = []
    for name, key, fmt, _ in MACROS:
        v = h[key]
        if fmt == "byx":
            for x, suf in zip(XS, XSUF):
                out.append((name + suf, _of(v[x])))
        else:
            out.append((name, FMT[fmt](v)))
    names = [n for n, _ in out]
    assert len(names) == len(set(names)), "duplicate macro name"
    return out


def quoted_names():
    """Macro names the manuscript must use at least once."""
    out = []
    for name, _, fmt, q in MACROS:
        if q:
            out += [name + s for s in XSUF] if fmt == "byx" else [name]
    return out


def numbers_tex(h):
    lines = ["% generated by docs/paper3/latex_tables.py from paper3/FROZEN.json -- do not edit",
             "% every number the manuscript quotes; rounding happens here only"]
    lines += [f"\\newcommand{{\\{n}}}{{{t}\\xspace}}" for n, t in macros(h)]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- tables
def _j(p):
    return json.loads(Path(p).read_text())


def _xkey(x):
    return f"{x:g}"


XLAB = {"0.001": r"$10^{-3}$", "0.01": r"$10^{-2}$", "0.1": r"$10^{-1}$"}


def _median(v):
    v = sorted(v)
    n = len(v)
    return float("nan") if n == 0 else (v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2]))


def verdict_rows(h, dest):
    """Table 1: one row per X_lan."""
    T0, T1 = _j(dest["sensitivity"]), _j(dest["sensitivity_T1"])
    g, ob = _j(dest["grid_table"]), _j(dest["observability"])
    rows = []
    for x in XS:
        pts = [p for p in g["points"] if _xkey(p["point"][2]) == x]
        worst = [p["C_both"]["worst_dcolor"] for p in pts if p["C_both"]["worst_dcolor"] == p["C_both"]["worst_dcolor"]]
        cells = [c for c in g["cells"] if _xkey(c["point"][2]) == x and c["ran"] and c["n_used"] >= g["floor"]["n_min"]
                 and c["floor"] == c["floor"]]
        floor = _median([c["floor"] for c in cells])
        chi_T1 = _median([r["legs"]["C_both"]["chi2_res_dof"] for k, r in T1["points"].items()
                          if _xkey(ast.literal_eval(k)[2]) == x and r["legs"]["C_both"].get("status") == "ok"
                          and not r["legs"]["C_both"].get("underdetermined")])
        s = ob["summary"]["C_both"]
        rows.append({
            "X": XLAB[x],
            "worst": f"{min(worst):.1f}--{max(worst):.1f}",
            "floor": f"{floor:.3f}",
            "T0": _slash(h["gate2.T0.cb_by_X"][x]),
            "T1": _slash(h["gate2.T1.cb_by_X"][x]),
            "T1chi": f"{chi_T1:.0f}",
            "dense": _slash(s["dense"]["T0"]["by_X"][x]),
            "sparse": _slash(s["sparse"]["T0"]["by_X"][x]),
            "optical": _slash(s["optical"]["T0"]["by_X"][x]),
            "denseT1": _slash(s["dense"]["T1"]["by_X"][x]),
        })
    return rows


def _check_columns(tex):
    """Every body row must carry exactly as many cells as the preamble declares."""
    import re
    spec = re.search(r"\\begin\{tabular\}\{([^}]*)\}", tex).group(1)
    ncol = len(re.sub(r"[^lcrp]", "", spec))
    for line in tex.splitlines():
        if line.endswith(r"\\") and r"\multicolumn" not in line:
            n = line.count("&") + 1
            assert n == ncol, f"{n} cells in a {ncol}-column table: {line[:60]}"
    return tex


def tab_verdict(h, dest):
    rows = verdict_rows(h, dest)
    head = (r"\begin{tabular}{lccccccccc}" "\n" r"\toprule" "\n"
            r"$X_{\rm lan}$ & $|\Delta({\rm colour})|_{\max}$ & floor & Gate 2 & \multicolumn{2}{c}{free $L(t)$} & "
            r"\multicolumn{3}{c}{Gate 3, $(M,v,X)$} & Gate 3, free $L(t)$ \\" "\n"
            r" & (mag) & (mag) & C-B & C-B & $\chi^2_{\rm res}/{\rm dof}$ & dense & sparse & optical & dense \\" "\n"
            r"\midrule" "\n")
    body = "".join(f"{r['X']} & {r['worst']} & {r['floor']} & {r['T0']} & {r['T1']} & {r['T1chi']} & "
                   f"{r['dense']} & {r['sparse']} & {r['optical']} & {r['denseT1']} \\\\\n" for r in rows)
    return _check_columns(head + body + r"\bottomrule" "\n" r"\end{tabular}" "\n")


def grid_rows(dest):
    """Extended Data Table 1: one row per grid point."""
    S = {tg: _j(dest[f"sensitivity{s}"]) for tg, s in (("T0", ""), ("T1", "_T1"), ("T2", "_T2"), ("T3", "_T3"))}
    g = _j(dest["grid_table"])
    rows = []
    for p in g["points"]:
        key = str(tuple(p["point"]))
        cls = {}
        for tg, s in S.items():
            r = s["points"][key]["legs"]["C_both"]
            cls[tg] = ("underdet." if r.get("underdetermined") else r.get("cls") or "--") if r.get("status") == "ok" else "--"
        cells = [c for c in g["cells"] if c["point"] == p["point"] and c["ran"]]
        trapped = max(c["trapped_frac"] for c in cells)
        rows.append({
            "M": f"{p['point'][0]:g}", "v": f"{p['point'][1]:g}", "X": XLAB[_xkey(p["point"][2])],
            "N": S["T0"]["points"][key]["legs"]["C_both"].get("N", 0),
            "n_used": f"{p['n_used_min'] / 1e3:.0f}--{p['n_used_max'] / 1e3:.0f}",
            "trapped": f"{100 * trapped:.1f}",
            "floor": f"{p['floor_max_well']:.3f}" if p["floor_max_well"] is not None and p["floor_max_well"] == p["floor_max_well"] else "--",
            "cboth": f"{p['C_both']['worst_dcolor']:.2f} ({p['C_both']['worst_key']})",
            "cbinned": f"{p['C_binned']['worst_dcolor']:.2f} ({p['C_binned']['worst_key']})",
            "redone": len(p["redone_epochs"]), "floored": len(p["floored_epochs"]),
            **cls,
        })
    return rows


def tab_grid(dest):
    rows = grid_rows(dest)
    head = (r"\begin{tabular}{cccrcrrcccccc}" "\n" r"\toprule" "\n"
            r"$M_{\rm ej}$ & $v_{\rm ej}$ & $X_{\rm lan}$ & $N$ & $n_{\rm used}$ & trapped & floor & "
            r"$|\Delta{\rm col}|_{\max}$ & $|\Delta{\rm col}|_{\max}$ & \multicolumn{4}{c}{class} \\" "\n"
            r"($M_\odot$) & ($c$) & & & ($10^3$) & (\%) & (mag) & C & C$_{\rm bin}$ & T0 & T1 & T2 & T3 \\" "\n"
            r"\midrule" "\n")
    body = "".join(f"{r['M']} & {r['v']} & {r['X']} & {r['N']} & {r['n_used']} & {r['trapped']} & {r['floor']} & "
                   f"{r['cboth']} & {r['cbinned']} & {r['T0']} & {r['T1']} & {r['T2']} & {r['T3']} \\\\\n" for r in rows)
    return _check_columns(head + body + r"\bottomrule" "\n" r"\end{tabular}" "\n")


def _epochs(obs, bands):
    ts = sorted({t for b in bands for t in obs.get(b, [])})
    return ", ".join(f"{t:g}" for t in ts) if ts else "--"


def tab_scenarios(dest):
    ob = _j(dest["observability"])
    head = (r"\begin{tabular}{llcccccc}" "\n" r"\toprule" "\n"
            r"scenario & epochs (d) & depth & $\sigma_{\rm sys}$ & \multicolumn{4}{c}{eligible / survives} \\" "\n"
            r" & optical; NIR & opt / NIR & opt / NIR & T0 & T1 & T2 & T3 \\" "\n" r"\midrule" "\n")
    body = ""
    for name, sc in ob["scenarios"].items():
        opt, nir = _epochs(sc["obs"], "griz"), _epochs(sc["obs"], "JHK")
        d = sc["depth"]; s = sc["sys"]
        depth = f"{d['g']:.1f} / {d.get('J', float('nan')):.1f}".replace("nan", "--")
        sy = f"{s['g']:.2f} / {s.get('J', float('nan')):.2f}".replace("nan", "--")
        cols = " & ".join(f"{ob['summary']['C_both'][name][tg]['eligible']} / {ob['summary']['C_both'][name][tg]['survives']}"
                          for tg in ("T0", "T1", "T2", "T3"))
        body += f"{name} & {opt}; {nir} & {depth} & {sy} & {cols} \\\\\n"
    return _check_columns(head + body + r"\bottomrule" "\n" r"\end{tabular}" "\n")


# ------------------------------------------------------------------ SI tables
LEGLAB = {"A_redist": "A", "B_opacity": "B", "C_both": "C", "C_binned": r"C$_{\rm bin}$"}


def _tab(spec, header_rows, body_rows):
    head = "\\begin{tabular}{" + spec + "}\n\\toprule\n" + "".join(r + " \\\\\n" for r in header_rows) + "\\midrule\n"
    body = "".join(" & ".join(r) + " \\\\\n" for r in body_rows)
    return _check_columns(head + body + "\\bottomrule\n\\end{tabular}\n")


def _f(v, nd=2):
    if v is None or v != v:
        return "--"
    s = f"{abs(v):.{nd}f}"
    return f"$-{s}$" if v < 0 and float(s) != 0 else s  # no "-0.00"


def tab_si_robustness(dest):
    """SI Table 1: Gate 2 summary under every rule and tangent space, all four legs."""
    variants = [("T0, conserving, floored excluded (baseline)", "sensitivity"),
                ("T0, floored epochs included", "sensitivity_floored_incl"),
                ("T0, absorbing core", "sensitivity_absorbing"),
                ("T0, chain cap 8000 at the four cells", "sensitivity_chain8000"),
                ("T1: + free $L$ per epoch", "sensitivity_T1"),
                ("T2: + free colour temperature", "sensitivity_T2"),
                ("T3: + free blue component", "sensitivity_T3")]
    rows = []
    for label, key in variants:
        sm = _j(dest[key])["summary"]
        for i, leg in enumerate(("C_both", "C_binned", "B_opacity", "A_redist")):
            m = sm[leg]
            rows.append([label if i == 0 else "", LEGLAB[leg], str(m["n"]), str(m["C-C"]), str(m["C-A"]), str(m["C-B"]),
                         str(m["underdetermined"]), _f(m["median_R"]), _f(m["median_chi2_res_dof"], 1)])
    return _tab("llrrrrrcc",
                [r"rule / tangent space & leg & $n$ & C-C & C-A & C-B & underdet. & median $R$ & median $\chi^2_{\rm res}/{\rm dof}$"],
                rows)


def tab_si_points(dest):
    """SI Table 2: every grid point under T0-T3: R, chi2_res/dof, class."""
    S = {tg: _j(dest[f"sensitivity{s}"]) for tg, s in (("T0", ""), ("T1", "_T1"), ("T2", "_T2"), ("T3", "_T3"))}
    rows = []
    for key, p0 in S["T0"]["points"].items():
        pt = ast.literal_eval(key)
        cells = [f"{pt[0]:g}", f"{pt[1]:g}", XLAB[_xkey(pt[2])], str(p0["legs"]["C_both"].get("N", 0)), _f(p0["noise_floor"], 3)]
        for tg in ("T0", "T1", "T2", "T3"):
            r = S[tg]["points"][key]["legs"]["C_both"]
            if r.get("status") != "ok":
                cells += ["--", "--", "--"]; continue
            cells += [_f(r["R"]), _f(r["chi2_res_dof"], 1), "underdet." if r.get("underdetermined") else r["cls"]]
        rows.append(cells)
    hdr1 = (r"$M_{\rm ej}$ & $v_{\rm ej}$ & $X_{\rm lan}$ & $N$ & floor & \multicolumn{3}{c}{T0} & \multicolumn{3}{c}{T1} & "
            r"\multicolumn{3}{c}{T2} & \multicolumn{3}{c}{T3}")
    sub = r"$R$ & $\chi^2_{\rm res}/{\rm dof}$ & class"
    hdr2 = r"($M_\odot$) & ($c$) & & & (mag) & " + " & ".join([sub] * 4)
    return _tab("ccc" + "r" * 2 + "ccc" * 4, [hdr1, hdr2], rows)


def tab_si_chain(dest):
    """SI Table 3: the chain-cap test at the four cells, per cap."""
    c = _j(dest["chain_table"])
    rows = []
    for cell in sorted(c["cells"], key=lambda x: (-x["runs"]["2000"]["trapped_frac"])):
        pt = cell["point"]
        for i, cap in enumerate(("2000", "4000", "8000")):
            run = cell["runs"][cap]
            leg = run["legs"]["C_both"]
            dc = leg["dcolor"]
            crit = "--" if cap == "2000" else f"{sum(leg['criterion_met'].values())}/{len(leg['criterion_met'])}"
            rows.append([f"({pt[0]:g}, {pt[1]:g}, {pt[2]:g}) at {cell['t_d']:g} d" if i == 0 else "",
                         cap, f"{100 * run['trapped_frac']:.1f}", _f(dc["g-r"]), _f(dc["i-J"]), _f(dc["J-K"]),
                         "--" if cap == "2000" else _f(run["max_dm_ref_change"]),
                         "--" if cap == "2000" else _f(leg["max_dm_change"]), crit])
    return _tab("llrrrrrrc",
                [r"cell & cap & trapped & $\Delta(g-r)$ & $\Delta(i-J)$ & $\Delta(J-K)$ & $\max|\delta m_{\rm ref}|$ & $\max|\delta d_{\rm RT}|$ & criterion",
                 r" & & (\%) & (mag) & (mag) & (mag) & (mag) & (mag) & met"],
                rows)


def tab_si_tscale(dest):
    """SI Table 4: the Monte Carlo temperature direction against the Planck proxy, per live observable."""
    t = _j(dest["tscale"])
    il, gas = t["variants"]["illumination_only"]["compare"]["C_both"], t["variants"]["with_T_gas"]["compare"]["C_both"]
    rows = []
    for key in il["proxy"]:
        b, td = key.split(",")
        rows.append([b, f"{float(td):g}", _f(il["proxy"][key]), _f(il["mc"].get(key)), _f(gas["mc"].get(key))])
    rows.append([r"\multicolumn{2}{l}{cosine with the proxy}", "", _f(il["cos"]), _f(gas["cos"])])
    rows.append([r"\multicolumn{2}{l}{norm ratio MC / proxy}", "", _f(il["norm_ratio"]), _f(gas["norm_ratio"])])
    return _tab("clrrr",
                [r"band & $t$ (d) & Planck proxy & MC, illumination only & MC, with $T_{\rm gas}$"],
                rows)


def tab_si_syserr(dest):
    """SI Table 5: the closure error per band-epoch key over the live points, and the one-mode shape."""
    y = _j(dest["syserr"])["legs"]
    L = y["C_both"]
    mode = dict(zip([f"{b},{t:g}" for b, t in L["keys"]], L["one_mode"]["mode"])) if isinstance(L["one_mode"].get("mode"), list) else {}
    rows = []
    for e, eb in zip(L["per_key"], y["C_binned"]["per_key"]):
        k = f"{e['band']},{e['t_d']:g}"
        rows.append([e["band"], f"{e['t_d']:g}", str(e["n"]), _f(e["median"]), f"{_f(e['p16'])} to {_f(e['p84'])}",
                     _f(eb["median"]), _f(mode.get(k))])
    return _tab("clrrcrr",
                [r"band & $t$ (d) & $n$ & median $d_{\rm RT}$ & 16--84\% & median, C$_{\rm bin}$ & mode shape",
                 r" & & & (mag) & (mag) & (mag) & "],
                rows)


SI_TABLES = {"si_tab_robustness.tex": tab_si_robustness, "si_tab_points.tex": tab_si_points,
             "si_tab_chain.tex": tab_si_chain, "si_tab_tscale.tex": tab_si_tscale, "si_tab_syserr.tex": tab_si_syserr}


def main(h, dest, out_dir=HERE):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    files = {"numbers.tex": numbers_tex(h), "tab_verdict.tex": tab_verdict(h, dest),
             "tab_grid.tex": tab_grid(dest), "tab_scenarios.tex": tab_scenarios(dest)}
    files.update({name: fn(dest) for name, fn in SI_TABLES.items()})
    written = []
    for name, text in files.items():
        (out_dir / name).write_text(text)
        written.append(out_dir / name)
    return written


if __name__ == "__main__":
    sys.path.insert(0, str(ROOT / "paper3"))
    import freeze
    fz = json.loads(freeze.FROZEN.read_text())
    for p in main(fz["headline"], freeze.canonical(), HERE):
        print("wrote", p)
