"""Paper III display items, from the frozen derived JSONs only.

Main text: Fig. 1 (the closure error on one kilonova), Fig. 2 (across the
grid), Fig. 3 (not a parameter shift), Fig. 4 (observability and structure).
Extended Data: ED Fig. 1 (source model and gate), ED Fig. 2 (T_eff
validation), ED Fig. 3 (chain cap), ED Fig. 4 (the A_redist floor).
Table 1 and ED Tables 1-2 are written by `latex_tables.py`.

Everything is read from the files `paper3/freeze.py` records in `dest`
(sensitivity*.json, grid_table.json, observability.json, tscale.json,
robustness/chain_table.json, syserr.json) plus the central-model transport
JSON (an *input* of the freeze, so its path is fixed here).  The one
computation that is not a file read is ED Fig. 1's `SourceModel(M, v)`,
which is deterministic, sub-second and part of the frozen inputs
(`sobolev/source.py`).

Output is byte-stable within one environment: PDF metadata dates are
stripped and `SOURCE_DATE_EPOCH` is set by the freeze.

    main(dest, out_dir) -> list of written paths
"""
import ast, json, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev import photometry as phot          # noqa: E402
from sobolev.source import SourceModel, DAY    # noqa: E402

GRID = ROOT / "paper3" / "phase12_grid" / "grid"
CENTRAL = GRID / "model_M0.01_v0.1_X0.01.json"
CENTRAL_POINT = "(0.01, 0.1, 0.01)"
BANDS = ("g", "r", "i", "z", "J", "H", "K")
COLS = ("g-r", "r-i", "i-z", "i-J", "J-K")
NIR_COLS = ("i-J", "J-K")
EPOCHS = (0.5, 1.0, 2.0, 3.0, 5.0, 7.0)
XS = (0.001, 0.01, 0.1)
MS = (0.003, 0.01, 0.03)
VS = (0.05, 0.1, 0.2)
TANGENTS = ("T0", "T1", "T2", "T3")
TG_LABEL = {"T0": r"$(M,v,X)$", "T1": r"$+L(t)$", "T2": r"$+T$", "T3": r"+blue comp."}
SCENARIOS = ("dense", "sparse", "optical")
N_GREY = 5          # per-key medians on fewer live points than this are greyed out

# Okabe-Ito, colourblind-safe
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "red": "#D55E00",
      "purple": "#CC79A7", "sky": "#56B4E9", "yellow": "#F0E442", "black": "#000000"}
LEGC = {"ref": OI["black"], "A_redist": OI["green"], "B_opacity": OI["sky"],
        "C_both": OI["red"], "C_binned": OI["orange"]}
LEGN = {"ref": "reference", "A_redist": "A (grouped $R_{ij}$)", "B_opacity": "B (expansion opacity)",
        "C_both": "C (both)", "C_binned": "C$_{\\rm bin}$ (binned $\\Sigma\\tau$)"}
XC = {0.001: OI["blue"], 0.01: OI["purple"], 0.1: OI["red"]}
XN = {0.001: r"$X_{\rm lan}=10^{-3}$", 0.01: r"$X_{\rm lan}=10^{-2}$", 0.1: r"$X_{\rm lan}=10^{-1}$"}
VMK = {0.05: "v", 0.1: "o", 0.2: "^"}
MSZ = {0.003: 10, 0.01: 18, 0.03: 30}
BANDC = dict(zip(BANDS, [OI["purple"], OI["blue"], OI["sky"], OI["green"], OI["yellow"], OI["orange"], OI["red"]]))

mm = 1 / 25.4
W1, W2 = 88 * mm, 180 * mm
STYLE = {"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7, "legend.fontsize": 6,
         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "lines.linewidth": 0.9,
         "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
         "xtick.major.size": 2.5, "ytick.major.size": 2.5, "legend.frameon": False,
         "font.family": "sans-serif", "pdf.fonttype": 42, "savefig.dpi": 300}
plt.rcParams.update(STYLE)


# ------------------------------------------------------------------ helpers
def _j(p):
    return json.loads(Path(p).read_text())


def _pt(key):
    return tuple(ast.literal_eval(key)) if isinstance(key, str) else tuple(key)


def letter(ax, s, dx=-0.18, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontweight="bold", fontsize=8, va="bottom", ha="left")


def logx(ax, ticks, labels=None):
    """Log x axis with only the named ticks (no minor labels)."""
    from matplotlib.ticker import NullFormatter, NullLocator
    ax.set_xscale("log")
    ax.set_xticks(list(ticks)); ax.set_xticklabels(labels or [f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_locator(NullLocator()); ax.xaxis.set_minor_formatter(NullFormatter())


def save(fig, name, out_dir):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}),
                      ("png", {"Software": None})):
        p = out_dir / f"{name}.{ext}"
        fig.savefig(p, metadata=meta, bbox_inches="tight", pad_inches=0.02)
        written.append(p)
    plt.close(fig)
    return written


def _grid_cells(g, point=None):
    return [c for c in g["cells"] if c["ran"] and (point is None or tuple(c["point"]) == point)]


def _point_order():
    """27 points sorted (X, M, v): X groups outermost."""
    return [(m, v, x) for x in XS for m in MS for v in VS]


# ------------------------------------------------------------------ Fig. 1
def fig1(dest, out_dir, epoch=2.0):
    """The closure error on one kilonova (central model)."""
    model, g = _j(CENTRAL), _j(dest["grid_table"])
    rows = model["rows"]
    row = next(r for r in rows if r["t_d"] == epoch)
    edges = phot.nu_edges(1000.0, 30000.0, 200)
    nu_c = 0.5 * (edges[:-1] + edges[1:])
    lam = 2.99792458e18 / nu_c                       # Angstrom, descending
    pbs = phot.load_passbands()

    fig = plt.figure(figsize=(W2, 62 * mm))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.42)
    ax = fig.add_subplot(gs[0]); bx = fig.add_subplot(gs[1]); cx = fig.add_subplot(gs[2])

    # a: spectra
    ymax = 0
    for leg in ("ref", "C_both", "C_binned", "A_redist"):
        lnu = np.array(row["ref"]["L_nu"] if leg == "ref" else row["legs"][leg]["L_nu"])
        y = nu_c * lnu / 1e40
        ymax = max(ymax, y.max())
        ax.plot(lam, y, color=LEGC[leg], label=LEGN[leg], lw=1.1 if leg == "ref" else 0.9,
                ls="-" if leg in ("ref", "C_both") else "--", zorder=3 if leg == "ref" else 2)
    for b, pb in pbs.items():
        T = np.array(pb.T) / max(pb.T)
        ax.fill_between(np.array(pb.lam), 0, T * 0.12 * ymax, color=BANDC[b], alpha=0.35, lw=0)
        ax.text(pb.lam_eff, 0.13 * ymax, b, ha="center", va="bottom", fontsize=6, color=BANDC[b])
    logx(ax, [3000, 5000, 10000, 20000], ["0.3", "0.5", "1", "2"]); ax.set_xlim(2500, 26000); ax.set_ylim(0, 1.15 * ymax)
    ax.set_xlabel(r"wavelength ($\mu$m)"); ax.set_ylabel(r"$\nu L_\nu$ ($10^{40}$ erg s$^{-1}$)")
    ax.legend(loc="upper right", handlelength=1.6)
    ax.text(0.03, 0.95, f"$t={epoch:g}$ d", transform=ax.transAxes, va="top")
    letter(ax, "a", dx=-0.13)

    # b: Delta colour vs t with the A_redist floor band
    t = np.array([r["t_d"] for r in rows])
    cells = {c["t_d"]: c for c in _grid_cells(g, _pt(CENTRAL_POINT))}
    CC = {"g-r": OI["blue"], "i-J": OI["green"], "J-K": OI["red"]}
    for col, c in CC.items():
        for leg, ls, mk in (("C_both", "-", "o"), ("C_binned", "--", "s"), ("A_redist", ":", None)):
            y = [cells[td]["legs"][leg]["dcolor"].get(col, np.nan) if td in cells else np.nan for td in t]
            bx.plot(t, y, ls=ls, color=c, marker=mk, ms=2.5, lw=0.9 if leg != "A_redist" else 0.6,
                    label=f"$\\Delta$({col})" if leg == "C_both" else None)
    fl = [cells[td]["floor"] if td in cells else np.nan for td in t]
    bx.fill_between(t, -np.array(fl), np.array(fl), color="0.5", alpha=0.25, lw=0, label="A floor")
    bx.axhline(0, color="k", lw=0.4)
    logx(bx, EPOCHS)
    bx.set_xlabel("time (d)"); bx.set_ylabel(r"$\Delta$colour $=$ closure $-$ reference (mag)")
    h = bx.get_legend_handles_labels()
    h[0].extend([Line2D([], [], color="k", ls="-", marker="o", ms=2.5, label="C (both)"),
                 Line2D([], [], color="k", ls="--", marker="s", ms=2.5, label="C$_{\\rm bin}$"),
                 Line2D([], [], color="k", ls=":", lw=0.6, label="A (control)")])
    bx.legend(handles=h[0], loc="lower left", handlelength=1.8, ncol=2, columnspacing=0.8, fontsize=5.5)
    letter(bx, "b")

    # c: Delta m matrix (C_both), masked cells hatched
    M = np.full((len(BANDS), len(EPOCHS)), np.nan)
    for j, td in enumerate(EPOCHS):
        if td in cells:
            for i, b in enumerate(BANDS):
                M[i, j] = cells[td]["legs"]["C_both"]["dm"].get(b, np.nan)
    vmax = np.nanmax(np.abs(M))
    im = cx.imshow(np.ma.masked_invalid(M), cmap="RdBu", vmin=-vmax, vmax=vmax, aspect="auto", origin="upper")
    for i in range(len(BANDS)):
        for j in range(len(EPOCHS)):
            if np.isnan(M[i, j]):
                cx.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=True, facecolor="0.9", hatch="////", edgecolor="0.6", lw=0))
            else:
                cx.text(j, i, f"{M[i, j]:+.1f}", ha="center", va="center", fontsize=5.2,
                        color="white" if abs(M[i, j]) > 0.6 * vmax else "black")
    cx.set_xticks(range(len(EPOCHS))); cx.set_xticklabels([f"{e:g}" for e in EPOCHS])
    cx.set_yticks(range(len(BANDS))); cx.set_yticklabels(BANDS)
    cx.set_xlabel("time (d)"); cx.set_ylabel("band")
    cb = fig.colorbar(im, ax=cx, fraction=0.05, pad=0.03); cb.set_label(r"$\Delta m$ (C $-$ reference, mag)")
    cx.text(0.5, 1.02, "hatched: below the observable mask", transform=cx.transAxes, ha="center", va="bottom", fontsize=5.5)
    letter(cx, "c")
    return save(fig, "fig1_one_kilonova", out_dir)


# ------------------------------------------------------------------ Fig. 2
def fig2(dest, out_dir):
    """Across the grid."""
    g = _j(dest["grid_table"])
    pts = {tuple(p["point"]): p for p in g["points"]}
    order = _point_order()

    fig = plt.figure(figsize=(W2, 58 * mm))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.35, 0.8], wspace=0.4)
    ax = fig.add_subplot(gs[0]); bx = fig.add_subplot(gs[1]); cx = fig.add_subplot(gs[2])

    # a: worst live colour error per point, by X
    for k, x in enumerate(XS):
        for m in MS:
            for v in VS:
                p = pts[(m, v, x)]
                w = p["C_both"]["worst_dcolor"]; fl = p["C_both"].get("floor_at_cell") or 0.0
                jit = (MS.index(m) - 1) * 0.22
                ax.errorbar(k + jit, w, yerr=fl, fmt=VMK[v], ms=np.sqrt(MSZ[m]) * 1.1, color=XC[x],
                            mfc=XC[x] if m != 0.003 else "white", mew=0.7, ecolor="0.4", elinewidth=0.6, capsize=1.5)
                ax.plot(k + jit, pts[(m, v, x)]["C_binned"]["worst_dcolor"], marker="_", color=LEGC["C_binned"], ms=5, mew=0.9)
    ax.set_xticks(range(3)); ax.set_xticklabels([r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$"])
    ax.set_xlabel(r"$X_{\rm lan}$"); ax.set_ylabel(r"largest live $|\Delta$colour$|$ (mag)")
    ax.set_ylim(0, None)
    h = [Line2D([], [], marker=VMK[v], color="k", ls="", ms=3.5, label=f"$v={v:g}c$") for v in VS]
    h += [Line2D([], [], marker="o", color="k", ls="", ms=np.sqrt(MSZ[m]) * 1.1, mfc="white" if m == 0.003 else "k", label=f"$M={m:g}$") for m in MS]
    h += [Line2D([], [], marker="_", color=LEGC["C_binned"], ls="", ms=5, label="C$_{\\rm bin}$")]
    ax.legend(handles=h, loc="lower right", ncol=2, handletextpad=0.3, columnspacing=0.8)
    ax.text(0.02, 0.97, "bars: A floor at that cell", transform=ax.transAxes, va="top", fontsize=5.5)
    letter(ax, "a")

    # b: Delta(g-r) at 1 d and Delta(J-K) at 2 d for the 27 points
    xs = np.arange(len(order))
    for col, td, c, mk in (("g-r", 1.0, OI["blue"], "o"), ("J-K", 2.0, OI["red"], "s")):
        ys, es = [], []
        for p in order:
            cell = next((cc for cc in _grid_cells(g, p) if cc["t_d"] == td), None)
            y = cell["legs"]["C_both"]["dcolor"].get(col, np.nan) if cell else np.nan
            ys.append(y); es.append(cell["floor"] if cell and np.isfinite(cell["floor"]) else 0.0)
        bx.errorbar(xs, ys, yerr=es, fmt=mk, ms=3, color=c, ecolor="0.5", elinewidth=0.5, capsize=1.2,
                    label=f"$\\Delta$({col}) at {td:g} d")
    for k in range(3):
        bx.axvspan(9 * k - 0.5, 9 * k + 8.5, color=XC[XS[k]], alpha=0.07, lw=0)
    bx.axhline(0, color="k", lw=0.4)
    bx.set_xticks([4, 13, 22]); bx.set_xticklabels([r"$X_{\rm lan}=10^{-3}$", r"$10^{-2}$", r"$10^{-1}$"])
    bx.set_xlim(-0.7, 26.7)
    bx.set_ylabel(r"$\Delta$colour (C $-$ reference, mag)")
    bx.set_xlabel(r"grid point (within each $X_{\rm lan}$: $M$ outer, $v$ inner, ascending)")
    bx.legend(loc="lower left", ncol=2)
    letter(bx, "b", dx=-0.13)

    # c: histogram of the live NIR colour errors
    vals = [c["legs"]["C_both"]["dcolor"][k] for c in _grid_cells(g) for k in NIR_COLS if k in c["legs"]["C_both"]["dcolor"]]
    valsA = [c["legs"]["A_redist"]["dcolor"][k] for c in _grid_cells(g) for k in NIR_COLS if k in c["legs"]["A_redist"]["dcolor"]]
    bins = np.linspace(-4, 1, 26)
    cx.hist(vals, bins=bins, color=LEGC["C_both"], alpha=0.8, label=f"C, $i-J$ and $J-K$ ({len(vals)})")
    cx.hist(valsA, bins=bins, color=LEGC["A_redist"], alpha=0.6, label="A (control)")
    cx.axvline(0, color="k", lw=0.5)
    cx.set_xlabel(r"live NIR $\Delta$colour (mag)"); cx.set_ylabel("count")
    neg = sum(v < 0 for v in vals)
    cx.text(0.03, 0.72, f"{neg}/{len(vals)} negative\n(too blue)", transform=cx.transAxes, va="top", fontsize=6)
    cx.legend(loc="upper left", fontsize=5.5)
    letter(cx, "c", dx=-0.25)
    return save(fig, "fig2_grid", out_dir)


# ------------------------------------------------------------------ Fig. 3
def fig3(dest, out_dir):
    """Not a parameter shift."""
    S = {tg: _j(dest["sensitivity" + ("" if tg == "T0" else f"_{tg}")]) for tg in TANGENTS}
    thr = S["T0"]["thresholds"]
    p0 = S["T0"]["points"][CENTRAL_POINT]["legs"]["C_both"]
    d, f, keys, sig = np.array(p0["d_rt"]), np.array(p0["fit"]), p0["keys"], np.array(p0["sigma"])

    fig = plt.figure(figsize=(W2, 58 * mm))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.1, 1.1], wspace=0.45)
    ax = fig.add_subplot(gs[0]); bx = fig.add_subplot(gs[1]); cx = fig.add_subplot(gs[2])

    # a: d_RT vs best fit at the central point
    seen = set()
    for (b, t), di, fi, si in zip(keys, d, f, sig):
        ax.errorbar(di, fi, xerr=si, fmt="o", ms=3.2, color=BANDC[b], ecolor="0.6", elinewidth=0.5,
                    label=b if b not in seen else None)
        seen.add(b)
    ax.legend(loc="lower right", ncol=2, fontsize=5, handletextpad=0.2, columnspacing=0.6)
    lim = max(np.abs(d).max(), np.abs(f).max()) * 1.15
    ax.plot([-lim, lim], [-lim, lim], color="0.6", lw=0.6, ls="--")
    ax.axhline(0, color="k", lw=0.3); ax.axvline(0, color="k", lw=0.3)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.set_xlabel(r"$d_{\rm RT}$ = C $-$ reference (mag)"); ax.set_ylabel(r"best $(M,v,X)$ tangent fit (mag)")
    ax.text(0.03, 0.97, f"central model\n$R={p0['R']:.2f}$, $\\chi^2_{{\\rm res}}/{{\\rm dof}}={p0['chi2_res_dof']:.0f}$",
            transform=ax.transAxes, va="top", fontsize=6)
    letter(ax, "a")

    # b, c: R and chi2_res/dof per point under T0..T3, by X
    def per_point(tg, field):
        out = {}
        for k, pt in S[tg]["points"].items():
            r = pt["legs"]["C_both"]
            if r.get("status") == "ok" and not r.get("underdetermined"):
                out[_pt(k)] = r[field]
        return out

    for axx, field, thrv, lab in ((bx, "R", thr["R_max"], r"$R = |{\rm residual}|/|d_{\rm RT}|$"),
                                  (cx, "chi2_res_dof", thr["chi2_small"], r"$\chi^2_{\rm res}/{\rm dof}$")):
        for j, tg in enumerate(TANGENTS):
            vals = per_point(tg, field)
            for k, x in enumerate(XS):
                ys = [vals[p] for p in vals if p[2] == x]
                jit = (k - 1) * 0.22
                axx.plot(np.full(len(ys), j + jit) + np.linspace(-0.05, 0.05, len(ys)), ys, ls="", marker="o", ms=2.6,
                         color=XC[x], alpha=0.85, label=XN[x] if j == 0 else None)
        axx.axhline(thrv, color="k", lw=0.6, ls="--")
        axx.set_xticks(range(4)); axx.set_xticklabels([f"{tg}\n{TG_LABEL[tg]}" for tg in TANGENTS])
        axx.set_ylabel(lab)
    bx.set_ylim(0, 1.05); bx.text(3.45, thr["R_max"] + 0.02, "C-A limit", ha="right", fontsize=5.5)
    cx.set_yscale("log"); cx.text(3.45, thr["chi2_small"] * 1.2, "C-C limit", ha="right", fontsize=5.5)
    bx.legend(loc="lower left", fontsize=5.5)
    bx.set_xlabel("tangent space"); cx.set_xlabel("tangent space")
    letter(bx, "b"); letter(cx, "c")
    return save(fig, "fig3_not_parameters", out_dir)


# ------------------------------------------------------------------ Fig. 4
def fig4(dest, out_dir):
    """Observability and structure."""
    ob, sy = _j(dest["observability"]), _j(dest["syserr"])
    fig = plt.figure(figsize=(W2, 60 * mm))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.5], wspace=0.3)
    ax = fig.add_subplot(gs[0]); bx = fig.add_subplot(gs[1])

    # a: Gate 3 survival per scenario x tangent, by X
    sm = ob["summary"]["C_both"]
    width = 0.26
    for i, scn in enumerate(SCENARIOS):
        for j, tg in enumerate(TANGENTS):
            x0 = i * 5 + j
            for k, x in enumerate(XS):
                bx_ = sm[scn][tg]["by_X"][f"{x:g}"]
                frac = bx_["survives"] / bx_["eligible"] if bx_["eligible"] else 0
                ax.bar(x0 + (k - 1) * width, frac, width, color=XC[x], alpha=0.9 if tg == "T0" else 0.55,
                       label=XN[x] if (i == 0 and j == 0) else None)
                ax.text(x0 + (k - 1) * width, frac + 0.02, f"{bx_['survives']}/{bx_['eligible']}", ha="center", va="bottom",
                        fontsize=3.8, rotation=90)
    ax.set_xticks([i * 5 + j for i in range(3) for j in range(4)])
    ax.set_xticklabels([tg for _ in range(3) for tg in TANGENTS], fontsize=5.5)
    for i, scn in enumerate(SCENARIOS):
        ax.text(i * 5 + 1.5, 1.22, scn, ha="center", fontsize=6.5)
    ax.set_ylim(0, 1.55); ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("fraction of eligible points\nwhere the residual survives")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=5.5, handlelength=1.0, columnspacing=0.7)
    letter(ax, "a", dx=-0.22)

    # b: sigma_sys comparison
    leg = sy["legs"]["C_both"]
    tt = np.array(EPOCHS)
    bx.axhspan(-1, 1, color="0.85", lw=0, label=r"$\sigma_{\rm sys}=1$ mag")
    bx.axhspan(-0.5, 0.5, color="0.7", lw=0, label=r"$\sigma_{\rm sys}=0.5$ mag")
    minn = N_GREY
    for b in BANDS:
        rows = {r["t_d"]: r for r in leg["per_key"] if r["band"] == b}
        for td in EPOCHS:
            r = rows.get(td)
            if r is None:
                continue
            grey = r["n"] < minn
            c = "0.5" if grey else BANDC[b]
            off = (BANDS.index(b) - 3) * 0.02
            bx.errorbar(td * (1 + off), r["median"], yerr=[[r["median"] - r["p16"]], [r["p84"] - r["median"]]],
                        fmt="o", ms=2.8, color=c, ecolor=c, elinewidth=0.6, capsize=1.2, alpha=0.5 if grey else 1)
        good = [(td, rows[td]["median"]) for td in EPOCHS if td in rows and rows[td]["n"] >= minn]
        if good:
            bx.plot(*zip(*good), color=BANDC[b], lw=0.7, label=b)
    bx.axhline(0, color="k", lw=0.4)
    logx(bx, EPOCHS); bx.set_ylim(-2.4, 2.4)
    bx.set_xlabel("time (d)"); bx.set_ylabel(r"$d_{\rm RT}$ (mag): grid median, 16$-$84%")
    om, nl = leg["one_mode"]["f1"], leg["null_sign_scramble"]
    nlm = nl["median"] if isinstance(nl, dict) else nl
    bx.text(0.02, 0.03, f"one coherent mode: {om:.2f} of $\\|d_{{\\rm RT}}\\|^2$\n(sign-scrambled null {nlm:.2f}; control {sy['legs']['A_redist']['one_mode']['f1']:.2f})\ngrey: fewer than {minn} live points",
            transform=bx.transAxes, va="bottom", fontsize=5.5)
    bx.legend(loc="upper right", ncol=5, fontsize=5.5, handlelength=1.4, columnspacing=0.8)
    letter(bx, "b", dx=-0.12)
    return save(fig, "fig4_observability_syserr", out_dir)


# --------------------------------------------------------------- ED Fig. 1
def edfig1(dest, out_dir):
    """Source model and gate: L, T_eff, v_ph/v_ej, and the line strength S."""
    g = _j(dest["grid_table"])
    fig, axs = plt.subplots(1, 4, figsize=(W2, 42 * mm))
    fig.subplots_adjust(wspace=0.5)
    t = np.geomspace(0.3, 10, 200) * DAY
    for m in MS:
        for v in VS:
            sm = SourceModel(m, v)
            lw = 1.2 if (m, v) == (0.01, 0.1) else 0.6
            c = {0.003: OI["blue"], 0.01: OI["purple"], 0.03: OI["red"]}[m]
            ls = {0.05: ":", 0.1: "-", 0.2: "--"}[v]
            axs[0].plot(t / DAY, sm.luminosity(t), color=c, ls=ls, lw=lw)
            axs[1].plot(t / DAY, sm.t_eff(t), color=c, ls=ls, lw=lw)
            axs[2].plot(t / DAY, sm.v_ph(t), color=c, ls=ls, lw=lw, label=f"$M={m:g}$, $v={v:g}c$")
    axs[2].axhline(0.5, color="k", lw=0.5); axs[2].text(0.35, 0.515, "floor (epochs excluded)", fontsize=5.5)
    for x in XS:
        cells = [c for c in g["cells"] if c["ran"] and c["point"][2] == x and c["S"] is not None]
        axs[3].plot([c["t_d"] for c in cells], [c["S"] for c in cells], ls="", marker="o", ms=2, color=XC[x], label=XN[x], alpha=0.8)
    for a, yl, log in ((axs[0], r"$L_{\rm bol}$ (erg s$^{-1}$)", True), (axs[1], r"$T_{\rm eff}$ (K)", True),
                       (axs[2], r"$v_{\rm ph}/v_{\rm ej}$", False), (axs[3], r"line strength $S$ (band max)", True)):
        logx(a, [0.5, 1, 2, 5]); a.set_xlabel("time (d)"); a.set_ylabel(yl)
        if log:
            a.set_yscale("log")
        for e in EPOCHS:
            a.axvline(e, color="0.85", lw=0.4, zorder=0)
    fig.legend(*axs[2].get_legend_handles_labels(), fontsize=5.5, ncol=5, loc="lower center",
               bbox_to_anchor=(0.5, -0.3), handlelength=1.8, columnspacing=0.9)
    axs[3].legend(fontsize=5, loc="upper right")
    for a, s in zip(axs, "abcd"):
        letter(a, s, dx=-0.3)
    return save(fig, "edfig1_source", out_dir)


# --------------------------------------------------------------- ED Fig. 2
def edfig2(dest, out_dir):
    """T_eff validation at the central point: MC direction vs Planck proxy."""
    ts = _j(dest["tscale"])
    fl = ts["noise_floor"]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 55 * mm))
    fig.subplots_adjust(wspace=0.35)
    for a, (tag, title) in zip(axs, (("illumination_only", "illumination only"), ("with_T_gas", r"with $T_{\rm gas}$"))):
        cc = ts["variants"][tag]["compare"]["C_both"]
        for key, mc in cc["mc"].items():
            b, td = key.split(","); pr = cc["proxy"][key]
            a.plot(pr, mc, ls="", marker="o", ms=3.5, color=BANDC[b])
        lim = max(max(abs(v) for v in cc["mc"].values()), max(abs(v) for v in cc["proxy"].values())) * 1.15
        a.plot([-lim, lim], [-lim, lim], color="0.6", lw=0.6, ls="--")
        a.axhspan(-fl, fl, color="0.85", lw=0)
        a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect("equal")
        a.set_xlabel(r"Planck proxy $\partial m/\partial\ln T$ at $T_{\rm eff}$ (mag)")
        a.set_ylabel(r"MC $d_T$ (mag per e-fold)")
        a.text(0.03, 0.97, f"{title}\ncos $= {cc['cos']:.2f}$, norm ratio $= {cc['norm_ratio']:.2f}$\n$N={cc['N']}$; grey: A floor",
               transform=a.transAxes, va="top", fontsize=6)
    axs[1].legend(handles=[Line2D([], [], marker="o", ls="", color=BANDC[b], ms=3, label=b) for b in BANDS],
                  ncol=2, fontsize=5, loc="lower right")
    for a, s in zip(axs, "ab"):
        letter(a, s, dx=-0.2)
    return save(fig, "edfig2_tscale", out_dir)


# --------------------------------------------------------------- ED Fig. 3
def edfig3(dest, out_dir):
    """Chain cap: per-band closure error vs cap at the four worst-trapped cells."""
    ct = _j(dest["chain_table"])
    cells = ct["cells"]
    fig, axs = plt.subplots(1, len(cells), figsize=(W2, 48 * mm), sharey=True)
    fig.subplots_adjust(wspace=0.15)
    for a, cell in zip(np.atleast_1d(axs), cells):
        caps = sorted(int(k) for k in cell["runs"])
        for b in BANDS:
            ys = [cell["runs"][str(c)]["legs"]["C_both"]["dm"].get(b, np.nan) for c in caps]
            a.plot(caps, ys, marker="o", ms=2.5, color=BANDC[b], label=b)
        a.axhline(0, color="k", lw=0.4)
        a.set_xscale("log", base=2); a.set_xticks(caps); a.set_xticklabels([str(c) for c in caps])
        a.set_xlabel("chain cap")
        tr = [cell["runs"][str(c)]["trapped_frac"] for c in caps]
        m, v, x = cell["point"]
        a.text(0.5, 0.97, f"$M={m:g}$, $v={v:g}c$, $X={x:g}$, {cell['t_d']:g} d\ntrapped {100*tr[0]:.1f}% $\\to$ {100*tr[-1]:.1f}%",
               transform=a.transAxes, ha="center", va="top", fontsize=5.5)
        a.set_ylim(-1.15, 1.25)
        crit = cell["runs"][str(caps[-1])]["legs"]["C_both"]
        a.text(0.03, 0.03, f"max change {crit['max_dm_change']:.2f} mag", transform=a.transAxes, fontsize=5.5)
    axs[0].set_ylabel(r"$\Delta m$ (C $-$ reference, mag)")
    axs[1].legend(loc="center right", ncol=2, fontsize=5, handlelength=1.2, columnspacing=0.6)
    for a, s in zip(axs, "abcd"):
        letter(a, s, dx=-0.12 if s != "a" else -0.3, dy=1.02)
    return save(fig, "edfig3_chain_cap", out_dir)


# --------------------------------------------------------------- ED Fig. 4
def edfig4(dest, out_dir):
    """The A_redist control across the grid: per-band error vs the cell floor, by n_used."""
    g = _j(dest["grid_table"])
    nmin = g["floor"]["n_min"]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 52 * mm), sharey=True)
    fig.subplots_adjust(wspace=0.12)
    a, b_ = axs
    for c in _grid_cells(g):
        well = c["n_used"] >= nmin
        for b in c["live"]:
            a.plot(c["n_used"], abs(c["legs"]["A_redist"]["dm"][b]), ls="", marker="o", ms=2.2, color=BANDC[b],
                   alpha=0.8 if well else 0.35)
            b_.plot(c["n_used"], abs(c["legs"]["C_both"]["dm"][b]), ls="", marker="o", ms=2.2, color=BANDC[b],
                    alpha=0.8 if well else 0.35)
    for ax_, lab in ((a, "A (grouped $R_{ij}$, the control)"), (b_, "C (both approximations)")):
        logx(ax_, [2e4, 5e4, 1e5, 2e5, 3e5], ["20k", "50k", "100k", "200k", "300k"]); ax_.set_yscale("log")
        ax_.axvline(nmin, color="k", lw=0.5, ls="--")
        ax_.axhline(g["floor"]["median"], color=LEGC["A_redist"], lw=0.7)
        ax_.axhline(g["floor"]["p90"], color=LEGC["A_redist"], lw=0.7, ls=":")
        ax_.set_xlabel(r"packets per seed $n_{\rm used}$"); ax_.set_ylabel(r"$|\Delta m|$ (mag)")
        ax_.text(0.03, 0.97, lab, transform=ax_.transAxes, va="top", fontsize=6)
    a.text(nmin * 1.1, g["floor"]["median"] * 1.15, "floor median / 90%", fontsize=5, color=LEGC["A_redist"])
    b_.set_ylabel("")
    b_.legend(handles=[Line2D([], [], marker="o", ls="", color=BANDC[b], ms=3, label=b) for b in BANDS], ncol=4, fontsize=5, loc="lower right")
    letter(a, "a", dx=-0.2); letter(b_, "b", dx=-0.2)
    return save(fig, "edfig4_control_floor", out_dir)


ITEMS = {"fig1": fig1, "fig2": fig2, "fig3": fig3, "fig4": fig4,
         "edfig1": edfig1, "edfig2": edfig2, "edfig3": edfig3, "edfig4": edfig4}
# file stem each item writes (the manuscript includes them without extension)
NAMES = {"fig1": "fig1_one_kilonova", "fig2": "fig2_grid", "fig3": "fig3_not_parameters",
         "fig4": "fig4_observability_syserr", "edfig1": "edfig1_source", "edfig2": "edfig2_tscale",
         "edfig3": "edfig3_chain_cap", "edfig4": "edfig4_control_floor"}


def main(dest, out_dir, which=None):
    out_dir = Path(out_dir)
    written = []
    for name, fn in ITEMS.items():
        if which and name not in which:
            continue
        written += fn(dest, out_dir)
        print("wrote", name)
    return written


if __name__ == "__main__":
    import argparse
    sys.path.insert(0, str(ROOT / "paper3"))
    import freeze
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "figures"))
    ap.add_argument("--which", nargs="*")
    a = ap.parse_args()
    main(freeze.canonical(), a.out, a.which)
