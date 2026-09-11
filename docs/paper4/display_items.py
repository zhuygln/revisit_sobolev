"""Paper IV display items, from the committed run records only (the same
files paper4/freeze.py reads): Fig. 1 the closure/redistribution matrix
(F62), Fig. 2 the Fontes snapshot (F64), Fig. 3 the Fontes light curve
(F65), Fig. 4 the macroatom test (F63), Fig. 5 the P1-xkn light curve (F66,
written when its record exists). Tables are written by latex_tables.py.

PDF + PNG into docs/paper4/figures/ with the metadata dates stripped, so the
output is byte-stable within one environment (SOURCE_DATE_EPOCH set by the
Makefile).

    main(out_dir, which=None) -> list of written paths
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BANDS = "grizJHK"
DAY = 86400.0
P1XKN = ROOT / "paper4/phase10_fontes/prod_record/p1xkn/summary.json"


def save(fig, name, out_dir):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}),
                      ("png", {"Software": None})):
        p = out_dir / f"{name}.{ext}"
        fig.savefig(p, metadata=meta, bbox_inches="tight", pad_inches=0.02, dpi=150)
        written.append(p)
    plt.close(fig)
    return written


def fig_matrix_p1(out_dir):
    """F62: closure - resolved per band on P1 (1-5 d), P2 and the robustness
    patterns, under eps = 1 and under fluorescence, from the corrected legs."""
    states = [("P1 1 d", "legs_P1_t1"), ("P1 2 d", "legs_P1_t2"), ("P1 3 d", "legs_P1_t3"), ("P1 5 d", "legs_P1_t5"),
              ("P1r1", "legs_P1r1_t2"), ("P1r2", "legs_P1r2_t2"), ("P2 3.4 d", "legs_P2_t3.4")]
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True, sharey="row")
    pairs = [(("Bth", "Rth"), "expansion, ε = 1"), (("Bbinth", "Rth"), "line-binned, ε = 1"),
             (("B2", "R2"), "expansion, fluorescence"), (("Bbin2", "R2"), "line-binned, fluorescence")]
    for ax, ((a, r), title) in zip(axes.T.flatten(), pairs):
        for name, f in states:
            d = json.load(open(ROOT / "paper4/phase9_final" / f"{f}.json"))["legs"]
            dm = [d[a]["mags"][b] - d[r]["mags"][b] for b in BANDS]
            ax.plot(range(7), dm, "o-", label=name, lw=1.2, ms=4)
        ax.axhline(0, color="k", lw=0.6); ax.set_title(title, fontsize=10)
        ax.set_xticks(range(7)); ax.set_xticklabels(list(BANDS))
    axes[0, 0].set_ylabel("closure − resolved [mag]"); axes[1, 0].set_ylabel("closure − resolved [mag]")
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.tight_layout(); return save(fig, "fig_f62_matrix", out_dir)


def fig_fontes_snapshot(out_dir):
    """F64: the Fontes 4-d snapshot, both brackets, both line lists."""
    runs = [("GSI, reemit", "fontes_t4_z45_reemit"), ("GSI, reflect", "fontes_t4_z45_reflect"),
            ("JPLT, reemit", "fontes_t4_z45_jplt_reemit"), ("JPLT, reflect", "fontes_t4_z45_jplt_reflect")]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, (pairs, title) in zip(axes, [([("Bth", "Rth", "expansion"), ("Bbinth", "Rth", "line-binned")], "ε = 1"),
                                         ([("B2", "R2", "expansion"), ("Bbin2", "R2", "line-binned")], "fluorescence")]):
        for name, f in runs:
            p = ROOT / "paper4/phase10_fontes" / f"{f}.json"
            if not p.exists():
                continue
            d = json.load(open(p))["legs"]
            for a, r, lab in pairs:
                if a in d and r in d:
                    dm = [d[a]["mags"][b] - d[r]["mags"][b] for b in BANDS]
                    ax.plot(range(7), dm, ("o-" if lab == "expansion" else "s--"), lw=1, ms=4, label=f"{lab}, {name}")
        ax.axhline(0, color="k", lw=0.6); ax.set_title(f"Fontes 4 d snapshot, {title}", fontsize=10)
        ax.set_xticks(range(7)); ax.set_xticklabels(list(BANDS)); ax.set_ylim(-1.6, 4.0)
    axes[0].set_ylabel("closure − resolved [mag]"); axes[1].legend(fontsize=6, ncol=2)
    fig.tight_layout(); return save(fig, "fig_f64_snapshot", out_dir)


def fig_lightcurve(out_dir):
    """F65: the six light curves, the band residuals and the colour residuals."""
    s = json.load(open(ROOT / "paper4/phase10_fontes/prod_record/merged/summary.json"))
    t = np.array(s["t_grid"]) / DAY; tm = np.sqrt(t[1:] * t[:-1])
    tp = np.array(s["t_phot"]) / DAY; tpm = np.sqrt(tp[1:] * tp[:-1])
    style = {"Rth": ("k", "-", "resolved, ε = 1"), "Bth": ("tab:blue", "-", "expansion, ε = 1"), "Bbinth": ("tab:red", "-", "line-binned, ε = 1"),
             "R2": ("k", "--", "resolved, fluorescence"), "B2": ("tab:blue", "--", "expansion, fluorescence"), "Bbin2": ("tab:red", "--", "line-binned, fluorescence")}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for leg, o in s["legs"].items():
        c, ls, lab = style[leg]
        axes[0].plot(tm, np.array(o["L_esc"]) / 1e40, color=c, ls=ls, label=lab)
    axes[0].set_xlabel("t [d]"); axes[0].set_ylabel("L_esc [10⁴⁰ erg s⁻¹]"); axes[0].legend(fontsize=7); axes[0].set_title("bolometric light curves", fontsize=10)
    for leg, o in s["legs"].items():
        if "dm_vs_ref" not in o:
            continue
        c, ls, lab = style[leg]
        axes[1].plot(tpm, [d["z"] for d in o["dm_vs_ref"]], color=c, ls=ls, marker="o", ms=3, label=lab + " (z)")
        axes[1].plot(tpm, [d["H"] for d in o["dm_vs_ref"]], color=c, ls=ls, marker="^", ms=3, alpha=0.6, label=lab + " (H)")
        axes[2].plot(tpm, [d["J-K"] for d in o["dcolour_vs_ref"]], color=c, ls=ls, marker="o", ms=3, label=lab)
    for ax in axes[1:]:
        ax.axhline(0, color="k", lw=0.6); ax.set_xlabel("t [d]")
    axes[1].set_ylabel("closure − resolved [mag]"); axes[1].set_title("band residuals (z circles, H triangles)", fontsize=10); axes[1].legend(fontsize=6, ncol=2)
    axes[2].set_ylabel("Δ(J−K) [mag]"); axes[2].set_title("colour residual", fontsize=10); axes[2].legend(fontsize=7)
    fig.tight_layout(); return save(fig, "fig_f65_lightcurve", out_dir)


def fig_macroatom(out_dir):
    """F63: the full macroatom check on P1 2 d."""
    d = json.load(open(ROOT / "paper4/phase10_fontes/macro_P1_t2_s28.json"))["legs"]
    fig, ax = plt.subplots(figsize=(6, 4))
    for a, r, lab, st in (("B2", "R2", "expansion, downward", "o-"), ("Bbin2", "R2", "line-binned, downward", "s-"),
                          ("B2M", "R2M", "expansion, full macroatom", "o--"), ("Bbin2M", "R2M", "line-binned, full macroatom", "s--")):
        ax.plot(range(7), [d[a]["mags"][b] - d[r]["mags"][b] for b in BANDS], st, ms=4, label=lab)
    ax.axhline(0, color="k", lw=0.6); ax.set_xticks(range(7)); ax.set_xticklabels(list(BANDS)); ax.set_ylabel("closure − resolved [mag]")
    ax.set_title("P1, 2 d", fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); return save(fig, "fig_f63_macroatom", out_dir)


def fig_p1xkn(out_dir):
    """F66: the P1-xkn transport replacement test, fluorescence trio, with
    xkn's own L_thick + L_thin for reference."""
    if not P1XKN.exists():
        return []
    s = json.load(open(P1XKN))
    t = np.array(s["t_grid"]) / DAY; tm = np.sqrt(t[1:] * t[:-1])
    tp = np.array(s["t_phot"]) / DAY; tpm = np.sqrt(tp[1:] * tp[:-1])
    style = {"R2": ("k", "resolved"), "B2": ("tab:blue", "expansion"), "Bbin2": ("tab:red", "line-binned")}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for leg, o in s["legs"].items():
        c, lab = style[leg]
        axes[0].plot(tm, np.array(o["L_esc"]) / 1e40, color=c, label=lab)
    if "L_xkn" in s:
        axes[0].plot(tm, np.array(s["L_xkn"]) / 1e40, color="grey", ls=":", label="xkn L_thick + L_thin")
    axes[0].set_xlabel("t [d]"); axes[0].set_ylabel("L_esc [10⁴⁰ erg s⁻¹]"); axes[0].legend(fontsize=7)
    axes[0].set_title("P1-xkn: bolometric, fluorescence", fontsize=10)
    for leg, o in s["legs"].items():
        if "dm_vs_ref" not in o:
            continue
        c, lab = style[leg]
        axes[1].plot(tpm, [d["z"] for d in o["dm_vs_ref"]], color=c, marker="o", ms=3, label=lab + " (z)")
        axes[1].plot(tpm, [d["K"] for d in o["dm_vs_ref"]], color=c, marker="^", ms=3, alpha=0.6, label=lab + " (K)")
        axes[2].plot(tpm, [d["J-K"] for d in o["dcolour_vs_ref"]], color=c, marker="o", ms=3, label=lab)
    for ax in axes[1:]:
        ax.axhline(0, color="k", lw=0.6); ax.set_xlabel("t [d]")
    axes[1].set_ylabel("closure − resolved [mag]"); axes[1].set_title("band residuals (z circles, K triangles)", fontsize=10); axes[1].legend(fontsize=6, ncol=2)
    axes[2].set_ylabel("Δ(J−K) [mag]"); axes[2].set_title("colour residual", fontsize=10); axes[2].legend(fontsize=7)
    fig.tight_layout(); return save(fig, "fig_f66_p1xkn", out_dir)


ITEMS = {"fig1": fig_matrix_p1, "fig2": fig_fontes_snapshot, "fig3": fig_lightcurve, "fig4": fig_macroatom, "fig5": fig_p1xkn}


def main(out_dir, which=None):
    written = []
    for name, fn in ITEMS.items():
        if which and name not in which:
            continue
        w = fn(out_dir)
        written += w
        print("wrote" if w else "skipped (no record)", name)
    return written


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "figures"))
    ap.add_argument("--which", nargs="*")
    a = ap.parse_args()
    main(a.out, a.which)
