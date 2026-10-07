#!/usr/bin/env python3
"""The Paper B display items, generated from paperB/FROZEN.json alone (the
same frozen record the manuscript's macros come from; no number is typed
here). The PI froze these two on 2026-09-22.

  (Revised 2026-10-05, the PI: the discovery is the first figure.)

  Figure 1 -- the flagship. (a) a scalar closure fails where a small
      operator succeeds (G1); (b-d) better on the microscopic events, worse
      on the light (G2). Figure 2 -- generality (G3). End Matter figure --
      the operator follows the fluorescence physics it is built from (the
      R2M robustness check). The former standalone figures are kept
      reproducible below.

  Former Figure 1 -- Compression survives the physics.  For La II, Ce II and Nd II:
      the optimally tuned scalar eps* and the group operators R_2 ... R_32
      against the energy-conserving downward macroatom, with the R2M
      robustness check (the same operator rebuilt from the
      radiation-field-driven macroatom's events) on the right.
      The message: eps* fails where a small R succeeds, and R follows the
      reference physics when that physics changes.

  Figure 3 -- Transfer is not compression.  Per ion and state axis, the
      anchor operator transported unchanged, the operator rebuilt at the
      state, and the whole-operator interpolation at the interior point,
      with the preregistered A/B/C/D reading; G3 (F70).

  Figure 2 -- Why the compression works.  Transport error and event-level
      error against the number of archetypal exit distributions, for local
      frequency coarsening and for the global non-negative factorisation;
      and the scatter of m_event against m_band, where the central
      phenomenon is visible directly: the global points move LEFT (they fit
      the microscopic events better) while staying HIGH (they reproduce the
      light worse).

PDF + PNG with the metadata dates stripped, so the output is byte-stable.

    .venv/bin/python docs/paperB/figures.py [--out DIR]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = ROOT / "paperB/FROZEN.json"
COL = {"57LaII": "#0072B2", "58CeII": "#E69F00", "60NdII": "#D55E00"}
MARK = {"57LaII": "o", "58CeII": "s", "60NdII": "^"}


def save(fig, out_dir, stem):
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}),
                      ("png", {"Software": None})):
        p = out_dir / f"{stem}.{ext}"
        fig.savefig(p, metadata=meta, dpi=200, bbox_inches="tight")
        written.append(p)
    plt.close(fig)
    return written


def _panel_existence(ax, h):
    names, dm_max = h["ion_names"], h["g1"]["dm_max"]
    for ion in h["ions"]:
        r = h["g1"]["per_ion"][ion]; c = COL[ion]
        ng = [t["ng"] for t in r["table"]]; band = [t["band_max"] for t in r["table"]]
        ax.plot(ng, band, MARK[ion] + "-", color=c, label=f"{names[ion]}: $R_{{N_g}}$", ms=5)
        ax.axhline(r["eps_max_dm"], color=c, ls="--", lw=1.1)
        ax.annotate(rf"$\epsilon^\star$, {names[ion]}", (ng[-1], r["eps_max_dm"]), textcoords="offset points",
                    xytext=(-4, 4), fontsize=7, color=c, ha="right")
    ax.axhline(dm_max, color="k", ls=":", lw=1.2)
    ax.annotate("criterion", (2, dm_max), textcoords="offset points", xytext=(2, 4), fontsize=7)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks([2, 4, 8, 16, 32]); ax.set_xticklabels(["2", "4", "8", "16", "32"])
    ax.set_xlabel("frequency groups $N_g$")
    ax.set_ylabel(r"max $|\Delta m|$ over live bands  [mag]")
    ax.legend(fontsize=7, loc="lower left")


def _panels_inversion(axes, h, legend=True):
    names, dm_max = h["ion_names"], h["g1"]["dm_max"]
    for ion in h["ions"]:
        r = h["g2"]["per_ion"][ion]; c = COL[ion]
        loc, glo = r["local"], r["glob"]
        axes[0].plot([e["k"] for e in loc], [e["band"] for e in loc], MARK[ion] + "-", color=c, ms=5, label=f"{names[ion]}: local")
        axes[0].plot([e["k"] for e in glo], [e["band"] for e in glo], MARK[ion] + "--", color=c, ms=5, mfc="none", label=f"{names[ion]}: global")
        axes[1].plot([e["k"] for e in loc], [e["event"] for e in loc], MARK[ion] + "-", color=c, ms=5)
        axes[1].plot([e["k"] for e in glo], [e["event"] for e in glo], MARK[ion] + "--", color=c, ms=5, mfc="none")
        axes[2].plot([e["event"] for e in loc], [e["band"] for e in loc], MARK[ion] + "-", color=c, ms=5)
        axes[2].plot([e["event"] for e in glo], [e["band"] for e in glo], MARK[ion] + "--", color=c, ms=5, mfc="none")
    for ax in (axes[0], axes[2]):
        ax.axhline(dm_max, color="k", ls=":", lw=1.2)
    for ax in axes[:2]:
        ax.set_xscale("log", base=2); ax.set_xticks([1, 2, 4, 8, 16, 32])
        ax.set_xticklabels(["1", "2", "4", "8", "16", "32"])
        ax.set_xlabel("archetypal exit distributions")
    axes[0].set_yscale("log"); axes[0].set_ylabel(r"max $|\Delta m|$  [mag]")
    if legend:
        axes[0].legend(fontsize=6, ncol=2)
    axes[1].set_ylabel("event-level total-variation loss")
    axes[2].set_yscale("log"); axes[2].set_xlabel("event-level total-variation loss")
    axes[2].set_ylabel(r"max $|\Delta m|$  [mag]")
    axes[2].annotate("", xy=(0.08, 0.88), xytext=(0.40, 0.88), xycoords="axes fraction", textcoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", color="0.4", lw=1.1))
    axes[2].annotate("the global family fits the events better", (0.10, 0.91), xycoords="axes fraction", fontsize=7, color="0.35")
    axes[2].annotate("passes", (0.02, 0.10), xycoords="axes fraction", fontsize=7, color="0.35")


def figure1(h, out_dir):
    """The flagship (the PI, 2026-10-05): the discovery is the first figure.
    (a) a scalar closure fails where a small operator succeeds (G1);
    (b-d) better on the microscopic events, worse on the light (G2)."""
    # two by two: at PRL full width (about 7 in) a 1 x 4 strip renders its
    # labels below 4 pt; this layout keeps them near 7 pt
    fig, ax2 = plt.subplots(2, 2, figsize=(9.6, 7.6))
    axes = [ax2[0][0], ax2[0][1], ax2[1][0], ax2[1][1]]
    _panel_existence(axes[0], h)
    _panels_inversion(axes[1:], h)
    axes[0].set_title("a. a scalar closure fails; a small operator succeeds", fontsize=9, loc="left")
    axes[1].set_title("b. band error (solid: local, open: global)", fontsize=9, loc="left")
    axes[2].set_title("c. microscopic error: the global family wins at high count", fontsize=9, loc="left")
    axes[3].set_title("d. at high count: better on the events, worse on the light", fontsize=9, loc="left")
    fig.tight_layout()
    return save(fig, out_dir, "fig1_flagship")


def figure_robustness(h, out_dir):
    """End Matter: the operator follows the fluorescence physics it is built
    from (the R2M robustness check, F68)."""
    names, dm_max = h["ion_names"], h["g1"]["dm_max"]
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    x = np.arange(len(h["ions"])); w = 0.28
    ref_shift = [abs(h["r2m"]["per_ion"][i]["shift_max"]) for i in h["ions"]]
    rebuilt = [h["r2m"]["per_ion"][i]["rebuilt_band"] for i in h["ions"]]
    downward = [h["r2m"]["per_ion"][i]["downward_band"] for i in h["ions"]]
    ax.bar(x - w, ref_shift, w, color="0.55", label="the reference itself moves")
    ax.bar(x, downward, w, color="#CC79A7", label=r"$R_8$ trained on the old reference")
    ax.bar(x + w, rebuilt, w, color="#009E73", label=r"$R_8$ rebuilt on the new reference")
    ax.axhline(dm_max, color="k", ls=":", lw=1.2)
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([names[i] for i in h["ions"]])
    ax.set_ylabel(r"max $|\Delta m|$ vs the new reference  [mag]")
    ax.set_ylim(top=max(ref_shift) * 6)
    ax.legend(fontsize=7, loc="lower left", framealpha=0.95)
    fig.tight_layout()
    return save(fig, out_dir, "figEM_robustness")


def figure2(h, out_dir):
    """Kept for the record: the three inversion panels alone (the former
    Figure 2); the Letter now carries them inside Figure 1."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.1))
    _panels_inversion(axes, h)
    fig.tight_layout()
    return save(fig, out_dir, "fig2_mechanism")


AXLAB = {"T": r"$\log T_{\rm gas}$", "D": r"$\log n_{\rm ion}$", "J": r"$\log T_{\rm core}$", "P": r"$\log t$"}
CLS = {"A": "#009E73", "B": "#0072B2", "C": "#E69F00", "D": "#D55E00", "GRAY": "0.6"}


def figure3(h, out_dir):
    """The Letter's Figure 2 -- transfer is not compression. Per ion (rows)
    and state axis (columns): the anchor operator transported unchanged
    (filled squares), the operator rebuilt at the state (open circles), and
    at the interior point the whole-operator interpolation (diamonds), each
    against the state's own reference; the letter is the preregistered
    reading (A transfers, B interpolates, C only a fresh fit passes, D none).
    Revised 2026-10-06 (the referee, the PI): the y value is the JOINT
    decision statistic max(max |dm|, max |dcolour|), so the letters and the
    dotted threshold refer to the same quantity (three Ce/Nd readings are
    decided by colour alone and were invisible on a band-only axis)."""
    g3 = h["g3"]; names = h["ion_names"]; dm_max = h["g1"]["dm_max"]; ng_t = g3["ng_t"]
    ions = [i for i in ("58CeII", "60NdII", "57LaII") if i in g3["per_ion"]]
    axes_names = g3["axes"]
    # sized for PRL full width: 12 panels at (3.0 x 2.5) in each render at ~0.6 scale,
    # which keeps tick labels near 6 pt (the earlier 3.6 x 2.9 fell below 5 pt)
    fig, axes = plt.subplots(len(ions), len(axes_names), figsize=(3.0 * len(axes_names), 2.5 * len(ions)), squeeze=False, sharey="row")
    for ri, ion in enumerate(ions):
        sts = g3["per_ion"][ion]["states"]
        for ci, ax_name in enumerate(axes_names):
            ax = axes[ri][ci]
            for st in sorted((x for x in sts if x["axis"] == ax_name), key=lambda x: x["coord_value"]):
                c = st["coord_value"]; col = CLS[st["cls"]]
                # the y value is the joint decision statistic max(band, colour): the one the letters are read on
                ax.plot([c], [st["fresh_joint"]], "o", color=col, ms=6, mfc="none", mew=1.4)
                if st["anchor_joint"] is not None:
                    ax.plot([c], [st["anchor_joint"]], "s", color=col, ms=6)
                if st["interp_whole_joint"] is not None:
                    ax.plot([c], [st["interp_whole_joint"]], "D", color=col, ms=6)
                    ax.plot([c], [st["interp_matrix_joint"]], "d", color=col, ms=5, mfc="none")
                # the paired-seed 95 % interval of G3U where the reading was near the threshold (annotates; the frozen point stays)
                for fam, dx in (("fresh", -0.03), ("anchor", 0.03), ("interp_whole", -0.03), ("interp_matrix", 0.03)):
                    uu = st.get(f"{fam}_u")
                    if uu:
                        ucol = "0.45" if uu["status"] == "GRAY" else col
                        ax.plot([c + dx, c + dx], [uu["lo"], uu["hi"]], "-", color=ucol, lw=1.0, alpha=0.9)
                        ax.plot([c + dx], [uu["e"]], "_", color=ucol, ms=5, mew=1.0)
                top = max(st["fresh_joint"], st["anchor_joint"] or 0.0, st["interp_whole_joint"] or 0.0,
                          *[st[f"{f}_u"]["hi"] for f in ("fresh", "anchor", "interp_whole", "interp_matrix") if st.get(f"{f}_u")])
                ax.annotate(st["cls"], (c, top), textcoords="offset points", xytext=(0, 5), ha="center", fontsize=9, color=col)
            ax.axhline(dm_max, color="k", ls=":", lw=1)
            ax.set_yscale("log"); ax.set_ylim(2e-3, 1.5)
            ax.tick_params(labelsize=9)
            if ri == len(ions) - 1:
                ax.set_xlabel(AXLAB[ax_name])
            if ci == 0:
                ax.set_ylabel(f"{names[ion]}\nmax(band, colour) error [mag]", fontsize=8)
            if ri == 0:
                ax.set_title({"T": "gas temperature", "D": "density", "J": "source spectrum", "P": "trajectory (all coordinates)"}[ax_name], fontsize=10)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker="s", color="0.3", ls="", ms=6, label=f"anchor $R_{{{ng_t}}}$ transported"),
          Line2D([], [], marker="o", color="0.3", ls="", ms=6, mfc="none", mew=1.4, label=f"fresh $R_{{{ng_t}}}$ built at the state"),
          Line2D([], [], marker="D", color="0.3", ls="", ms=6, label="whole-operator interpolation"),
          Line2D([], [], marker="d", color="0.3", ls="", ms=5, mfc="none", label="matrix-only interpolation (diagnostic)"),
          Line2D([], [], marker="_", color="0.3", ls="-", lw=1.0, ms=5, label="paired-seed 95 % interval, 12 seeds (near-threshold cases)")]
    fig.legend(handles=hs, loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return save(fig, out_dir, "fig2_generality")


def figure_adjacency(h, out_dir):
    """The Letter's Figure 3 (2026-10-07, G2R): per decisive ion the joint
    transport error of the physical frequency ordering (diamond) against
    the 31 scrambled orderings (grey) at every block count -- same blocks,
    rank bound, parameter count, exit tables, events, seeds and transport;
    only the adjacency of the frequencies in a block differs. Single column."""
    g = h["g2r"]; names = h["ion_names"]; dm_max = h["g1"]["dm_max"]
    ions = [i for i in h["decisive"] if i in g["per_ion"]]
    fig, axes = plt.subplots(len(ions), 1, figsize=(3.4, 2.6 * len(ions)), squeeze=False, sharex=True)
    for ri, ion in enumerate(ions):
        ax = axes[ri][0]; r = g["per_ion"][ion]
        for c in r["cells"]:
            k = c["k"]
            jit = np.random.default_rng(k).uniform(-0.12, 0.12, len(c["scrambled"]))
            ax.plot(k * 2 ** jit, c["scrambled"], ".", color="0.6", ms=3.5, alpha=0.8)
            ax.plot([k], [c["e_physical"]], "D", color=COL[ion], ms=6, mec="k")
        ax.axhline(dm_max, color="k", ls=":", lw=1.1)
        ax.set_yscale("log"); ax.set_xscale("log", base=2)
        ax.set_xticks(g["k_grid"]); ax.set_xticklabels([str(k) for k in g["k_grid"]])
        ax.set_ylabel("max(band, colour) error  [mag]", fontsize=8)
        ax.set_title(f"{names[ion]}: physical order (diamond) vs {g['n_scrambled']} scrambled", fontsize=8, loc="left")
        ax.tick_params(labelsize=8)
    axes[-1][0].set_xlabel("contiguous blocks $N_g$ of the ordering", fontsize=8)
    fig.tight_layout()
    return save(fig, out_dir, "fig3_adjacency")


def main(out_dir=None, which=None):
    h = json.loads(FROZEN.read_text())
    out_dir = Path(out_dir) if out_dir else HERE / "figures"
    written = []
    if which in (None, 1):
        written += figure1(h, out_dir)
    if which in (None, 2):
        written += figure3(h, out_dir)            # the Letter's Figure 2: generality
    if which in (None, 3):
        written += figure_adjacency(h, out_dir)   # the Letter's Figure 3: adjacency (G2R)
    if which in (None, 5):
        written += figure_robustness(h, out_dir)  # End Matter
    if which == 4:
        written += figure2(h, out_dir)            # the former standalone mechanism figure, for the record
    return written


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--which", type=int, default=None)
    a = ap.parse_args()
    for p in main(a.out, a.which):
        print("wrote", p)
