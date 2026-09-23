#!/usr/bin/env python3
"""The Paper B display items, generated from paperB/FROZEN.json alone (the
same frozen record the manuscript's macros come from; no number is typed
here). The PI froze these two on 2026-09-22.

  Figure 1 -- Compression survives the physics.  For La II, Ce II and Nd II:
      the optimally tuned scalar eps* and the group operators R_2 ... R_32
      against the energy-conserving downward macroatom, with the R2M
      robustness check (the same operator rebuilt from the
      radiation-field-driven macroatom's events) on the right.
      The message: eps* fails where a small R succeeds, and R follows the
      reference physics when that physics changes.

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


def figure1(h, out_dir):
    """eps* fails where a small R succeeds; and R follows the physics."""
    names, dm_max = h["ion_names"], h["g1"]["dm_max"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.1), gridspec_kw=dict(width_ratios=[1.45, 1]))
    ax = axes[0]
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
    ax.set_title("a. the scalar fails where a small operator succeeds", fontsize=9, loc="left")
    ax.legend(fontsize=7, loc="lower left")

    ax = axes[1]
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
    ax.set_title("b. when the fluorescence physics changes, the operator follows it", fontsize=9, loc="left")
    ax.legend(fontsize=7, loc="lower left", framealpha=0.95)
    ax.set_ylim(top=max(ref_shift) * 6)
    fig.tight_layout()
    return save(fig, out_dir, "fig1_compression")


def figure2(h, out_dir):
    """Locality, not rank; and events are not observables."""
    names, dm_max = h["ion_names"], h["g1"]["dm_max"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.1))
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
    axes[0].set_title("a. transport error (solid: local, open: global)", fontsize=9, loc="left")
    axes[0].legend(fontsize=6, ncol=2)
    axes[1].set_ylabel("event-level total-variation loss")
    axes[1].set_title("b. microscopic error: the global family wins", fontsize=9, loc="left")
    axes[2].set_yscale("log"); axes[2].set_xlabel("event-level total-variation loss")
    axes[2].set_ylabel(r"max $|\Delta m|$  [mag]")
    axes[2].set_title("c. better on events, worse on the light", fontsize=9, loc="left")
    axes[2].annotate("", xy=(0.08, 0.88), xytext=(0.40, 0.88), xycoords="axes fraction", textcoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", color="0.4", lw=1.1))
    axes[2].annotate("the global family fits the events better", (0.10, 0.91), xycoords="axes fraction", fontsize=7, color="0.35")
    axes[2].annotate("passes", (0.02, 0.10), xycoords="axes fraction", fontsize=7, color="0.35")
    fig.tight_layout()
    return save(fig, out_dir, "fig2_mechanism")


def main(out_dir=None, which=None):
    h = json.loads(FROZEN.read_text())
    out_dir = Path(out_dir) if out_dir else HERE / "figures"
    written = []
    if which in (None, 1):
        written += figure1(h, out_dir)
    if which in (None, 2):
        written += figure2(h, out_dir)
    return written


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--which", type=int, default=None)
    a = ap.parse_args()
    for p in main(a.out, a.which):
        print("wrote", p)
