#!/usr/bin/env python3
"""The G3 figure: per ion, the anchor's transfer error, the fresh R_16's
error and (at the interior) the whole-operator interpolation's error along
each axis against the frozen coordinate, with the A/B/C/D reading as the
marker colour and the 0.10 mag criterion drawn. Reads gate3_verdict.json;
writes docs/figures/paperB/gate3_transfer_vs_existence."""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "docs/figures/paperB"
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}
CLS = {"A": "#009E73", "B": "#0072B2", "C": "#E69F00", "D": "#D55E00", "GRAY": "0.6", "REF": "k"}
AXLAB = {"T": r"$\log T_{\rm gas}$", "D": r"$\log n_{\rm ion}$", "J": r"$\log T_{\rm core}$", "P": r"$\log t$"}


def main(verdict_path=HERE / "gate3_verdict.json", out_dir=OUT):
    v = json.loads(Path(verdict_path).read_text())
    ions = [i for i in ("58CeII", "60NdII", "57LaII") if i in v["per_ion"]]
    axes_names = list(v["prereg"]["axes"])
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(ions), len(axes_names), figsize=(4.2 * len(axes_names), 3.3 * len(ions)), squeeze=False)
    ref = {s["ion"]: s for s in v["states"] if s["axis"] == "ref"}
    for ri, ion in enumerate(ions):
        for ci, ax_name in enumerate(axes_names):
            ax = axes[ri][ci]
            sts = sorted([s for s in v["states"] if s["ion"] == ion and s["axis"] == ax_name], key=lambda s: s["coord_value"])
            # the reference sits on every axis at its own coordinate
            r = ref.get(ion)
            if r:
                a = v["prereg"]["axes"][ax_name]
                cref = dict(T=np.log(3401.0), J=np.log(3401.0), P=np.log(2 * 86400.0)).get(ax_name)
                if cref is None:                                  # D: the reference density from the neighbours' coordinates
                    cref = next((s["coord_value"] - np.log(s["value"]) for s in sts), None)
                if cref is not None:
                    ax.plot([cref], [r["recomputed"].get("16", r["recomputed"].get(16, {})).get("band", np.nan)], "k*", ms=9, label="θ₀ (anchor)")
            for s in sts:
                c = s["coord_value"]; col = CLS[s["classification"]]
                r16 = s["recomputed"].get("16", s["recomputed"].get(16, {})).get("band", np.nan)
                ax.plot([c], [r16], "o", color=col, ms=7, mfc="none", label="fresh $R_{16}$" if s is sts[0] else None)
                if s["transfer"]:
                    ax.plot([c], [s["transfer"]["band"]], "s", color=col, ms=6, label="anchor $R_{16}$" if s is sts[0] else None)
                if s["interpolation"]:
                    ax.plot([c], [s["interpolation"]["whole"]["band"]], "D", color=col, ms=6, label="whole-operator interp.")
                    ax.plot([c], [s["interpolation"]["matrix_only"]["band"]], "d", color=col, ms=5, mfc="none", label="matrix-only interp.")
                ax.annotate(s["classification"], (c, max(r16, s["transfer"]["band"] if s["transfer"] else r16)), textcoords="offset points",
                            xytext=(0, 6), ha="center", fontsize=8, color=col)
            ax.axhline(v["prereg"]["dm_max"], color="k", ls=":", lw=1)
            ax.set_yscale("log"); ax.set_xlabel(AXLAB[ax_name])
            if ci == 0:
                ax.set_ylabel(f"{NAME[ion]}\nmax |Δm| vs R2(θ) [mag]")
            if ri == 0:
                ax.set_title(f"axis {ax_name}", fontsize=9)
            if ri == 0 and ci == 0:
                ax.legend(fontsize=6, loc="lower left")
    fig.suptitle(f"Gate G3 — C1 {v['C1']}, C2 {v['C2']}, C3 {v['C3']}; A transfers, B tabulable, C fresh operator passes, D fails", fontsize=10)
    fig.tight_layout()
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}), ("png", {"Software": None})):
        p = out_dir / f"gate3_transfer_vs_existence.{ext}"; fig.savefig(p, metadata=meta, dpi=150, bbox_inches="tight"); written.append(p)
    plt.close(fig)
    return written


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
