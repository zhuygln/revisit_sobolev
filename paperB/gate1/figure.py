#!/usr/bin/env python3
"""The central figure of gate G1: transport error against effective-model
complexity for the three ions, with eps* as a line and the seed noise as a
band. Reads gate1_verdict.json; writes docs/figures/paperB/."""
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
COL = {"57LaII": "#0072B2", "58CeII": "#E69F00", "60NdII": "#D55E00"}
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}


def main(verdict_path=HERE / "gate1_verdict.json", out_dir=OUT):
    v = json.loads(Path(verdict_path).read_text())
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.2))
    for ion, r in v["per_ion"].items():
        x = [t["n_matrix_dof"] for t in r["table"]]; c = COL[ion]; lab = NAME[ion] + (" (gray)" if r["gray"] else "")
        axes[0].plot(x, [t["band_max"] for t in r["table"]], "o-", color=c, label=lab)
        axes[0].axhline(r["eps_star"]["max_dm"], color=c, ls="--", lw=0.9)
        noise = r["eps_star"]["noise_on_dm"]
        axes[0].axhspan(0, noise, color="grey", alpha=0.08)
        axes[1].plot(x, [t["sed"] for t in r["table"]], "o-", color=c, label=lab); axes[1].axhline(r["eps_star"]["sed"], color=c, ls="--", lw=0.9)
        axes[2].plot(x, [t["event"] for t in r["table"]], "o-", color=c, label=lab)
        axes[3].plot([t["ng"] for t in r["table"]], [t["table_kb"] for t in r["table"]], "o-", color=c, label=lab + " (matrix + exit tables)")
    axes[0].axhline(v["prereg"]["dm_max"], color="k", ls=":", lw=1, label="B1 threshold")
    axes[0].set_ylabel("max |Δm| over live bands vs macroatom [mag]"); axes[0].set_title("band error (dashed: ε*)", fontsize=9)
    axes[1].set_ylabel("SED L1 error"); axes[1].set_title("SED error (dashed: ε*)", fontsize=9)
    axes[2].set_ylabel("event-level total-variation loss vs the independent 128-group matrix"); axes[2].set_title("information discarded by the binning", fontsize=9)
    for ax in axes[:3]:
        ax.set_xscale("log"); ax.set_xlabel("coarse transition-matrix degrees of freedom $N_g^2$"); ax.legend(fontsize=7)
    axes[3].set_xscale("log"); axes[3].set_yscale("log"); axes[3].set_xlabel("$N_g$"); axes[3].set_ylabel("serialized kernel size [kB]"); axes[3].set_title("what the whole representation stores", fontsize=9); axes[3].legend(fontsize=7)
    fig.suptitle(f"Gate G1: transport error against the coarse matrix size — B1 {v['B1']}, B2 {v['B2']}, B3 {v['B3']} → {v['decision']}", fontsize=10)
    fig.tight_layout()
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}), ("png", {"Software": None})):
        p = out_dir / f"gate1_error_vs_complexity.{ext}"; fig.savefig(p, metadata=meta, dpi=150, bbox_inches="tight"); written.append(p)
    plt.close(fig)
    return written


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
