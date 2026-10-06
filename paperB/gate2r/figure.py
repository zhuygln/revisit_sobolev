#!/usr/bin/env python3
"""G2R figure: per decisive ion, the joint transport error of the physical
ordering against the 31 scrambled orderings at every block count (top),
and the event-level loss of the same orderings (bottom).

    .venv/bin/python paperB/gate2r/figure.py      # docs/figures/paperB/g2r_adjacency.{pdf,png}
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "docs/figures/paperB"
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}
COL = {"57LaII": "#0072B2", "58CeII": "#E69F00", "60NdII": "#D55E00"}


def draw(out, axes=None):
    ions = [i for i in out["decisive"] if i in out["per_ion"]]
    fig = None
    if axes is None:
        fig, axes = plt.subplots(2, len(ions), figsize=(4.2 * len(ions), 6.4), squeeze=False)
    dm_max = out["prereg"]["dm_max"]
    for ci, ion in enumerate(ions):
        r = out["per_ion"][ion]
        for ri, key in enumerate(("e_joint", "event")):
            ax = axes[ri][ci]
            for k, c in sorted(r["cells"].items(), key=lambda kv: int(kv[0])):
                k = int(k)
                scr = [o[key] for o in c["orderings"] if o["perm_id"] != 0]
                phys = next(o[key] for o in c["orderings"] if o["perm_id"] == 0)
                jit = np.random.default_rng(k).uniform(-0.12, 0.12, len(scr))
                ax.plot(k * 2 ** jit, scr, ".", color="0.6", ms=4, alpha=0.8)
                ax.plot([k], [phys], "D", color=COL[ion], ms=7, mec="k")
            if key == "e_joint":
                ax.axhline(dm_max, color="k", ls=":", lw=1.1)
                ax.set_yscale("log"); ax.set_ylabel("max(band, colour) error vs $R_2$  [mag]")
                ax.set_title(f"{NAME[ion]}: physical order (diamond) vs 31 scrambled", fontsize=9, loc="left")
            else:
                ax.set_ylabel("event-level total-variation loss")
            ax.set_xscale("log", base=2); ax.set_xticks([2, 4, 8, 16, 32]); ax.set_xticklabels(["2", "4", "8", "16", "32"])
            ax.set_xlabel("contiguous blocks $N_g$ (the ordering's)")
    return fig


def main(out_dir=OUT):
    out = json.loads((HERE / "g2r_verdict.json").read_text())
    fig = draw(out); fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}), ("png", {"Software": None})):
        p = out_dir / f"g2r_adjacency.{ext}"; fig.savefig(p, metadata=meta, dpi=200, bbox_inches="tight"); paths.append(p)
    plt.close(fig)
    return paths


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
