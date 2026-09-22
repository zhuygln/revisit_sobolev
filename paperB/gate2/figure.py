#!/usr/bin/env python3
"""The G2 figure: transport error against the archetype count for the local
and the global family, the event-level loss alongside, the truncation
curve, and every operator on the (m_event, max|dm|) plane. Reads
gate2_verdict.json; writes docs/figures/paperB/gate2_locality_vs_rank."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "docs/figures/paperB"
COL = {"57LaII": "#0072B2", "58CeII": "#E69F00", "60NdII": "#D55E00"}
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}


def main(verdict_path=HERE / "gate2_verdict.json", out_dir=OUT):
    v = json.loads(Path(verdict_path).read_text())
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.2))
    for ion, r in v["per_ion"].items():
        c = COL[ion]; lab = NAME[ion] + (" (gray)" if r["gray"] else "")
        loc = sorted((int(k), m) for k, m in r["local"].items()); glo = sorted((int(k), m) for k, m in r["global_nmf"].items())
        tru = sorted((float(k), m) for k, m in r["truncation"].items())
        axes[0].plot([k for k, _ in loc], [m["band_max"] for _, m in loc], "o-", color=c, label=f"{lab} local (N_g)")
        axes[0].plot([k for k, _ in glo], [m["band_max"] for _, m in glo], "s--", color=c, mfc="none", label=f"{lab} global (rank k)")
        axes[1].plot([k for k, _ in loc], [m["event"] for _, m in loc], "o-", color=c, label=f"{lab} local")
        axes[1].plot([k for k, _ in glo], [m["event"] for _, m in glo], "s--", color=c, mfc="none", label=f"{lab} global")
        axes[2].plot([m["rho_exit"] for _, m in tru], [m["band_max"] for _, m in tru], "^-", color=c, label=lab)
        for fam, mk in ((loc, "o"), (glo, "s"), (tru, "^")):
            axes[3].scatter([m["event"] for _, m in fam], [m["band_max"] for _, m in fam], marker=mk, color=c, facecolors="none" if mk == "s" else c, s=28)
    for ax in (axes[0], axes[2], axes[3]):
        ax.axhline(v["prereg"]["dm_max"], color="k", ls=":", lw=1)
    axes[0].set_xscale("log", base=2); axes[0].set_xlabel("archetypal exit distributions ($N_g$ local, $k$ global)"); axes[0].set_ylabel("max |Δm| over live bands vs macroatom [mag]")
    axes[0].set_title("H1: locality vs rank", fontsize=9); axes[0].legend(fontsize=6)
    axes[1].set_xscale("log", base=2); axes[1].set_xlabel("archetypes"); axes[1].set_ylabel("event-level TV loss vs the independent 128-group matrix"); axes[1].set_title("what each family discards", fontsize=9); axes[1].legend(fontsize=6)
    axes[2].set_xscale("log"); axes[2].set_xlabel("fraction of exit lines kept (per-group energy truncation)"); axes[2].set_ylabel("max |Δm| [mag]"); axes[2].set_title("H2: the exit lines the light needs", fontsize=9); axes[2].legend(fontsize=7)
    axes[3].set_xlabel("event-level TV loss"); axes[3].set_ylabel("max |Δm| [mag]"); axes[3].set_title("every operator: events discarded vs observable error\n(o local, □ global, △ truncation)", fontsize=9)
    fig.suptitle(f"Gate G2 — H1 {v['H1']}, H2 {v['H2']} → {v['decision']}", fontsize=10)
    fig.tight_layout()
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}), ("png", {"Software": None})):
        p = out_dir / f"gate2_locality_vs_rank.{ext}"; fig.savefig(p, metadata=meta, dpi=150, bbox_inches="tight"); written.append(p)
    plt.close(fig)
    return written


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
