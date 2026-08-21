"""Summarize the single-line ladder (ladder.py) into a table and a figure.

Targets: frozen law 0.1371 (what a steady-iterate code solves), classical
e^-tau_S = 0.1353 (the beta -> 0 limit). Each rung: seed mean, sem, Poisson
expectation, and the offset from the frozen target in percent. The anchor
(3.2e7 packets, 4.17e-5 grid, 5 seeds) is the converged value the paper
quotes; the zero-opacity control must be 1.000.
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN, CLASSICAL = 0.1371, 0.1353
res = json.loads((HERE / "ladder_results.json").read_text())

rows = []
print(f"{'rung':20s} {'n':>2s} {'mean':>8s} {'sem':>7s} {'poisson':>7s} {'vs frozen':>9s}")
for tag in sorted(res):
    vals = [v["trough"] for v in res[tag].values() if v["trough"] is not None]
    errs = [v["poisson_err"] for v in res[tag].values() if v["poisson_err"] is not None]
    if not vals:
        continue
    a = np.array(vals)
    row = dict(tag=tag, n=len(a), mean=float(a.mean()), sem=float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else None,
               poisson=float(np.mean(errs)) if errs else None,
               vs_frozen=float(a.mean() / FROZEN - 1.0) if not tag.startswith("Z_") else None)
    rows.append(row)
    sem_s = "   n/a " if row["sem"] is None else f"{row['sem']:7.5f}"
    poi_s = "   n/a " if row["poisson"] is None else f"{row['poisson']:7.4f}"
    vsf_s = "" if row["vs_frozen"] is None else f"{100*row['vs_frozen']:+8.2f}%"
    print(f"{tag:20s} {len(a):2d} {a.mean():8.5f} {sem_s:>7s} {poi_s:>7s} {vsf_s:>9s}")
(HERE / "ladder_summary.json").write_text(json.dumps(rows, indent=1))

# figure: one panel per axis
fig, axs = plt.subplots(1, 4, figsize=(13, 3.4), sharey=True)
axes = {"A": ("packets", lambda t: float(t.split("emit")[1]), "log"),
        "B": (r"transport d$\nu/\nu$", lambda t: float(t.split("dnut")[1]), "log"),
        "C": (r"spectrum d$\nu/\nu$", lambda t: float(t.split("dnus")[1]), "log"),
        "D": ("zones", lambda t: float(t.split("nz")[1]), "linear")}
for ax, (key, (xl, fx, sc)) in zip(axs, axes.items()):
    pts = [(fx(r["tag"]), r["mean"], r["sem"] or 0.0) for r in rows if r["tag"].startswith(key + "_")]
    if pts:
        x, y, e = map(np.array, zip(*sorted(pts)))
        ax.errorbar(x, y, yerr=e, fmt="o-", color="C1", capsize=3)
    ax.axhline(FROZEN, color="k", ls="--", lw=0.9, label="frozen law 0.1371")
    ax.axhline(CLASSICAL, color="r", ls=":", lw=0.9, label=r"$e^{-\tau_S}$ 0.1353")
    ax.set_xscale(sc); ax.set_xlabel(xl); ax.grid(alpha=.3)
anchor = [r for r in rows if r["tag"] == "E_anchor_bb"]
if anchor:
    for ax in axs:
        ax.axhspan(anchor[0]["mean"] - (anchor[0]["sem"] or 0), anchor[0]["mean"] + (anchor[0]["sem"] or 0), color="C1", alpha=0.15)
axs[0].set_ylabel("trough depth (1400-2600 km/s)"); axs[0].legend(fontsize=7)
fig.suptitle("single-line benchmark: SEDONA resolved, one axis at a time (fixed seeds)", fontsize=9)
fig.tight_layout()
for out in (ROOT / "outputs/fig_ladder.png", ROOT / "docs/figures/fig_ladder.png"):
    out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=170)
print("wrote ladder_summary.json, fig_ladder.png")
