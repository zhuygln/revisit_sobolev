"""Figure 11: is the separation universal? 36 conditions across
windows x epochs x ion mixes, plotted against realized tau_max."""

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
rows = [r for r in json.loads((HERE / "breadth_results_fixed.json").read_text())
        if r.get("d_exp") is not None]

WCOL = {4300: "C0", 4900: "C1", 7000: "C2", 9100: "C3"}
MMK = {"La": "o", "LaCe": "s", "LaCeCe3": "^"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True)
for r in rows:
    kw = dict(color=WCOL[int(r["window"])], marker=MMK[r["mix"]], ms=6,
              alpha=0.85, ls="none")
    t = max(r["tau_max"], 3e-3)
    ax1.plot(t, 100 * r["d_exp"], **kw)
    ax2.plot(t, 100 * r["d_sob"], **kw)

for ax, ttl in ((ax1, "expansion opacity"), (ax2, "Sobolev proper")):
    ax.set_xscale("log")
    ax.set_xlabel(r"realized $\tau_{\max}$")
    ax.set_ylabel(r"$\Delta$ [%]")
    ax.set_title(ttl)
    ax.axhline(0, color="k", lw=0.6)
    ax.grid(alpha=0.3)
ax1.set_ylim(-5, 68)
ax2.set_ylim(-5, 68)  # same scale: the point is the size difference

wh = [plt.Line2D([], [], color=c, marker="o", ls="none", label=f"{w} $\\AA$")
      for w, c in WCOL.items()]
mh = [plt.Line2D([], [], color="gray", marker=m, ls="none", label=k)
      for k, m in MMK.items()]
ax1.legend(handles=wh, fontsize=7, loc="upper left", title="window",
           title_fontsize=7)
ax2.legend(handles=mh, fontsize=7, loc="upper left", title="ion mix",
           title_fontsize=7)

fig.suptitle(
    "Breadth sweep: 4 windows $\\times$ 3 epochs $\\times$ 3 ion mixes "
    "(36 conditions, same y-scale)", y=1.0
)
fig.tight_layout()
out = HERE.parents[1] / "outputs" / "fig11_breadth.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("saved", out)

# quantitative summary
for name, key in (("expansion", "d_exp"), ("Sobolev", "d_sob")):
    v = np.array([r[key] for r in rows]) * 100
    t = np.array([r["tau_max"] for r in rows])
    print(f"{name:10s}: median {np.median(v):+5.1f}%  max {v.max():+5.1f}%  "
          f"| tau>3: median {np.median(v[t>3]):+5.1f}%  "
          f"| tau<0.5: median {np.median(v[t<0.5]):+5.1f}%")
ratio = np.array([r["d_exp"] / r["d_sob"] for r in rows
                  if r["d_sob"] > 0.01 and r["tau_max"] > 1])
print(f"expansion/Sobolev error ratio for tau_max>1: median {np.median(ratio):.1f}x")
