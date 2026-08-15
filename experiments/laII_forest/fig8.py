"""Figure 8: temperature axis + thermal-width frontier of the validity map."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
rows = json.loads((HERE / "tsweep_results.json").read_text())
sweep = json.loads((HERE / "sweep_results.json").read_text())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.2))

# --- left: Delta_Sob vs T (tau_max pinned at 5, v_D = 100 km/s) ---
trows = sorted((r for r in rows if r["axis"] == "T"), key=lambda r: r["T"])
T = [r["T"] for r in trows]
ax1.plot(T, [100 * r["delta_sob"] for r in trows], "o-", color="C1")
for r in trows:
    ax1.annotate(
        f"{r['n_tau_gt1']}+{r['n_tau_01_1']}",
        (r["T"], 100 * r["delta_sob"]),
        textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8,
    )
ax1.set_xlabel("shell temperature T [K]")
ax1.set_ylabel(r"$\Delta_{\rm Sob}$ [%]")
ax1.set_title(
    r"T axis ($\tau_{\max}$ pinned at 5); labels: lines with $\tau>1$ + $\tau$ 0.1-1"
)
ax1.grid(alpha=0.3)

# --- right: v_D axis at tau_max = 5, T = 3000 K, extended to the frontier ---
v_old = sorted(
    (r for r in sweep if r["tau_max"] == 5.0), key=lambda r: r["v_d_kms"]
)
v_new = sorted((r for r in rows if r["axis"] == "vd"), key=lambda r: r["v_d_kms"])
pts = [(r["v_d_kms"], 100 * r["delta_sob"]) for r in v_new] + [
    (r["v_d_kms"], 100 * r["delta_sob"]) for r in v_old
]
pts.sort()
ax2.semilogx([p[0] for p in pts], [p[1] for p in pts], "o-", color="C1")
ax2.axvline(0.6, color="gray", ls="--", lw=1)
ax2.text(0.65, ax2.get_ylim()[0] + 2, "La thermal\n0.6 km/s", fontsize=8,
         color="gray")
for r in v_new:
    ax2.annotate(
        f"{r['bb_wall_s']:.0f}s", (r["v_d_kms"], 100 * r["delta_sob"]),
        textcoords="offset points", xytext=(0, -14), ha="center", fontsize=8,
        color="C3",
    )
ax2.set_xlabel(r"$v_D$ [km/s]")
ax2.set_ylabel(r"$\Delta_{\rm Sob}$ [%]")
ax2.set_title(
    r"$v_D$ axis at $\tau_{\max}$ = 5 (red: resolved-run wall time)"
)
ax2.grid(alpha=0.3)

fig.suptitle(
    "La II 3850-3950 $\\AA$ forest, day 1 -- validity-map rows 2-3", y=1.0
)
fig.tight_layout()
out = HERE.parents[1] / "outputs" / "fig8_T_and_frontier.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("saved", out)
