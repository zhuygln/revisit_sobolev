"""Figure 7: first validity-map slice -- Delta_Sob(tau_max, v_D) for the
La II 3850-3950 A window at T = 3000 K, day 1."""

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
rows = json.loads((HERE / "sweep_results.json").read_text())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.2))
for tau_max, color in [(0.5, "C2"), (5.0, "C1"), (50.0, "C3")]:
    sel = sorted(
        (r for r in rows if r["tau_max"] == tau_max), key=lambda r: r["v_d_kms"]
    )
    v = [r["v_d_kms"] for r in sel]
    ax1.semilogx(v, [100 * r["delta_sob"] for r in sel], "o-", color=color,
                 label=rf"$\tau_{{\max}}$ = {tau_max:g}")
    ax2.semilogx(v, [r["bb"] for r in sel], "o-", color=color)
    ax2.semilogx(v, [r["exp"] for r in sel], "o--", color=color, alpha=0.55)

ax1.set_xlabel(r"$v_D$ [km/s]")
ax1.set_ylabel(r"$\Delta_{\rm Sob}$ = (F$_{\rm exp}$ − F$_{\rm bb}$)/F$_{\rm bb}$  [%]")
ax1.set_title("Expansion-opacity band-flux error")
ax1.legend()
ax1.grid(alpha=0.3)

ax2.set_xlabel(r"$v_D$ [km/s]")
ax2.set_ylabel(r"band-averaged $F/F_{\rm cont}$")
ax2.set_title("Resolved (solid) vs expansion (dashed)")
ax2.grid(alpha=0.3)

fig.suptitle(
    "La II 3850-3950 $\\AA$ forest, T = 3000 K, day 1 -- validity-map slice",
    y=1.0,
)
fig.tight_layout()
out = HERE.parents[1] / "outputs" / "fig7_validity_slice.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("saved", out)
