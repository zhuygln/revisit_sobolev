"""Figure 13: is the expansion-opacity error a bin-width artifact? No."""
import json
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = Path(__file__).parent; ROOT = HERE.parents[1]
V_D = 100.0
la = [r for r in json.loads((HERE/"binwidth_results.json").read_text()) if r.get("d_exp") is not None]
bl = [r for r in json.loads((HERE/"blend_binwidth.json").read_text()) if r.get("d_exp") is not None]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
ok  = [r for r in la if r["bin_kms"] <= 1.5*V_D]
bad = [r for r in la if r["bin_kms"] >  1.5*V_D]

# left: the reference's own validity
ax1.semilogx([r["bin_kms"] for r in ok],  [r["bb"] for r in ok],  "o-", color="C1",
             label="resolved reference (valid)")
ax1.semilogx([r["bin_kms"] for r in bad], [r["bb"] for r in bad], "x--", color="C3",
             label="reference unconverged")
ax1.axvline(V_D, color="gray", ls=":", lw=1.2)
ax1.text(V_D*1.15, 0.60, "bin = $v_D$", fontsize=8, color="gray")
ax1.set_xlabel("transport bin width [km/s]"); ax1.set_ylabel(r"$F_{\rm resolved}$")
ax1.set_title("The reference fails once bins exceed the profile")
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

# right: the actual test, against lines per bin
ax2.semilogx([r["lines_per_bin"] for r in ok], [100*r["d_exp"] for r in ok],
             "o-", color="C0", label="La II (153 lines)")
ax2.semilogx([r["lines_per_bin"] for r in bl], [100*r["d_exp"] for r in bl],
             "s-", color="C2", label="La II + Ce II (2529 lines)")
ax2.semilogx([r["lines_per_bin"] for r in bad], [100*r["d_exp"] for r in bad],
             "x", color="C3", label="excluded: reference unconverged")
ax2.axvline(1.0, color="gray", ls=":", lw=1.2)
ax2.text(1.15, -25, "1 line per bin", fontsize=8, color="gray", rotation=90)
ax2.set_xlabel("lines per transport bin"); ax2.set_ylabel(r"$\Delta_{\rm expansion}$ [%]")
ax2.set_title("The error is flat across the design regime")
ax2.legend(fontsize=7); ax2.grid(alpha=0.3); ax2.axhline(0, color="k", lw=0.6)

fig.tight_layout()
out = ROOT/"outputs"/"fig13_binwidth.png"; fig.savefig(out, dpi=200)
print("saved", out)
d1 = np.array([r["d_exp"] for r in ok]); d2 = np.array([r["d_exp"] for r in bl])
print(f"La II  : {100*d1.mean():+.1f}% +- {100*d1.std():.1f}, "
      f"{ok[0]['lines_per_bin']:.3f}-{ok[-1]['lines_per_bin']:.2f} lines/bin")
print(f"blend  : {100*d2.mean():+.1f}% +- {100*d2.std():.1f}, "
      f"{bl[0]['lines_per_bin']:.1f}-{bl[-1]['lines_per_bin']:.1f} lines/bin")
