"""Figure for E4/E5: F_b(eps) on both TLA legs against the direct-branching
targets, per band, and chi^2/dof(eps). Reads e4_eps_sweep.json."""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
d = json.load(open(HERE / "e4_eps_sweep.json")); EPS = np.array(d["eps"])
bands = [b for b in d["bands"] if b != "band3800"] + ["band3800"]
fig, axs = plt.subplots(2, 3, figsize=(12, 7)); axs = axs.ravel()
for ax, b in zip(axs, bands):
    for leg, c in (("sobolev_tla", "C0"), ("expansion_tla", "C3")):
        F = [np.mean(d["legs"][f"{leg}_eps{e:g}"]["bands"][b]) for e in EPS]
        Fs = [np.std(d["legs"][f"{leg}_eps{e:g}"]["bands"][b], ddof=1) for e in EPS]
        ax.errorbar(EPS, F, Fs, color=c, marker="o", ms=3, label=leg.replace("_", "+"))
    for ref, c, ls in (("sobolev_branch", "C0", "--"), ("expansion_branch", "C3", ":")):
        m = np.mean(d["legs"][ref]["bands"][b]); s = np.std(d["legs"][ref]["bands"][b], ddof=1)
        ax.axhline(m, color=c, ls=ls, label=ref.replace("_", "+")); ax.axhspan(m - s, m + s, color=c, alpha=0.12)
    lo, hi = d["bands"][b]; ax.set_title(f"{b}: {lo:.0f}-{hi:.0f} A", fontsize=10); ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel("escaped / launched energy")
axs[0].legend(fontsize=7)
fig.suptitle("La II, 6000 K Planck launch: two-level-atom closure vs direct A-branching fluorescence", fontsize=11)
fig.tight_layout()
for out in (ROOT / "outputs/fig_p2_eps_sweep.png", ROOT / "docs/figures/fig_p2_eps_sweep.png"):
    out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=150)
fig2, ax = plt.subplots(figsize=(5, 3.5))
for leg, c in (("sobolev_tla", "C0"), ("expansion_tla", "C3")):
    chi = [d["chi2"][leg][f"{e:g}"]["total"] / d["chi2"][leg][f"{e:g}"]["dof"] for e in EPS]
    ax.plot(EPS, chi, color=c, marker="o", label=leg.replace("_", "+"))
ax.set_yscale("log"); ax.set_xlabel(r"$\epsilon$"); ax.set_ylabel(r"$\chi^2$/dof vs Sobolev+branch (200 bins)"); ax.legend(fontsize=8)
fig2.tight_layout()
for out in (ROOT / "outputs/fig_p2_eps_chi2.png", ROOT / "docs/figures/fig_p2_eps_chi2.png"):
    fig2.savefig(out, dpi=150)
print("wrote fig_p2_eps_sweep.png, fig_p2_eps_chi2.png")
