"""Paper II Phase 1 figure: emergent spectra under the five treatments.

Left: Paper I's window atom, flat launch, with SEDONA's RE N=1 spectra
overlaid (the calibration). Right: the full La II atom with a Planck 6000 K
launch, 3700-4100 A zoom of the emergent/incident ratio (the physics).
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.formal_transfer import planck_bnu

z = np.load(HERE / "forest_spectra_s1.npz")
lamA = C / np.sqrt(z["edges_A"][1:] * z["edges_A"][:-1]) * 1e8
lamB = C / np.sqrt(z["edges_B"][1:] * z["edges_B"][:-1]) * 1e8


def sed(path, scale=None):
    s = np.loadtxt(path, comments="#"); nu, l = s[:, 0], s[:, 1]; g = l > 0; nu, l = nu[g], l[g]
    lam = C / nu * 1e8; r = l / (4 * np.pi**2 * 8.64e12**2 * planck_bnu(nu, 6000.0))
    o = np.argsort(lam); return lam[o], r[o]


lam0, r0 = sed(ROOT / "experiments/laII_forest/run_bb/spectrum_1.dat")
scale = r0[(lam0 > 3952) & (lam0 < 3970)].mean()
sed_bb = sed(ROOT / "experiments/laII_forest/re_prod_N1_bb_s1/spectrum_1.dat"); sed_ex = sed(ROOT / "experiments/laII_forest/re_prod_N1_exp_s1/spectrum_1.dat")


def smooth(y, k=5):
    y = np.nan_to_num(y); return np.convolve(y, np.ones(k) / k, mode="same")


fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axs[0]
ax.plot(sed_bb[0], sed_bb[1] / scale, color="0.6", lw=0.8, label="SEDONA RE N=1, resolved")
ax.plot(sed_ex[0], sed_ex[1] / scale, color="0.8", lw=0.8, label="SEDONA RE N=1, expansion")
for tag, c, lab in (("A_sobolev_absorb", "C2", "Sobolev, absorb"), ("A_expansion_absorb", "C0", "expansion, absorb"),
                    ("A_sobolev_thermal_window", "C2", "Sobolev + thermal (window)"), ("A_expansion_thermal_window", "C0", "expansion + thermal (window)"),
                    ("A_sobolev_branch", "C3", "Sobolev + fluorescence")):
    ls = "--" if "absorb" in tag else "-"
    ax.plot(lamA, smooth(z[tag]), color=c, ls=ls, lw=1.2, label=lab)
ax.set_xlim(3770, 3997); ax.set_ylim(0, 1.9); ax.set_xlabel("wavelength [A]"); ax.set_ylabel("emergent / incident")
ax.set_title("Paper I's atom (3850-3950 A lines), flat launch: calibration vs SEDONA", fontsize=9); ax.legend(fontsize=7, ncol=2); ax.grid(alpha=.3)
ax = axs[1]
for tag, c, lab in (("B_sobolev_absorb", "C2", "Sobolev, absorb"), ("B_expansion_absorb", "C0", "expansion, absorb"),
                    ("B_sobolev_thermal", "C2", "Sobolev + thermal (whole atom)"), ("B_expansion_thermal", "C0", "expansion + thermal"),
                    ("B_sobolev_branch", "C3", "Sobolev + fluorescence")):
    ls = "--" if "absorb" in tag else "-"
    ax.plot(lamB, smooth(z[tag], 3), color=c, ls=ls, lw=1.2, label=lab)
ax.axvspan(3800, 3955, color="k", alpha=0.06, label="band 3800-3955")
ax.set_xlim(3600, 4200); ax.set_ylim(0, 1.5); ax.set_xlabel("wavelength [A]"); ax.set_title("full La II atom, Planck 6000 K launch: the physics", fontsize=9)
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=.3)
fig.tight_layout()
for out in (ROOT / "outputs/fig_p2_phase1.png", ROOT / "docs/figures/fig_p2_phase1.png"):
    out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=170)
print("wrote fig_p2_phase1.png")
