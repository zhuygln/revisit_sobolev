"""The controlling statistic (referee Comment 8).

ln(1 + Delta_exp) against three abscissae: the exact predictor
ln <e^D>_w (w ~ e^-S), the transmission-weighted mean deficit <D>_w (a Jensen
lower bound), and tau_max (the empirical proxy). Points: the 12-point
tau_max x v_D grid (always) and the 36 breadth conditions (when
breadth_results_v2.json exists). The analytic pair sits on the 1:1 line by
identity; the SEDONA pair sits above it by its expansion-leg bin systematic.
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.rays import RaySet
from sobolev.sobolev_leg import crossing_depths
from sobolev.optical_depth import stimulated_emission_factor

pts = []  # dict(label, tau_max, x_exact, x_dw, d_exp_sedona, d_exp_ana, d_sob_det, marker)

# ---- 12-grid + frontier: recompute predictors from the forest
FOREST = ROOT / "experiments/laII_forest"
T_EXP, R_CORE, R_OUT = 86400.0, 8.64e12, 2.592e13
BAND = (3800.0, 3955.0)
d = np.load(FOREST / "forest_lines.npz")
N_ION = float(d["n_ion"])
LINES = [(C/(l*1e-8), f, p*stimulated_emission_factor(C/(l*1e-8), 3000.0)) for l, f, p in zip(d["lam"], d["f_lu"], d["pop"])]
nu = np.geomspace(7.50e14, 7.95e14, 4000); lam = C/nu*1e8; m = (lam > BAND[0]) & (lam < BAND[1])
rays = RaySet.midpoint(R_CORE, R_OUT, 400)
sep = json.loads((ROOT / "experiments/sobolev_proper/separation_results.json").read_text())
for r in sep["grid"] + sep["frontier"]:
    S, E, p, w = crossing_depths(nu, LINES, R_CORE, R_OUT, T_EXP, N_ION * r["tau_max"]/5.0 if "tau_max" in r else N_ION, rays=rays)
    wt = (w[:, None]*np.exp(-S))[:, m]
    x_exact = float(np.log(np.sum(wt*np.exp((S-E)[:, m]))/wt.sum()))
    x_dw = float(np.sum(wt*(S-E)[:, m])/wt.sum())
    # y = ln(F_exp / F_Sob): the closure's departure from the Sobolev leg, which
    # is what D controls (the identity is for F_exp/F_Sob, not F_exp/F_res --
    # dividing by the resolved reference would add the Sobolev localization
    # error, +33% at tau_max=50, v_D=300).
    pts.append(dict(label="grid", tau_max=r.get("tau_max", 5.0), v_d=r["v_d"], x_exact=x_exact, x_dw=x_dw,
                    y_sed=np.log(r["f_exp"] / r["f_sob"]),
                    y_ana=np.log(r.get("f_exp_ana", r["f_res_det"] * (1 + r["d_exp_det"])) / r["f_sob"]), marker="o"))

# ---- breadth, if available
bp = ROOT / "experiments/breadth/breadth_results_v2.json"
if bp.exists():
    for r in json.loads(bp.read_text()):
        if not r["stim"] or r.get("d_exp") is None: continue
        pts.append(dict(label="breadth", tau_max=r["tau_max"], v_d=100.0, x_exact=r["ln_ew_expD"], x_dw=r["d_w"],
                        y_sed=np.log(r["f_exp"] / r["f_sob"]), y_ana=np.log(r["f_exp_ana"] / r["f_sob"]),
                        marker={"La": "s", "LaCe": "^", "LaCeCe3": "D"}[r["mix"]]))
else:
    print("breadth_results_v2.json not yet available -- grid points only")

fig, axs = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
tm = np.array([q["tau_max"] for q in pts]); y = np.array([q["y_sed"] for q in pts]); ya = np.array([q["y_ana"] for q in pts])
cmap = plt.cm.viridis; norm = plt.matplotlib.colors.LogNorm(vmin=max(tm.min(), 1e-2), vmax=tm.max())
for ax, key, xl in ((axs[0], "x_exact", r"$\ln\langle e^{D}\rangle_{\tilde w}$ (exact)"),
                    (axs[1], "x_dw", r"$\langle D\rangle_{\tilde w}$ (Jensen bound)"),
                    (axs[2], "tau_max", r"$\tau_{\max}$ (proxy)")):
    x = np.array([q[key] for q in pts])
    for q, xi, yi, yai in zip(pts, x, y, ya):
        ax.scatter(xi, yi, marker=q["marker"], c=[cmap(norm(max(q["tau_max"], 1e-2)))], edgecolor="k", lw=0.4, s=36, zorder=3)
        if np.isfinite(yai): ax.scatter(xi, yai, marker=q["marker"], facecolor="none", edgecolor="gray", s=36, zorder=2)
    if key != "tau_max":
        lim = [0, max(x.max(), y.max())*1.05]; ax.plot(lim, lim, "k--", lw=0.8); ax.set_xlim(lim)
    else:
        ax.set_xscale("log")
    ax.set_xlabel(xl); ax.grid(alpha=.3)
axs[0].set_ylabel(r"$\ln(F_{\rm exp}/F_{\rm Sob})$")
axs[0].set_title("filled: SEDONA pair   open: analytic closure", fontsize=8)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
fig.colorbar(sm, ax=axs, label=r"$\tau_{\max}$", pad=0.01)
for out in (ROOT/"outputs/fig_predictor.png", ROOT/"docs/figures/fig_predictor.png"):
    out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=170, bbox_inches="tight")
print(f"{len(pts)} points; wrote fig_predictor.png")
(HERE/"predictor_points.json").write_text(json.dumps(pts, indent=1))
