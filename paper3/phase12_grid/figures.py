"""Figures 2-4 of the action plan, from committed JSON only.

  fig2_bol_vs_colour   central model: dm_bol(t) (conserving core, and the
                       absorbing-core bookkeeping) against d(g-r), d(r-i), d(i-J)
                       for C_both and C_binned, plus the spectra at the epoch
                       where |d(g-r)| is smallest
  fig3_velocity        from phase11's velocity_*.json: worst |dm| and |dm_bol|
                       against the ejecta velocity
  fig4_vectors         sensitivity.json: cosine matrix over the grid, and
                       d_RT against sum_theta a_theta d_theta for the central point

Usage: python figures.py [--central 0.01,0.1,0.01] [--out ../figures]
"""
import argparse, json, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
P11 = HERE.parents[0] / "phase11_observables"
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1]))
from grid import model_name          # noqa: E402
from sobolev import photometry as phot   # noqa: E402

COLS = ("g-r", "r-i", "i-J")
STYLE = {"C_both": dict(color="C3", marker="o"), "C_binned": dict(color="C0", marker="s"),
         "B_opacity": dict(color="C1", marker="^"), "A_redist": dict(color="0.5", marker="x")}


def save(fig, out, name):
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out / name}.png/.pdf")


def fig2(central, out):
    d = json.loads(central.read_text())
    rows = [r for r in d["rows"] if r.get("status") in ("ok", "reduced_n")]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4), gridspec_kw={"width_ratios": [1.1, 1.1, 1.4]})
    t = np.array([r["t_d"] for r in rows])
    for leg in ("C_both", "C_binned"):
        st = STYLE[leg]
        ax[0].plot(t, [r["legs"][leg]["dm_bol_absorbing"] for r in rows], "--", mfc="none",
                   label=f"{leg}: absorbing core", **st)
        ax[0].plot(t, [r["legs"][leg]["dm_bol"] for r in rows], "-", label=f"{leg}: conserving core", **st)
        for c, ls in zip(COLS, ("--", ":", "-.")):
            ax[1].plot(t, [r["legs"][leg]["dcolor"][c] for r in rows], ls, color=st["color"],
                       marker=st["marker"], label=f"{leg}: Δ({c})")
    for a in ax[:2]:
        a.axhline(0, color="k", lw=0.5); a.set_xlabel("t (d)"); a.set_xscale("log")
        a.set_xticks(t); a.set_xticklabels([f"{x:g}" for x in t]); a.minorticks_off()
    ax[0].set_ylabel("Δm_bol (mag)"); ax[0].legend(fontsize=7)
    ax[0].set_title("bolometric error (window 1000–30000 Å)\nabsorbing = inner-boundary bookkeeping", fontsize=9)
    ax[1].set_ylabel("Δcolour (mag)"); ax[1].set_title("colour error, DECam/2MASS"); ax[1].legend(fontsize=7, ncol=2)
    # spectra at the epoch where C_both's optical colour error |d(g-r)| is smallest
    k = int(np.argmin([abs(r["legs"]["C_both"]["dcolor"]["g-r"]) for r in rows]))
    r = rows[k]
    edges = phot.nu_edges(*d["lam_window"], d["n_spec"]); lam = phot.C / np.sqrt(edges[1:] * edges[:-1]) * 1e8
    for tag, lnu, st in (("reference", r["ref"]["L_nu"], dict(color="k")),
                         ("C_both", r["legs"]["C_both"]["L_nu"], STYLE["C_both"]),
                         ("C_binned", r["legs"]["C_binned"]["L_nu"], STYLE["C_binned"])):
        ax[2].plot(lam, np.asarray(lnu) * phot.C / (lam * 1e-8) ** 2 * 1e-8, lw=1, color=st["color"], label=tag)
    ax[2].set_xscale("log"); ax[2].set_yscale("log"); ax[2].set_xlabel("λ (Å)"); ax[2].set_ylabel("L_λ (erg/s/Å)")
    ax[2].set_xlim(3000, 25000); ax[2].set_xticks([3000, 5000, 10000, 20000]); ax[2].set_xticklabels(["3000", "5000", "10000", "20000"])
    ax[2].minorticks_off()
    lo = max(np.asarray(r["ref"]["L_nu"]).max() * phot.C / (3000e-8) ** 2 * 1e-8 * 1e-4, 1.0); ax[2].set_ylim(bottom=lo)
    cb = r["legs"]["C_both"]
    ax[2].set_title(f"t = {r['t_d']:g} d, C_both: Δm_bol = {cb['dm_bol']:+.2f}, Δ(g−r) = {cb['dcolor']['g-r']:+.2f}, "
                    f"Δ(i−J) = {cb['dcolor']['i-J']:+.2f}, Δ(J−K) = {cb['dcolor']['J-K']:+.2f}", fontsize=9)
    ax[2].legend(fontsize=8)
    fig.suptitle(f"M = {d['m_ej_msun']} M☉, v = {d['v_ej_c']} c, X_lan = {d['x_lan']}  (κ_src = {d['kappa_src']})", fontsize=10, y=1.04)
    save(fig, out, "fig2_bol_vs_colour")


def fig3(out):
    files = sorted(P11.glob("velocity_*.json"))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for f in files:
        d = json.loads(f.read_text())
        lab = f"{d['ion']} t={d['t_d']:g} d"
        rows = [r for r in d["rows"] if "legs" in r]
        b = [r["v_out"] for r in rows]
        for leg in ("C_both", "C_binned"):
            st = STYLE[leg]
            ax[0].plot(b, [max(abs(x) for x in r["legs"][leg]["dm"].values() if np.isfinite(x)) for r in rows],
                       "-" if leg == "C_both" else "--", marker=st["marker"], ms=4, label=f"{lab} {leg}")
            ax[1].plot(b, [abs(r["legs"][leg]["dm_bol"]) for r in rows],
                       "-" if leg == "C_both" else "--", marker=st["marker"], ms=4)
    for a, yl in zip(ax, ("worst |Δm| (mag)", "|Δm_bol| (mag)")):
        a.set_xscale("log"); a.set_xlabel("v_out / c"); a.set_ylabel(yl); a.grid(alpha=0.3)
    ax[0].legend(fontsize=6, ncol=2)
    fig.suptitle("closure error against ejecta velocity (phase 11 velocity scans, top-hat bands)", fontsize=10)
    save(fig, out, "fig3_velocity")


def fig4(sens, central_key, out):
    d = json.loads(sens.read_text())
    pts = list(d["points"])
    legs = ("C_both", "C_binned")
    fig, ax = plt.subplots(1, 3, figsize=(16, 6.5), gridspec_kw={"width_ratios": [1.3, 1.5, 1.4], "wspace": 0.12})
    for j, leg in enumerate(legs):
        M = np.full((len(pts), 3), np.nan)
        for i, p in enumerate(pts):
            r = d["points"][p]["legs"][leg]
            if r.get("status") == "ok":
                M[i] = r["cos"]
        im = ax[j].imshow(M, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
        ax[j].set_xticks(range(3)); ax[j].set_xticklabels(["cos(d_RT, d_M)", "cos(d_RT, d_v)", "cos(d_RT, d_X)"], fontsize=8)
        ax[j].set_yticks(range(len(pts)))
        ax[j].set_yticklabels(pts if j == 0 else [], fontsize=6)
        ax[j].set_xlim(-0.5, 3.6)
        for i, p in enumerate(pts):
            r = d["points"][p]["legs"][leg]
            txt = (r["cls"] + ("*" if r.get("unstable") else "")) if r.get("status") == "ok" else "—"
            ax[j].text(2.6, i, txt, fontsize=6, va="center")
        ax[j].set_title(f"{leg}: cosines and class (* = unstable)", fontsize=9)
    fig.colorbar(im, ax=ax[1], shrink=0.7, pad=0.02)
    r = d["points"].get(central_key, {}).get("legs", {}).get("C_both", {})
    if r.get("status") == "ok":
        # d_RT against the fitted combination sum_theta a_theta d_theta, one point per (band, epoch)
        d_rt, fit = np.array(r["d_rt"]), np.array(r["fit"])
        bands = [k[0] for k in r["keys"]]
        cmap = {b: f"C{i}" for i, b in enumerate(("g", "r", "i", "z", "J", "H", "K"))}
        for b in dict.fromkeys(bands):
            m = np.array([bb == b for bb in bands])
            ax[2].scatter(fit[m], d_rt[m], s=18, color=cmap[b], label=b)
        lim = 1.05 * max(np.max(np.abs(d_rt)), np.max(np.abs(fit)), 1e-3)
        ax[2].plot([-lim, lim], [-lim, lim], "k-", lw=0.5)
        ax[2].axhline(0, color="k", lw=0.3); ax[2].axvline(0, color="k", lw=0.3)
        ax[2].set_xlim(-lim, lim); ax[2].set_ylim(-lim, lim)
        ax[2].set_xlabel("Σ a_θ ∂m/∂lnθ (mag)"); ax[2].set_ylabel("d_RT = m_leg − m_ref (mag)")
        ax[2].legend(fontsize=7, ncol=2)
        ax[2].set_title(f"central point, C_both: a = ({', '.join(f'{a:+.2f}' for a in r['a'])})\n"
                        f"R = {r['R']:.2f}, χ²_RT/N = {r['chi2_RT_N']:.1f}, χ²_res/N = {r['chi2_res_N']:.1f}, class {r['cls']}",
                        fontsize=8)
    else:
        ax[2].axis("off")
    save(fig, out, "fig4_vectors")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--central", default="0.01,0.1,0.01")
    ap.add_argument("--out", default=str(HERE.parents[0] / "figures"))
    ap.add_argument("--which", default="2,3,4")
    ap.add_argument("--sens", default=str(HERE / "sensitivity.json"))
    a = ap.parse_args()
    m, v, x = (float(s) for s in a.central.split(","))
    out = Path(a.out)
    which = a.which.split(",")
    if "2" in which:
        fig2(HERE / "grid" / f"{model_name(m, v, x)}.json", out)
    if "3" in which:
        fig3(out)
    if "4" in which:
        fig4(Path(a.sens), str((m, v, x)), out)
