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
    n_phys = 3   # cosines of the physical columns; nuisance amplitudes go in the text column
    for j, leg in enumerate(legs):
        M = np.full((len(pts), n_phys), np.nan)
        for i, p in enumerate(pts):
            r = d["points"][p]["legs"][leg]
            if r.get("status") == "ok":
                M[i] = r["cos"][:n_phys]
        im = ax[j].imshow(M, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
        ax[j].set_xticks(range(n_phys)); ax[j].set_xticklabels(["cos(d_RT, d_M)", "cos(d_RT, d_v)", "cos(d_RT, d_X)"], fontsize=8)
        ax[j].set_yticks(range(len(pts)))
        ax[j].set_yticklabels(pts if j == 0 else [], fontsize=6)
        ax[j].set_xlim(-0.5, n_phys + 0.6 + (0.6 if d.get("nuisance") else 0.0))
        for i, p in enumerate(pts):
            r = d["points"][p]["legs"][leg]
            txt = (r["cls"] + ("*" if r.get("unstable") else "")) if r.get("status") == "ok" else "—"
            if r.get("status") == "ok" and r.get("a_nuisance"):
                nz = {n: a for n, a in r["a_nuisance"].items() if not n.startswith("L_")}
                txt += " " + " ".join(f"{n}={a:+.1f}" for n, a in nz.items())
            ax[j].text(n_phys - 0.4, i, txt, fontsize=6, va="center")
        title = f"{leg}: cosines and class (* = unstable)"
        if d.get("nuisance"):
            title += f"\ntangent {d.get('tangent')}: + {', '.join(d['nuisance'])}"
        ax[j].set_title(title, fontsize=9)
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
        ax[2].set_title(f"central point, C_both: a = ({', '.join(f'{a:+.2f}' for a in r['a'][:n_phys])})\n"
                        f"R = {r['R']:.2f}, χ²_RT/N = {r['chi2_RT_N']:.1f}, χ²_res/N = {r['chi2_res_N']:.1f}, class {r['cls']}",
                        fontsize=8)
    else:
        ax[2].axis("off")
    save(fig, out, "fig4_vectors")


def fig6(sens_dir, central_key, out, leg="C_both"):
    """R and chi2_res/dof per point against the tangent space T0..T3 (§4.41),
    with the nuisance amplitudes the fits demand."""
    files = {"T0": "sensitivity.json", "T1": "sensitivity_T1.json", "T2": "sensitivity_T2.json", "T3": "sensitivity_T3.json"}
    D = {t: json.loads((Path(sens_dir) / f).read_text()) for t, f in files.items() if (Path(sens_dir) / f).exists()}
    tgs = list(D)
    pts = list(D[tgs[0]]["points"])
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for i, (key, lab, log) in enumerate((("R", "R = |residual| / |d_RT|", False), ("chi2_res_dof", "χ²_res / dof", True))):
        for pt in pts:
            y = []
            for t in tgs:
                q = D[t]["points"][pt]["legs"][leg]
                y.append(q[key] if q.get("status") == "ok" and not q.get("underdetermined") else np.nan)
            style = dict(color="C3", lw=2, marker="o", zorder=5) if pt == central_key else dict(color="0.5", lw=0.7, marker=".", alpha=0.7)
            ax[i].plot(range(len(tgs)), y, **style)
        med = [D[t]["summary"][leg].get("median_" + key) for t in tgs]
        ax[i].plot(range(len(tgs)), med, "k-", lw=2, marker="s", label="median")
        ax[i].axhline(0.3 if key == "R" else 4.0, color="C0", ls="--", lw=0.8, label="threshold")
        ax[i].set_xticks(range(len(tgs))); ax[i].set_xticklabels([f"{t}\n+{','.join(D[t]['nuisance']) or '—'}" for t in tgs], fontsize=8)
        ax[i].set_ylabel(lab)
        if log:
            ax[i].set_yscale("log")
        ax[i].legend(fontsize=7)
        ax[i].set_title(f"{leg}: {lab} vs tangent space\n(red = {central_key}; absent where underdetermined)", fontsize=8)
    # amplitudes demanded: |a_L| max per point (T1), a_T (T2), a_2c (T3)
    a_L = [max(abs(a) for a in q["a_nuisance"].values()) for q in (D["T1"]["points"][pt]["legs"][leg] for pt in pts)
           if q.get("status") == "ok" and not q.get("underdetermined") and q["a_nuisance"]] if "T1" in D else []
    a_T = [q["a_nuisance"]["T_bb"] for q in (D["T2"]["points"][pt]["legs"][leg] for pt in pts)
           if q.get("status") == "ok" and not q.get("underdetermined") and "T_bb" in q["a_nuisance"]] if "T2" in D else []
    a_2 = [q["a_nuisance"]["2c"] for q in (D["T3"]["points"][pt]["legs"][leg] for pt in pts)
           if q.get("status") == "ok" and not q.get("underdetermined") and "2c" in q["a_nuisance"]] if "T3" in D else []
    data = [v for v in (a_L, a_T, a_2) if v]
    labels = [l for v, l in zip((a_L, a_T, a_2), ("max |a_L| (mag, T1)", "Δln T (T2)", "a_2c (blue-model flux, T3)")) if v]
    ax[2].boxplot(data, widths=0.5, showfliers=True)
    ax[2].set_xticks(range(1, len(data) + 1)); ax[2].set_xticklabels(labels, fontsize=8)
    ax[2].axhline(0, color="k", lw=0.5)
    ax[2].set_title("nuisance amplitudes the fits demand (determined points)", fontsize=9)
    save(fig, out, "fig6_tangent")


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
    if "6" in which:
        fig6(Path(a.sens).parent, str((m, v, x)), out)
