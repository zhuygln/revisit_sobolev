"""Phase 3A (§4.42): does the closure residual survive a real observation?

Three pre-declared observing scenarios at 40 Mpc on the grid's own epochs,
a magnitude-dependent noise model, and the Gate 2 projection of
`sensitivity.py` re-run with the scenario's observables and sigmas:

    sigma(m) = sqrt(sigma_sys^2 + (1.0857 / SNR)^2),  SNR = 5 * 10^(-0.4 (m - m_5sigma))

An observation is used only if the *reference* magnitude is brighter than
the 5-sigma depth and the (band, epoch) passes the Gate 2 masks (>= 1 % of
L_bol, floor). The reference sets sigma for every leg, so all legs carry the
same weights. chi2_RT,obs = sum (d_RT / sigma)^2 is the expected chi-square
of the closure error against the scenario's errors (analytic, Gaussian, no
noise realizations; the within-epoch correlation caveat of §4.40 stands).

Gate 3, pre-declared: at points with N_obs >= 8, the residual survives if
chi2_RT,obs / N_obs >= 4 and the class is C-B under the tangent space.

    python observe.py                 # writes observability.json, prints the table
    python observe.py --table         # table from the JSON
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
P12 = HERE.parents[0] / "phase12_grid"
sys.path.insert(0, str(P12))
import sensitivity as S      # noqa: E402

SCENARIOS = {
    # name: {"obs": {band: [epochs]}, "depth": {band: m_5sigma}, "sys": {band: sigma_sys}}
    "dense": {"obs": {b: [0.5, 1.0, 2.0, 3.0, 5.0, 7.0] for b in S.BANDS},
              "depth": {**{b: 23.5 for b in "griz"}, **{b: 21.5 for b in "JHK"}},
              "sys": {**{b: 0.03 for b in "griz"}, **{b: 0.05 for b in "JHK"}},
              "note": "AT2017gfo-like: griz + JHK at every grid epoch"},
    "sparse": {"obs": {**{b: [1.0, 3.0, 7.0] for b in "griz"}, **{b: [2.0, 5.0] for b in "JHK"}},
               "depth": {**{b: 22.5 for b in "griz"}, **{b: 20.5 for b in "JHK"}},
               "sys": {**{b: 0.03 for b in "griz"}, **{b: 0.05 for b in "JHK"}},
               "note": "typical follow-up: optical at 1, 3, 7 d, NIR at 2, 5 d, 1 mag shallower"},
    "optical": {"obs": {b: [0.5, 1.0, 2.0, 3.0, 5.0, 7.0] for b in "griz"},
                "depth": {b: 23.5 for b in "griz"}, "sys": {b: 0.03 for b in "griz"},
                "note": "no NIR"},
}
TANGENTS = ("T0", "T1", "T2", "T3")
LEGS = ("C_both", "C_binned")
N_OBS_MIN = 8
GATE = {"chi2_RT_N_min": 4.0, "cls": "C-B"}
NIR = ("J", "H", "K")


def sigma_of(m, m5, sys_):
    """Photometric error at magnitude m for a 5-sigma depth m5 and floor sys_."""
    snr = 5.0 * 10 ** (-0.4 * (m - m5))
    return float(np.sqrt(sys_ ** 2 + (1.0857 / snr) ** 2))


def scenario_obs(scn, keys, m_ref):
    """(band, t) the scenario provides with the reference above its depth, and their sigmas."""
    sc = SCENARIOS[scn]
    mask, sig = set(), {}
    for k in keys:
        b, t = k
        if b in sc["obs"] and t in sc["obs"][b] and np.isfinite(m_ref[k]) and m_ref[k] <= sc["depth"][b]:
            mask.add(k); sig[k] = sigma_of(m_ref[k], sc["depth"][b], sc["sys"][b])
    return mask, sig


def analyse(vecs, point, leg, scn, tangent, noise_floor):
    keys, m_ref = vecs[point][0], vecs[point][1]
    mask, sig = scenario_obs(scn, keys, m_ref)
    p = S.analyse_point(vecs, point, leg, noise_floor=noise_floor, sigma=sig, obs_mask=mask,
                        nuisance=S.TANGENT[tangent])
    out = {"scenario": scn, "tangent": tangent, "n_scenario": len(mask)}
    if p.get("status") != "ok":
        out.update(status=p["status"], N_obs=p.get("N", 0))
        return out
    d, s = np.array(p["d_rt"]), np.array(p["sigma"])
    chi2 = float(np.sum((d / s) ** 2))
    nir = np.array([k[0] in NIR for k in p["keys"]])
    out.update(status="ok", N_obs=p["N"], N_nir=int(nir.sum()), chi2_RT_obs=chi2, chi2_RT_N=p["chi2_RT_N"],
               chi2_res_dof=p["chi2_res_dof"], R=p["R"], dof=p["dof"], cls=p["cls"],
               underdetermined=p["underdetermined"], a=p["a"][:S.N_PHYS], a_nuisance=p["a_nuisance"],
               nir_share=float(np.sum((d[nir] / s[nir]) ** 2) / chi2) if chi2 > 0 else np.nan,
               sigma_median=float(np.median(s)), keys=[list(k) for k in p["keys"]], sigma=p["sigma"],
               d_rt=p["d_rt"], eligible=p["N"] >= N_OBS_MIN)
    out["detectable"] = out["chi2_RT_N"] >= GATE["chi2_RT_N_min"]
    out["survives"] = bool(out["eligible"] and out["detectable"] and out["cls"] == GATE["cls"])
    out["leftover_detectable"] = bool(np.isfinite(p["chi2_res_dof"]) and p["chi2_res_dof"] > S.THRESH["chi2_small"])
    return out


def summarize(results):
    summary = {}
    for leg in LEGS:
        for scn in SCENARIOS:
            for tg in TANGENTS:
                rs = [r["cells"][leg][scn][tg] for r in results.values()]
                ok = [r for r in rs if r["status"] == "ok"]
                el = [r for r in ok if r["eligible"]]
                det = [r for r in el if not r["underdetermined"]]
                s = {"analysable": len(ok), "eligible": len(el), "determined": len(det),
                     "detectable": sum(r["detectable"] for r in el),
                     "survives": sum(r["survives"] for r in el),
                     "leftover_detectable": sum(r["leftover_detectable"] for r in det),
                     "underdetermined": sum(r["underdetermined"] for r in el),
                     "median_N_obs": float(np.median([r["N_obs"] for r in el])) if el else None,
                     "median_chi2_RT_N": float(np.median([r["chi2_RT_N"] for r in el])) if el else None,
                     "median_chi2_res_dof": float(np.median([r["chi2_res_dof"] for r in det])) if det else None,
                     "median_R": float(np.median([r["R"] for r in det])) if det else None,
                     "median_nir_share": float(np.median([r["nir_share"] for r in el])) if el else None,
                     "by_X": {}}
                for X in S.X_GRID:
                    sel = [r for k, r in ((k, r["cells"][leg][scn][tg]) for k, r in results.items())
                           if eval(k)[2] == X and r["status"] == "ok" and r["eligible"]]
                    s["by_X"][str(X)] = {"eligible": len(sel), "survives": sum(r["survives"] for r in sel),
                                         "detectable": sum(r["detectable"] for r in sel)}
                summary.setdefault(leg, {}).setdefault(scn, {})[tg] = s
    return summary


def main(grid_dir, out, sens_json=None):
    models = S.load_grid(grid_dir)
    vecs = {p: S.vectors(m, S.FRAC_MIN) for p, m in models.items()}
    results = {}
    for point in vecs:
        floor = S.noise_floor_of(vecs, point)
        results[str(point)] = {"noise_floor": floor, "cells": {
            leg: {scn: {tg: analyse(vecs, point, leg, scn, tg, floor) for tg in TANGENTS}
                  for scn in SCENARIOS} for leg in LEGS}}
    d = {"scenarios": SCENARIOS, "tangents": TANGENTS, "n_obs_min": N_OBS_MIN, "gate": GATE,
         "frac_min": S.FRAC_MIN, "floored": "exclude", "core": "conserving",
         "points": results, "summary": summarize(results)}
    Path(out).write_text(json.dumps(d, indent=1))
    print_table(d)
    print(f"wrote {out}")
    return d


def print_table(d, leg="C_both"):
    print(f"\nGate 3 ({leg}): eligible = N_obs >= {d['n_obs_min']}; survives = chi2_RT,obs/N >= "
          f"{d['gate']['chi2_RT_N_min']} and class {d['gate']['cls']}")
    print("| scenario | tangent | analysable | eligible | detectable | survives | underdet. | leftover > 4 (of determined) | "
          "median N_obs | median χ²_RT,obs/N | median R | median χ²_res/dof | median NIR share | survives by X (10⁻³/10⁻²/10⁻¹) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for scn in d["scenarios"]:
        for tg in d["tangents"]:
            s = d["summary"][leg][scn][tg]
            f = lambda v, fmt: (fmt % v) if v is not None else "—"
            bx = " / ".join(f"{s['by_X'][str(X)]['survives']}/{s['by_X'][str(X)]['eligible']}" for X in S.X_GRID)
            print(f"| {scn} | {tg} | {s['analysable']} | {s['eligible']} | {s['detectable']} | **{s['survives']}** | "
                  f"{s['underdetermined']} | {s['leftover_detectable']}/{s['determined']} | {f(s['median_N_obs'], '%.0f')} | "
                  f"{f(s['median_chi2_RT_N'], '%.0f')} | {f(s['median_R'], '%.2f')} | {f(s['median_chi2_res_dof'], '%.1f')} | "
                  f"{f(s['median_nir_share'], '%.2f')} | {bx} |")
    print(f"\nper point ({leg}, dense / sparse / optical; T0 and T1): N_obs, χ²_RT,obs/N, R, class")
    for pt, r in d["points"].items():
        cells = r["cells"][leg]
        row = []
        for scn in d["scenarios"]:
            for tg in ("T0", "T1"):
                q = cells[scn][tg]
                row.append(f"{q['N_obs']}: {q['chi2_RT_N']:.0f} / {q['R']:.2f} / {q['cls']}" if q["status"] == "ok"
                           else f"{q.get('N_obs', 0)}: —")
        print(f"| {pt} | " + " | ".join(row) + " |")


def fig5(d, out, leg="C_both", central="(0.01, 0.1, 0.01)"):
    """Survival / detectability per scenario x tangent space, the per-band share
    of chi2_RT,obs under `dense`, and the central point's d_RT/sigma."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(P12))
    from figures import save
    scns, tgs = list(d["scenarios"]), list(d["tangents"])
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), gridspec_kw={"width_ratios": [1.4, 1.0, 1.2]})
    w = 0.8 / len(tgs)
    for j, tg in enumerate(tgs):
        x = np.arange(len(scns)) + (j - (len(tgs) - 1) / 2) * w
        surv = [d["summary"][leg][s][tg]["survives"] / max(d["summary"][leg][s][tg]["eligible"], 1) for s in scns]
        left = [d["summary"][leg][s][tg]["leftover_detectable"] / max(d["summary"][leg][s][tg]["determined"], 1) for s in scns]
        ax[0].bar(x, surv, w, color=f"C{j}", label=f"{tg}: survives (C-B, χ²/N ≥ 4)")
        ax[0].plot(x, left, "k_", ms=12, mew=1.5)
    ax[0].plot([], [], "k_", ms=12, mew=1.5, label="leftover χ²_res/dof > 4 (of determined)")
    ax[0].set_xticks(range(len(scns))); ax[0].set_xticklabels([f"{s}\n(N_obs ≥ {d['n_obs_min']})" for s in scns])
    ax[0].set_ylim(0, 1.05); ax[0].set_ylabel(f"fraction of eligible points ({leg})")
    ax[0].legend(fontsize=7, loc="upper right"); ax[0].set_title("Gate 3 by scenario and tangent space", fontsize=9)
    # per-band share of chi2_RT,obs, dense, eligible points
    shares = {b: [] for b in S.BANDS}
    for pt, r in d["points"].items():
        q = r["cells"][leg]["dense"]["T0"]
        if q["status"] != "ok" or not q["eligible"]:
            continue
        c = np.array([(dd / ss) ** 2 for dd, ss in zip(q["d_rt"], q["sigma"])])
        for b in S.BANDS:
            shares[b].append(sum(c[i] for i, k in enumerate(q["keys"]) if k[0] == b) / c.sum())
    pos = np.arange(len(S.BANDS))
    ax[1].boxplot([shares[b] for b in S.BANDS], positions=pos, widths=0.6, showfliers=True)
    ax[1].set_xticks(pos); ax[1].set_xticklabels(S.BANDS)
    ax[1].set_ylabel("share of χ²_RT,obs (dense, T0)"); ax[1].set_title("which band carries the χ²", fontsize=9)
    ax[1].axhline(1 / len(S.BANDS), color="0.6", lw=0.5, ls="--")
    # central point: d_RT / sigma per observable, dense
    q = d["points"].get(central, {}).get("cells", {}).get(leg, {}).get("dense", {}).get("T0", {})
    if q.get("status") == "ok":
        keys = [tuple(k) for k in q["keys"]]
        cmap = {b: f"C{i}" for i, b in enumerate(S.BANDS)}
        for b in S.BANDS:
            idx = [i for i, k in enumerate(keys) if k[0] == b]
            if idx:
                ax[2].plot([keys[i][1] for i in idx], [q["d_rt"][i] / q["sigma"][i] for i in idx], "o-", color=cmap[b], label=b, ms=4)
        ax[2].axhline(0, color="k", lw=0.5); ax[2].set_xlabel("t (d)"); ax[2].set_ylabel("d_RT / σ")
        ax[2].legend(fontsize=7, ncol=2)
        ax[2].set_title(f"{central}, dense: N_obs = {q['N_obs']}, χ²/N = {q['chi2_RT_N']:.0f}, {q['cls']}", fontsize=9)
    else:
        ax[2].axis("off")
    save(fig, Path(out), "fig5_observability")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid-dir", default=str(P12 / "grid"))
    ap.add_argument("--out", default=str(HERE / "observability.json"))
    ap.add_argument("--table", action="store_true")
    ap.add_argument("--leg", default="C_both")
    ap.add_argument("--fig", action="store_true", help="fig5 from the JSON")
    ap.add_argument("--fig-dir", default=str(HERE.parents[0] / "figures"))
    a = ap.parse_args()
    if a.table:
        print_table(json.loads(Path(a.out).read_text()), a.leg)
    elif a.fig:
        fig5(json.loads(Path(a.out).read_text()), a.fig_dir, a.leg)
    else:
        main(a.grid_dir, a.out)
