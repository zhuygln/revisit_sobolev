"""§4.43: the measured T_eff direction against the Planck proxy (F47).

Four model runs of the central model with the launch temperature scaled
(`grid.py --t-scale 0.8 / 1.25`, each with and without `--t-scale-gas`)
give the MC direction at fixed L per (band, epoch):

    d_T^MC(b, t) = [m_ref(1.25) - m_ref(0.8) + 2.5 log10(L_bol(1.25) / L_bol(0.8))] / ln(1.25 / 0.8)

(the grey luminosity term removed, so the direction is at fixed L like the
proxy `sensitivity.tbb_derivative`). Reports the cosine and norm ratio with
the proxy over the live observables, and reclassifies the central point
under T2/T3 with d_T^MC in place of the proxy column.

    python tscale.py            # reads grid/tscale/*.json, writes tscale.json
    python tscale.py --fig      # also paper3/figures/fig7_tscale.{png,pdf}
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sensitivity as S              # noqa: E402
from grid import model_name          # noqa: E402

CENTRAL = (0.01, 0.1, 0.01)
SCALES = (0.8, 1.25)


def load_pair(tdir, m, v, x, gas):
    out = {}
    for sc in SCALES:
        name = model_name(m, v, x, sc) + ("_gas" if gas else "")
        p = Path(tdir) / f"{name}.json"
        if not p.exists():
            return None
        out[sc] = json.loads(p.read_text())
    return out


def mc_direction(pair):
    """{(band, t): dm/dlnT at fixed L} from the two scaled runs, plus checks."""
    lo, hi = pair[SCALES[0]], pair[SCALES[1]]
    dl = np.log(SCALES[1] / SCALES[0])
    d, checks = {}, {}
    rows_lo = {r["t_d"]: r for r in lo["rows"] if r.get("status") in ("ok", "reduced_n")}
    rows_hi = {r["t_d"]: r for r in hi["rows"] if r.get("status") in ("ok", "reduced_n")}
    for t in sorted(set(rows_lo) & set(rows_hi)):
        a, b = rows_lo[t], rows_hi[t]
        grey = 2.5 * np.log10(b["ref"]["L_bol"] / a["ref"]["L_bol"])
        checks[t] = {"L_ratio": b["ref"]["L_bol"] / a["ref"]["L_bol"], "grey_mag": grey,
                     "T_eff": (a["source"]["T_eff"], b["source"]["T_eff"]),
                     "T_gas": (a.get("T_gas"), b.get("T_gas")),
                     "same_L_R_v": all(abs(a["source"][k] - b["source"][k]) < 1e-9 * abs(a["source"][k]) for k in ("L", "R_ph", "v_ph")),
                     "n_used": (a.get("n_used"), b.get("n_used"))}
        for band in S.BANDS:
            ma, mb = a["ref"]["mags"].get(band, np.nan), b["ref"]["mags"].get(band, np.nan)
            d[(band, t)] = (mb - ma + grey) / dl
    return d, checks


def compare(vecs, point, d_mc, leg="C_both"):
    """Cosine / norm ratio of d_mc against the proxy on the live observables, per epoch and overall."""
    keys, m0, d_rt, mask, info = vecs[point]
    live = [k for k in keys if mask[k] and k in d_mc and np.isfinite(d_mc[k])]
    proxy = {k: S.tbb_derivative(info["T_eff"][k[1]])[k[0]] for k in live}
    w = np.array([1.0 / S.SIGMA[k[0]] ** 2 for k in live])
    a = np.array([d_mc[k] for k in live]); b = np.array([proxy[k] for k in live])
    def cos(a, b, w):
        den = np.sqrt(np.sum(w * a * a) * np.sum(w * b * b))
        return float(np.sum(w * a * b) / den) if den > 0 else np.nan
    out = {"N": len(live), "cos": cos(a, b, w),
           "norm_ratio": float(np.sqrt(np.sum(w * a * a) / np.sum(w * b * b))) if np.sum(w * b * b) > 0 else np.nan,
           "per_epoch": {}, "mc": {f"{k[0]},{k[1]:g}": d_mc[k] for k in live},
           "proxy": {f"{k[0]},{k[1]:g}": proxy[k] for k in live}}
    for t in sorted({k[1] for k in live}):
        idx = [i for i, k in enumerate(live) if k[1] == t]
        out["per_epoch"][f"{t:g}"] = {"N": len(idx), "cos": cos(a[idx], b[idx], w[idx]),
                                       "norm_ratio": float(np.sqrt(np.sum(w[idx] * a[idx] ** 2) / np.sum(w[idx] * b[idx] ** 2)))}
    # cosine of d_RT with each direction
    r = np.array([d_rt[leg][k] for k in live])
    out["cos_dRT_mc"] = cos(r, a, w); out["cos_dRT_proxy"] = cos(r, b, w)
    return out


def main(grid_dir, tdir, out):
    models = S.load_grid(grid_dir)
    vecs = {p: S.vectors(m, S.FRAC_MIN) for p, m in models.items()}
    floor = S.noise_floor_of(vecs, CENTRAL)
    res = {"point": CENTRAL, "scales": SCALES, "noise_floor": floor, "variants": {}}
    for gas in (False, True):
        pair = load_pair(tdir, *CENTRAL, gas)
        tag = "with_T_gas" if gas else "illumination_only"
        if pair is None:
            res["variants"][tag] = {"status": "runs not present"}; continue
        n_done = min(len(pair[sc]["rows"]) for sc in SCALES)
        if n_done < len(pair[SCALES[0]]["epochs"]):
            res["variants"][tag] = {"status": f"runs incomplete ({n_done} epochs)"}; continue
        d_mc, checks = mc_direction(pair)
        v = {"status": "ok", "checks": {f"{t:g}": c for t, c in checks.items()},
             "compare": {leg: compare(vecs, CENTRAL, d_mc, leg) for leg in ("C_both", "C_binned")}, "reclass": {}}
        for tg in ("T2", "T3"):
            for leg in ("C_both", "C_binned"):
                base = S.analyse_point(vecs, CENTRAL, leg, noise_floor=floor, nuisance=S.TANGENT[tg])
                alt = S.analyse_point(vecs, CENTRAL, leg, noise_floor=floor, nuisance=S.TANGENT[tg],
                                      column_override={"T_bb": d_mc})
                pick = lambda q: ({k: q[k] for k in ("N", "dof", "R", "chi2_res_dof", "cls")} | {"a_T": q["a_nuisance"].get("T_bb")}
                                  if q.get("status") == "ok" else {"status": q["status"]})
                v["reclass"][f"{tg}/{leg}"] = {"proxy": pick(base), "mc": pick(alt)}
        res["variants"][tag] = v
    Path(out).write_text(json.dumps(res, indent=1))
    print_table(res)
    print(f"wrote {out}")
    return res


def fig7(res, out_dir, noise_floor=None):
    """d_T per band and epoch: MC (illumination only / with T_gas) against the Planck proxy."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    tags = [t for t in ("illumination_only", "with_T_gas") if res["variants"][t]["status"] == "ok"]
    fig, axes = plt.subplots(1, len(tags), figsize=(4.6 * len(tags), 3.8), sharey=True, squeeze=False)
    for ax, tag in zip(axes[0], tags):
        c = res["variants"][tag]["compare"]["C_both"]
        epochs = sorted({float(k.split(",")[1]) for k in c["mc"]})
        for j, t in enumerate(epochs):
            bands = [b for b in S.BANDS if f"{b},{t:g}" in c["mc"]]
            x = [S.BANDS.index(b) for b in bands]
            col = plt.cm.viridis(j / max(len(epochs) - 1, 1))
            ax.plot(x, [c["proxy"][f"{b},{t:g}"] for b in bands], "-", color=col, lw=1.2, alpha=0.8)
            ax.plot(x, [c["mc"][f"{b},{t:g}"] for b in bands], "o", color=col, ms=5, label=f"{t:g} d")
        if noise_floor:
            ax.axhspan(-noise_floor / np.log(SCALES[1] / SCALES[0]), noise_floor / np.log(SCALES[1] / SCALES[0]),
                       color="0.85", zorder=0)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(S.BANDS))); ax.set_xticklabels(S.BANDS)
        ax.set_title(f"{tag.replace('_', ' ')}: cos {c['cos']:.2f}, |MC|/|proxy| {c['norm_ratio']:.2f}", fontsize=9)
    axes[0][0].set_ylabel(r"$\partial m_b / \partial \ln T$ at fixed $L$ (mag)")
    axes[0][-1].legend(fontsize=7, title="MC (points), proxy (lines)", title_fontsize=7, loc="best")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"fig7_tscale.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)


def print_table(res):
    for tag, v in res["variants"].items():
        print(f"\n== {tag}: {v['status']}")
        if v["status"] != "ok":
            continue
        print("| t (d) | T_eff (0.8 / 1.25) | T_gas | L ratio | grey term (mag) | same L, R_ph, v_ph | n_used |")
        print("|---|---|---|---|---|---|---|")
        for t, c in v["checks"].items():
            print(f"| {t} | {c['T_eff'][0]:.0f} / {c['T_eff'][1]:.0f} | {c['T_gas'][0]:.0f} / {c['T_gas'][1]:.0f} | {c['L_ratio']:.3f} | "
                  f"{c['grey_mag']:+.3f} | {c['same_L_R_v']} | {c['n_used'][0]} / {c['n_used'][1]} |")
        for leg, c in v["compare"].items():
            print(f"{leg}: N = {c['N']}, cos(MC, proxy) = {c['cos']:.3f}, |MC|/|proxy| = {c['norm_ratio']:.2f}, "
                  f"cos(d_RT, MC) = {c['cos_dRT_mc']:.2f}, cos(d_RT, proxy) = {c['cos_dRT_proxy']:.2f}")
            print("   per epoch: " + "; ".join(f"{t} d: N={e['N']} cos={e['cos']:.2f} ratio={e['norm_ratio']:.2f}" for t, e in c["per_epoch"].items()))
        print("| space / leg | proxy: N / dof / R / χ²_res/dof / class / a_T | MC direction: same |")
        print("|---|---|---|")
        for k, r in v["reclass"].items():
            f = lambda q: (f"{q['N']} / {q['dof']} / {q['R']:.2f} / {q['chi2_res_dof']:.1f} / {q['cls']} / {q['a_T']:+.2f}"
                           if "cls" in q else q["status"])
            print(f"| {k} | {f(r['proxy'])} | {f(r['mc'])} |")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid-dir", default=str(HERE / "grid"))
    ap.add_argument("--tscale-dir", default=str(HERE / "grid" / "tscale"))
    ap.add_argument("--out", default=str(HERE / "tscale.json"))
    ap.add_argument("--fig", action="store_true", help="also draw paper3/figures/fig7_tscale")
    a = ap.parse_args()
    res = main(a.grid_dir, a.tscale_dir, a.out)
    if a.fig:
        fig7(res, HERE.parent / "figures", noise_floor=res.get("noise_floor"))
