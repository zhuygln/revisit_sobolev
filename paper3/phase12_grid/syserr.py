"""Paper III §4.45: the closure error against the model-error allowance.

Kilonova inference commonly adds a generic, uncorrelated "systematic" of
σ_sys ≈ 0.5–1 mag per photometric point to absorb model error. This module
asks what the *specific* closure error d_RT(band, epoch) of the grouped-opacity
transport looks like against that allowance, using only the observable vectors
Gate 2 already stored (sensitivity.json, T0): no transport is run.

Everything is unweighted and in magnitudes, because the allowance is; the
sensitivity σ's play no role here.

Per leg (C_both, C_binned, and the redistribution-only control A_redist):
  1. per (band, epoch) key: n live points, median, 16–84 % of d_RT; the fraction
     of live observables with |d_RT| > 0.5 and > 1 mag;
  2. coherence: (a) the fraction of coepochal live (g, K) pairs with
     d_RT(g) < 0 and d_RT(K) > 0 (the too-blue signature: dm = m_leg − m_ref);
     (b) the one-mode fraction f1 of a masked rank-1 fit X ≈ u vᵀ over all live
     entries (fraction of the squared Frobenius norm the first mode carries;
     the matrix is not centred -- the mean *is* the coherent mode), with the
     mode shape v; a median-filled SVD on keys live at >= `min_n` points as a
     cross-check; (c) two nulls: A_redist through the identical construction,
     and a sign-scrambled copy of the leg (each entry × ±1, `n_draws` draws)
     that keeps the amplitudes and the mask and destroys the sign coherence;
  3. against the allowance: Σ(d_RT/σ_sys)²/N per point at σ_sys = 1 and 0.5 mag
     (grid median and range), and the per-point fraction of ‖d_RT‖² along the
     first mode.
Optionally the same on the residual d_RT − fit of a tangent-space file (T1).

Usage: python syserr.py [--sens sensitivity.json] [--out syserr.json] [--residual sensitivity_T1.json]
"""
import argparse, ast, json, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from sensitivity import BANDS   # noqa: E402

LEGS = ("C_both", "C_binned", "A_redist")
SIGMAS = (1.0, 0.5)
THRESHOLDS = (0.5, 1.0)
MIN_N = 10


def matrix(sens, leg, field="d_rt"):
    """Points × keys matrix of the observable vector, NaN where a key is not live.
    field = "d_rt" (the closure error) or "residual" (d_rt − fit)."""
    pts = sorted(sens["points"], key=lambda s: ast.literal_eval(s))
    keys = sorted({tuple(k) for p in pts for k in sens["points"][p]["legs"][leg].get("keys", [])},
                  key=lambda k: (k[1], BANDS.index(k[0])))
    col = {k: j for j, k in enumerate(keys)}
    X = np.full((len(pts), len(keys)), np.nan)
    for i, p in enumerate(pts):
        q = sens["points"][p]["legs"][leg]
        if q.get("status") != "ok":
            continue
        vals = np.array(q["d_rt"]) if field == "d_rt" else np.array(q["d_rt"]) - np.array(q["fit"])
        for k, v in zip(q["keys"], vals):
            X[i, col[tuple(k)]] = v
    return X, pts, keys


def rank1_masked(X, iters=200, tol=1e-13):
    """Least-squares rank-1 fit X ≈ u vᵀ on the finite entries (alternating
    normal equations). Returns f1 = 1 − Σ_obs(X − uvᵀ)²/Σ_obs X², u, v with v unit-norm."""
    M = np.isfinite(X)
    X0 = np.where(M, X, 0.0)
    v = np.nanmedian(np.where(M, X, np.nan), axis=0)
    v = np.where(np.isfinite(v), v, 0.0)
    if not np.any(v):
        v = np.ones(X.shape[1])
    tot = float(np.sum(X0 ** 2))
    prev = None
    for _ in range(iters):
        den = M @ (v ** 2)
        u = np.where(den > 0, (X0 @ v) / np.where(den > 0, den, 1), 0.0)
        den = M.T @ (u ** 2)
        v = np.where(den > 0, (X0.T @ u) / np.where(den > 0, den, 1), 0.0)
        res = float(np.sum(((X0 - np.outer(u, v)) * M) ** 2))
        if prev is not None and abs(prev - res) <= tol * max(tot, 1e-300):
            break
        prev = res
    f1 = 1.0 - res / tot if tot > 0 else float("nan")
    nv = np.linalg.norm(v)
    if nv > 0:
        u, v = u * nv, v / nv
    if np.sum(v) < 0:          # sign convention: the mode's mean is positive
        u, v = -u, -v
    return f1, u, v


def svd_filled(X, min_n=MIN_N):
    """Cross-check: keep keys live at >= min_n points, fill gaps with the key median,
    return the first-singular-value fraction of the squared Frobenius norm."""
    M = np.isfinite(X)
    cols = np.where(M.sum(axis=0) >= min_n)[0]
    if len(cols) == 0:
        return float("nan"), 0
    Y = X[:, cols].copy()
    med = np.nanmedian(Y, axis=0)
    Y = np.where(np.isfinite(Y), Y, med)
    s = np.linalg.svd(Y, compute_uv=False)
    return float(s[0] ** 2 / np.sum(s ** 2)), int(len(cols))


def sign_scramble_null(X, n_draws=1000, seed=0):
    rng = np.random.default_rng(seed)
    f = []
    for _ in range(n_draws):
        f.append(rank1_masked(X * rng.choice([-1.0, 1.0], size=X.shape))[0])
    f = np.array(f)
    return {"n_draws": n_draws, "seed": seed, "median": float(np.median(f)), "p95": float(np.percentile(f, 95)),
            "max": float(np.max(f))}


def mp_scale(m, n):
    """Marchenko–Pastur: the first singular value of an iid m×n matrix carries
    ≈ (√m+√n)²/(mn) of the squared Frobenius norm."""
    return (np.sqrt(m) + np.sqrt(n)) ** 2 / (m * n)


def threshold_fractions(X, thresholds=THRESHOLDS):
    v = np.abs(X[np.isfinite(X)])
    return {f"{t:g}": [int(np.sum(v > t)), int(v.size)] for t in thresholds}


def sign_pattern(X, keys, blue="g", red="K"):
    """(count, total) of coepochal live (blue, red) pairs with d_RT(blue) < 0 and d_RT(red) > 0."""
    col = {k: j for j, k in enumerate(keys)}
    epochs = sorted({k[1] for k in keys})
    n = tot = 0
    for t in epochs:
        if (blue, t) not in col or (red, t) not in col:
            continue
        a, b = X[:, col[(blue, t)]], X[:, col[(red, t)]]
        ok = np.isfinite(a) & np.isfinite(b)
        tot += int(ok.sum())
        n += int(np.sum((a[ok] < 0) & (b[ok] > 0)))
    return [n, tot]


def chi2_equiv(X, sigma):
    """Σ(d_RT/σ_sys)²/N per point (rows with no live entry are skipped)."""
    out = []
    for row in X:
        v = row[np.isfinite(row)]
        if v.size:
            out.append(float(np.sum((v / sigma) ** 2) / v.size))
    return out


def per_key(X, keys):
    out = []
    for j, k in enumerate(keys):
        v = X[:, j][np.isfinite(X[:, j])]
        out.append({"band": k[0], "t_d": k[1], "n": int(v.size),
                    "median": float(np.median(v)) if v.size else None,
                    "p16": float(np.percentile(v, 16)) if v.size else None,
                    "p84": float(np.percentile(v, 84)) if v.size else None})
    return out


def analyse(X, pts, keys, n_draws=1000, seed=0):
    f1, u, v = rank1_masked(X)
    f_svd, n_cols = svd_filled(X)
    M = np.isfinite(X)
    # per-point share of ‖d_RT‖² along the first mode
    share = []
    for i in range(X.shape[0]):
        m = M[i]
        if m.sum() == 0:
            continue
        x = X[i, m]; vv = v[m]
        c = float(np.dot(x, vv) ** 2 / (np.dot(x, x) * np.dot(vv, vv))) if np.dot(vv, vv) > 0 else 0.0
        share.append(c)
    chi2 = {f"{s:g}": chi2_equiv(X, s) for s in SIGMAS}
    return {"n_points": int(np.sum(M.any(axis=1))), "n_keys": len(keys), "n_live": int(M.sum()),
            "keys": [list(k) for k in keys], "per_key": per_key(X, keys),
            "frac_gt": threshold_fractions(X), "sign_pattern_gK": sign_pattern(X, keys),
            "one_mode": {"f1": f1, "mode": v.tolist(), "amplitude": u.tolist(), "points": pts,
                         "svd_filled": f_svd, "svd_n_keys": n_cols, "min_n": MIN_N,
                         "mp_scale": float(mp_scale(X.shape[0], n_cols)) if n_cols else None,
                         "point_share_median": float(np.median(share)) if share else None},
            "null_sign_scramble": sign_scramble_null(X, n_draws, seed),
            "chi2_equiv": {s: {"median": float(np.median(c)), "min": float(np.min(c)), "max": float(np.max(c)),
                               "per_point": c} for s, c in chi2.items()}}


def main(sens_path, out=None, residual=None, n_draws=1000, seed=0, legs=LEGS):
    sens = json.loads(Path(sens_path).read_text())
    d = {"source": Path(sens_path).name, "tangent": sens.get("tangent"), "sigmas": list(SIGMAS),
         "thresholds": list(THRESHOLDS), "note": "unweighted, magnitudes; dm = m_leg - m_ref", "legs": {}}
    for leg in legs:
        X, pts, keys = matrix(sens, leg)
        d["legs"][leg] = analyse(X, pts, keys, n_draws, seed)
    if residual:
        r = json.loads(Path(residual).read_text())
        d["residual"] = {"source": Path(residual).name, "tangent": r.get("tangent"), "legs": {}}
        for leg in legs:
            X, pts, keys = matrix(r, leg, field="residual")
            d["residual"]["legs"][leg] = analyse(X, pts, keys, n_draws, seed)
    if out:
        Path(out).write_text(json.dumps(d, indent=1))
    return d


def print_table(d):
    print("| leg | live | > 0.5 mag | > 1 mag | (g<0, K>0) pairs | one-mode f1 | SVD (≥10 pts) | scrambled null median / 95 % | χ²/N at 1 mag (median, range) | at 0.5 mag |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for leg, a in d["legs"].items():
        f5, f1 = a["frac_gt"]["0.5"], a["frac_gt"]["1"]
        c1, c5 = a["chi2_equiv"]["1"], a["chi2_equiv"]["0.5"]
        print(f"| {leg} | {a['n_live']} | {f5[0]/f5[1]:.2f} | {f1[0]/f1[1]:.2f} | {a['sign_pattern_gK'][0]}/{a['sign_pattern_gK'][1]} | "
              f"{a['one_mode']['f1']:.2f} | {a['one_mode']['svd_filled']:.2f} ({a['one_mode']['svd_n_keys']} keys) | "
              f"{a['null_sign_scramble']['median']:.2f} / {a['null_sign_scramble']['p95']:.2f} | "
              f"{c1['median']:.2f} ({c1['min']:.2f}–{c1['max']:.2f}) | {c5['median']:.2f} |")
    a = d["legs"]["C_both"]
    print("\nC_both per-key median (n ≥ 5):")
    for k in a["per_key"]:
        if k["n"] >= 5:
            print(f"  {k['band']} {k['t_d']:g} d  n={k['n']:2d}  median {k['median']:+.2f}  16–84 % {k['p16']:+.2f}…{k['p84']:+.2f}")
    if "residual" in d:
        r = d["residual"]["legs"]["C_both"]
        print(f"\nresidual after {d['residual']['tangent']} (C_both): one-mode {r['one_mode']['f1']:.2f}, "
              f"> 0.5 mag {r['frac_gt']['0.5'][0]}/{r['frac_gt']['0.5'][1]}, χ²/N at 1 mag median {r['chi2_equiv']['1']['median']:.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sens", default=str(HERE / "sensitivity.json"))
    ap.add_argument("--out", default=str(HERE / "syserr.json"))
    ap.add_argument("--residual", default=str(HERE / "sensitivity_T1.json"))
    ap.add_argument("--draws", type=int, default=1000)
    a = ap.parse_args()
    d = main(a.sens, a.out, residual=a.residual if Path(a.residual).exists() else None, n_draws=a.draws)
    print_table(d)
    print("wrote", a.out)
