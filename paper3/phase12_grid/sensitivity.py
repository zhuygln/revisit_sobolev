"""Workstream C: does the closure error point along a physical-parameter
direction? (Gate 2 of the 2026-09-02 plan.)

For every grid point and closure leg, the observable vector is the masked set
of (filter, epoch) reference magnitudes. `d_RT` is the closure's error in
that vector; `d_theta` are the finite-difference derivatives dm/dln theta for
theta in (M_ej, v_ej, X_lan) on the grid. The weighted least-squares
projection  a = (D'WD)^-1 D'W d_RT  is the linearized inference bias -- the
parameter shift a fit would absorb the closure error into -- with covariance
(D'WD)^-1, residual fraction R, chi^2 of the residual and of d_RT itself, and
the cosines between d_RT and each d_theta.

Pre-declared classification per (point, leg):
  C-C small      chi2_RT/N <= 4
  C-A degenerate not C-C, R <= 0.3 and max |a|/sqrt(Cov) >= 3
  C-B distinct   otherwise
Robustness: reclassified with sigma doubled, one-sided derivatives, and a 3 %
live-band cut; a point whose class changes is flagged unstable.

Noise floor: the A_redist leg's |dm| is the empirical MC floor; finite
differences amplify it by sqrt(2)/dln theta, and a derivative component below
3x that is zeroed.

Usage: python sensitivity.py [--grid-dir grid] [--out sensitivity.json]
"""
import argparse, json, sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_grid import M_GRID, V_GRID, X_GRID   # noqa: E402

BANDS = ("g", "r", "i", "z", "J", "H", "K")
LEGS = ("A_redist", "B_opacity", "C_both", "C_binned")
SIGMA = {"g": 0.05, "r": 0.05, "i": 0.05, "z": 0.05, "J": 0.10, "H": 0.10, "K": 0.10}
MAG_LIMIT = {"g": 23.5, "r": 23.5, "i": 23.5, "z": 23.5, "J": 21.5, "H": 21.5, "K": 21.5}
FRAC_MIN = 0.01
THETA = ("M", "v", "X")
GRIDS = {"M": M_GRID, "v": V_GRID, "X": X_GRID}
THRESH = {"chi2_small": 4.0, "R_max": 0.3, "signif_min": 3.0}
N_LOW = 8      # fewer observables than this (3 parameters) is flagged low_N, not excluded


# --- loading ----------------------------------------------------------------

def load_grid(grid_dir):
    from grid import model_name
    out = {}
    for m, v, x in product(M_GRID, V_GRID, X_GRID):
        p = Path(grid_dir) / f"{model_name(m, v, x)}.json"
        if p.exists():
            out[(m, v, x)] = json.loads(p.read_text())
    return out


def band_fraction(row, lam_window, n_spec):
    """Each band's share of the reference L_bol -- `summary.band_fractions`, on
    the real-filter extent (bin weights)."""
    sys.path.insert(0, str(HERE.parents[1]))
    from sobolev import photometry as phot
    edges = phot.nu_edges(*lam_window, n_spec)
    lnu = np.asarray(row["ref"]["L_nu"], float)
    tot = float(np.sum(lnu * np.diff(edges)))
    out = {}
    for b, pb in phot.load_passbands().items():
        w = pb.bin_weights(edges) > 0
        out[b] = float(np.sum((lnu * np.diff(edges))[w])) / tot if tot > 0 else 0.0
    return out


def vectors(model, frac_min=FRAC_MIN):
    """Reference magnitudes m[(band, epoch)], per-leg errors, and the mask."""
    epochs = model["epochs"]
    keys = [(b, t) for t in epochs for b in BANDS]
    m_ref = {k: np.nan for k in keys}
    d_rt = {leg: {k: np.nan for k in keys} for leg in LEGS}
    mask = {k: False for k in keys}
    for r in model["rows"]:
        if r.get("status") not in ("ok", "reduced_n"):
            continue
        t = r["t_d"]
        frac = band_fraction(r, model["lam_window"], model["n_spec"])
        for b in BANDS:
            mag = r["ref"]["mags"].get(b, np.nan)
            m_ref[(b, t)] = mag
            mask[(b, t)] = bool(np.isfinite(mag) and frac[b] >= frac_min and mag <= MAG_LIMIT[b]
                                and not r["source"].get("v_ph_floored", False))
            for leg in LEGS:
                d_rt[leg][(b, t)] = r["legs"][leg]["dm"].get(b, np.nan)
    return keys, m_ref, d_rt, mask


# --- derivatives --------------------------------------------------------------

def neighbours(point, theta):
    """Grid neighbours along theta: (lower, upper) or None at the edge."""
    i = THETA.index(theta)
    g = GRIDS[theta]
    j = g.index(point[i])
    lo = list(point); hi = list(point)
    lo[i] = g[j - 1] if j > 0 else None
    hi[i] = g[j + 1] if j < len(g) - 1 else None
    return (tuple(lo) if lo[i] is not None else None, tuple(hi) if hi[i] is not None else None)


def derivative(vecs, point, theta, scheme="central"):
    """dm/dln theta at `point` by finite difference on the grid.

    Returns (d, mask, info): central where both neighbours exist (the actual
    dln theta of the two-sided secant), one-sided otherwise, flagged.
    """
    lo, hi = neighbours(point, theta)
    i = THETA.index(theta)
    have = {p for p in (lo, hi, point) if p is not None and p in vecs}
    keys, m0, _, mask0 = vecs[point]
    if scheme == "central" and lo in have and hi in have:
        _, m_lo, _, mk_lo = vecs[lo]; _, m_hi, _, mk_hi = vecs[hi]
        dl = np.log(hi[i] / lo[i])
        d = {k: (m_hi[k] - m_lo[k]) / dl for k in keys}
        mask = {k: mask0[k] and mk_lo[k] and mk_hi[k] for k in keys}
        return d, mask, {"scheme": "central", "dln": dl, "one_sided": False}
    # one-sided (edge, or requested): prefer the upper neighbour
    for nb, sign in ((hi, 1.0), (lo, -1.0)):
        if nb in have:
            _, m_nb, _, mk_nb = vecs[nb]
            dl = sign * np.log(nb[i] / point[i])
            d = {k: (m_nb[k] - m0[k]) / dl for k in keys}
            mask = {k: mask0[k] and mk_nb[k] for k in keys}
            return d, mask, {"scheme": "one_sided", "dln": dl, "one_sided": True,
                             "side": "upper" if sign > 0 else "lower"}
    return None, None, {"scheme": "none"}


def secant_disagreement(vecs, point, theta):
    """The two one-sided secants and their normalized disagreement: the
    nonlinearity metric the report must quote for the x10 X spacing."""
    lo, hi = neighbours(point, theta)
    if lo not in vecs or hi not in vecs:
        return None
    i = THETA.index(theta)
    keys, m0, _, mask0 = vecs[point]
    _, m_lo, _, mk_lo = vecs[lo]; _, m_hi, _, mk_hi = vecs[hi]
    ks = [k for k in keys if mask0[k] and mk_lo[k] and mk_hi[k]]
    if not ks:
        return None
    du = np.array([(m_hi[k] - m0[k]) / np.log(hi[i] / point[i]) for k in ks])
    dd = np.array([(m0[k] - m_lo[k]) / np.log(point[i] / lo[i]) for k in ks])
    n = np.linalg.norm(du) + np.linalg.norm(dd)
    return {"n": len(ks), "norm_upper": float(np.linalg.norm(du)), "norm_lower": float(np.linalg.norm(dd)),
            "disagreement": float(2 * np.linalg.norm(du - dd) / n) if n > 0 else None}


# --- projection ----------------------------------------------------------------

def project(d_rt, D, sigma):
    """WLS of d_rt on the columns of D (N x p) with W = diag(1/sigma^2)."""
    W = 1.0 / sigma ** 2
    G = D.T @ (W[:, None] * D)
    cond = float(np.linalg.cond(G)) if np.all(np.isfinite(G)) else np.inf
    cov = np.linalg.pinv(G)
    a = cov @ (D.T @ (W * d_rt))
    res = d_rt - D @ a
    n_rt = np.sqrt(np.sum(W * d_rt ** 2))
    R = float(np.sqrt(np.sum(W * res ** 2)) / n_rt) if n_rt > 0 else np.nan
    N = d_rt.size
    cos = []
    for j in range(D.shape[1]):
        dj = D[:, j]
        den = np.sqrt(np.sum(W * dj ** 2)) * n_rt
        cos.append(float(np.sum(W * dj * d_rt) / den) if den > 0 else np.nan)
    signif = np.abs(a) / np.sqrt(np.maximum(np.diag(cov), 1e-300))
    return {"a": a.tolist(), "cov": cov.tolist(), "signif": signif.tolist(), "R": R,
            "chi2_res_N": float(np.sum(W * res ** 2) / N), "chi2_RT_N": float(np.sum(W * d_rt ** 2) / N),
            "cos": cos, "cond": cond, "N": int(N)}


def classify(p, thresh=THRESH):
    if not np.isfinite(p["chi2_RT_N"]):
        return "undefined"
    if p["chi2_RT_N"] <= thresh["chi2_small"]:
        return "C-C"
    if p["R"] <= thresh["R_max"] and max(p["signif"]) >= thresh["signif_min"]:
        return "C-A"
    return "C-B"


def analyse_point(vecs, point, leg, sigma_scale=1.0, scheme="central", frac_min=None,
                  noise_floor=None):
    keys, m0, d_rt_all, mask0 = vecs[point]
    ders, masks, infos = [], [], []
    for th in THETA:
        d, mk, info = derivative(vecs, point, th, scheme)
        if d is None:
            return {"status": f"no derivative along {th}"}
        ders.append(d); masks.append(mk); infos.append(info)
    ks = [k for k in keys if mask0[k] and all(mk[k] for mk in masks)
          and np.isfinite(d_rt_all[leg][k]) and all(np.isfinite(d[k]) for d in ders)]
    if len(ks) < 4:
        return {"status": f"only {len(ks)} usable observables"}
    d_rt = np.array([d_rt_all[leg][k] for k in ks])
    D = np.array([[d[k] for d in ders] for k in ks])
    if noise_floor is not None:
        # zero derivative components below 3x the amplified MC floor
        for j, info in enumerate(infos):
            amp = np.sqrt(2.0) / abs(info["dln"])
            D[np.abs(D[:, j]) < 3.0 * noise_floor * amp, j] = 0.0
    sigma = np.array([SIGMA[k[0]] for k in ks]) * sigma_scale
    p = project(d_rt, D, sigma)
    p.update(status="ok", keys=[list(k) for k in ks], derivative_info=infos,
             d_rt_norm=float(np.sqrt(np.mean(d_rt ** 2))),
             d_rt=d_rt.tolist(), fit=(D @ np.array(p["a"])).tolist(),
             D=D.tolist(), cls=classify(p), low_N=len(ks) < N_LOW,
             live_params=[bool(np.any(D[:, j] != 0)) for j in range(D.shape[1])])
    return p


def noise_floor_of(vecs, point):
    keys, _, d_rt_all, mask = vecs[point]
    v = [abs(d_rt_all["A_redist"][k]) for k in keys if mask[k] and np.isfinite(d_rt_all["A_redist"][k])]
    return float(np.max(v)) if v else 0.0


def main(grid_dir, out):
    models = load_grid(grid_dir)
    print(f"{len(models)} of 27 models present")
    vecs1 = {p: vectors(m, FRAC_MIN) for p, m in models.items()}
    vecs3 = {p: vectors(m, 0.03) for p, m in models.items()}
    results = {}
    for point in vecs1:
        floor = noise_floor_of(vecs1, point)
        rec = {"noise_floor": floor, "secants": {th: secant_disagreement(vecs1, point, th) for th in THETA},
               "legs": {}}
        for leg in LEGS:
            base = analyse_point(vecs1, point, leg, noise_floor=floor)
            if base["status"] != "ok":
                rec["legs"][leg] = base; continue
            alt = {"sigma_x2": analyse_point(vecs1, point, leg, sigma_scale=2.0, noise_floor=floor),
                   "one_sided": analyse_point(vecs1, point, leg, scheme="one_sided", noise_floor=floor),
                   "frac_3pct": analyse_point(vecs3, point, leg, noise_floor=floor)}
            base["robustness"] = {k: v.get("cls", v["status"]) for k, v in alt.items()}
            base["unstable"] = any(v.get("cls") not in (None, base["cls"]) for v in alt.values())
            rec["legs"][leg] = base
        results[str(point)] = rec
    summary = {}
    for leg in LEGS:
        cls = [r["legs"][leg].get("cls") for r in results.values() if r["legs"][leg].get("status") == "ok"]
        summary[leg] = {c: cls.count(c) for c in ("C-A", "C-B", "C-C")}
        summary[leg]["n"] = len(cls)
        summary[leg]["unstable"] = sum(r["legs"][leg].get("unstable", False) for r in results.values())
        ok = [r["legs"][leg] for r in results.values() if r["legs"][leg].get("status") == "ok"]
        summary[leg]["low_N"] = sum(q["low_N"] for q in ok)
        summary[leg]["well_sampled"] = {c: sum(q["cls"] == c and not q["low_N"] for q in ok)
                                        for c in ("C-A", "C-B", "C-C")}
        summary[leg]["median_chi2_RT_N"] = float(np.median([q["chi2_RT_N"] for q in ok])) if ok else None
        summary[leg]["median_R"] = float(np.median([q["R"] for q in ok])) if ok else None
    d = {"thresholds": THRESH, "n_low": N_LOW, "sigma": SIGMA, "mag_limit": MAG_LIMIT, "frac_min": FRAC_MIN,
         "points": results, "summary": summary}
    Path(out).write_text(json.dumps(d, indent=1))
    print_table(d)
    print(f"wrote {out}")
    return d


def print_table(d):
    print(f"\n{'point':22s}{'leg':11s}{'N':>4}{'chi2RT/N':>10}{'R':>7}{'chi2res/N':>10}"
          f"{'a_M':>7}{'a_v':>7}{'a_X':>7}{'cosM':>6}{'cosv':>6}{'cosX':>6}  class")
    for pt, r in d["points"].items():
        for leg, p in r["legs"].items():
            if p.get("status") != "ok":
                print(f"{pt:22s}{leg:11s}  {p['status']}"); continue
            print(f"{pt:22s}{leg:11s}{p['N']:4d}{p['chi2_RT_N']:10.1f}{p['R']:7.2f}{p['chi2_res_N']:10.1f}"
                  + "".join(f"{a:7.2f}" for a in p["a"]) + "".join(f"{c:6.2f}" for c in p["cos"])
                  + f"  {p['cls']}{' *' if p.get('unstable') else ''}{' lowN' if p.get('low_N') else ''}")
    print("\nsummary:", json.dumps(d["summary"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid-dir", default=str(HERE / "grid"))
    ap.add_argument("--out", default=str(HERE / "sensitivity.json"))
    a = ap.parse_args()
    main(a.grid_dir, a.out)
