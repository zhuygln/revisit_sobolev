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

§4.41-§4.43 generalization (the design matrix is explicit):
  * `core` -- the magnitudes of every leg re-derived under the absorbing core
    from the stored L_bol / L_bol_absorbing (a grey per-leg factor, so colours
    are unchanged and the shift is exact up to seed averaging);
  * `floored` -- rows whose grey photosphere sits on the v_ej/2 floor are
    excluded (the declared Gate-2 rule; the committed run lacked the flag and
    included them) or included;
  * `nuisance` columns appended after (M, v, X): `L_t` one grey column per
    epoch (any luminosity history; also makes the class core-independent),
    `T_bb` dm/dlnT of a Planck spectrum at the row's T_eff at fixed L through
    the real passbands, `2c` a linearized second component with the same-(M, v)
    X = 1e-3 reference's light curve. Tangent spaces T0 = {M, v, X},
    T1 = T0 + L_t, T2 = T1 + T_bb, T3 = T2 + 2c. The class is decided on the
    physical columns' significance only; dof = N - rank, and a point with
    dof < 4 is `underdetermined`;
  * `sigma` / `obs_mask` -- per-observable weights and an observing mask
    (phase 13's scenarios).

Usage: python sensitivity.py [--grid-dir grid] [--out sensitivity.json]
       python sensitivity.py --core absorbing | --floored include | --tangent T1
       python sensitivity.py --override robustness/chain_<model>_t<t>.json --override-chain 8000
"""
import argparse, json, sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_grid import M_GRID, V_GRID, X_GRID   # noqa: E402
sys.path.insert(0, str(HERE.parents[1]))
from sobolev import photometry as phot        # noqa: E402
from sobolev.source import V_PH_MIN           # noqa: E402

BANDS = ("g", "r", "i", "z", "J", "H", "K")
LEGS = ("A_redist", "B_opacity", "C_both", "C_binned")
SIGMA = {"g": 0.05, "r": 0.05, "i": 0.05, "z": 0.05, "J": 0.10, "H": 0.10, "K": 0.10}
MAG_LIMIT = {"g": 23.5, "r": 23.5, "i": 23.5, "z": 23.5, "J": 21.5, "H": 21.5, "K": 21.5}
FRAC_MIN = 0.01
THETA = ("M", "v", "X")
GRIDS = {"M": M_GRID, "v": V_GRID, "X": X_GRID}
THRESH = {"chi2_small": 4.0, "R_max": 0.3, "signif_min": 3.0}
N_LOW = 8      # fewer observables than this (3 parameters) is flagged low_N, not excluded
N_PHYS = 3
DOF_MIN = 4    # N - rank below this: the fit is not a test of anything
NUISANCE = ("L_t", "T_bb", "2c")
TANGENT = {"T0": (), "T1": ("L_t",), "T2": ("L_t", "T_bb"), "T3": ("L_t", "T_bb", "2c")}
LAM_WIN, N_SPEC = (1000.0, 30000.0), 200


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


def row_floored(row, model):
    """Is the grey photosphere on the v_ej/2 floor? The pre-§4.41 JSONs lack
    the flag; v_core / v_ej == V_PH_MIN identifies the floored rows exactly."""
    src = row.get("source", {})
    if "v_ph_floored" in src:
        return bool(src["v_ph_floored"])
    return abs(row["v_core"] / model["v_ej_c"] - V_PH_MIN) < 1e-6


def core_shift(o, core):
    """m^core - m^conserving for a leg/ref dict: a grey magnitude shift."""
    if core == "conserving":
        return 0.0
    if core == "absorbing":
        return 2.5 * np.log10(o["L_bol"] / o["L_bol_absorbing"])
    raise ValueError(core)


def vectors(model, frac_min=FRAC_MIN, core="conserving", floored="exclude"):
    """Reference magnitudes m[(band, epoch)], per-leg errors, the mask and
    per-epoch info (T_eff, floored flag) for the nuisance columns."""
    epochs = model["epochs"]
    keys = [(b, t) for t in epochs for b in BANDS]
    m_ref = {k: np.nan for k in keys}
    d_rt = {leg: {k: np.nan for k in keys} for leg in LEGS}
    mask = {k: False for k in keys}
    info = {"core": core, "floored": floored, "T_eff": {}, "floored_rows": {}, "point":
            (model["m_ej_msun"], model["v_ej_c"], model["x_lan"])}
    for r in model["rows"]:
        if r.get("status") not in ("ok", "reduced_n"):
            continue
        t = r["t_d"]
        fl = row_floored(r, model)
        info["T_eff"][t] = r["source"]["T_eff"]; info["floored_rows"][t] = fl
        frac = band_fraction(r, model["lam_window"], model["n_spec"])
        s_ref = core_shift(r["ref"], core)
        for b in BANDS:
            mag = r["ref"]["mags"].get(b, np.nan) + s_ref
            m_ref[(b, t)] = mag
            mask[(b, t)] = bool(np.isfinite(mag) and frac[b] >= frac_min and mag <= MAG_LIMIT[b]
                                and not (fl and floored == "exclude"))
            for leg in LEGS:
                o = r["legs"][leg]
                d_rt[leg][(b, t)] = o["dm"].get(b, np.nan) + core_shift(o, core) - s_ref
    return keys, m_ref, d_rt, mask, info


# --- nuisance columns -------------------------------------------------------------

_TBB = {}


def tbb_derivative(T, eps=0.05):
    """dm_b/dlnT of a Planck spectrum at fixed L (R^2 ~ T^-4) through the real
    passbands on the grid's 200-bin window, central difference at T(1 +- eps)."""
    key = (round(float(T), 3), eps)
    if key not in _TBB:
        from sobolev.formal_transfer import planck_bnu
        edges = phot.nu_edges(*LAM_WIN, N_SPEC)
        nu_c = np.sqrt(edges[1:] * edges[:-1])
        pb = phot.load_passbands()

        def mags(TT):
            lnu = 4 * np.pi ** 2 * (1e15 / TT ** 2) ** 2 * planck_bnu(nu_c, TT)   # R^2 = const / T^4
            return phot.magnitudes(nu_c, lnu, pb, phot.D_40MPC, edges)
        hi, lo = mags(T * (1 + eps)), mags(T * (1 - eps))
        dl = np.log((1 + eps) / (1 - eps))
        _TBB[key] = {b: (hi[b] - lo[b]) / dl for b in hi}
    return _TBB[key]


def nuisance_columns(vecs, point, ks, nuisance, override=None):
    """Extra design columns for the keys `ks`; returns (names, columns, notes).
    `override` = {name: {(band, t): value}} replaces a column by a measured
    direction (the MC T_eff direction of §4.43 in place of the Planck proxy)."""
    keys, m0, _, _, info = vecs[point]
    names, cols, notes = [], [], {}
    override = override or {}
    if "L_t" in nuisance:
        for t in sorted({k[1] for k in ks}):
            names.append(f"L_{t:g}"); cols.append([1.0 if k[1] == t else 0.0 for k in ks])
    if "T_bb" in nuisance:
        names.append("T_bb")
        if "T_bb" in override:
            cols.append([override["T_bb"][k] for k in ks]); notes["T_bb"] = "override"
        else:
            cols.append([tbb_derivative(info["T_eff"][k[1]])[k[0]] for k in ks])
    if "2c" in nuisance:
        blue = (point[0], point[1], X_GRID[0])
        if point[2] == X_GRID[0] or blue not in vecs:
            notes["2c"] = "skipped: the point is the blue component itself" if point[2] == X_GRID[0] \
                else "skipped: no X = 1e-3 model"
        else:
            m_blue = vecs[blue][1]
            ratio = np.array([10 ** (-0.4 * (m_blue[k] - m0[k])) for k in ks])
            names.append("2c"); cols.append((-(2.5 / np.log(10)) * ratio).tolist())
            notes["2c_flux_ratio"] = ratio.tolist()
    return names, cols, notes


def restrict_keys(vecs, point, ks, nuisance):
    """Keys a tangent space can use: under L_t an epoch with a single live
    observable carries no information (drop it); under 2c the blue model's
    magnitude must exist."""
    if "2c" in nuisance:
        blue = (point[0], point[1], X_GRID[0])
        if point[2] != X_GRID[0] and blue in vecs:
            m_blue = vecs[blue][1]
            ks = [k for k in ks if np.isfinite(m_blue[k])]
    if "L_t" in nuisance:
        count = {}
        for k in ks:
            count[k[1]] = count.get(k[1], 0) + 1
        ks = [k for k in ks if count[k[1]] >= 2]
    return ks


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
    keys, m0, _, mask0 = vecs[point][:4]
    if scheme == "central" and lo in have and hi in have:
        _, m_lo, _, mk_lo = vecs[lo][:4]; _, m_hi, _, mk_hi = vecs[hi][:4]
        dl = np.log(hi[i] / lo[i])
        d = {k: (m_hi[k] - m_lo[k]) / dl for k in keys}
        mask = {k: mask0[k] and mk_lo[k] and mk_hi[k] for k in keys}
        return d, mask, {"scheme": "central", "dln": dl, "one_sided": False}
    # one-sided (edge, or requested): prefer the upper neighbour
    for nb, sign in ((hi, 1.0), (lo, -1.0)):
        if nb in have:
            _, m_nb, _, mk_nb = vecs[nb][:4]
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
    keys, m0, _, mask0 = vecs[point][:4]
    _, m_lo, _, mk_lo = vecs[lo][:4]; _, m_hi, _, mk_hi = vecs[hi][:4]
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
    signif[~np.any(D != 0, axis=0)] = 0.0          # a dead (all-zero) column has no amplitude
    rank = int(np.linalg.matrix_rank(D * np.sqrt(W)[:, None])) if np.all(np.isfinite(D)) else 0
    dof = N - rank
    chi2_res = float(np.sum(W * res ** 2))
    return {"a": a.tolist(), "cov": cov.tolist(), "signif": signif.tolist(), "R": R,
            "chi2_res_N": chi2_res / N, "chi2_RT_N": float(np.sum(W * d_rt ** 2) / N),
            "chi2_res_dof": chi2_res / dof if dof > 0 else np.nan, "rank": rank, "dof": int(dof),
            "cos": cos, "cond": cond, "N": int(N)}


def classify(p, thresh=THRESH, n_phys=None):
    """Decided on the physical columns' significance only (`n_phys`); a
    significant nuisance amplitude never makes C-A."""
    if not np.isfinite(p["chi2_RT_N"]):
        return "undefined"
    if p.get("dof", DOF_MIN) < DOF_MIN:
        return "underdetermined"
    if p["chi2_RT_N"] <= thresh["chi2_small"]:
        return "C-C"
    signif = p["signif"][:n_phys] if n_phys else p["signif"]
    if p["R"] <= thresh["R_max"] and max(signif) >= thresh["signif_min"]:
        return "C-A"
    return "C-B"


def analyse_point(vecs, point, leg, sigma_scale=1.0, scheme="central", frac_min=None,
                  noise_floor=None, sigma=None, obs_mask=None, nuisance=(), column_override=None):
    """WLS projection of d_RT on [d_M, d_v, d_X | nuisance columns].

    `sigma`: {(band, t): sigma} overriding SIGMA; `obs_mask`: the (band, t)
    an observing scenario provides, intersected with the physical mask;
    `nuisance`: a subset of NUISANCE (see `TANGENT`).
    """
    keys, m0, d_rt_all, mask0 = vecs[point][:4]
    ders, masks, infos = [], [], []
    for th in THETA:
        d, mk, info = derivative(vecs, point, th, scheme)
        if d is None:
            return {"status": f"no derivative along {th}"}
        ders.append(d); masks.append(mk); infos.append(info)
    ks = [k for k in keys if mask0[k] and all(mk[k] for mk in masks)
          and np.isfinite(d_rt_all[leg][k]) and all(np.isfinite(d[k]) for d in ders)
          and (obs_mask is None or k in obs_mask)
          and all(np.isfinite(c.get(k, np.nan)) for c in (column_override or {}).values())]
    ks = restrict_keys(vecs, point, ks, nuisance)
    if len(ks) < 4:
        return {"status": f"only {len(ks)} usable observables", "N": len(ks)}
    d_rt = np.array([d_rt_all[leg][k] for k in ks])
    D = np.array([[d[k] for d in ders] for k in ks])
    if noise_floor is not None:
        # zero derivative components below 3x the amplified MC floor
        for j, info in enumerate(infos):
            amp = np.sqrt(2.0) / abs(info["dln"])
            D[np.abs(D[:, j]) < 3.0 * noise_floor * amp, j] = 0.0
    names, cols, notes = nuisance_columns(vecs, point, ks, nuisance, column_override)
    columns = list(THETA) + names
    if cols:
        D = np.hstack([D, np.array(cols).T])
    sig = np.array([(sigma[k] if sigma is not None else SIGMA[k[0]]) for k in ks]) * sigma_scale
    p = project(d_rt, D, sig)
    p.update(status="ok", keys=[list(k) for k in ks], derivative_info=infos,
             d_rt_norm=float(np.sqrt(np.mean(d_rt ** 2))),
             d_rt=d_rt.tolist(), fit=(D @ np.array(p["a"])).tolist(), sigma=sig.tolist(),
             D=D.tolist(), columns=columns, n_phys=N_PHYS, p=len(columns),
             nuisance=list(nuisance), nuisance_used=names, nuisance_notes=notes,
             a_nuisance={n: p["a"][N_PHYS + j] for j, n in enumerate(names)},
             signif_nuisance={n: p["signif"][N_PHYS + j] for j, n in enumerate(names)},
             low_N=len(ks) < N_LOW,
             live_params=[bool(np.any(D[:, j] != 0)) for j in range(D.shape[1])])
    p["cls"] = classify(p, n_phys=N_PHYS)
    p["underdetermined"] = p["cls"] == "underdetermined"
    # is the fitted shift inside the stencil the derivative was taken on? |a| > dln
    # means the projection extrapolates beyond the neighbouring grid model
    p["a_over_dln"] = [abs(p["a"][j]) / abs(info["dln"]) for j, info in enumerate(infos)]
    p["extrapolated"] = bool(max(p["a_over_dln"]) > 1.0)
    if names:
        # what the nuisance columns absorb on their own, and what the three
        # physical parameters add beyond them
        q = project(d_rt, D[:, N_PHYS:], sig)
        p["R_nuisance_only"] = q["R"]; p["chi2_res_dof_nuisance_only"] = q["chi2_res_dof"]
        p["a_nuisance_only"] = {n: q["a"][j] for j, n in enumerate(names)}
    if "2c" in names:
        j = columns.index("2c")
        p["lin_2c"] = float(np.max(np.abs(p["a"][j] * np.array(notes["2c_flux_ratio"]))))
    return p


def noise_floor_of(vecs, point):
    keys, _, d_rt_all, mask = vecs[point][:4]
    v = [abs(d_rt_all["A_redist"][k]) for k in keys if mask[k] and np.isfinite(d_rt_all["A_redist"][k])]
    return float(np.max(v)) if v else 0.0


def apply_override(models, chain_file, chain_max):
    """Substitute one cell's reference and legs by a `robustness.py chain` run."""
    rec = json.loads(Path(chain_file).read_text())
    run = rec["runs"][str(chain_max)]
    point = (rec["m_ej_msun"], rec["v_ej_c"], rec["x_lan"])
    row = next(r for r in models[point]["rows"] if r["t_d"] == rec["t_d"])
    row["ref"].update({k: run["ref"][k] for k in ("mags", "colors", "L_bol", "L_bol_absorbing",
                                                    "n_trapped", "f_return", "f_dep")})
    for leg, o in run["legs"].items():
        row["legs"][leg].update({k: o[k] for k in ("mags", "colors", "dm", "dcolor", "L_bol",
                                                     "L_bol_absorbing", "dm_bol_absorbing")})
    row["override"] = {"file": str(chain_file), "chain_max": chain_max}
    return point, rec["t_d"]


def summarize(results, nuisance=()):
    summary = {}
    for leg in LEGS:
        ok = [r["legs"][leg] for r in results.values() if r["legs"][leg].get("status") == "ok"]
        cls = [q["cls"] for q in ok]
        s = {c: cls.count(c) for c in ("C-A", "C-B", "C-C", "underdetermined")}
        s["n"] = len(cls)
        s["unstable"] = sum(q.get("unstable", False) for q in ok)
        s["low_N"] = sum(q["low_N"] for q in ok)
        s["well_sampled"] = {c: sum(q["cls"] == c and not q["low_N"] for q in ok)
                             for c in ("C-A", "C-B", "C-C", "underdetermined")}
        det = [q for q in ok if not q["underdetermined"]]
        s["median_chi2_RT_N"] = float(np.median([q["chi2_RT_N"] for q in det])) if det else None
        s["median_R"] = float(np.median([q["R"] for q in det])) if det else None
        s["median_chi2_res_dof"] = float(np.median([q["chi2_res_dof"] for q in det])) if det else None
        s["median_dof"] = float(np.median([q["dof"] for q in ok])) if ok else None
        # the two quantities behind the class: how much of the norm the expanded
        # model absorbs, and whether what is left is still detectable at SIGMA
        s["absorbed_R_le_0.3"] = sum(q["R"] <= THRESH["R_max"] for q in det)
        s["leftover_chi2_res_dof_gt_4"] = sum(q["chi2_res_dof"] > THRESH["chi2_small"] for q in det)
        s["absorbed_and_leftover"] = sum(q["R"] <= THRESH["R_max"] and q["chi2_res_dof"] > THRESH["chi2_small"]
                                         for q in det)
        s["n_determined"] = len(det)
        s["C-A_extrapolated"] = sum(q["cls"] == "C-A" and q["extrapolated"] for q in det)
        s["max_a_over_dln"] = {"median": float(np.median([max(q["a_over_dln"]) for q in det])),
                               "max": float(np.max([max(q["a_over_dln"]) for q in det]))} if det else None
        if "T_bb" in nuisance:
            v = [q["a_nuisance"]["T_bb"] for q in det if "T_bb" in q["a_nuisance"]]
            s["a_T_bb"] = {"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v))} if v else None
        if "2c" in nuisance:
            v = [q["a_nuisance"]["2c"] for q in det if "2c" in q["a_nuisance"]]
            s["a_2c"] = {"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v)),
                         "max_lin": float(np.max([q["lin_2c"] for q in det if "lin_2c" in q]))} if v else None
        if "L_t" in nuisance:
            v = [abs(a) for q in det for n, a in q["a_nuisance"].items() if n.startswith("L_")]
            s["abs_a_L"] = {"median": float(np.median(v)), "max": float(np.max(v))} if v else None
        summary[leg] = s
    return summary


def main(grid_dir, out=None, core="conserving", floored="exclude", tangent="T0",
         override=None, override_chain=None):
    models = load_grid(grid_dir)
    print(f"{len(models)} of 27 models present")
    if override:
        for f in override:
            print("override:", apply_override(models, f, override_chain), f)
    nuisance = TANGENT[tangent]
    vecs1 = {p: vectors(m, FRAC_MIN, core, floored) for p, m in models.items()}
    vecs3 = {p: vectors(m, 0.03, core, floored) for p, m in models.items()}
    results = {}
    for point in vecs1:
        floor = noise_floor_of(vecs1, point)
        rec = {"noise_floor": floor, "secants": {th: secant_disagreement(vecs1, point, th) for th in THETA},
               "floored_rows": [t for t, f in vecs1[point][4]["floored_rows"].items() if f], "legs": {}}
        kw = {"noise_floor": floor, "nuisance": nuisance}
        for leg in LEGS:
            base = analyse_point(vecs1, point, leg, **kw)
            if base["status"] != "ok":
                rec["legs"][leg] = base; continue
            alt = {"sigma_x2": analyse_point(vecs1, point, leg, sigma_scale=2.0, **kw),
                   "one_sided": analyse_point(vecs1, point, leg, scheme="one_sided", **kw),
                   "frac_3pct": analyse_point(vecs3, point, leg, **kw)}
            base["robustness"] = {k: v.get("cls", v["status"]) for k, v in alt.items()}
            base["unstable"] = any(v.get("cls") not in (None, base["cls"]) for v in alt.values())
            rec["legs"][leg] = base
        results[str(point)] = rec
    summary = summarize(results, nuisance)
    d = {"thresholds": THRESH, "n_low": N_LOW, "dof_min": DOF_MIN, "sigma": SIGMA, "mag_limit": MAG_LIMIT,
         "frac_min": FRAC_MIN, "core": core, "floored": floored, "tangent": tangent,
         "nuisance": list(nuisance), "override": override, "override_chain": override_chain,
         "points": results, "summary": summary}
    if out is None:
        tags = ([core] if core != "conserving" else []) + (["floored_incl"] if floored == "include" else []) \
            + ([tangent] if tangent != "T0" else [])
        out = HERE / ("sensitivity" + "".join(f"_{t}" for t in tags) + ".json")
    Path(out).write_text(json.dumps(d, indent=1))
    print_table(d)
    print(f"wrote {out}")
    return d


def print_table(d):
    print(f"\n{'point':22s}{'leg':11s}{'N':>4}{'dof':>4}{'chi2RT/N':>10}{'R':>7}{'chi2res/dof':>12}"
          f"{'a_M':>7}{'a_v':>7}{'a_X':>7}{'cosM':>6}{'cosv':>6}{'cosX':>6}  class  nuisance")
    for pt, r in d["points"].items():
        for leg, p in r["legs"].items():
            if p.get("status") != "ok":
                print(f"{pt:22s}{leg:11s}  {p['status']}"); continue
            nuis = " ".join(f"{n}={a:+.2f}" for n, a in p.get("a_nuisance", {}).items()
                            if not n.startswith("L_"))
            aL = [abs(a) for n, a in p.get("a_nuisance", {}).items() if n.startswith("L_")]
            if aL:
                nuis += f" |a_L|max={max(aL):.2f}"
            print(f"{pt:22s}{leg:11s}{p['N']:4d}{p['dof']:4d}{p['chi2_RT_N']:10.1f}{p['R']:7.2f}"
                  f"{p['chi2_res_dof']:12.1f}"
                  + "".join(f"{a:7.2f}" for a in p["a"][:N_PHYS]) + "".join(f"{c:6.2f}" for c in p["cos"][:N_PHYS])
                  + f"  {p['cls']}{' *' if p.get('unstable') else ''}{' lowN' if p.get('low_N') else ''}  {nuis}")
    print("\nsummary:", json.dumps(d["summary"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid-dir", default=str(HERE / "grid"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--core", default="conserving", choices=("conserving", "absorbing"))
    ap.add_argument("--floored", default="exclude", choices=("exclude", "include"))
    ap.add_argument("--tangent", default="T0", choices=tuple(TANGENT))
    ap.add_argument("--override", nargs="*", default=None, help="robustness chain JSON(s)")
    ap.add_argument("--override-chain", type=int, default=8000)
    a = ap.parse_args()
    main(a.grid_dir, a.out, a.core, a.floored, a.tangent, a.override, a.override_chain)
