"""Markdown tables from the grid JSONs (report section 4.40).

  --which cells    one row per (model, epoch): status, n_used, S, f_ret, f_dep,
                   the A_redist noise floor, C_both / C_binned worst |dcolour|
                   and the absorbing-core dm_bol the conserving normalization
                   removes
  --which models   one row per model: epochs run / over budget / failed,
                   worst colour error over the epochs that ran
  --which central  the central model's per-epoch colour table (fig 2's numbers)
  --which summary  the numbers the report quotes (NIR sign count, band extremes, floor)
  --which json     write everything (`summary`) to grid_table.json for freeze.py

Every amplitude is printed as value ± floor, the floor being that cell's worst
live |Δm| of A_redist (the MC-noise statement); the chain-cap systematic is a
separate statement (robustness.chain_summary).

Usage: python grid_table.py [--which cells|models|central|summary|json|all] [--grid-dir grid]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_grid import M_GRID, V_GRID, X_GRID   # noqa: E402
from grid import model_name                    # noqa: E402
import sensitivity as sens                     # noqa: E402  (mask rule shared with Gate 2)

COLS = ("g-r", "r-i", "i-z", "i-J", "J-K")
RAN = ("ok", "reduced_n")
N_WELL = 100_000   # the A_redist floor is a MC-noise statement; quote it from well-sampled cells only


def load(grid_dir):
    out = []
    for m in M_GRID:
        for v in V_GRID:
            for x in X_GRID:
                p = Path(grid_dir) / f"{model_name(m, v, x)}.json"
                if p.exists():
                    out.append(((m, v, x), json.loads(p.read_text())))
    return out


def worst(d):
    vals = [abs(v) for v in d.values() if v is not None and np.isfinite(v)]
    if not vals:
        return np.nan, ""
    k = max(d, key=lambda k: abs(d[k]) if d[k] is not None and np.isfinite(d[k]) else -1)
    return abs(d[k]), k


def live_bands(r, d):
    """The sensitivity mask for one row: band >= FRAC_MIN of L_bol, brighter than
    the magnitude limit, photosphere not floored (same rule as sensitivity.py)."""
    frac = sens.band_fraction(r, d["lam_window"], d["n_spec"])
    return {b for b in sens.BANDS
            if np.isfinite(r["ref"]["mags"].get(b, np.nan)) and frac[b] >= sens.FRAC_MIN
            and r["ref"]["mags"][b] <= sens.MAG_LIMIT[b] and not sens.row_floored(r, d)}


def masked(d, live):
    """Keep dm entries whose band is live, and dcolor entries whose two bands are live."""
    return {k: v for k, v in d.items() if all(b in live for b in k.split("-"))}


NIR_COLS = ("i-J", "J-K")


def trapped_fraction(r):
    """Reference packets thermalized by the chain cap, as a fraction of the 3-seed total."""
    return r["ref"]["n_trapped"] / (3 * r["n_used"]) if r.get("n_used") else float("nan")


def rows(models):
    """One record per (model, epoch): the cell table as data (report §4.40/§4.44,
    Paper III ED Table 1). Colour and magnitude errors are the live (masked) ones;
    `floor` is the cell's A_redist worst live |Δm| -- the MC-noise floor every
    quoted amplitude is printed against."""
    out = []
    for (m, v, x), d in models:
        for r in d["rows"]:
            ran = r["status"] in RAN
            rec = {"point": [m, v, x], "t_d": r["t_d"], "status": r["status"], "ran": ran,
                   "floored": bool(ran and sens.row_floored(r, d)), "redo_budget_s": r.get("redo", {}).get("budget_s"),
                   "n_used": r.get("n_used"), "S": r.get("band_S_band"), "T_gas": r.get("T_gas")}
            if not ran:
                rec.update(projected_s=r.get("projected_s"), error=r.get("error"),
                           f_return=r.get("probe_f_return"), f_dep=r.get("probe_f_dep"))
                out.append(rec); continue
            lv = live_bands(r, d)
            rec.update(f_return=r["ref"]["f_return"], f_dep=r["ref"]["f_dep"], trapped_frac=trapped_fraction(r),
                       live=sorted(lv, key=sens.BANDS.index), floor=worst(masked(r["legs"]["A_redist"]["dm"], lv))[0],
                       dm_bol_absorbing=r["legs"]["C_both"]["dm_bol_absorbing"], legs={})
            for leg in ("A_redist", "B_opacity", "C_both", "C_binned"):
                l = r["legs"][leg]
                w, k = worst(masked(l["dcolor"], lv))
                rec["legs"][leg] = {"dm": masked(l["dm"], lv), "dcolor": masked(l["dcolor"], lv),
                                    "worst_dcolor": w, "worst_key": k}
            out.append(rec)
    return out


def per_point(cells):
    """Aggregate `rows` per grid point (the per-model table)."""
    pts = {}
    for c in cells:
        pts.setdefault(tuple(c["point"]), []).append(c)
    out = []
    for pt, cs in pts.items():
        ran = [c for c in cs if c["ran"]]
        rec = {"point": list(pt), "n_epochs": len(cs), "ran": len(ran),
               "over_budget": sum(c["status"] == "over_budget" for c in cs),
               "failed": sum(not c["ran"] and c["status"] != "over_budget" for c in cs),
               "epochs_run": [c["t_d"] for c in ran], "floored_epochs": [c["t_d"] for c in ran if c["floored"]],
               "redone_epochs": [c["t_d"] for c in ran if c["redo_budget_s"]],
               "n_used_min": min((c["n_used"] for c in ran), default=None),
               "n_used_max": max((c["n_used"] for c in ran), default=None),
               "n_live": sum(len(c["live"]) for c in ran),
               "trapped_frac_max": max((c["trapped_frac"] for c in ran), default=None),
               # the A_redist floor is quoted from well-sampled cells only (n_used >= N_WELL)
               "floor_max_well": max((c["floor"] for c in ran if c["n_used"] >= N_WELL and np.isfinite(c["floor"])),
                                     default=None),
               "floor_max_all": max((c["floor"] for c in ran if np.isfinite(c["floor"])), default=None)}
        for leg in ("C_both", "C_binned"):
            best = max((c for c in ran if np.isfinite(c["legs"][leg]["worst_dcolor"])),
                       key=lambda c: c["legs"][leg]["worst_dcolor"], default=None)
            rec[leg] = ({"worst_dcolor": best["legs"][leg]["worst_dcolor"], "worst_key": best["legs"][leg]["worst_key"],
                         "t_d": best["t_d"], "floor_at_cell": best["floor"]} if best else None)
        out.append(rec)
    return out


def nir_sign_count(cells, leg="C_both", cols=NIR_COLS):
    """(negative, total) over the live NIR colour errors of `leg` -- the F48
    '195 of 199' statement (negative = the closure leg is too blue in the NIR)."""
    vals = [c["legs"][leg]["dcolor"][k] for c in cells if c["ran"] for k in cols if k in c["legs"][leg]["dcolor"]]
    return sum(v < 0 for v in vals), len(vals)


def band_extremes(cells, leg="C_both"):
    """Per band: n live, min, max, median of the live per-band closure error of `leg`."""
    out = {}
    for b in sens.BANDS:
        v = [c["legs"][leg]["dm"][b] for c in cells if c["ran"] and b in c["legs"][leg]["dm"]]
        out[b] = {"n": len(v), "min": min(v), "max": max(v), "median": float(np.median(v))} if v else None
    return out


def floor_stats(cells, n_min=N_WELL):
    """The A_redist noise floor: median and 90th percentile over well-sampled cells
    (n_used >= n_min) with a live band, and the range over the redone cells."""
    well = [c["floor"] for c in cells if c["ran"] and c["n_used"] >= n_min and np.isfinite(c["floor"])]
    redo = [c["floor"] for c in cells if c["ran"] and c["redo_budget_s"] and np.isfinite(c["floor"])]
    return {"n_min": n_min, "n_well": len(well), "median": float(np.median(well)), "p90": float(np.percentile(well, 90)),
            "max": float(np.max(well)), "n_redone": len(redo),
            "redone_range": [float(min(redo)), float(max(redo))] if redo else None}


def summary(models):
    """Everything freeze.py writes to grid_table.json."""
    cs = rows(models)
    neg, tot = nir_sign_count(cs)
    pp = per_point(cs)
    wc = [p["C_both"]["worst_dcolor"] for p in pp if p["C_both"]]
    return {"n_models": len(models), "n_well": N_WELL, "cells": cs, "points": pp,
            "nir_negative": [neg, tot], "band_extremes": {leg: band_extremes(cs, leg) for leg in ("C_both", "C_binned", "A_redist")},
            "floor": floor_stats(cs), "worst_dcolor_range": [float(min(wc)), float(max(wc))],
            "n_redone_cells": sum(bool(c["redo_budget_s"]) for c in cs)}


def _pm(v, f):
    return f"{v:.2f} ± {f:.2f}" if np.isfinite(v) else "—"


def cells(models):
    lines = ["| M | v | X | t (d) | status | n_used | S | f_ret | f_dep | trapped | floor (A) | C_both worst \\|Δcol\\| ± floor | C_binned worst \\|Δcol\\| ± floor | Δm_bol (absorbing) |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in rows(models):
        m, v, x = c["point"]
        s = c["status"] + (" (floored)" if c["floored"] else "") + (f" (redo {c['redo_budget_s']:g} s)" if c["redo_budget_s"] else "")
        base = f"| {m:g} | {v:g} | {x:g} | {c['t_d']:g} | {s} | {c['n_used'] or ''} | {c['S'] if c['S'] is not None else float('nan'):.1e} |"
        if c["ran"]:
            cb, cn = c["legs"]["C_both"], c["legs"]["C_binned"]
            lines.append(base + f" {c['f_return']:.2f} | {c['f_dep']:+.2f} | {100*c['trapped_frac']:.1f} % | {c['floor']:.3f} | "
                         f"{_pm(cb['worst_dcolor'], c['floor'])} ({cb['worst_key']}) | {_pm(cn['worst_dcolor'], c['floor'])} ({cn['worst_key']}) | "
                         f"{c['dm_bol_absorbing']:+.2f} |")
        else:
            extra = f"projected {c['projected_s']/3600:.1f} h" if c["status"] == "over_budget" else (c.get("error") or "")[:40]
            fr = c["f_return"] if c["f_return"] is not None else float("nan"); fd = c["f_dep"] if c["f_dep"] is not None else float("nan")
            lines.append(base + f" {fr:.2f} | {fd:+.2f} | | | {extra} | | |")
    return "\n".join(lines)


def per_model(models):
    lines = ["| M | v | X | ran | over budget | failed | epochs run | n_used (min–max) | floored epochs | redone | live (band, epoch) | max trapped | worst C_both \\|Δcol\\| ± floor | worst C_binned \\|Δcol\\| ± floor | max floor (A), n_used ≥ 1e5 |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for p in per_point(rows(models)):
        m, v, x = p["point"]
        ep = ", ".join(f"{t:g}" for t in p["epochs_run"])
        n_used = f"{p['n_used_min']}–{p['n_used_max']}" if p["n_used_min"] is not None else ""
        fl_ep = ", ".join(f"{t:g}" for t in p["floored_epochs"]) or "—"
        rd = ", ".join(f"{t:g}" for t in p["redone_epochs"]) or "—"
        w = []
        for leg in ("C_both", "C_binned"):
            q = p[leg]
            w.append(f"{_pm(q['worst_dcolor'], q['floor_at_cell'])} ({q['worst_key']}, {q['t_d']:g} d)" if q else "—")
        fl = f"{p['floor_max_well']:.3f}" if p["floor_max_well"] is not None else "—"
        tr = f"{100*p['trapped_frac_max']:.1f} %" if p["trapped_frac_max"] is not None else "—"
        lines.append(f"| {m:g} | {v:g} | {x:g} | {p['ran']} | {p['over_budget']} | {p['failed']} | {ep} | {n_used} | {fl_ep} | {rd} | "
                     f"{p['n_live']} | {tr} | {w[0]} | {w[1]} | {fl} |")
    return "\n".join(lines)


def central(models, key=(0.01, 0.1, 0.01)):
    d = dict(models).get(key)
    if d is None:
        return "(central model not present)"
    lines = ["| t (d) | T_eff | S | n_used | f_ret | f_dep | trapped | floor (A) | leg | Δ(g−r) | Δ(r−i) | Δ(i−z) | Δ(i−J) | Δ(J−K) | Δm_bol abs |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in d["rows"]:
        if r["status"] not in RAN:
            lines.append(f"| {r['t_d']:g} | {r['T_gas']:.0f} | {r.get('band_S_band', float('nan')):.1e} | — | | | | | {r['status']} | | | | | | |")
            continue
        lv = live_bands(r, d)
        fa, _ = worst(masked(r["legs"]["A_redist"]["dm"], lv))
        for leg in ("C_both", "C_binned"):
            l = r["legs"][leg]
            lines.append(f"| {r['t_d']:g} | {r['T_gas']:.0f} | {r['band_S_band']:.1e} | {r['n_used']} | {r['ref']['f_return']:.2f} | "
                         f"{r['ref']['f_dep']:+.2f} | {100*trapped_fraction(r):.1f} % | {fa:.3f} | {leg} | "
                         + " | ".join(f"{l['dcolor'][c]:+.2f}" + ("" if c.split("-")[0] in lv and c.split("-")[1] in lv else "†") for c in COLS)
                         + f" | {l['dm_bol_absorbing']:+.2f} |")
    lines.append("\n† colour with a masked band (below 1 % of L_bol, fainter than the limit, or a floored photosphere); not used in any claim.")
    return "\n".join(lines)


def summary_text(models):
    s = summary(models)
    be = s["band_extremes"]["C_both"]
    return (f"{s['n_models']} models, {len(s['cells'])} cells, {s['n_redone_cells']} redone; live NIR colour errors of C_both negative "
            f"at {s['nir_negative'][0]} of {s['nir_negative'][1]}; per-point worst C_both |Δcolour| "
            f"{s['worst_dcolor_range'][0]:.2f}–{s['worst_dcolor_range'][1]:.2f} mag; C_both Δm(g) {be['g']['min']:+.2f}…{be['g']['max']:+.2f} "
            f"(n = {be['g']['n']}), Δm(K) up to {be['K']['max']:+.2f} (n = {be['K']['n']}); A_redist floor on n_used ≥ {s['n_well']} "
            f"({s['floor']['n_well']} cells): median {s['floor']['median']:.3f}, 90 % {s['floor']['p90']:.3f}, max {s['floor']['max']:.3f}; "
            f"redone cells {s['floor']['redone_range'][0]:.2f}–{s['floor']['redone_range'][1]:.2f}.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="all", help="cells|models|central|summary|all, or json")
    ap.add_argument("--grid-dir", default=str(HERE / "grid"))
    ap.add_argument("--out", default=str(HERE / "grid_table.json"))
    a = ap.parse_args()
    models = load(a.grid_dir)
    if a.which == "json":
        Path(a.out).write_text(json.dumps(summary(models), indent=1))
        print("wrote", a.out); sys.exit()
    print(f"{len(models)} of 27 models present\n")
    if a.which in ("cells", "all"):
        print(cells(models), "\n")
    if a.which in ("models", "all"):
        print(per_model(models), "\n")
    if a.which in ("central", "all"):
        print(central(models), "\n")
    if a.which in ("summary", "all"):
        print(summary_text(models))
