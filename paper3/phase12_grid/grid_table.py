"""Markdown tables from the grid JSONs (report section 4.40).

  --which cells    one row per (model, epoch): status, n_used, S, f_ret, f_dep,
                   the A_redist noise floor, C_both / C_binned worst |dcolour|
                   and the absorbing-core dm_bol the conserving normalization
                   removes
  --which models   one row per model: epochs run / over budget / failed,
                   worst colour error over the epochs that ran
  --which central  the central model's per-epoch colour table (fig 2's numbers)

Usage: python grid_table.py [--which cells|models|central|all] [--grid-dir grid]
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
            and r["ref"]["mags"][b] <= sens.MAG_LIMIT[b] and not r["source"].get("v_ph_floored", False)}


def masked(d, live):
    """Keep dm entries whose band is live, and dcolor entries whose two bands are live."""
    return {k: v for k, v in d.items() if all(b in live for b in k.split("-"))}


def cells(models):
    lines = ["| M | v | X | t (d) | status | n_used | S | f_ret | f_dep | floor (A) | C_both worst \\|Δcol\\| | C_binned worst \\|Δcol\\| | Δm_bol (absorbing) |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (m, v, x), d in models:
        for r in d["rows"]:
            s = r["status"]
            base = f"| {m:g} | {v:g} | {x:g} | {r['t_d']:g} | {s} | {r.get('n_used', '')} | {r.get('band_S_band', float('nan')):.1e} |"
            if s in RAN:
                lv = live_bands(r, d)
                fa, _ = worst(masked(r["legs"]["A_redist"]["dm"], lv))
                wb, kb = worst(masked(r["legs"]["C_both"]["dcolor"], lv)); wn, kn = worst(masked(r["legs"]["C_binned"]["dcolor"], lv))
                lines.append(base + f" {r['ref']['f_return']:.2f} | {r['ref']['f_dep']:+.2f} | {fa:.3f} | "
                             f"{wb:.2f} ({kb}) | {wn:.2f} ({kn}) | {r['legs']['C_both']['dm_bol_absorbing']:+.2f} |")
            else:
                extra = f"projected {r['projected_s']/3600:.1f} h" if s == "over_budget" else r.get("error", "")[:40]
                lines.append(base + f" {r.get('probe_f_return', float('nan')):.2f} | {r.get('probe_f_dep', float('nan')):+.2f} | | {extra} | | |")
    return "\n".join(lines)


def per_model(models):
    lines = ["| M | v | X | ran | over budget | failed | epochs run | live (band, epoch) | worst C_both \\|Δcol\\| | worst C_binned \\|Δcol\\| | max floor (A) |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for (m, v, x), d in models:
        rows = d["rows"]
        ran = [r for r in rows if r["status"] in RAN]
        ob = sum(r["status"] == "over_budget" for r in rows)
        fail = sum(r["status"] not in RAN + ("over_budget",) for r in rows)
        lv = {r["t_d"]: live_bands(r, d) for r in ran}
        wb = max((worst(masked(r["legs"]["C_both"]["dcolor"], lv[r["t_d"]])) for r in ran), default=(np.nan, ""))
        wn = max((worst(masked(r["legs"]["C_binned"]["dcolor"], lv[r["t_d"]])) for r in ran), default=(np.nan, ""))
        fl = max((worst(masked(r["legs"]["A_redist"]["dm"], lv[r["t_d"]]))[0] for r in ran), default=np.nan)
        ep = ", ".join(f"{r['t_d']:g}" for r in ran)
        n_live = sum(len(v) for v in lv.values())
        lines.append(f"| {m:g} | {v:g} | {x:g} | {len(ran)} | {ob} | {fail} | {ep} | {n_live} | {wb[0]:.2f} ({wb[1]}) | {wn[0]:.2f} ({wn[1]}) | {fl:.3f} |")
    return "\n".join(lines)


def central(models, key=(0.01, 0.1, 0.01)):
    d = dict(models).get(key)
    if d is None:
        return "(central model not present)"
    lines = ["| t (d) | T_eff | S | n_used | f_ret | f_dep | floor (A) | leg | Δ(g−r) | Δ(r−i) | Δ(i−z) | Δ(i−J) | Δ(J−K) | Δm_bol abs |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in d["rows"]:
        if r["status"] not in RAN:
            lines.append(f"| {r['t_d']:g} | {r['T_gas']:.0f} | {r.get('band_S_band', float('nan')):.1e} | — | | | | {r['status']} | | | | | | |")
            continue
        lv = live_bands(r, d)
        fa, _ = worst(masked(r["legs"]["A_redist"]["dm"], lv))
        for leg in ("C_both", "C_binned"):
            l = r["legs"][leg]
            lines.append(f"| {r['t_d']:g} | {r['T_gas']:.0f} | {r['band_S_band']:.1e} | {r['n_used']} | {r['ref']['f_return']:.2f} | "
                         f"{r['ref']['f_dep']:+.2f} | {fa:.3f} | {leg} | "
                         + " | ".join(f"{l['dcolor'][c]:+.2f}" for c in COLS) + f" | {l['dm_bol_absorbing']:+.2f} |")
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="all")
    ap.add_argument("--grid-dir", default=str(HERE / "grid"))
    a = ap.parse_args()
    models = load(a.grid_dir)
    print(f"{len(models)} of 27 models present\n")
    if a.which in ("cells", "all"):
        print(cells(models), "\n")
    if a.which in ("models", "all"):
        print(per_model(models), "\n")
    if a.which in ("central", "all"):
        print(central(models))
