"""Paper III, phase 12, §4.41: is the closure error an artefact of the chain cap?

The grid thermalizes a re-absorption chain that has not left its line after
`CHAIN_MAX = 2000` draws (`chain_overflow="absorb"`); the reference leg
thermalizes up to 14 % of its packets that way at the deepest X_lan = 0.1
cells, the closure legs none. `chain` re-runs one cell -- the reference AND
every leg, because the three group legs build their kernel from the
reference's event stream -- at its stored `n_used` and seeds for a set of chain
caps, through the same `grid.run_epoch` code path, and stores the per-leg
magnitudes and colour errors per cap. The 2000 run is the provenance check
(it must reproduce the committed row exactly); B_opacity does not depend on
the reference and must be identical at every cap.

Usage: python robustness.py chain --model model_M0.03_v0.05_X0.1 --t 3 --chain-max 2000 4000 8000
       python robustness.py chain ... --probe        (5000-packet reference timing only)
       python robustness.py table          (Markdown from chain_summary)
       python robustness.py summary        (chain_summary -> robustness/chain_table.json)
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grid                                                     # noqa: E402
from grid import DAY, SourceModel, git_sha                      # noqa: E402

OUT = HERE / "robustness"
CELLS = (("model_M0.03_v0.05_X0.1", 3.0), ("model_M0.03_v0.05_X0.1", 2.0),
         ("model_M0.01_v0.05_X0.1", 2.0), ("model_M0.03_v0.1_X0.1", 1.0))
COLS = ("g-r", "i-J", "J-K")
LEGS = ("C_both", "C_binned")


def load_row(model, t_d):
    d = json.loads((HERE / "grid" / f"{model}.json").read_text())
    row = next(r for r in d["rows"] if r["t_d"] == t_d)
    return d, row


def leg_summary(o):
    return {"mags": o["mags"], "colors": o["colors"], "dm": o.get("dm"), "dcolor": o.get("dcolor"),
            "n_trapped": o["n_trapped"], "f_return": o["f_return"], "f_dep": o["f_dep"],
            "L_bol": o["L_bol"], "L_bol_absorbing": o["L_bol_absorbing"],
            "dm_bol_absorbing": o.get("dm_bol_absorbing"), "t_wall": o.get("t_wall")}


def chain(model, t_d, chain_maxes, probe=False, wall_budget_s=3600.0, verbose=True):
    d, row = load_row(model, t_d)
    n_used = 5000 if probe else row["n_used"]
    src = SourceModel(d["m_ej_msun"], d["v_ej_c"], kappa=d["kappa_src"])
    st = src.state(t_d)
    t0 = time.time()
    atom = grid.build_state_atom(st, d["x_lan"])
    out = OUT / (f"chain_{model}_t{t_d:g}" + ("_probe" if probe else "") + ".json")
    rec = json.loads(out.read_text()) if out.exists() else {
        "model": model, "t_d": t_d, "m_ej_msun": d["m_ej_msun"], "v_ej_c": d["v_ej_c"],
        "x_lan": d["x_lan"], "n_used": n_used, "seeds": d["seeds"], "probe": probe,
        "stored": {"chain_max": d["chain_max"], "ref": leg_summary(row["ref"]),
                   "legs": {k: leg_summary(v) for k, v in row["legs"].items()},
                   "git": d.get("git")},
        "git": git_sha(), "runs": {}}
    if verbose:
        print(f"{model} t={t_d:g} d: n_used={n_used} atom {time.time()-t0:.0f}s  "
              f"stored n_trapped(ref)={row['ref']['n_trapped']} of {3*row['n_used']}", flush=True)
    for c in chain_maxes:
        if str(c) in rec["runs"]:
            continue
        r = grid.run_epoch(st, d["x_lan"], grid.N_DEFAULT, budget_s=wall_budget_s,
                           n_override=n_used, chain_max=c, atom=atom)
        run = {"status": r["status"], "t_wall": r["t_wall"], "timing": r.get("timing")}
        if r["status"] in ("ok", "reduced_n"):
            run.update(ref=leg_summary(r["ref"]), legs={k: leg_summary(v) for k, v in r["legs"].items()})
            if not probe:
                # provenance: same seeds and n through the same code path
                run["max_dmag_vs_stored"] = {
                    "ref": max(abs(r["ref"]["mags"][b] - row["ref"]["mags"][b]) for b in row["ref"]["mags"]),
                    **{k: max(abs(r["legs"][k]["mags"][b] - row["legs"][k]["mags"][b])
                              for b in row["ref"]["mags"]) for k in row["legs"]}}
        rec["runs"][str(c)] = run
        OUT.mkdir(exist_ok=True)
        out.write_text(json.dumps(rec, indent=1))
        if verbose:
            line = f"  chain_max={c:5d}  [{r['status']}, {r['t_wall']:.0f}s]"
            if "ref" in run:
                cb = r["legs"]["C_both"]
                line += (f"  ref {r['timing']['reference']:.0f}s  n_trapped={r['ref']['n_trapped']}  "
                         f"f_ret={r['ref']['f_return']:.2f}  d(g-r)={cb['dcolor']['g-r']:+.3f}  "
                         f"d(i-J)={cb['dcolor']['i-J']:+.3f}  d(J-K)={cb['dcolor']['J-K']:+.3f}")
                if "max_dmag_vs_stored" in run:
                    line += f"  |dm vs stored| ref {run['max_dmag_vs_stored']['ref']:.1e}" \
                            f" B {run['max_dmag_vs_stored']['B_opacity']:.1e}"
            print(line, flush=True)
    return rec


def _files(files=None):
    files = files or sorted(OUT.glob("chain_model_*.json"))
    return [Path(f) for f in files if not Path(f).stem.endswith("_probe")]


def chain_summary(files=None, base_cap="2000", rel_max=0.25):
    """The chain-cap test as data (§4.44.3): per cell and cap, the reference's
    trapped fraction, each leg's colour errors, the largest per-band change of the
    reference magnitude and of the C_both/C_binned closure error against the grid's
    cap, the relative change of each colour error and whether its sign is kept;
    plus the pre-declared criterion (each colour error changes by < `rel_max` and
    keeps its sign) tallied at the largest cap. Written to robustness/chain_table.json
    by freeze.py; `table` renders it."""
    cells = []
    for f in _files(files):
        rec = json.loads(f.read_text())
        n3 = 3 * rec["n_used"]
        base = rec["runs"].get(base_cap)
        bands = list(rec["stored"]["ref"]["mags"])
        cell = {"file": f.name, "model": rec["model"], "t_d": rec["t_d"],
                "point": [rec["m_ej_msun"], rec["v_ej_c"], rec["x_lan"]], "n_used": rec["n_used"],
                "floor_A_redist": max(abs(v) for v in rec["stored"]["legs"]["A_redist"]["dcolor"].values()),
                "stored_matches_base": None, "runs": {}}
        if base and "ref" in base:
            cell["stored_matches_base"] = max(abs(base["ref"]["mags"][b] - rec["stored"]["ref"]["mags"][b])
                                              for b in bands) < 1e-9
        for c, run in sorted(rec["runs"].items(), key=lambda kv: int(kv[0])):
            e = {"status": run["status"]}
            if "ref" in run:
                e.update(trapped_frac=run["ref"]["n_trapped"] / n3, t_ref_s=run["timing"]["reference"],
                         legs={leg: {"dcolor": {k: run["legs"][leg]["dcolor"][k] for k in COLS},
                                     "dm": dict(run["legs"][leg]["dm"])} for leg in LEGS})
                if base and "ref" in base and c != base_cap:
                    e["max_dm_ref_change"] = max(abs(run["ref"]["mags"][b] - base["ref"]["mags"][b]) for b in bands)
                    e["B_opacity_mags_identical"] = all(run["legs"]["B_opacity"]["mags"][b] == base["legs"]["B_opacity"]["mags"][b]
                                                        for b in bands)
                    for leg in LEGS:
                        dc, dc0 = run["legs"][leg]["dcolor"], base["legs"][leg]["dcolor"]
                        el = e["legs"][leg]
                        el["max_dm_change"] = max(abs(run["legs"][leg]["dm"][b] - base["legs"][leg]["dm"][b]) for b in bands)
                        el["rel"] = {k: abs(dc[k] - dc0[k]) / max(abs(dc0[k]), 1e-9) for k in COLS}
                        el["sign_kept"] = {k: bool(np.sign(dc[k]) == np.sign(dc0[k])) for k in COLS}
                        el["criterion_met"] = {k: bool(el["rel"][k] < rel_max and el["sign_kept"][k]) for k in COLS}
            cell["runs"][c] = e
        cells.append(cell)
    caps = sorted({c for cell in cells for c in cell["runs"]}, key=int)
    top = caps[-1] if caps else None
    summ = {"n_cells": len(cells), "base_cap": base_cap, "top_cap": top, "rel_max": rel_max,
            "trapped_range_base": _rng([cell["runs"][base_cap]["trapped_frac"] for cell in cells
                                        if "trapped_frac" in cell["runs"].get(base_cap, {})]),
            "trapped_range_top": _rng([cell["runs"][top]["trapped_frac"] for cell in cells
                                       if top and "trapped_frac" in cell["runs"].get(top, {})]),
            "floor_range": _rng([cell["floor_A_redist"] for cell in cells]),
            "stored_reproduced": all(cell["stored_matches_base"] for cell in cells)}
    for leg in LEGS:
        at_top = [cell["runs"][top]["legs"][leg] for cell in cells if top and "legs" in cell["runs"].get(top, {})]
        every = [cell["runs"][c]["legs"][leg] for cell in cells for c in cell["runs"]
                 if c != base_cap and "legs" in cell["runs"][c]]
        summ[leg] = {"dm_change_range_top": _rng([e["max_dm_change"] for e in at_top]),
                     "dm_change_range_all_caps": _rng([e["max_dm_change"] for e in every]),
                     "signs_kept_top": [sum(sum(e["sign_kept"].values()) for e in at_top), len(COLS) * len(at_top)],
                     "criterion_met_top": [sum(sum(e["criterion_met"].values()) for e in at_top), len(COLS) * len(at_top)],
                     "worst_rel_top": max((max(e["rel"].values()) for e in at_top), default=None)}
    summ["dm_ref_change_range_top"] = _rng([cell["runs"][top]["max_dm_ref_change"] for cell in cells
                                            if top and "max_dm_ref_change" in cell["runs"].get(top, {})])
    return {"cols": list(COLS), "legs": list(LEGS), "cells": cells, "summary": summ}


def _rng(v):
    v = [float(x) for x in v if x is not None]
    return [min(v), max(v)] if v else None


def table(files=None, data=None):
    """Markdown rendering of `chain_summary` (report §4.44.3)."""
    data = data or chain_summary(files)
    base_cap = data["summary"]["base_cap"]
    print("| cell | chain_max | trapped (ref) | ref t (s) | leg | Δ(g−r) | Δ(i−J) | Δ(J−K) | max \\|Δm_ref\\| vs base | max \\|Δ(Δm_b)\\| vs base | max Δcol change (% of residual) | floor (A) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for cell in data["cells"]:
        name = f"({cell['point'][0]:g}, {cell['point'][1]:g}, {cell['point'][2]:g}) @ {cell['t_d']:g} d"
        for c, run in cell["runs"].items():
            if "legs" not in run:
                print(f"| {name} | {c} | — | — | — | {run['status']} | | | | | | |"); continue
            for leg in data["legs"]:
                el = run["legs"][leg]
                chg = f"{100*max(el['rel'].values()):.0f}" if "rel" in el else "–"
                dref = f"{run['max_dm_ref_change']:.2f}" if "max_dm_ref_change" in run else "–"
                dleg = f"{el['max_dm_change']:.2f}" if "max_dm_change" in el else "–"
                print(f"| {name} | {c} | {100*run['trapped_frac']:.1f} % | {run['t_ref_s']:.0f} | {leg} | "
                      + " | ".join(f"{el['dcolor'][k]:+.3f}" for k in data["cols"])
                      + f" | {dref} | {dleg} | {chg} | {cell['floor_A_redist']:.2f} |")
    s = data["summary"]; top = s["top_cap"]
    for leg in data["legs"]:
        t = s[leg]
        print(f"\n{leg} at cap {top}: per-band change {t['dm_change_range_top'][0]:.2f}–{t['dm_change_range_top'][1]:.2f} mag, "
              f"signs kept {t['signs_kept_top'][0]}/{t['signs_kept_top'][1]}, < {100*s['rel_max']:.0f} % criterion met "
              f"{t['criterion_met_top'][0]}/{t['criterion_met_top'][1]}, worst change {100*t['worst_rel_top']:.0f} %")
    print(f"reference trapped {100*s['trapped_range_base'][0]:.1f}–{100*s['trapped_range_base'][1]:.1f} % at cap {base_cap}, "
          f"{100*s['trapped_range_top'][0]:.1f}–{100*s['trapped_range_top'][1]:.1f} % at cap {top}; "
          f"stored rows reproduced at cap {base_cap}: {s['stored_reproduced']}")
    return data


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("chain", "table", "summary"))
    ap.add_argument("--out", default=str(OUT / "chain_table.json"))
    ap.add_argument("--model", default=None)
    ap.add_argument("--t", type=float, default=None)
    ap.add_argument("--chain-max", type=int, nargs="+", default=[2000, 4000, 8000])
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--wall-budget", type=float, default=3600.0)
    a = ap.parse_args()
    if a.cmd == "chain":
        chain(a.model, a.t, a.chain_max, probe=a.probe, wall_budget_s=a.wall_budget)
    elif a.cmd == "table":
        table()
    else:
        Path(a.out).write_text(json.dumps(chain_summary(), indent=1))
        print("wrote", a.out)
