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
       python robustness.py table
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


def table(files=None):
    files = files or sorted(OUT.glob("chain_model_*.json"))
    files = [f for f in files if not f.stem.endswith("_probe")]
    print("| cell | chain_max | trapped (ref) | ref t (s) | leg | Δ(g−r) | Δ(i−J) | Δ(J−K) | max Δcol change vs 2000 (% of residual) |")
    print("|---|---|---|---|---|---|---|---|---|")
    summary = []
    for f in files:
        rec = json.loads(f.read_text())
        n3 = 3 * rec["n_used"]
        cell = f"({rec['m_ej_msun']:g}, {rec['v_ej_c']:g}, {rec['x_lan']:g}) @ {rec['t_d']:g} d"
        base = rec["runs"].get("2000")
        for c, run in sorted(rec["runs"].items(), key=lambda kv: int(kv[0])):
            if "ref" not in run:
                print(f"| {cell} | {c} | — | — | — | {run['status']} | | | |"); continue
            for leg in LEGS:
                dc = run["legs"][leg]["dcolor"]
                chg = ""
                if base and "ref" in base and c != "2000":
                    rel = [abs(dc[k] - base["legs"][leg]["dcolor"][k]) / max(abs(base["legs"][leg]["dcolor"][k]), 1e-9)
                           for k in COLS]
                    chg = f"{100*max(rel):.0f}"
                    summary.append((cell, int(c), leg, max(rel),
                                    all(np.sign(dc[k]) == np.sign(base["legs"][leg]["dcolor"][k]) for k in COLS)))
                print(f"| {cell} | {c} | {run['ref']['n_trapped']/n3*100:.1f} % | "
                      f"{run['timing']['reference']:.0f} | {leg} | "
                      + " | ".join(f"{dc[k]:+.3f}" for k in COLS) + f" | {chg} |")
    if summary:
        worst = max(summary, key=lambda s: s[3])
        print(f"\nworst change of a colour error vs chain 2000: {100*worst[3]:.0f} % "
              f"({worst[0]}, chain {worst[1]}, {worst[2]}); signs preserved in "
              f"{sum(s[4] for s in summary)}/{len(summary)} (cell, cap, leg) triples")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("chain", "table"))
    ap.add_argument("--model", default=None)
    ap.add_argument("--t", type=float, default=None)
    ap.add_argument("--chain-max", type=int, nargs="+", default=[2000, 4000, 8000])
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--wall-budget", type=float, default=3600.0)
    a = ap.parse_args()
    if a.cmd == "chain":
        chain(a.model, a.t, a.chain_max, probe=a.probe, wall_budget_s=a.wall_budget)
    else:
        table()
