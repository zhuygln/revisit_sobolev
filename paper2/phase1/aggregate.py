"""Combine the production seeds of run_forest.py into means and scatters.

Reads forest_results_s{1,2,3}.json; writes forest_results.json with, per part
and leg, the seed mean and standard deviation of every band quantity and of
every differential, plus the SEDONA comparison rows for Part A.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
files = sorted(HERE.glob("forest_results_s[0-9]*.json"))
runs = [json.loads(f.read_text()) for f in files]
if not runs:
    sys.exit("no forest_results_s*.json yet")
out = {"n_seeds": len(runs), "seeds": [r["seed"] for r in runs], "n_packets": runs[0]["n_packets"],
       "geometry": runs[0]["geometry"]}
for part in ("part_A", "part_B"):
    legs = {}
    for tag in runs[0][part]["legs"]:
        vals = {k: np.array([r[part]["legs"][tag][k] for r in runs]) for k in ("band", "red", "bluewing", "forest")}
        dest = {k: np.array([r[part]["legs"][tag]["destination_of_band_launches"].get(k, np.nan) for r in runs])
                for k in ("escaped", "in_band", "redward", "blueward", "beyond_1um")}
        legs[tag] = {**{k: float(v.mean()) for k, v in vals.items()},
                     **{k + "_std": float(v.std(ddof=1)) if len(v) > 1 else None for k, v in vals.items()},
                     "destination": {k: float(np.nanmean(v)) for k, v in dest.items()},
                     "n_interactions": int(np.mean([r[part]["legs"][tag]["n_interactions"] for r in runs]))}
    diffs = {}
    for k in runs[0][part]["differentials"]:
        v = np.array([r[part]["differentials"][k] for r in runs])
        diffs[k] = dict(mean=float(v.mean()), std=float(v.std(ddof=1)) if len(v) > 1 else None)
    out[part] = dict(atom=runs[0][part]["atom"], legs=legs, differentials=diffs)
(HERE / "forest_results.json").write_text(json.dumps(out, indent=1))

print(f"{len(runs)} seeds x {out['n_packets']:.0e} packets")
for part in ("part_A", "part_B"):
    print(f"\n== {part} ==")
    for tag, v in out[part]["legs"].items():
        d = v["destination"]
        print(f"  {tag:28s} bluewing {v['bluewing']:.3f} forest {v['forest']:.3f} BAND {v['band']:.4f}+-{(v['band_std'] or 0):.4f} red {v['red']:.3f} | "
              f"esc {d['escaped']:.3f} in {d['in_band']:.2f} red {d['redward']:.2f} blue {d['blueward']:.2f} >1um {d['beyond_1um']:.2f}")
    for k, v in out[part]["differentials"].items():
        print(f"    {k:60s} {100*v['mean']:+7.2f}% +- {100*(v['std'] or 0):.2f}")
print("\nwrote forest_results.json")
