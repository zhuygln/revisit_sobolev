"""Print the light-curve residual tables from a merged summary.json."""
import json, sys
import numpy as np
s = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "prod/merged/summary.json"))
t = np.array(s["t_grid"]) / 86400.0; tm = np.sqrt(t[1:] * t[:-1])
for leg, o in s["legs"].items():
    if "dm_vs_ref" not in o:
        continue
    dm = o["dm_vs_ref"]; dc = o["dcolour_vs_ref"]
    bands = list(dm[0].keys())
    mx = {b: np.nanmax(np.abs([d[b] for d in dm])) if np.isfinite([d[b] for d in dm]).any() else np.nan for b in bands}
    i_pk = int(np.nanargmax(np.where(np.isfinite(o["L_esc"]), o["L_esc"], -np.inf)))
    at_pk = {b: dm[i_pk][b] for b in bands}
    col_pk = {k: dc[i_pk][k] for k in dc[i_pk]}
    print(f"{leg:6s} dL_peak {o['dL_peak']:+.3f} dt_peak {o['dt_peak']:+.3f} dE_rad {o['dE_rad']:+.3f} max|dcolour| {o['max_abs_dcolour']:.2f} NaN cells {o['nan_band_cells']}")
    print("       max|dm| per band: " + " ".join(f"{b} {mx[b]:.2f}" for b in bands))
    print("       dm at the peak:   " + " ".join(f"{b} {at_pk[b]:+.2f}" for b in bands))
    print("       dcolour at peak:  " + " ".join(f"{k} {v:+.2f}" for k, v in col_pk.items()))
    # the slab of the maximum colour residual
    allc = [(abs(v), k, j) for j, d in enumerate(dc) for k, v in d.items() if np.isfinite(v)]
    if allc:
        m = max(allc); print(f"       max residual {m[1]} at t = {tm[m[2]]:.1f} d")
print("readings:", {k: (v.get("outcome") if isinstance(v, dict) and "outcome" in v else v) for k, v in s["readings"].items()})
