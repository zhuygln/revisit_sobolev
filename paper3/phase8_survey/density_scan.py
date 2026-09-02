"""Paper III: the Ce II / La II density scan behind F35 -- as a driver.

§4.32 reports the closure error changing sign along a density scan of Ce II
(binned: -33.4 -> -15.3 -> +21.4 -> +124.6 -> +94.7 -> +129.8 % at n_ion =
2910 x {1, 2, 3, 4, 6, 10}) and La II staying negative. Until 2026-09-02 those
numbers existed only in the lab notebook (§9af): no script produced them and no
JSON stored them. This is the driver. It reuses `survey.measure` at explicit
n_ion, so the geometry (Paper I), temperature (T_SHELL), kernel (N_g = 32) and
band (3800-3955 A) are exactly the survey's. The base density is the F36
global normalization (tau_max = 5 on the ion's strongest line), which for Ce II
is the third point of the scan.

The notebook numbers were a single seed; this runs the survey's three, so the
target is agreement to a few points, not bit equality.

Usage: python density_scan.py [--n 500000] [--ion 58CeII]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import survey as S

# multipliers of the global-normalization n_ion; Ce II's base is the scan's
# third point (n_ion = 8731 in §4.32, 8675 here), La II's is its second
SCAN = {"58CeII": (1 / 3, 2 / 3, 1.0, 4 / 3, 2.0, 10 / 3),
        "57LaII": (0.25, 1.0, 2.5)}


def crossing(rows, key="binned"):
    """S at which the signed band-3800 error crosses zero (either direction)."""
    r = sorted(rows, key=lambda x: x["band_S_band"])
    for a, b in zip(r, r[1:]):
        va, vb = a[key], b[key]
        if va * vb < 0:
            f = -va / (vb - va)
            return {"S": float(np.exp(np.log(a["band_S_band"]) +
                                      f * (np.log(b["band_S_band"]) - np.log(a["band_S_band"])))),
                    "direction": "neg_to_pos" if va < 0 else "pos_to_neg"}
    return None


def main(ion, n):
    lev, tr = S.DATA / f"{ion}_levels_calib.txt", S.DATA / f"{ion}_transitions_calib.txt"
    base, _ = S.n_ion_for(lev, tr)
    print(f"{ion}: global n_ion = {base:.0f}, {n} packets x {len(S.SEEDS)} seeds", flush=True)
    rows, t0 = [], time.time()
    for f in SCAN[ion]:
        r = S.measure(ion, n, n_ion=base * f)
        row = {"factor": f, "n_ion": r["n_ion"], "tau_max": r["tau_max"],
               "band_S_band": r["band_S_band"], "band_n_sat_band": r["band_n_sat_band"],
               "binned": r["legs"]["binned_group"]["band3800"],
               "expansion": r["legs"]["expansion_group"]["band3800"],
               "redist": r["legs"]["sobolev_group"]["band3800"],
               "opacity": r["legs"]["expansion_branch"]["band3800"],
               "wall_s": r["wall_s"]}
        rows.append(row)
        print(f"  n_ion {row['n_ion']:8.0f} tau_max {row['tau_max']:6.1f} S {row['band_S_band']:7.1f} | "
              f"binned {100*row['binned']:+7.1f}%  exp {100*row['expansion']:+7.1f}%  "
              f"[{row['wall_s']:.0f}s]", flush=True)
    out = {"ion": ion, "n": n, "ng": S.NG, "n_ion_global": float(base),
           "band": S.BAND, "rows": rows,
           "crossing": {k: crossing(rows, k) for k in ("binned", "expansion")}}
    (HERE / f"density_scan_{ion}.json").write_text(json.dumps(out, indent=1))
    print(f"crossing: {out['crossing']}\nwrote density_scan_{ion}.json  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="58CeII", choices=tuple(SCAN))
    ap.add_argument("--n", type=float, default=5e5)
    a = ap.parse_args()
    main(a.ion, int(a.n))
