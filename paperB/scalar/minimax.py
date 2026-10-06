#!/usr/bin/env python3
"""The fair scalar comparator (the referee's point, the PI's 2026-10-06
decision): the operator is judged by max |dm| over the live bands AND max
|dcolour| at one threshold, while the preregistered eps* minimises the
MEAN band error. From the existing eps grids -- no new transport -- this
recomputes, per record, the preregistered eps* with its max band and
colour errors, the band-only minimax eps, and the joint minimax

    eps*_mm = argmin_eps max( max_b |dm_b|, max_c |dcolour_c| )

with its value. The preregistered eps* stays the gate definition (End
Matter); the Letter's main comparison quotes eps*_mm.

    .venv/bin/python paperB/scalar/minimax.py            # -> paperB/scalar/minimax.json; --markdown
"""
import importlib.util as _ilu
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
_spec = _ilu.spec_from_file_location("paperB_gate1_analyse", ROOT / "paperB/gate1/analyse.py")
A = _ilu.module_from_spec(_spec); _spec.loader.exec_module(A)
IONS = ("57LaII", "58CeII", "60NdII")
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II", "blend13": "13-ion blend"}


def g1_record(ion):
    """G1's highest-packet-count record (the one its analyse.py reads)."""
    cands = sorted((ROOT / "paperB/gate1").glob(f"gate1_{ion}*.json"), key=lambda p: json.loads(p.read_text())["n"])
    return (cands[-1].name, json.loads(cands[-1].read_text())) if cands else (None, None)


def curve(row):
    """[(eps, mean band, max band, max colour, joint)] over the record's eps legs, on the frozen live-band rule."""
    M = A.metrics(row)
    pts = []
    for tag, m in sorted(M["legs"].items(), key=lambda kv: kv[1].get("eps", -1)):
        if "eps" in m:
            pts.append(dict(eps=float(m["eps"]), mean_band=m["band"]["mean"], max_band=m["band"]["max"], max_colour=m["colour"]["max"],
                            joint=float(max(m["band"]["max"], m["colour"]["max"]))))
    return pts, M["live_bands"]


def read(row):
    pts, live = curve(row)
    i_pre = int(np.nanargmin([p["mean_band"] for p in pts]))          # the preregistered rule: mean |dm| over the live bands
    i_band = int(np.nanargmin([p["max_band"] for p in pts]))
    i_joint = int(np.nanargmin([p["joint"] for p in pts]))
    n = len(pts)
    return dict(live_bands=live, n_eps=n, curve=pts,
                prereg=dict(pts[i_pre], interior=0 < i_pre < n - 1),
                band_minimax=dict(pts[i_band], interior=0 < i_band < n - 1),
                joint_minimax=dict(pts[i_joint], interior=0 < i_joint < n - 1),
                joint_gain=float(pts[i_pre]["joint"] - pts[i_joint]["joint"]))


def main():
    out = dict(rule="eps*_mm = argmin_eps max(max_b |dm_b|, max_c |dcolour_c|) over the live bands and colours of the record's R2; "
                    "the preregistered eps* = argmin_eps mean_b |dm_b| (paperB/prl_gate.md, G1)", per_record={})
    for ion in IONS:
        name, row = g1_record(ion)
        if row is None:
            continue
        out["per_record"][ion] = dict(record=f"paperB/gate1/{name}", n=row["n"], **read(row))
    pc = ROOT / "paperB/gate3/gate3_partc_p1blend.json"
    if pc.exists():
        row = json.loads(pc.read_text())
        out["per_record"]["blend13"] = dict(record="paperB/gate3/gate3_partc_p1blend.json", n=row["n"], **read(row))
    (HERE / "minimax.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for key, r in out["per_record"].items():
        p, j = r["prereg"], r["joint_minimax"]
        print(f"{NAME[key]}: prereg eps* {p['eps']:.2f} (mean {p['mean_band']:.3f}, max band {p['max_band']:.3f}, colour {p['max_colour']:.3f}); "
              f"joint minimax eps {j['eps']:.2f} -> {j['joint']:.3f} (band {j['max_band']:.3f}, colour {j['max_colour']:.3f})")
    return out


if __name__ == "__main__" and "--markdown" not in sys.argv:
    main()


def markdown(out=None):
    out = out or json.loads((HERE / "minimax.json").read_text())
    L = ["| record | live bands | preregistered ε* (mean-optimal) | its max band / colour | band-only minimax ε | its max band / colour | joint minimax ε*_mm | its max band / colour | E_joint at ε*_mm |",
         "|---|---|---|---|---|---|---|---|---|"]
    for key, r in out["per_record"].items():
        p, b, j = r["prereg"], r["band_minimax"], r["joint_minimax"]
        L.append(f"| {NAME[key]} | {' '.join(r['live_bands'])} | {p['eps']:.2f} | {p['max_band']:.3f} / {p['max_colour']:.3f} | {b['eps']:.2f} | "
                 f"{b['max_band']:.3f} / {b['max_colour']:.3f} | {j['eps']:.2f} | {j['max_band']:.3f} / {j['max_colour']:.3f} | {j['joint']:.3f} |")
    return "\n".join(L)


if __name__ == "__main__" and "--markdown" in sys.argv:
    print(markdown())
