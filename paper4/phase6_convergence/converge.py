"""Paper IV Phase 6: numerical convergence of R2 and B2 on a benchmark zone.

Sweeps, one axis at a time about the production settings
(tau_min = 1e-3, dnu/nu = 4.17e-5, n = 3e5, no table cut):

    tau_min     1e-2, 1e-3, 1e-4          (the opacity-line cut; atom rebuilt)
    dnu_over_nu 4.17e-6, 4.17e-5, 4.17e-4 (1.25, 12.5, 125 km/s bins; B2 only)
    n           1e5, 2e5, 4e5              (packets per seed)
    a_cut       1e-3, 1e-4, 1e-5           (macroatom table cut)

and reports |dm| of every band between each setting and the finest one on
that axis. The plan's criterion is |dm| < 0.03-0.05 mag in the bands used
for a headline; where that is not met the measured spread is the quoted
numerical uncertainty and Gate 2 turns gray.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper4/phase2_energy"):
    sys.path.insert(0, str(p))

from sobolev.ejecta import EjectaState          # noqa: E402
import legs                                     # noqa: E402

AXES = dict(tau_min=(1e-2, 1e-3, 1e-4), dnu_over_nu=(4.17e-6, 4.17e-5, 4.17e-4),
            n=(100000, 200000, 400000), a_cut=(1e-3, 1e-4, 1e-5))
FINEST = dict(tau_min=1e-4, dnu_over_nu=4.17e-6, n=400000, a_cut=1e-5)
BASE = dict(tau_min=1e-3, dnu_over_nu=4.17e-5, n=300000, a_cut=None)


def one(state, shell, tau_min, dnu, n, a_cut, seeds, which=("R2", "B2")):
    atom, n_ion = legs.atom_for_zone(state, shell, tau_min=tau_min)
    zone = state.local_zone(shell)
    row = legs.run_legs(zone, atom, n, which, seeds, dnu_over_nu=dnu, a_cut=a_cut, verbose=False)
    return {t: dict(mags=row["legs"][t]["mags"], std=row["legs"][t]["mags_seed_std"],
                    dm_vs_R2=row["legs"][t]["dm_vs_R2"], t_wall=row["timing"][t]) for t in which} | \
        dict(n_opacity=row["n_opacity"], tau_min=tau_min, dnu_over_nu=dnu, n=n, a_cut=a_cut)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state")
    ap.add_argument("--shell", type=int, default=None)
    ap.add_argument("--axes", default="tau_min,dnu_over_nu,n,a_cut")
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    st = EjectaState.from_json(a.state)
    shell = st.meta.get("photospheric_shell", 0) if a.shell is None else a.shell
    seeds = tuple(int(s) for s in a.seeds.split(","))
    out = {"state": a.state, "shell": shell, "base": BASE, "axes": {}}
    t0 = time.time()
    for axis in a.axes.split(","):
        rows = []
        for val in AXES[axis]:
            kw = dict(BASE); kw[axis] = val
            r = one(st, shell, kw["tau_min"], kw["dnu_over_nu"], kw["n"], kw["a_cut"], seeds)
            rows.append(r)
            print(f"{axis}={val:g}: n_op={r['n_opacity']} R2 g={r['R2']['mags']['g']:.2f} K={r['R2']['mags']['K']:.2f} "
                  f"B2-R2 z={r['B2']['dm_vs_R2']['z']:+.2f} K={r['B2']['dm_vs_R2']['K']:+.2f} "
                  f"({r['R2']['t_wall']:.0f}+{r['B2']['t_wall']:.0f}s)", flush=True)
        fine = rows[[r[axis] for r in rows].index(FINEST[axis])]
        for r in rows:
            r["dm_to_finest"] = {t: {b: r[t]["mags"][b] - fine[t]["mags"][b] for b in r[t]["mags"]} for t in ("R2", "B2")}
            r["dclosure_to_finest"] = {b: r["B2"]["dm_vs_R2"][b] - fine["B2"]["dm_vs_R2"][b] for b in r["B2"]["dm_vs_R2"]}
        out["axes"][axis] = rows
    out["t_wall"] = time.time() - t0
    path = a.out or HERE / f"converge_{st.meta.get('name')}_t{st.t / 86400:g}_s{shell}.json"
    Path(path).write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(f"wrote {path} in {out['t_wall']:.0f}s")


if __name__ == "__main__":
    main()
