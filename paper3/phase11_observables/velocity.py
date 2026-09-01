"""Paper III §10b: how the photometric closure error scales with ejecta speed.

`trajectory.py` inherits Paper I's shell, 1000-3000 km/s. Sobolev optical depth
does not care -- tau ~ n t is velocity-free in homologous flow -- but the
wavelength interval a packet sweeps before escaping does: d ln lambda = dv/c.
At 0.003-0.01c that is 0.7%; real kilonova ejecta sweep ~20%. Since the swept
interval is exactly what a grouped closure coarse-grains, the F40 geometry
should be a LOWER bound on the closure error, and §4.37's magnitudes with it.

Rather than repeat the whole trajectory at kilonova speed -- which is expensive
and confounds velocity with epoch -- this holds the epoch, the density, the
composition and the temperature fixed and moves only the shell's velocity. The
shell shape is held fixed too (v_core = v_out / 3, as in Paper I), so the only
thing changing is how far in frequency a packet travels.

`relativity="worldline"` throughout: F23 showed the time-frozen classical
transport is not trustworthy at 0.1c, and the difference here is larger than
the effect being measured.

Usage: python velocity.py [--ion 57LaII] [--t 1.0] [--n 150000]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for p in (HERE,):
    sys.path.insert(0, str(p))
from observables import LEGS, epoch

BETAS = (0.01, 0.03, 0.06, 0.10, 0.20, 0.30)


def main(ion, t_d, n, out_name, core_law="cool"):
    print(f"{ion} at t = {t_d} d, worldline, v_out = {BETAS} c "
          f"(v_core = v_out/3), {n} packets", flush=True)
    print(f"{'v_out/c':>8}{'dlnlam':>8}{'S_band':>9}{'dF3800':>9}{'dm_bol':>8}"
          f"{'|dm|max':>9}{'|dcol|max':>10}  worst", flush=True)
    rows, t0 = [], time.time()
    for b in BETAS:
        r = epoch(ion, t_d, n, core_law, vel=(b / 3.0, b), relativity="worldline")
        if "skipped" in r:
            print(f"{b:8.3f}  skipped -- {r['skipped']}", flush=True); continue
        L = r["legs"]["C_both"]
        dm = {k: v for k, v in L["dm"].items() if np.isfinite(v)}
        dc = {k: v for k, v in L["dcolor"].items() if np.isfinite(v)}
        wm = max(dm, key=lambda k: abs(dm[k])); wc = max(dc, key=lambda k: abs(dc[k]))
        r["worst"] = {"band": wm, "dm": dm[wm], "color": wc, "dcolor": dc[wc]}
        rows.append(r)
        f = L["dF_b3800"]
        fs = "    --  " if f is None else f"{100*f:+8.1f}%"
        print(f"{b:8.3f}{b*(1-1/3.0):8.3f}{r['band_S_band']:9.1f}{fs}"
              f"{L['dm_bol']:+8.3f}{abs(dm[wm]):9.3f}{abs(dc[wc]):10.3f}"
              f"  {wm} / {wc}", flush=True)
    (HERE / out_name).write_text(json.dumps(
        {"ion": ion, "t_d": t_d, "n": n, "core_law": core_law,
         "relativity": "worldline", "betas": list(BETAS), "rows": rows}, indent=1))
    print(f"wrote {out_name}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="57LaII")
    ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--n", type=float, default=1.5e5)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    main(a.ion, a.t, int(a.n), a.out or f"velocity_{a.ion}_t{a.t:g}.json")
