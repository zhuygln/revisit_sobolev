"""One energy-conserving check: the La II forest with radiative equilibrium on.

The paper's +45% is the attenuation of an externally supplied continuum, not
an emergent spectrum. The referee (Comment 4) is right that nothing measured
there bounds what happens once absorbed energy is re-emitted, and that the
earlier "upper bound" and "sign survives" arguments were not tested. This
runs the test: the same resolved vs expansion pair with
`transport_radiative_equilibrium = 1`, and reports the emergent-band
differential.

What RE does in SEDONA (verified in source): every absorption becomes an
effective scatter -- isotropic re-emission at a frequency drawn from the
zone emissivity (resolved: A_ul n_u phi over the lines; expansion: kappa_exp
B_nu(T) per bin) -- and after each iteration T_gas is re-solved from the
absorbed-energy tally. Iteration 1 re-emits at the INPUT T = 3000 K with the
populations unchanged, which isolates redistribution from population
feedback; later iterations let T move, and the two modes converge to
DIFFERENT temperatures, so that differential includes feedback. Both are
reported.

Runs (fixed, matched seeds so the differential is paired):
  N=1   bb, exp     seeds 1-3     redistribution only
  N=15  bb, exp     seeds 1-3     self-consistent T
  null: tau_max=0.05 (rho x 0.01), N=15, bb + exp, seed 1 -- both modes must
        agree with each other and with the RE-off continuum within noise.

Normalization: the red margin 3952-3970 A is contaminated by re-emission
redshifted up to 3000 km/s (the 3950 A lines land at ~3990 A), so the band
ratio uses a BLUE margin 3785-3805 A; re-emission cannot reach blueward of
3850 x (1 - 0.01) = 3811 A. On the RE-off production runs the blue- and
red-margin ratios must agree within noise (the control on this choice).

Usage:  nohup <abs venv python> -u re_run.py [--workers 3] > re_run.log 2>&1 &
"""

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.spectra import band_ratio

SEDONA = "/home/yozhuz_223/personal/pubsed/src/sedona6.ex"
SEDONA_HOME = "/home/yozhuz_223/personal/pubsed"
R_CORE, T_CORE = 8.64e12, 6000.0
BAND = (3800.0, 3955.0)
BLUE_MARGIN, RED_MARGIN = (3785.0, 3805.0), (3952.0, 3970.0)

PARAM = """sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "{here}/atom_laII.hdf5"
grid_type    = "grid_1D_sphere"
model_file   = "{model}"
hydro_module = "homologous"
transport_nu_grid  = {{7.50e14, 7.95e14, 4.17e-5, 1}}
spectrum_nu_grid   = {{7.50e14, 7.95e14, 1.0e-4, 1}}
transport_radiative_equilibrium = 1
transport_steady_iterate        = {n_iter}
transport_fix_rng_seed          = 1
transport_rng_seed              = {seed}
texp             = 86400.0
tstep_time_start = 86400.0
core_n_emit      = 2e6
core_radius      = 8.64e12
core_temperature = 6000.0
core_luminosity  = 6.8937e37
opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = {bb}
opacity_line_expansion      = {exp}
opacity_epsilon             = 1
line_velocity_width         = 1.0e7
output_write_radiation = 0
"""
# plt_0000N.dat is written each iteration (default output_write_plt_file_time
# = 1); column 4 is T_gas, which is how the T convergence is read off.


def null_model():
    """tau_max = 0.05 version of the production model: rho x 0.01."""
    src, dst = HERE / "laII.mod", HERE / "laII_null.mod"
    if dst.exists():
        return dst
    out = []
    for line in src.read_text().splitlines():
        parts = line.split()
        if len(parts) == 5:
            try:
                vals = [float(x) for x in parts]
                vals[2] *= 0.01
                line = " ".join(f"{v:.6e}" for v in vals)
            except ValueError:
                pass
        out.append(line)
    dst.write_text("\n".join(out) + "\n")
    return dst


def run_one(tag, model, n_iter, mode, seed):
    run = HERE / f"re_{tag}_N{n_iter}_{mode}_s{seed}"
    run.mkdir(parents=True, exist_ok=True)
    last = run / f"spectrum_{n_iter}.dat"
    if not last.exists():
        (run / "param.lua").write_text(PARAM.format(
            here=HERE, model=model, n_iter=n_iter, seed=seed,
            bb=1 if mode == "bb" else 0, exp=0 if mode == "bb" else 1,
        ))
        r = subprocess.run([SEDONA, "param.lua"], cwd=run, capture_output=True,
                           text=True,
                           env={**os.environ, "SEDONA_HOME": SEDONA_HOME},
                           timeout=60000)
        (run / "run.log").write_text(r.stdout[-6000:] + "\n" + r.stderr[-2000:])
        if r.returncode != 0 or not last.exists():
            return tag, n_iter, mode, seed, None
    per_iter = []
    for it in range(1, n_iter + 1):
        p = run / f"spectrum_{it}.dat"
        if p.exists():
            per_iter.append(band_ratio(p, BAND, BLUE_MARGIN, R_CORE, T_CORE))
    return tag, n_iter, mode, seed, per_iter


def main(workers):
    prod = HERE / "laII.mod"
    jobs = []
    for n_iter in (1, 15):
        for seed in (1, 2, 3):
            for mode in ("bb", "exp"):
                jobs.append(("prod", prod, n_iter, mode, seed))
    nm = null_model()
    for mode in ("bb", "exp"):
        jobs.append(("null", nm, 15, mode, 1))
    jobs.sort(key=lambda j: -j[2])
    print(f"{len(jobs)} runs, {workers} workers", flush=True)

    # the margin control, on the RE-off production spectra
    ctrl = {}
    for mode in ("bb", "exp"):
        spec = HERE / f"run_{mode}" / "spectrum_1.dat"
        if spec.exists():
            ctrl[mode] = {
                "blue_margin": band_ratio(spec, BAND, BLUE_MARGIN, R_CORE, T_CORE),
                "red_margin": band_ratio(spec, BAND, RED_MARGIN, R_CORE, T_CORE),
            }
            print(f"RE-off {mode}: blue-margin {ctrl[mode]['blue_margin']:.5f}  "
                  f"red-margin {ctrl[mode]['red_margin']:.5f}", flush=True)

    out_path = HERE / "re_results.json"
    results = {"margin_control": ctrl, "runs": {}}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(run_one, *j) for j in jobs]
        for f in as_completed(futs):
            tag, n_iter, mode, seed, vals = f.result()
            key = f"{tag}_N{n_iter}_{mode}_s{seed}"
            results["runs"][key] = vals
            shown = "FAILED" if vals is None else f"{vals[-1]:.5f} (per-iter {len(vals)})"
            print(f"{key:22s}: {shown}", flush=True)
            out_path.write_text(json.dumps(results, indent=1, sort_keys=True))
    print("done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    main(ap.parse_args().workers)
