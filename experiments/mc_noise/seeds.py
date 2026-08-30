"""Seed-matched SEDONA pairs for every quoted band flux (referee Comment 7).

The production runs were seeded from time(NULL) with the resolved and
expansion legs launched minutes apart, so their sampling errors are
independent and the manuscript's "partially correlated" was unjustified as
written. With `transport_fix_rng_seed = 1` and the SAME seed integer in both
legs the core-emission stream is identical -- same packets, positions,
directions, frequencies -- and the streams diverge only at the first
interaction. From the five pairs already in `run.py`: corr(F_bb, F_exp) =
+0.97 and the paired Delta_exp scatter is 0.11%, against 0.40% from adding
the legs in quadrature. So matched seeds turn the claim into a measurement.

This script extends the pairs to everything the paper quotes:
  headline  tau_max=5, v_D=100     seeds 1-10 (1-5 exist), bb + exp
  12-grid   tau {0.5,5,50} x v_D {300,100,30,10}   seeds 1-3, bb + exp
  frontier  tau 5, v_D 3 km/s      seeds 2-3 (seed 1 is the production run)
The headline becomes the seed mean +- sem; the time-seed production runs are
kept as unpaired extra samples.

dnu/nu follows production: 8 transport bins per Doppler width. Runs that
already have a spectrum are not repeated, so this is safe to re-launch.
Results go to mc_noise_seeds.json; `analyze.py` reduces them.

Usage:  nohup <abs venv python> -u seeds.py [--workers 6] > seeds.log 2>&1 &
"""

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.spectra import band_ratio

from sobolev.sedona import sedona_cmd, sedona_home, sedona_timeout

SEDONA_HOME = sedona_home()
FOREST = ROOT / "experiments/laII_forest"
R_CORE, T_CORE = 8.64e12, 6000.0
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)

PARAM = """sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "{forest}/atom_laII.hdf5"
grid_type    = "grid_1D_sphere"
model_file   = "{model}"
hydro_module = "homologous"
transport_nu_grid  = {{7.50e14, 7.95e14, {dnu:.4e}, 1}}
spectrum_nu_grid   = {{7.50e14, 7.95e14, 1.0e-4, 1}}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1
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
line_velocity_width         = {vd:.3e}
output_write_radiation = 0
"""


def dnu_for(vd_kms, bins_per_width=8):
    return vd_kms * 1e5 / C / bins_per_width


# (label, model, v_D km/s, seeds)
CASES = [("tau5_vd100", FOREST / "sweep_tau5.mod", 100.0, range(1, 11))]
for tau in (0.5, 5.0, 50.0):
    for vd in (300.0, 100.0, 30.0, 10.0):
        if tau == 5.0 and vd == 100.0:
            continue  # the headline case above
        CASES.append((f"tau{tau:g}_vd{vd:g}", FOREST / f"sweep_tau{tau:g}.mod",
                      vd, range(1, 4)))
CASES.append(("tau5_vd3", FOREST / "sweep_tau5.mod", 3.0, range(2, 4)))

# The five existing tau5_vd10 bb runs are reused under their old directory
# names; only their exp partners are new.
ALIASES = {("tau5_vd10", "bb"): "tau5_vd10"}


def run_dir(label, mode, seed):
    return HERE / f"run_{label}_{mode}_s{seed}"


def run_one(label, model, vd, mode, seed):
    run = run_dir(label, mode, seed)
    run.mkdir(parents=True, exist_ok=True)
    spec = run / "spectrum_1.dat"
    if not spec.exists():
        (run / "param.lua").write_text(PARAM.format(
            forest=FOREST, model=model, dnu=dnu_for(vd), seed=seed, vd=vd * 1e5,
            bb=1 if mode == "bb" else 0, exp=0 if mode == "bb" else 1,
        ))
        r = subprocess.run(sedona_cmd(), cwd=run, capture_output=True,
                           text=True,
                           env={**os.environ, "SEDONA_HOME": SEDONA_HOME},
                           timeout=sedona_timeout(40000))
        if r.returncode != 0 or not spec.exists():
            (run / "run.log").write_text(r.stdout[-4000:] + "\n" + r.stderr[-2000:])
            return label, mode, seed, None
    return label, mode, seed, band_ratio(spec, BAND, MARGIN, R_CORE, T_CORE)


def main(workers):
    jobs = [(label, model, vd, mode, seed)
            for label, model, vd, seeds in CASES
            for seed in seeds for mode in ("bb", "exp")]
    # longest first so the pool tail is short
    cost = {300.0: 1, 100.0: 1, 30.0: 3, 10.0: 10, 3.0: 30}
    jobs.sort(key=lambda j: -cost[j[2]])
    print(f"{len(jobs)} runs, {workers} workers", flush=True)

    out_path = HERE / "mc_noise_seeds.json"
    results = json.loads(out_path.read_text()) if out_path.exists() else {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(run_one, *j) for j in jobs]
        for f in as_completed(futs):
            label, mode, seed, val = f.result()
            print(f"{label:14s} {mode:3s} seed={seed:2d}: {val}", flush=True)
            results.setdefault(label, {}).setdefault(mode, {})[str(seed)] = val
            out_path.write_text(json.dumps(results, indent=1, sort_keys=True))
    print("done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    main(ap.parse_args().workers)
