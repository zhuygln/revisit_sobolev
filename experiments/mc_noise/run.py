"""P1-D: measure the Monte Carlo noise floor instead of asserting it.

THE PROBLEM. The paper claims the Sobolev approximation is accurate to
"<=0.5%" while separately stating a Monte Carlo noise floor of "~1-2% per band
flux". Those two statements cannot both be load-bearing: the Sobolev leg is
analytic, but the SEDONA run it is differenced against is not, so a quoted
-0.3% is meaningless if the reference carries +-1.5%. A referee will notice,
and the fix is not to reword the claim but to measure the uncertainty.

THE MEASUREMENT. Rerun identical configurations with different RNG seeds
(SEDONA's transport_fix_rng_seed / transport_rng_seed) and take the spread of
the band flux. Everything else -- model, atom, grid, packet count -- is held
fixed, so the spread is pure sampling noise.

TRANSFERABILITY. Noise is measured at v_D = 100 and v_D = 10 km/s. The band
flux is a sum over all bins in the window, so its sampling error is set by the
number of packets landing in the band, not by how finely the band is
subdivided; it should therefore be v_D-independent at fixed core_n_emit.
Measuring at two widths an order of magnitude apart TESTS that expectation
rather than assuming it, which is what licenses quoting the same uncertainty
at the v_D = 1 km/s frontier point (a run ~60x more expensive to repeat).
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.spectra import band_ratio

SEDONA_HOME = os.environ.get("SEDONA_HOME", os.path.expanduser("~/personal/pubsed"))
SEDONA = os.environ.get("SEDONA_EXE", f"{SEDONA_HOME}/src/sedona6.ex")
FOREST = ROOT / "experiments/laII_forest"
R_CORE, T_CORE = 8.64e12, 6000.0
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)

N_SEEDS = 5

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

# (label, model file, transport dnu/nu, v_D km/s, which legs)
# dnu/nu follows the production runs: 8 transport bins per Doppler width.
CASES = [
    ("tau5_vd100", FOREST / "sweep_tau5.mod", 4.17e-5, 100.0, ("bb", "exp")),
    ("tau5_vd10", FOREST / "sweep_tau5.mod", 4.17e-6, 10.0, ("bb",)),
]


def run_one(case, mode, seed):
    label, model, dnu, vd, _ = case
    run = HERE / f"run_{label}_{mode}_s{seed}"
    run.mkdir(parents=True, exist_ok=True)
    spec = run / "spectrum_1.dat"
    if spec.exists():
        return band_ratio(spec, BAND, MARGIN, R_CORE, T_CORE)
    (run / "param.lua").write_text(PARAM.format(
        forest=FOREST, model=model, dnu=dnu, seed=seed, vd=vd * 1e5,
        bb=1 if mode == "bb" else 0, exp=0 if mode == "bb" else 1,
    ))
    r = subprocess.run([SEDONA, "param.lua"], cwd=run, capture_output=True,
                       text=True, env={**os.environ, "SEDONA_HOME": SEDONA_HOME},
                       timeout=20000)
    if r.returncode != 0:
        print(f"  FAIL {label} {mode} seed={seed} rc={r.returncode}", flush=True)
        print(r.stdout[-1500:], flush=True)
        return None
    return band_ratio(spec, BAND, MARGIN, R_CORE, T_CORE)


results = {}
for case in CASES:
    label, _, _, vd, modes = case
    for mode in modes:
        vals = []
        for seed in range(1, N_SEEDS + 1):
            v = run_one(case, mode, seed)
            if v is not None:
                vals.append(v)
            print(f"{label} {mode} seed={seed}: {v}", flush=True)
        if len(vals) < 2:
            continue
        a = np.array(vals)
        results[f"{label}_{mode}"] = {
            "v_d_kms": vd, "mode": mode, "values": vals,
            "mean": float(a.mean()), "std": float(a.std(ddof=1)),
            "frac_std": float(a.std(ddof=1) / a.mean()),
            "n_seeds": len(vals),
        }
        print(f"  -> mean {a.mean():.5f}  std {a.std(ddof=1):.5f} "
              f"({a.std(ddof=1)/a.mean():.3%})", flush=True)
        (HERE / "mc_noise_results.json").write_text(json.dumps(results, indent=1))

print("\n=== MEASURED MONTE CARLO NOISE (1 sigma, fractional) ===")
for k, v in results.items():
    print(f"  {k:20s} v_D={v['v_d_kms']:5g}  {v['frac_std']:.3%} "
          f"({v['n_seeds']} seeds)")

# Does the noise depend on v_D? If not, the value measured here transfers to
# the frontier runs that are too expensive to repeat.
a, b = results.get("tau5_vd100_bb"), results.get("tau5_vd10_bb")
if a and b:
    ratio = b["frac_std"] / a["frac_std"]
    print(f"\n  v_D=10 / v_D=100 noise ratio: {ratio:.2f}")
    print("  => band-flux noise is v_D-independent as expected; the value "
          "above may be quoted at the 1 km/s frontier point"
          if 0.4 < ratio < 2.5 else
          "  => NOT v_D-independent -- the frontier points need their own "
          "seed repeats before any sub-percent claim is made")

# What the differentials inherit.
r, e = results.get("tau5_vd100_bb"), results.get("tau5_vd100_exp")
if r and e:
    # d = (F_x - F_res)/F_res. For the analytic Sobolev leg F_x is exact, so
    # only the reference fluctuates: sigma_d = (1+d) * sigma_res/F_res.
    print(f"\n  sigma on Delta_Sobolev  (analytic vs MC reference): "
          f"~{r['frac_std']:.2%}")
    both = np.hypot(r["frac_std"], e["frac_std"])
    print(f"  sigma on Delta_expansion (MC vs MC, uncorrelated bound): "
          f"~{both:.2%}")
    print("  NOTE: the two SEDONA legs share a model and differ only in "
          "opacity treatment, so their errors are partially correlated and "
          "the expansion figure is an upper bound.")
