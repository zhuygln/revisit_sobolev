"""Single-line convergence ladder for SEDONA (referee Comment 6).

The published single-line numbers are analytic 0.1353 / solver 0.1372 /
SEDONA 0.1420. Two things are now known about them. First, SEDONA's
steady-iterate mode is a frozen-snapshot calculation (F11), whose law for one
line is e^{-tau_S (1-beta)/gamma}: over the 1400-2600 km/s trough window that
is 0.1371, and the solver in frozen mode reproduces it -- so 0.1372 was never
a discrepancy, and 0.1353 is the beta -> 0 limit, not the target. Second,
the one thing left is SEDONA's +3.5% above that, and the production grid
puts dnu/nu = 2e-4 (~60 km/s bins) against a 100 km/s Doppler width: under
two bins per profile, against eight in every forest run, and the forest
bin-width study already shows the resolved leg degrading past that.

So the ladder converges SEDONA toward 0.1371 one axis at a time, with fixed
matched seeds, and reports the residual whatever it is:

  A  packets      5e5 2e6 8e6 3.2e7          at dnu_t 2e-4
  B  transport    4e-4 2e-4 1e-4 4.17e-5 2e-5 at 8e6
  C  spectrum     5e-4 2e-4 1e-4              at 8e6, dnu_t 4.17e-5
  D  zones        25 100 101 400              (101 = half-zone shift: line vs cell edge)
  E  anchor       3.2e7, 4.17e-5, 2e-4, 100 zones, 5 seeds, bb AND exp
  Z  zero opacity (rho x 1e-6) at the anchor grid: trough must be 1.000

Seeds 1-3 per rung (5 at E). particles_max_total raised for the 3.2e7 rung
or SEDONA silently emits nothing. Common random numbers break when the
transport grid changes (frequency sampling uses the grid CDF), so rung B's
seeds are not paired across its steps -- say so.
"""
import argparse, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.formal_transfer import planck_bnu

from sobolev.sedona import sedona_cmd, sedona_home, sedona_timeout

SEDONA_HOME = sedona_home()
M_P = 1.67262192e-24; SIGMA_CLASSICAL = 0.026540083433884684
T_EXP = 20 * 86400.0; V_CORE, V_MAX = 1.0e8, 3.0e8; T_SHELL = 2000.0
F_LU = 0.6647; LAMBDA0_CM = 12398.42 / 10.2 * 1e-8; NU0 = C / LAMBDA0_CM
TAU_S = 2.0; R_CORE = V_CORE * T_EXP; T_CORE = 2.0e4
N_H = TAU_S / (SIGMA_CLASSICAL * F_LU * LAMBDA0_CM * T_EXP)

PARAM = """sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/2level_atomdata.hdf5"
grid_type    = "grid_1D_sphere"
model_file   = "{model}"
hydro_module = "homologous"
transport_nu_grid  = {{2.30e15, 2.62e15, {dnu_t:.3e}, 1}}
spectrum_nu_grid   = {{2.30e15, 2.62e15, {dnu_s:.3e}, 1}}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1
transport_fix_rng_seed          = 1
transport_rng_seed              = {seed}
particles_max_total = {pmax:.0e}
texp             = {texp:.6e}
tstep_time_start = texp
core_n_emit      = {n_emit:.0e}
core_radius      = {rcore:.6e}
core_temperature = 2.0e4
core_luminosity  = 3.40e42
opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = {bb}
opacity_line_expansion      = {exp}
opacity_epsilon             = 1
line_velocity_width         = 1.0e7
output_write_radiation = 0
"""


def model(n_zones, rho_scale=1.0):
    p = HERE / f"ladder_model_n{n_zones}_r{rho_scale:g}.mod"
    if p.exists():
        return p
    rho = N_H * M_P * rho_scale
    v_edges = np.linspace(V_CORE, V_MAX, n_zones + 1)[1:]
    with open(p, "w") as fh:
        fh.write("1D_sphere standard\n")
        fh.write(f"{n_zones}\t{R_CORE:.6e}\t{T_EXP:.6e} 1 \n1.1\n")
        for v in v_edges:
            fh.write(f"{v*T_EXP:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} 1.0\n")
    return p


def trough(spec):
    s = np.loadtxt(spec, comments="#"); nu, lum = s[:, 0], s[:, 1]
    good = lum > 0; nu, lum = nu[good], lum[good]
    cont = 4 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    v = C * (nu / NU0 - 1.0) / 1e5
    red = (v > -25000) & (v < -500)
    ratio = lum / np.mean(lum[red] / cont[red]) / cont
    sel = (v > 1400) & (v < 2600)
    counts = s[good, 2][sel] if s.shape[1] > 2 else None
    return float(np.mean(ratio[sel])), (float(1 / np.sqrt(counts.sum())) if counts is not None and counts.sum() > 0 else None)


RUNGS = []  # (tag, dnu_t, dnu_s, n_emit, n_zones, rho_scale, mode, seeds)
for n in (5e5, 2e6, 8e6, 3.2e7):
    RUNGS.append((f"A_emit{n:.0e}", 2e-4, 5e-4, n, 100, 1.0, "bb", (1, 2, 3)))
for dt in (4e-4, 2e-4, 1e-4, 4.17e-5, 2e-5):
    RUNGS.append((f"B_dnut{dt:.2e}", dt, 5e-4, 8e6, 100, 1.0, "bb", (1, 2, 3)))
for ds in (5e-4, 2e-4, 1e-4):
    RUNGS.append((f"C_dnus{ds:.0e}", 4.17e-5, ds, 8e6, 100, 1.0, "bb", (1, 2, 3)))
for nz in (25, 100, 101, 400):
    RUNGS.append((f"D_nz{nz}", 4.17e-5, 2e-4, 8e6, nz, 1.0, "bb", (1, 2, 3)))
for mode in ("bb", "exp"):
    RUNGS.append((f"E_anchor_{mode}", 4.17e-5, 2e-4, 3.2e7, 100, 1.0, mode, (1, 2, 3, 4, 5)))
RUNGS.append(("Z_zero", 4.17e-5, 2e-4, 8e6, 100, 1e-6, "bb", (1,)))


def run_one(tag, dnu_t, dnu_s, n_emit, n_zones, rho_scale, mode, seed):
    run = HERE / f"ladder_{tag}_s{seed}"
    run.mkdir(exist_ok=True)
    spec = run / "spectrum_1.dat"
    if not spec.exists():
        (run / "param.lua").write_text(PARAM.format(
            model=model(n_zones, rho_scale), dnu_t=dnu_t, dnu_s=dnu_s, seed=seed,
            pmax=max(4e7, 2 * n_emit), texp=T_EXP, n_emit=n_emit, rcore=R_CORE,
            bb=1 if mode == "bb" else 0, exp=0 if mode == "bb" else 1))
        r = subprocess.run(sedona_cmd(), cwd=run, capture_output=True, text=True,
                           env={**os.environ, "SEDONA_HOME": SEDONA_HOME}, timeout=sedona_timeout(60000))
        (run / "run.log").write_text(r.stdout[-4000:] + "\n" + r.stderr[-2000:])
        if r.returncode != 0 or not spec.exists():
            return tag, seed, None, None
    t, err = trough(spec)
    return tag, seed, t, err


def main(workers):
    jobs = [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], s) for r in RUNGS for s in r[7]]
    jobs.sort(key=lambda j: -j[3])
    print(f"{len(jobs)} runs, {workers} workers; targets: frozen 0.1371, classical 0.1353", flush=True)
    out_path = HERE / "ladder_results.json"
    results = json.loads(out_path.read_text()) if out_path.exists() else {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(run_one, *j) for j in jobs]
        for f in as_completed(futs):
            tag, seed, t, err = f.result()
            results.setdefault(tag, {})[str(seed)] = dict(trough=t, poisson_err=err)
            print(f"  {tag:18s} seed {seed}: {t if t is None else f'{t:.5f}'}  (poisson {err if err is None else f'{err:.4f}'})", flush=True)
            out_path.write_text(json.dumps(results, indent=1, sort_keys=True))
    print("done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4)
    main(ap.parse_args().workers)
