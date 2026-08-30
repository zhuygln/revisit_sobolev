"""Validity-map rows 2-3: temperature sweep + thermal-width frontier.

Temperature axis: T in {2500, 3000, 4000, 5000} K at v_D = 100 km/s. At each
T the Boltzmann populations are recomputed and n_ion is RESCALED so the
strongest window line keeps tau_S = 5 -- this isolates population
redistribution (which lines are active, how they overlap) from the overall
strength scale already mapped by sweep.py. Note: at 5000 K the shell's
thermal emission at 3900 A reaches ~30% of the core surface brightness, so
absolute fluxes include fill-in; Delta_Sob remains a same-code differential.

Frontier axis: v_D in {3, 1} km/s at T = 3000 K, tau_max = 5, approaching
the 0.6 km/s La thermal width. Transport grids keep 8 bins per Doppler width
(46k / 140k bins) -- this also measures the resolved-mode cost curve.

Run from experiments/laII_forest/ after setup.py.
"""

import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.atomic_data import load_gsi
from sobolev.constants import SIGMA_CLASSICAL
from sobolev.populations import boltzmann_fractions_from_levels

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
from sobolev.sedona import sedona_cmd, sedona_home, sedona_timeout

SEDONA_HOME = sedona_home()

M_P = 1.67262192e-24
A_LA = 139.0
C_CM = 2.99792458e10
T_EXP = 86400.0
V_CORE, V_MAX = 1.0e8, 3.0e8
N_ZONES = 100
WINDOW = (3850.0, 3950.0)
TAU_MAX = 5.0

levels = load_gsi(ROOT / "data/57LaII_levels_calib.txt")
lines = load_gsi(ROOT / "data/57LaII_transitions_calib.txt")
lam_all = lines["WV_Transition"].to_numpy()
sel = (lam_all >= WINDOW[0]) & (lam_all < WINDOW[1])
win = lines[sel].reset_index(drop=True)
g_l = 2 * win["J_Lower"].to_numpy() + 1
f_lu = 10 ** win["Log(gf)"].to_numpy() / g_l
lam_cm = win["WV_Transition"].to_numpy() * 1e-8

PARAM = """sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "../atom_laII.hdf5"
grid_type    = "grid_1D_sphere"
model_file   = "../{mod}"
hydro_module = "homologous"
transport_nu_grid  = {{7.50e14, 7.95e14, {dnu:.3e}, 1}}
spectrum_nu_grid   = {{7.50e14, 7.95e14, 1.0e-4, 1}}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1
texp             = 86400.0
tstep_time_start = texp
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


def band_flux(spec_path):
    s = np.loadtxt(spec_path, comments="#")
    nu, lum = s[:, 0], s[:, 1]
    lam = C_CM / nu * 1e8
    red = (lam > 3952) & (lam < 3978) & (lum > 0)
    m = (lam > 3800) & (lam < 3955)
    order = np.argsort(lam[m])
    return np.trapezoid((lum[m] / np.mean(lum[red]))[order], lam[m][order]) / 155.0


def run_pair(tag, mod, vd):
    dnu = (vd / C_CM) / 8.0
    out = {}
    for mode, bb, ex in [("bb", 1, 0), ("exp", 0, 1)]:
        run = HERE / f"tsweep_{tag}_{mode}"
        run.mkdir(exist_ok=True)
        (run / "param.lua").write_text(
            PARAM.format(mod=mod, dnu=dnu, bb=bb, exp=ex, vd=vd)
        )
        t0 = time.time()
        r = subprocess.run(
            sedona_cmd(), cwd=run, capture_output=True, text=True,
            env={**os.environ, "SEDONA_HOME": SEDONA_HOME}, timeout=sedona_timeout(3000),
        )
        wall = time.time() - t0
        if r.returncode != 0:
            print(f"FAIL {tag} {mode} rc={r.returncode}")
            out[mode] = None
            continue
        out[mode] = band_flux(run / "spectrum_1.dat")
        out[f"{mode}_wall_s"] = round(wall, 1)
    return out


results = []

# ---- Temperature axis ----
for T in [2500.0, 3000.0, 4000.0, 5000.0]:
    frac = boltzmann_fractions_from_levels(levels, T)
    pop = frac[win["Lower"].to_numpy()]
    tau_per_n = SIGMA_CLASSICAL * f_lu * pop * lam_cm * T_EXP
    n_ion = TAU_MAX / tau_per_n.max()
    taus = tau_per_n * n_ion
    rho = n_ion * A_LA * M_P
    mod = f"tsweep_T{T:g}.mod"
    v_edges = np.linspace(V_CORE, V_MAX, N_ZONES + 1)[1:]
    with open(HERE / mod, "w") as fh:
        fh.write("1D_sphere standard\n")
        fh.write(f"{N_ZONES}\t{V_CORE * T_EXP:.6e}\t{T_EXP:.6e} 1 \n")
        fh.write("57.139\n")
        for v in v_edges:
            fh.write(f"{v * T_EXP:.6e} {v:.6e} {rho:.6e} {T:.6e} 1.0\n")
    row = {
        "axis": "T", "T": T, "v_d_kms": 100.0, "n_ion": n_ion,
        "n_tau_gt1": int((taus > 1).sum()),
        "n_tau_01_1": int(((taus > 0.1) & (taus <= 1)).sum()),
    }
    row.update(run_pair(f"T{T:g}", mod, 1.0e7))
    if row.get("bb") and row.get("exp"):
        row["delta_sob"] = (row["exp"] - row["bb"]) / row["bb"]
        print(
            f"T={T:5g} K  n_ion={n_ion:8.1f}  lines(tau>1)={row['n_tau_gt1']:2d} "
            f"(0.1-1)={row['n_tau_01_1']:2d}  F_bb={row['bb']:.4f} "
            f"F_exp={row['exp']:.4f}  Delta_Sob={row['delta_sob']:+.1%}"
        )
    results.append(row)

# ---- Thermal-width frontier at T = 3000 K, tau_max = 5 ----
for vd in [3.0e5, 1.0e5]:  # 3 km/s, 1 km/s
    row = {"axis": "vd", "T": 3000.0, "v_d_kms": vd / 1e5}
    row.update(run_pair(f"T3000_vd{vd/1e5:g}", "tsweep_T3000.mod", vd))
    if row.get("bb") and row.get("exp"):
        row["delta_sob"] = (row["exp"] - row["bb"]) / row["bb"]
        print(
            f"v_D={vd/1e5:3g} km/s  F_bb={row['bb']:.4f}  F_exp={row['exp']:.4f}  "
            f"Delta_Sob={row['delta_sob']:+.1%}  "
            f"walls bb/exp = {row['bb_wall_s']}/{row['exp_wall_s']} s"
        )
    results.append(row)

(HERE / "tsweep_results.json").write_text(json.dumps(results, indent=1))
print("wrote tsweep_results.json")
