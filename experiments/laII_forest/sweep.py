"""Validity-map sweep over the La II window: tau_max x v_D.

For each condition, rho rescales all line strengths together (populations are
fixed Boltzmann at 3000 K) and line_velocity_width sets the resolved profile
width; the transport grid keeps 8 bins per Doppler width. Delta_Sob is the
band-averaged flux error of the expansion-opacity run relative to the
resolved run from the SAME code -- the cleanest differential, immune to
cross-code systematics.

Run from experiments/laII_forest/ after setup.py. Sequential: ~24 x 30 s.
"""

import json
import os
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SEDONA_HOME = os.environ.get("SEDONA_HOME", os.path.expanduser("~/personal/pubsed"))
SEDONA = os.environ.get("SEDONA_EXE", f"{SEDONA_HOME}/src/sedona6.ex")

M_P = 1.67262192e-24
A_LA = 139.0
T_EXP = 86400.0
V_CORE, V_MAX = 1.0e8, 3.0e8
T_SHELL = 3000.0
N_ZONES = 100

d = np.load(HERE / "forest_lines.npz")
N_ION_TAU5 = float(d["n_ion"])  # n_ion giving tau_max = 5 (setup.py)

TAU_MAXES = [0.5, 5.0, 50.0]
V_DS = [3.0e7, 1.0e7, 3.0e6, 1.0e6]  # 300, 100, 30, 10 km/s

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


def write_model(tau_max, name):
    rho = N_ION_TAU5 * (tau_max / 5.0) * A_LA * M_P
    v_edges = np.linspace(V_CORE, V_MAX, N_ZONES + 1)[1:]
    with open(HERE / name, "w") as fh:
        fh.write("1D_sphere standard\n")
        fh.write(f"{N_ZONES}\t{V_CORE * T_EXP:.6e}\t{T_EXP:.6e} 1 \n")
        fh.write("57.139\n")
        for v in v_edges:
            fh.write(f"{v * T_EXP:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} 1.0\n")


def band_flux(spec_path):
    s = np.loadtxt(spec_path, comments="#")
    nu, lum = s[:, 0], s[:, 1]
    lam = 2.99792458e10 / nu * 1e8
    red = (lam > 3952) & (lam < 3978) & (lum > 0)
    scale = np.mean(lum[red])  # continuum level in the line-free margin
    m = (lam > 3800) & (lam < 3955)
    order = np.argsort(lam[m])
    # band-average of L/L_cont; continuum is flat enough over the band that a
    # scalar red-side normalization suffices for a DIFFERENTIAL measure.
    return np.trapezoid((lum[m] / scale)[order], lam[m][order]) / (3955 - 3800)


results = []
for tau_max in TAU_MAXES:
    mod = f"sweep_tau{tau_max:g}.mod"
    write_model(tau_max, mod)
    for vd in V_DS:
        dnu = (vd / 2.99792458e10) / 8.0
        row = {"tau_max": tau_max, "v_d_kms": vd / 1e5}
        for mode, bb, ex in [("bb", 1, 0), ("exp", 0, 1)]:
            run = HERE / f"sweep_tau{tau_max:g}_vd{vd/1e5:g}_{mode}"
            run.mkdir(exist_ok=True)
            (run / "param.lua").write_text(
                PARAM.format(mod=mod, dnu=dnu, bb=bb, exp=ex, vd=vd)
            )
            r = subprocess.run(
                [SEDONA, "param.lua"], cwd=run, capture_output=True, text=True,
                env={**os.environ, "SEDONA_HOME": SEDONA_HOME},
                timeout=560,
            )
            if r.returncode != 0:
                row[mode] = None
                print(f"FAIL tau={tau_max} vd={vd/1e5} {mode}: rc={r.returncode}")
                continue
            row[mode] = band_flux(run / "spectrum_1.dat")
        results.append(row)
        if row["bb"] and row["exp"]:
            row["delta_sob"] = (row["exp"] - row["bb"]) / row["bb"]
            print(
                f"tau_max={tau_max:5g} v_D={vd/1e5:4g} km/s  "
                f"F_bb={row['bb']:.4f}  F_exp={row['exp']:.4f}  "
                f"Delta_Sob={row['delta_sob']:+.1%}"
            )

(HERE / "sweep_results.json").write_text(json.dumps(results, indent=1))
print("wrote sweep_results.json")
