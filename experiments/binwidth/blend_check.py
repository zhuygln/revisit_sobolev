"""Reach the many-lines-per-bin regime with a reference that is still valid.

The La II window gives at most ~2.5 lines/bin before the resolved leg stops
resolving profiles (bin > v_D). The La+Ce blend has 2529 lines in the same
window -- 16x denser -- so at the same bin widths it reaches tens of lines per
bin while the reference remains converged. That is the regime the
expansion-opacity construction was actually designed for.
"""
import json, os, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.spectra import band_ratio
SEDONA_HOME = os.environ.get("SEDONA_HOME", os.path.expanduser("~/personal/pubsed"))
SEDONA = os.environ.get("SEDONA_EXE", f"{SEDONA_HOME}/src/sedona6.ex")
BLEND = ROOT / "experiments/multiion"
BAND, MARGIN, RC, TC = (3800.0, 3955.0), (3952.0, 3970.0), 8.64e12, 6000.0
C = 2.99792458e10
PARAM = """sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "{b}/atom_multiion.hdf5"
grid_type    = "grid_1D_sphere"
model_file   = "{b}/multiion.mod"
hydro_module = "homologous"
transport_nu_grid  = {{7.50e14, 7.95e14, {dnu:.4e}, 1}}
spectrum_nu_grid   = {{7.50e14, 7.95e14, 1.0e-4, 1}}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1
texp             = 86400.0
tstep_time_start = 86400.0
core_n_emit      = 4e6
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
out = []
print(" bin [km/s]  lines/bin |  F_res   F_exp   Delta_exp")
for dnu in (4.17e-5, 1.25e-4, 4.17e-4):
    bk = dnu*C/1e5; lpb = 2529*bk/7700.0
    row = {"bin_kms": bk, "lines_per_bin": lpb}
    for mode, bb, ex in (("bb",1,0),("exp",0,1)):
        run = HERE/f"blend_{dnu:.2e}_{mode}"; run.mkdir(exist_ok=True)
        (run/"param.lua").write_text(PARAM.format(b=BLEND,dnu=dnu,bb=bb,exp=ex))
        r = subprocess.run([SEDONA,"param.lua"],cwd=run,capture_output=True,
            text=True,env={**os.environ,"SEDONA_HOME":SEDONA_HOME},timeout=6000)
        row[mode] = band_ratio(run/"spectrum_1.dat",BAND,MARGIN,RC,TC) if r.returncode==0 else None
    if row.get("bb") and row.get("exp"):
        row["d_exp"]=(row["exp"]-row["bb"])/row["bb"]
        print(f"{bk:10.2f} {lpb:10.1f} | {row['bb']:.4f} {row['exp']:.4f} {row['d_exp']:+9.1%}",flush=True)
    out.append(row)
    (HERE/"blend_binwidth.json").write_text(json.dumps(out,indent=1))
