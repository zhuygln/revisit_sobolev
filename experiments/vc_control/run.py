"""Step 1: the v/c control experiment.

A single isolated line in a shell of CONSTANT lower-level density, so the
nonrelativistic Sobolev depth tau_S is the same at every radius. Any
dependence of the measured trough on the resonance velocity is therefore the
O(v/c) effect and nothing else.

One spectrum supplies the whole sweep: an observer frequency nu resonates
where nu' = nu0, so scanning nu scans beta. That makes this experiment one
SEDONA run plus two solver calls, not eight of each.

Legs
----
  solver "first"  -- historical first-order mode
  solver "exact"  -- full SR (Doppler on opacity, gamma, source transform)
  SEDONA resolved bound-bound
  analytic        -- tau_S; tau_S(1-b); tau_S(1+b); tau_rel(b)

Geometry note: relativistically the resonance surface is not a plane
(gamma depends on r, not z) -- Jeffery's CD/CP surfaces. The core is made
small (beta_core = 1e-4) so all core rays are nearly radial and the plane
picture is a good approximation; this keeps the analytic legs interpretable.
"""

import subprocess
import sys
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C, SIGMA_CLASSICAL
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.optical_depth import tau_sobolev
from sobolev.sobolev_leg import tau_sobolev_relativistic

SEDONA = "/home/yozhuz_223/personal/pubsed/src/sedona6.ex"
SEDONA_HOME = "/home/yozhuz_223/personal/pubsed"
M_P = 1.67262192e-24

# --- the shipped 2-level atom: one line, f from A ---
LAMBDA0 = 12398.42 / 10.2 * 1e-8      # cm
NU0 = C / LAMBDA0
F_OSC = 0.6647
T_EXP = 86400.0
BETA_CORE, BETA_OUT = 1.0e-4, 0.35
R_CORE, R_OUT = BETA_CORE * C * T_EXP, BETA_OUT * C * T_EXP
T_CORE, T_SHELL = 2.0e4, 2000.0
V_D = 1.0e7                            # 100 km/s
TAU_TARGET = 1.0                       # sensitive, unsaturated
N_ZONES = 200

BETAS = np.array([0.001, 0.003, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3])

n_l = TAU_TARGET / (SIGMA_CLASSICAL * F_OSC * LAMBDA0 * T_EXP)
rho = n_l * M_P
print(f"n_l = {n_l:.4f} cm^-3 -> tau_S = {TAU_TARGET} everywhere")
print(f"shell beta {BETA_CORE:g} .. {BETA_OUT:g}, r = {R_CORE:.3e} .. {R_OUT:.3e} cm")

# Observer frequency resonating at beta, first-order convention nu = nu0/(1-b).
# (The exact condition differs at O(b^2); the legs are compared at the SAME
# observer frequency, so the labelling convention cancels.)
nu_probe = NU0 / (1.0 - BETAS)


def solver_transmission(mode):
    lum = emergent_luminosity(
        nu_probe, [(NU0, F_OSC)],
        lambda r: np.full_like(r, n_l), lambda r: np.full_like(r, T_SHELL),
        T_EXP, R_CORE, R_OUT, T_CORE, V_D,
        n_impact=120, relativity=mode,
    )
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu_probe, T_CORE)
    return lum / cont


# --------------------------------------------------------------- SEDONA
v_edges = np.geomspace(BETA_CORE * C, BETA_OUT * C, N_ZONES + 1)[1:]
with open(HERE / "vc.mod", "w") as fh:
    fh.write("1D_sphere standard\n")
    fh.write(f"{N_ZONES}\t{R_CORE:.6e}\t{T_EXP:.6e} 1 \n")
    fh.write("1.1\n")
    for v in v_edges:
        fh.write(f"{v*T_EXP:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} 1.0\n")

SB = 5.670374e-5
core_lum = 4 * np.pi * R_CORE**2 * SB * T_CORE**4
nu_lo, nu_hi = NU0 * 0.97, NU0 / (1.0 - 0.34)
PARAM = f"""sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/2level_atomdata.hdf5"
grid_type    = "grid_1D_sphere"
model_file   = "../vc.mod"
hydro_module = "homologous"
transport_nu_grid  = {{{nu_lo:.6e}, {nu_hi:.6e}, 4.17e-5, 1}}
spectrum_nu_grid   = {{{nu_lo:.6e}, {nu_hi:.6e}, 1.0e-4, 1}}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1
texp             = {T_EXP:.6e}
tstep_time_start = {T_EXP:.6e}
core_n_emit      = 1e7
core_radius      = {R_CORE:.6e}
core_temperature = {T_CORE:.1f}
core_luminosity  = {core_lum:.6e}
opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = 1
opacity_line_expansion      = 0
opacity_epsilon             = 1
line_velocity_width         = {V_D:.3e}
output_write_radiation = 0
"""
run = HERE / "run_bb"
run.mkdir(exist_ok=True)
(run / "param.lua").write_text(PARAM)
print("running SEDONA ...", flush=True)
r = subprocess.run([SEDONA, "param.lua"], cwd=run, capture_output=True, text=True,
                   env={**os.environ, "SEDONA_HOME": SEDONA_HOME}, timeout=5000)
print(f"  rc={r.returncode}")

s = np.loadtxt(run / "spectrum_1.dat", comments="#")
nu_s, lum_s = s[:, 0], s[:, 1]
cont_s = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu_s, T_CORE)
# normalize on the red side, well clear of the line and of both grid edges
red = (nu_s > NU0 * 0.975) & (nu_s < NU0 * 0.99) & (lum_s > 0)
ratio_s = lum_s / np.mean(lum_s[red] / cont_s[red]) / cont_s
sed = np.interp(nu_probe, nu_s, ratio_s)

# --------------------------------------------------------------- assemble
first = solver_transmission("first")
exact = solver_transmission("exact")
tau_S = tau_sobolev(F_OSC, n_l, LAMBDA0, T_EXP)
ana = {
    "tau_S": np.exp(-tau_S * np.ones_like(BETAS)),
    "tau_S(1-b)": np.exp(-tau_S * (1.0 - BETAS)),
    "tau_S(1+b)": np.exp(-tau_S * (1.0 + BETAS)),
    "tau_rel": np.exp(-tau_sobolev_relativistic(F_OSC, n_l, LAMBDA0, T_EXP, BETAS)),
}

print(f"\ntau_S = {tau_S:.4f}, exp(-tau_S) = {np.exp(-tau_S):.4f}\n")
hdr = f"{'beta':>7} {'solv_1st':>9} {'solv_ex':>9} {'SEDONA':>9} | " + \
      " ".join(f"{k:>11}" for k in ana)
print(hdr)
for i, b in enumerate(BETAS):
    print(f"{b:7.3f} {first[i]:9.4f} {exact[i]:9.4f} {sed[i]:9.4f} | " +
          " ".join(f"{ana[k][i]:11.4f}" for k in ana))

np.savez(HERE / "vc_results.npz", betas=BETAS, first=first, exact=exact,
         sedona=sed, tau_S=tau_S, **{f"ana_{k}": v for k, v in ana.items()})
print("\nwrote vc_results.npz")
