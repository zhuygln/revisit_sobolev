"""Real GSI La II forest, 3850-3950 A window (Weeks 3-4 first experiment).

Population control (research_requirements.md section 17: same populations in
both codes): the SEDONA atom file is built FROM the GSI La II data -- all 472
levels, window lines only -- with the ionization potential set to 1e5 eV so
Saha keeps 100% of the element in this stage. SEDONA's LTE Boltzmann then
runs over exactly the level list my sobolev.populations module uses.

Conditions: T = 3000 K, t = 1 day, shell 1000-3000 km/s, v_D = 100 km/s
(controlled; the 0.6 km/s thermal width is the expensive frontier), n_ion
chosen so the strongest window line has tau_S = 5.
"""

import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.atomic_data import load_gsi
from sobolev.constants import SIGMA_CLASSICAL
from sobolev.populations import boltzmann_fractions_from_levels

M_P = 1.67262192e-24
CM1_TO_EV = 1.239841984e-4
A_LA = 139.0

WINDOW = (3850.0, 3950.0)
T_SHELL = 3000.0
T_EXP = 86400.0
V_CORE, V_MAX = 1.0e8, 3.0e8
T_CORE = 6000.0
N_ZONES = 100
TAU_MAX_TARGET = 5.0

root = Path(__file__).resolve().parents[2]
levels = load_gsi(root / "data/57LaII_levels_calib.txt")
lines = load_gsi(root / "data/57LaII_transitions_calib.txt")

lam = lines["WV_Transition"].to_numpy()
sel = (lam >= WINDOW[0]) & (lam < WINDOW[1])
win = lines[sel].reset_index(drop=True)
print(f"{len(win)} lines in {WINDOW[0]:.0f}-{WINDOW[1]:.0f} A")

# Boltzmann populations over the full GSI level list.
frac = boltzmann_fractions_from_levels(levels, T_SHELL)
pop = frac[win["Lower"].to_numpy()]
g_l = 2 * win["J_Lower"].to_numpy() + 1
f_lu = 10 ** win["Log(gf)"].to_numpy() / g_l
tau_per_n = SIGMA_CLASSICAL * f_lu * pop * win["WV_Transition"].to_numpy() * 1e-8 * T_EXP
n_ion = TAU_MAX_TARGET / tau_per_n.max()
rho = n_ion * A_LA * M_P
taus = tau_per_n * n_ion
print(f"n_ion = {n_ion:.1f} cm^-3, rho = {rho:.3e} g/cm^3")
print(f"tau > 1: {(taus > 1).sum()},  0.1-1: {((taus > 0.1) & (taus <= 1)).sum()},  <0.1: {(taus <= 0.1).sum()}")

# Cross-check f(log gf) against f(A): the dataset should be self-consistent.
f_from_a = win["A"].to_numpy() * win["WV_Transition"].to_numpy() ** 2 \
    * (2 * win["J_Upper"].to_numpy() + 1) / (6.6702e15 * g_l)
worst = np.max(np.abs(f_from_a / f_lu - 1))
print(f"max |f(A)/f(gf) - 1| in window = {worst:.2e}")

# --- SEDONA atom file: all 472 levels, window lines only, no ionization ---
n_lev = len(levels)
with h5py.File(Path(__file__).parent / "atom_laII.hdf5", "w") as f:
    g = f.create_group("57")
    g.attrs["n_ions"] = np.int64(2)
    g.attrs["n_levels"] = np.int64(n_lev + 1)
    g.attrs["n_lines"] = np.int64(len(win))
    g.create_dataset("ion_chi", data=np.array([9.9999e4, 9.9999e4]))
    g.create_dataset("ion_ground", data=np.array([0, n_lev], dtype=np.int64))
    g.create_dataset(
        "level_E",
        data=np.concatenate([levels["Energy"].to_numpy() * CM1_TO_EV, [0.0]]),
    )
    g.create_dataset(
        "level_g",
        data=np.concatenate(
            [(2 * levels["J"].to_numpy() + 1).astype(np.int64), [1]]
        ),
    )
    g.create_dataset(
        "level_i", data=np.array([0] * n_lev + [1], dtype=np.int64)
    )
    g.create_dataset("line_A", data=win["A"].to_numpy())
    g.create_dataset("line_l", data=win["Lower"].to_numpy().astype(np.int64))
    g.create_dataset("line_u", data=win["Upper"].to_numpy().astype(np.int64))

# --- model file ---
v_edges = np.linspace(V_CORE, V_MAX, N_ZONES + 1)[1:]
with open(Path(__file__).parent / "laII.mod", "w") as fh:
    fh.write("1D_sphere standard\n")
    fh.write(f"{N_ZONES}\t{V_CORE * T_EXP:.6e}\t{T_EXP:.6e} 1 \n")
    fh.write("57.139\n")
    for v in v_edges:
        fh.write(f"{v * T_EXP:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} 1.0\n")

# core luminosity consistent with its blackbody temperature
SB = 5.670374e-5
r_core = V_CORE * T_EXP
lum = 4 * np.pi * r_core**2 * SB * T_CORE**4
print(f"r_core = {r_core:.4e} cm, core L = {lum:.4e} erg/s")

np.savez(
    Path(__file__).parent / "forest_lines.npz",
    lam=win["WV_Transition"].to_numpy(),
    f_lu=f_lu,
    pop=pop,
    tau=taus,
    n_ion=n_ion,
    rho=rho,
    core_lum=lum,
)
