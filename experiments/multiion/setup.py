"""Multi-ion forest: La II + Ce II blended in the 3850-3950 A window.

Population control as in the single-ion experiment, per element: the SEDONA
atom file carries TWO element groups (Z=57, Z=58), each holding the full GSI
level list of its singly-ionized stage as "ion 0" with ionization disabled
(chi = 1e5 eV). Equal mass fractions X = 0.5 each; the total density is set
so the strongest line of EITHER species has tau_S = 5.

Ce III is extracted alongside for later but deliberately excluded here:
including two ionization stages of one element would hand the II/III split
to Saha, breaking exact population control. One stage per element keeps the
match to the Python side exact.
"""

import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.atomic_data import load_gsi
from sobolev.constants import SIGMA_CLASSICAL
from sobolev.populations import boltzmann_fractions_from_levels, statistical_weight

M_P = 1.67262192e-24
CM1_TO_EV = 1.239841984e-4
WINDOW = (3850.0, 3950.0)
T_SHELL = 3000.0
T_EXP = 86400.0
V_CORE, V_MAX = 1.0e8, 3.0e8
N_ZONES = 100
TAU_MAX_TARGET = 5.0

SPECIES = [
    # (Z, mass number, levels file, transitions file, mass fraction)
    (57, 139, "57LaII_levels_calib.txt", "57LaII_transitions_calib.txt", 0.5),
    (58, 140, "58CeII_levels_calib.txt", "58CeII_transitions_calib.txt", 0.5),
]

root = Path(__file__).resolve().parents[2]
here = Path(__file__).parent

# --- load, window-select, and compute tau per unit total mass density ---
per_species = []
for z, a, lev_f, tr_f, x_frac in SPECIES:
    levels = load_gsi(root / "data" / lev_f)
    lines = load_gsi(root / "data" / tr_f)
    lam = lines["WV_Transition"].to_numpy()
    win = lines[(lam >= WINDOW[0]) & (lam < WINDOW[1])].reset_index(drop=True)
    frac = boltzmann_fractions_from_levels(levels, T_SHELL)
    pop = frac[win["Lower"].to_numpy()]
    g_l = statistical_weight(win["J_Lower"].to_numpy())
    f_lu = 10 ** win["Log(gf)"].to_numpy() / g_l
    # n_ion per unit rho: X * rho / (A m_p) -> tau per unit rho:
    tau_per_rho = (
        SIGMA_CLASSICAL * f_lu * pop
        * win["WV_Transition"].to_numpy() * 1e-8 * T_EXP
        * x_frac / (a * M_P)
    )
    per_species.append(
        dict(z=z, a=a, x=x_frac, levels=levels, win=win, pop=pop,
             f_lu=f_lu, tau_per_rho=tau_per_rho)
    )
    print(f"Z={z}: {len(levels)} levels, {len(win)} window lines, "
          f"max tau/rho = {tau_per_rho.max():.3e}")

rho = TAU_MAX_TARGET / max(s["tau_per_rho"].max() for s in per_species)
print(f"rho = {rho:.3e} g/cm^3")
for s in per_species:
    taus = s["tau_per_rho"] * rho
    s["taus"] = taus
    print(f"  Z={s['z']}: tau>1: {(taus>1).sum():3d}   0.1-1: "
          f"{((taus>0.1)&(taus<=1)).sum():3d}   <=0.1: {(taus<=0.1).sum():4d}")

# --- SEDONA atom file: one group per element ---
with h5py.File(here / "atom_multiion.hdf5", "w") as f:
    for s in per_species:
        n_lev = len(s["levels"])
        g = f.create_group(str(s["z"]))
        g.attrs["n_ions"] = np.int64(2)
        g.attrs["n_levels"] = np.int64(n_lev + 1)
        g.attrs["n_lines"] = np.int64(len(s["win"]))
        g.create_dataset("ion_chi", data=np.array([9.9999e4, 9.9999e4]))
        g.create_dataset("ion_ground", data=np.array([0, n_lev], dtype=np.int64))
        g.create_dataset(
            "level_E",
            data=np.concatenate([s["levels"]["Energy"].to_numpy() * CM1_TO_EV, [0.0]]),
        )
        g.create_dataset(
            "level_g",
            data=np.concatenate(
                [statistical_weight(s["levels"]["J"].to_numpy()).astype(np.int64), [1]]
            ),
        )
        g.create_dataset("level_i", data=np.array([0] * n_lev + [1], dtype=np.int64))
        g.create_dataset("line_A", data=s["win"]["A"].to_numpy())
        g.create_dataset("line_l", data=s["win"]["Lower"].to_numpy().astype(np.int64))
        g.create_dataset("line_u", data=s["win"]["Upper"].to_numpy().astype(np.int64))

# --- model: two elements ---
v_edges = np.linspace(V_CORE, V_MAX, N_ZONES + 1)[1:]
with open(here / "multiion.mod", "w") as fh:
    fh.write("1D_sphere standard\n")
    fh.write(f"{N_ZONES}\t{V_CORE * T_EXP:.6e}\t{T_EXP:.6e} 2 \n")
    fh.write("57.139 58.140\n")
    for v in v_edges:
        fh.write(
            f"{v * T_EXP:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} 0.5 0.5\n"
        )

# --- line list for the Python solver: per-line pop_frac folds in the
#     species' ion density relative to the reference density rho/m_p ---
lam_out, f_out, pf_out, tau_out, z_out = [], [], [], [], []
for s in per_species:
    n_species_per_rho = s["x"] / (s["a"] * M_P)  # ions per gram
    for lam, f_lu, pop, tau in zip(
        s["win"]["WV_Transition"], s["f_lu"], s["pop"], s["taus"]
    ):
        lam_out.append(lam)
        f_out.append(f_lu)
        pf_out.append(pop * n_species_per_rho)  # multiply by rho later
        tau_out.append(tau)
        z_out.append(s["z"])

np.savez(
    here / "multiion_lines.npz",
    lam=np.array(lam_out), f_lu=np.array(f_out),
    popfrac_per_rho=np.array(pf_out), tau=np.array(tau_out),
    z=np.array(z_out), rho=rho,
)
print(f"total window lines: {len(lam_out)}")

# cross-species crowding: nearest-neighbour spacing of tau > 0.1 lines
lam_strong = np.sort(np.array(lam_out)[np.array(tau_out) > 0.1])
if len(lam_strong) > 1:
    dv = 2.99792458e10 * np.diff(lam_strong) / lam_strong[:-1] / 1e5
    print(f"strong lines (tau>0.1): {len(lam_strong)}, "
          f"min spacing = {dv.min():.0f} km/s, median = {np.median(dv):.0f} km/s")
