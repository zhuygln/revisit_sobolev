"""Generate the 2-line and 20-line SEDONA atoms + shared model (Week 2 ladder).

Atoms mirror the structure of pubsed's 2level_atomdata.hdf5: a hydrogen-like
element (Z=1) whose ion 0 has a ground level plus N excited levels, each with
one E1 line to ground (same A = 1e9 s^-1), and a bare ion stage. Line k sits
at E_k = 10.2 eV * (1 + k dv/c): equally spaced in velocity.

- 2-line atom:  dv = 1500 km/s -- troughs overlap (separation < v_max = 3000
  km/s), so the blend tests multiplicative attenuation exp(-(tau_1 + tau_2)).
- 20-line atom: dv = 750 km/s -- a mini forest spanning 15,000 km/s; a photon
  crosses ~4 resonances, the regime where expansion opacity is supposed to be
  a statistical approximation.

The shell model is the minimal_1line one at lower density: tau_S = 0.5 per
line, so the 4-line blend floor exp(-2) stays measurable above MC noise and
the per-line expansion-opacity error (1 - e^-tau vs tau) is ~20%.
"""

import h5py
import numpy as np

M_P = 1.67262192e-24
C = 2.99792458e10
SIGMA_CLASSICAL = 0.026540083433884684
HC_EV_A = 12398.42

T_EXP = 20 * 86400.0
V_CORE = 1.0e8
V_MAX = 3.0e8
T_SHELL = 2000.0
N_ZONES = 100
E0_EV = 10.2
A_UL = 1.0e9
G_L, G_U = 2, 6
TAU_TARGET = 0.5


def f_lu(lambda_angstrom):
    return A_UL * lambda_angstrom**2 * G_U / (6.6702e15 * G_L)


def write_atom(path, dv_cms, n_lines):
    energies = E0_EV * (1.0 + np.arange(n_lines) * dv_cms / C)
    n_lev = n_lines + 1  # ground + excited
    with h5py.File(path, "w") as f:
        g = f.create_group("1")
        # SEDONA sizes its arrays from these group attributes; without them
        # the reader segfaults.
        g.attrs["n_ions"] = np.int64(2)
        g.attrs["n_levels"] = np.int64(n_lev + 1)
        g.attrs["n_lines"] = np.int64(n_lines)
        g.create_dataset("ion_chi", data=np.array([13.6, 9.9999e4]))
        g.create_dataset("ion_ground", data=np.array([0, n_lev], dtype=np.int64))
        g.create_dataset(
            "level_E", data=np.concatenate([[0.0], energies, [0.0]])
        )
        g.create_dataset(
            "level_g",
            data=np.array([G_L] + [G_U] * n_lines + [1], dtype=np.int64),
        )
        g.create_dataset(
            "level_i", data=np.array([0] * n_lev + [1], dtype=np.int64)
        )
        g.create_dataset("line_A", data=np.full(n_lines, A_UL))
        g.create_dataset(
            "line_l", data=np.zeros(n_lines, dtype=np.int64)
        )
        g.create_dataset(
            "line_u", data=np.arange(1, n_lines + 1, dtype=np.int64)
        )
    lam = HC_EV_A / energies
    return lam, np.array([f_lu(l) for l in lam])


lam2, f2 = write_atom("atom_2line.hdf5", 1.5e8, 2)
lam20, f20 = write_atom("atom_20line.hdf5", 7.5e7, 20)

# Shared model: density for tau_S = TAU_TARGET on the reddest (10.2 eV) line.
lam0_cm = HC_EV_A / E0_EV * 1e-8
n_h = TAU_TARGET / (SIGMA_CLASSICAL * f_lu(HC_EV_A / E0_EV) * lam0_cm * T_EXP)
rho = n_h * M_P

v_edges = np.linspace(V_CORE, V_MAX, N_ZONES + 1)[1:]
with open("ladder.mod", "w") as fh:
    fh.write("1D_sphere standard\n")
    fh.write(f"{N_ZONES}\t{V_CORE * T_EXP:.6e}\t{T_EXP:.6e} 1 \n")
    fh.write("1.1\n")
    for v in v_edges:
        fh.write(f"{v * T_EXP:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} 1.0\n")

np.savez("ladder_lines.npz", lam2=lam2, f2=f2, lam20=lam20, f20=f20, n_h=n_h)
print(f"n_H = {n_h:.4f} cm^-3 (tau_S = {TAU_TARGET} per line)")
print(f"2-line: {lam2.round(2)} A")
print(f"20-line: {lam20[0]:.2f} .. {lam20[-1]:.2f} A")
