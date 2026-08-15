"""Generate the minimal 1-element / 1-ion / 1-line SEDONA model (Week 2).

Design (babystep_plan.md section 19):
- SEDONA's shipped 2-level "hydrogen" atom: one line at hc/10.2 eV = 1215.5 A,
  A = 1e9 s^-1  ->  f_lu = A lambda^2 g_u / (6.6702e15 g_l) = 0.6647.
- Shell at T = 2000 K: cold enough that (a) Saha keeps H neutral even at the
  very low density below (at 5000 K it would be FULLY ionized: the Saha RHS
  ~ 3e7 cm^-3 dwarfs n_H ~ 5), and (b) B_nu(T_shell) ~ 0 at 1215 A, so the
  line is pure absorption against the core -- exactly the regime the Python
  formal solver treats.
- Density chosen for tau_Sobolev ~ 2 at t = 20 days.
"""

import numpy as np

M_P = 1.67262192e-24
SIGMA_CLASSICAL = 0.026540083433884684

T_EXP = 20 * 86400.0  # 20 days, s
V_CORE = 1.0e8  # cm/s -- core boundary
V_MAX = 3.0e8  # cm/s -- outer edge
T_SHELL = 2000.0
N_ZONES = 100
F_LU = 0.6647
LAMBDA0_CM = 12398.42 / 10.2 * 1e-8  # hc/E in Angstrom -> cm

TAU_TARGET = 2.0
n_h = TAU_TARGET / (SIGMA_CLASSICAL * F_LU * LAMBDA0_CM * T_EXP)
rho = n_h * M_P

v_edges = np.linspace(V_CORE, V_MAX, N_ZONES + 1)[1:]
r_edges = v_edges * T_EXP

with open("minimal_1line.mod", "w") as fh:
    fh.write("1D_sphere standard\n")
    fh.write(f"{N_ZONES}\t{V_CORE * T_EXP:.6e}\t{T_EXP:.6e} 1 \n")
    fh.write("1.1\n")
    for r, v in zip(r_edges, v_edges):
        fh.write(f"{r:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} 1.0\n")

print(f"n_H = {n_h:.4f} cm^-3, rho = {rho:.4e} g/cm^3")
print(f"tau_S(target) = {TAU_TARGET}, lambda0 = {LAMBDA0_CM*1e8:.3f} A")
print(f"r_core = {V_CORE*T_EXP:.4e}, r_out = {V_MAX*T_EXP:.4e} cm")
