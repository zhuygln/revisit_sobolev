"""Is the closed-form resolved leg right? Brute force on IDENTICAL rays.

`resolved_attenuation` (erf) must agree with `emergent_luminosity` in
matching mode to <= 2e-4 -- it passed that on a 3-line toy at 300 km/s, but
the forest at 100 km/s gave erf 0.3405 against the cached solver_cold 0.3439
(default rays, worldline mode). Settle it: forest, same RaySet, every mode.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference import *  # noqa

lines, n_ion = forest_lines(stim=True)
nu = nu_grid(1600)
rs = RaySet.midpoint(R_CORE, R_OUT, 200, n_env=0)
L = legs_erf(lines, n_ion, rs, nu=nu)
print(f"erf classical {L['res']:.5f}   erf first {L['res_first']:.5f}   (v_D = 100 km/s, 200 core rays)\n")
print(" solver mode            F_res      vs erf-first   vs erf-classical   [s]")
for mode, dil in (("first", None), ("exact", None), ("worldline", False), ("worldline", True)):
    t0 = time.time()
    f = leg_solver(lines, n_ion, rs, relativity=mode, dilution=dil, nu=nu)
    tag = mode + ("" if dil is None else f"(dilution={dil})")
    print(f"  {tag:22s} {f:.5f}    {100*delta(f, L['res_first']):+7.3f}%      {100*delta(f, L['res']):+7.3f}%       [{time.time()-t0:.0f}]", flush=True)
# and the legacy call that produced forest_table.json's solver_cold = 0.3439
t0 = time.time()
from sobolev.formal_transfer import emergent_luminosity
lum = emergent_luminosity(nu, lines, const(n_ion), const(0.0), T_EXP, R_CORE, R_OUT, T_CORE, V_D, n_impact=150)
cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
print(f"  legacy(n_impact=150,worldline) {band_average(C/nu*1e8, lum/cont, BAND):.5f}   [{time.time()-t0:.0f}]")
