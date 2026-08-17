"""Why does the Sobolev residual grow with v_D? A boundary experiment.

Hypothesis: the Sobolev staircase applies a hard crossed/not-crossed STEP at
the shell edges, while the resolved calculation has a profile of finite width
that is CLIPPED there. The disagreement is confined to within a few Doppler
widths of a boundary, so its band-averaged size scales as v_D / (velocity
span of the shell) -- which is exactly the observed scaling.

For a central ray the prediction is analytic. With the profile centred at
z_res and the shell running z_lo..z_hi,

    tau/tau_S = 1/2 [ erf((z_hi - z_res)/Delta) - erf((z_lo - z_res)/Delta) ],
    Delta = v_D t,

going to 1 deep inside, 1/2 exactly at an edge, and 0 outside -- against the
Sobolev step of 1 inside and 0 outside. No second line is involved: this is a
single-line, isolation-irrelevant effect.

Also emits the null overlap test: two identical lines far from any boundary,
separation swept through Delta v / v_D = 20 ... 0. In pure absorption with
fixed populations the optical depths add exactly, so the prediction is
Delta_overlap = 0 at every separation.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.special import erf
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C, SIGMA_CLASSICAL

T_EXP = 86400.0
V_CORE, V_MAX = 1.0e8, 3.0e8           # the shell used throughout the project
Z_LO, Z_HI = V_CORE * T_EXP, V_MAX * T_EXP
LAMBDA0 = 4000e-8
NU0 = C / LAMBDA0
F_OSC, N_L = 0.5, 20.0
TAU_S = SIGMA_CLASSICAL * F_OSC * N_L * LAMBDA0 * T_EXP


def tau_resolved(z_res, v_d, n=400001, half_widths=60.0):
    """Direct integration of the resolved profile, clipped by the shell."""
    delta = v_d * T_EXP
    lo = max(Z_LO, z_res - half_widths * delta)
    hi = min(Z_HI, z_res + half_widths * delta)
    if hi <= lo:
        return 0.0
    z = np.linspace(lo, hi, n)
    dnu_d = NU0 * v_d / C
    # comoving offset for a photon resonant at z_res
    dnu = NU0 * (z_res - z) / (C * T_EXP)
    phi = np.exp(-((dnu / dnu_d) ** 2)) / (np.sqrt(np.pi) * dnu_d)
    # int(phi dz) = c t / nu0, so int(alpha dz) = sigma f n lambda0 t = tau_S
    return np.trapezoid(SIGMA_CLASSICAL * F_OSC * N_L * phi, z)


print("=== single-line boundary sweep ===")
print("Resonance walked out through the OUTER edge; d = distance in Doppler widths")
print(" v_D    d/v_D   tau_res/tau_S   erf prediction   Sobolev step")
rows = []
for v_d in (10.0e5, 100.0e5, 300.0e5):
    delta = v_d * T_EXP
    for d_over in (-6, -3, -1, 0, 1, 3, 6):
        z_res = Z_HI + d_over * delta          # negative = inside the shell
        got = tau_resolved(z_res, v_d) / TAU_S
        pred = 0.5 * (erf((Z_HI - z_res) / delta) - erf((Z_LO - z_res) / delta))
        step = 1.0 if Z_LO < z_res < Z_HI else 0.0
        rows.append((v_d / 1e5, d_over, got, pred, step))
        print(f"{v_d/1e5:5.0f} {d_over:8.0f} {got:14.4f} {pred:16.4f} {step:13.0f}")

print("\n=== null overlap test: two identical lines, far from boundaries ===")
print("Resolved transmission vs Sobolev exp(-sum tau_S) at the blended centre")
print("  dv/v_D   resolved      Sobolev      ratio")
v_d = 100.0e5
delta = v_d * T_EXP
z1 = 0.5 * (Z_LO + Z_HI)
for sep in (20, 10, 5, 2, 1, 0.5, 0.0):
    z2 = z1 + sep * delta
    if not (Z_LO < z2 < Z_HI):
        continue
    # A photon whose resonance for line 1 sits at z1 also samples line 2's
    # profile; integrate BOTH profiles along the ray and exponentiate.
    lo, hi = Z_LO, Z_HI
    z = np.linspace(lo, hi, 600001)
    dnu_d = NU0 * v_d / C
    tot = np.zeros_like(z)
    for zc in (z1, z2):
        dnu = NU0 * (zc - z) / (C * T_EXP)
        tot += np.exp(-((dnu / dnu_d) ** 2)) / (np.sqrt(np.pi) * dnu_d)
    tau_res = np.trapezoid(SIGMA_CLASSICAL * F_OSC * N_L * tot, z)
    print(f"{sep:8.1f} {np.exp(-tau_res):10.6f} {np.exp(-2*TAU_S):12.6f}"
          f" {np.exp(-tau_res)/np.exp(-2*TAU_S):10.6f}")

# ---------------- figure ----------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
d_fine = np.linspace(-8, 8, 400)
for v_d, col in ((10.0e5, "C0"), (100.0e5, "C1"), (300.0e5, "C3")):
    delta = v_d * T_EXP
    got = [tau_resolved(Z_HI + x * delta, v_d) / TAU_S for x in d_fine]
    ax1.plot(d_fine, got, color=col, label=f"$v_D$ = {v_d/1e5:.0f} km/s")
ax1.plot(d_fine, 0.5 * (1 - erf(d_fine)), "k--", lw=1, label="erf prediction")
ax1.step([-8, 0, 0, 8], [1, 1, 0, 0], where="post", color="gray", lw=1.2,
         label="Sobolev step")
ax1.set_xlabel(r"distance past the outer edge  $d/v_D$")
ax1.set_ylabel(r"$\tau_{\rm resolved}/\tau_S$")
ax1.set_title("A single line at a shell boundary")
ax1.legend(fontsize=7)
ax1.grid(alpha=0.3)

# band-averaged consequence: fraction of the shell within a few widths of an edge
vds = np.geomspace(1, 300, 40)
frac = 2 * vds / ((V_MAX - V_CORE) / 1e5)
meas = {1: 0.3, 3: 0.3, 10: 0.2, 30: 0.5, 100: 2.5, 300: 9.2}
ax2.loglog(vds, 100 * frac, "k--", lw=1.2,
           label=r"$\propto v_D/\Delta v_{\rm shell}$")
ax2.loglog(list(meas), [abs(v) for v in meas.values()], "o", color="C1",
           label=r"measured $|\Delta_{\rm Sobolev}|$")
ax2.set_xlabel(r"$v_D$ [km/s]")
ax2.set_ylabel(r"$|\Delta_{\rm Sobolev}|$  [%]")
ax2.set_title("Band-averaged effect scales with the edge fraction")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3, which="both")

fig.tight_layout()
out = ROOT / "outputs" / "fig12_boundary.png"
fig.savefig(out, dpi=200)
print(f"\nsaved {out}")
