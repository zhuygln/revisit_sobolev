"""PILOT: does the Sobolev/expansion separation survive at 0.1-0.3c?

NOT the Stage 6 measurement. A synthetic 20-line forest at uniform density,
one epoch, one tau_max -- enough to establish the sign and rough size of the
beta dependence and to expose the design traps, not enough to quote as a
result. The real version needs the La II forest and a tau_max range.

Two false starts are preserved as `--confounded`, because both produced
plausible tables that meant nothing:

  1. Scaling the shell but not the line list. At beta = 0.3 each line's
     resonance region covers 30x more velocity space than at 0.01, so the
     forest goes from sparse to fully blanketed and one measures F7's density
     axis instead of relativity. Delta_expansion appeared to COLLAPSE, +51% to
     +1.6%. Entirely an artifact.
  2. Fixing (1) but leaving v_D fixed. Then v_D/delta_v_shell varies 30x and
     the low-beta rows carry F12's finite-region boundary effect.

The self-similar family holds fixed every dimensionless ratio known to move
Delta -- tau_max (F10), lines per crossing (F7), v_D/delta_v_shell (F12), and
beta_core/beta_out -- and varies only beta. Anything left is relativity plus
light-travel dilution.

CONVERGENCE. The resolved leg is the one that has to be converged, and it is
not converged at the defaults. `emergent_luminosity` averages over
n_impact//2 core rays while `sobolev_attenuation` averages over n_p = 200, so
Delta drifts as the solver's p-sampling settles: Delta_Sobolev ran -0.93% ->
+0.59% over n_impact = 12 -> 96 and CHANGED SIGN. Frequency resolution, by
contrast, was already converged at n_nu = 200 (Delta_exp stable to +-0.12
points out to 3200), because the structure scale in the band is the shell
span, not the Doppler width.

The beta DIFFERENCE is far better behaved than either endpoint: the ray error
is beta-independent and cancels, so d(Delta_exp) held at +4.40/+4.57/+4.51/
+4.58 across an 8x ray refinement while the endpoints moved by 1.8 points.

Usage:
    python pilot.py                # the self-similar family, converged settings
    python pilot.py --converge     # the n_nu and n_impact ladders
    python pilot.py --control      # the tau_max-collapse test
    python pilot.py --confounded   # the two false starts, for the record

Runtime warning: the converged points are ~10-20 min each at n_impact = 96.
Matching the two legs' p-sampling would remove most of that cost and is the
first thing to fix in a committed Stage 6 experiment.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sobolev.constants import C, SIGMA_CLASSICAL
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.sobolev_leg import expansion_damp, sobolev_attenuation

T_EXP, T_CORE = 86400.0, 2.0e4
N_LINES, F_OSC, LAM_MID = 20, 0.5, 4000e-8
TAU = 5.0
# anchored so that beta_out = 0.3 keeps the project's fiducial v_D = 100 km/s
VD_OVER_SPAN = 1.0e7 / (0.2 * C)


def _const(v):
    return lambda r: np.full_like(np.asarray(r, dtype=float), v)


def legs(b_out, tau_max_target, n_nu=800, n_impact=96, self_similar=True,
         fixed_v_d=None, fixed_window=None):
    """Return (resolved, sobolev, expansion) band-averaged transmitted fractions.

    All three legs are worldline-consistent: the resolved leg runs
    relativity="worldline", dilution=True with T_shell -> 0 (the F8 emission
    convention), and both analytic legs run relativity="worldline", differing
    only by `damp`. Delta_expansion is therefore a same-code differential.
    """
    ct = C * T_EXP
    b_core = b_out / 3.0
    r_core, r_out = b_core * ct, b_out * ct
    d_beta = b_out - b_core

    v_d = fixed_v_d if fixed_v_d is not None else VD_OVER_SPAN * d_beta * C
    span = fixed_window if fixed_window is not None else d_beta
    if not self_similar:
        span = fixed_window if fixed_window is not None else 0.025

    lam0 = LAM_MID * (1.0 + np.linspace(-0.5, 0.5, N_LINES) * span)
    nu_lines = C / lam0
    lines = [(nu, F_OSC) for nu in nu_lines]
    n_ref = tau_max_target / (SIGMA_CLASSICAL * F_OSC * lam0.max() * T_EXP)

    nu = np.linspace(
        nu_lines.min(), nu_lines.max() * np.sqrt((1 + b_out) / (1 - b_out)), n_nu
    )
    width = nu[-1] - nu[0]

    lum = emergent_luminosity(
        nu, lines, _const(n_ref), _const(0.0), T_EXP, r_core, r_out, T_CORE,
        v_d, n_impact=n_impact, relativity="worldline", dilution=True,
    )
    cont = 4.0 * np.pi**2 * r_core**2 * planck_bnu(nu, T_CORE)
    resolved = np.trapezoid(lum / cont, nu) / width

    kw = dict(r_core=r_core, r_out=r_out, t_exp=T_EXP, n_ref=n_ref,
              relativity="worldline")
    sob = np.trapezoid(sobolev_attenuation(nu, lines, **kw), nu) / width
    exp_ = np.trapezoid(
        sobolev_attenuation(nu, lines, damp=expansion_damp, **kw), nu
    ) / width
    return resolved, sob, exp_


def _delta(res, sob, exp_):
    return 100 * (sob - res) / res, 100 * (exp_ - res) / res


def run_family(n_impact=96, n_nu=800):
    print("Self-similar family: tau_max, lines/crossing, v_D/span and")
    print("beta_core/beta_out all held fixed; only beta varies.")
    print(f"n_impact = {n_impact}, n_nu = {n_nu}, tau_max(snapshot) = {TAU}\n")
    print(" beta_out   v_D km/s   D_Sobolev   D_expansion   R=(1-b)^2/g")
    out = {}
    for b_out in (0.01, 0.30):
        t0 = time.time()
        d_s, d_e = _delta(*legs(b_out, TAU, n_nu=n_nu, n_impact=n_impact))
        g = 1.0 / np.sqrt(1.0 - b_out**2)
        v_d = VD_OVER_SPAN * (b_out - b_out / 3.0) * C
        out[b_out] = (d_s, d_e)
        print(f"  {b_out:5.2f}     {v_d/1e5:7.1f}    {d_s:+7.2f}%   {d_e:+7.2f}%"
              f"      {(1-b_out)**2/g:.4f}   [{time.time()-t0:.0f}s]")
        sys.stdout.flush()
    ds = out[0.30][0] - out[0.01][0]
    de = out[0.30][1] - out[0.01][1]
    print(f"\n  d(D_Sobolev)   = {ds:+.2f} points   -- no beta dependence")
    print(f"  d(D_expansion) = {de:+.2f} points   -- expansion opacity gets WORSE")
    return out


def run_converge():
    b = 0.30
    print("Frequency resolution (n_impact = 12):")
    print("  n_nu    D_Sobolev   D_expansion")
    for n_nu in (200, 400, 800, 1600, 3200):
        d_s, d_e = _delta(*legs(b, TAU, n_nu=n_nu, n_impact=12))
        print(f"  {n_nu:5d}   {d_s:+7.2f}%   {d_e:+7.2f}%")
        sys.stdout.flush()

    print("\nRay count, BOTH endpoints (n_nu = 800). The endpoints drift; the")
    print("difference does not, because the ray error is beta-independent.")
    print(" n_impact | b=0.01  D_Sob   D_exp  | b=0.30  D_Sob   D_exp  | dD_Sob  dD_exp")
    for n_imp in (12, 24, 48, 96):
        dl_s, dl_e = _delta(*legs(0.01, TAU, n_nu=800, n_impact=n_imp))
        dh_s, dh_e = _delta(*legs(0.30, TAU, n_nu=800, n_impact=n_imp))
        print(f"   {n_imp:4d}   |        {dl_s:+6.2f}% {dl_e:+7.2f}% |"
              f"        {dh_s:+6.2f}% {dh_e:+7.2f}% | {dh_s-dl_s:+6.2f}  {dh_e-dl_e:+6.2f}")
        sys.stdout.flush()


def run_control(n_impact=96):
    """Delta(tau_max) collapse: does high beta just move you along the curve?"""
    b = 0.30
    g = 1.0 / np.sqrt(1.0 - b**2)
    R = (1.0 - b) ** 2 / g
    d_ctrl = _delta(*legs(0.01, TAU * R, n_nu=800, n_impact=n_impact))[1]
    d_hi = _delta(*legs(b, TAU, n_nu=800, n_impact=n_impact))[1]
    print(f"beta=0.01 at tau_max = {TAU*R:.3f} (= {TAU} x R, R = {R:.4f}):"
          f"  D_exp = {d_ctrl:+.2f}%")
    print(f"beta=0.30 at tau_max = {TAU:.3f}:"
          f"                          D_exp = {d_hi:+.2f}%")
    print(f"ratio = {d_hi/d_ctrl:.3f}   (the collapse hypothesis predicts 1.000)")
    print("\nRefuted. Only about two thirds of the high-beta growth is the")
    print("tau_max shift; the rest is an unidentified mechanism.")


def run_confounded():
    print("The two false starts, kept because both looked plausible.\n")
    print("(1) shell scaled, line list NOT -- measures F7's blanketing axis:")
    print("  beta_out   D_expansion")
    for b_out in (0.01, 0.1, 0.2, 0.3):
        d_e = _delta(*legs(b_out, TAU, n_nu=200, n_impact=12,
                           self_similar=False, fixed_window=0.025))[1]
        print(f"    {b_out:5.2f}     {d_e:+7.2f}%     <- apparent collapse, an artifact")
        sys.stdout.flush()
    print("\n(2) window scaled but v_D fixed -- low beta carries F12's boundary effect:")
    print("  beta_out   D_Sobolev   D_expansion")
    for b_out in (0.01, 0.3):
        d_s, d_e = _delta(*legs(b_out, TAU, n_nu=200, n_impact=12, fixed_v_d=1.0e7))
        print(f"    {b_out:5.2f}    {d_s:+7.2f}%    {d_e:+7.2f}%")
        sys.stdout.flush()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--converge", action="store_true")
    ap.add_argument("--control", action="store_true")
    ap.add_argument("--confounded", action="store_true")
    ap.add_argument("--n-impact", type=int, default=96)
    a = ap.parse_args()
    if a.converge:
        run_converge()
    elif a.control:
        run_control(n_impact=a.n_impact)
    elif a.confounded:
        run_confounded()
    else:
        run_family(n_impact=a.n_impact)
