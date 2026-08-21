"""The deterministic reference for Delta_Sobolev (referee Comment 5).

Delta_Sob used to compare an analytic Sobolev leg against SEDONA's resolved
Monte Carlo, which exposed it to three conventions that had to be matched by
hand -- thermal emission, normalization, transport treatment -- and each had
already produced a spurious few-percent result. The referee's point stands:
the cleanest quantity is the finite-profile vs delta-resonance difference on
IDENTICAL rays with the same source convention, populations and transport
law, with SEDONA as the independent validation of the reference rather than
the reference itself.

This module builds that. Four matched pairs, each a same-code differential:

    (Sob None,      erf classical)   the classical plane, no nu0/nu factor
    (Sob first,     erf first)       first-order Doppler, the nu0/nu factor kept
    (Sob exact,     solver exact)    frozen snapshot, exact SR  [brute force]
    (Sob worldline, solver worldline + dilution)  the physical law [brute force]

The erf pairs use `sobolev_leg.resolved_attenuation`: closed form, seconds,
and independent of v_D, so the 1 km/s frontier is free. The solver pairs use
`formal_transfer.emergent_luminosity` on the same RaySet and are the check
that the closed form is right.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.formal_transfer import emergent_luminosity, planck_bnu
from sobolev.optical_depth import stimulated_emission_factor
from sobolev.rays import RaySet
from sobolev.sobolev_leg import expansion_damp, resolved_attenuation, sobolev_attenuation
from sobolev.spectra import band_average

FOREST = ROOT / "experiments/laII_forest"
T_EXP, R_CORE, R_OUT, T_CORE, T_SHELL = 86400.0, 8.64e12, 2.592e13, 6000.0, 3000.0
V_D = 1.0e7
BAND = (3800.0, 3955.0)


def forest_lines(stim=True, t_shell=T_SHELL):
    """The La II forest as (nu0, f, pop) tuples; `stim` folds the LTE
    stimulated-emission factor into pop (what SEDONA does)."""
    d = np.load(FOREST / "forest_lines.npz")
    lines = []
    for lam, f, p in zip(d["lam"], d["f_lu"], d["pop"]):
        nu0 = C / (lam * 1e-8)
        if stim:
            p = p * stimulated_emission_factor(nu0, t_shell)
        lines.append((nu0, f, p))
    return lines, float(d["n_ion"])


def const(v):
    return lambda r: np.full_like(np.asarray(r, dtype=float), v)


def nu_grid(n=1600):
    return np.geomspace(7.50e14, 7.95e14, n)


def legs_erf(lines, n_ion, rays, v_d=V_D, nu=None, band=BAND):
    """Band-averaged {sob, sob_first, res, res_first, exp} on shared rays."""
    nu = nu_grid() if nu is None else nu
    lam = C / nu * 1e8
    kw = dict(r_core=R_CORE, r_out=R_OUT, t_exp=T_EXP, n_ref=n_ion, rays=rays)
    out = {
        "sob": sobolev_attenuation(nu, lines, **kw),
        "sob_first": sobolev_attenuation(nu, lines, relativity="first", **kw),
        "exp": sobolev_attenuation(nu, lines, damp=expansion_damp, **kw),
        "res": resolved_attenuation(nu, lines, v_doppler=v_d, sweep="classical", **kw),
        "res_first": resolved_attenuation(nu, lines, v_doppler=v_d, sweep="first", **kw),
    }
    return {k: band_average(lam, v, band) for k, v in out.items()}


def leg_solver(lines, n_ion, rays, relativity, dilution=None, v_d=V_D, nu=None,
               t_shell=0.0, cutoff_widths=None, band=BAND):
    """Brute-force resolved leg on the same rays, T_shell -> 0 (F8)."""
    nu = nu_grid() if nu is None else nu
    lam = C / nu * 1e8
    lum = emergent_luminosity(
        nu, lines, const(n_ion), const(t_shell), T_EXP, R_CORE, R_OUT, T_CORE, v_d,
        relativity=relativity, dilution=dilution, rays=rays, cutoff_widths=cutoff_widths,
    )
    cont = 4.0 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)
    return band_average(lam, lum / cont, band)


def delta(a, b):
    return (a - b) / b
