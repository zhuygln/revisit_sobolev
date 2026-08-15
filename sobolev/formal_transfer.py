"""Deterministic formal radiative transfer through a 1-D homologous sphere.

The independent reference solver of babystep_plan.md section 13 / research
requirements section 13: solve

    dI_nu/ds = -alpha_nu I_nu + j_nu,      S_nu = j_nu/alpha_nu = B_nu(T)

along parallel rays (impact parameter p, path coordinate z toward the
observer) with pure-absorption LTE line opacity. In homologous flow the
line-of-sight velocity is simply v_z = z/t for every ray -- resonance
surfaces are planes of constant z -- and the comoving frequency to first
order in v/c is nu' = nu (1 - z/(c t)).

Geometry: an emitting core (blackbody at T_core, radius r_core) inside a
line-forming shell extending to r_out. Rays with p < r_core start at the
core surface with I = B_nu(T_core); rays with p >= r_core traverse the
whole shell starting with I = 0.

Deliberately NOT here: electron scattering, non-LTE source functions,
relativistic terms beyond first order, time dependence. This solver's job
is to be small enough to trust, not complete.

KNOWN FIRST-ORDER LIMITATION: this observer-frame integration and the
comoving-frame Sobolev formula tau_S = sigma f n lambda0 t disagree at
O(v_bulk/c), where v_bulk is the flow velocity at the resonance point --
the observer-frame sweep rate is nu/(ct) while the comoving derivation
uses nu0/(ct), and the path-length transformation adds another such
factor. Measured empirically: a resonance at v/c = 7.7% produced a 7%
shallower trough than exp(-tau_S). Comparisons against Sobolev
predictions are only clean when the resonance lies at v_bulk << c, or
when both sides adopt the same frame convention. This is a systematic
that the Week 3+ SEDONA comparisons must control for explicitly.
"""

import numpy as np

from .constants import C, H, K_B, SIGMA_CLASSICAL
from .profiles import doppler_width_hz, gaussian, voigt


def planck_bnu(nu, temperature):
    """B_nu(T) in erg s^-1 cm^-2 Hz^-1 sr^-1.

    Deep in the Wien tail (h nu >> kT) expm1 overflows and the division
    correctly returns 0; the overflow itself is benign and suppressed.
    """
    nu = np.asarray(nu, dtype=float)
    x = H * nu / (K_B * temperature)
    with np.errstate(over="ignore"):
        return 2.0 * H * nu**3 / C**2 / np.expm1(x)


def emergent_luminosity(
    nu_grid,
    lines,
    n_l_of_r,
    temp_of_r,
    t_exp,
    r_core,
    r_out,
    t_core,
    v_doppler,
    n_impact=100,
    n_z_per_doppler=8.0,
    cutoff_widths=None,
):
    """Emergent line spectrum L_nu (erg s^-1 Hz^-1) of the core+shell system.

    Parameters
    ----------
    nu_grid : observer-frame frequencies (Hz)
    lines : sequence of (nu0, f_osc), (nu0, f_osc, pop_frac) or
        (nu0, f_osc, pop_frac, gamma) tuples -- one entry per transition.
        pop_frac (default 1) multiplies n_l_of_r for that line: for a real
        forest it is the Boltzmann fraction of the line's own lower level,
        with n_l_of_r the total ion density. gamma (Hz, default 0) is the
        Lorentzian FWHM of a Voigt profile; for natural broadening
        gamma = A_ul / (2 pi), which reproduces SEDONA's damping parameter
        a = A_ul / (4 pi dnu_D). gamma = 0 keeps the pure Gaussian.
    n_l_of_r : callable r -> lower-level number density (cm^-3)
    temp_of_r : callable r -> gas temperature (K), used for the LTE source
    t_exp : ejecta age (s); v = r / t_exp
    r_core, r_out : core and outer shell radii (cm)
    t_core : core blackbody temperature (K)
    v_doppler : intrinsic Doppler width (cm/s)
    n_impact : number of impact-parameter rays
    n_z_per_doppler : z-resolution in units of the Doppler length v_D * t_exp.
        The resonance region of every line must be resolved; unresolved grids
        underestimate tau exactly like an unconverged tau_exact.
    cutoff_widths : if not None, each line's profile is zeroed beyond this
        many Doppler widths from its centre -- mimicking SEDONA's hard-coded
        +-5-width truncation (AtomicSpecies_opacities.cpp). None (default)
        evaluates profiles everywhere.

    Returns L_nu on nu_grid. Luminosity from intensity: L = 8 pi^2 int I p dp.
    """
    nu_grid = np.asarray(nu_grid, dtype=float)
    # Normalize entries to (nu0, f_osc, pop_frac, gamma) and precompute widths.
    lines = [
        (
            entry[0],
            entry[1],
            entry[2] if len(entry) > 2 else 1.0,
            entry[3] if len(entry) > 3 else 0.0,
        )
        for entry in lines
    ]
    dnu_d = {nu0: doppler_width_hz(nu0, v_doppler) for nu0, _, _, _ in lines}

    # z step from the Doppler length; the same uniform grid serves all rays.
    dz = v_doppler * t_exp / n_z_per_doppler
    z_max = np.sqrt(max(r_out**2 - 0.0, 0.0))

    # Impact parameters: p-midpoint rule on [0, r_out]. Concentrate half the
    # rays on the core (p < r_core) since they carry the absorption trough.
    p_core = np.linspace(0.0, r_core, n_impact // 2, endpoint=False)
    p_core += 0.5 * (p_core[1] - p_core[0])
    p_env = np.linspace(r_core, r_out, n_impact - n_impact // 2, endpoint=False)
    p_env += 0.5 * (p_env[1] - p_env[0])
    p_all = np.concatenate([p_core, p_env])

    intensity = np.zeros((p_all.size, nu_grid.size))
    for i, p in enumerate(p_all):
        if p < r_core:
            z0 = np.sqrt(r_core**2 - p**2)
            i0 = planck_bnu(nu_grid, t_core)
        else:
            z0 = -np.sqrt(r_out**2 - p**2)
            i0 = np.zeros_like(nu_grid)
        z1 = np.sqrt(r_out**2 - p**2)
        n_z = max(int(np.ceil((z1 - z0) / dz)), 4)
        z = np.linspace(z0, z1, n_z)
        r = np.sqrt(p**2 + z**2)
        n_l = n_l_of_r(r)
        temp = temp_of_r(r)

        # alpha(z, nu): sum the resolved profiles of all lines.
        nu_com = nu_grid[None, :] * (1.0 - z[:, None] / (C * t_exp))
        alpha = np.zeros((n_z, nu_grid.size))
        for nu0, f_osc, pop, gamma in lines:
            dnu = nu_com - nu0
            if gamma > 0.0:
                phi = voigt(dnu, dnu_d[nu0], gamma)
            else:
                phi = gaussian(dnu, dnu_d[nu0])
            if cutoff_widths is not None:
                phi = np.where(
                    np.abs(dnu) <= cutoff_widths * dnu_d[nu0], phi, 0.0
                )
            alpha += SIGMA_CLASSICAL * f_osc * pop * n_l[:, None] * phi

        # Formal solution via cumulative optical depth (trapezoidal steps):
        #   I_out = I_0 e^{-tau_tot} + sum_k S_k (e^{-t_above_k+1} - e^{-t_above_k})
        dtau = 0.5 * (alpha[1:, :] + alpha[:-1, :]) * np.diff(z)[:, None]
        tau_below = np.vstack([np.zeros_like(nu_grid), np.cumsum(dtau, axis=0)])
        tau_tot = tau_below[-1, :]
        tau_above = tau_tot[None, :] - tau_below
        # Source at segment midpoints, at the local comoving frequency.
        s_mid = 0.5 * (
            planck_bnu(nu_com[1:, :], 0.5 * (temp[1:] + temp[:-1])[:, None])
            + planck_bnu(nu_com[:-1, :], 0.5 * (temp[1:] + temp[:-1])[:, None])
        )
        emis = np.sum(
            s_mid * (np.exp(-tau_above[1:, :]) - np.exp(-tau_above[:-1, :])), axis=0
        )
        intensity[i, :] = i0 * np.exp(-tau_tot) + emis

    # L_nu = 8 pi^2 int I(p) p dp, midpoint rule with the two p sections.
    dp_core = p_core[1] - p_core[0] if p_core.size > 1 else r_core
    dp_env = p_env[1] - p_env[0] if p_env.size > 1 else r_out - r_core
    weights = np.concatenate(
        [p_core * dp_core, p_env * dp_env]
    )
    return 8.0 * np.pi**2 * np.sum(weights[:, None] * intensity, axis=0)
