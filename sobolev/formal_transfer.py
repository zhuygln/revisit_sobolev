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

FRAME TREATMENT (`relativity` argument):

  "first"  -- nu' = nu (1 - z/ct), opacity used as sigma f n phi with NO
              Doppler factor. This is the historical behaviour of this
              module and the default, so that previously reported numbers
              are reproducible. It is NOT the correct relativistic result.
  "exact"  -- full special relativity: nu' = gamma nu (1 - beta_z) with
              beta = r/(ct) and beta_z = z/(ct); the comoving opacity is
              transformed to the lab frame as chi_lab = D chi', with the
              Doppler factor D = nu'/nu (from the invariance of nu*chi_nu);
              and the source function as S_lab = S'(nu')/D^3 (from the
              invariance of I_nu/nu^3).

The two modes agree at O(beta) and separate only at O(beta^2): the D factor
the "first" mode omits is cancelled, to leading order, by the sweep-rate
factor it also gets wrong (see sobolev_leg.tau_sobolev_relativistic for the
algebra). Both give an effective optical depth tau_S (1 - beta) at first
order. Consequently the systematic recorded as Finding F3 -- a resonance at
beta = 7.7% giving a 7% shallower trough than exp(-tau_S) -- is NOT an error
in this module: it is the correct relativistic behaviour, and it is the
nonrelativistic tau_S that is wrong at O(v/c).

Which mode is physically right is settled by experiment: see
experiments/vc_control/, which sweeps beta with both solver modes, SEDONA,
and every candidate Sobolev expression on one axis.
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
    relativity="first",
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
        # Frame transformation. beta_z is the line-of-sight component of the
        # homologous velocity; beta is its magnitude (both < 1 by construction
        # for a shell with v_max < c).
        beta_z = z / (C * t_exp)
        if relativity == "exact":
            beta = r / (C * t_exp)
            lorentz = 1.0 / np.sqrt(np.maximum(1.0 - beta**2, 1e-300))
            doppler = lorentz * (1.0 - beta_z)  # D = nu'/nu
        elif relativity == "first":
            doppler = 1.0 - beta_z
        else:
            raise ValueError(f"relativity must be 'first' or 'exact', got {relativity!r}")
        nu_com = nu_grid[None, :] * doppler[:, None]

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
            chi_comoving = SIGMA_CLASSICAL * f_osc * pop * n_l[:, None] * phi
            if relativity == "exact":
                # nu * chi_nu is a Lorentz invariant, so chi_lab = D chi'.
                alpha += doppler[:, None] * chi_comoving
            else:
                alpha += chi_comoving

        # Formal solution via cumulative optical depth (trapezoidal steps):
        #   I_out = I_0 e^{-tau_tot} + sum_k S_k (e^{-t_above_k+1} - e^{-t_above_k})
        dtau = 0.5 * (alpha[1:, :] + alpha[:-1, :]) * np.diff(z)[:, None]
        tau_below = np.vstack([np.zeros_like(nu_grid), np.cumsum(dtau, axis=0)])
        tau_tot = tau_below[-1, :]
        tau_above = tau_tot[None, :] - tau_below
        # Source at segment midpoints, at the local comoving frequency. In
        # exact mode the comoving source must be carried back to the lab
        # frame: I_nu/nu^3 is invariant, so S_lab(nu) = S'(nu')/D^3.
        t_mid = 0.5 * (temp[1:] + temp[:-1])[:, None]
        s_mid = 0.5 * (
            planck_bnu(nu_com[1:, :], t_mid) + planck_bnu(nu_com[:-1, :], t_mid)
        )
        if relativity == "exact":
            d_mid = 0.5 * (doppler[1:] + doppler[:-1])[:, None]
            s_mid = s_mid / d_mid**3
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
