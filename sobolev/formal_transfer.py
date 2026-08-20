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

FRAME TREATMENT (`relativity` argument). The important distinction here is
not first-order versus exact special relativity -- it is whether the ejecta
is treated as a FROZEN SNAPSHOT or as a genuinely time-dependent flow.

  "worldline" -- DEFAULT and the physical one. The photon takes time to
              travel, so t advances along the path: t(s) = t_exp + s/c,
              with s measured from the ray's starting point. In homologous
              flow beta = r/(c t) then changes because BOTH r and t change,
              giving db/dt = (1-b)/t rather than the frozen db/dz = 1/(ct).
  "exact"   -- frozen snapshot, exact special relativity. t is held fixed
              while integrating over z. Exact Doppler nu' = gamma nu (1-b),
              chi_lab = D chi' from the invariance of nu*chi_nu, and
              S_lab = S'/D^3 from the invariance of I_nu/nu^3.
  "first"   -- frozen snapshot, first order, no opacity transformation.
              The module's original behaviour; retained to reproduce
              previously published numbers.

For a single line of comoving depth tau_S the three give, relative to
tau_S evaluated with the local density and the local ejecta age at the
resonance,

    worldline : tau/tau_S = 1/gamma        = 1 - beta^2/2 + ...
    exact     : tau/tau_S = (1-beta)/gamma = 1 - beta + ...
    first     : tau/tau_S = (1-beta)          (to first order)

verified against a first-principles integration in tests/test_relativistic.py
to five decimals. The physical law therefore has NO first-order term: at
beta = 0.1/0.2/0.3 the correction is 0.5%/2.0%/4.6%, not 10%/20%/30%.

This supersedes two earlier readings of Finding F3. The offset measured at
beta = 7.7% is neither a frame ambiguity nor "the leading relativistic
correction": it is an artifact of integrating a frozen snapshot, which is
not the problem the photon solves.

THE MEDIUM (`dilution` argument). Worldline kinematics alone are not a
self-consistent flow, and the inconsistency is not small. With `dilution`
False the ejecta age advances along the ray while the density and the outer
boundary stay at their t_exp snapshot; with it True the medium co-evolves,
homology making beta the Lagrangian label, so the element now at
(r, t_local) is read off `n_l_of_r` at r * t_exp/t_local and diluted by
(t_exp/t_local)^3, and the outer boundary recedes as r_out(t) = beta_out c t.

This matters more than the relativity does. A photon leaving the origin at
t0 meets the material moving at beta at t_res = t0/(1-beta) -- later than the
frozen-geometry guess t0(1+beta), because the fluid element is itself moving
outward while the photon chases it. Normalized against tau_S at the emission
epoch, which is the epoch a spectrum is labelled by,

    frozen medium, worldline kinematics : tau/tau_S = 1/[(1-beta) gamma]
    self-consistent flow                : tau/tau_S = (1-beta)^2 / gamma

At beta = 0.3 that is 1.363 against 0.467 -- a factor 2.9, against the 4.6%
that 1/gamma contributes. The dilution term is roughly 11x the relativistic
one, and it points the other way. Both laws are confirmed against direct
integration of the resolved opacity in tests/test_relativistic.py.

A frozen outer boundary also caps the line-of-sight velocity a photon can
ever reach at beta_out/(1+beta_out) -- 0.23 for a 0.3c shell -- so high-beta
resonances simply go missing. `dilution` defaults to None, meaning "not
specified", and worldline mode raises above beta_out = 0.02 rather than
silently picking one; pass False explicitly to hold the medium fixed on
purpose, as the transport-law tests do.

Still not modelled: the core is a fixed-radius emitting surface (correct as
a starting condition, since the photon leaves it at t_exp, but it does not
co-move), `temp_of_r` is remapped to the fluid element without any adiabatic
cooling law, and there is no time-dependent transport of the source itself.
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
    relativity="worldline",
    dilution=None,
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
    dilution : None (default, "not specified"), True or False. True lets the
        medium co-evolve with the flow -- density diluted as t^-3 off the
        Lagrangian radius, outer boundary receding as beta_out c t. See the
        module docstring; at beta_out > 0.02 worldline mode requires an
        explicit choice rather than defaulting.
    cutoff_widths : if not None, each line's profile is zeroed beyond this
        many Doppler widths from its centre -- mimicking SEDONA's hard-coded
        +-5-width truncation (AtomicSpecies_opacities.cpp). None (default)
        evaluates profiles everywhere.

    Returns L_nu on nu_grid. Luminosity from intensity: L = 8 pi^2 int I p dp.
    """
    nu_grid = np.asarray(nu_grid, dtype=float)

    if relativity not in ("worldline", "exact", "first"):
        raise ValueError(
            f"relativity must be 'worldline', 'exact' or 'first', got {relativity!r}"
        )
    # dilution=None means "not specified"; the guard below fires on that but
    # not on an explicit False, so a caller can still hold the medium frozen
    # on purpose (to isolate the transport law) without tripping it.
    undecided = dilution is None
    dilution = bool(dilution)
    if dilution and relativity != "worldline":
        raise ValueError(
            "dilution=True is only meaningful with relativity='worldline'; the "
            f"frozen modes hold t fixed, so nothing dilutes (got {relativity!r})"
        )
    beta_out = r_out / (C * t_exp)
    if relativity == "worldline" and undecided and beta_out > 0.02:
        raise ValueError(
            f"worldline transport through a frozen medium at beta_out = {beta_out:.3g}. "
            "The ejecta age advances along the ray but the density and the outer "
            "boundary do not, which is not a small inconsistency: at beta = 0.3 the "
            "density is wrong by (1-beta)^3 = 0.34 and tau by (1-beta)^2 = 0.49, an "
            "error 11x the 1/gamma term this mode exists to capture. It also caps the "
            "reachable line-of-sight velocity at beta_out/(1+beta_out). Pass "
            "dilution=True for the self-consistent flow, relativity='exact' to "
            "integrate a frozen snapshot deliberately, or dilution=False to state "
            "that you want worldline kinematics on a frozen medium anyway."
        )

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
        if dilution:
            # The outer boundary is a fluid surface, so it recedes with the
            # flow: r_out(t) = beta_out c t. Exit where sqrt(p^2+z^2) =
            # beta_out (z + Z0), with c t(z) = z + Z0 along the worldline.
            big_z = C * t_exp - z0
            b_out = r_out / (C * t_exp)
            disc = b_out**2 * big_z**2 - (1.0 - b_out**2) * p**2
            z1 = (b_out**2 * big_z + np.sqrt(max(disc, 0.0))) / (1.0 - b_out**2)
        else:
            z1 = np.sqrt(r_out**2 - p**2)
        n_z = max(int(np.ceil((z1 - z0) / dz)), 4)
        z = np.linspace(z0, z1, n_z)
        r = np.sqrt(p**2 + z**2)

        # alpha(z, nu): sum the resolved profiles of all lines.
        # Frame transformation. beta_z is the line-of-sight component of the
        # homologous velocity; beta is its magnitude (both < 1 by construction
        # for a shell with v_max < c). In "worldline" mode the ejecta age
        # advances along the photon path, so beta = r/(c t) changes through
        # BOTH r and t -- that time term is the whole difference between the
        # frozen and physical laws, and it is O(beta), not O(beta^2).
        if relativity == "worldline":
            t_local = t_exp + (z - z[0]) / C
        else:  # "exact" / "first" -- validated above
            t_local = np.full_like(z, t_exp)

        if dilution:
            # Homology makes beta the Lagrangian label: the element now at
            # (r, t_local) sat at r * t_exp/t_local at t_exp, and its density
            # has fallen by (t_exp/t_local)^3 since. So the caller's n_l_of_r
            # is read as the INITIAL profile and diluted, rather than being
            # re-read as a fresh snapshot at every epoch.
            shrink = t_exp / t_local
            n_l = n_l_of_r(r * shrink) * shrink**3
            temp = temp_of_r(r * shrink)
        else:
            n_l = n_l_of_r(r)
            temp = temp_of_r(r)

        beta_z = z / (C * t_local)
        if relativity in ("exact", "worldline"):
            beta = r / (C * t_local)
            lorentz = 1.0 / np.sqrt(np.maximum(1.0 - beta**2, 1e-300))
            doppler = lorentz * (1.0 - beta_z)  # D = nu'/nu
        else:
            doppler = 1.0 - beta_z
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
            if relativity in ("exact", "worldline"):
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
        if relativity in ("exact", "worldline"):
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
