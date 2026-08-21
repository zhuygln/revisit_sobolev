"""Analytic Sobolev attenuation of core light through an expanding shell.

This is the *Sobolev approximation proper*: each line is a delta-function
resonance at one plane, attenuating by exp(-tau_S), with no profile width and
no interaction between lines. In the pure-absorption LTE setup used
throughout this project (no scattering, no thermal emission from the shell)
that prediction is exact analytics -- a Monte Carlo implementation of
per-line Sobolev interactions would add nothing but noise.

The same machinery yields the expansion-opacity prediction by damping each
crossing: pass damp = lambda t: 1 - np.exp(-t) (see optical_depth and
docs/paper Appendix A.4 for why that substitution is what the binned
formalism actually applies).

CONVENTION (Finding F8): there is no emission term here, so results must be
compared against the deterministic solver run with T_shell -> 0, or against
SEDONA's resolved mode, which likewise deposits absorbed energy without
re-emitting it at fixed temperature. Comparing against a solver run that
includes S = B_nu(T_shell) reintroduces a 3-6% offset that is a convention
mismatch, not physics.
"""

import numpy as np
from scipy.special import erf

from .constants import C
from .optical_depth import tau_sobolev
from .rays import RaySet, core_rays_midpoint


def sobolev_attenuation(
    nu_grid,
    lines,
    r_core,
    r_out,
    t_exp,
    n_ref,
    damp=None,
    n_p=200,
    relativity=None,
    n_of_beta=None,
    rays=None,
):
    """Fraction of core continuum transmitted at each observer frequency.

    A ray at impact parameter p leaving the core crosses line k's resonance
    plane iff sqrt(r_core^2 - p^2) < z_res,k < sqrt(r_out^2 - p^2), and is
    attenuated by exp(-sum_k damp(tau_k)) over the lines it crosses.
    Averaging exp(...) over p with weight p is essential: planes closer than
    r_core are crossed by only the outer rays, so plane-counting without the
    average predicts visibly too-shallow troughs.

    Parameters
    ----------
    nu_grid : observer-frame frequencies (Hz)
    lines : sequence of (nu0, f_osc) or (nu0, f_osc, pop_frac) tuples -- the
        same form `formal_transfer.emergent_luminosity` accepts, so both legs
        consume one line list.
    r_core, r_out : core and outer shell radii (cm)
    t_exp : ejecta age (s); homologous v = r / t_exp
    n_ref : reference density multiplying each line's pop_frac, i.e. the same
        quantity the solver's `n_l_of_r` returns (ion density, or total mass
        density when pop_frac carries ions-per-gram).
    damp : None for Sobolev proper (tau used as is), or a callable applied to
        each tau_k -- e.g. lambda t: 1 - np.exp(-t) for expansion opacity.
    n_p : number of impact-parameter rays across the core disk.
    relativity : None (default) keeps the classical first-order plane
        `z_res = c t (1 - nu0/nu)` and the beta-free tau_S -- the path every
        result in Paper I was computed with, left bit-for-bit intact.
        "first" keeps the classical plane but carries the factor nu0/nu that
        the resolved first-order integral has and the classical leg drops
        (a 0.3-1% effect across a 1000-3000 km/s shell), so it is the exact
        delta-resonance limit of `emergent_luminosity(relativity="first")`
        and of `resolved_attenuation(sweep="first")`. "worldline" and
        "exact" name the same two treatments as `emergent_luminosity`; see
        below.
    rays : optional `sobolev.rays.RaySet`. Replaces the n_p midpoint rule with
        explicit impact parameters and weights, so that this leg and the
        resolved legs can be evaluated on IDENTICAL rays. Only the core rays
        of the set are used (envelope rays see no continuum).
    n_of_beta : callable beta -> n_l/n_ref at the resonance, i.e. the shape of
        the INITIAL (t_exp) profile as a function of the fluid element's
        velocity. None means uniform. Ignored when relativity is None, where
        n_ref alone sets the depth.

    Returns the transmitted fraction on nu_grid (1.0 = untouched continuum).

    RELATIVISTIC MODES. Both make tau and the resonance locus functions of the
    impact parameter as well as the frequency, so `tau_k` and `z_res` become
    (n_p, n_nu) arrays. The loop structure is otherwise unchanged and `damp` is
    still applied to tau per crossing before accumulation, so the
    expansion-opacity leg remains the same code path with one keyword flipped.

    "worldline" -- the physical law. With ct(z) = z + Z0 along the photon path
        (Z0 = c t_exp - z0, anchored where the ray leaves the core),

            D = gamma (1 - beta_z) = Z0 / sqrt(Z0^2 + 2 Z0 z - p^2)

        so D = nu0/nu is LINEAR in z and the locus has a single root:

            z_res(p, nu) = (Z0/2)[(nu/nu0)^2 - 1] + p^2/(2 Z0).

        Since |dD/dz| = D^3/Z0 and D(z_res + Z0) = gamma Z0, the impact
        parameter cancels identically and tau/tau_S(t_res) = 1/gamma for EVERY
        ray. tau_S is evaluated at the resonance's own epoch t_res = (z_res +
        Z0)/c and its own diluted density, n_of_beta(beta_res) (t_exp/t_res)^3.

    "exact" -- frozen snapshot, exact SR, kept because it is what an
        observer-frame integration at fixed t computes. Here the locus really
        is the two-root quadratic behind Jeffery's CD/CP surfaces,

            u+- = [1 +- x sqrt(x^2 (1-a^2) - a^2)] / (1 + x^2),

        with u = beta_z, a = p/(c t_exp), x = nu0/nu. Both roots are carried;
        u+ sits above beta_z = 1 - a^2 and so never falls inside a shell with
        beta_out <~ 0.35, which is *why* the plane picture was adequate at the
        velocities of Paper I. tau comes from `tau_sobolev_frozen` called with
        beta_z and beta distinct -- the general-ray law
        (1-beta_z)^2/[gamma(1 - beta_z - beta_p^2)] -- rather than letting beta
        default to |beta_z|, which is the radial-ray special case.

    Note the asymmetry: the CD/CP geometry is a property of the FROZEN problem.
    In the mode that is physically correct the locus is linear and the geometry
    objection does not arise.
    """
    nu_grid = np.asarray(nu_grid, dtype=float)
    lines = [
        (entry[0], entry[1], entry[2] if len(entry) > 2 else 1.0)
        for entry in lines
    ]

    # Midpoint rays across the core disk (or the caller's shared RaySet).
    if rays is None:
        p = core_rays_midpoint(r_core, n_p)
        w_ray = p
    else:
        p = np.asarray(rays.p, dtype=float)[rays.is_core]
        w_ray = np.asarray(rays.w, dtype=float)[rays.is_core]
    z_lo = np.sqrt(np.maximum(r_core**2 - p**2, 0.0))
    z_hi = np.sqrt(np.maximum(r_out**2 - p**2, 0.0))

    if relativity not in (None, "first", "worldline", "exact"):
        raise ValueError(
            "relativity must be None, 'first', 'worldline' or 'exact', "
            f"got {relativity!r}"
        )

    att = np.zeros((p.size, nu_grid.size))
    for nu0, f_osc, pop in lines:
        if relativity == "first":
            # Classical plane, but the resolved first-order integral gives
            # tau_S * (nu0/nu) at resonance: the Doppler sweep rate is
            # dnu'/dz = nu/(ct), not nu0/(ct). Kept separate so the None
            # path stays untouched.
            tau_k = tau_sobolev(f_osc, pop * n_ref, C / nu0, t_exp) * (nu0 / nu_grid)
            if damp is not None:
                tau_k = damp(tau_k)
            z_res = C * t_exp * (1.0 - nu0 / nu_grid)
            crossed = (z_res[None, :] > z_lo[:, None]) & (z_res[None, :] < z_hi[:, None])
            att += np.where(crossed, tau_k[None, :], 0.0)
        elif relativity is None:
            # Untouched: the classical path every Paper I number came from.
            tau_k = tau_sobolev(f_osc, pop * n_ref, C / nu0, t_exp)
            if damp is not None:
                tau_k = damp(tau_k)
            if tau_k <= 0.0:
                continue
            # Resonance plane for this line at each observer frequency.
            z_res = C * t_exp * (1.0 - nu0 / nu_grid)
            crossed = (z_res[None, :] > z_lo[:, None]) & (z_res[None, :] < z_hi[:, None])
            att += np.where(crossed, tau_k, 0.0)
        else:
            if f_osc * pop <= 0.0:
                continue
            crossed, tau_k = _relativistic_crossings(
                p, z_lo, z_hi, nu_grid, nu0, f_osc, pop * n_ref, t_exp,
                r_core, r_out, relativity, n_of_beta,
            )
            # Zero the non-crossing entries before damping so that a locus
            # outside the shell cannot leak a NaN into the sum; damp(0) = 0,
            # so this does not change what a real crossing contributes.
            tau_k = np.where(crossed, tau_k, 0.0)
            if damp is not None:
                tau_k = damp(tau_k)
            att += np.where(crossed, tau_k, 0.0)

    return np.sum(w_ray[:, None] * np.exp(-att), axis=0) / np.sum(w_ray)


def resolved_attenuation(
    nu_grid,
    lines,
    r_core,
    r_out,
    t_exp,
    n_ref,
    v_doppler,
    n_p=200,
    rays=None,
    cutoff_widths=None,
    sweep="classical",
):
    """Closed-form RESOLVED attenuation: Gaussian profiles, uniform n_l,
    pure absorption, first-order Doppler.

    For a line of rest frequency nu0 on a ray crossing the shell from z_lo to
    z_hi, the resolved optical depth at observer frequency nu is the profile
    integrated over the part of the resonance that lies inside the shell:

        tau(p, nu) = tau_S [nu0/nu] * 1/2 [erf((z_hi - z_k)/Delta)
                                           - erf((z_lo - z_k)/Delta)],
        z_k = c t (1 - nu0/nu),     Delta = v_D t (nu0/nu),

    the formula `experiments/boundary/run.py` verified against a brute-force z
    integration to four decimals, applied per line and per ray with no z loop
    at all. Cost O(n_p n_nu n_lines) and INDEPENDENT of v_D, so the 1 km/s
    frontier that costs the brute-force solver hours is free here.

    This is the deterministic reference the Sobolev leg is differenced
    against (referee Comment 5): same rays, same source convention (no
    emission), same populations, same transport treatment. It reproduces
    `sobolev_attenuation` as v_D -> 0 (the erf bracket -> the crossing test)
    and the F12 boundary law (a plane at an edge is half-counted).

    Assumptions to keep visible: n_l uniform in r (true in every experiment
    of this project -- the tau_max sweeps scale a constant density); pure
    Gaussian profile (Voigt wings beyond the cutoff are ~1e-5 of tau at these
    damping parameters); first-order Doppler only. `sweep="first"` keeps the
    [nu0/nu] factor that the resolved integral has; `sweep="classical"` drops
    it, matching the classical Sobolev plane exactly (for the v_D -> 0 test).
    `cutoff_widths` clips the erf arguments at +-cutoff, mimicking SEDONA's
    hard-coded +-5 width truncation; numerically irrelevant (erf(5) = 1 -
    1.5e-12) but set it for like-for-like and say so.
    """
    if sweep not in ("classical", "first"):
        raise ValueError(f"sweep must be 'classical' or 'first', got {sweep!r}")
    nu_grid = np.asarray(nu_grid, dtype=float)
    lines = [
        (entry[0], entry[1], entry[2] if len(entry) > 2 else 1.0)
        for entry in lines
    ]
    if rays is None:
        p = core_rays_midpoint(r_core, n_p)
        w_ray = p
    else:
        p = np.asarray(rays.p, dtype=float)[rays.is_core]
        w_ray = np.asarray(rays.w, dtype=float)[rays.is_core]
    z_lo = np.sqrt(np.maximum(r_core**2 - p**2, 0.0))[:, None]
    z_hi = np.sqrt(np.maximum(r_out**2 - p**2, 0.0))[:, None]

    att = np.zeros((p.size, nu_grid.size))
    for nu0, f_osc, pop in lines:
        tau_k = tau_sobolev(f_osc, pop * n_ref, C / nu0, t_exp)
        if tau_k <= 0.0:
            continue
        x = nu0 / nu_grid
        z_k = C * t_exp * (1.0 - x)
        delta = v_doppler * t_exp * x
        a_hi = (z_hi - z_k[None, :]) / delta[None, :]
        a_lo = (z_lo - z_k[None, :]) / delta[None, :]
        if cutoff_widths is not None:
            a_hi = np.clip(a_hi, -cutoff_widths, cutoff_widths)
            a_lo = np.clip(a_lo, -cutoff_widths, cutoff_widths)
        frac = 0.5 * (erf(a_hi) - erf(a_lo))
        scale = x if sweep == "first" else 1.0
        att += tau_k * scale * frac
    return np.sum(w_ray[:, None] * np.exp(-att), axis=0) / np.sum(w_ray)


def crossing_depths(nu_grid, lines, r_core, r_out, t_exp, n_ref, n_p=200, rays=None):
    """Per-ray, per-frequency Sobolev sums S and E over the crossed resonances.

        S(p, nu) = sum_k tau_k                 -- what Sobolev transfer exponentiates
        E(p, nu) = sum_k (1 - exp(-tau_k))     -- what expansion opacity integrates to

    E is the expected number of line interactions per crossing for a photon
    that is counted but not removed (Karp's mean-free-path statistic), and
    expansion opacity preserves it exactly. Transmission is exp(-S) (product
    of Bernoulli survivals) against exp(-E) (Poisson with the same mean), so
    F_exp / F_Sob = E_w[exp(S - E)] with weights w ~ e^{-S} -- an identity,
    tested as one. Classical plane, same rays as `sobolev_attenuation(None)`.

    Returns (S, E, p, w) with S, E of shape (n_p, n_nu).
    """
    nu_grid = np.asarray(nu_grid, dtype=float)
    lines = [
        (entry[0], entry[1], entry[2] if len(entry) > 2 else 1.0)
        for entry in lines
    ]
    if rays is None:
        p = core_rays_midpoint(r_core, n_p)
        w_ray = p
    else:
        p = np.asarray(rays.p, dtype=float)[rays.is_core]
        w_ray = np.asarray(rays.w, dtype=float)[rays.is_core]
    z_lo = np.sqrt(np.maximum(r_core**2 - p**2, 0.0))
    z_hi = np.sqrt(np.maximum(r_out**2 - p**2, 0.0))
    S = np.zeros((p.size, nu_grid.size))
    E = np.zeros_like(S)
    for nu0, f_osc, pop in lines:
        tau_k = tau_sobolev(f_osc, pop * n_ref, C / nu0, t_exp)
        if tau_k <= 0.0:
            continue
        z_res = C * t_exp * (1.0 - nu0 / nu_grid)
        crossed = (z_res[None, :] > z_lo[:, None]) & (z_res[None, :] < z_hi[:, None])
        S += np.where(crossed, tau_k, 0.0)
        E += np.where(crossed, 1.0 - np.exp(-tau_k), 0.0)
    return S, E, p, w_ray


def _relativistic_crossings(
    p, z_lo, z_hi, nu_grid, nu0, f_osc, n_line, t_exp,
    r_core, r_out, relativity, n_of_beta,
):
    """Resonance locus and optical depth per (impact parameter, frequency).

    Returns (crossed, tau), both (n_p, n_nu). See `sobolev_attenuation`'s
    docstring for the two laws; the algebra is derived there and pinned
    against a numerical root-find in tests/test_relativistic.py.
    """
    shape = lambda b: 1.0 if n_of_beta is None else n_of_beta(b)
    lambda0 = C / nu0
    ct = C * t_exp
    b_core, b_out = r_core / ct, r_out / ct

    if relativity == "worldline":
        # Anchor the clock where the ray leaves the core, as the solver does.
        big_z = (ct - z_lo)[:, None]
        y2 = (nu_grid / nu0)[None, :] ** 2
        z_res = 0.5 * big_z * (y2 - 1.0) + (p[:, None] ** 2) / (2.0 * big_z)
        ct_res = z_res + big_z
        beta_res = np.sqrt(p[:, None] ** 2 + z_res**2) / ct_res
        lorentz = 1.0 / np.sqrt(np.maximum(1.0 - beta_res**2, 1e-300))

        # tau_S at the resonance's OWN epoch and its own diluted density.
        t_res = ct_res / C
        n_res = n_line * shape(beta_res) * (t_exp / t_res) ** 3
        tau = tau_sobolev(f_osc, n_res, lambda0, t_res) / lorentz

        # Shell membership on beta, not z: homology fixes the velocity span at
        # every epoch, so this needs no moving-boundary root.
        crossed = (z_res > z_lo[:, None]) & (beta_res > b_core) & (beta_res < b_out)
        return crossed, tau

    # "exact": frozen snapshot, two-root quadratic, z bounds frozen too.
    a = (p / ct)[:, None]
    x = (nu0 / nu_grid)[None, :]
    disc = x**2 * (1.0 - a**2) - a**2
    root = np.sqrt(np.maximum(disc, 0.0))
    tau = np.zeros((p.size, nu_grid.size))
    crossed = np.zeros((p.size, nu_grid.size), dtype=bool)
    for sign in (-1.0, +1.0):
        u = (1.0 + sign * x * root) / (1.0 + x**2)
        z_res = u * ct
        hit = (
            (disc > 0.0)
            & (z_res > z_lo[:, None])
            & (z_res < z_hi[:, None])
        )
        if not hit.any():
            continue
        beta = np.sqrt(a**2 + u**2)
        tau_root = tau_sobolev_frozen(
            f_osc, n_line * shape(beta), lambda0, t_exp, beta_z=u, beta=beta
        )
        # A ray can pierce a CD/CP surface twice; optical depths add.
        tau = tau + np.where(hit, tau_root, 0.0)
        crossed = crossed | hit
    return crossed, tau


def tau_sobolev_relativistic(f_osc, n_l, lambda0, t_exp, beta):
    """Sobolev optical depth for homologous flow along the photon worldline.

    THIS is the physical law. The photon takes time to cross the resonance,
    so the ejecta age advances as it goes. For homologous flow the photon
    worldline gives

        dbeta_vec/dt = (n_hat - beta_vec)/t   =>   db/dt = (1 - b)/t,

    and with D = gamma(1-b), dD/dt = -D gamma^2 (1-b)/t. Using the opacity
    transformation chi_lab = D chi' and the resonance condition D nu = nu0,
    the profile integral collapses to

        tau = D chi' c / |dnu'/dt| = sigma f n_l c t D / (nu0 gamma^2 (1-b))

    and since D = gamma(1-b) the bracket reduces to a single 1/gamma:

        tau_rel / tau_S = 1/gamma = sqrt(1 - beta^2).

    There is NO first-order term: 1/gamma = 1 - beta^2/2 + O(beta^4). At
    beta = 0.1/0.2/0.3 the correction is 0.5%/2.0%/4.6%.

    tau_S must be evaluated with the density and the ejecta age local to the
    resonance; with that convention the result is anchoring-independent.

    Contrast `tau_sobolev_frozen`, which holds t fixed and yields
    (1-beta)/gamma -- a first-order effect that is an artifact of the frozen
    snapshot, not physics. Both laws are confirmed against a first-principles
    integration to five decimals in tests/test_relativistic.py.
    """
    beta = np.asarray(beta, dtype=float)
    lorentz = 1.0 / np.sqrt(np.maximum(1.0 - beta**2, 1e-300))
    return tau_sobolev(f_osc, n_l, lambda0, t_exp) / lorentz


def tau_sobolev_frozen(f_osc, n_l, lambda0, t_exp, beta_z, beta=None):
    """Sobolev optical depth for a FROZEN SNAPSHOT of homologous flow.

    Retained because it is what an observer-frame integration at fixed t
    computes -- including this project's solver before the worldline mode
    existed -- and because the difference from the physical law is itself a
    result. It is NOT the correct relativistic Sobolev depth; see
    `tau_sobolev_relativistic`.

    Original derivation follows.

    Derivation (the relativistic counterpart of Appendix A.2). Along a ray
    toward the observer the gas at coordinate z has line-of-sight velocity
    beta_z = z/(ct) and speed beta = r/(ct), so the comoving frequency is

        nu' = D nu,      D = gamma (1 - beta_z),   gamma = 1/sqrt(1-beta^2).

    Two separate Doppler factors then enter the optical depth:

    1. The opacity itself. nu*chi_nu is a Lorentz invariant, so the lab-frame
       opacity is chi_lab = D chi'_comoving -- absorption is WEAKER in the lab
       frame by D.
    2. The sweep rate. tau = integral of chi dz collapses to chi / |dnu'/dz|,
       and the resonance condition nu' = nu0 fixes nu = nu0/D, so the sweep
       rate carries a further factor.

    Keeping gamma exact and using db/dz = 1/(ct),

        dnu'/dz = (nu/(ct)) [ gamma^3 b (1-b) - gamma ],      b = beta_z,

    and the resonance condition nu' = nu0 fixes nu = nu0/D. Since the profile
    integrates to unity in nu', tau = D chi' / |dnu'/dz| carries only ONE
    explicit D, but substituting nu = nu0/D puts a second one in the sweep
    rate:

        tau_rel / tau_S  =  D^2 / |gamma^3 b (1-b) - gamma| .

    Expanding, D^2 ~ 1 - 2b while the denominator ~ 1 - b, so

        tau_rel  ~  tau_S (1 - beta_z)      to first order,

    NOT (1 - beta_z)^2: the denominator cancels one of the two Doppler
    factors. This matters for interpreting Finding F3. The solver's "first"
    mode, which omits the opacity transformation entirely, also yields
    tau_S (1 - beta_z) at first order -- so the empirically measured
    exp[-tau_S(1-beta)] of F3 is reproduced by this frozen law. For radial
    rays the expression above equals (1-beta)/gamma exactly.

    What that means was misread twice. It is not a frame ambiguity, and it is
    not the leading relativistic correction either: holding t fixed while a
    photon crosses the resonance is simply the wrong problem, and the O(beta)
    term is the price of that choice. See `tau_sobolev_relativistic`.

    Parameters
    ----------
    beta_z : line-of-sight velocity / c at the resonance point.
    beta : total speed / c there. Defaults to |beta_z| (exact for a ray
        through the centre, where the velocity is purely radial along z).

    Returns the frozen-snapshot Sobolev optical depth. Reduces to
    `tau_sobolev` as beta_z -> 0.
    """
    beta_z = np.asarray(beta_z, dtype=float)
    if beta is None:
        beta = np.abs(beta_z)
    beta = np.asarray(beta, dtype=float)

    lorentz = 1.0 / np.sqrt(np.maximum(1.0 - beta**2, 1e-300))
    d_factor = lorentz * (1.0 - beta_z)  # D = nu'/nu

    # |dnu'/dz| in units of nu0/(ct), with nu = nu0/D at resonance:
    #   dnu'/dz = nu [ gamma^3 beta_z (1-beta_z)/(ct) - gamma/(ct) ]
    # so   |dnu'/dz| (ct/nu0) = |gamma^3 beta_z (1-beta_z) - gamma| / D.
    sweep = np.abs(lorentz**3 * beta_z * (1.0 - beta_z) - lorentz) / d_factor

    tau_nonrel = tau_sobolev(f_osc, n_l, lambda0, t_exp)
    return tau_nonrel * d_factor / sweep


def expansion_damp(tau):
    """The per-crossing substitution the binned expansion-opacity formalism
    applies: tau -> 1 - exp(-tau), capped at one effective optical depth."""
    return 1.0 - np.exp(-tau)
