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

from .constants import C
from .optical_depth import tau_sobolev


def sobolev_attenuation(
    nu_grid,
    lines,
    r_core,
    r_out,
    t_exp,
    n_ref,
    damp=None,
    n_p=200,
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

    Returns the transmitted fraction on nu_grid (1.0 = untouched continuum).
    """
    nu_grid = np.asarray(nu_grid, dtype=float)
    lines = [
        (entry[0], entry[1], entry[2] if len(entry) > 2 else 1.0)
        for entry in lines
    ]

    # Midpoint rays across the core disk.
    p = np.linspace(0.0, r_core, n_p, endpoint=False) + r_core / (2 * n_p)
    z_lo = np.sqrt(np.maximum(r_core**2 - p**2, 0.0))
    z_hi = np.sqrt(np.maximum(r_out**2 - p**2, 0.0))

    att = np.zeros((p.size, nu_grid.size))
    for nu0, f_osc, pop in lines:
        tau_k = tau_sobolev(f_osc, pop * n_ref, C / nu0, t_exp)
        if damp is not None:
            tau_k = damp(tau_k)
        if tau_k <= 0.0:
            continue
        # Resonance plane for this line at each observer frequency.
        z_res = C * t_exp * (1.0 - nu0 / nu_grid)
        crossed = (z_res[None, :] > z_lo[:, None]) & (z_res[None, :] < z_hi[:, None])
        att += np.where(crossed, tau_k, 0.0)

    return np.sum(p[:, None] * np.exp(-att), axis=0) / np.sum(p)


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
