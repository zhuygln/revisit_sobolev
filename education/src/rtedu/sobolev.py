"""Homologous expansion and the Sobolev resonance (chapter 4, the E1 capstone).

A packet of fixed lab-frame frequency nu_lab moving radially outward through
homologous flow v = r/t sees the comoving frequency

    nu_com(r) = nu_lab (1 - v/c) = nu_lab (1 - r / (c t))      (first order),

which falls as it travels: it redshifts into every line below its launch
frequency in turn, at the radius where nu_com(r_res) = nu_line. The line's
Sobolev depth in homologous flow is

    tau_S = (pi e^2 / m_e c) f n_l lambda_0 t,

the same formula shape as sobolev/optical_depth.py::tau_sobolev in the
production code (written out here so the chapter needs no import). Two
different quantities carry the name "Sobolev" (PI amendment 8):

    P_int(tau) = 1 - e^{-tau}          the interaction probability,
                                        the coin the transport flips;
    beta(tau)  = (1 - e^{-tau}) / tau   the escape probability of a photon
                                        emitted in the line; beta -> 1 as
                                        tau -> 0 and beta -> 1/tau as
                                        tau -> infinity. It enters in E2.
"""
import numpy as np

from . import C, SIGMA_CLASSICAL


def comoving_frequency(nu_lab, r, t, mu=1.0, c=C):
    """First-order Doppler shift for a packet at radius r in homologous flow
    v = r/t, direction cosine mu with respect to the radial direction."""
    return nu_lab * (1.0 - mu * np.asarray(r, float) / (c * t))


def resonance_position(nu_lab, nu_line, t, mu=1.0, c=C):
    """The radius at which nu_com(r) = nu_line, or nan if the packet never
    reaches that line (nu_line above nu_lab, or a receding geometry)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        r = c * t * (1.0 - np.asarray(nu_line, float) / nu_lab) / mu
    return np.where(r >= 0.0, r, np.nan)


def tau_sobolev_toy(n_lower, f_osc, lambda0_cm, t):
    """The Sobolev depth of one line in homologous flow (no stimulated
    emission, no level degeneracy: the toy)."""
    return SIGMA_CLASSICAL * f_osc * n_lower * lambda0_cm * t


def interaction_probability(tau):
    """P_int = 1 - e^{-tau}: the coin the transport flips at a resonance."""
    return -np.expm1(-np.asarray(tau, float))


def escape_probability(tau):
    """beta = (1 - e^{-tau}) / tau, with beta(0) = 1 (a different quantity
    from P_int: the fraction of photons emitted in the line that escape it)."""
    tau = np.asarray(tau, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = -np.expm1(-tau) / tau
    return np.where(tau > 1e-12, beta, 1.0)


def resonance_crossings(nu_lab, nu_lines, tau_lines, r0, r_out, t, rng=None, mu=1.0):
    """Walk one packet outward from r0 to r_out and list every line it
    resonates with, in the order met: (line index, r_res, tau, p_int,
    interacted). With rng=None no coin is flipped (interacted is False
    everywhere: the deterministic list of crossings); with an rng the walk
    stops at the first interaction.

    A line is met when nu_com(r) sweeps down through nu_line, i.e. when
    nu_line < nu_com(r0) and r_res < r_out. Lines are met in decreasing
    frequency order because nu_com decreases monotonically outward.
    """
    nu_lines = np.asarray(nu_lines, float); tau_lines = np.asarray(tau_lines, float)
    nu_start = comoving_frequency(nu_lab, r0, t, mu)
    r_res = resonance_position(nu_lab, nu_lines, t, mu)
    ok = (nu_lines < nu_start) & (r_res < r_out)
    order = np.argsort(-nu_lines[ok])
    out = []
    for k in np.flatnonzero(ok)[order]:
        p = float(interaction_probability(tau_lines[k]))
        hit = bool(rng.random() < p) if rng is not None else False
        out.append(dict(line=int(k), r_res=float(r_res[k]), tau=float(tau_lines[k]), p_int=p, interacted=hit))
        if hit:
            break
    return out


def first_interaction(nu_lab, nu_lines, tau_lines, r0, r_out, t, rng, mu=1.0):
    """The crossing at which the packet interacts, or None if it escapes."""
    cr = resonance_crossings(nu_lab, nu_lines, tau_lines, r0, r_out, t, rng, mu)
    return cr[-1] if cr and cr[-1]["interacted"] else None


def toy_line_list(t=2.0 * 86400.0, nu_lab=6.0e14, v_max_c=0.15, f_osc=0.05):
    """The twelve-line toy forest of chapter 4 (one source for the notebook,
    the video and the tests): frequencies as fractions of the launch
    frequency, toy lower-level densities chosen so that the Sobolev depths
    span 0.01 to 2, and an outer edge at v_max_c that leaves the four
    reddest lines out of reach. Returns a dict."""
    frac = np.array([0.985, 0.97, 0.955, 0.94, 0.925, 0.905, 0.885, 0.86, 0.84, 0.815, 0.79, 0.77])
    n_l = np.array([4.0, 12.0, 25.0, 0.6, 60.0, 20.0, 2.5, 120.0, 30.0, 6.0, 80.0, 1.5])   # cm^-3
    nu_lines = frac * nu_lab
    tau = tau_sobolev_toy(n_l, f_osc, C / nu_lines, t)
    return dict(t=t, nu_lab=nu_lab, v_max_c=v_max_c, r_out=v_max_c * C * t, frac=frac, n_l=n_l, f_osc=f_osc,
                nu_lines=nu_lines, tau=tau)
