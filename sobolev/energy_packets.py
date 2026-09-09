"""The Paper IV packet contract: indivisible energy packets in homologous flow.

Paper II/III's `run_mc` propagates *photons*: each packet is one photon
(weight `w` photons), its frequency changes at every fluorescence and the
level-energy difference is booked as a deposit in the gas that is never
re-radiated. Paper IV makes the packets *indivisible energy packets*
(Lucy 2002, 2003, 2005). The contract, stated once here and implemented as
hooks in `paper2/phase1/forest_mc.run_mc(packets="energy")`:

A packet is (E_cm, nu_cm, nu_lab, mu, r, t) with

    E_cm  = w h nu_cm            its energy in the comoving frame of the gas
                                 at its position,
    E_lab = w h nu_lab = E_cm D  its energy in the lab frame, D the local
                                 Doppler factor nu_lab / nu_cm.

The two invariants:

  1. E_cm is conserved at every atomic interaction. Absorption, internal
     macroatom jumps and re-emission happen in the gas frame; whatever line
     the packet leaves through, it leaves with the energy it arrived with.
     `w` is therefore rescaled at re-emission by nu_abs,cm / nu_rest so that
     w h nu is unchanged in the comoving frame.
  2. E_lab is conserved in free flight (no interaction changes nu_lab), so
     between interactions E_cm = E_lab / D drifts with D: the O(v/c)
     adiabatic loss of an expanding flow. Its running sum over a packet's
     interactions is the Doppler work term W that `run_mc` already reports
     (E_dep_lab - E_dep_cm); under energy packets E_dep_cm = 0 by
     construction and E_dep_lab = W.

Doppler factors as `run_mc` computes them (first-order classical, or exact
worldline with the packet's own clock c t_p):

    classical : nu_cm = nu_lab (1 - r mu / c t)               D = 1/(1 - r mu/ct)
    worldline : nu_cm = nu_lab gamma (1 - beta mu),           D = 1/(gamma (1 - beta mu))
                beta = r / (c t_p), gamma = (1 - beta^2)^-1/2

Re-emission is isotropic in the comoving frame; the lab direction and
frequency follow by aberration, nu_lab = nu_rest gamma (1 + beta mu_c)
(worldline) or nu_rest / (1 - r mu / c t) (classical) -- the same D at the
new direction. So a coherent scattering changes E_lab by D_out / D_in and
E_cm not at all.

What the switch does NOT change: packet histories under photon-number
probabilities. `run_mc` never reads `w` inside its step loop, so
`packets="energy"` with the Paper III outcomes (branch, thermal, tla,
group) reproduces the photon-mode fates and frequencies bit for bit and
differs only in `w` -- this is the R1^E rung of the Paper IV ladder, the
bookkeeping change isolated from any change of transition probabilities.
Where probabilities must follow the energy language (the thermal sampler,
the bin emissivity, the kernel rows, the downward macroatom) the callers
ask for energy weights explicitly.
"""
import numpy as np

from .constants import H


def doppler_factor(r, mu, ct, relativity=None, ctime=None):
    """D = nu_lab / nu_cm at lab position r, direction cosine mu.

    `ct` is c t_exp (classical); `ctime` the packet's own c t_p (worldline)."""
    r = np.asarray(r, float); mu = np.asarray(mu, float)
    if relativity is None:
        return 1.0 / (1.0 - r * mu / ct)
    if relativity == "worldline":
        cti = ct if ctime is None else np.asarray(ctime, float)
        beta = r / cti
        gamma = 1.0 / np.sqrt(1.0 - beta * beta)
        return 1.0 / (gamma * (1.0 - beta * mu))
    raise ValueError(f"relativity must be None or 'worldline', got {relativity!r}")


def comoving_frequency(nu_lab, r, mu, ct, relativity=None, ctime=None):
    """nu_cm from nu_lab, with exactly `run_mc`'s operation order."""
    z = np.asarray(r, float) * np.asarray(mu, float)
    if relativity is None:
        return np.asarray(nu_lab, float) * (1.0 - z / ct)
    cti = ct if ctime is None else np.asarray(ctime, float)
    beta = np.asarray(r, float) / cti
    gamma = 1.0 / np.sqrt(1.0 - beta * beta)
    return np.asarray(nu_lab, float) * gamma * (1.0 - z / cti)


def lab_frequency(nu_rest, r, mu_lab, ct, relativity=None, ctime=None):
    """nu_lab of a packet re-emitted at comoving frequency nu_rest into lab
    direction mu_lab (the inverse of `comoving_frequency`)."""
    return np.asarray(nu_rest, float) * doppler_factor(r, mu_lab, ct, relativity, ctime)


def energies(w, nu_lab, r, mu, ct, relativity=None, ctime=None):
    """(E_lab, E_cm) of packets with weight w and lab frequency nu_lab."""
    e_lab = np.asarray(w, float) * H * np.asarray(nu_lab, float)
    return e_lab, e_lab / doppler_factor(r, mu, ct, relativity, ctime)


def reweight(w, nu_abs_cm, nu_rest):
    """Weight after re-emission that keeps E_cm = w h nu_cm invariant.

    Returns w * nu_abs_cm / nu_rest; coherent scattering (nu_rest ==
    nu_abs_cm) returns w unchanged exactly (x / x == 1.0 in IEEE 754)."""
    return np.asarray(w, float) * (np.asarray(nu_abs_cm, float) / np.asarray(nu_rest, float))


def scattered_lab_energy(e_lab_in, r, mu_in, mu_out, ct, relativity=None, ctime=None):
    """Lab energy after one coherent scattering at (r, mu_in -> mu_out):
    E_cm conserved, so E_lab,out = E_lab,in * D_out / D_in."""
    d_in = doppler_factor(r, mu_in, ct, relativity, ctime)
    d_out = doppler_factor(r, mu_out, ct, relativity, ctime)
    return np.asarray(e_lab_in, float) * d_out / d_in


def energy_accounting(res):
    """Paper IV's view of `run_mc`'s accounting: the identity residual, the
    comoving deposit fraction (== 0 for a radiative downward macroatom), the
    adiabatic work fraction, and the two normalisation ratios.

    Every quantity is a ratio to E_inj so it can be quoted per leg."""
    a = res["accounting"]
    e_inj = a["E_inj"]
    e_exact = e_inj - a["E_core"] - a["E_abs"]       # what leaves the core for good
    return dict(
        packets=res.get("packets", "photon"),
        identity_residual=a["identity_residual"],
        dep_cm_frac=a["E_dep_cm"] / e_inj,
        dep_lab_frac=a["E_dep_lab"] / e_inj,
        work_frac=a["W"] / e_inj,
        core_frac=a["E_core"] / e_inj,
        abs_frac=a["E_abs"] / e_inj,
        thermal_frac=a.get("E_thermal", 0.0) / e_inj,
        esc_frac=a["E_esc"] / e_inj,
        # conserving / equilibrium scale - 1: the grey renormalisation the
        # Paper III harness needed; ~0 when transport itself conserves energy
        renorm_ratio=(e_exact / a["E_esc"] - 1.0) if a["E_esc"] > 0 else np.nan,
        n_core_passes=int(res.get("n_core_passes_total", 0)),
        n_dead_end=int(res.get("n_dead_end", 0)),
        n_kpackets=int(res.get("n_kpackets", 0)),
    )
