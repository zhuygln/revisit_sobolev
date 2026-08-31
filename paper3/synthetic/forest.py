"""Synthetic line forests with independently dialled structure.

F30/F31 read the grouped-opacity failure as a two-regime story separated by
line spacing, but that reading rests on two lanthanides. To make it a law we
need forests whose crowding, saturation, spacing and *redistribution range* can
be varied one at a time -- which no real atom lets you do.

`ForestAtom.__init__` already accepts arbitrary arrays and builds
`branch_lines` / `branch_cum` / `exit_cum` / `beta_all` from `upper` and `A`
alone, so no new class is needed. Two constructors:

  synthetic_forest  -- the exit rule is DIALLED. Each absorbing line's upper
      level gets `n_exit` downward channels placed at a controlled |d ln lambda|,
      so redistribution range is an axis orthogonal to crowding and tau.
  synthetic_ladder  -- the exit rule EMERGES from a parameterized level ladder
      and its A values, as in a real atom. Slower to control, and its job is to
      validate that the dial reproduces what real cascade structure does.

Both follow `tests/test_forest_mc.py:three_level`: exit channels are given
f_osc = 0 and n_lower = 0, so tau = 0 excludes them from the opacity set while
they remain in the branching tables, and `beta_all` is forced to 1 for them
(they escape freely). Line placement follows
`experiments/highbeta/pilot.py:legs` -- self-similar in ln lambda, with n_lower
inverted from the target tau.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper2/phase1"))
from forest_mc import ForestAtom

from sobolev.constants import C
from sobolev.optical_depth import tau_sobolev

T_EXP = 86400.0
LAM_MID = 5000e-8       # cm; the forest is centred here
F_OSC = 0.1


def _tau_draw(n, tau, spread, rng):
    """Optical depths: a scalar gives a uniform forest, spread > 0 gives a
    lognormal spread about it (the realistic case -- real forests are not
    monodisperse, and saturation statistics depend on the tail)."""
    tau = np.asarray(tau, float)
    if tau.ndim:                      # explicit per-line array
        return tau
    if spread <= 0:
        return np.full(n, float(tau))
    return float(tau) * np.exp(rng.normal(0.0, spread, n))


def synthetic_forest(n_lines=200, tau=5.0, tau_spread=0.0, span=0.2,
                     n_exit=2, dlnlam=0.1, dlnlam_spread=0.0, f_return=0.5,
                     jitter=0.0, seed=0, t_exp=T_EXP, lam_mid=LAM_MID,
                     delocalize=0.0, exit_tau=0.0, n_sink=3):
    """A forest whose redistribution range is a dial.

    Parameters
    ----------
    n_lines : absorbing lines (each with its own upper level).
    tau, tau_spread : Sobolev depth per line; lognormal sigma if spread > 0.
    span : total extent of the forest in ln lambda. With n_lines fixed this is
        the spacing axis; crowding is n_lines / span.
    n_exit : downward channels per upper level. **n_exit = 1 is the star
        topology** -- the only way down is the absorbing line itself, so the
        photon re-emerges where it was absorbed and redistribution range is
        identically zero (coherent scattering).
    dlnlam, dlnlam_spread : characteristic |d ln lambda| of the exit channels,
        and its lognormal spread. THE REDISTRIBUTION-RANGE AXIS.
    f_return : share of the upper level's A that goes back down the absorbing
        line. Sets how often a photon returns to its own resonance versus
        leaving it, independently of how far the exits reach.
    jitter : if > 0, line positions are perturbed by this fraction of the mean
        spacing (a regular comb is an unphysical special case -- it makes every
        group identical).
    delocalize : probability that an exit channel is placed ANYWHERE in the
        forest rather than at +-dlnlam from its own absorbing line.

        F37: with purely local exits the model redistributes energy near where
        it was absorbed and never delivers a NET INFLOW to a sub-band, so the
        measured band only ever darkens and the closure can only ever be too
        opaque. Real forests do the opposite -- Tm II's 3800-3955 A band
        transmits 1.049, MORE than the continuum entering it, because energy
        absorbed elsewhere is re-emitted into it. An upper level in a real atom
        decays to lower levels spread across the whole term structure, so its
        exit wavelengths are spread across the whole forest; that is what
        delocalize models. Without it the too-bright branch of the sign
        boundary (§4.32) does not exist in the model.
    exit_tau, n_sink : optical depth carried by the exit lines, and how many
        shared lower levels they terminate on.

        With exit_tau = 0 (the original) each channel ends on its own
        unpopulated sink, so a photon leaving through one can never be
        re-absorbed: the exit is terminal and the cascade stops after one step.
        Real forests refill a band MORE as saturation rises, because more
        absorption elsewhere feeds more re-emission and the cascade continues;
        a terminal exit refills LESS, which is why the delocalized model
        produced a sign change running the wrong way (§4.34). Giving the exit
        lines their own opacity, on a small set of SHARED populated lower
        levels as a real atom has, makes absorption on an exit line return the
        photon to that upper level and branch again -- the recurrent
        fluorescence cycle.

    Returns (atom, info) with info carrying the requested parameters.
    """
    rng = np.random.default_rng(seed)
    if n_exit < 1:
        raise ValueError("n_exit must be >= 1 (1 = star topology)")

    # --- absorbing lines: self-similar placement in ln lambda
    u = np.linspace(-0.5, 0.5, n_lines) if n_lines > 1 else np.zeros(1)
    if jitter > 0 and n_lines > 1:
        u = u + rng.normal(0.0, jitter / n_lines, n_lines)
    lam_abs = lam_mid * np.exp(u * span)
    order = np.argsort(lam_abs)
    lam_abs = lam_abs[order]
    nu_abs = C / lam_abs
    tau_abs = _tau_draw(n_lines, tau, tau_spread, rng)[order] if np.ndim(tau) \
        else _tau_draw(n_lines, tau, tau_spread, rng)

    # n_lower inverted from the target tau, per line
    n_low = np.array([t / tau_sobolev(F_OSC, 1.0, lam, t_exp)
                      for t, lam in zip(tau_abs, lam_abs)])

    # --- assemble the line table. Upper level u_i belongs to absorbing line i.
    nu0 = list(nu_abs)
    f_osc = [F_OSC] * n_lines
    n_lower = list(n_low)
    A = [f_return if n_exit > 1 else 1.0] * n_lines
    lower = [0] * n_lines                      # ground
    upper = list(range(1, n_lines + 1))

    n_ch = n_exit - 1
    if n_ch:
        a_each = (1.0 - f_return) / n_ch
        sink = n_lines + 1                     # first sink level index
        n_sink_eff = max(1, n_sink) if exit_tau > 0 else n_ch
        lo_ln, hi_ln = np.log(lam_abs.min()), np.log(lam_abs.max())
        for i in range(n_lines):
            d = dlnlam * np.exp(rng.normal(0.0, dlnlam_spread, n_ch)) \
                if dlnlam_spread > 0 else np.full(n_ch, dlnlam)
            sgn = rng.choice([-1.0, 1.0], n_ch)   # blueward and redward exits
            lam_x = lam_abs[i] * np.exp(sgn * d)
            if delocalize > 0:
                # a fraction of channels land anywhere in the forest, so a
                # photon absorbed outside a sub-band can be re-emitted into it
                far = rng.uniform(size=n_ch) < delocalize
                if far.any():
                    lam_x = lam_x.copy()
                    lam_x[far] = np.exp(rng.uniform(lo_ln, hi_ln, int(far.sum())))
            for k in range(n_ch):
                nu0.append(C / lam_x[k])
                if exit_tau > 0:
                    # a populated shared lower level: the exit line absorbs, so
                    # the photon returns to this upper level and branches again
                    f_osc.append(F_OSC)
                    n_lower.append(exit_tau / tau_sobolev(F_OSC, 1.0, lam_x[k], t_exp))
                    lower.append(sink + int(rng.integers(n_sink_eff)))
                else:
                    f_osc.append(0.0)          # carries no opacity ...
                    n_lower.append(0.0)        # ... so tau = 0, beta = 1
                    lower.append(sink + k)
                A.append(a_each)
                upper.append(i + 1)            # SAME upper level -> branching
        n_levels = sink + n_sink_eff
    else:
        n_levels = n_lines + 1

    atom = ForestAtom(nu0=np.array(nu0), f_osc=np.array(f_osc),
                      n_lower=np.array(n_lower), n_upper=np.zeros(len(nu0)),
                      A=np.array(A), lower=np.array(lower, int),
                      upper=np.array(upper, int), t_exp=t_exp,
                      tau_min=1e-3, stim=False)
    atom.emis_w = np.ones(len(nu0))
    atom.temperature = 3000.0
    info = dict(kind="dialled", n_lines=n_lines, tau=float(np.median(tau_abs)),
                tau_spread=tau_spread, span=span, n_exit=n_exit,
                dlnlam=dlnlam, f_return=f_return, jitter=jitter, seed=seed,
                n_levels=n_levels, delocalize=delocalize,
                exit_tau=exit_tau, n_sink=n_sink)
    return atom, info


def synthetic_ladder(n_lines=200, tau=5.0, tau_spread=0.0, span=0.2,
                     n_rungs=4, seed=0, t_exp=T_EXP, lam_mid=LAM_MID):
    """A forest whose redistribution EMERGES from a level ladder.

    Levels are a ladder of `n_rungs` energies; every absorbing line pumps from
    the ground level to a randomly chosen upper rung, and each upper rung
    decays to every lower rung with A ~ nu^3 (the Einstein scaling), so the exit
    wavelengths are set by the ladder geometry rather than prescribed. This is
    the control for `synthetic_forest`: matched forests should land on the same
    point of the phase diagram once their *measured* redistribution range
    agrees, and if they do not, the dial is not physical.
    """
    rng = np.random.default_rng(seed)
    u = np.linspace(-0.5, 0.5, n_lines) if n_lines > 1 else np.zeros(1)
    lam_abs = np.sort(lam_mid * np.exp(u * span))
    nu_abs = C / lam_abs
    tau_abs = _tau_draw(n_lines, tau, tau_spread, rng)

    # ladder energies in cm^-1, spread over the forest's own photon energy
    e_top = float(np.max(nu_abs) / C)
    rungs = np.linspace(0.0, e_top, n_rungs + 1)      # rung 0 = ground

    nu0, f_osc, n_lower, A, lower, upper = [], [], [], [], [], []
    for i in range(n_lines):
        up = int(rng.integers(1, n_rungs + 1))
        nu0.append(nu_abs[i]); f_osc.append(F_OSC)
        n_lower.append(tau_abs[i] / tau_sobolev(F_OSC, 1.0, lam_abs[i], t_exp))
        A.append(1.0); lower.append(0); upper.append(up)
        for lo in range(up):                           # every lower rung
            dnu = (rungs[up] - rungs[lo]) * C
            if dnu <= 0:
                continue
            nu0.append(dnu); f_osc.append(0.0); n_lower.append(0.0)
            A.append((dnu / C) ** 3 * 1e-30 + 1e-12)   # Einstein nu^3
            lower.append(lo); upper.append(up)

    atom = ForestAtom(nu0=np.array(nu0), f_osc=np.array(f_osc),
                      n_lower=np.array(n_lower), n_upper=np.zeros(len(nu0)),
                      A=np.array(A), lower=np.array(lower, int),
                      upper=np.array(upper, int), t_exp=t_exp,
                      tau_min=1e-3, stim=False)
    atom.emis_w = np.ones(len(nu0))
    atom.temperature = 3000.0
    info = dict(kind="ladder", n_lines=n_lines, tau=float(np.median(tau_abs)),
                span=span, n_rungs=n_rungs, seed=seed)
    return atom, info
