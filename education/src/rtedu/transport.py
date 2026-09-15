"""The radial Sobolev sweep with a redistribution hook: the toy transport of
chapters 5-12. A packet launched at r0 moves radially outward through the
five-level atom's line forest; at each resonance it flips chapter 4's coin;
at an interaction the redistribution model names the exit line and the
packet continues outward from that radius with the exit line's comoving
frequency (its lab frequency is nu_com / (1 - r/ct)). The energy of an
energy packet is unchanged by the interaction (the first-order toy ignores
the expansion work that the production code books as W).

It is one dimensional on purpose: angular redistribution and the trapping
it causes are chapter 13's.
"""
import numpy as np

from . import C


def sweep(rng, nu_launch, nu_lines, tau, r_out, t, redistribute, max_interactions=200):
    """One packet. Returns (nu_lab at escape, number of interactions, list
    of (line absorbed, line emitted))."""
    r = 0.0; nu_lab = float(nu_launch); history = []
    for _ in range(max_interactions):
        nu_now = nu_lab * (1.0 - r / (C * t))
        r_res = C * t * (1.0 - nu_lines / nu_lab)
        ok = (nu_lines < nu_now) & (r_res < r_out)
        if not ok.any():
            return nu_lab, len(history), history
        order = np.argsort(-nu_lines[ok]); hit = None
        for k in np.flatnonzero(ok)[order]:
            if rng.random() < 1.0 - np.exp(-tau[k]):
                hit = int(k); break
        if hit is None:
            return nu_lab, len(history), history
        j = redistribute(rng, hit)
        history.append((hit, j))
        r = float(r_res[hit])
        nu_lab = nu_lines[j] / (1.0 - r / (C * t))       # the exit line's comoving frequency, seen in the lab
    return nu_lab, len(history), history


def run(rng, nu_launch, n, nu_lines, tau, r_out, t, redistribute):
    """n packets of one launch frequency: escape frequencies, the line each
    packet last emitted in (-1 if never interacted), interaction counts."""
    nu = np.empty(n); last = np.full(n, -1); n_int = np.empty(n, int)
    for i in range(n):
        nu[i], n_int[i], h = sweep(rng, nu_launch, nu_lines, tau, r_out, t, redistribute)
        if h:
            last[i] = h[-1][1]
    return nu, last, n_int


def emergent_by_line(last, n_lines):
    """Fraction of packets that escaped after last emitting in each line
    (index n_lines = never interacted)."""
    out = np.zeros(n_lines + 1)
    np.add.at(out, np.where(last < 0, n_lines, last), 1.0)
    return out / last.size
