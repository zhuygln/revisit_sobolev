"""The effective redistribution matrix R_ij (chapters 9-11) and the models
that use it, built from the macroatom's own statistics.

Groups: the atom's lines are sorted by frequency and cut into N_g
contiguous groups (with ten lines, N_g = 10 resolves every line and
N_g = 1 is a single group). From many interactions of the macroatom
transport, each an (absorbed line, emitted line) event, the matrix

    R_ij = N(i -> j) / sum_k N(i -> k)

is the fraction of the energy absorbed in group i that leaves in group j
(energy packets: every event carries the same energy). A transport that
uses R alone draws the exit group from row i and then a line within the
group from a fixed within-group weight, the thermal emissivity: that
within-group choice is exactly the information the binning discarded.
"""
import numpy as np

from .macroatom import ToyMacroAtom
from .sobolev import escape_probability


def group_of_line(nu, n_g):
    """Contiguous frequency groups: the index of each line's group, and the
    groups' frequency edges (for plotting)."""
    order = np.argsort(nu)
    n_g = int(min(n_g, nu.size))
    g = np.empty(nu.size, int)
    for j, chunk in enumerate(np.array_split(order, n_g)):
        g[chunk] = j
    return g, n_g


class MacroatomRedistribution:
    """The reference: an interaction in line k activates the macroatom at
    the line's upper level and the walk names the exit line; beta from the
    line list."""

    def __init__(self, atom, tau):
        self.atom = atom
        self.m = ToyMacroAtom(atom, beta=escape_probability(tau))
        self.events = []

    def __call__(self, rng, k):
        j, _ = self.m.walk(rng, int(self.atom.upper[k]))
        self.events.append((int(k), int(j)))
        return int(j)

    @property
    def n_parameters(self):
        return int(self.atom.n_levels + self.atom.n_lines)


def build_R(events, group, n_g):
    """R from (absorbed line, emitted line) events; rows without events are
    identity rows (nothing known: coherent)."""
    R = np.zeros((n_g, n_g))
    for k, j in events:
        R[group[k], group[j]] += 1.0
    rows = R.sum(axis=1)
    for i in range(n_g):
        if rows[i] > 0:
            R[i] /= rows[i]
        else:
            R[i, i] = 1.0
    return R


class MatrixRedistribution:
    """Transport with R alone: exit group from row i, exit line within the
    group from `within` (the thermal emissivity restricted to the group)."""

    def __init__(self, R, group, within):
        self.R = np.asarray(R, float); self.group = np.asarray(group, int)
        self.cum = np.cumsum(self.R, axis=1)
        self.within = np.asarray(within, float)
        self.lines_in = [np.flatnonzero(self.group == j) for j in range(self.R.shape[0])]
        self.within_cum = []
        for lines in self.lines_in:
            w = self.within[lines]
            w = w / w.sum() if w.sum() > 0 else np.ones(lines.size) / max(lines.size, 1)
            self.within_cum.append(np.cumsum(w))

    def __call__(self, rng, k):
        i = self.group[k]
        j = int(np.searchsorted(self.cum[i], rng.random())); j = min(j, self.R.shape[0] - 1)
        lines = self.lines_in[j]
        if lines.size == 0:
            return int(k)
        return int(lines[min(int(np.searchsorted(self.within_cum[j], rng.random())), lines.size - 1)])

    @property
    def n_parameters(self):
        return int(self.R.size)


def low_rank(R, rank):
    """A rank-`rank` approximation of R by truncated SVD, clipped to be
    non-negative and renormalised row by row (a stand-in for a constrained
    factorisation; the teaching point is rank against frequency locality)."""
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    A = (U[:, :rank] * s[:rank]) @ Vt[:rank]
    A = np.clip(A, 0.0, None)
    rows = A.sum(axis=1, keepdims=True)
    return np.where(rows > 0, A / np.where(rows > 0, rows, 1.0), np.eye(R.shape[0]))


def interpolate_R(T, Ts, Rs):
    """Row-wise linear interpolation in log T between tabulated matrices,
    renormalised; clamped at the ends of the table."""
    Ts = np.asarray(Ts, float); x = np.log(T); xs = np.log(Ts)
    if x <= xs[0]:
        return np.array(Rs[0])
    if x >= xs[-1]:
        return np.array(Rs[-1])
    k = int(np.searchsorted(xs, x)) - 1
    w = (x - xs[k]) / (xs[k + 1] - xs[k])
    R = (1 - w) * np.asarray(Rs[k]) + w * np.asarray(Rs[k + 1])
    return R / R.sum(axis=1, keepdims=True)


def row_error(R_model, R_true, weights=None):
    """Mean absolute row error, weighted by how often each row is used."""
    d = np.abs(np.asarray(R_model) - np.asarray(R_true)).sum(axis=1)
    if weights is None:
        return float(d.mean())
    w = np.asarray(weights, float); return float((d * w).sum() / w.sum())
