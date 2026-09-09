"""The downward macroatom: energy-conserving downward fluorescence with fixed
populations (Paper IV, Phase 2).

WHAT IT IS. Lucy's (2002, 2003) macroatom formalism makes indivisible energy
packets flow through an atom's level structure: a packet absorbed in line
l -> u activates the macroatom in level u, which then either de-activates
radiatively through one of its downward lines, carrying the WHOLE packet
energy out at that line's frequency, or makes an internal jump to another
level and tries again. The probabilities are chosen so that, in statistical
equilibrium, the emitted energy in every line equals its net emissivity.
Here the populations are frozen (LTE Boltzmann at T_gas, as in every Paper
II/III leg), the only processes are radiative bound-bound, and the Sobolev
escape probability beta_ul is the net radiative rate factor. Activated in
level i with energy eps_i above the ion's ground:

    de-activate via i -> j :  p  ~  A_ij beta_ij (eps_i - eps_j)   emit E at nu_ij
    internal jump  i -> j  :  p  ~  A_ij beta_ij  eps_j            continue from j

normalised over all downward lines of i. The walk terminates because the
ground has eps = 0 (no internal jump lands there) and every level's energy
is strictly positive.

WHAT IT IS NOT. A full Lucy macroatom also has internal UPWARD transitions,
driven by the radiation field J_nu through B_ji J_bar (and collisions).
Those need a mean-intensity estimator that a fixed-atmosphere study does not
have; an imposed J_bar would contradict the program's premise. They arrive
with the radiative-equilibrium iteration (Paper IV WP8-10), when this class
grows into a `MacroAtom`. Until then every leg built on this table is
"energy-conserving downward fluorescence" -- TARDIS's `downbranch` with
internal downward jumps, closer to but not the same as its `macroatom` --
and is named `dmacro` in `forest_mc.MODES` for that reason.

BETA ONCE. The escape probability is inside these weights, exactly as in the
closed-form exit kernel A beta / sum A beta of `expansion_branch` (Paper II
E8). A packet that de-activates has already paid its re-absorption odds;
it leaves the resonance and must NOT be sent through the beta escape /
re-absorption chain again (that would give exits ~ A beta^2). The chain was
rejection sampling of the A table with acceptance beta; this table is its
closed form, with zero chain trapping.

OPTIONAL CHANNELS (both off by default):
  * thermal de-activation with probability `eps_k` per activation -- the
    packet becomes a k-packet (Test C of the plan): booked as deposit, or
    re-emitted from the energy-weighted net line emissivity A n_u h nu beta
    with the same energy;
  * a table cut `a_cut`: per level, downward lines are kept in decreasing
    weight until 1 - a_cut of the level's total A beta eps_i is covered (the
    rest is renormalised away). Cut lines still carry opacity if they have
    it; only their role as exits is dropped. Measured, not assumed.

Dead ends -- a level with eps > 0 and no downward line in the data -- cannot
de-activate radiatively; they are k-packets and counted.
"""
import numpy as np

from .constants import C, H

HC = H * C


class DownwardMacroAtom:
    """CSR tables of the downward macroatom, one segment per level.

    Parameters
    ----------
    nu0, A, lower, upper : per-line arrays (all lines of the atom).
    beta : per-line Sobolev escape probability (1 where a line carries no
        opacity).
    level_energy_cm : per-level excitation energy above the ion's ground,
        cm^-1 (a blend concatenates ions; each ion's ground is 0).
    a_cut : optional table cut, see the module docstring.
    """

    def __init__(self, nu0, A, lower, upper, beta, level_energy_cm, a_cut=None):
        nu0 = np.asarray(nu0, float); A = np.asarray(A, float)
        lower = np.asarray(lower, int); upper = np.asarray(upper, int)
        beta = np.asarray(beta, float)
        eps = HC * np.asarray(level_energy_cm, float)
        n_lev = eps.size
        ab = A * beta
        good = (ab > 0) & (nu0 > 0)
        # two entries per line: de-activation (kind 0) and internal jump (kind 1)
        w_deact = np.where(good, ab * H * nu0, 0.0)
        w_int = np.where(good, ab * eps[lower], 0.0)
        line2 = np.concatenate([np.flatnonzero(good), np.flatnonzero(good)])
        kind2 = np.concatenate([np.zeros(good.sum(), np.int8), np.ones(good.sum(), np.int8)])
        w2 = np.concatenate([w_deact[good], w_int[good]])
        lev2 = upper[line2]
        # group by level (stable so entries stay in a deterministic order)
        order = np.argsort(lev2, kind="stable")
        line2, kind2, w2, lev2 = line2[order], kind2[order], w2[order], lev2[order]
        off = np.searchsorted(lev2, np.arange(n_lev + 1))
        keep = np.ones(line2.size, bool)
        if a_cut is not None and a_cut > 0:
            for i in range(n_lev):
                a, b = off[i], off[i + 1]
                if b <= a:
                    continue
                ws = w2[a:b]
                o = np.argsort(-ws, kind="stable")
                cum = np.cumsum(ws[o]) / ws.sum()
                n_keep = int(np.searchsorted(cum, 1.0 - a_cut, side="left")) + 1
                drop = o[n_keep:]
                keep[a + drop] = False
            line2, kind2, w2, lev2 = line2[keep], kind2[keep], w2[keep], lev2[keep]
            off = np.searchsorted(lev2, np.arange(n_lev + 1))
        cum = np.empty(w2.size)
        for i in range(n_lev):
            a, b = off[i], off[i + 1]
            if b > a:
                c = np.cumsum(w2[a:b]); cum[a:b] = c / c[-1]
        self.n_levels = n_lev
        self.off = off
        self.line = line2
        self.kind = kind2
        self.target = lower[line2]           # level reached by an internal jump
        self.cum = cum
        self.key = np.repeat(np.arange(n_lev, dtype=float), np.diff(off)) + cum
        self.eps = eps
        self.has_exit = np.diff(off) > 0
        self.a_cut = a_cut
        self.n_entries = int(line2.size)
        # levels with energy above ground but no downward line: dead ends
        self.dead_end = (~self.has_exit) & (eps > 0)

    # ---- sampling ------------------------------------------------------
    def sample(self, levels, v):
        """One table entry per packet: the first entry in the level's segment
        whose cumulative weight is >= v (searchsorted 'left' semantics, exact
        -- see `ForestAtom.sample_branch` for the rounding argument).
        Levels without a segment must be masked by the caller (`has_exit`)."""
        levels = np.asarray(levels, int); v = np.asarray(v, float)
        lo = self.off[levels]; hi = self.off[levels + 1] - 1
        i = np.clip(np.searchsorted(self.key, levels + v, side="left"), lo, hi)
        for _ in range(8):
            down = (i > lo) & (self.cum[np.maximum(i - 1, 0)] >= v)
            i = np.where(down, i - 1, i)
            up = (i < hi) & (self.cum[i] < v)
            i = np.where(up, i + 1, i)
            if not (down.any() or up.any()):
                break
        return i

    def walk(self, levels, rng, max_jumps=10000):
        """Run the downward walk for every packet from `levels`. Returns
        (exit_line, n_jumps, dead) with exit_line = -1 where the walk ended in
        a dead end. One uniform per step; the number of draws depends on the
        walk, exactly as the branching chain's did."""
        lev = np.array(levels, int, copy=True)
        exit_line = np.full(lev.size, -1, np.int64)
        n_jumps = np.zeros(lev.size, np.int32)
        dead = np.zeros(lev.size, bool)
        active = np.flatnonzero(np.ones(lev.size, bool))
        for _ in range(max_jumps):
            if active.size == 0:
                break
            ok = self.has_exit[lev[active]]
            dead[active[~ok]] = True
            active = active[ok]
            if active.size == 0:
                break
            e = self.sample(lev[active], rng.uniform(size=active.size))
            de = self.kind[e] == 0
            exit_line[active[de]] = self.line[e[de]]
            jump = active[~de]
            lev[jump] = self.target[e[~de]]
            n_jumps[jump] += 1
            active = jump
        else:
            raise RuntimeError("downward macroatom walk did not terminate")
        return exit_line, n_jumps, dead

    # ---- diagnostics ---------------------------------------------------
    def probabilities(self, level):
        """(lines, p_deactivate, p_internal) for one level -- analytic
        fractions for the toy tests and for a TARDIS table comparison."""
        a, b = self.off[level], self.off[level + 1]
        seg_cum = self.cum[a:b]
        p = np.diff(np.concatenate([[0.0], seg_cum]))
        lines, kinds = self.line[a:b], self.kind[a:b]
        ul = np.unique(lines)
        p_de = np.array([p[(lines == l) & (kinds == 0)].sum() for l in ul])
        p_in = np.array([p[(lines == l) & (kinds == 1)].sum() for l in ul])
        return ul, p_de, p_in
