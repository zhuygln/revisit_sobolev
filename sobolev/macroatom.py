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


class MacroAtom(DownwardMacroAtom):
    """Lucy's (2003) macroatom with internal UPWARD transitions driven by an
    imposed radiation field, on top of the downward tables (Paper IV,
    Phase 10's "one stronger macroatom check").

    Activated in level i (energy eps_i above the ion's ground), with the
    radiative rates per atom in level i:

        de-activate via i -> j (j < i):  p ~ A_ij beta_ij (eps_i - eps_j)
        internal down  i -> j (j < i):   p ~ A_ij beta_ij  eps_j
        internal up    i -> k (k > i):   p ~ B_ik Jbar_ik  eps_i

    B_ik Jbar_ik = (g_k / g_i) A_ki W / (exp(h nu / k T_rad) - 1) for the
    imposed field Jbar = W B_nu(T_rad): a diluted Planck field at the local
    radiation temperature, the fixed-state stand-in for the estimator that a
    radiative-equilibrium iteration would provide (W = 1/2 at a photosphere,
    smaller above it). Collisions are not included. Downward weights are the
    net radiative rates with the Sobolev beta, as in `DownwardMacroAtom`;
    the upward rate uses the imposed field without a beta factor (the
    absorbed field is the external one, not the line's own trapped
    radiation -- stated, not derived). The walk is a Markov chain with
    de-activation as the absorbing state; it terminates with probability
    one and is capped at `max_jumps` (an overflow is booked as a dead end
    and counted, never raised). A level with no downward line but with
    upward lines is no longer a dead end: it climbs and cascades elsewhere,
    which is the physical answer to the 7-10 % dead-end fraction of the
    downward table.

    Parameters as `DownwardMacroAtom` plus `g_lev` (statistical weights),
    `T_rad` (K) and `W` (dilution). `a_cut` is not supported here.
    """

    def __init__(self, nu0, A, lower, upper, beta, level_energy_cm, g_lev, T_rad, W=0.5):
        nu0 = np.asarray(nu0, float); A = np.asarray(A, float)
        lower = np.asarray(lower, int); upper = np.asarray(upper, int)
        beta = np.asarray(beta, float); g = np.asarray(g_lev, float)
        eps = HC * np.asarray(level_energy_cm, float)
        n_lev = eps.size
        ab = A * beta
        good = (ab > 0) & (nu0 > 0)
        gi = np.flatnonzero(good)
        gu = np.flatnonzero((A > 0) & (nu0 > 0))
        from .constants import K_B
        x = H * nu0[gu] / (K_B * float(T_rad))
        with np.errstate(over="ignore"):
            occ = np.where(x < 700.0, 1.0 / np.expm1(np.minimum(x, 700.0)), 0.0)
        w_up = (g[upper[gu]] / g[lower[gu]]) * A[gu] * float(W) * occ * eps[lower[gu]]
        line3 = np.concatenate([gi, gi, gu])
        kind3 = np.concatenate([np.zeros(gi.size, np.int8), np.ones(gi.size, np.int8), np.full(gu.size, 2, np.int8)])
        w3 = np.concatenate([ab[gi] * H * nu0[gi], ab[gi] * eps[lower[gi]], w_up])
        lev3 = np.concatenate([upper[gi], upper[gi], lower[gu]])
        tgt3 = np.concatenate([lower[gi], lower[gi], upper[gu]])
        keep = (kind3 < 2) | (w3 > 0)          # downward entries exactly as the downward table keeps them
        line3, kind3, w3, lev3, tgt3 = line3[keep], kind3[keep], w3[keep], lev3[keep], tgt3[keep]
        order = np.argsort(lev3, kind="stable")
        line3, kind3, w3, lev3, tgt3 = line3[order], kind3[order], w3[order], lev3[order], tgt3[order]
        off = np.searchsorted(lev3, np.arange(n_lev + 1))
        cum = np.empty(w3.size)
        for i in range(n_lev):
            a, b = off[i], off[i + 1]
            if b > a:
                c = np.cumsum(w3[a:b]); cum[a:b] = c / c[-1]
        self.n_levels = n_lev; self.off = off; self.line = line3; self.kind = kind3
        self.target = tgt3; self.cum = cum
        self.key = np.repeat(np.arange(n_lev, dtype=float), np.diff(off)) + cum
        self.eps = eps; self.has_exit = np.diff(off) > 0; self.a_cut = None
        self.n_entries = int(line3.size)
        self.dead_end = (~self.has_exit) & (eps > 0)
        self.T_rad, self.W = float(T_rad), float(W)
        self.n_up_entries = int((kind3 == 2).sum())
        self.n_overflow = 0

    def walk(self, levels, rng, max_jumps=100000):
        """The macroatom walk: de-activation ends it; internal jumps go down
        (kind 1) or up (kind 2). Overflowing `max_jumps` is booked as a dead
        end and counted in `n_overflow`."""
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
            dead[active] = True
            self.n_overflow += int(active.size)
        return exit_line, n_jumps, dead

    def probabilities(self, level):
        """(lines, p_deactivate, p_internal_down, p_internal_up) for one level."""
        a, b = self.off[level], self.off[level + 1]
        p = np.diff(np.concatenate([[0.0], self.cum[a:b]]))
        lines, kinds = self.line[a:b], self.kind[a:b]
        ul = np.unique(lines)
        return (ul, np.array([p[(lines == l) & (kinds == 0)].sum() for l in ul]),
                np.array([p[(lines == l) & (kinds == 1)].sum() for l in ul]),
                np.array([p[(lines == l) & (kinds == 2)].sum() for l in ul]))
