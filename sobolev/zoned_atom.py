"""Radially zoned atoms for multi-shell transport (Paper IV, Phases 8-10).

A single-zone `ForestAtom` holds one set of optical depths, one set of
populations and one macroatom table. A zoned state has N transported shells
with their own (rho, T, X_Z, f_ion); the line list is the same in every
shell, and what changes per shell is exactly what depends on populations:
the Sobolev optical depths (hence the interaction probabilities and the
escape probabilities beta), the emissivities, and the downward-macroatom
weights A beta eps through beta. This module holds those per-shell arrays
behind one object so that `run_mc` can be zoned with thin hooks, and so
that N = 1 reproduces `ForestAtom` bit for bit.

Measured on P1 at 2 d (13 La-Yb II ions, tau_min = 1e-3): the opacity-line
sets of the shells NEST (the union over shells 26-31 is shell 26's set), so
the union is bounded by the densest transported shell; 71 % of the 39.7 M
macroatom entries sit in levels with at least one union opacity line and
are shell-dependent through beta, 29 % are shell-independent and shared.
A per-shell k-packet sampler needs 2.9 M lines for 1 - 1e-6 of its weight.

MEMORY. The first version kept every per-shell array over all 20 M lines
(populations, tau, beta, emissivity) and was killed at 23 GB on twelve
shells. This one streams: shell by shell it computes the all-line
transients, keeps only the union-subset opacity arrays, the shell's
macroatom block and its two prebuilt emissivity samplers, and frees the
rest. The macroatom's shell-dependent levels are those with at least one
union opacity line (decided before any beta exists), so identical shells
hold identical blocks rather than a shared one.

Two classes:

`ZonedMacroAtom`  the downward-macroatom tables of `sobolev/macroatom.py`
                  with a shared entry order and per-shell cumulative
                  weights for the beta-dependent level segments only;
                  sampling by vectorised bisection, which is exactly the
                  searchsorted-left semantics `DownwardMacroAtom.sample`
                  reproduces.
`ZonedAtom`       the line list, the union opacity subset with per-shell
                  tau / p / beta and a per-shell skip table, per-shell
                  emissivity samplers with a weight cut, per-shell bins.
"""
import numpy as np

from .constants import C, H, SIGMA_CLASSICAL, K_B
from .macroatom import HC

# --------------------------------------------------------------------------
# the macroatom tables
# --------------------------------------------------------------------------


def _segment_cum(w, off, n_lev, out=None):
    """cum = cumsum(w)/sum(w) per level segment, with the same per-level loop
    as `DownwardMacroAtom` (a vectorised cumsum rounds differently)."""
    cum = np.empty(w.size) if out is None else out
    for i in range(n_lev):
        a, b = off[i], off[i + 1]
        if b > a:
            c = np.cumsum(w[a:b]); cum[a:b] = c / c[-1]
    return cum


class ZonedMacroAtom:
    """Per-shell downward-macroatom tables with a shared entry order.

    Parameters
    ----------
    nu0, A, lower, upper : per-line arrays (all lines).
    level_energy_cm : per-level energies above the ion's ground.
    n_shell : number of shells.
    dep_lines : bool per line -- lines whose beta may differ between shells
        (the union opacity lines); every other line has beta = 1 in every
        shell. Levels with any such line get one block per shell, filled
        by `add_shell(s, beta_s)`; the others share the base block.
    """

    def __init__(self, nu0, A, lower, upper, level_energy_cm, n_shell, dep_lines):
        nu0 = np.asarray(nu0, float); A = np.asarray(A, float)
        lower = np.asarray(lower, int); upper = np.asarray(upper, int)
        eps = HC * np.asarray(level_energy_cm, float)
        n_lev = eps.size
        good = (A > 0) & (nu0 > 0)
        gi = np.flatnonzero(good)
        line2 = np.concatenate([gi, gi])
        kind2 = np.concatenate([np.zeros(gi.size, np.int8), np.ones(gi.size, np.int8)])
        lev2 = upper[line2]
        order = np.argsort(lev2, kind="stable")
        line2, kind2, lev2 = line2[order], kind2[order], lev2[order]
        off = np.searchsorted(lev2, np.arange(n_lev + 1))
        self.n_levels = n_lev; self.off = off; self.line = line2.astype(np.int32)
        self.kind = kind2; self.target = lower[line2].astype(np.int32)
        self.eps = eps
        seg_len = np.diff(off)
        self.seg_len = seg_len
        self.has_exit = seg_len > 0
        self.dead_end = (~self.has_exit) & (eps > 0)
        self.n_entries = int(line2.size)
        self.n_shell = int(n_shell)
        # weights exactly as DownwardMacroAtom forms them: (A beta) * H * nu0 for
        # de-activation, (A beta) * eps_lower for internal jumps -- the
        # operation order matters for bit-identity at one shell
        self._nu_line = nu0[line2]; self._eps_low = eps[lower[line2]]; self._a_line = A[line2]
        self._deact = kind2 == 0
        dep = np.zeros(n_lev, bool)
        dl = np.asarray(dep_lines, bool)
        np.logical_or.at(dep, lev2[dl[line2]], True)
        self.dep = dep
        self._dep_entries = dep[lev2]
        # base block: beta = 1 for every line of an independent level
        base = _segment_cum(self._weights(np.ones(nu0.size)), off, n_lev)[~self._dep_entries]
        base_off = np.zeros(n_lev + 1, np.int64); base_off[1:] = np.cumsum(np.where(dep, 0, seg_len))
        dep_off = np.zeros(n_lev + 1, np.int64); dep_off[1:] = np.cumsum(np.where(dep, seg_len, 0))
        n_dep = int(dep_off[-1])
        self.n_dep_entries = n_dep
        self._base_size = int(base.size)
        self.cum = np.empty(base.size + self.n_shell * n_dep)
        self.cum[:base.size] = base
        self.start = np.empty((self.n_shell, n_lev), np.int64)
        for s in range(self.n_shell):
            self.start[s] = np.where(dep, base.size + s * n_dep + dep_off[:-1], base_off[:-1])
        self.filled = np.zeros(self.n_shell, bool)

    def _weights(self, beta_all):
        ab = self._a_line * beta_all[self.line]
        return np.where(self._deact, ab * H * self._nu_line, ab * self._eps_low)

    def add_shell(self, s, beta_all):
        """Fill shell s's block from its beta over ALL lines (a transient
        of the caller; nothing over all lines is kept here)."""
        cum = _segment_cum(self._weights(np.asarray(beta_all, float)), self.off, self.n_levels)
        a = self._base_size + s * self.n_dep_entries
        self.cum[a:a + self.n_dep_entries] = cum[self._dep_entries]
        self.filled[s] = True

    @classmethod
    def from_beta(cls, nu0, A, lower, upper, beta, level_energy_cm):
        """Toys and tests: beta as (n_shell, n_lines) (1-D = one shell);
        lines with beta < 1 in any shell are the dependent ones."""
        beta = np.asarray(beta, float)
        if beta.ndim == 1:
            beta = beta[None, :]
        zm = cls(nu0, A, lower, upper, level_energy_cm, beta.shape[0], np.any(beta < 1.0, axis=0))
        for s in range(beta.shape[0]):
            zm.add_shell(s, beta[s])
        return zm

    def level_cum(self, level, shell=0):
        """The level's cumulative segment in that shell (the slice of
        `DownwardMacroAtom.cum` for the same level at one shell)."""
        st = self.start[shell, level]
        return self.cum[st:st + self.seg_len[level]]

    # ---- sampling ------------------------------------------------------
    def sample(self, levels, v, shells):
        """Entry index (into line/kind/target) of the first entry of the
        level's segment in that shell whose cumulative weight is >= v:
        searchsorted-left by bisection, exact."""
        levels = np.asarray(levels, int); v = np.asarray(v, float); shells = np.asarray(shells, int)
        n = self.seg_len[levels]
        base = self.start[shells, levels]
        lo = np.zeros(levels.size, np.int64); hi = n.copy()          # search in [lo, hi)
        for _ in range(int(np.ceil(np.log2(max(int(n.max()), 1)))) + 1):
            act = lo < hi
            if not act.any():
                break
            mid = (lo + hi) // 2
            less = np.zeros(levels.size, bool)
            less[act] = self.cum[base[act] + mid[act]] < v[act]
            lo = np.where(act & less, mid + 1, lo)
            hi = np.where(act & ~less, mid, hi)
        j = np.minimum(lo, n - 1)              # v beyond the last cum (roundoff): last entry
        return self.off[levels] + j

    def walk(self, levels, rng, shells, max_jumps=10000):
        """The downward walk per packet from `levels` in `shells`; one uniform
        per active walker per step (exactly `DownwardMacroAtom.walk`)."""
        lev = np.array(levels, int, copy=True); shells = np.asarray(shells, int)
        exit_line = np.full(lev.size, -1, np.int64)
        n_jumps = np.zeros(lev.size, np.int32)
        dead = np.zeros(lev.size, bool)
        active = np.arange(lev.size)
        for _ in range(max_jumps):
            if active.size == 0:
                break
            ok = self.has_exit[lev[active]]
            dead[active[~ok]] = True
            active = active[ok]
            if active.size == 0:
                break
            e = self.sample(lev[active], rng.uniform(size=active.size), shells[active])
            de = self.kind[e] == 0
            exit_line[active[de]] = self.line[e[de]]
            jump = active[~de]
            lev[jump] = self.target[e[~de]]
            n_jumps[jump] += 1
            active = jump
        else:
            raise RuntimeError("downward macroatom walk did not terminate")
        return exit_line, n_jumps, dead

    def probabilities(self, level, shell=0):
        a, b = self.off[level], self.off[level + 1]
        st = self.start[shell, level]
        seg = self.cum[st:st + (b - a)]
        p = np.diff(np.concatenate([[0.0], seg]))
        lines, kinds = self.line[a:b], self.kind[a:b]
        ul = np.unique(lines)
        return (ul, np.array([p[(lines == l) & (kinds == 0)].sum() for l in ul]),
                np.array([p[(lines == l) & (kinds == 1)].sum() for l in ul]))


# --------------------------------------------------------------------------
# the zoned atom
# --------------------------------------------------------------------------


def _beta_of(tau):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(tau > 1e-12, -np.expm1(-tau) / np.where(tau > 1e-12, tau, 1.0), 1.0)


def _cut_sampler(w, emis_cut):
    """u -> line index from the weights w (all lines); lines below the
    weight cut are dropped from the table (None = full, bit-identical)."""
    tot = w.sum()
    if tot <= 0:
        return None                                   # no emissivity (toys): error on use
    cum_full = np.cumsum(w / tot)
    if emis_cut is None:
        idx = np.arange(w.size); cum = cum_full
    else:
        # keep the lines that carry all but emis_cut of the weight, in order
        order = np.argsort(-w, kind="stable")
        c = np.cumsum(w[order]) / tot
        n_keep = int(np.searchsorted(c, 1.0 - emis_cut, side="left")) + 1
        idx = np.sort(order[:n_keep])
        cum = cum_full[idx]
    def sample(u):
        return idx[np.minimum(np.searchsorted(cum, u), idx.size - 1)]
    sample.n_lines = int(idx.size)
    return sample


class ZonedAtom:
    """A line list with per-shell populations for `run_mc`.

    Built from per-shell lower/upper populations of every line
    (`from_arrays`: (n_shell, n_lines) arrays, or a callable
    `populations(s) -> (n_lower_s, n_upper_s)` with `n_shell`, which is how
    `from_state` streams one shell at a time) or from an `EjectaState` and
    the compact cache (`from_state`). Public arrays mirror `ForestAtom`'s
    where the transport reads them, with a leading shell axis where they
    depend on the shell. Nothing over all lines is kept per shell.
    """
    is_zoned = True

    def __init__(self, nu0, f_osc, n_lower, n_upper, A, lower, upper, r_edges, t_exp,
                 tau_min=1e-3, stim=True, temperature=None, level_energy_cm=None,
                 emis_cut=1e-6, ions=None, ion_of_line=None, ion_of_level=None, n_shell=None):
        nu0 = np.asarray(nu0, float); f_osc = np.asarray(f_osc, float)
        A = np.asarray(A, float); lower = np.asarray(lower, int); upper = np.asarray(upper, int)
        if callable(n_lower):
            if n_shell is None:
                raise ValueError("n_shell is needed with a populations callable")
            pops = n_lower
        else:
            nl = np.atleast_2d(np.asarray(n_lower, float)); nu_ = np.atleast_2d(np.asarray(n_upper, float))
            n_shell = nl.shape[0]
            pops = lambda s: (nl[s], nu_[s])
        self.n_shell = int(n_shell)
        self.r_edges = np.asarray(r_edges, float)
        if self.r_edges.size != self.n_shell + 1:
            raise ValueError(f"r_edges has {self.r_edges.size} entries for {self.n_shell} shells")
        self.t_exp = float(t_exp)
        self.v_edges = self.r_edges / self.t_exp
        T = None if temperature is None else np.broadcast_to(np.asarray(temperature, float), (self.n_shell,)).copy()
        self.T = T
        self.temperature = None if T is None else float(T[0])     # ForestAtom compatibility (bins)
        self.n_lines_total = nu0.size
        self.nu0_all, self.A_all, self.lower_all, self.upper_all = nu0, A, lower, upper
        self.level_energy_cm = None if level_energy_cm is None else np.asarray(level_energy_cm, float)
        self.ions, self.ion_of_line, self.ion_of_level = ions, ion_of_line, ion_of_level
        self.tau_min, self.emis_cut = tau_min, emis_cut
        self._pops = pops; self._f_osc = f_osc
        self._stim = bool(stim and T is not None)

        # pass 1: the union of the shells' opacity lines (one transient at a time)
        keep = np.zeros(nu0.size, bool)
        for s in range(self.n_shell):
            keep |= self._tau_of(s) > tau_min
        op = np.flatnonzero(keep)
        order = op[np.argsort(nu0[op])]
        self.op_idx = order
        self.op_nu = nu0[order]
        self.op_upper = upper[order]
        self.n_opacity = order.size
        U = order.size
        self.line_to_union = np.full(nu0.size, -1, np.int64); self.line_to_union[order] = np.arange(U)
        self.op_tau = np.empty((self.n_shell, U)); self.op_p = np.empty((self.n_shell, U))
        self.op_beta = np.empty((self.n_shell, U))
        self.op_nxt = np.empty((self.n_shell, U), np.int64)
        self.n_opacity_shell = np.zeros(self.n_shell, np.int64)
        self.macro = None
        if self.level_energy_cm is not None:
            self.macro = ZonedMacroAtom(nu0, A, lower, upper, self.level_energy_cm, self.n_shell, keep)
        self._samplers = {}
        self.tau_all = None; self.beta_all = None                # single-shell compatibility
        # pass 2: shell by shell; tau, beta and the emissivity are transients
        for s in range(self.n_shell):
            tau = self._tau_of(s)
            t_op = tau[order]
            live = t_op > tau_min
            self.op_tau[s] = np.where(live, t_op, 0.0)
            self.op_p[s] = np.where(live, -np.expm1(-t_op), 0.0)
            beta = _beta_of(tau)
            free = np.ones(nu0.size, bool); free[order[live]] = False
            beta = np.where(free, 1.0, beta)
            self.op_beta[s] = beta[order]
            self.n_opacity_shell[s] = int(live.sum())
            # skip table: largest k' <= k with a live line in this shell, else -1
            idx = np.where(live, np.arange(U), -1)
            self.op_nxt[s] = np.maximum.accumulate(idx) if U else idx
            if self.macro is not None:
                self.macro.add_shell(s, beta)
            w = A * pops(s)[1] * nu0                           # A n_u nu (h dropped, as ForestAtom)
            self._samplers[(s, "energy", None)] = _cut_sampler(w, emis_cut)
            self._samplers[(s, "energy_beta", None)] = _cut_sampler(w * beta, emis_cut)
            if self.n_shell == 1:
                self.tau_all = tau; self.beta_all = beta
            del tau, beta, w

    def _tau_of(self, s):
        """Sobolev optical depth of every line in shell s (transient)."""
        nu0 = self.nu0_all
        tau = SIGMA_CLASSICAL * self._f_osc * self._pops(s)[0] * (C / nu0) * self.t_exp
        if self._stim:
            tau = tau * (1.0 - np.exp(-H * nu0 / (K_B * self.T[s])))
        return tau

    def beta_of_lines(self, shells, lines):
        """Escape probability of `lines` for packets in `shells`: 1 outside
        the union opacity set, the shell's beta inside it."""
        lu = self.line_to_union[np.asarray(lines, int)]
        return np.where(lu >= 0, self.op_beta[np.asarray(shells, int), np.maximum(lu, 0)], 1.0)

    # ---- builders ------------------------------------------------------
    @classmethod
    def from_arrays(cls, nu0, f_osc, n_lower, n_upper, A, lower, upper, r_edges, t_exp, **kw):
        return cls(nu0, f_osc, n_lower, n_upper, A, lower, upper, r_edges, t_exp, **kw)

    @classmethod
    def from_state(cls, state, shells, stages=("II",), tau_min=1e-3, emis_cut=1e-6, cache_dir=None,
                   n_ion_min=1e-30, f_min=None):
        """The zoned blend of the transported `shells` of an `EjectaState`,
        from the compact cache; populations Boltzmann at each shell's T,
        computed per shell on demand (nothing over all lines is held).
        `f_min`: drop lines with oscillator strength below it before anything
        else (Fontes et al. 2020's f_c cut; None keeps every line)."""
        from .abundances import ATOMIC_MASS, Z_OF, GSI_IONS
        from .atomic_cache import load_cached
        from .populations import boltzmann_fractions
        shells = list(shells)
        kw = {} if cache_dir is None else {"cache_dir": cache_dir}
        elements = [el for el in state.X if el != "bulk" and el in Z_OF]
        specs = []
        for el in elements:
            for st in stages:
                if f"{Z_OF[el]}{el}{st}" not in GSI_IONS:
                    continue
                if f"{el} {st}" in state.f_ion or (st == "II" and not state.f_ion):
                    n = np.array([state.n_ion(el, st, s, ATOMIC_MASS[el]) for s in shells])
                    if n.max() > n_ion_min:
                        specs.append((f"{Z_OF[el]}{el}{st}", n))
        T = np.array([float(state.T_gas[s]) for s in shells])
        ions_d = [load_cached(ion, **kw) for ion, _ in specs]
        if f_min is not None:
            cut = []
            for d in ions_d:
                keep = d["f_lu"] > f_min
                cut.append({**d, "nu0": d["nu0"][keep], "f_lu": d["f_lu"][keep], "A": d["A"][keep],
                            "lower": d["lower"][keep], "upper": d["upper"][keep], "n_lines": int(keep.sum())})
            ions_d = cut
        offs = np.cumsum([0] + [d["n_levels"] for d in ions_d])
        cat = dict(nu0=np.concatenate([d["nu0"] for d in ions_d]), f=np.concatenate([d["f_lu"] for d in ions_d]),
                   A=np.concatenate([d["A"] for d in ions_d]), E=np.concatenate([d["E_lev"] for d in ions_d]),
                   low=np.concatenate([d["lower"].astype(int) + offs[i] for i, d in enumerate(ions_d)]),
                   up=np.concatenate([d["upper"].astype(int) + offs[i] for i, d in enumerate(ions_d)]))
        ion_line = np.concatenate([np.full(d["n_lines"], i) for i, d in enumerate(ions_d)])
        ion_lev = np.concatenate([np.full(d["n_levels"], i) for i, d in enumerate(ions_d)])
        n_per = [n for _, n in specs]

        def populations(k):
            """Boltzmann populations of every line's lower and upper level
            in transported shell k (a transient of the zoned builder)."""
            nl, nu_ = [], []
            for i, d in enumerate(ions_d):
                fr = boltzmann_fractions(d["g_lev"], d["E_lev"], T[k])
                nl.append(fr[d["lower"].astype(int)] * n_per[i][k]); nu_.append(fr[d["upper"].astype(int)] * n_per[i][k])
            return np.concatenate(nl), np.concatenate(nu_)

        r_edges = np.array([state.r_edges[s] for s in shells] + [state.r_edges[shells[-1] + 1]])
        atom = cls(cat["nu0"], cat["f"], populations, None, cat["A"], cat["low"], cat["up"], r_edges,
                   float(state.t), tau_min=tau_min, stim=True, temperature=T, level_energy_cm=cat["E"],
                   emis_cut=emis_cut, ions=[sp[0] for sp in specs], ion_of_line=ion_line, ion_of_level=ion_lev,
                   n_shell=len(shells))
        atom.shells = shells
        atom.f_min = f_min
        atom.n_ion = {sp[0]: sp[1].tolist() for sp in specs}
        atom.rho = np.array([float(state.rho[s]) for s in shells])
        return atom

    # ---- samplers ------------------------------------------------------
    def thermal_sampler(self, shell, emit_window=None, weight="energy_beta"):
        """u -> line index from the shell's LTE line emissivity: "energy"
        (A n_u h nu) or "energy_beta" (A n_u h nu beta); lines below the
        weight cut `emis_cut` are dropped from the table (None = full).
        The windowless samplers are prebuilt; a window rebuilds the shell's
        transients once."""
        if weight not in ("energy", "energy_beta"):
            raise ValueError(f"zoned samplers are energy-weighted: 'energy' or 'energy_beta', got {weight!r}")
        key = (int(shell), weight, None if emit_window is None else tuple(emit_window))
        if key in self._samplers:
            if self._samplers[key] is None:
                raise ValueError("thermal emissivity is empty")
            return self._samplers[key]
        s = int(shell)
        w = self.A_all * self._pops(s)[1] * self.nu0_all
        if weight == "energy_beta":
            beta = _beta_of(self._tau_of(s))
            free = np.ones(w.size, bool); free[self.op_idx[self.op_p[s] > 0]] = False
            w = w * np.where(free, 1.0, beta)
        lo, hi = emit_window
        w = np.where((self.nu0_all < lo) | (self.nu0_all > hi), 0.0, w)
        self._samplers[key] = _cut_sampler(w, self.emis_cut)
        if self._samplers[key] is None:
            raise ValueError("thermal emissivity is empty in the requested window")
        return self._samplers[key]

    def expansion_bins(self, dnu_over_nu=4.17e-5, nu_lo=None, nu_hi=None, weight="poisson"):
        """Log bins (shared edges) with per-shell E[s, b]; sets the per-shell
        within-bin cumulatives for `sample_line_in_bin`."""
        if weight not in ("poisson", "exact", "dual"):
            raise ValueError(f"weight must be 'poisson', 'exact' or 'dual', got {weight!r}")
        lo = self.op_nu.min() * (1 - 10 * dnu_over_nu) if nu_lo is None else nu_lo
        hi = self.op_nu.max() * (1 + 10 * dnu_over_nu) if nu_hi is None else nu_hi
        n = int(np.ceil(np.log(hi / lo) / np.log1p(dnu_over_nu))) + 1
        edges = lo * (1 + dnu_over_nu) ** np.arange(n + 1)
        b = np.clip(np.searchsorted(edges, self.op_nu, side="right") - 1, 0, n - 1)
        self._bin_of_line = b
        self._bin_start = np.searchsorted(b, np.arange(n + 1), side="left")
        E = np.zeros((self.n_shell, n))
        self._cum_in_bin = np.empty((self.n_shell, self.n_opacity))
        for s in range(self.n_shell):
            wgt = self.op_p[s] if weight == "poisson" else self.op_tau[s]
            np.add.at(E[s], b, wgt)
            sel = self.op_p[s] if weight in ("poisson", "dual") else self.op_tau[s]
            self._cum_in_bin[s] = np.cumsum(sel) - np.repeat(
                np.concatenate([[0.0], np.cumsum(sel)])[self._bin_start[:-1]], np.diff(self._bin_start))
        return edges, E

    def sample_line_in_bin(self, shells, b, u):
        """Absorbing line within bin b for packets in `shells` (the
        `ForestAtom` loop, with the shell's cumulative)."""
        b = np.array(b, copy=True); shells = np.asarray(shells, int)
        out = np.empty(b.size, int)
        for i in range(b.size):
            while self._bin_start[b[i]] == self._bin_start[b[i] + 1] and b[i] > 0:
                b[i] -= 1
            start, stop = self._bin_start[b[i]], self._bin_start[b[i] + 1]
            c = self._cum_in_bin[shells[i], start:stop]
            out[i] = min(start + np.searchsorted(c, u[i] * c[-1], side="right"), stop - 1)
        return out

    def shell_view(self, s):
        """The shell's opacity lines as a `band_saturation`-compatible object."""
        from types import SimpleNamespace
        live = self.op_p[s] > 0
        return SimpleNamespace(op_nu=self.op_nu[live], op_tau=self.op_tau[s][live], n_opacity=int(live.sum()))
