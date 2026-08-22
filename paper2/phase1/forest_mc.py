"""Paper II Phase 1: a vectorized Sobolev / expansion-opacity Monte Carlo with
line fluorescence, for a whole ion.

Phase 0 built and calibrated `paper2/phase0/three_level_atom/branching_mc.py`
-- per-packet Python loops, a three-level atom, ~250 lines -- and showed it
reproduces known branching ratios, interaction probabilities and Paper I's
analytic Sobolev leg. It has measured nothing, because O(lines x packets)
Python cannot carry a real atom: La II at 3000 K has ~950 lines with
tau_S > 1e-3 from 1148 to 17,600 A, and a pumped upper level may decay through
any of its 17,743 transitions. This module is the successor: same physics,
same conventions, packets advanced in lockstep with numpy, and three
interaction treatments in ONE code so that every comparison is same-code.

THE THREE LEGS (Phase 1's measurement):

  sobolev+thermal     Sobolev point interactions; on interaction, re-emit
                      isotropically at a line drawn from the LTE line
                      emissivity, photon-number weighted: P(k) ~ A_k n_u,k.
                      This is what SEDONA's resolved mode does with radiative
                      equilibrium on, at fixed T -- the leg Paper I's RE check
                      stood on -- except that SEDONA can only re-emit inside
                      its transport window. `emit_window` reproduces that.
  expansion+thermal   The expansion-opacity closure, alpha_exp per comoving
                      bin, as a continuous absorber; on interaction, the same
                      thermal sampler (the closure has no line identity).
  sobolev+branch      Sobolev point interactions; on interaction, sample the
                      upper level's downward channels by A and re-emit at that
                      line. The physics Paper II is about.

Plus two pure-absorption controls (sobolev, expansion) that must reproduce
Paper I's analytic legs -- the calibration against a calculation written
independently.

PHYSICS AND CONVENTIONS (unchanged from Phase 0, restated because they are
where cross-code comparisons go wrong):

* Homologous flow; first-order Doppler, nu_cm = nu_lab (1 - z/ct), strictly
  decreasing along any ray, so the next resonance a packet can meet is the
  opacity line with the LARGEST rest frequency below its current comoving
  frequency -- one searchsorted, not a loop. F11's worldline correction is
  O(beta^2) and irrelevant at 0.003-0.02 c.
* Sobolev point interaction with probability 1 - e^{-tau_S}; a re-emitted
  packet sits at s = 0 of its own resonance and moves away, so the resonance
  is consumed by the monotonicity, not by special-casing.
* ESCAPE PROBABILITY. A photon re-emitted inside a resonance zone has still
  to traverse the rest of that line's profile, and escapes the emitting line
  only with Sobolev's beta = (1 - e^{-tau})/tau; otherwise it is re-absorbed
  by the SAME line, at the same place, and re-emitted again with a fresh
  thermal or branching draw. A resolved calculation has this automatically;
  a point-interaction code must impose it. For pure resonant scattering the
  omission is invisible (same place, same frequency, isotropic either way),
  which is why Phase 0 could not catch it; for thermal and branching
  redistribution it decides which line the photon finally leaves through --
  trapped in strong lines until a draw lands on a weak one -- and it is what
  the SEDONA comparison in run_forest.py exposed.
  It applies to the SOBOLEV legs. The expansion closure is a continuous
  absorber whose E_bin already contains the emitting line's (1 - e^-tau): the
  next interaction in the same bin IS the re-absorption, with a fresh outcome
  draw -- SEDONA's expansion-mode semantics, which carry no beta. Imposing
  beta there as well counts the self-absorption twice; `beta_on_expansion`
  (default False) exists so the double count can be measured, not used. The
  expansion legs re-emit from the closure's OWN Kirchhoff emissivity,
  kappa_exp B_nu per bin -- which saturates at (1 - e^-tau) per strong line,
  unlike the Sobolev line emissivity A n_u -- placed uniformly within the bin
  (`exp_emit="bin"`); with that, and no beta, the leg reproduces SEDONA's
  expansion-mode radiative-equilibrium spectrum sub-band by sub-band (E0).
* ENERGY BOOKKEEPING. Packets are photons; each carries h nu, which changes at
  every re-emission. `accounting` reports E_inj = E_esc + E_core + E_abs +
  E_dep_lab (an identity, asserted to roundoff), with E_dep split into the
  comoving exchange with the gas, sum h(nu_cm,abs - nu_rest,emit), and the
  O(v/c) Doppler work term W = E_dep_lab - E_dep_cm. The branch leg's
  comoving exchange equals the level-energy difference hc(E_l,exit - E_l,pump)
  per chain (one emission per absorption event; the final lower level's
  excitation returns to the pool, a modelling choice). Per-packet weights `w`
  make photon- and energy-weighted spectra available from one run.
* Expansion opacity as a continuous absorber. Along a ray the comoving
  frequency falls linearly with path length, so the optical depth accumulated
  between two comoving frequencies is the sum over bins of (fraction of bin
  swept) x E_bin, E_bin = sum_{k in bin}(1 - e^{-tau_k}) -- the F13
  cancellation, which is exactly what a code integrating alpha_exp computes.
  Distance to interaction is therefore found by inverting a piecewise-linear
  cumulative, not by stepping. The bin width enters only through which lines
  share a bin; 4.17e-5 (SEDONA's production grid) is the default.
* Populations frozen (LTE at T, Boltzmann over the ion), opaque core,
  photon number conserved as launched = escaped + core-absorbed (+ absorbed,
  in the control modes). Packets are photons, not energy: a fluorescent
  packet keeps its count and changes its frequency.
* Flat launch continuum over [nu_min, nu_max], so escaped/launched per unit
  frequency IS the transmitted fraction, with nothing to divide out.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from sobolev.atomic_data import load_gsi
from sobolev.constants import C, SIGMA_CLASSICAL, H, K_B
from sobolev.populations import boltzmann_fractions_from_levels, statistical_weight

MODES = ("sobolev_absorb", "expansion_absorb", "sobolev_thermal",
         "expansion_thermal", "sobolev_branch",
         # E4: two-level-atom thermalisation parameter eps -- with probability
         # eps a thermal re-emission, otherwise coherent (resonant) scattering;
         # every event, including re-absorptions, draws again (SEDONA's
         # do_scatter semantics)
         "sobolev_tla", "expansion_tla",
         # E8: the closure with line identity restored at the interaction --
         # sample the absorbing line within the bin, exit by the A*beta kernel
         "expansion_branch")


# --------------------------------------------------------------------------
# Atom
# --------------------------------------------------------------------------


class ForestAtom:
    """A whole ion: opacity lines, branching tables, and an LTE emissivity.

    Parameters
    ----------
    nu0, f_osc, n_lower, A : per-line arrays (all lines of the ion).
    lower, upper : per-line level indices.
    n_upper : per-line upper-level number density (for the thermal sampler).
    t_exp : s.
    tau_min : lines with tau_S below this carry no opacity (they still
        appear in the branching tables and the emissivity). 1e-3 keeps ~950
        of La II's 17,743 lines at 3000 K.
    """

    def __init__(self, nu0, f_osc, n_lower, n_upper, A, lower, upper, t_exp,
                 tau_min=1e-3, stim=True, temperature=None, opacity_window=None):
        nu0 = np.asarray(nu0, float); f_osc = np.asarray(f_osc, float)
        n_lower = np.asarray(n_lower, float); n_upper = np.asarray(n_upper, float)
        A = np.asarray(A, float); lower = np.asarray(lower, int); upper = np.asarray(upper, int)
        self.n_lines_total = nu0.size
        tau = SIGMA_CLASSICAL * f_osc * n_lower * (C / nu0) * t_exp
        if stim and temperature is not None:
            # SEDONA's (1 - n_u g_l / n_l g_u); for Boltzmann populations
            # exactly 1 - exp(-h nu / kT). Paper I's convention, stated.
            tau = tau * (1.0 - np.exp(-H * nu0 / (K_B * temperature)))
        self.tau_all = tau

        # opacity lines, sorted by frequency ascending. `opacity_window`
        # (nu_lo, nu_hi) keeps opacity only from lines inside it -- Paper I's
        # window atom, whose branching tables and emissivity still span the
        # whole ion -- so the Paper I controls can be reproduced exactly.
        keep = tau > tau_min
        if opacity_window is not None:
            keep &= (nu0 >= opacity_window[0]) & (nu0 <= opacity_window[1])
        op = np.flatnonzero(keep)
        order = op[np.argsort(nu0[op])]
        self.op_idx = order                       # line index of each opacity line
        self.op_nu = nu0[order]
        self.op_tau = tau[order]
        self.op_upper = upper[order]
        self.op_p = -np.expm1(-self.op_tau)       # interaction probability

        # branching tables: for every upper level that any opacity line feeds,
        # the downward lines (ALL of them) and cumulative A
        self.nu0_all = nu0; self.A_all = A; self.upper_all = upper; self.lower_all = lower
        self.level_energy_cm = None   # set by from_gsi (or by hand) for the level-energy identity
        self.branch_lines = {}
        self.branch_cum = {}
        for u in np.unique(self.op_upper):
            idx = np.flatnonzero(upper == u)
            a = A[idx]
            if a.sum() <= 0:
                raise ValueError(f"upper level {u} has no downward rate")
            self.branch_lines[u] = idx
            self.branch_cum[u] = np.cumsum(a / a.sum())

        # thermal (LTE) line emissivity, photon-number weighted: A n_u
        w = A * n_upper
        self.emis_w = w
        self.n_opacity = order.size
        # Sobolev escape probability of every line, for photons re-emitted in it
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.where(tau > 1e-12, -np.expm1(-tau) / np.where(tau > 1e-12, tau, 1.0), 1.0)
        self.beta_all = beta
        # lines that carry no opacity here (below tau_min or outside the
        # opacity window) are treated as freely escaping
        free = np.ones(nu0.size, bool); free[order] = False
        self.beta_all = np.where(free, 1.0, beta)
        # Exit kernel (E8): the chain of A-branching draws with re-absorption
        # by the emitting line exits through j with probability
        # A_uj beta_uj / sum_m A_um beta_um -- the closed form of the loop.
        self.exit_cum = {}
        for u, idx in self.branch_lines.items():
            w_exit = A[idx] * self.beta_all[idx]
            self.exit_cum[u] = np.cumsum(w_exit / w_exit.sum())

    @classmethod
    def from_gsi(cls, levels_path, transitions_path, temperature, n_ion, t_exp,
                 tau_min=1e-3, stim=True, opacity_window=None):
        lev = load_gsi(levels_path); tr = load_gsi(transitions_path)
        frac = boltzmann_fractions_from_levels(lev, temperature)
        low = tr["Lower"].to_numpy(); up = tr["Upper"].to_numpy()
        g_l = statistical_weight(tr["J_Lower"].to_numpy())
        f_lu = 10 ** tr["Log(gf)"].to_numpy() / g_l
        nu0 = C / (tr["WV_Transition"].to_numpy() * 1e-8)
        atom = cls(nu0, f_lu, frac[low] * n_ion, frac[up] * n_ion, tr["A"].to_numpy(),
                   low, up, t_exp, tau_min=tau_min, stim=stim, temperature=temperature,
                   opacity_window=opacity_window)
        atom.levels, atom.transitions, atom.temperature, atom.n_ion = lev, tr, temperature, n_ion
        atom.level_energy_cm = lev["Energy"].to_numpy(dtype=float)
        return atom

    # ---- samplers ------------------------------------------------------
    def thermal_sampler(self, emit_window=None):
        """Return a function u -> line index drawing from A n_u, optionally
        restricted to lines with rest frequency inside `emit_window`
        (nu_lo, nu_hi) -- SEDONA's window-confined re-emission."""
        w = self.emis_w.copy()
        if emit_window is not None:
            lo, hi = emit_window
            w[(self.nu0_all < lo) | (self.nu0_all > hi)] = 0.0
        if w.sum() <= 0:
            raise ValueError("thermal emissivity is empty in the requested window")
        cum = np.cumsum(w / w.sum())
        def sample(u):
            return np.searchsorted(cum, u)
        return sample

    def expansion_bins(self, dnu_over_nu=4.17e-5, nu_lo=None, nu_hi=None):
        """Log-spaced comoving bins and E_bin = sum (1 - e^-tau) per bin, over
        the opacity lines. Returns (edges ascending, E per bin, cumulative E
        from the TOP down)."""
        lo = self.op_nu.min() * (1 - 10 * dnu_over_nu) if nu_lo is None else nu_lo
        hi = self.op_nu.max() * (1 + 10 * dnu_over_nu) if nu_hi is None else nu_hi
        n = int(np.ceil(np.log(hi / lo) / np.log1p(dnu_over_nu))) + 1
        edges = lo * (1 + dnu_over_nu) ** np.arange(n + 1)
        E = np.zeros(n)
        b = np.clip(np.searchsorted(edges, self.op_nu, side="right") - 1, 0, n - 1)
        np.add.at(E, b, self.op_p)
        # CSR by bin for sampling the absorbing line within a bin (E8): lines
        # are already sorted by frequency, so each bin's lines are contiguous
        self._bin_of_line = b
        self._bin_start = np.searchsorted(b, np.arange(n + 1), side="left")
        self._cum_in_bin = np.cumsum(self.op_p) - np.repeat(np.concatenate([[0.0], np.cumsum(self.op_p)])[self._bin_start[:-1]], np.diff(self._bin_start))
        return edges, E

    def sample_line_in_bin(self, b, u):
        """Absorbing line index (into op_*) within bin b, with probability
        op_p[k]/E_b; u uniform. Vectorized over packets."""
        b = np.array(b, copy=True)
        out = np.empty(b.size, int)
        for i in range(b.size):   # small loops only over interacting packets per step
            # an interaction point can land exactly on a bin edge and be
            # assigned to the (empty) bin above; the optical depth was
            # accumulated in the populated bin just below -- step down to it
            while self._bin_start[b[i]] == self._bin_start[b[i] + 1] and b[i] > 0:
                b[i] -= 1
            start, stop = self._bin_start[b[i]], self._bin_start[b[i] + 1]
            c = self._cum_in_bin[start:stop]
            out[i] = min(start + np.searchsorted(c, u[i] * c[-1], side="right"), stop - 1)
        return out


# --------------------------------------------------------------------------
# Geometry (vectorized versions of Phase 0's helpers)
# --------------------------------------------------------------------------


def distance_to_boundary(r, mu, r_core, r_out):
    """(distance, hits_core) for arrays r, mu."""
    perp2 = r * r * (1.0 - mu * mu)
    s_out = -r * mu + np.sqrt(np.maximum(r_out * r_out - perp2, 0.0))
    hits_core = (mu < 0.0) & (perp2 < r_core * r_core)
    s_core = np.where(hits_core, -r * mu - np.sqrt(np.maximum(r_core * r_core - perp2, 0.0)), np.inf)
    core_first = hits_core & (s_core < s_out)
    return np.where(core_first, s_core, s_out), core_first


def advance(r, mu, s):
    r_new = np.sqrt(np.maximum(r * r + s * s + 2.0 * r * s * mu, 0.0))
    mu_new = np.where(r_new > 0.0, (r * mu + s) / np.where(r_new > 0.0, r_new, 1.0), 1.0)
    return r_new, mu_new


# --------------------------------------------------------------------------
# Transport
# --------------------------------------------------------------------------

_S_MIN = 1.0  # cm; a resonance must lie strictly ahead (see Phase 0)


def sample_launch(rng, nu_min, nu_max, n, t_core=None):
    """Launch frequencies: flat in nu, or photon-number weighted Planck
    B_nu(T)/(h nu) over [nu_min, nu_max] (rejection from a log grid CDF)."""
    if t_core is None:
        return rng.uniform(nu_min, nu_max, n)
    grid = np.geomspace(nu_min, nu_max, 20001)
    x = H * grid / (K_B * t_core)
    w = grid**2 / np.expm1(np.minimum(x, 700.0))       # photon-number emissivity ~ nu^2/(e^x-1)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    return np.interp(rng.uniform(0.0, 1.0, n), cdf, grid)


def run_mc(atom, r_core, r_out, t_exp, nu_min, nu_max, n_packets, mode,
           seed=0, emit_window=None, dnu_over_nu=4.17e-5, max_steps=100000,
           t_core=None, beta_on_expansion=False, exp_emit="bin",
           launch_weight="photon", eps=1.0):
    """exp_emit : how the expansion legs place a thermally re-emitted packet
    in frequency. "bin" (default) -- uniformly within the sampled line's
    comoving bin, which is what an emissivity kappa_exp B_nu per bin means and
    what SEDONA does; the packet then sweeps the rest of the bin at the bin's
    opacity, and an escape probability (1 - e^-E_b)/E_b emerges from the
    geometry. "line" -- at the line's exact rest frequency (Phase-1 behaviour,
    kept for the E0 record)."""
    """Propagate packets from the core through the shell, in lockstep.

    mode : one of MODES.
    emit_window : (nu_lo, nu_hi) to confine thermal re-emission to a window
        (SEDONA-like); None re-emits over the whole atom.
    Returns dict with nu_launch, nu_out, counts, and n_interactions (total
    interaction events, a diagnostic).
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    rng = np.random.default_rng(seed)
    ct = C * t_exp
    sobolev = mode.startswith("sobolev")
    outcome = mode.split("_")[1]
    needs_thermal = outcome in ("thermal", "tla")
    thermal = atom.thermal_sampler(emit_window) if (needs_thermal and (sobolev or exp_emit == "line")) else None
    if not sobolev:
        edges, E = atom.expansion_bins(dnu_over_nu, nu_lo=min(nu_min, atom.op_nu.min()) * 0.5,
                                       nu_hi=max(nu_max, atom.op_nu.max()) * 1.01)
        # cumulative E from the top (high frequency) downward; G(nu) = optical
        # depth accumulated sweeping from the top edge down to nu.
        # G is piecewise linear in nu inside each bin.
        G_edges = np.concatenate([[0.0], np.cumsum(E[::-1])])[::-1]   # G at each edge, top = 0
        width = np.diff(edges)
        if needs_thermal and exp_emit == "bin":
            # Kirchhoff for the closure's OWN opacity: emissivity per bin is
            # kappa_exp(b) B_nu(T) -- photon-number weight E_b * B_nu/(h nu) --
            # which saturates at (1 - e^-tau) per strong line, unlike the
            # Sobolev line emissivity A n_u (~ tau). This is what SEDONA's
            # expansion mode re-emits from; frequency uniform within the bin.
            nu_b = np.sqrt(edges[1:] * edges[:-1])
            T_em = getattr(atom, "temperature", None)
            if T_em is None:
                # toy atoms without a temperature: flat B_nu (weights E_b only)
                w_b = E * width
            else:
                xb = H * nu_b / (K_B * T_em)
                w_b = E * nu_b**2 / np.expm1(np.minimum(xb, 700.0)) * width
            if emit_window is not None:
                w_b[(nu_b < emit_window[0]) | (nu_b > emit_window[1])] = 0.0
            if w_b.sum() <= 0:
                raise ValueError("expansion thermal emissivity is empty in the window")
            cum_b = np.cumsum(w_b / w_b.sum())
            def thermal_bin(u):
                return np.searchsorted(cum_b, u)
        def G_of(nu):
            b = np.clip(np.searchsorted(edges, nu, side="right") - 1, 0, E.size - 1)
            frac = (edges[b + 1] - nu) / width[b]
            return G_edges[b + 1] + frac * E[b]
        def nu_of_G(g):
            # invert: find bin where G_edges[b+1] <= g < G_edges[b]
            b = np.clip(E.size - 1 - np.searchsorted(G_edges[::-1], g, side="right"), 0, E.size - 1)
            # G_edges descending in b; G_edges[::-1] ascending
            # within bin b: g = G_edges[b+1] + frac*E[b] -> frac
            frac = np.where(E[b] > 0, (g - G_edges[b + 1]) / np.where(E[b] > 0, E[b], 1.0), 0.0)
            return edges[b + 1] - frac * width[b]

    if launch_weight == "photon" or t_core is None:
        nu_launch = sample_launch(rng, nu_min, nu_max, n_packets, t_core)
        w = np.ones(n_packets)
    elif launch_weight == "energy":
        # launch in proportion to B_nu (energy); each packet then stands for
        # w ~ 1/nu photons, normalized to n_packets photons in total
        grid = np.geomspace(nu_min, nu_max, 20001)
        x = H * grid / (K_B * t_core)
        wgt = grid**3 / np.expm1(np.minimum(x, 700.0))
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * (wgt[1:] + wgt[:-1]) * np.diff(grid))]); cdf /= cdf[-1]
        nu_launch = np.interp(rng.uniform(0.0, 1.0, n_packets), cdf, grid)
        w = 1.0 / nu_launch; w *= n_packets / w.sum()
    else:
        raise ValueError("launch_weight must be 'photon' or 'energy'")
    r = np.full(n_packets, float(r_core))
    mu = np.sqrt(rng.uniform(0.0, 1.0, n_packets))
    nu = nu_launch.copy()
    alive = np.ones(n_packets, bool)
    nu_out = np.full(n_packets, np.nan)
    fate = np.zeros(n_packets, np.int8)   # 0 alive, 1 escaped, 2 core, 3 absorbed
    n_inter = np.zeros(n_packets, np.int32)
    # branch chosen at each packet's FIRST interaction (branch mode), tallied
    # per line index: subsequent interactions are conditioned on the first,
    # so pooling them would bias a branching-ratio estimate (Phase 0's rule).
    first_branch = np.zeros(atom.n_lines_total, np.int64)
    # line through which each packet's FIRST interaction chain finally
    # escaped (after any re-absorptions in the emitting lines)
    chain_exit = np.zeros(atom.n_lines_total, np.int64)
    # energy bookkeeping (per packet, in erg per photon; weights applied at the end)
    e_dep_lab = np.zeros(n_packets)   # sum h (nu_lab before - nu_lab after) over re-emissions
    e_dep_cm = np.zeros(n_packets)    # sum h (nu_cm at absorption - nu_rest emitted)
    e_lev = np.zeros(n_packets)       # branch mode: sum hc (E_l,exit - E_l,pump) per chain
    n_events = np.zeros(n_packets, np.int32)   # interaction events (chains)
    n_reabs = np.zeros(n_packets, np.int32)    # re-absorptions inside emitting lines
    nu_final = np.full(n_packets, np.nan)      # lab frequency at death (escape/core/absorb)
    first_line = np.full(n_packets, -1, np.int64)   # E7: first absorbing line (-1 expansion)
    last_line = np.full(n_packets, -1, np.int64)    # E7: last emitting line
    # expansion mode: optical depth still to travel before the next interaction
    tau_r = rng.exponential(1.0, n_packets) if not sobolev else None

    for step in range(max_steps):
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        ri, mi, ni = r[idx], mu[idx], nu[idx]
        s_b, core_first = distance_to_boundary(ri, mi, r_core, r_out)
        z = ri * mi
        nu_cm = ni * (1.0 - z / ct)

        if sobolev:
            # next opacity line strictly below the current comoving frequency
            k = np.searchsorted(atom.op_nu, nu_cm, side="left") - 1
            has = k >= 0
            kk = np.where(has, k, 0)
            s_res = np.where(has, ct * (1.0 - atom.op_nu[kk] / ni) - z, np.inf)
            # a resonance at (numerically) zero distance is the one just used
            s_res = np.where(s_res > _S_MIN, s_res, np.inf)
            interact_candidate = s_res < s_b
        else:
            g_now = G_of(nu_cm)
            g_target = g_now + tau_r[idx]
            nu_target = nu_of_G(g_target)
            # can't go below the bottom edge
            reachable = g_target <= G_edges[0]
            s_res = np.where(reachable, ct * (1.0 - nu_target / ni) - z, np.inf)
            s_res = np.where(s_res > _S_MIN, s_res, np.inf)
            interact_candidate = s_res < s_b

        # --- packets that reach a boundary first
        esc = ~interact_candidate
        e_idx = idx[esc]
        nu_out[e_idx[~core_first[esc]]] = ni[esc][~core_first[esc]]
        fate[e_idx[~core_first[esc]]] = 1
        fate[e_idx[core_first[esc]]] = 2
        nu_final[e_idx] = ni[esc]
        alive[e_idx] = False

        # --- packets that reach a resonance / interaction point
        c = interact_candidate
        if not c.any():
            continue
        ci = idx[c]
        rn, mn = advance(ri[c], mi[c], s_res[c])
        r[ci], mu[ci] = rn, mn
        if sobolev:
            kc = kk[c]
            hit = rng.uniform(size=ci.size) < atom.op_p[kc]
        else:
            hit = np.ones(ci.size, bool)
            # expansion mode: those that don't interact don't exist; every
            # arrival at nu_target is an interaction. Redraw tau_r for survivors later.
        # non-interacting: continue from the resonance (Sobolev) -- nothing to do
        hi = ci[hit]
        if hi.size == 0:
            continue
        first = n_inter[hi] == 0
        n_inter[hi] += 1
        n_events[hi] += 1
        # comoving frequency at absorption: the line's rest frequency in the
        # Sobolev legs, the interaction point's comoving frequency otherwise
        nu_abs_cm = (atom.op_nu[kc] if sobolev else nu_target[c])[hit]
        nu_lab_before = nu[hi].copy()
        if sobolev:
            first_line[hi[first]] = atom.op_idx[kc[hit][first]]
        if outcome == "absorb":
            fate[hi] = 3; alive[hi] = False
            nu_final[hi] = nu_lab_before
            continue
        # re-emit isotropically, then apply the emitting line's escape
        # probability: re-absorbed packets draw again (same place, new line,
        # new direction) until one escapes. `cur_up` tracks the upper level
        # the packet currently sits in (branch mode); a re-absorption in the
        # line just emitted returns it to that line's upper level.
        todo = np.arange(hi.size)
        bin_emit = not sobolev and exp_emit == "bin"   # expansion legs: thermal draws are bins
        is_bin = np.zeros(hi.size, bool)               # which entries of new_line are bin ids
        coherent = np.zeros(hi.size, bool)             # expansion tla: coherent scatter
        if outcome == "branch" and sobolev:
            cur_up = atom.op_upper[kc[hit]].copy()
        if outcome == "tla" and sobolev:
            cur_line = atom.op_idx[kc[hit]].copy()     # the line just interacted with
        if outcome == "branch" and not sobolev:
            # E8: restore line identity at the interaction point -- the
            # absorbing line within the bin, with probability op_p[k]/E_b
            b_hit = np.clip(np.searchsorted(edges, nu_abs_cm, side="right") - 1, 0, E.size - 1)
            k_abs = atom.sample_line_in_bin(b_hit, rng.uniform(size=hi.size))
            first_line[hi[first]] = atom.op_idx[k_abs[first]]
            cur_up = atom.op_upper[k_abs].copy()
        new_line = np.empty(hi.size, int)
        n_chain = 0
        while todo.size:
            n_chain += 1
            if outcome == "branch":
                if sobolev:
                    u = rng.uniform(size=todo.size)
                    for uval in np.unique(cur_up[todo]):
                        m = cur_up[todo] == uval; sel = todo[m]
                        lines_u = atom.branch_lines[uval]; cum_u = atom.branch_cum[uval]
                        new_line[sel] = lines_u[np.searchsorted(cum_u, u[m])]
                else:
                    # exit by the A*beta kernel: the chain in closed form
                    u = rng.uniform(size=todo.size)
                    for uval in np.unique(cur_up[todo]):
                        m = cur_up[todo] == uval; sel = todo[m]
                        lines_u = atom.branch_lines[uval]; cum_u = atom.exit_cum[uval]
                        new_line[sel] = lines_u[np.searchsorted(cum_u, u[m])]
            elif outcome == "tla":
                th = rng.uniform(size=todo.size) < eps
                if sobolev:
                    new_line[todo[~th]] = cur_line[todo[~th]]          # resonant: same line
                    if th.any():
                        new_line[todo[th]] = thermal(rng.uniform(size=th.sum()))
                else:
                    coherent[todo[~th]] = True
                    if th.any():
                        if bin_emit:
                            new_line[todo[th]] = thermal_bin(rng.uniform(size=th.sum())); is_bin[todo[th]] = True
                        else:
                            new_line[todo[th]] = thermal(rng.uniform(size=th.sum()))
            elif thermal is not None:  # thermal, line-based (Sobolev legs, or exp_emit="line")
                new_line[todo] = thermal(rng.uniform(size=todo.size))
            else:  # expansion + thermal, bin-based
                new_line[todo] = thermal_bin(rng.uniform(size=todo.size)); is_bin[todo] = True
            if n_chain == 1:
                ok = first & ~is_bin & ~coherent
                np.add.at(first_branch, new_line[ok], 1)
            # escape the emitting line? (Sobolev legs only; the continuous
            # absorber re-absorbs through its own bins; the E8 kernel already
            # contains its beta)
            if sobolev or (beta_on_expansion and exp_emit == "line" and outcome == "thermal"):
                esc = rng.uniform(size=todo.size) < atom.beta_all[new_line[todo]]
            else:
                esc = np.ones(todo.size, bool)
            stay = todo[~esc]
            if stay.size:
                n_inter[hi[stay]] += 1
                n_reabs[hi[stay]] += 1
                if outcome == "branch" and sobolev:
                    cur_up[stay] = atom.upper_all[new_line[stay]]
                if outcome == "tla" and sobolev:
                    cur_line[stay] = new_line[stay]
            todo = stay
            if n_chain > 10000:
                raise RuntimeError("re-absorption chain did not terminate")
        nu_rest = np.empty(hi.size)
        nu_rest[~is_bin & ~coherent] = atom.nu0_all[new_line[~is_bin & ~coherent]]
        if is_bin.any():
            nu_rest[is_bin] = edges[new_line[is_bin]] + rng.uniform(size=is_bin.sum()) * width[new_line[is_bin]]
        if coherent.any():
            nu_rest[coherent] = nu_abs_cm[coherent]
        ok = first & ~is_bin & ~coherent
        np.add.at(chain_exit, new_line[ok], 1)
        mu_new = rng.uniform(-1.0, 1.0, hi.size)
        zn = r[hi] * mu_new
        mu[hi] = mu_new
        nu[hi] = nu_rest / (1.0 - zn / ct)
        e_dep_lab[hi] += H * (nu_lab_before - nu[hi])
        e_dep_cm[hi] += H * (nu_abs_cm - nu_rest)
        lined = ~is_bin & ~coherent
        last_line[hi[lined]] = new_line[lined]
        if outcome == "branch" and atom.level_energy_cm is not None:
            E_cm = atom.level_energy_cm
            pump = atom.op_idx[kc[hit]] if sobolev else atom.op_idx[k_abs]
            e_lev[hi] += H * C * (E_cm[atom.lower_all[new_line]] - E_cm[atom.lower_all[pump]])
        if not sobolev:
            tau_r[hi] = rng.exponential(1.0, hi.size)
    else:
        raise RuntimeError(f"packets still alive after {max_steps} steps")

    E_inj = float(np.sum(w * H * nu_launch))
    E_esc = float(np.sum(w[fate == 1] * H * nu_final[fate == 1]))
    E_core = float(np.sum(w[fate == 2] * H * nu_final[fate == 2]))
    E_abs = float(np.sum(w[fate == 3] * H * nu_final[fate == 3]))
    E_dep_lab = float(np.sum(w * e_dep_lab)); E_dep_cm = float(np.sum(w * e_dep_cm))
    accounting = dict(E_inj=E_inj, E_esc=E_esc, E_core=E_core, E_abs=E_abs,
                      E_dep_lab=E_dep_lab, E_dep_cm=E_dep_cm, W=E_dep_lab - E_dep_cm,
                      E_interacting=float(np.sum(w[n_events > 0] * H * nu_launch[n_events > 0])),
                      identity_residual=(E_esc + E_core + E_abs + E_dep_lab - E_inj) / E_inj,
                      N_inj=float(w.sum()), N_esc=float(w[fate == 1].sum()),
                      N_core=float(w[fate == 2].sum()), N_abs=float(w[fate == 3].sum()))
    return dict(nu_launch=nu_launch, nu_out=nu_out[fate == 1], nu_out_all=nu_out, fate=fate, w=w,
                n_packets=n_packets, n_escaped=int((fate == 1).sum()),
                n_core=int((fate == 2).sum()), n_absorbed=int((fate == 3).sum()),
                n_interactions=int(n_inter.sum()),
                n_interacted=int((n_inter > 0).sum()), first_branch=first_branch,
                chain_exit=chain_exit, steps=step + 1, accounting=accounting,
                e_dep_lab=e_dep_lab, e_dep_cm=e_dep_cm, e_lev=e_lev, n_events=n_events,
                n_reabs=n_reabs, nu_final=nu_final, first_line=first_line, last_line=last_line)


def _weights(res, weight):
    """Per-packet weights for launched and escaped packets: photon number, or
    energy (h nu) on top of the photon weight."""
    w = res.get("w", np.ones(res["nu_launch"].size))
    esc = res["fate"] == 1
    if weight == "photon":
        return w, w[esc]
    if weight == "energy":
        return w * res["nu_launch"], (w * res["nu_out_all"])[esc]
    raise ValueError("weight must be 'photon' or 'energy'")


def band_ratio(res, nu_lo, nu_hi, weight="photon"):
    """Escaped/launched in [nu_lo, nu_hi], photon- or energy-weighted; with a
    flat launch this is the transmitted (or emergent) fraction, with a Planck
    launch the emergent-to-incident ratio. Returns (value, poisson_error)."""
    w_in, w_out = _weights(res, weight)
    sel_in = (res["nu_launch"] >= nu_lo) & (res["nu_launch"] < nu_hi)
    sel_out = (res["nu_out"] >= nu_lo) & (res["nu_out"] < nu_hi)
    s_in, s_out = w_in[sel_in].sum(), w_out[sel_out].sum()
    if s_in <= 0:
        return np.nan, np.nan
    n_out = max(np.count_nonzero(sel_out), 1)
    return s_out / s_in, (s_out / s_in) / np.sqrt(n_out)


def spectrum(res, edges, weight="photon"):
    """Emergent / launched per frequency bin (photon- or energy-weighted)."""
    w_in, w_out = _weights(res, weight)
    out, _ = np.histogram(res["nu_out"], edges, weights=w_out)
    inn, _ = np.histogram(res["nu_launch"], edges, weights=w_in)
    cnt, _ = np.histogram(res["nu_out"], edges)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (np.where(inn > 0, out / inn, np.nan),
                np.where(inn > 0, (out / inn) / np.sqrt(np.maximum(cnt, 1)), np.nan))
