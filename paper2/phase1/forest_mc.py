"""Paper II Phase 1: a vectorized Sobolev / expansion-opacity Monte Carlo with
line fluorescence, for a whole ion.

Phase 0 built and calibrated `paper2/phase0/three_level_atom/branching_mc.py`
-- per-packet Python loops, a three-level atom, ~250 lines -- and showed it
reproduces known branching ratios, interaction probabilities and Paper I's
analytic Sobolev leg. It has measured nothing, because O(lines x packets)
Python cannot carry a real atom: La II at 3000 K has ~950 lines with
tau_S > 1e-3 from 1148 to 23,000 A, and a pumped upper level may decay through
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
         "expansion_thermal", "sobolev_branch")


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
        self.nu0_all = nu0; self.A_all = A; self.upper_all = upper
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
        return edges, E


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
           t_core=None):
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

    thermal = atom.thermal_sampler(emit_window) if outcome == "thermal" else None
    if not sobolev:
        edges, E = atom.expansion_bins(dnu_over_nu, nu_lo=min(nu_min, atom.op_nu.min()) * 0.5,
                                       nu_hi=max(nu_max, atom.op_nu.max()) * 1.01)
        # cumulative E from the top (high frequency) downward; G(nu) = optical
        # depth accumulated sweeping from the top edge down to nu.
        # G is piecewise linear in nu inside each bin.
        G_edges = np.concatenate([[0.0], np.cumsum(E[::-1])])[::-1]   # G at each edge, top = 0
        width = np.diff(edges)
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

    nu_launch = sample_launch(rng, nu_min, nu_max, n_packets, t_core)
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
        if outcome == "absorb":
            fate[hi] = 3; alive[hi] = False
            continue
        # re-emit isotropically, then apply the emitting line's escape
        # probability: re-absorbed packets draw again (same place, new line,
        # new direction) until one escapes. `cur_up` tracks the upper level
        # the packet currently sits in (branch mode); a re-absorption in the
        # line just emitted returns it to that line's upper level.
        todo = np.arange(hi.size)
        if outcome == "branch":
            cur_up = atom.op_upper[kc[hit]].copy()
        new_line = np.empty(hi.size, int)
        n_chain = 0
        while todo.size:
            n_chain += 1
            if outcome == "branch":
                u = rng.uniform(size=todo.size)
                for uval in np.unique(cur_up[todo]):
                    sel = todo[cur_up[todo] == uval]
                    lines_u = atom.branch_lines[uval]; cum_u = atom.branch_cum[uval]
                    new_line[sel] = lines_u[np.searchsorted(cum_u, u[cur_up[todo] == uval])]
            else:  # thermal
                new_line[todo] = thermal(rng.uniform(size=todo.size))
            if n_chain == 1:
                np.add.at(first_branch, new_line[first], 1)
            # escape the emitting line?
            esc = rng.uniform(size=todo.size) < atom.beta_all[new_line[todo]]
            stay = todo[~esc]
            if stay.size:
                n_inter[hi[stay]] += 1
                if outcome == "branch":
                    cur_up[stay] = atom.upper_all[new_line[stay]]
            todo = stay
            if n_chain > 10000:
                raise RuntimeError("re-absorption chain did not terminate")
        np.add.at(chain_exit, new_line[first], 1)
        nu_rest = atom.nu0_all[new_line]
        mu_new = rng.uniform(-1.0, 1.0, hi.size)
        zn = r[hi] * mu_new
        mu[hi] = mu_new
        nu[hi] = nu_rest / (1.0 - zn / ct)
        if not sobolev:
            tau_r[hi] = rng.exponential(1.0, hi.size)
    else:
        raise RuntimeError(f"packets still alive after {max_steps} steps")

    return dict(nu_launch=nu_launch, nu_out=nu_out[fate == 1], nu_out_all=nu_out, fate=fate,
                n_packets=n_packets, n_escaped=int((fate == 1).sum()),
                n_core=int((fate == 2).sum()), n_absorbed=int((fate == 3).sum()),
                n_interactions=int(n_inter.sum()),
                n_interacted=int((n_inter > 0).sum()), first_branch=first_branch,
                chain_exit=chain_exit, steps=step + 1)


def band_ratio(res, nu_lo, nu_hi, n_bins=1):
    """Escaped/launched per unit frequency in [nu_lo, nu_hi]; with a flat
    launch continuum this is the transmitted (or emergent) fraction. Returns
    (value, poisson_error)."""
    n_out = np.count_nonzero((res["nu_out"] >= nu_lo) & (res["nu_out"] < nu_hi))
    n_in = np.count_nonzero((res["nu_launch"] >= nu_lo) & (res["nu_launch"] < nu_hi))
    if n_in == 0:
        return np.nan, np.nan
    return n_out / n_in, np.sqrt(max(n_out, 1)) / n_in


def spectrum(res, edges):
    """Emergent / launched per frequency bin."""
    out, _ = np.histogram(res["nu_out"], edges)
    inn, _ = np.histogram(res["nu_launch"], edges)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(inn > 0, out / inn, np.nan), np.where(inn > 0, np.sqrt(np.maximum(out, 1)) / inn, np.nan)
