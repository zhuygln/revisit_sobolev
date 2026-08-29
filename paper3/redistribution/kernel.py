"""The reduced redistribution operator (Paper III, plan section 4).

Event-level: one row draw per absorption event, where an event is the
branch leg's full re-absorption chain collapsed to (nu_absorbed_cm ->
nu_exit_rest) -- exactly what `run_mc(collect_events=True)` records.
Photon-number rows drive the sampling (each event re-emits one photon); the
energy matrix R^E with the net comoving deposit fraction q_dep (negative
for net blueward fluorescence) is the conservation object, and
sum_j R^E_ij + q_dep_i = 1 holds exactly by construction per event.
Core loss is transport's job, not the kernel's.
"""
import json
import numpy as np


class RedistributionKernel:
    """Group-to-group redistribution: edges (n_g+1, ascending in frequency),
    R (energy, rows = absorbed group), N_cum (photon-count row cumulatives,
    the sampling object), q_dep, and per-output-group within-group photon
    sub-histograms (`sub_cum`, n_g x n_sub) marginalized over input groups.
    Empty rows sample to nan (transport scatters those coherently)."""

    def __init__(self, edges, R, N_cum, q_dep, sub_cum, counts, metadata=None,
                 disc_vals=None, disc_cum=None, disc_off=None):
        self.edges = np.asarray(edges, float)
        self.R = np.asarray(R, float)
        self.N_cum = np.asarray(N_cum, float)
        self.q_dep = np.asarray(q_dep, float)
        self.sub_cum = np.asarray(sub_cum, float)
        self.counts = np.asarray(counts, float)
        self.metadata = dict(metadata or {})
        self.n_groups = self.edges.size - 1
        self.empty_rows = self.counts <= 0
        # discrete within-group emission tables: exits are exact line rest
        # frequencies, and sampling them exactly lets the transport's
        # at-resonance convention skip the just-emitted line -- a continuous
        # PDF lands half its draws just ABOVE that line and double-counts the
        # self-absorption the kernel's training already resolved (the error
        # GROWS with n_groups; measured in compression_laII.json v1)
        self.disc_vals = None if disc_vals is None else np.asarray(disc_vals, float)
        self.disc_cum = None if disc_cum is None else np.asarray(disc_cum, float)
        self.disc_off = None if disc_off is None else np.asarray(disc_off, int)

    # ---- construction --------------------------------------------------
    @classmethod
    def from_branching_mc(cls, nu_in, nu_out, w, n_groups, nu_lo=None, nu_hi=None,
                          n_sub=16, metadata=None):
        """Build from event arrays (absorbed comoving frequency, exit rest
        frequency, packet weight). Groups are log-spaced over [nu_lo, nu_hi]
        (default: the events' own extent with a 0.1% margin); events outside
        are clipped into the edge groups."""
        nu_in = np.asarray(nu_in, float); nu_out = np.asarray(nu_out, float)
        w = np.asarray(w, float)
        lo = (min(nu_in.min(), nu_out.min()) * 0.999) if nu_lo is None else nu_lo
        hi = (max(nu_in.max(), nu_out.max()) * 1.001) if nu_hi is None else nu_hi
        edges = np.geomspace(lo, hi, n_groups + 1)
        gi = np.clip(np.searchsorted(edges, nu_in, side="right") - 1, 0, n_groups - 1)
        gj = np.clip(np.searchsorted(edges, nu_out, side="right") - 1, 0, n_groups - 1)
        E_in = np.bincount(gi, weights=w * nu_in, minlength=n_groups)
        flow = np.zeros((n_groups, n_groups))
        np.add.at(flow, (gi, gj), w * nu_out)
        Ncnt = np.zeros((n_groups, n_groups))
        np.add.at(Ncnt, (gi, gj), w)
        counts = Ncnt.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            R = np.where(E_in[:, None] > 0, flow / np.where(E_in[:, None] > 0, E_in[:, None], 1.0), 0.0)
            q_dep = np.where(E_in > 0, 1.0 - flow.sum(axis=1) / np.where(E_in > 0, E_in, 1.0), 0.0)
            N_cum = np.where(counts[:, None] > 0,
                             np.cumsum(Ncnt, axis=1) / np.where(counts[:, None] > 0, counts[:, None], 1.0), 0.0)
        # within-group emission PDF per OUTPUT group (photon-weighted,
        # marginal over input groups), on log-uniform sub-bins
        sub_cum = np.zeros((n_groups, n_sub))
        for j in range(n_groups):
            m = gj == j
            if not m.any():
                continue
            se = np.geomspace(edges[j], edges[j + 1], n_sub + 1)
            h, _ = np.histogram(nu_out[m], se, weights=w[m])
            tot = h.sum()
            sub_cum[j] = np.cumsum(h / tot) if tot > 0 else np.linspace(1 / n_sub, 1, n_sub)
        # discrete tables: the distinct exit frequencies per output group with
        # cumulative photon weights
        vals_u, inv = np.unique(nu_out, return_inverse=True)
        w_u = np.bincount(inv, weights=w)
        g_u = np.clip(np.searchsorted(edges, vals_u, side="right") - 1, 0, n_groups - 1)
        order = np.argsort(g_u, kind="stable")
        vals_s, w_s, g_s = vals_u[order], w_u[order], g_u[order]
        disc_off = np.searchsorted(g_s, np.arange(n_groups + 1))
        disc_cum = np.empty_like(w_s)
        for j in range(n_groups):
            a, b = disc_off[j], disc_off[j + 1]
            if b > a:
                c = np.cumsum(w_s[a:b]); disc_cum[a:b] = c / c[-1]
        md = dict(metadata or {}); md.setdefault("n_events", int(nu_in.size)); md["n_sub"] = n_sub
        return cls(edges, R, N_cum, q_dep, sub_cum, counts, md,
                   disc_vals=vals_s, disc_cum=disc_cum, disc_off=disc_off)

    # ---- use -----------------------------------------------------------
    def group_index(self, nu):
        return np.clip(np.searchsorted(self.edges, nu, side="right") - 1, 0, self.n_groups - 1)

    def sample_nu_out(self, nu_abs, rng, within="discrete"):
        """One output frequency per absorbed frequency; nan where the row was
        never populated in the reference (caller scatters those coherently)."""
        gi = self.group_index(np.asarray(nu_abs, float))
        out = np.full(gi.size, np.nan)
        u = rng.uniform(size=gi.size)
        for i in np.unique(gi):
            m = gi == i
            if self.empty_rows[i]:
                continue
            gj = np.searchsorted(self.N_cum[i], u[m])
            gj = np.minimum(gj, self.n_groups - 1)
            v = rng.uniform(size=gj.size)
            if within == "discrete" and self.disc_vals is not None:
                nu_j = np.empty(gj.size)
                for j in np.unique(gj):
                    mm = gj == j
                    a, b = self.disc_off[j], self.disc_off[j + 1]
                    if b <= a:      # no exits recorded in this output group
                        nu_j[mm] = np.nan
                        continue
                    k = np.minimum(np.searchsorted(self.disc_cum[a:b], v[mm]), b - a - 1)
                    nu_j[mm] = self.disc_vals[a + k]
                out[m] = nu_j
            elif within == "pdf":
                nu_j = np.empty(gj.size)
                for j in np.unique(gj):
                    mm = gj == j
                    se = np.geomspace(self.edges[j], self.edges[j + 1], self.sub_cum.shape[1] + 1)
                    k = np.searchsorted(self.sub_cum[j], v[mm])
                    k = np.minimum(k, self.sub_cum.shape[1] - 1)
                    lo = np.where(k > 0, self.sub_cum[j][k - 1], 0.0)
                    f = (v[mm] - lo) / np.maximum(self.sub_cum[j][k] - lo, 1e-300)
                    nu_j[mm] = se[k] * (se[k + 1] / se[k]) ** f
                out[m] = nu_j
            else:   # uniform in log within the group
                out[m] = self.edges[gj] * (self.edges[gj + 1] / self.edges[gj]) ** v
        return out

    # ---- validation / io ----------------------------------------------
    def validate_energy(self):
        """max |sum_j R_ij + q_dep_i - 1| over populated rows (exact up to
        float roundoff by construction)."""
        rows = ~self.empty_rows
        return float(np.abs(self.R[rows].sum(axis=1) + self.q_dep[rows] - 1.0).max()) if rows.any() else 0.0

    def block_flow(self, lam_in, lam_out, c_cgs=2.99792458e10):
        """Energy-flow fraction from wavelength band lam_in (A) into lam_out,
        relative to the energy absorbed in lam_in (uses R^E)."""
        mid = np.sqrt(self.edges[1:] * self.edges[:-1]); lam = c_cgs / mid * 1e8
        mi = (lam >= lam_in[0]) & (lam < lam_in[1]); mj = (lam >= lam_out[0]) & (lam < lam_out[1])
        E_in = self.counts[mi]  # photon-count proxy weights for the row average
        if not mi.any() or E_in.sum() <= 0:
            return np.nan
        return float((E_in[:, None] * self.R[mi][:, mj]).sum() / E_in.sum())

    def save(self, path):
        extra = {}
        if self.disc_vals is not None:
            extra = dict(disc_vals=self.disc_vals, disc_cum=self.disc_cum, disc_off=self.disc_off)
        np.savez(path, edges=self.edges, R=self.R, N_cum=self.N_cum, q_dep=self.q_dep,
                 sub_cum=self.sub_cum, counts=self.counts, metadata=json.dumps(self.metadata), **extra)

    @classmethod
    def load(cls, path):
        d = np.load(path, allow_pickle=False)
        disc = {k: d[k] for k in ("disc_vals", "disc_cum", "disc_off") if k in d.files}
        return cls(d["edges"], d["R"], d["N_cum"], d["q_dep"], d["sub_cum"], d["counts"],
                   json.loads(str(d["metadata"])),
                   disc_vals=disc.get("disc_vals"), disc_cum=disc.get("disc_cum"),
                   disc_off=disc.get("disc_off"))
