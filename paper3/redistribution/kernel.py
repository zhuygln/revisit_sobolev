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

    # ---- composition (P9) ----------------------------------------------
    @classmethod
    def mix(cls, kernels, w, metadata=None):
        """Opacity-weighted composition rule: R_mix[i] = sum_s w[i,s] R_s[i].

        `kernels` are per-species kernels on IDENTICAL edges (build them with
        the same nu_lo/nu_hi/n_groups); `w` is (n_groups, n_species) with rows
        summing to 1 -- the fraction of group-i absorption belonging to each
        species. The point of the rule is that w comes from the mixture's
        opacity alone, so a blend needs no blend-specific training.

        Convex mixing preserves the conservation object exactly: sum_j R_mix
        + q_mix = sum_s w[i,s] (sum_j R_s + q_s) = sum_s w[i,s] = 1.

        Rows a species never populated contribute nothing and its weight is
        redistributed over the species that did; a row no species populated
        stays empty (transport scatters those coherently).

        The within-group exit tables are marginal over the input group -- the
        same approximation the single-species kernel already makes -- so they
        are merged per output group, weighting each species by the flow the
        mixture sends it there.
        """
        ks = list(kernels)
        edges = ks[0].edges
        for k in ks[1:]:
            if k.edges.shape != edges.shape or not np.allclose(k.edges, edges):
                raise ValueError("mix() needs kernels on identical edges")
        w = np.asarray(w, float)
        ng, ns = ks[0].n_groups, len(ks)
        if w.shape != (ng, ns):
            raise ValueError(f"w must be (n_groups, n_species) = {(ng, ns)}, got {w.shape}")

        # per-species row-normalized photon rows
        rows = np.stack([np.diff(np.hstack([np.zeros((ng, 1)), k.N_cum]), axis=1) for k in ks])
        live = np.stack([~k.empty_rows for k in ks])                  # (ns, ng)
        wl = w.T * live                                                # zero out unpopulated
        norm = wl.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            wn = np.where(norm > 0, wl / np.where(norm > 0, norm, 1.0), 0.0)

        N_row = np.einsum("si,sij->ij", wn, rows)
        R = np.einsum("si,sij->ij", wn, np.stack([k.R for k in ks]))
        q_dep = np.einsum("si,si->i", wn, np.stack([k.q_dep for k in ks]))
        counts = np.einsum("si,si->i", wn, np.stack([k.counts for k in ks]))
        counts = np.where(norm > 0, counts, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            tot = N_row.sum(axis=1, keepdims=True)
            N_cum = np.where(tot > 0, np.cumsum(N_row, axis=1) / np.where(tot > 0, tot, 1.0), 0.0)

        # Flow the mixture sends into each output group, per species. The
        # per-group absorption budget must be SPECIES-INDEPENDENT (counts,
        # the mixture's own row occupancy), or each species is weighted by
        # its absolute run size as well as by w and the composition is
        # counted twice -- which shows up as a mixture worse than its own
        # dominant component.
        T = np.einsum("i,si,sij->sj", counts, wn, rows)

        sub_cum = np.zeros((ng, ks[0].sub_cum.shape[1]))
        for j in range(ng):
            if T[:, j].sum() <= 0:
                sub_cum[j] = np.linspace(1 / sub_cum.shape[1], 1, sub_cum.shape[1])
                continue
            h = sum(T[si, j] * np.diff(np.hstack([0.0, k.sub_cum[j]])) for si, k in enumerate(ks))
            sub_cum[j] = np.cumsum(h / h.sum())

        vals_out, cum_out, off_out = [], [], [0]
        for j in range(ng):
            v, ww = [], []
            for si, k in enumerate(ks):
                if k.disc_vals is None or T[si, j] <= 0:
                    continue
                a, b = k.disc_off[j], k.disc_off[j + 1]
                if b <= a:
                    continue
                c = k.disc_cum[a:b]
                v.append(k.disc_vals[a:b])
                ww.append(T[si, j] * np.diff(np.hstack([0.0, c])))
            if v:
                vv = np.concatenate(v); wwv = np.concatenate(ww)
                u, inv = np.unique(vv, return_inverse=True)
                wu = np.bincount(inv, weights=wwv)
                vals_out.append(u); cum_out.append(np.cumsum(wu / wu.sum()))
            off_out.append(off_out[-1] + (vals_out[-1].size if v else 0))
        disc_vals = np.concatenate(vals_out) if vals_out else np.zeros(0)
        disc_cum = np.concatenate(cum_out) if cum_out else np.zeros(0)

        md = dict(metadata or {}); md["mixed_from"] = len(ks)
        return cls(edges, R, N_cum, q_dep, sub_cum, counts, md,
                   disc_vals=disc_vals, disc_cum=disc_cum,
                   disc_off=np.array(off_out, int))

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
