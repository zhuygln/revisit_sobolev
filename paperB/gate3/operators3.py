"""G3's derived operators, all on the frozen per-ion support.

  interpolate_whole(ka, kb, lam)  the whole effective operator interpolated
      linearly in the frozen coordinate between two endpoint kernels on
      identical edges: the energy and photon transition rows (a row populated
      at both endpoints is interpolated; at one endpoint only, that endpoint's
      row is used and counted; at neither, it stays empty), the deposition
      channel q_dep, and the conditional exit-frequency weights per output
      group on the UNION of the two endpoints' exit lines, a line absent at
      one endpoint zero-filled there. Built from the two endpoint states'
      build events only -- never from the interior state's.
  interpolate_matrix(ka, kb, lam, tables_from)  the diagnostic: the
      interpolated transition matrix carrying one endpoint's exit tables.
  support_diag(...)  a transform for the fine kernel leg that records how many
      of the state's absorbed frequencies fall outside the frozen support
      ("outside fixed support", clipped into an edge group) -- reported
      separately from "row never trained" (an endpoint-empty row hit at the
      state, which transport counts as coherent fallback).
"""
import numpy as np

from redistribution import RedistributionKernel


def _rows_from_cum(cum):
    return np.diff(np.hstack([np.zeros((cum.shape[0], 1)), cum]), axis=1)


def _group_weights(k, j, which):
    """(vals, relative weights) of output group j's exit table."""
    a, b = k.disc_off[j], k.disc_off[j + 1]
    if b <= a:
        return np.zeros(0), np.zeros(0)
    cum = k.disc_cum if which == "photon" else (k.disc_cum_E if k.disc_cum_E is not None else k.disc_cum)
    return k.disc_vals[a:b], np.diff(np.hstack([0.0, cum[a:b]]))


def _check_edges(ka, kb):
    if ka.edges.shape != kb.edges.shape or not np.allclose(ka.edges, kb.edges, rtol=0, atol=0):
        raise ValueError("interpolation needs kernels on the identical frozen support")


def interpolate_whole(ka, kb, lam, coord=None, endpoints=None):
    _check_edges(ka, kb)
    lam = float(lam); ng = ka.n_groups
    live_a, live_b = ~ka.empty_rows, ~kb.empty_rows
    both = live_a & live_b; only_a = live_a & ~live_b; only_b = live_b & ~live_a
    R = np.zeros_like(ka.R); q = np.zeros(ng); counts = np.zeros(ng)
    Na, Nb = _rows_from_cum(ka.N_cum), _rows_from_cum(kb.N_cum)
    N = np.zeros_like(Na)
    for src, w_a, w_b in ((both, 1.0 - lam, lam), (only_a, 1.0, 0.0), (only_b, 0.0, 1.0)):
        R[src] = w_a * ka.R[src] + w_b * kb.R[src]
        q[src] = w_a * ka.q_dep[src] + w_b * kb.q_dep[src]
        N[src] = w_a * Na[src] + w_b * Nb[src]
        counts[src] = w_a * ka.counts[src] + w_b * kb.counts[src]
    with np.errstate(invalid="ignore", divide="ignore"):
        tot = N.sum(axis=1, keepdims=True)
        N_cum = np.where(tot > 0, np.cumsum(N, axis=1) / np.where(tot > 0, tot, 1.0), 0.0)
    sub_cum = (1.0 - lam) * ka.sub_cum + lam * kb.sub_cum
    # the conditional exit spectrum on the union of exit lines, zero-filled
    vals, cum_p, cum_e, off = [], [], [], [0]
    n_a = n_b = n_union = 0
    for j in range(ng):
        va, wa = _group_weights(ka, j, "photon"); vb, wb = _group_weights(kb, j, "photon")
        _, ea = _group_weights(ka, j, "energy"); _, eb = _group_weights(kb, j, "energy")
        n_a += va.size; n_b += vb.size
        u = np.unique(np.concatenate([va, vb]))
        if u.size == 0:
            off.append(off[-1]); continue
        def fill(v, w):
            out = np.zeros(u.size); out[np.searchsorted(u, v)] = w; return out
        wp = (1.0 - lam) * fill(va, wa) + lam * fill(vb, wb)
        we = (1.0 - lam) * fill(va, ea) + lam * fill(vb, eb)
        if wp.sum() <= 0 or we.sum() <= 0:                    # one side empty: take the other whole
            wp = fill(va, wa) if wa.size else fill(vb, wb); we = fill(va, ea) if ea.size else fill(vb, eb)
        vals.append(u); cum_p.append(np.cumsum(wp) / wp.sum()); cum_e.append(np.cumsum(we) / we.sum())
        off.append(off[-1] + u.size); n_union += u.size
    md = dict(ka.metadata, transform=dict(kind="interp_whole", lam=lam, coord=coord, endpoints=list(endpoints or []),
                                          rows_both=int(both.sum()), rows_one_endpoint=int((only_a | only_b).sum()),
                                          n_exit_a=int(n_a), n_exit_b=int(n_b), n_exit_union=int(n_union),
                                          archetypes=int(ng), n_params=int(ng * ng)))
    return RedistributionKernel(ka.edges, R, N_cum, q, sub_cum, counts, md,
                                disc_vals=np.concatenate(vals) if vals else np.zeros(0),
                                disc_cum=np.concatenate(cum_p) if cum_p else np.zeros(0),
                                disc_off=np.array(off, int),
                                disc_cum_E=np.concatenate(cum_e) if cum_e else np.zeros(0))


def interpolate_matrix(ka, kb, lam, tables_from, coord=None, endpoints=None):
    """The diagnostic: interpolated transition rows, one endpoint's tables."""
    whole = interpolate_whole(ka, kb, lam, coord, endpoints)
    src = ka if tables_from == "a" else kb
    md = dict(whole.metadata); md["transform"] = dict(whole.metadata["transform"], kind="interp_matrix", tables_from=tables_from)
    k = RedistributionKernel(src.edges, whole.R, whole.N_cum, whole.q_dep, src.sub_cum, whole.counts, md,
                             disc_vals=src.disc_vals, disc_cum=src.disc_cum, disc_off=src.disc_off, disc_cum_E=src.disc_cum_E)
    return k


def support_diag():
    """Transform for a kernel-only leg: records the state's absorbed
    frequencies outside the frozen support (clipped into an edge group)."""
    def f(kern, ev):
        nu = ev["nu_in"]; e = ev["w_in"] * nu
        lo, hi = kern.edges[0], kern.edges[-1]
        below, above = nu < lo, nu > hi
        kern.metadata["transform"] = dict(kind="support", clipped_frac=float((below | above).mean()),
                                          clipped_energy_frac=float(e[below | above].sum() / e.sum()) if e.sum() > 0 else 0.0,
                                          below_frac=float(below.mean()), above_frac=float(above.mean()),
                                          nu_lo=float(lo), nu_hi=float(hi), n_events=int(nu.size))
        return kern
    return f
