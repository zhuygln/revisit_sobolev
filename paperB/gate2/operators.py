"""The G2 operator families, every one a 128-group energy matrix on G1's
edges with the SAME 128-group discrete exit tables (from the build
reference), so the families differ only in the group-to-group matrix.

  local(N)     the N-group coarse matrix expanded onto the fine groups, each
               coarse column's mass spread by the fine exit-energy marginal
               within it -- sampling-equivalent to G1's A2_ng{N} kernel
  nmf_rank(k)  the rank-k non-negative factorisation W H of the fine energy
               matrix's populated rows (Lee-Seung multiplicative updates, the
               rule of paper3/phase7_rank/rank.py::nmf, copied here so that
               module's Paper III imports are not needed), rows rescaled to
               the original row sums; k archetypal exit distributions with
               free mixing weights per fine input group
  truncate(f)  the full matrix with each output group's exit table cut to
               the highest-weight lines carrying fraction f of its exit energy

Each returns a callable (kernel, events) -> kernel for `run_legs`' spec key
`transform`; the derived kernel's metadata["transform"] records what was done.
"""
import numpy as np

from redistribution import RedistributionKernel


def nmf(V, k, iters=600, seed=0, tail=100):
    """V ~ W H, non-negative, Lee-Seung. Returns W, H and the convergence
    record: the relative Frobenius error at the end and its relative change
    over the last `tail` iterations."""
    rng = np.random.default_rng(seed)
    m, n = V.shape
    scale = np.sqrt(V.mean() / k) if V.mean() > 0 else 1.0
    W = rng.uniform(0.5, 1.5, (m, k)) * scale
    H = rng.uniform(0.5, 1.5, (k, n)) * scale
    eps = 1e-12
    norm = np.linalg.norm(V)
    err_tail = None
    for it in range(iters):
        H *= (W.T @ V) / (W.T @ W @ H + eps)
        W *= (V @ H.T) / (W @ H @ H.T + eps)
        if it == iters - tail - 1:
            err_tail = np.linalg.norm(V - W @ H) / norm
    err = float(np.linalg.norm(V - W @ H) / norm)
    change = float(abs(err_tail - err) / max(err, 1e-300)) if err_tail is not None else float("nan")
    return W, H, dict(rel_frobenius=err, rel_change_tail=change, iters=iters, seed=seed)


def energy_in(kern, ev):
    """Energy absorbed per fine row from the events (h omitted)."""
    gi = kern.group_index(ev["nu_in"])
    return np.bincount(gi, weights=ev["w_in"] * ev["nu_in"], minlength=kern.n_groups)


def _tv_rows(Ra, Rb, W):
    rows = W > 0
    d = 0.5 * np.abs(Ra[rows] - Rb[rows]).sum(axis=1)
    return float((d * W[rows]).sum() / W[rows].sum())


def local(n_coarse):
    def f(kern, ev):
        nf = kern.n_groups
        if nf % n_coarse:
            raise ValueError(f"{nf} fine groups are not a multiple of {n_coarse}")
        coarse = RedistributionKernel.from_branching_mc(ev["nu_in"], ev["nu_out"], ev["w_in"], n_coarse,
                                                        nu_lo=kern.edges[0], nu_hi=kern.edges[-1], w_out=ev["w_out"])
        E = energy_in(kern, ev)
        flow = E[:, None] * kern.R                            # energy flow fine i -> fine j
        blk = nf // n_coarse
        R_loc = np.zeros_like(kern.R)
        for J in range(n_coarse):
            cols = slice(J * blk, (J + 1) * blk)
            fj = flow[:, cols].sum(axis=0)
            share = fj / fj.sum() if fj.sum() > 0 else np.full(blk, 1.0 / blk)
            for I in range(n_coarse):
                rows = slice(I * blk, (I + 1) * blk)
                R_loc[rows, cols] = coarse.R[I, J] * share[None, :]
        md = dict(transform=dict(kind="local", n_coarse=int(n_coarse), archetypes=int(n_coarse), n_params=int(n_coarse ** 2),
                                 in_sample_tv=_tv_rows(R_loc, kern.R, E)))
        return kern.with_matrix(R_loc, md)
    return f


def nmf_rank(k, iters=600, seed=0):
    def f(kern, ev):
        live = ~kern.empty_rows
        V = kern.R[live]
        kk = min(int(k), min(V.shape))
        W, H, conv = nmf(V, kk, iters=iters, seed=seed)
        Rk = W @ H
        rs_old = V.sum(axis=1, keepdims=True); rs_new = Rk.sum(axis=1, keepdims=True)
        Rk = np.where(rs_new > 0, Rk * rs_old / np.where(rs_new > 0, rs_new, 1.0), 0.0)
        R_full = np.zeros_like(kern.R); R_full[live] = Rk
        E = energy_in(kern, ev)
        md = dict(transform=dict(kind="nmf", k=kk, archetypes=kk, n_params=int(2 * kern.n_groups * kk),
                                 in_sample_tv=_tv_rows(R_full, kern.R, E), **conv))
        return kern.with_matrix(R_full, md)
    return f


def truncate(frac):
    def f(kern, ev):
        n_total = int(kern.disc_vals.size)
        k_new, kept = kern.truncate_exits(frac)
        k_new.metadata["transform"] = dict(kind="truncate", f=float(frac), n_exit_kept=kept, n_exit_total=n_total,
                                           archetypes=int(kern.n_groups), n_params=int(kern.n_groups ** 2))
        return k_new
    return f
