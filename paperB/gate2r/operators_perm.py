"""G2R -- the frequency-adjacency ablation (paperB/prl_gate.md "G2R").

The G2 local family is the 128-group fine energy matrix block-coarsened
into N_g CONTIGUOUS index blocks (gate2/operators.py::local); its blocks
are contiguous in frequency because the 128 groups are. The ablation keeps
everything -- the fine operator, the build events, the block count, the
parameter count N_g^2, the 128-group exit tables, the evaluation seeds,
the transport -- and destroys only the physical adjacency: the 128 indices
are permuted, the SAME block coarse-graining is applied, and the
permutation is undone. Ordering 0 is the physical order and must reproduce
`local(N_g)`; orderings 1..31 are fixed random permutations.

The block coarse-graining is written here directly on the fine matrix and
its energy-in vector. It equals `local(N_g)`'s route through
`from_branching_mc` on the coarse edges because the coarse edges are every
blk-th fine edge (both log-spaced over the same support), so
coarse.R[I,J] = sum_{i in I, j in J} E_i R_ij / sum_{i in I} E_i -- the
identity `tests/test_paperB_gate2r.py` pins against the frozen G2 records.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT / "paperB/gate2", ROOT / "paper3"):
    sys.path.insert(0, str(p))
from operators import energy_in, _tv_rows      # noqa: E402  (G2's, unchanged)

N_FINE = 128
N_SCRAMBLED = 31
PERM_SEED0 = 1000


def orderings(n_scrambled=N_SCRAMBLED, n_fine=N_FINE, seed0=PERM_SEED0):
    """[(perm_id, perm, seed)]: id 0 the physical order (seed None), then
    the fixed random permutations default_rng(seed0 + m).permutation(n_fine)."""
    out = [(0, np.arange(n_fine), None)]
    for m in range(1, n_scrambled + 1):
        out.append((m, np.random.default_rng(seed0 + m).permutation(n_fine), seed0 + m))
    return out


def block_coarsen(R, E, n_coarse, perm=None):
    """The fine energy matrix R (nf x nf) with energy-in E (nf) coarsened into
    n_coarse equal index blocks of the order `perm` (None: the physical
    order), expanded back onto the fine groups exactly as
    gate2/operators.py::local does: every fine row of block I carries the
    block's row, and each block column's mass is spread over its fine groups
    by the fine exit-energy marginal within it (over all input rows)."""
    R = np.asarray(R, float); E = np.asarray(E, float); nf = R.shape[0]
    if nf % n_coarse:
        raise ValueError(f"{nf} fine groups are not a multiple of {n_coarse}")
    perm = np.arange(nf) if perm is None else np.asarray(perm, int)
    if perm.shape != (nf,) or not np.array_equal(np.sort(perm), np.arange(nf)):
        raise ValueError("perm must be a permutation of the fine group indices")
    flow = E[:, None] * R
    Fp = flow[np.ix_(perm, perm)]; Ep = E[perm]
    blk = nf // n_coarse
    Rp = np.zeros_like(R)
    for J in range(n_coarse):
        cols = slice(J * blk, (J + 1) * blk)
        fj = Fp[:, cols].sum(axis=0)
        share = fj / fj.sum() if fj.sum() > 0 else np.full(blk, 1.0 / blk)
        for I in range(n_coarse):
            rows = slice(I * blk, (I + 1) * blk)
            EI = Ep[rows].sum()
            cIJ = Fp[rows, cols].sum() / EI if EI > 0 else 0.0
            Rp[rows, cols] = cIJ * share[None, :]
    inv = np.argsort(perm)
    return Rp[np.ix_(inv, inv)]


def permuted_local(n_coarse, perm, perm_id, perm_seed=None):
    """The `transform` for run_legs: (kernel, events) -> the operator of
    ordering `perm_id` at block count n_coarse, on the kernel's own tables."""
    perm = np.asarray(perm, int)

    def f(kern, ev):
        E = energy_in(kern, ev)
        R_loc = block_coarsen(kern.R, E, n_coarse, perm)
        md = dict(transform=dict(kind="local_perm", n_coarse=int(n_coarse), perm_id=int(perm_id),
                                 perm_seed=None if perm_seed is None else int(perm_seed), physical=bool(perm_id == 0),
                                 archetypes=int(n_coarse), n_params=int(n_coarse ** 2),
                                 in_sample_tv=_tv_rows(R_loc, kern.R, E)))
        return kern.with_matrix(R_loc, md)
    return f
