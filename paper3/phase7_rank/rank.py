"""Paper III E1 (was P8): how many independent redistribution modes are there?

Reframed from an optimization question to a physics one. Nobody needs to save
14 kB. The question is whether the redistribution operator of a lanthanide
forest -- built from 10^3-10^6 atomic pathways -- actually has only a few
macroscopic modes, because that would explain WHY these enormous atomic
networks compress at all (F25, F27), rather than merely recording that they do.

Two operators, and the distinction matters:

  R      the ENERGY matrix. Rows sum to 1 - q_dep, and q_dep goes NEGATIVE
         under blueward fluorescence, so R is not row-stochastic.
  N_row  the PHOTON matrix, recovered by differencing N_cum. Row-stochastic on
         populated rows, and the object transport actually samples.

They have different spectra, so "the rank of the kernel" is ambiguous until you
say which. Both are reported.

Empty rows are excluded throughout: they are identically zero and would inject
spurious zero singular values (a material fraction of the small-N_g kernels).

Usage: python rank.py [--ion laII|ceII|ndII|all] [--groups 32] [--transport]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase0_reference"):
    sys.path.insert(0, str(p))
from forest_mc import band_ratio, run_mc
from run_forest import R_CORE, R_OUT, T_EXP, T_CORE, nu_of
from redistribution import RedistributionKernel
from reference import BANDS, SEEDS, build_atom

KDIR = ROOT / "paper3/phase1_groups"


def photon_matrix(k):
    """Row-stochastic photon operator from the stored row cumulatives."""
    ng = k.n_groups
    return np.diff(np.hstack([np.zeros((ng, 1)), k.N_cum]), axis=1)


def spectrum(M, live):
    """Singular values of the populated block, plus effective-dimension measures.

    Two ways of saying "how many modes", because they answer different
    questions: the participation ratio is a smooth count that a single dominant
    mode drives to ~1, while k90/k99 say how many you must actually keep.
    """
    A = M[live][:, ~np.all(M[live] == 0, axis=0)] if live.any() else M[:0]
    if A.size == 0 or min(A.shape) == 0:
        return dict(n_rows=0)
    s = np.linalg.svd(A, compute_uv=False)
    s2 = s ** 2
    tot = s2.sum()
    if tot <= 0:
        return dict(n_rows=int(A.shape[0]))
    p = s2 / tot
    cum = np.cumsum(p)
    nz = p[p > 0]
    return dict(n_rows=int(A.shape[0]), n_cols=int(A.shape[1]),
                sigma=[float(x) for x in s[:16]],
                sigma1_frac=float(p[0]),
                participation_ratio=float(tot ** 2 / (s2 ** 2).sum()),
                spectral_entropy=float(-(nz * np.log(nz)).sum()),
                exp_entropy=float(np.exp(-(nz * np.log(nz)).sum())),
                k90=int(np.searchsorted(cum, 0.90) + 1),
                k99=int(np.searchsorted(cum, 0.99) + 1),
                rank_numeric=int(np.linalg.matrix_rank(A)))


def nmf(V, k, iters=600, seed=0):
    """Non-negative factorization V ~ W H by Lee-Seung multiplicative updates.

    SVD is the WRONG decomposition here and using it was a mistake worth
    recording: truncating a stochastic matrix produces negative entries, and
    clipping them to zero then renormalizing is a violent nonlinear operation
    that destroys the distribution -- it gave non-monotone transport errors of
    several hundred per cent, which is a property of the truncation and not of
    the physics.

    NMF is both well-posed and the physically meaningful model: with H's rows
    normalized, each input group's exit distribution becomes a MIXTURE OF k
    ARCHETYPAL DISTRIBUTIONS, which is exactly the "are there a few macroscopic
    redistribution modes" question. Non-negativity holds by construction, so no
    clipping is needed.
    """
    rng = np.random.default_rng(seed)
    m, n = V.shape
    scale = np.sqrt(V.mean() / k) if V.mean() > 0 else 1.0
    W = rng.uniform(0.5, 1.5, (m, k)) * scale
    H = rng.uniform(0.5, 1.5, (k, n)) * scale
    eps = 1e-12
    for _ in range(iters):
        H *= (W.T @ V) / (W.T @ W @ H + eps)
        W *= (V @ H.T) / (W @ H @ H.T + eps)
    return W, H


def truncate(kern, k, seed=0):
    """Rank-k kernel that transport can still sample, via NMF.

    The within-group discrete exit tables are untouched -- this experiment asks
    about the GROUP-TO-GROUP structure, not the within-group placement.
    """
    N = photon_matrix(kern)
    live = ~kern.empty_rows
    V = N[live]
    kk = min(k, min(V.shape))
    W, H = nmf(V, kk, seed=seed)
    Nk = W @ H
    rs = Nk.sum(axis=1, keepdims=True)
    Nk = np.where(rs > 0, Nk / np.where(rs > 0, rs, 1.0), 0.0)
    # how well the factorization reproduced the operator, so a transport error
    # can be attributed to the approximation rather than to sensitivity
    rec = float(np.abs(Nk - V).sum(axis=1).mean())
    Nfull = np.zeros_like(N); Nfull[live] = Nk
    kern_k = RedistributionKernel(kern.edges, kern.R, np.cumsum(Nfull, axis=1),
                                  kern.q_dep, kern.sub_cum, kern.counts,
                                  dict(kern.metadata, truncated_to=int(kk)),
                                  disc_vals=kern.disc_vals, disc_cum=kern.disc_cum,
                                  disc_off=kern.disc_off)
    return kern_k, rec


def bootstrap_entropy(ion, ng, n_boot, rng):
    """Resample the event log and rebuild -- error bars with no transport."""
    f = ROOT / f"paper3/phase0_reference/reference_events_{ion}.npz"
    if not f.exists():
        return None
    d = np.load(f); nu_in, nu_out = d["nu_in"], d["nu_out"]
    lo, hi = None, None
    ref = RedistributionKernel.load(KDIR / f"kernel_{ion}_ng{ng}.npz")
    lo, hi = ref.edges[0], ref.edges[-1]
    out = []
    for _ in range(n_boot):
        i = rng.integers(0, nu_in.size, nu_in.size)
        k = RedistributionKernel.from_branching_mc(nu_in[i], nu_out[i],
                                                   np.ones(i.size), ng,
                                                   nu_lo=lo, nu_hi=hi)
        out.append(spectrum(photon_matrix(k), ~k.empty_rows)["exp_entropy"])
    return dict(mean=float(np.mean(out)), std=float(np.std(out, ddof=1)))


def main(ions, ng, do_transport, n_boot):
    rng = np.random.default_rng(0)
    out = {"ng": ng, "ions": {}}
    for ion in ions:
        path = KDIR / f"kernel_{ion}_ng{ng}.npz"
        if not path.exists():
            print(f"  {ion}: no kernel at N_g = {ng}, skipping"); continue
        kern = RedistributionKernel.load(path)
        live = ~kern.empty_rows
        e_spec = spectrum(kern.R, live)
        p_spec = spectrum(photon_matrix(kern), live)
        row = {"n_groups": ng, "n_live_rows": int(live.sum()),
               "n_empty_rows": int((~live).sum()),
               "energy": e_spec, "photon": p_spec}
        print(f"{ion}  N_g={ng}  live rows {int(live.sum())}/{ng}")
        for tag, sp in (("energy R", e_spec), ("photon N", p_spec)):
            print(f"    {tag:9s} sigma1 {100*sp['sigma1_frac']:5.1f}%  "
                  f"PR {sp['participation_ratio']:5.2f}  expH {sp['exp_entropy']:5.2f}  "
                  f"k90 {sp['k90']:3d}  k99 {sp['k99']:3d}  rank {sp['rank_numeric']:3d}")
        if n_boot:
            b = bootstrap_entropy(ion, ng, n_boot, rng)
            if b:
                row["photon"]["exp_entropy_boot"] = b
                print(f"    bootstrap expH  {b['mean']:.2f} +- {b['std']:.2f} ({n_boot} resamples)")

        if do_transport:
            atom, _, _ = build_atom(ion)
            lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
            ref = []
            for s in SEEDS:
                r = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 300000,
                           "sobolev_branch", seed=s, t_core=T_CORE)
                ref.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0]
                            for b, w in BANDS.items()})
            refm = {b: float(np.mean([x[b] for x in ref])) for b in BANDS}
            row["ref_bands"] = refm; row["truncated"] = {}
            for k in (1, 2, 3, 4, 6, 8, 12, 16):
                if k > ng:
                    break
                kt, rec = truncate(kern, k)
                rows = []
                for s in SEEDS:
                    r = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 300000,
                               "sobolev_group", seed=s, t_core=T_CORE, kernel=kt)
                    rows.append({b: band_ratio(r, *nu_of(*w), weight="energy")[0]
                                 for b, w in BANDS.items()})
                dF = {b: float(np.mean([x[b] for x in rows]) / refm[b] - 1) for b in BANDS}
                worst = max(abs(v) for v in dF.values())
                row["truncated"][str(k)] = {"dF": dF, "worst": worst,
                                            "row_L1_reconstruction": rec}
                print(f"    rank {k:2d}: worst |dF_b| {100*worst:6.2f}%   "
                      f"row-L1 reconstruction {rec:.4f}", flush=True)
        out["ions"][ion] = row
    (HERE / f"rank_ng{ng}.json").write_text(json.dumps(out, indent=1))
    print(f"wrote rank_ng{ng}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="all")
    ap.add_argument("--groups", type=int, default=32)
    ap.add_argument("--transport", action="store_true")
    ap.add_argument("--boot", type=int, default=0)
    a = ap.parse_args()
    ions = ["laII", "ceII", "ndII"] if a.ion == "all" else [a.ion]
    main(ions, a.groups, a.transport, a.boot)
