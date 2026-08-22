"""E6: the redistribution matrix P(lambda_out | lambda_in) -- direct branching
vs the thermal closure vs TLA at eps_best -- and a difference map.

run_mc returns nu_launch and nu_out_all aligned by packet, so the packet-level
matrix is a 2-D histogram; rows are normalized by the launched energy in the
row so each row reads "where the energy launched at lambda_in emerged".
Escaped-only (core-absorbed packets do not appear). The closure is expected
to return energy toward the LTE emissivity peak (infrared at 3000 K) while
branching carries it through the cascade into the optical.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C, H
from forest_mc import ForestAtom, run_mc
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE


def matrix(res, edges):
    """Row-normalized energy redistribution matrix: M[i, j] = escaped energy in
    out-bin j from launches in in-bin i / launched energy in in-bin i."""
    esc = res["fate"] == 1
    w_in = H * res["nu_launch"]; w_out = H * res["nu_out_all"]
    Hm, _, _ = np.histogram2d(res["nu_launch"][esc], res["nu_out_all"][esc], bins=[edges, edges], weights=w_out[esc])
    row, _ = np.histogram(res["nu_launch"], edges, weights=w_in)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(row[:, None] > 0, Hm / row[:, None], np.nan)


def main(n, seed, eps_best):
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz"); n_ion = float(d["n_ion"])
    atom = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    edges = np.geomspace(lo, hi, 61)
    lam = C / np.sqrt(edges[1:] * edges[:-1]) * 1e8
    legs = {"sobolev_branch": ("sobolev_branch", {}), "expansion_thermal": ("expansion_thermal", {}),
            "sobolev_thermal": ("sobolev_thermal", {}),
            f"expansion_tla_eps{eps_best:g}": ("expansion_tla", {"eps": eps_best}),
            f"sobolev_tla_eps{eps_best:g}": ("sobolev_tla", {"eps": eps_best})}
    M = {}
    for tag, (mode, kw) in legs.items():
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=seed, t_core=T_CORE, **kw)
        M[tag] = matrix(res, edges)
        # diagonal share: energy that emerged within its own launch bin
        diag = np.nanmean(np.diag(M[tag])); print(f"  {tag:26s} mean diagonal (own-bin) share {diag:.3f}")
    np.savez(HERE / "e6_redistribution.npz", edges=edges, **{k: v for k, v in M.items()})

    fig, axs = plt.subplots(1, 4, figsize=(16, 4.2))
    ext = [lam[-1], lam[0], lam[-1], lam[0]]
    for ax, tag in zip(axs[:3], ("sobolev_branch", "expansion_thermal", f"expansion_tla_eps{eps_best:g}")):
        im = ax.imshow(np.log10(np.clip(M[tag], 1e-4, None))[::-1, ::-1], origin="lower", extent=ext, aspect="auto", vmin=-4, vmax=0, cmap="viridis")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(tag, fontsize=9); ax.set_xlabel(r"$\lambda_{\rm out}$ [A]")
        ax.axvspan(3800, 3955, color="w", alpha=0.15)
    axs[0].set_ylabel(r"$\lambda_{\rm in}$ [A]")
    fig.colorbar(im, ax=axs[:3], label=r"$\log_{10}$ P(out | in), energy, row-normalized", pad=0.01)
    ax = axs[3]
    D = np.nan_to_num(M[f"expansion_tla_eps{eps_best:g}"]) - np.nan_to_num(M["sobolev_branch"])
    im2 = ax.imshow(D[::-1, ::-1], origin="lower", extent=ext, aspect="auto", vmin=-0.2, vmax=0.2, cmap="RdBu_r")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(f"TLA(eps={eps_best:g}) - branch", fontsize=9); ax.set_xlabel(r"$\lambda_{\rm out}$ [A]")
    fig.colorbar(im2, ax=ax, pad=0.02)
    for out in (ROOT / "outputs/fig_p2_redistribution.png", ROOT / "docs/figures/fig_p2_redistribution.png"):
        out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote e6_redistribution.npz, fig_p2_redistribution.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e6); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--eps", type=float, default=0.3, help="eps_best from e4_eps_sweep.json")
    a = ap.parse_args(); main(int(a.n), a.seed, a.eps)
