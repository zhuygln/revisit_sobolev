"""Paper III P9: is composition cheap?

The question the plan asks: does an opacity-weighted composition rule

    R_mix[i] = sum_s w[i,s] R_s[i]

reproduce a kernel trained explicitly on the blend? If it does, a kernel
library is per-ION -- every mixture is assembled from it for free -- and the
table stops growing with the number of compositions.

The weights come from the blend's OPACITY ALONE, never from a blend run:
w[i,s] is species s's share of the interaction probability sum(1 - e^-tau) in
group i. That is what makes it a prediction rather than a fit.

Legs, all against the SAME blend branch reference (3 seeds x 2e6,
energy-weighted bands), all on identical group edges spanning the blend's
opacity extent:

  explicit   kernel trained on the blend's own event log      (the ceiling)
  mixed      RedistributionKernel.mix(per-ion kernels, w)     (the rule)
  la_only    the La II kernel alone                           (control)
  ce_only    the Ce II kernel alone                           (control)

The controls matter: if the blend is simply dominated by Ce (E10 found
eps_best tracks the dominant forest), then "ce_only" already works and the
mixture rule is untested. The rule earns its keep only by beating them.

Usage: python mixture.py [--groups 32] [--n 2000000]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase0_reference"):
    sys.path.insert(0, str(p))
from sobolev.constants import C
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of
from e9_ceII import CE_LEV, CE_TR, ce_n_ion
from redistribution import RedistributionKernel
from reference import BANDS, SEEDS


def composition_weights(atom, edges):
    """w[i,s]: species s's share of sum(1 - e^-tau) among the opacity lines of
    group i. Groups with no opacity line get an even split (nothing is
    absorbed there, so the row is never used)."""
    ion = atom.ion_of_line[atom.op_idx]
    ns = int(ion.max()) + 1
    gi = np.clip(np.searchsorted(edges, atom.op_nu, side="right") - 1, 0, edges.size - 2)
    W = np.zeros((edges.size - 1, ns))
    for s in range(ns):
        m = ion == s
        W[:, s] = np.bincount(gi[m], weights=atom.op_p[m], minlength=edges.size - 1)
    tot = W.sum(axis=1, keepdims=True)
    return np.where(tot > 0, W / np.where(tot > 0, tot, 1.0), 1.0 / ns)


def branch_run(atom, lo, hi, n, collect):
    """Branch reference on one atom; returns per-seed bands and the event log."""
    bands, ev_in, ev_out = [], [], []
    for s in SEEDS:
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_branch",
                     seed=s, t_core=T_CORE, collect_events=collect)
        bands.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
        if collect:
            e = res["events"]; ev_in.append(e[0]); ev_out.append(e[1])
    if not collect:
        return bands, None, None
    return bands, np.concatenate(ev_in), np.concatenate(ev_out)


def group_run(atom, lo, hi, n, kernel):
    rows = []
    for s in SEEDS:
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_group",
                     seed=s, t_core=T_CORE, kernel=kernel)
        rows.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
    return rows


def main(ng, n):
    n_la = float(np.load(ROOT / "experiments/laII_forest/forest_lines.npz")["n_ion"])
    n_ce, _ = ce_n_ion()
    t0 = time.time()
    blend = ForestAtom.from_gsi_blend([(LEV, TR, n_la), (CE_LEV, CE_TR, n_ce)],
                                      T_SHELL, T_EXP, tau_min=1e-3)
    # every kernel lives on the blend's grid, so the rows are mixable
    lo, hi = blend.op_nu.min() * 0.995, blend.op_nu.max() * 1.005
    la = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_la, T_EXP, tau_min=1e-3)
    ce = ForestAtom.from_gsi(CE_LEV, CE_TR, T_SHELL, n_ce, T_EXP, tau_min=1e-3)
    print(f"blend {blend.n_opacity} opacity lines (La {la.n_opacity}, Ce {ce.n_opacity}), "
          f"built in {time.time()-t0:.0f}s", flush=True)

    edges = np.geomspace(lo, hi, ng + 1)
    w = composition_weights(blend, edges)
    occupied = np.bincount(np.clip(np.searchsorted(edges, blend.op_nu, side="right") - 1,
                                   0, ng - 1), minlength=ng) > 0
    print(f"N_g = {ng}: {occupied.sum()} groups carry opacity; "
          f"mean Ce share over those {w[occupied, 1].mean():.3f}", flush=True)

    # per-ion kernels, each trained on its own forest but on the SHARED grid
    kern = {}
    for tag, atom in (("laII", la), ("ceII", ce)):
        _, nin, nout = branch_run(atom, lo, hi, n, True)
        kern[tag] = RedistributionKernel.from_branching_mc(
            nin, nout, np.ones(nin.size), ng, nu_lo=lo, nu_hi=hi,
            metadata={"ion": tag, "shared_grid": True})
        print(f"  {tag} kernel: {nin.size} events", flush=True)

    ref_bands, nin_b, nout_b = branch_run(blend, lo, hi, n, True)
    explicit = RedistributionKernel.from_branching_mc(
        nin_b, nout_b, np.ones(nin_b.size), ng, nu_lo=lo, nu_hi=hi,
        metadata={"ion": "blend", "shared_grid": True})
    print(f"  blend reference: {nin_b.size} events", flush=True)

    mixed = RedistributionKernel.mix([kern["laII"], kern["ceII"]], w,
                                     metadata={"rule": "opacity-weighted"})
    assert mixed.validate_energy() < 1e-10

    ref = {b: float(np.mean([r[b] for r in ref_bands])) for b in BANDS}
    out = {"ng": ng, "n": n, "seeds": list(SEEDS), "n_ion": {"LaII": n_la, "CeII": n_ce},
           "n_opacity": {"blend": int(blend.n_opacity), "laII": int(la.n_opacity),
                         "ceII": int(ce.n_opacity)},
           "ref_bands": ref, "mean_ce_share": float(w[occupied, 1].mean()),
           "row_l1": {}, "legs": {}}
    for tag, k in (("explicit", explicit), ("mixed", mixed),
                   ("la_only", kern["laII"]), ("ce_only", kern["ceII"])):
        rows = group_run(blend, lo, hi, n, k)
        dF = {b: float(np.mean([r[b] for r in rows]) / ref[b] - 1) for b in BANDS}
        worst = max(abs(v) for v in dF.values())
        live = ~explicit.empty_rows & ~k.empty_rows
        out["row_l1"][tag] = float(np.abs(k.R[live] - explicit.R[live]).sum(axis=1).mean())
        out["legs"][tag] = {"dF": dF, "worst": worst}
        print(f"  {tag:9s} worst |dF_b| {100*worst:6.2f}%   rowL1 vs explicit "
              f"{out['row_l1'][tag]:.4f}   " +
              " ".join(f"{b}={100*dF[b]:+.1f}%" for b in BANDS), flush=True)

    # Two independent questions, reported separately -- collapsing them into
    # one pass/fail hides which half failed.
    m, e = out["legs"]["mixed"]["worst"], out["legs"]["explicit"]["worst"]
    best_ctl = min(out["legs"]["la_only"]["worst"], out["legs"]["ce_only"]["worst"])
    out["verdict"] = {
        "beats_controls": bool(m < best_ctl),
        "gain_over_best_control": float(best_ctl / m) if m > 0 else float("inf"),
        "matches_explicit": bool(m < max(1.5 * e, e + 0.01)),
        "cost_over_explicit": float(m / e) if e > 0 else float("inf"),
        "meets_strong_bar": bool(m < 0.05),
    }
    v = out["verdict"]
    print(f"GATE P9: beats controls {v['beats_controls']} "
          f"({v['gain_over_best_control']:.1f}x better than the best single ion); "
          f"matches explicit {v['matches_explicit']} "
          f"({v['cost_over_explicit']:.1f}x its error); "
          f"|dF_b| < 5% {v['meets_strong_bar']}")
    (HERE / f"mixture_ng{ng}.json").write_text(json.dumps(out, indent=1))
    print(f"wrote mixture_ng{ng}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--groups", type=int, default=32)
    ap.add_argument("--n", type=float, default=2e6)
    a = ap.parse_args(); main(a.groups, int(a.n))
