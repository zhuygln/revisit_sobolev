"""Paper III P11: the target architecture, kappa_grouped + R_ij.

Everything up to here held the opacity FIXED -- Sobolev line-by-line on both
sides of every comparison -- so that any error was redistribution compression
alone. That was deliberate (the plan: "do NOT combine with expansion opacity
until the redistribution approximation is validated"). It is validated now
(F25, F27), so this is where the two halves meet.

A production code cannot afford line-by-line Sobolev opacity; that is the
entire reason expansion opacity exists. So the question is whether the kernel
survives being put on top of a *grouped* opacity, and which grouping.

Three opacity treatments, ALL carrying the same R_ij, so the opacity is the
only thing that differs:

  sobolev_group     exact per-line tau, exact resonance frequency
  binned_group      bin optical depth sum(tau)        -- exact attenuation
                    (F12: optical depths add), frequency at bin resolution
  expansion_group   bin optical depth sum(1 - e^-tau) -- the Poisson
                    substitution of F15, frequency at bin resolution

sobolev -> binned isolates the cost of losing the resonance frequency.
binned -> expansion isolates the cost of the substitution itself.

Reference is the explicit physics, sobolev_branch, on the same atom; and
expansion_branch is carried as the F24 comparison -- the branching-aware
Poisson closure that is +21% on La II and +113% on Ce II. The live question
is whether kappa_grouped + R_ij inherits that density limit. If it does, the
architecture does not work on dense forests and P12 must not be built on it.

Usage: python opacity.py [--ion laII|ceII|ndII] [--groups 32] [--n 2000000]
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

# SEDONA's production transport bin, the grid Paper I's F13 bin-width study
# swept; the expansion legs are the closure as codes actually run it.
DNU = 4.17e-5
# E2: memory depth. How much of the recent ORDERED resonance history must group
# transport keep? m = 1 is the at-resonance skip; m > 1 is the open question.
MEM = (1, 2, 4, 8, 16)


def legs(atom, lo, hi, n, kernel, modes):
    out = {}
    for tag in modes:
        mode = tag
        t0 = time.time()
        rows, inter = [], []
        for s in SEEDS:
            kw = {"kernel": kernel} if "_group" in mode else {}
            if "+m" in mode:
                mode, m = mode.split("+m")
                kw["line_memory"] = int(m)
            res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode,
                         seed=s, t_core=T_CORE, dnu_over_nu=DNU, **kw)
            rows.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
            inter.append(float(res["n_events"].mean()))
        out[tag] = {"bands": {b: float(np.mean([r[b] for r in rows])) for b in BANDS},
                     "events_per_packet": float(np.mean(inter)),
                     "wall_s": time.time() - t0}
        print(f"  {tag:20s} " + " ".join(f"{b}={out[tag]['bands'][b]:.4f}" for b in BANDS)
              + f"  [{out[tag]['wall_s']:.0f}s, {out[tag]['events_per_packet']:.3f} ev/pkt]",
              flush=True)
    return out


def main(ion, ng, n):
    atom, n_ion, _ = build_atom(ion)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    print(f"{ion}: {atom.n_opacity} opacity lines, n_ion {n_ion:.1f}, "
          f"tau_max {atom.op_tau.max():.1f}, N_g = {ng}, dnu/nu = {DNU:g}", flush=True)

    # the kernel: trained once, from the explicit branch physics, and reused
    # unchanged by all three opacity treatments
    ev_in, ev_out = [], []
    for s in SEEDS:
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_branch",
                     seed=s, t_core=T_CORE, collect_events=True)
        e = res["events"]; ev_in.append(e[0]); ev_out.append(e[1])
    nu_in, nu_out = np.concatenate(ev_in), np.concatenate(ev_out)
    kernel = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size),
                                                    ng, nu_lo=lo, nu_hi=hi)
    print(f"  kernel trained on {nu_in.size} events", flush=True)

    res = legs(atom, lo, hi, n, kernel,
               ["sobolev_branch", "expansion_branch", "dual_branch",
                "sobolev_group", "binned_group", "expansion_group", "dual_group"]
               + [f"binned_group+m{m}" for m in MEM]
               + [f"expansion_group+m{m}" for m in MEM])
    ref = res["sobolev_branch"]["bands"]

    out = {"ion": ion, "ng": ng, "n": n, "dnu_over_nu": DNU, "n_ion": n_ion,
           "n_opacity": int(atom.n_opacity), "tau_max": float(atom.op_tau.max()),
           "n_events": int(nu_in.size), "ref_bands": ref, "legs": {}}
    print(f"\n  vs sobolev_branch (the explicit physics):")
    for mode, r in res.items():
        if mode == "sobolev_branch":
            continue
        dF = {b: r["bands"][b] / ref[b] - 1 for b in BANDS}
        worst = max(abs(v) for v in dF.values())
        out["legs"][mode] = {"dF": dF, "worst": worst,
                             "events_per_packet": r["events_per_packet"]}
        print(f"    {mode:18s} worst {100*worst:7.2f}%   "
              + " ".join(f"{b}={100*dF[b]:+7.1f}%" for b in BANDS), flush=True)

    g = out["legs"]
    out["decomposition"] = {
        "redistribution_only": g["sobolev_group"]["worst"],
        "plus_bin_resolution": g["binned_group"]["worst"],
        "plus_poisson": g["expansion_group"]["worst"],
        "dual_bin_group": g["dual_group"]["worst"],
        "memory_sweep": {f"m{m}": {"binned": g[f"binned_group+m{m}"]["worst"],
                                   "expansion": g[f"expansion_group+m{m}"]["worst"]}
                         for m in MEM},
        "poisson_with_exact_exit": g["expansion_branch"]["worst"],
        "dual_with_exact_exit": g["dual_branch"]["worst"],
    }
    d = out["decomposition"]
    print(f"\n  P11 decomposition (worst band):"
          f"\n    R_ij alone                  {100*d['redistribution_only']:7.2f}%"
          f"\n    + bin resolution            {100*d['plus_bin_resolution']:7.2f}%"
          f"\n    + Poisson substitution      {100*d['plus_poisson']:7.2f}%"
          f"\n    + two-quantity bin          {100*d['dual_bin_group']:7.2f}%"
          f"\n  memory depth m (worst band):"
          + "".join(f"\n    m = {m:2d}   binned {100*d['memory_sweep'][f'm{m}']['binned']:7.2f}%"
                    f"   expansion {100*d['memory_sweep'][f'm{m}']['expansion']:7.2f}%"
                    for m in MEM)
          + f"\n  with the exact A*beta exit (opacity alone):"
          f"\n    Poisson  (expansion_branch) {100*d['poisson_with_exact_exit']:7.2f}%  [F24]"
          f"\n    two-quantity (dual_branch)  {100*d['dual_with_exact_exit']:7.2f}%")
    (HERE / f"opacity_{ion}_ng{ng}.json").write_text(json.dumps(out, indent=1))
    print(f"wrote opacity_{ion}_ng{ng}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="laII", choices=["laII", "ceII", "ndII"])
    ap.add_argument("--groups", type=int, default=32)
    ap.add_argument("--n", type=float, default=2e6)
    a = ap.parse_args(); main(a.ion, a.groups, int(a.n))
