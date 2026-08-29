"""Paper III P5: thermodynamic robustness of the redistribution kernel.

--which src : T_src = 4000/5000/6000/8000 K, material state fixed. The
  6000 K kernel is reused at every T_src (transferability) and compared
  against a kernel recomputed at that T_src (representation vs
  state-transfer error). Prediction, in advance: near-exact transfer -- the
  kernel's rows depend on the radiation field only through the absorbing-
  line mix within a group.
--which gas : T_gas = 2500/3000/4000/5000 K at fixed n_ion; LTE populations
  and every tau recomputed, so the branch chains themselves change. The
  genuine dependence, if there is one.

Every closure run is scored against the SAME-configuration branch
reference (3 seeds x 2e6, energy-weighted bands); kernels train on that
configuration's own event log except the "fixed" kernel, which is
phase1_groups/kernel_laII_ng32.npz (T_src 6000, T_gas 3000, t 1 d).
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase0_reference"):
    sys.path.insert(0, str(p))
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, nu_of
from redistribution import RedistributionKernel
from reference import BANDS, SEEDS, N

NG = 32
FIXED = ROOT / "paper3/phase1_groups/kernel_laII_ng32.npz"


def run_config(atom, lo, hi, t_src, kernels):
    """Branch reference + one closure run per kernel; returns bands."""
    ev_in, ev_out = [], []
    ref = []
    for s in SEEDS:
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, N, "sobolev_branch",
                     seed=s, t_core=t_src, collect_events=True)
        ref.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
        e = res["events"]; ev_in.append(e[0]); ev_out.append(e[1])
    nu_in, nu_out = np.concatenate(ev_in), np.concatenate(ev_out)
    own = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size), NG)
    out = {"ref": ref, "own_kernel_events": int(nu_in.size)}
    for name, kern in list(kernels.items()) + [("own", own)]:
        rows = []
        for s in SEEDS:
            res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, N, "sobolev_group",
                         seed=s, t_core=t_src, kernel=kern)
            rows.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
        out[name] = rows
    return out, own


def summarize(tag, out):
    ref = {b: np.mean([r[b] for r in out["ref"]]) for b in BANDS}
    line = f"  {tag:14s}"
    res = {}
    for name in ("fixed", "own"):
        dF = {b: float(np.mean([r[b] for r in out[name]]) / ref[b] - 1) for b in BANDS}
        worst = max(abs(v) for v in dF.values())
        res[name] = {"dF": dF, "worst": worst}
        line += f"  {name}: worst {100*worst:5.2f}% (band3800 {100*dF['band3800']:+5.2f}%, blue {100*dF['blue']:+5.2f}%)"
    print(line, flush=True)
    return {"ref_bands": ref, **res}


def main(which):
    fixed = RedistributionKernel.load(FIXED)
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz"); n_ion = float(d["n_ion"])
    results = {"which": which, "ng": NG, "runs": {}}
    if which == "src":
        atom = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
        lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
        for t_src in (4000.0, 5000.0, 6000.0, 8000.0):
            out, own = run_config(atom, lo, hi, t_src, {"fixed": fixed})
            s = summarize(f"T_src={t_src:.0f}", out)
            s["kernel_row_l1"] = float(np.abs(own.R - fixed.R).sum(axis=1).mean())
            results["runs"][f"{t_src:.0f}"] = s
    else:
        for t_gas in (2500.0, 3000.0, 4000.0, 5000.0):
            atom = ForestAtom.from_gsi(LEV, TR, t_gas, n_ion, T_EXP, tau_min=1e-3)
            lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
            out, own = run_config(atom, lo, hi, 6000.0, {"fixed": fixed})
            s = summarize(f"T_gas={t_gas:.0f}", out)
            s["n_opacity"] = int(atom.n_opacity); s["tau_max"] = float(atom.op_tau.max())
            # R-difference on the shared row support (edges differ per config;
            # compare via the fixed kernel's edges by rebuilding own on them)
            results["runs"][f"{t_gas:.0f}"] = s
            print(f"      (opacity lines {atom.n_opacity}, tau_max {atom.op_tau.max():.1f})", flush=True)
    (HERE / f"tsweep_{which}.json").write_text(json.dumps(results, indent=1))
    print(f"wrote tsweep_{which}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--which", choices=["src", "gas"], required=True)
    main(ap.parse_args().which)
