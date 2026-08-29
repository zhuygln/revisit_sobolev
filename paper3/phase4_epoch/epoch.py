"""Paper III P6: epoch dependence, and the tau-collapse test.

t = 0.5/1/2/4 d, homologous: fixed velocity shell (r = v t), n = n_ref
(t/1d)^-3, so tau propto n t = tau_ref (t/1d)^-2. Three kernels per epoch:
"fixed" (the 1 d kernel), "own" (trained at that epoch), and "tau_matched"
(trained at 1 d geometry but with n_ion scaled by (t/1d)^-2, i.e. the SAME
tau set as the target epoch). Prediction, in advance: tau_matched == own to
within noise -- beta and the branch chains depend only on {tau_j}, and
geometry never enters the kernel -- so R(t) collapses onto R(tau_scale).
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase0_reference"):
    sys.path.insert(0, str(p))
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of
from redistribution import RedistributionKernel
from reference import BANDS, SEEDS, N

NG = 32
FIXED = ROOT / "paper3/phase1_groups/kernel_laII_ng32.npz"


def train(atom, r_core, r_out, t_exp, lo, hi):
    ev_in, ev_out = [], []
    for s in SEEDS:
        res = run_mc(atom, r_core, r_out, t_exp, lo, hi, N, "sobolev_branch",
                     seed=s, t_core=T_CORE, collect_events=True)
        e = res["events"]; ev_in.append(e[0]); ev_out.append(e[1])
    nu_in, nu_out = np.concatenate(ev_in), np.concatenate(ev_out)
    return RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size), NG)


def main():
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz"); n_ref = float(d["n_ion"])
    fixed = RedistributionKernel.load(FIXED)
    results = {"ng": NG, "runs": {}}
    for t_d in (0.5, 1.0, 2.0, 4.0):
        t_exp = T_EXP * t_d; rc, ro = R_CORE * t_d, R_OUT * t_d
        atom = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ref * t_d ** -3, t_exp, tau_min=1e-3)
        lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
        # tau-matched auxiliary: 1 d geometry, n scaled to the target tau set
        aux = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ref * t_d ** -2, T_EXP, tau_min=1e-3)
        alo, ahi = aux.op_nu.min() * 0.995, aux.op_nu.max() * 1.005
        k_tau = train(aux, R_CORE, R_OUT, T_EXP, alo, ahi)
        # reference + own kernel at the target epoch
        ref, ev_in, ev_out = [], [], []
        for s in SEEDS:
            res = run_mc(atom, rc, ro, t_exp, lo, hi, N, "sobolev_branch",
                         seed=s, t_core=T_CORE, collect_events=True)
            ref.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
            e = res["events"]; ev_in.append(e[0]); ev_out.append(e[1])
        nu_in, nu_out = np.concatenate(ev_in), np.concatenate(ev_out)
        k_own = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size), NG)
        refm = {b: np.mean([r[b] for r in ref]) for b in BANDS}
        row = {"ref_bands": refm, "tau_max": float(atom.op_tau.max()), "n_opacity": int(atom.n_opacity)}
        line = f"  t={t_d:3.1f}d (tau_max {atom.op_tau.max():5.1f})"
        for name, kern in (("fixed", fixed), ("tau_matched", k_tau), ("own", k_own)):
            rows = []
            for s in SEEDS:
                res = run_mc(atom, rc, ro, t_exp, lo, hi, N, "sobolev_group",
                             seed=s, t_core=T_CORE, kernel=kern)
                rows.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
            dF = {b: float(np.mean([r[b] for r in rows]) / refm[b] - 1) for b in BANDS}
            row[name] = {"dF": dF, "worst": max(abs(v) for v in dF.values())}
            line += f"  {name}: {100*row[name]['worst']:5.2f}%"
        print(line, flush=True)
        results["runs"][f"{t_d:g}"] = row
    (HERE / "epoch.json").write_text(json.dumps(results, indent=1))
    print("wrote epoch.json")


if __name__ == "__main__":
    main()
