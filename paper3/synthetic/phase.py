"""Paper III E3: when does a grouped-opacity closure fail?

F30/F31 read the failure as a two-regime story separated by line spacing, but
that reading rests on two lanthanides and a referee can fairly say "two peculiar
atoms". This sweeps synthetic forests whose crowding, saturation, spacing and
redistribution range vary independently (`forest.py`), measures the closure
error on each, and asks whether one dimensionless combination collapses them.

The precedent is F15/§4.18: `experiments/r8_saturation/plot.py` already collapses
54 conditions of the pure-absorption problem onto ln<e^D>_w with D = S - E. That
worked because absorption has one control parameter. The scattering problem
plainly needs at least a second, because F31 showed the dense-forest failure is
fluorescent refill rather than saturation.

Candidate axes, deliberately few:

    E/S                     saturation (F15's own dimensionless ratio)
    N_sat per unit ln(lam)  crowding
    <|d ln lam|> / w_group  redistribution range in GROUP widths
    same-group fraction     how often a photon returns to its own group

Each condition: reference `sobolev_branch`, closure `binned_group` with and
without one remembered line, all on the same forest, same kernel, same seeds.

Usage: python phase.py [--n 200000] [--groups 32] [--quick]
"""
import argparse, itertools, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", HERE):
    sys.path.insert(0, str(p))
from forest import synthetic_forest, synthetic_ladder
from forest_mc import band_ratio, run_mc
from redistribution import RedistributionKernel

from sobolev.constants import C
from sobolev.forest_stats import forest_summary

R_CORE, R_OUT, T_EXP, T_CORE = 8.64e13, 2.592e14, 86400.0, 6000.0
SEEDS = (1, 2, 3)

# Bands are defined relative to the forest, not to a fixed astronomical window:
# a synthetic forest has no "optical". Five equal slices in ln lambda plus the
# forest core, which is the analogue of the 3800-3955 A band that carries every
# real-atom failure.
def bands_of(lo, hi):
    e = np.geomspace(lo, hi, 6)
    b = {f"s{i}": (e[i], e[i + 1]) for i in range(5)}
    mid = np.sqrt(lo * hi)
    b["core"] = (mid / 1.02, mid * 1.02)
    return b


def measure(atom, lo, hi, n, mode, kernel=None, memory=False):
    bands = bands_of(lo, hi)
    rows, ev_in, ev_out, inter = [], [], [], []
    for s in SEEDS:
        kw = {}
        if kernel is not None:
            kw["kernel"] = kernel
        if memory:
            kw["line_memory"] = True
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode,
                     seed=s, t_core=T_CORE, collect_events=(kernel is None), **kw)
        rows.append({b: band_ratio(res, w[0], w[1], weight="energy")[0]
                     for b, w in bands.items()})
        inter.append(float(res["n_events"].mean()))
        if kernel is None:
            e = res["events"]; ev_in.append(e[0]); ev_out.append(e[1])
    out = {b: float(np.mean([r[b] for r in rows])) for b in bands}
    ev = (np.concatenate(ev_in), np.concatenate(ev_out)) if kernel is None else None
    return out, ev, float(np.mean(inter))


def condition(params, n, ng):
    """One forest: reference, closure, closure+memory, and every statistic."""
    atom, info = synthetic_forest(**params)
    lo, hi = atom.op_nu.min() * 0.99, atom.op_nu.max() * 1.01
    ref, ev, ev_ref = measure(atom, lo, hi, n, "sobolev_branch")
    if ev is None or ev[0].size < 1000:
        return None
    kern = RedistributionKernel.from_branching_mc(ev[0], ev[1], np.ones(ev[0].size),
                                                  ng, nu_lo=lo, nu_hi=hi)
    row = {"params": info, "ref_bands": ref, "events_per_packet_ref": ev_ref}
    row.update(forest_summary(atom, events=(ev[0], ev[1]), edges=kern.edges))
    # group width in ln lambda -- the scale the closure resolves
    w_group = float(np.log(hi / lo) / ng)
    row["w_group"] = w_group
    row["range_in_groups"] = row["mean_abs_dlnlam"] / w_group

    for tag, mem in (("group", False), ("group_mem", True)):
        b, _, ipp = measure(atom, lo, hi, n, "binned_group", kernel=kern, memory=mem)
        dF = {k: b[k] / ref[k] - 1 for k in ref}
        row[tag] = {"dF": dF, "worst": max(abs(v) for v in dF.values()),
                    "core": dF["core"], "events_per_packet": ipp}
    # the redistribution-only control: same kernel, exact line opacity
    b, _, ipp = measure(atom, lo, hi, n, "sobolev_group", kernel=kern)
    dF = {k: b[k] / ref[k] - 1 for k in ref}
    row["redistribution_only"] = {"worst": max(abs(v) for v in dF.values()),
                                  "core": dF["core"], "events_per_packet": ipp}
    return row


def grid(quick):
    if quick:
        return dict(n_lines=[50, 200], tau=[1.0, 8.0], dlnlam=[0.01, 0.2],
                    f_return=[0.5])
    return dict(n_lines=[25, 50, 100, 200, 400],
                tau=[0.5, 2.0, 8.0, 30.0],
                dlnlam=[0.004, 0.015, 0.06, 0.25],
                f_return=[0.2, 0.8])


def main(n, ng, quick, out_name):
    g = grid(quick)
    keys = list(g)
    combos = list(itertools.product(*(g[k] for k in keys)))
    print(f"{len(combos)} conditions, {n} packets x {len(SEEDS)} seeds, N_g = {ng}", flush=True)
    rows, t0 = [], time.time()
    for i, vals in enumerate(combos):
        p = dict(zip(keys, vals))
        p.update(span=0.3, n_exit=2, tau_spread=0.4, jitter=0.5, seed=7)
        r = condition(p, n, ng)
        if r is None:
            print(f"  [{i+1}/{len(combos)}] {p} -- too few events, skipped", flush=True)
            continue
        rows.append(r)
        print(f"  [{i+1}/{len(combos)}] N={p['n_lines']:3d} tau={p['tau']:4.1f} "
              f"dlnlam={p['dlnlam']:.3f} fret={p['f_return']:.1f} | "
              f"E/S {r['E_over_S']:.3f} Nsat {r['n_sat_per_lnlam']:6.0f} "
              f"rng/grp {r['range_in_groups']:6.2f} same {r['same_group_frac']:.2f} | "
              f"Rij {100*r['redistribution_only']['worst']:6.2f}% "
              f"grp {100*r['group']['worst']:7.2f}% "
              f"+mem {100*r['group_mem']['worst']:7.2f}%", flush=True)
    (HERE / out_name).write_text(json.dumps({"n": n, "ng": ng, "rows": rows}, indent=1))
    print(f"wrote {out_name}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e5)
    ap.add_argument("--groups", type=int, default=32)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    main(int(a.n), a.groups, a.quick,
         a.out or ("phase_quick.json" if a.quick else "phase.json"))
