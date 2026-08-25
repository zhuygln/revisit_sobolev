"""E10: the La II + Ce II mixture -- does dense blanketing suppress the
redistribution error, move it, or change eps_best?

Each ion at its own reference density: La II's n_ion from forest_lines.npz
(tau_max = 5 in 3850-3950 A), Ce II's from the same recipe (e9_ceII.ce_n_ion).
The blend is the union of the two reference forests -- denser than either,
with ion-internal branching (ions share only the radiation field). Slow
shell, classical transport, Planck 6000 K photon launch over the union
opacity extent, energy-weighted bands, 3 seeds x 2e6. Legs: both absorption
controls (the Paper I F7 blanketing connection), both branching legs, and
the TLA sweep on both opacities. eps_best per band vs the blend's own
sobolev_branch, compared against La-only (e4_eps_sweep.json) and Ce-only
(e9_ceII.json).

Usage: nohup <venv python> -u e10_blend.py > e10.log 2>&1 &
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of
from e9_ceII import CE_LEV, CE_TR, ce_n_ion, BANDS, EPS


def main(n, seeds):
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz")
    n_la = float(d["n_ion"]); n_ce, _ = ce_n_ion()
    print(f"blend: La II n_ion = {n_la:.2f}, Ce II n_ion = {n_ce:.2f} (each its reference)")
    t0 = time.time()
    atom = ForestAtom.from_gsi_blend([(LEV, TR, n_la), (CE_LEV, CE_TR, n_ce)],
                                     T_SHELL, T_EXP, tau_min=1e-3)
    lam_op = C / atom.op_nu * 1e8
    print(f"atom built in {time.time()-t0:.0f}s: {atom.n_lines_total} lines, "
          f"{atom.n_opacity} opacity lines {lam_op.min():.0f}-{lam_op.max():.0f} A, "
          f"tau_max {atom.op_tau.max():.1f}, N(tau>1) {(atom.op_tau>1).sum()}", flush=True)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005

    def measure(mode, seed, eps=1.0):
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=seed, t_core=T_CORE, eps=eps)
        return {b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()}

    out = {"n": n, "seeds": list(seeds), "bands": BANDS, "eps": EPS,
           "n_ion": {"LaII": n_la, "CeII": n_ce}, "n_opacity": int(atom.n_opacity), "legs": {}}
    for tag, mode, eps_list in (("sobolev_absorb", "sobolev_absorb", [None]),
                                ("expansion_absorb", "expansion_absorb", [None]),
                                ("sobolev_branch", "sobolev_branch", [None]),
                                ("expansion_branch", "expansion_branch", [None]),
                                ("sobolev_tla", "sobolev_tla", EPS),
                                ("expansion_tla", "expansion_tla", EPS)):
        for eps in eps_list:
            key = tag if eps is None else f"{tag}_eps{eps:g}"
            t0 = time.time()
            rows = [measure(mode, s, eps if eps is not None else 1.0) for s in seeds]
            out["legs"][key] = {b: [r[b] for r in rows] for b in BANDS}
            print(f"  {key:24s} " + " ".join(f"{b}={np.mean([r[b] for r in rows]):.4f}" for b in BANDS)
                  + f"  [{time.time()-t0:.0f}s]", flush=True)

    ref = {b: np.mean(out["legs"]["sobolev_branch"][b]) for b in BANDS}
    out["summary"] = {}
    for leg in ("sobolev_tla", "expansion_tla"):
        out["summary"][leg] = {}
        for b in BANDS:
            F = np.array([np.mean(out["legs"][f"{leg}_eps{e:g}"][b]) for e in EPS])
            eg = np.linspace(0, 1, 1001); Fi = np.interp(eg, EPS, F)
            i = int(np.argmin(np.abs(Fi - ref[b])))
            out["summary"][leg][b] = dict(eps_best=float(eg[i]), residual=float(Fi[i] - ref[b]),
                                          target=float(ref[b]), F_min=float(F.min()), F_max=float(F.max()),
                                          reachable=bool(F.min() - 0.01 <= ref[b] <= F.max() + 0.01))

    # comparisons: single-ion eps_best (La from E4, Ce from E9) vs blend
    la = json.load(open(HERE / "e4_eps_sweep.json"))["summary"]
    ce = json.load(open(HERE / "e9_ceII.json"))["summary"]
    print("\n=== eps_best per band: La-only vs Ce-only vs blend ===")
    for leg in ("sobolev_tla", "expansion_tla"):
        print(f"  {leg}:")
        for b in BANDS:
            def show(v):
                return f"{v['eps_best']:.2f}" if v["reachable"] else "unreach"
            print(f"    {b:9s} La {show(la[leg][b])}   Ce {show(ce[leg][b])}   blend {show(out['summary'][leg][b])}")
    # F7 connection: blanketing and the two errors in the 3800-3955 band
    print("\n=== blanketing (blend vs single forests), key differentials ===")
    sa, ea = np.mean(out["legs"]["sobolev_absorb"]["band3800"]), np.mean(out["legs"]["expansion_absorb"]["band3800"])
    sb, eb = ref["band3800"], np.mean(out["legs"]["expansion_branch"]["band3800"])
    et1 = np.mean(out["legs"]["expansion_tla_eps1"]["band3800"])
    da = f"{100*(ea/sa-1):+.1f}%" if sa > 0 else "band black under Sobolev"
    print(f"  3800-3955: absorb Sob {sa:.4f} exp {ea:.4f} (Delta_absorb {da})")
    print(f"             branch Sob {sb:.4f} exp {eb:.4f} (opacity error under fluorescence {100*(eb/sb-1):+.1f}%)")
    print(f"             exp TLA eps=1 {et1:.4f} (redistribution error {100*(et1/sb-1):+.1f}%)")
    (HERE / "e10_blend.json").write_text(json.dumps(out, indent=1))
    print("wrote e10_blend.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e6); ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    a = ap.parse_args(); main(int(a.n), tuple(a.seeds))
