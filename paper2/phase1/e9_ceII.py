"""E9: the same minimum suite on Ce II -- is outcome B ion-dependent (outcome C)?

The Ce II atom is normalized by the SAME recipe as the La II reference
(setup.py): n_ion such that the strongest classical Sobolev depth in
3850-3950 A at T = 3000 K equals 5. Everything else identical: slow shell,
6000 K Planck photon launch over the ion's own opacity extent, classical
transport (the Phase-2 reference convention), energy-weighted bands,
3 seeds x 2e6. Legs: both pure-absorption controls, both branching legs
(physics and E8), and the TLA sweep on both opacities. The test:
eps_best^La ?= eps_best^Ce per band, and whether reachability changes.
La II values are read from e4_eps_sweep.json.

Usage: nohup <venv python> -u e9_ceII.py > e9.log 2>&1 &
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C
from sobolev.optical_depth import SIGMA_CLASSICAL
from forest_mc import ForestAtom, band_ratio, run_mc, spectrum
from run_forest import R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of
from sobolev.atomic_data import load_gsi
from sobolev.populations import boltzmann_fractions_from_levels, statistical_weight

CE_LEV = ROOT / "data/58CeII_levels_calib.txt"
CE_TR = ROOT / "data/58CeII_transitions_calib.txt"
WINDOW = (3850.0, 3950.0); TAU_MAX_TARGET = 5.0
EPS = (0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
BANDS = {"UV": (1142.0, 3300.0), "blue": (3300.0, 4500.0), "optical": (4500.0, 6000.0),
         "red": (6000.0, 9000.0), "NIR": (9000.0, 17697.0), "band3800": (3800.0, 3955.0)}


def ce_n_ion():
    """setup.py's recipe verbatim, on the Ce II files."""
    levels = load_gsi(CE_LEV); lines = load_gsi(CE_TR)
    lam = lines["WV_Transition"].to_numpy()
    win = lines[(lam >= WINDOW[0]) & (lam < WINDOW[1])].reset_index(drop=True)
    frac = boltzmann_fractions_from_levels(levels, T_SHELL)
    pop = frac[win["Lower"].to_numpy()]
    g_l = statistical_weight(win["J_Lower"].to_numpy())   # Ce II carries '7/2'-style J strings
    f_lu = 10 ** win["Log(gf)"].to_numpy() / g_l
    tau_per_n = SIGMA_CLASSICAL * f_lu * pop * win["WV_Transition"].to_numpy() * 1e-8 * T_EXP
    return TAU_MAX_TARGET / tau_per_n.max(), len(win)


def main(n, seeds):
    n_ion, n_win = ce_n_ion()
    print(f"Ce II: {n_win} lines in {WINDOW[0]:.0f}-{WINDOW[1]:.0f} A, n_ion = {n_ion:.2f} cm^-3")
    t0 = time.time()
    atom = ForestAtom.from_gsi(CE_LEV, CE_TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    lam_op = C / atom.op_nu * 1e8
    print(f"atom built in {time.time()-t0:.0f}s: {atom.n_lines_total} lines total, "
          f"{atom.n_opacity} opacity lines {lam_op.min():.0f}-{lam_op.max():.0f} A, "
          f"tau_max {atom.op_tau.max():.1f}, N(tau>1) {(atom.op_tau>1).sum()}, "
          f"N(tau>0.1) {(atom.op_tau>0.1).sum()}", flush=True)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    edges = np.geomspace(lo, hi, 201)

    def measure(mode, seed, eps=1.0):
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=seed, t_core=T_CORE, eps=eps)
        bands = {b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()}
        sp, se = spectrum(res, edges, weight="energy")
        return bands, sp, se

    out = {"n": n, "seeds": list(seeds), "bands": BANDS, "eps": EPS, "n_ion": float(n_ion),
           "n_opacity": int(atom.n_opacity), "legs": {}}
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
            out["legs"][key] = {"bands": {b: [r[0][b] for r in rows] for b in BANDS},
                                "spec": np.mean([r[1] for r in rows], axis=0).tolist(),
                                "spec_err": (np.mean([r[2] for r in rows], axis=0) / np.sqrt(len(seeds))).tolist()}
            print(f"  {key:24s} " + " ".join(f"{b}={np.mean([r[0][b] for r in rows]):.4f}" for b in BANDS)
                  + f"  [{time.time()-t0:.0f}s]", flush=True)

    # eps_best per band vs Ce sobolev_branch; then the La comparison
    ref = {b: np.mean(out["legs"]["sobolev_branch"]["bands"][b]) for b in BANDS}
    out["summary"] = {}
    for leg in ("sobolev_tla", "expansion_tla"):
        out["summary"][leg] = {}
        for b in BANDS:
            F = np.array([np.mean(out["legs"][f"{leg}_eps{e:g}"]["bands"][b]) for e in EPS])
            eg = np.linspace(0, 1, 1001); Fi = np.interp(eg, EPS, F)
            i = int(np.argmin(np.abs(Fi - ref[b])))
            out["summary"][leg][b] = dict(eps_best=float(eg[i]), residual=float(Fi[i] - ref[b]),
                                          target=float(ref[b]), F_min=float(F.min()), F_max=float(F.max()),
                                          reachable=bool(F.min() - 0.01 <= ref[b] <= F.max() + 0.01))
    la = json.load(open(HERE / "e4_eps_sweep.json"))["summary"]
    print("\n=== eps_best per band: La II (E4) vs Ce II (E9) ===")
    for leg in ("sobolev_tla", "expansion_tla"):
        print(f"  {leg}:")
        for b in BANDS:
            lv, cv = la[leg][b], out["summary"][leg][b]
            def show(v):
                return f"{v['eps_best']:.2f}" if v["reachable"] else "unreach"
            print(f"    {b:9s} La {show(lv)} (target {lv['target']:.4f})   "
                  f"Ce {show(cv)} (target {cv['target']:.4f})   d_eps={cv['eps_best']-lv['eps_best']:+.2f}")
    # chi2 over the Ce spectrum
    sref = np.array(out["legs"]["sobolev_branch"]["spec"]); eref = np.array(out["legs"]["sobolev_branch"]["spec_err"])
    chi = {}
    for leg in ("sobolev_tla", "expansion_tla"):
        chi[leg] = {}
        for e in EPS:
            s_ = np.array(out["legs"][f"{leg}_eps{e:g}"]["spec"]); se_ = np.array(out["legs"][f"{leg}_eps{e:g}"]["spec_err"])
            ok = np.isfinite(s_) & np.isfinite(sref) & (eref > 0)
            z2 = (s_ - sref) ** 2 / (se_ ** 2 + eref ** 2)
            chi[leg][f"{e:g}"] = dict(total=float(np.nansum(z2[ok])), dof=int(ok.sum()))
    out["chi2"] = chi
    print("\n=== chi^2/dof vs Ce sobolev_branch ===")
    for leg in ("sobolev_tla", "expansion_tla"):
        print(f"  {leg}: " + "  ".join(f"eps={e:g}: {chi[leg][f'{e:g}']['total']/chi[leg][f'{e:g}']['dof']:.1f}" for e in EPS))
    (HERE / "e9_ceII.json").write_text(json.dumps(out, indent=1))
    print("wrote e9_ceII.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e6); ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    a = ap.parse_args(); main(int(a.n), tuple(a.seeds))
