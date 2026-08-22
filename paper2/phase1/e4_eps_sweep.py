"""E4/E5: the scalar thermalisation parameter sweep.

The decisive Paper II question: is there one eps in [0, 1] for which the
two-level-atom closure F(eps) reproduces direct A-branching fluorescence
across the spectrum? Phase 1 showed eps = 1 does not. This sweeps
eps in {0, .1, .2, .3, .5, .7, .9, 1} on BOTH legs -- expansion + TLA (the
literature object, SEDONA's opacity_epsilon) and Sobolev + TLA (same opacity
as the branch leg, so only the redistribution closure differs) -- against
sobolev_branch and expansion_branch (E8), on the full La II atom with the
Planck 6000 K launch; 3 seeds; energy-weighted bands.

Per band b: eps_best^(b) = argmin |F_b(eps) - F_b^branch|, by linear
interpolation between grid values, with the 3-seed scatter; and whether
F_b^branch lies outside [min_eps F_b, max_eps F_b] at all ("no eps reproduces
it"). E5: chi^2(eps) over ~200 log bins, energy-weighted, sigma^2 = sigma_TLA^2
+ sigma_branch^2, per band; it ranks eps and locates the residual, it is not a
goodness-of-fit.

Bands (energy-weighted F_b = escaped/launched energy): UV 1142-3300 (the pump
reservoir), blue 3300-4500 (contains 3800-3955), optical 4500-6000, red
6000-9000, NIR 9000-17,697 (transparent: sum(1-e^-tau) < 1, a null control
excluded from the same-eps test). Edges are checked against strong lines
within +-20 A and printed.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C
from forest_mc import ForestAtom, band_ratio, run_mc, spectrum
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of

EPS = (0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0)
BANDS = {"UV": (1142.0, 3300.0), "blue": (3300.0, 4500.0), "optical": (4500.0, 6000.0),
         "red": (6000.0, 9000.0), "NIR": (9000.0, 17697.0), "band3800": (3800.0, 3955.0)}


def check_edges(atom):
    lam = C / atom.op_nu * 1e8
    for name, (lo, hi) in BANDS.items():
        for e in (lo, hi):
            near = np.abs(lam - e) < 20
            strong = near & (atom.op_tau > 1.0)
            if strong.any():
                print(f"  edge {e:.0f} A of {name}: strong lines within 20 A: {np.round(lam[strong],1)} tau {np.round(atom.op_tau[strong],2)}")


def main(n, seeds, tag):
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz"); n_ion = float(d["n_ion"])
    atom = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    print(f"La II full atom, {atom.n_opacity} opacity lines, Planck {T_CORE:.0f} K launch, {n:.0e} packets, seeds {seeds}")
    check_edges(atom)
    edges = np.geomspace(lo, hi, 201)

    def measure(mode, seed, eps=1.0):
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=seed, t_core=T_CORE, eps=eps)
        bands = {b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()}
        sp, se = spectrum(res, edges, weight="energy")
        acc = res["accounting"]
        return bands, sp, se, dict(E_dep_cm=acc["E_dep_cm"] / acc["E_inj"], E_esc=acc["E_esc"] / acc["E_inj"])

    out = {"eps": EPS, "bands": BANDS, "n": n, "seeds": list(seeds), "legs": {}}
    refs = {}
    for ref in ("sobolev_branch", "expansion_branch"):
        rows = [measure(ref, s) for s in seeds]
        refs[ref] = rows
        out["legs"][ref] = {"bands": {b: [r[0][b] for r in rows] for b in BANDS},
                            "spec": np.mean([r[1] for r in rows], axis=0).tolist(),
                            "spec_err": (np.mean([r[2] for r in rows], axis=0) / np.sqrt(len(seeds))).tolist(),
                            "acc": [r[3] for r in rows]}
        print(f"  {ref:18s} " + " ".join(f"{b}={np.mean([r[0][b] for r in rows]):.4f}" for b in BANDS))
    for leg in ("sobolev_tla", "expansion_tla"):
        for eps in EPS:
            t0 = time.time()
            rows = [measure(leg, s, eps) for s in seeds]
            key = f"{leg}_eps{eps:g}"
            out["legs"][key] = {"eps": eps, "bands": {b: [r[0][b] for r in rows] for b in BANDS},
                                "spec": np.mean([r[1] for r in rows], axis=0).tolist(),
                                "spec_err": (np.mean([r[2] for r in rows], axis=0) / np.sqrt(len(seeds))).tolist(),
                                "acc": [r[3] for r in rows]}
            print(f"  {key:22s} " + " ".join(f"{b}={np.mean([r[0][b] for r in rows]):.4f}" for b in BANDS) + f"  [{time.time()-t0:.0f}s]")

    # ---- eps_best per band, per leg, vs sobolev_branch
    ref = out["legs"]["sobolev_branch"]["bands"]
    summary = {}
    for leg in ("sobolev_tla", "expansion_tla"):
        summary[leg] = {}
        for b in BANDS:
            F = np.array([np.mean(out["legs"][f"{leg}_eps{e:g}"]["bands"][b]) for e in EPS])
            Fs = np.array([np.std(out["legs"][f"{leg}_eps{e:g}"]["bands"][b], ddof=1) for e in EPS])
            target = np.mean(ref[b]); tstd = np.std(ref[b], ddof=1)
            # fine interpolation of |F(eps) - target|
            eg = np.linspace(0, 1, 1001); Fi = np.interp(eg, EPS, F)
            i = int(np.argmin(np.abs(Fi - target)))
            inside = (target >= F.min() - 2 * Fs.max() - 2 * tstd) and (target <= F.max() + 2 * Fs.max() + 2 * tstd)
            summary[leg][b] = dict(eps_best=float(eg[i]), residual=float(Fi[i] - target), target=float(target),
                                   F_min=float(F.min()), F_max=float(F.max()), reachable=bool(inside))
    out["summary"] = summary

    # ---- E5: chi^2(eps) over the spectrum, per band contributions
    sref = np.array(out["legs"]["sobolev_branch"]["spec"]); eref = np.array(out["legs"]["sobolev_branch"]["spec_err"])
    lam_mid = C / np.sqrt(edges[1:] * edges[:-1]) * 1e8
    chi = {}
    for leg in ("sobolev_tla", "expansion_tla"):
        chi[leg] = {}
        for e in EPS:
            s_ = np.array(out["legs"][f"{leg}_eps{e:g}"]["spec"]); se_ = np.array(out["legs"][f"{leg}_eps{e:g}"]["spec_err"])
            ok = np.isfinite(s_) & np.isfinite(sref) & (eref > 0)
            z2 = (s_ - sref)**2 / (se_**2 + eref**2)
            per_band = {b: float(np.nansum(z2[ok & (lam_mid >= w[0]) & (lam_mid < w[1])])) for b, w in BANDS.items() if b != "band3800"}
            chi[leg][f"{e:g}"] = dict(total=float(np.nansum(z2[ok])), dof=int(ok.sum()), per_band=per_band)
    out["chi2"] = chi

    print("\n=== eps_best per band (vs sobolev_branch), and whether the branch value is reachable at all ===")
    for leg in ("sobolev_tla", "expansion_tla"):
        print(f"  {leg}:")
        for b in BANDS:
            v = summary[leg][b]
            print(f"    {b:9s} target {v['target']:.4f}  range over eps [{v['F_min']:.4f}, {v['F_max']:.4f}]  eps_best {v['eps_best']:.2f}  residual {v['residual']:+.4f}  {'reachable' if v['reachable'] else 'NOT reachable'}")
    print("\n=== chi^2 / dof over the spectrum ===")
    for leg in ("sobolev_tla", "expansion_tla"):
        print(f"  {leg}: " + "  ".join(f"eps={e:g}: {chi[leg][f'{e:g}']['total']/chi[leg][f'{e:g}']['dof']:.1f}" for e in EPS))
    (HERE / f"e4_eps_sweep{tag}.json").write_text(json.dumps(out, indent=1))
    np.savez(HERE / f"e4_spectra{tag}.npz", edges=edges, **{k: np.array(v["spec"]) for k, v in out["legs"].items()})
    print(f"wrote e4_eps_sweep{tag}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e6); ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3]); ap.add_argument("--tag", default="")
    a = ap.parse_args(); main(int(a.n), tuple(a.seeds), a.tag)
