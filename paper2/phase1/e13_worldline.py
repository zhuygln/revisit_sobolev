"""E13: does the closure verdict survive v_bulk ~ 0.1c with worldline-
consistent transport?

Part A -- E13.3 single-line control at beta_out = 0.01/0.05/0.10/0.20: the MC
against `sobolev_attenuation` in both conventions, plus the analytic
frozen-first-order-vs-worldline gap (the trap: a frozen snapshot manufactures
an O(beta) artifact; the physical worldline correction is ~beta^2/2).

Part B -- E13.2 matched line strength: tau_S = sigma f n lambda t has no
shell-velocity dependence, so keeping n_ion, T, t_exp from the slow reference
makes the fast shell's tau distribution IDENTICAL by construction; the
distribution is printed for the record. What changes is the sweep:
Delta v_shell = 0.1c crosses ~15x more resonances per unit ln nu.

Part C -- E13.4/13.5 full La II: slow shell (0.0033-0.01c, the Phase-2
reference) and fast shell (0.05-0.15c), worldline transport on both, plus the
frozen-first-order control on the fast shell; sobolev_branch as the physics,
TLA legs bracketing eps_best^slow; energy-weighted bands; eps_best per band by
interpolation; P(lambda_out | lambda_in) for sobolev_branch on both shells.
v_D never enters: the MC's resonances are delta functions (E2).

Usage: nohup <venv python> -u e13_worldline.py [--n 2e6] > e13.log 2>&1 &
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C, H
from sobolev.optical_depth import tau_sobolev
from sobolev.sobolev_leg import sobolev_attenuation
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of

CT = C * T_EXP
BANDS = {"UV": (1142.0, 3300.0), "blue": (3300.0, 4500.0), "optical": (4500.0, 6000.0),
         "red": (6000.0, 9000.0), "NIR": (9000.0, 17697.0), "band3800": (3800.0, 3955.0)}
SHELLS = {"slow": (R_CORE, R_OUT), "fast": (0.05 * CT, 0.15 * CT)}
EPS_SOB = (0.0, 0.065, 0.3, 1.0)
EPS_EXP = (0.0, 0.2, 0.3, 0.5, 1.0)


def part_a():
    print("=== E13.3 single-line control (tau_S = 5) ===")
    nu0, f = 7.6e14, 0.01
    n_ion = 5.0 / tau_sobolev(f, 1.0, C / nu0, T_EXP)
    atom = ForestAtom(nu0=np.array([nu0]), f_osc=np.array([f]), n_lower=np.array([n_ion]),
                      n_upper=np.array([0.0]), A=np.array([1e8]), lower=np.array([0]),
                      upper=np.array([1]), t_exp=T_EXP, stim=False)
    rows = []
    for b_out in (0.01, 0.05, 0.10, 0.20):
        r_core, r_out = b_out / 3 * CT, b_out * CT
        nu_lo, nu_hi = nu0 * 0.98, nu0 * (1 + 1.6 * b_out)
        edges = np.linspace(nu_lo, nu_hi, 121); mid = 0.5 * (edges[1:] + edges[:-1])
        ana = {rel: sobolev_attenuation(mid, [(nu0, f, 1.0)], r_core, r_out, T_EXP, n_ion,
                                        relativity=rel, n_p=400) for rel in (None, "worldline")}
        row = {"beta_out": b_out}
        for rel in (None, "worldline"):
            res = run_mc(atom, r_core, r_out, T_EXP, nu_lo, nu_hi, 2_000_000,
                         "sobolev_absorb", seed=1, relativity=rel)
            Nl, _ = np.histogram(res["nu_launch"], edges)
            Ne, _ = np.histogram(res["nu_launch"][res["fate"] == 1], edges)
            Nc, _ = np.histogram(res["nu_launch"][res["fate"] == 2], edges)
            tr = Ne / np.maximum(Nl - Nc, 1); m = Nl > 2000
            key = "frozen_first" if rel is None else "worldline"
            row[key] = dict(rms=float(np.sqrt(np.mean((tr[m] - ana[rel][m]) ** 2))),
                            mean_mc=float(tr[m].mean()), mean_ana=float(ana[rel][m].mean()))
        g = ana[None] - ana["worldline"]
        row["gap_mean"] = float(g.mean()); row["gap_max"] = float(np.abs(g).max())
        rows.append(row)
        print(f"  beta_out={b_out:4.2f}  worldline MC/ana {row['worldline']['mean_mc']:.4f}/"
              f"{row['worldline']['mean_ana']:.4f} (rms {row['worldline']['rms']:.4f})   "
              f"frozen MC/ana {row['frozen_first']['mean_mc']:.4f}/{row['frozen_first']['mean_ana']:.4f}   "
              f"convention gap mean {row['gap_mean']:+.4f} max {row['gap_max']:.4f}")
    return rows


def part_b(atom):
    t = atom.op_tau
    stats = dict(tau_max=float(t.max()), n_tau_gt1=int((t > 1).sum()),
                 n_tau_gt01=int((t > 0.1).sum()), n_opacity=int(t.size))
    print(f"=== E13.2 line strength (both shells, by construction identical): "
          f"tau_max={stats['tau_max']:.1f}, N(tau>1)={stats['n_tau_gt1']}, "
          f"N(tau>0.1)={stats['n_tau_gt01']}, opacity lines {stats['n_opacity']} ===")
    return stats


def matrix(res, edges):
    esc = res["fate"] == 1
    w_in = H * res["nu_launch"]; w_out = H * res["nu_out_all"]
    Hm, _, _ = np.histogram2d(res["nu_launch"][esc], res["nu_out_all"][esc],
                              bins=[edges, edges], weights=w_out[esc])
    row, _ = np.histogram(res["nu_launch"], edges, weights=w_in)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(row[:, None] > 0, Hm / row[:, None], np.nan)


def main(n, seeds):
    out = {"n": n, "seeds": list(seeds), "bands": BANDS,
           "shells": {k: [v[0] / CT, v[1] / CT] for k, v in SHELLS.items()}}
    out["control"] = part_a()
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz"); n_ion = float(d["n_ion"])
    atom = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    out["tau_stats"] = part_b(atom)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    m_edges = np.geomspace(lo, hi, 61)

    legs = [("sobolev_branch", {}, "worldline")]
    legs += [(f"sobolev_tla_eps{e:g}", {"mode": "sobolev_tla", "eps": e}, "worldline") for e in EPS_SOB]
    legs += [(f"expansion_tla_eps{e:g}", {"mode": "expansion_tla", "eps": e}, "worldline") for e in EPS_EXP]
    frozen_legs = [("sobolev_branch", {}, None),
                   ("expansion_tla_eps0.3", {"mode": "expansion_tla", "eps": 0.3}, None),
                   ("expansion_tla_eps1", {"mode": "expansion_tla", "eps": 1.0}, None)]

    out["legs"] = {}
    for shell, (r_core, r_out) in SHELLS.items():
        todo = legs + (frozen_legs if shell == "fast" else [])
        for tag, kw, rel in todo:
            mode = kw.get("mode", tag); eps = kw.get("eps", 1.0)
            key = f"{shell}_{tag}" + ("_frozen" if rel is None else "")
            t0 = time.time(); rows = []
            for s in seeds:
                res = run_mc(atom, r_core, r_out, T_EXP, lo, hi, n, mode, seed=s,
                             t_core=T_CORE, eps=eps, relativity=rel)
                rows.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0]
                             for b, w in BANDS.items()})
                if tag == "sobolev_branch" and rel == "worldline" and s == seeds[0]:
                    np.savez(HERE / f"e13_matrix_{shell}.npz", edges=m_edges,
                             M=matrix(res, m_edges))
            out["legs"][key] = {b: [r[b] for r in rows] for b in BANDS}
            print(f"  {key:34s} " + " ".join(
                f"{b}={np.mean([r[b] for r in rows]):.4f}" for b in BANDS)
                + f"  [{time.time() - t0:.0f}s]", flush=True)

    # eps_best per band per shell (worldline), vs that shell's own branch leg
    out["summary"] = {}
    for shell in SHELLS:
        ref = {b: np.mean(out["legs"][f"{shell}_sobolev_branch"][b]) for b in BANDS}
        out["summary"][shell] = {}
        for leg, grid in (("sobolev_tla", EPS_SOB), ("expansion_tla", EPS_EXP)):
            res = {}
            for b in BANDS:
                F = np.array([np.mean(out["legs"][f"{shell}_{leg}_eps{e:g}"][b]) for e in grid])
                eg = np.linspace(0, 1, 1001); Fi = np.interp(eg, grid, F)
                i = int(np.argmin(np.abs(Fi - ref[b])))
                res[b] = dict(eps_best=float(eg[i]), residual=float(Fi[i] - ref[b]),
                              target=float(ref[b]), F_min=float(F.min()), F_max=float(F.max()),
                              reachable=bool(F.min() - 0.01 <= ref[b] <= F.max() + 0.01))
            out["summary"][shell][leg] = res

    print("\n=== eps_best, slow vs fast (worldline) ===")
    for leg in ("sobolev_tla", "expansion_tla"):
        print(f"  {leg}:")
        for b in BANDS:
            s_, f_ = out["summary"]["slow"][leg][b], out["summary"]["fast"][leg][b]
            def show(v):
                return f"{v['eps_best']:.2f}" if v["reachable"] else "unreach"
            print(f"    {b:9s} slow {show(s_)} (target {s_['target']:.4f})   "
                  f"fast {show(f_)} (target {f_['target']:.4f})   "
                  f"d_eps={f_['eps_best'] - s_['eps_best']:+.2f}")
    print("\n=== frozen-first-order control on the fast shell ===")
    for tag in ("sobolev_branch", "expansion_tla_eps0.3", "expansion_tla_eps1"):
        w_ = {b: np.mean(out["legs"][f"fast_{tag}"][b]) for b in BANDS}
        f_ = {b: np.mean(out["legs"][f"fast_{tag}_frozen"][b]) for b in BANDS}
        print(f"  {tag:22s} " + " ".join(f"{b}: wl {w_[b]:.4f} fr {f_[b]:.4f}" for b in ("blue", "optical", "band3800")))
    (HERE / "e13_worldline.json").write_text(json.dumps(out, indent=1))
    print("wrote e13_worldline.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e6); ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    a = ap.parse_args(); main(int(a.n), tuple(a.seeds))
