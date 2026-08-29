"""Paper III Phase 0: freeze the reference problem (plan section 3).

Wraps the Paper II instrument unchanged: Sobolev + A*beta branching on the
frozen configuration, 3 seeds x 2e6, with the per-event (nu_abs -> nu_exit)
log the kernel trains on. Gate 0: bands must reproduce the Paper II values
(e4_eps_sweep.json for La II, e9_ceII.json for Ce II) within 3 sigma.

reference_events.npz (the raw event pairs, ~tens of MB) is gitignored --
it regenerates deterministically from the recorded seeds; the json and the
60-bin end-to-end matrix are committed.

Usage: python reference.py [--ion laII|ceII]
"""
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3"):
    sys.path.insert(0, str(p))
from sobolev.constants import C, H
from forest_mc import ForestAtom, band_ratio, run_mc, spectrum
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of
from e9_ceII import CE_LEV, CE_TR, ce_n_ion

BANDS = {"UV": (1142.0, 3300.0), "blue": (3300.0, 4500.0), "optical": (4500.0, 6000.0),
         "red": (6000.0, 9000.0), "NIR": (9000.0, 17697.0), "band3800": (3800.0, 3955.0)}
SEEDS = (1, 2, 3); N = 2_000_000


def build_atom(ion):
    if ion == "laII":
        d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz")
        n_ion, lev, tr = float(d["n_ion"]), LEV, TR
    else:
        n_ion, _ = ce_n_ion(); lev, tr = CE_LEV, CE_TR
    atom = ForestAtom.from_gsi(lev, tr, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in (lev, tr)}
    return atom, n_ion, hashes


def end_to_end_matrix(res, edges):
    esc = res["fate"] == 1
    M, _, _ = np.histogram2d(res["nu_launch"][esc], res["nu_out_all"][esc], bins=[edges, edges],
                             weights=H * res["nu_out_all"][esc])
    row, _ = np.histogram(res["nu_launch"], edges, weights=H * res["nu_launch"])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(row[:, None] > 0, M / row[:, None], 0.0)


def main(ion):
    atom, n_ion, hashes = build_atom(ion)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    sp_edges = np.geomspace(lo, hi, 201); m_edges = np.geomspace(lo, hi, 61)
    bands, specs, errs, acc, ev = [], [], [], [], []
    Msum = 0; t0 = time.time()
    for s in SEEDS:
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, N, "sobolev_branch",
                     seed=s, t_core=T_CORE, collect_events=True)
        bands.append({b: band_ratio(res, *nu_of(*w), weight="energy")[0] for b, w in BANDS.items()})
        sp, se = spectrum(res, sp_edges, weight="energy"); specs.append(sp); errs.append(se)
        a = res["accounting"]
        acc.append(dict(bol=a["E_esc"] / a["E_inj"], core=a["E_core"] / a["E_inj"],
                        dep_cm=a["E_dep_cm"] / a["E_inj"], events_per_packet=float(res["n_events"].mean())))
        ev.append(res["events"]); Msum = Msum + end_to_end_matrix(res, m_edges)
    wall = time.time() - t0
    # float64: exits are exact line rest frequencies, and the kernel's
    # discrete within-group tables must reproduce them bit-for-bit so a
    # re-emitted packet does not re-sweep its own line (see notebook 9t)
    nu_in = np.concatenate([e[0] for e in ev])
    nu_out = np.concatenate([e[1] for e in ev])
    out = dict(ion=ion, n_ion=n_ion, hashes=hashes, seeds=list(SEEDS), n_packets=N,
               t_gas=T_SHELL, t_src=T_CORE, r_core=R_CORE, r_out=R_OUT, t_exp=T_EXP,
               nu_lo=lo, nu_hi=hi, n_opacity=int(atom.n_opacity), wall_s=wall,
               n_events=int(nu_in.size),
               bands={b: [r[b] for r in bands] for b in BANDS},
               accounting=acc)
    (HERE / f"reference_{ion}.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(HERE / f"reference_spectrum_{ion}.npz", edges=sp_edges,
                        spec=np.mean(specs, 0), err=np.mean(errs, 0) / np.sqrt(len(SEEDS)))
    np.savez_compressed(HERE / f"reference_matrix_{ion}.npz", edges=m_edges, M=Msum / len(SEEDS))
    np.savez_compressed(HERE / f"reference_events_{ion}.npz", nu_in=nu_in, nu_out=nu_out)
    print(f"{ion}: {nu_in.size} events, {wall:.0f}s wall; bands " +
          " ".join(f"{b}={np.mean(out['bands'][b]):.4f}" for b in BANDS))

    # ---- Gate 0
    ref_file = HERE.parents[0] / ("phase1_groups" if False else ".")  # keep flat
    src = ROOT / ("paper2/phase1/e4_eps_sweep.json" if ion == "laII" else "paper2/phase1/e9_ceII.json")
    prev = json.load(open(src))["legs"]["sobolev_branch"]["bands"] if ion == "laII" else \
           json.load(open(src))["legs"]["sobolev_branch"]["bands"]
    ok = True
    for b in BANDS:
        m_new, s_new = np.mean(out["bands"][b]), np.std(out["bands"][b], ddof=1)
        m_old, s_old = np.mean(prev[b]), np.std(prev[b], ddof=1)
        sig = max(np.hypot(s_new, s_old), 1e-4)
        pull = abs(m_new - m_old) / sig
        flag = "ok" if pull < 3 else "FAIL"
        ok &= pull < 3
        print(f"  Gate 0 {b:9s}: new {m_new:.4f} vs paper2 {m_old:.4f}  ({pull:.1f} sigma) {flag}")
    print("GATE 0 " + ("PASSED" if ok else "FAILED -- do not proceed"))
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ion", default="laII", choices=["laII", "ceII"])
    main(ap.parse_args().ion)
