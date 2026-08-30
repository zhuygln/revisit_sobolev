"""Paper III Phase 0: freeze the reference problem (plan section 3).

Wraps the Paper II instrument unchanged: Sobolev + A*beta branching on the
frozen configuration, 3 seeds x 2e6, with the per-event (nu_abs -> nu_exit)
log the kernel trains on. Gate 0: bands must reproduce the Paper II values
(e4_eps_sweep.json for La II, e9_ceII.json for Ce II) within 3 sigma.

reference_events.npz (the raw event pairs, ~tens of MB) is gitignored --
it regenerates deterministically from the recorded seeds; the json and the
60-bin end-to-end matrix are committed.

Nd II (R5) has no Paper II counterpart, so its Gate 0 is skipped -- see the
note at the check itself.

Usage: python reference.py [--ion laII|ceII|ndII]
"""
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3"):
    sys.path.insert(0, str(p))
from sobolev.constants import C, H
from sobolev.optical_depth import SIGMA_CLASSICAL
from sobolev.atomic_data import load_gsi
from sobolev.populations import boltzmann_fractions_from_levels, statistical_weight
from forest_mc import ForestAtom, band_ratio, run_mc, spectrum
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of
from e9_ceII import CE_LEV, CE_TR, ce_n_ion, WINDOW, TAU_MAX_TARGET

ND_LEV = ROOT / "data/60NdII_levels_calib.txt"
ND_TR = ROOT / "data/60NdII_transitions_calib.txt"

BANDS = {"UV": (1142.0, 3300.0), "blue": (3300.0, 4500.0), "optical": (4500.0, 6000.0),
         "red": (6000.0, 9000.0), "NIR": (9000.0, 17697.0), "band3800": (3800.0, 3955.0)}
SEEDS = (1, 2, 3); N = 2_000_000


def nd_n_ion():
    """setup.py's recipe verbatim on the Nd II files (as ce_n_ion does for Ce).

    Nd II is the densest forest in the GSI set: 57,916 lines land in the
    3850-3950 A window against Ce II's 2,376, so pinning the strongest line
    at tau = 5 puts n_ion an order of magnitude below Ce's and leaves FEWER
    lines above the tau > 1e-3 opacity cut (4,496 vs 22,960). The recipe is
    held fixed for comparability across the three ions; the regime it selects
    is deliberately the blanketed corner of F7.
    """
    levels = load_gsi(ND_LEV); lines = load_gsi(ND_TR)
    lam = lines["WV_Transition"].to_numpy()
    win = lines[(lam >= WINDOW[0]) & (lam < WINDOW[1])].reset_index(drop=True)
    frac = boltzmann_fractions_from_levels(levels, T_SHELL)
    pop = frac[win["Lower"].to_numpy()]
    g_l = statistical_weight(win["J_Lower"].to_numpy())   # Nd II carries '9/2'-style J strings
    f_lu = 10 ** win["Log(gf)"].to_numpy() / g_l
    tau_per_n = SIGMA_CLASSICAL * f_lu * pop * win["WV_Transition"].to_numpy() * 1e-8 * T_EXP
    return TAU_MAX_TARGET / tau_per_n.max(), len(win)


def ion_inputs(ion):
    """(levels path, transitions path, reference n_ion) for an ion.

    Split out of build_atom so the P5/P6 sweeps, which rebuild the atom at
    other temperatures and densities, can reach the same inputs without
    duplicating each ion's normalization.
    """
    if ion == "laII":
        d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz")
        return LEV, TR, float(d["n_ion"])
    if ion == "ceII":
        n_ion, _ = ce_n_ion(); return CE_LEV, CE_TR, n_ion
    n_ion, _ = nd_n_ion(); return ND_LEV, ND_TR, n_ion


def build_atom(ion):
    lev, tr, n_ion = ion_inputs(ion)
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
    if ion == "ndII":
        # Gate 0 checks this wrapper against the Paper II result for the same
        # ion. There is no Paper II Nd II run, so there is nothing to
        # reproduce: the check is skipped rather than faked. The wrapper is
        # already validated -- it reproduced La II and Ce II bit-for-bit
        # (identical event counts, every band 0.0 sigma), and Nd changes only
        # which files build_atom reads.
        print("  Gate 0 skipped: no Paper II Nd II baseline exists to check against")
        print("GATE 0 N/A")
        return True
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
    ap = argparse.ArgumentParser(); ap.add_argument("--ion", default="laII", choices=["laII", "ceII", "ndII"])
    main(ap.parse_args().ion)
