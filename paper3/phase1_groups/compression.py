"""Paper III R1-R3: the compression sweep (plan sections 6-7).

For N_g in {4, 8, 16, 32, 64, 128}: build the kernel from the reference's
event log, validate energy conservation, run the sobolev_group leg on the
identical Sobolev opacity, and score against Sobolev + full branching:
band residuals dF_b (incl. bolometric), E5-style spectral chi^2, end-to-end
matrix row-L1 / energy-weighted error / block flows, and cost (wall time,
table size, interactions per packet). Gate 1 thresholds: strong = |dF_b| <
5% in every major band at N_g <= 32; excellent = < 2%.

Usage: python compression.py [--ion laII|ceII] [--groups 4 8 ...]
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3"):
    sys.path.insert(0, str(p))
from sobolev.constants import C, H
from forest_mc import band_ratio, run_mc, spectrum
from run_forest import R_CORE, R_OUT, T_EXP, T_CORE, nu_of
from redistribution import RedistributionKernel
sys.path.insert(0, str(ROOT / "paper3/phase0_reference"))
from reference import BANDS, SEEDS, N, build_atom, end_to_end_matrix

P0 = ROOT / "paper3/phase0_reference"
BLOCKS = [("blue", "blue"), ("blue", "optical"), ("blue", "red"), ("optical", "red")]


def main(ion, groups):
    ref = json.load(open(P0 / f"reference_{ion}.json"))
    ev = np.load(P0 / f"reference_events_{ion}.npz")
    nu_in, nu_out = ev["nu_in"].astype(float), ev["nu_out"].astype(float)
    w_ev = np.ones(nu_in.size)
    rm = np.load(P0 / f"reference_matrix_{ion}.npz"); M_ref = rm["M"]; m_edges = rm["edges"]
    sp_ref = np.load(P0 / f"reference_spectrum_{ion}.npz")
    atom, _, _ = build_atom(ion)
    lo, hi = ref["nu_lo"], ref["nu_hi"]
    sp_edges = sp_ref["edges"]
    ref_bands = {b: np.mean(v) for b, v in ref["bands"].items()}
    ref_bands_std = {b: np.std(v, ddof=1) for b, v in ref["bands"].items()}
    ref_bol = np.mean([a["bol"] for a in ref["accounting"]])
    lam_mid = C / np.sqrt(m_edges[1:] * m_edges[:-1]) * 1e8

    def band_rows(M, name):
        sel = lambda lo_, hi_: (lam_mid >= lo_) & (lam_mid < hi_)
        return {f"{a}->{b}": float(M[sel(*BANDS[a])][:, sel(*BANDS[b])].sum(axis=1).mean())
                for a, b in BLOCKS}

    out = {"ion": ion, "groups": groups, "reference": {"bands": ref_bands, "bol": ref_bol,
           "blocks": band_rows(M_ref, "ref")}, "runs": {}}
    print(f"{ion}: reference bol {ref_bol:.4f}; blocks {out['reference']['blocks']}")
    for ng in groups:
        kern = RedistributionKernel.from_branching_mc(nu_in, nu_out, w_ev, ng,
                                                      metadata={"ion": ion, "n_groups": ng})
        eres = kern.validate_energy()
        kern.save(HERE / f"kernel_{ion}_ng{ng}.npz")
        table_kb = (HERE / f"kernel_{ion}_ng{ng}.npz").stat().st_size / 1024
        t0 = time.time(); rows = []; Msum = 0; specs = []; errs = []; bol = []; epp = []
        for s in SEEDS:
            res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, N, "sobolev_group",
                         seed=s, t_core=T_CORE, kernel=kern)
            rows.append({b: band_ratio(res, *nu_of(*wd), weight="energy")[0] for b, wd in BANDS.items()})
            sp, se = spectrum(res, sp_edges, weight="energy"); specs.append(sp); errs.append(se)
            a = res["accounting"]; bol.append(a["E_esc"] / a["E_inj"])
            epp.append(float(res["n_events"].mean())); Msum = Msum + end_to_end_matrix(res, m_edges)
        wall = time.time() - t0
        M = Msum / len(SEEDS)
        dF = {b: float(np.mean([r[b] for r in rows]) / ref_bands[b] - 1) for b in BANDS}
        dbol = float(np.mean(bol) / ref_bol - 1)
        sp = np.mean(specs, 0); se = np.mean(errs, 0) / np.sqrt(len(SEEDS))
        ok = (sp_ref["err"] > 0) & np.isfinite(sp) & np.isfinite(sp_ref["spec"])
        chi2 = float(np.nansum((sp[ok] - sp_ref["spec"][ok]) ** 2 / (se[ok] ** 2 + sp_ref["err"][ok] ** 2)) / ok.sum())
        row_l1 = float(np.abs(M - M_ref).sum(axis=1).mean())
        w_row = M_ref.sum(axis=1)
        ew_err = float((np.abs(M - M_ref).sum(axis=1) * w_row).sum() / max(w_row.sum(), 1e-300))
        out["runs"][str(ng)] = dict(energy_residual=eres, dF=dF, dbol=dbol, chi2_dof=chi2,
                                    row_l1=row_l1, energy_weighted_matrix_err=ew_err,
                                    blocks=band_rows(M, str(ng)), wall_s=wall,
                                    table_kb=table_kb, events_per_packet=float(np.mean(epp)),
                                    bands={b: [r[b] for r in rows] for b in BANDS},
                                    empty_rows=int(kern.empty_rows.sum()))
        worst = max(abs(v) for v in dF.values())
        print(f"  Ng={ng:4d}: worst |dF_b| {100*worst:5.2f}%  dbol {100*dbol:+.2f}%  "
              f"chi2/dof {chi2:6.1f}  rowL1 {row_l1:.4f}  table {table_kb:.0f} kB  "
              f"[{wall:.0f}s]  " + " ".join(f"{b}={100*v:+.1f}%" for b, v in dF.items()), flush=True)
    # Gate 1 verdict
    for thresh, label in ((0.02, "excellent"), (0.05, "strong")):
        ngs = [ng for ng in groups if ng <= 32 and
               max(abs(v) for v in out["runs"][str(ng)]["dF"].values()) < thresh]
        if ngs:
            print(f"GATE 1: {label} -- Ng={min(ngs)} achieves |dF_b| < {100*thresh:.0f}% in every band")
            out["gate1"] = {"verdict": label, "ng": min(ngs)}
            break
    else:
        fine = [ng for ng in groups if max(abs(v) for v in out["runs"][str(ng)]["dF"].values()) < 0.05]
        out["gate1"] = {"verdict": "weak" if fine else "failure", "ng": (min(fine) if fine else None)}
        print(f"GATE 1: {out['gate1']}")
    (HERE / f"compression_{ion}.json").write_text(json.dumps(out, indent=1))
    print(f"wrote compression_{ion}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="laII", choices=["laII", "ceII", "ndII"])
    ap.add_argument("--groups", type=int, nargs="+", default=[4, 8, 16, 32, 64, 128])
    a = ap.parse_args(); main(a.ion, a.groups)
