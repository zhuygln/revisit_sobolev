#!/usr/bin/env python3
"""G2R: the frequency-adjacency ablation read as preregistered in
paperB/prl_gate.md "G2R", gray first. Per (ion, N_g) the joint transport
error E_joint = max(max |dm| over the live bands, max |dcolour|) of the
physical ordering against the 31 scrambled orderings: its rank, the
permutation p-value, the scrambled distribution, and the event-level loss
of every ordering alongside.

    .venv/bin/python paperB/gate2r/analyse.py            # g2r_<ion>_k<k>.json -> g2r_verdict.json; --markdown
"""
import importlib.util as _ilu
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
_spec = _ilu.spec_from_file_location("paperB_gate1_analyse", ROOT / "paperB/gate1/analyse.py")
A = _ilu.module_from_spec(_spec); _spec.loader.exec_module(A)          # G1's live bands, band/colour metrics, gray conditions

PREREG = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, seeds=[1, 2, 3], build_seeds=[101, 102, 103],
              n={"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000}, ng_fine=128,
              k_grid=[2, 4, 8, 16, 32], n_scrambled=31, perm_seed0=1000,
              perm_seeds=[None] + [1000 + m for m in range(1, 32)],
              dm_max=0.10, dcolour_max=0.10, ref_match=1e-6, physical_match=1e-6,
              green_rank=1, green_min_cells=4, red_rank=8, decisive=["58CeII", "60NdII"])
IONS = A.IONS
GATE2 = ROOT / "paperB/gate2"
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}


def check_prereg(row, strict=True):
    bad = []
    for k in ("state", "shell", "seeds", "build_seeds", "ng_fine", "n_scrambled", "perm_seed0", "perm_seeds"):
        if row.get(k) != PREREG[k]:
            bad.append(f"{k}: {row.get(k)!r} != {PREREG[k]!r}")
    if row.get("k") not in PREREG["k_grid"]:
        bad.append(f"k: {row.get('k')} not on the grid")
    if row.get("n") != PREREG["n"][row.get("ion")]:
        bad.append(f"n: {row.get('n')} != {PREREG['n'].get(row.get('ion'))}")
    if bad and strict:
        raise ValueError("record is not the preregistered G2R experiment: " + "; ".join(bad))
    return bad


def metrics(row):
    """G1's metrics on a G2R record. The scrambled orderings' 128x128
    matrices were dropped from the record after their event-level loss was
    computed (run_g2r.py); the stored loss stands in for m_event there."""
    ref = row["legs"]["R2"]
    live, dropped, detectable = A.live_bands(row, "R2")
    noise = A.seed_noise(row, live)
    fine = row["kernels"][f"K{row['ng_fine']}"]; fine_build = row["kernels"].get(f"K{row['ng_fine']}build")
    edges = A.V.nu_edges(*row["lam_window"], row["n_spec"]); dnu = np.diff(edges)
    out = dict(live_bands=live, dropped_for_precision=dropped, detectable_40mpc=detectable, seed_std_R2=noise, legs={},
               fine_in_vs_out_of_sample=A.m_event(fine_build, fine) if fine_build else None)
    for tag, leg in row["legs"].items():
        m = dict(mode=leg["mode"], identity=abs(leg["energy"]["identity_residual"]),
                 trapped_frac=leg["n_trapped"] / (row["n"] * len(leg.get("seeds", row["seeds"]))),
                 fallback_frac=leg.get("n_coherent_fallback", 0) / max(leg.get("n_interactions", 0), 1),
                 events_per_packet=leg["events_per_packet"], t_wall=leg["t_wall"], seeds=leg.get("seeds", row["seeds"]))
        if tag not in ("R2", "R2build"):
            m["band"] = A.m_band(leg, ref, live); m["sed"] = A.m_sed(leg, ref, dnu); m["colour"] = A.m_colour(leg, ref, live)
        if tag in row["kernels"]:
            k = row["kernels"][tag]
            m.update(ng=k["ng"], n_exit_samples=k.get("n_exit_samples"), table_kb=k.get("table_kb"), kernel_energy=k["validate_energy"],
                     empty_rows=k["empty_rows"], kernel_source=k["source"],
                     event=k["event"] if "event" in k else (A.m_event(k, fine) if "R" in k else float("nan")))
        out["legs"][tag] = m
    return out


def e_joint(m):
    return float(max(m["band"]["max"], m["colour"]["max"]))


def passes(m):
    return bool(m["band"]["max"] <= PREREG["dm_max"] and m["colour"]["max"] <= PREREG["dcolour_max"])


def g2_record(ion):
    p = GATE2 / f"gate2_{ion}.json"
    return json.loads(p.read_text()) if p.exists() else None


def read_cell(row, g2):
    """One (ion, N_g) record -> the statistic of its 32 orderings."""
    P = PREREG; k = row["k"]
    M = metrics(row); gray = A.gray_checks(row, M); live = M["live_bands"]
    ref_dev = max((abs(row["legs"]["R2"]["mags"][b] - g2["legs"]["R2"]["mags"][b]) for b in live), default=0.0) if g2 else float("nan")
    if not g2 or not np.isfinite(ref_dev) or ref_dev > P["ref_match"]:
        gray.append(f"6: R2 differs from G2's record by {ref_dev:.2e} mag (or no G2 record)")
    phys = f"P000_ng{k}"; l128 = f"L128_ng{k}"
    phys_dev = (max((abs(row["legs"][phys]["mags"][b] - g2["legs"][l128]["mags"][b]) for b in live), default=0.0)
                if g2 and l128 in g2["legs"] and phys in row["legs"] else float("nan"))
    if not np.isfinite(phys_dev) or phys_dev > P["physical_match"]:
        gray.append(f"7: the physical ordering differs from G2's {l128} by {phys_dev:.2e} mag (or no such leg)")
    orders = {}
    for tag, m in M["legs"].items():
        tr = row["kernels"].get(tag, {}).get("transform")
        if not tr or tr.get("kind") != "local_perm" or "band" not in m:
            continue
        orders[int(tr["perm_id"])] = dict(perm_id=int(tr["perm_id"]), perm_seed=tr.get("perm_seed"), e_joint=e_joint(m),
                                          band_max=m["band"]["max"], colour_max=m["colour"]["max"], sed=m["sed"], event=m["event"],
                                          in_sample_tv=tr.get("in_sample_tv"), n_params=tr["n_params"], passes=passes(m),
                                          fallback_frac=m["fallback_frac"], empty_rows=m["empty_rows"])
    if 0 not in orders or len(orders) != P["n_scrambled"] + 1:
        gray.append(f"10: {len(orders)} orderings in the record, {P['n_scrambled'] + 1} preregistered")
    e0 = orders[0]["e_joint"] if 0 in orders else float("nan")
    scr = np.array([orders[i]["e_joint"] for i in sorted(orders) if i != 0]) if orders else np.array([])
    ev0 = orders[0]["event"] if 0 in orders else float("nan")
    ev_scr = np.array([orders[i]["event"] for i in sorted(orders) if i != 0]) if orders else np.array([])
    rank = int(1 + np.sum(scr <= e0)) if scr.size else None            # ties count against the physical ordering
    return dict(k=k, gray=gray, live_bands=live, seed_std_R2=M["seed_std_R2"], ref_dev_from_G2=ref_dev, physical_dev_from_G2=phys_dev,
                e_physical=e0, band_physical=orders[0]["band_max"] if 0 in orders else None,
                colour_physical=orders[0]["colour_max"] if 0 in orders else None,
                rank=rank, p_value=None if rank is None else rank / (scr.size + 1),
                scrambled=dict(n=int(scr.size), min=float(scr.min()) if scr.size else None, median=float(np.median(scr)) if scr.size else None,
                               max=float(scr.max()) if scr.size else None, frac_pass=float(np.mean([orders[i]["passes"] for i in orders if i != 0])) if scr.size else None),
                physical_passes=orders[0]["passes"] if 0 in orders else None,
                event_physical=ev0, event_scrambled=dict(min=float(ev_scr.min()) if ev_scr.size else None, median=float(np.median(ev_scr)) if ev_scr.size else None,
                                                         max=float(ev_scr.max()) if ev_scr.size else None),
                event_rank=int(1 + np.sum(ev_scr <= ev0)) if ev_scr.size else None,
                fine_in_vs_out_of_sample=M["fine_in_vs_out_of_sample"], orderings=[orders[i] for i in sorted(orders)])


def read_ion(cells, k_local):
    """cells: {k: read_cell}. The preregistered ladder on the rank of the
    physical ordering at K*_local and across the N_g grid."""
    P = PREREG
    gray = [f"N_g = {k}: {g}" for k, c in sorted(cells.items()) for g in c["gray"]]
    missing = [k for k in P["k_grid"] if k not in cells]
    if missing:
        gray.append(f"11: no record for N_g = {missing}")
    ranks = {k: c["rank"] for k, c in cells.items()}
    at = ranks.get(k_local)
    n_rank1 = sum(1 for r in ranks.values() if r == P["green_rank"])
    if gray or k_local is None or at is None:
        reading = "GRAY"
        if k_local is None:
            gray.append("12: no K*_local in G2's verdict")
    elif at > P["red_rank"]:
        reading = "RED"
    elif at == P["green_rank"] and n_rank1 >= P["green_min_cells"]:
        reading = "GREEN"
    else:
        reading = "YELLOW"
    return dict(gray=gray, k_local=k_local, rank_at_k_local=at, n_rank1_cells=n_rank1, ranks=ranks, reading=reading,
                cells={k: cells[k] for k in sorted(cells)})


def readings(per_ion_cells, k_locals):
    P = PREREG
    per_ion = {ion: read_ion(cells, k_locals.get(ion)) for ion, cells in per_ion_cells.items()}
    dec = [i for i in P["decisive"] if i in per_ion]
    out = dict(prereg=PREREG, per_ion=per_ion, decisive=dec, ions_gray=[i for i in per_ion if per_ion[i]["gray"]])
    readable = [i for i in dec if per_ion[i]["reading"] != "GRAY"]
    if len(readable) < len(P["decisive"]):
        out.update(reading="GRAY", decision="GRAY", reason=f"decisive ions readable: {readable}")
        return out
    rs = [per_ion[i]["reading"] for i in readable]
    R = "RED" if "RED" in rs else "GREEN" if all(r == "GREEN" for r in rs) else "YELLOW"
    out.update(reading=R, decision={"GREEN": "ADJACENCY", "RED": "REFRAME", "YELLOW": "PI"}[R])
    return out


def main():
    g2v = json.loads((GATE2 / "gate2_verdict.json").read_text())
    k_locals = {ion: g2v["per_ion"][ion]["k_local"] for ion in g2v["per_ion"]}
    cells = {}
    for p in sorted(HERE.glob("g2r_*_k*.json")):
        row = json.loads(p.read_text()); check_prereg(row)
        cells.setdefault(row["ion"], {})[int(row["k"])] = read_cell(row, g2_record(row["ion"]))
    if not cells:
        print("no G2R records"); return None
    out = readings(cells, k_locals)
    (HERE / "g2r_verdict.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for ion, r in out["per_ion"].items():
        print(f"\n{ion}: K*_local {r['k_local']}  rank at K*_local {r['rank_at_k_local']}  rank-1 cells {r['n_rank1_cells']}/{len(r['cells'])}  "
              f"-> {r['reading']}" + (f"  gray: {r['gray']}" if r["gray"] else ""))
        for k, c in r["cells"].items():
            s = c["scrambled"]
            print(f"  N_g {k:2d}: physical E_joint {c['e_physical']:.3f} (rank {c['rank']}/{s['n'] + 1}, p {c['p_value']:.3f}); "
                  f"scrambled min/median/max {s['min']:.3f}/{s['median']:.3f}/{s['max']:.3f}, {s['frac_pass']:.2f} pass; "
                  f"m_event physical {c['event_physical']:.3f} vs scrambled median {c['event_scrambled']['median']:.3f} (rank {c['event_rank']})")
    print(f"\nG2R {out['reading']} -> {out['decision']}")
    return out


if __name__ == "__main__" and "--markdown" not in sys.argv:
    main()


def markdown(out=None):
    """The report's tables, machine-generated from g2r_verdict.json."""
    out = out or json.loads((HERE / "g2r_verdict.json").read_text())
    L = []
    for ion, r in out["per_ion"].items():
        L += [f"{NAME[ion]}: K*_local = {r['k_local']}, the physical ordering ranks {r['rank_at_k_local']} of 32 there and first in "
              f"{r['n_rank1_cells']} of {len(r['cells'])} block counts; **{r['reading']}**" + (f"; gray: {'; '.join(r['gray'])}" if r["gray"] else "") + ":", "",
              "| N_g | E_joint physical | rank / 32 | p | scrambled min | median | max | scrambled passing | m_event physical | m_event scrambled median | m_event rank |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        for k, c in sorted(r["cells"].items(), key=lambda kv: int(kv[0])):
            s = c["scrambled"]; e = c["event_scrambled"]
            L.append(f"| {k} | {c['e_physical']:.3f} | {c['rank']} | {c['p_value']:.3f} | {s['min']:.3f} | {s['median']:.3f} | {s['max']:.3f} | "
                     f"{s['frac_pass']:.2f} | {c['event_physical']:.3f} | {e['median']:.3f} | {c['event_rank']} |")
        L.append("")
    L.append(f"G2R reading **{out['reading']}** -> {out['decision']}.")
    return "\n".join(L)


if __name__ == "__main__" and "--markdown" in sys.argv:
    print(markdown())
