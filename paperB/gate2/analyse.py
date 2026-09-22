#!/usr/bin/env python3
"""Paper B gate G2: the readings H1 (locality, not global low rank) and H2
(the observable-relevant dimension), gray first, exactly as preregistered
in paperB/prl_gate.md "G2". The reference is G1's R2 (rerun seed for seed
and checked against G1's record); the local family is G1's A2_ng{N} legs
(read from the gate1 records) with the L128 control in the G2 run.

    .venv/bin/python paperB/gate2/analyse.py        # gate2_<ion>.json + gate1 records -> gate2_verdict.json
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1"):
    sys.path.insert(0, str(p))
import analyse as A                                          # noqa: E402  (G1's metrics, live bands, gray conditions)

PREREG = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, seeds=[1, 2, 3], build_seeds=[101, 102, 103],
              n={"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000},
              k_grid=[1, 2, 4, 8, 16, 32], f_grid=[0.1, 0.2, 0.5, 0.9, 0.99, 0.999], ng_control=8, ng_fine=128,
              dm_max=0.10, dcolour_max=0.10, h1_red_ratio=4.0, h2_rho_green=0.10,
              nmf_conv_max=1e-4, control_sigma=2.0, ref_match=1e-6, decisive=["58CeII", "60NdII"])
IONS = A.IONS
GATE1 = ROOT / "paperB/gate1"


def check_prereg(row, strict=True):
    bad = []
    for k in ("state", "shell", "seeds", "build_seeds", "k_grid", "f_grid", "ng_control", "ng_fine"):
        if row.get(k) != PREREG[k]:
            bad.append(f"{k}: {row.get(k)!r} != {PREREG[k]!r}")
    if row.get("n") != PREREG["n"][row.get("ion")]:
        bad.append(f"n: {row.get('n')} != {PREREG['n'].get(row.get('ion'))}")
    if bad and strict:
        raise ValueError("record is not the preregistered G2 experiment: " + "; ".join(bad))
    return bad


def g1_record(ion):
    cands = sorted(GATE1.glob(f"gate1_{ion}*.json"), key=lambda p: json.loads(p.read_text())["n"])
    return json.loads(cands[-1].read_text()) if cands else None


def passes(m):
    return bool(m["band"]["max"] <= PREREG["dm_max"] and m["colour"]["max"] <= PREREG["dcolour_max"])


def read_ion(row, g1):
    P = PREREG
    M = A.metrics(row)                                     # live bands on R2, band/colour/sed vs R2, event vs the independent K128
    gray = A.gray_checks(row, M)                           # G1's conditions 1-4
    live = M["live_bands"]; noise = M["seed_std_R2"]
    # the reference must be G1's, seed for seed
    ref_dev = max((abs(row["legs"]["R2"]["mags"][b] - g1["legs"]["R2"]["mags"][b]) for b in live), default=0.0) if g1 else float("nan")
    if not g1 or not np.isfinite(ref_dev) or ref_dev > P["ref_match"]:
        gray.append(f"6: R2 differs from G1's record by {ref_dev:.2e} mag (or no G1 record)")
    # the local control: L128 must equal A2 within the noise on a difference of two legs
    ctrl = f"L128_ng{row['ng_control']}"; a2 = f"A2_ng{row['ng_control']}"
    dev = max((abs(row["legs"][ctrl]["mags"][b] - row["legs"][a2]["mags"][b]) / (np.sqrt(2.0) * max(noise[b], 1e-9)) for b in live), default=0.0)
    if dev > P["control_sigma"]:
        gray.append(f"7: the L128 control differs from A2 by {dev:.1f} sigma of the difference")
    # the families
    G1M = A.metrics(g1) if g1 else None
    local = {}
    if G1M:
        for tag, m in G1M["legs"].items():
            if tag.startswith("A2_ng") and "band" in m:
                local[int(m["ng"])] = dict(archetypes=int(m["ng"]), n_params=int(m["ng"]) ** 2, band_max=m["band"]["max"], band_mean=m["band"]["mean"],
                                           colour_max=m["colour"]["max"], sed=m["sed"], event=m["event"], n_exit=m["n_exit_samples"], table_kb=m["table_kb"],
                                           passes=passes(m))
    glob, trunc, nmf_gray = {}, {}, []
    for tag, m in M["legs"].items():
        tr = row["kernels"].get(tag, {}).get("transform")
        if not tr or "band" not in m:
            continue
        entry = dict(archetypes=tr["archetypes"], n_params=tr["n_params"], band_max=m["band"]["max"], band_mean=m["band"]["mean"],
                     colour_max=m["colour"]["max"], sed=m["sed"], event=m["event"], n_exit=m["n_exit_samples"], table_kb=m["table_kb"],
                     in_sample_tv=tr.get("in_sample_tv"), passes=passes(m), fallback_frac=m["fallback_frac"])
        if tr["kind"] == "nmf":
            entry.update(rel_frobenius=tr["rel_frobenius"], rel_change_tail=tr["rel_change_tail"])
            if tr["rel_change_tail"] > P["nmf_conv_max"]:
                nmf_gray.append(tr["k"]); entry["gray"] = "nmf not converged"
            glob[int(tr["k"])] = entry
        elif tr["kind"] == "truncate":
            entry.update(f=tr["f"], n_exit_total=tr["n_exit_total"], rho_exit=tr["n_exit_kept"] / tr["n_exit_total"])
            trunc[float(tr["f"])] = entry
    for k in nmf_gray:
        gray.append(f"8: G_k{k} NMF not converged (relative change over the last 100 iterations > {P['nmf_conv_max']})")
    # H1
    ng_star = next((n for n in sorted(local) if local[n]["passes"]), None)
    k_star = next((k for k in sorted(glob) if glob[k]["passes"] and "gray" not in glob[k]), None)
    matched = None
    if ng_star is not None and ng_star in glob and ng_star in local:
        matched = dict(k=ng_star, event_global=glob[ng_star]["event"], event_local=local[ng_star]["event"],
                       global_fits_events_better=bool(glob[ng_star]["event"] < local[ng_star]["event"]))
    if ng_star is None:
        h1 = "GRAY"; gray.append("9: no N_g* in G1's record")
    elif k_star is not None and k_star <= ng_star / P["h1_red_ratio"]:
        h1 = "RED"
    elif (k_star is None or k_star >= ng_star) and matched and matched["global_fits_events_better"]:
        h1 = "GREEN"
    else:
        h1 = "YELLOW"
    # H2
    f_star = next((f for f in sorted(trunc) if trunc[f]["passes"]), None)
    rho = trunc[f_star]["rho_exit"] if f_star is not None else None
    h2 = "RED" if f_star is None else "GREEN" if rho <= P["h2_rho_green"] else "YELLOW"
    event_at_ng_star = local[ng_star]["event"] if ng_star in local else None
    return dict(gray=gray, live_bands=live, seed_std_R2=noise, ref_dev_from_G1=ref_dev, control_sigma=float(dev),
                local=local, global_nmf=glob, truncation=trunc, ng_star=ng_star, k_star=k_star, matched=matched,
                H1=h1 if not gray else "GRAY", f_star=f_star, rho_exit=rho, event_at_ng_star=event_at_ng_star,
                H2=h2 if not gray else "GRAY", fine_in_vs_out_of_sample=M["fine_in_vs_out_of_sample"])


def readings(records, g1_records):
    P = PREREG
    per_ion = {ion: read_ion(row, g1_records.get(ion)) for ion, row in records.items()}
    dec = [i for i in P["decisive"] if i in per_ion]
    readable = [i for i in dec if not per_ion[i]["gray"]]
    out = dict(prereg=PREREG, per_ion=per_ion, decisive=dec, ions_gray=[i for i in per_ion if per_ion[i]["gray"]])
    if len(readable) < len(P["decisive"]):
        out.update(H1="GRAY", H2="GRAY", decision="GRAY", reason=f"decisive ions readable: {readable}")
        return out
    h1s = [per_ion[i]["H1"] for i in readable]; h2s = [per_ion[i]["H2"] for i in readable]
    H1 = "RED" if "RED" in h1s else "GREEN" if all(h == "GREEN" for h in h1s) else "YELLOW"
    H2 = "RED" if "RED" in h2s else "GREEN" if all(h == "GREEN" for h in h2s) else "YELLOW"
    decision = "WRITE" if H1 == "GREEN" else "REFRAME" if H1 == "RED" else "PI"
    out.update(H1=H1, H2=H2, decision=decision)
    return out


def main():
    records, g1s = {}, {}
    for ion in IONS:
        p = HERE / f"gate2_{ion}.json"
        if not p.exists():
            print(f"{ion}: no record"); continue
        row = json.loads(p.read_text()); check_prereg(row); records[ion] = row; g1s[ion] = g1_record(ion)
    out = readings(records, g1s)
    (HERE / "gate2_verdict.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for ion, r in out["per_ion"].items():
        print(f"\n{ion}: live {r['live_bands']} gray {r['gray'] or 'none'}; N_g* {r['ng_star']}  k* {r['k_star']}  f* {r['f_star']}  rho_exit {r['rho_exit']}")
        print("  local  N : " + "  ".join(f"{n}:{m['band_max']:.3f}/{m['event']:.2f}" for n, m in sorted(r["local"].items())))
        print("  global k : " + "  ".join(f"{k}:{m['band_max']:.3f}/{m['event']:.2f}" for k, m in sorted(r["global_nmf"].items())) + "   (max|dm| / m_event)")
        print("  trunc  f : " + "  ".join(f"{f:g}:{m['band_max']:.3f}/{m['rho_exit']:.3f}" for f, m in sorted(r["truncation"].items())) + "   (max|dm| / rho_exit)")
        print(f"  H1 {r['H1']}  H2 {r['H2']}  control {r['control_sigma']:.1f} sigma  ref dev {r['ref_dev_from_G1']:.1e}")
    print(f"\nH1 {out['H1']}  H2 {out['H2']}  -> {out['decision']}")
    return out


if __name__ == "__main__":
    main()
