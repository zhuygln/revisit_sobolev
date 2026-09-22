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
sys.path.insert(0, str(ROOT))
# G1's analyse.py by explicit path: both gates call their analysis module
# `analyse`, so a plain import resolves to whichever directory comes first
# on the path -- this file itself when it is imported from elsewhere.
import importlib.util as _ilu                                 # noqa: E402
_spec = _ilu.spec_from_file_location("paperB_gate1_analyse", ROOT / "paperB/gate1/analyse.py")
A = _ilu.module_from_spec(_spec); _spec.loader.exec_module(A)   # G1's metrics, live bands, gray conditions

PREREG = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, seeds=[1, 2, 3], build_seeds=[101, 102, 103],
              n={"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000},
              k_grid=[1, 2, 4, 8, 16, 32], f_grid=[0.1, 0.2, 0.5, 0.9, 0.99, 0.999], ng_control=8, ng_fine=128,
              dm_max=0.10, dcolour_max=0.10, h1_red_ratio=4.0,
              h2_green=0.10, h2_yellow=0.25, h2_red=0.50,
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
    # ---- H1: archetype counts only (the PI's amendment of 2026-09-22) ----
    # K*_local = G1's N_g*, K*_global = the smallest passing rank; the
    # event-level fit is a DIAGNOSTIC and never enters this logic.
    k_local = next((n for n in sorted(local) if local[n]["passes"]), None)
    k_global = next((k for k in sorted(glob) if glob[k]["passes"] and "gray" not in glob[k]), None)
    diag = None
    if k_local is not None and k_local in glob and k_local in local:
        g, l = glob[k_local], local[k_local]
        diag = dict(k=k_local, event_global=g["event"], event_local=l["event"], band_global=g["band_max"], band_local=l["band_max"],
                    global_fits_events_better=bool(g["event"] < l["event"]), local_reproduces_transport_better=bool(g["band_max"] > l["band_max"]),
                    events_vs_observables=bool(g["event"] < l["event"] and g["band_max"] > l["band_max"]))
    if k_local is None:
        h1 = "GRAY"; gray.append("9: no K*_local in G1's record (the local family never passes)")
    elif k_global is not None and k_global <= k_local / P["h1_red_ratio"]:
        h1 = "RED"
    elif k_global is None or k_global >= k_local:
        h1 = "GREEN"
    else:
        h1 = "YELLOW"
    # ---- H2: the retained fraction of distinct exit lines ----
    # L* = the smallest retained fraction that passes; the ladder is read
    # Red first, then Green, then Yellow. La II is a control, never decisive.
    passing = [(m["rho_exit"], f) for f, m in trunc.items() if m["passes"]]
    L_star, f_star = min(passing) if passing else (None, None)
    if L_star is None or L_star > P["h2_red"]:
        h2 = "RED"
    elif L_star <= P["h2_green"]:
        h2 = "GREEN"
    else:
        h2 = "YELLOW"
    return dict(gray=gray, live_bands=live, seed_std_R2=noise, ref_dev_from_G1=ref_dev, control_sigma=float(dev),
                local=local, global_nmf=glob, truncation=trunc, k_local=k_local, k_global=k_global, ng_star=k_local,
                k_star=k_global, diagnostic=diag, H1=h1 if not gray else "GRAY",
                f_star=f_star, L_star=L_star, rho_exit=L_star, event_at_k_local=local[k_local]["event"] if k_local in local else None,
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
    # both ladders: Red if either decisive ion is Red, Green if both are, else Yellow.
    # The decisive ions are Ce II and Nd II; La II is reported and never counted.
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
        L = f"{r['L_star']:.3f}" if r["L_star"] is not None else "—"
        print(f"\n{ion}: live {r['live_bands']} gray {r['gray'] or 'none'}; K*_local {r['k_local']}  K*_global {r['k_global']}  f* {r['f_star']}  L* {L}")
        print("  local  N : " + "  ".join(f"{n}:{m['band_max']:.3f}/{m['event']:.2f}" for n, m in sorted(r["local"].items())))
        print("  global k : " + "  ".join(f"{k}:{m['band_max']:.3f}/{m['event']:.2f}" for k, m in sorted(r["global_nmf"].items())) + "   (max|dm| / m_event)")
        print("  trunc  f : " + "  ".join(f"{f:g}:{m['band_max']:.3f}/{m['rho_exit']:.3f}" for f, m in sorted(r["truncation"].items())) + "   (max|dm| / retained fraction)")
        d = r["diagnostic"]
        if d:
            print(f"  diagnostic at k = {d['k']}: m_event global {d['event_global']:.3f} vs local {d['event_local']:.3f}; "
                  f"max|dm| global {d['band_global']:.3f} vs local {d['band_local']:.3f} -> events_vs_observables {d['events_vs_observables']}")
        print(f"  H1 {r['H1']}  H2 {r['H2']}  control {r['control_sigma']:.1f} sigma  ref dev {r['ref_dev_from_G1']:.1e}")
    print(f"\nH1 {out['H1']}  H2 {out['H2']}  -> {out['decision']}")
    return out


if __name__ == "__main__":
    main()
