#!/usr/bin/env python3
"""The R2M robustness readings (paperB/prl_gate.md "R2M robustness"), gray
first, with R2M as the reference leg: the reference shift R2M - R2, whether
the 8-group operator rebuilt from R2M's events survives G1's thresholds,
the reference dependence of G1's downward-trained operator, and the
event-level distances against the R2M fine matrix.

    .venv/bin/python paperB/r2m/analyse_r2m.py          # all r2m_<ion>.json present -> r2m_verdict.json
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1"):
    sys.path.insert(0, str(p))
import analyse as A                                          # noqa: E402  (G1's metrics and live-band rule)

PREREG = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, n=300_000, seeds=[1, 2, 3], build_seeds=[101, 102, 103],
              ng_check=8, ng_fine=128, macro_W=0.5, dm_max=0.10, dcolour_max=0.10,
              gray_min_live=2, gray_identity=1e-10, gray_kernel_energy=1e-12, gray_trapped_frac=0.01, gray_fallback_frac=0.01)
IONS = A.IONS
REF = "R2M"


def check_prereg(row, strict=True):
    bad = []
    for k in ("state", "shell", "ng_check", "ng_fine", "macro_W", "build_seeds", "seeds"):
        if row.get(k) != PREREG[k]:
            bad.append(f"{k}: {row.get(k)!r} != {PREREG[k]!r}")
    if row.get("n", 0) < PREREG["n"]:
        bad.append(f"n: {row.get('n')} < {PREREG['n']}")
    if bad and strict:
        raise ValueError("record is not the preregistered R2M check: " + "; ".join(bad))
    return bad


def gray_checks(row, live, dropped):
    P = PREREG; fired = []
    if len(live) < P["gray_min_live"]:
        fired.append(f"1: {len(live)} live bands on {REF}")
    if dropped:
        fired.append(f"2: band(s) {dropped} carry >= 1 % of L_bol but exceed the precision rule on {REF}: raise the packet count")
    for tag, leg in row["legs"].items():
        ident = abs(leg["energy"]["identity_residual"])
        if ident > P["gray_identity"]:
            fired.append(f"3: {tag} identity residual {ident:.1e}")
        n_run = row["n"] * len(leg.get("seeds", row["seeds"]))
        if leg["n_trapped"] / n_run > P["gray_trapped_frac"]:
            fired.append(f"4: {tag} chain-capped {leg['n_trapped'] / n_run:.3f}")
        fb = leg.get("n_coherent_fallback", 0) / max(leg.get("n_interactions", 0), 1)
        if fb > P["gray_fallback_frac"]:
            fired.append(f"4: {tag} empty-row fallback {fb:.3f}")
    for tag, k in row["kernels"].items():
        if k["validate_energy"] > P["gray_kernel_energy"]:
            fired.append(f"3: {tag} kernel energy {k['validate_energy']:.1e}")
    return fired


def read_ion(row):
    P = PREREG; ng = row["ng_check"]; nf = row["ng_fine"]
    ref = row["legs"][REF]; r2 = row["legs"]["R2"]
    live, dropped, detectable = A.live_bands(row, REF)
    gray = gray_checks(row, live, dropped)
    edges = A.V.nu_edges(*row["lam_window"], row["n_spec"]); dnu = np.diff(edges)
    shift = {b: float(ref["mags"][b] - r2["mags"][b]) for b in live}
    fine = row["kernels"][f"K{nf}M"]; fine_build = row["kernels"][f"K{nf}Mbuild"]
    legs = {}
    for tag in (f"A2M_ng{ng}", f"A2_ng{ng}", "R2"):
        leg = row["legs"][tag]
        m = dict(band=A.m_band(leg, ref, live), colour=A.m_colour(leg, ref, live), sed=A.m_sed(leg, ref, dnu),
                 events_per_packet=leg["events_per_packet"], t_wall=leg["t_wall"])
        if tag in row["kernels"]:
            k = row["kernels"][tag]
            m.update(event_vs_K128M=A.m_event(k, fine), kernel_source=k["source"], n_exit_samples=k.get("n_exit_samples"),
                     table_kb=k.get("table_kb"), empty_rows=k["empty_rows"])
        legs[tag] = m
    a2m = legs[f"A2M_ng{ng}"]
    survives = None if gray else bool(a2m["band"]["max"] <= P["dm_max"] and a2m["colour"]["max"] <= P["dcolour_max"])
    return dict(gray=gray, live_bands=live, dropped_for_precision=dropped, detectable_40mpc=detectable,
                seed_std_R2M={b: float(ref["mags_seed_std"][b]) for b in live},
                seed_std_R2={b: float(r2["mags_seed_std"][b]) for b in live},
                shift_R2M_minus_R2=shift, shift_max=float(max(abs(v) for v in shift.values())) if shift else float("nan"),
                events_per_packet=dict(R2M=ref["events_per_packet"], R2=r2["events_per_packet"]),
                ledger=dict(R2M=ref["ledger"], R2=r2["ledger"]),
                survives="GRAY" if gray else ("YES" if survives else "NO"), legs=legs,
                fine_in_vs_out_of_sample=A.m_event(fine_build, fine),
                n_exit_K128M=fine.get("n_exit_samples"), table_kb_K128M=fine.get("table_kb"))


def readings(records):
    per_ion = {ion: read_ion(row) for ion, row in records.items()}
    return dict(prereg=PREREG, reference=REF, per_ion=per_ion,
                survives={ion: r["survives"] for ion, r in per_ion.items()})


def main():
    records = {}
    for ion in IONS:
        cands = sorted(HERE.glob(f"r2m_{ion}*.json"), key=lambda p: json.loads(p.read_text())["n"])
        if not cands:
            print(f"{ion}: no record"); continue
        row = json.loads(cands[-1].read_text()); check_prereg(row); records[ion] = row
        print(f"{ion}: {cands[-1].name} (n = {row['n']})")
    out = readings(records)
    (HERE / "r2m_verdict.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for ion, r in out["per_ion"].items():
        print(f"\n{ion}: live {r['live_bands']} dropped {r['dropped_for_precision']} gray {r['gray'] or 'none'}")
        print(f"  R2M - R2: " + " ".join(f"{b} {v:+.2f}" for b, v in r["shift_R2M_minus_R2"].items()) + f"  (max {r['shift_max']:.2f}); events/pkt R2M {r['events_per_packet']['R2M']:.1f} vs R2 {r['events_per_packet']['R2']:.1f}")
        for tag, m in r["legs"].items():
            ev = f" m_event {m['event_vs_K128M']:.3f}" if "event_vs_K128M" in m else ""
            print(f"  {tag:8s} vs R2M: max|dm| {m['band']['max']:.3f} mean {m['band']['mean']:.3f} colour {m['colour']['max']:.3f} sed {m['sed']:.3f}{ev}")
        print(f"  survives: {r['survives']}   (fine in/out-of-sample floor {r['fine_in_vs_out_of_sample']:.3f})")
    return out


if __name__ == "__main__":
    main()
