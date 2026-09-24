#!/usr/bin/env python3
"""Paper B gate G3: per state the four-way reading A/B/C/D and the support
diagnostics, then C1-C5 over the domain, gray first, as preregistered in
paperB/prl_gate.md "G3". Every comparison is against R2 at the SAME state.

    .venv/bin/python paperB/gate3/analyse.py        # every gate3_*.json present -> gate3_verdict.json; --markdown
"""
import importlib.util as _ilu
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
_spec = _ilu.spec_from_file_location("paperB_gate1_analyse", ROOT / "paperB/gate1/analyse.py")
A = _ilu.module_from_spec(_spec); _spec.loader.exec_module(A)          # G1's metrics, live bands, gray conditions
import states as S                                                     # noqa: E402

PREREG = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, seeds=[1, 2, 3], build_seeds=[101, 102, 103],
              n={"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000}, n_blend=1_000_000,
              ng_grid=[2, 4, 8, 16, 32], ng_t=16, ng_fine=128,
              axes={ax: dict(grid=a["grid"], interior=a["interior"], bracket=list(a["bracket"]), coord=a["coord"]) for ax, a in S.AXES.items()},
              dm_max=0.10, dcolour_max=0.10, c1_small_ng=8, c1_small_frac=2 / 3, coverage_flag=0.05,
              decisive=["58CeII", "60NdII"], controlled_axes=["T", "D", "J"])
IONS = A.IONS


def check_prereg(row, strict=True):
    bad = []
    for k in ("seeds", "build_seeds", "ng_grid", "ng_t", "ng_fine"):
        if row.get(k) != PREREG[k]:
            bad.append(f"{k}: {row.get(k)!r} != {PREREG[k]!r}")
    if row.get("n") != PREREG["n"][row["ion"]]:
        bad.append(f"n: {row.get('n')} != {PREREG['n'][row['ion']]}")
    if row["axis"] != "ref" and row["value"] not in PREREG["axes"][row["axis"]]["grid"]:
        bad.append(f"state {row['axis']}={row['value']} not on the grid")
    if bad and strict:
        raise ValueError("record is not the preregistered G3 experiment: " + "; ".join(bad))
    return bad


def passes(m):
    return bool(m["band"]["max"] <= PREREG["dm_max"] and m["colour"]["max"] <= PREREG["dcolour_max"])


def read_state(row):
    """One (axis, state, ion) record -> its reading against R2 at that state."""
    P = PREREG; ng_t = row["ng_t"]
    M = A.metrics(row)
    live = M["live_bands"]
    # gray: G1's conditions on R2 and the recomputed legs only; the transfer and
    # interpolation legs' fallback IS the failure mechanism and is reported below
    M_gray = dict(M, legs={t: m for t, m in M["legs"].items() if not t.startswith(("Afix", "Aint"))})
    gray = A.gray_checks(row, M_gray)
    legs = M["legs"]
    rec = {int(m["ng"]): m for t, m in legs.items() if t.startswith("Arec_ng") and "band" in m}
    k_rec = next((n for n in sorted(rec) if passes(rec[n])), None)
    exists = bool(ng_t in rec and passes(rec[ng_t]))
    fix = legs.get(f"Afix_ng{ng_t}"); whole = legs.get(f"Aint_ng{ng_t}"); mat = legs.get(f"AintM_ng{ng_t}")
    transfers = None if fix is None else passes(fix)
    interpolates = None if whole is None else passes(whole)
    if row["axis"] == "ref":
        cls = "REF"
    elif transfers:
        cls = "A"
    elif interpolates:
        cls = "B"
    elif exists:
        cls = "C"
    else:
        cls = "D"
    sup = row["kernels"].get(f"K{row['ng_fine']}", {}).get("transform", {})
    fixk = row["kernels"].get(f"Afix_ng{ng_t}", {})
    out = dict(axis=row["axis"], value=row["value"], label=row["label"], ion=row["ion"], coord=row.get("coord"),
               coord_value=row.get("coord_value"), gray=gray, live_bands=live, seed_std_R2=M["seed_std_R2"],
               classification=cls, k_rec=k_rec, exists_at_ng_t=exists, transfers=transfers, interpolates=interpolates,
               recomputed={n: dict(band=m["band"]["max"], colour=m["colour"]["max"], sed=m["sed"], event=m["event"]) for n, m in sorted(rec.items())},
               transfer=None if fix is None else dict(band=fix["band"]["max"], colour=fix["colour"]["max"], sed=fix["sed"], event=fix["event"],
                                                        rows_never_trained_frac=fix["fallback_frac"], anchor_empty_rows=fixk.get("empty_rows")),
               interpolation=None if whole is None else dict(
                   whole=dict(band=whole["band"]["max"], colour=whole["colour"]["max"], sed=whole["sed"], event=whole["event"], fallback=whole["fallback_frac"]),
                   matrix_only=None if mat is None else dict(band=mat["band"]["max"], colour=mat["colour"]["max"], sed=mat["sed"], event=mat["event"]),
                   lam=row["interpolation"]["lam"], endpoints=row["interpolation"]["endpoints"],
                   transform=row["kernels"].get(f"Aint_ng{ng_t}", {}).get("transform")),
               outside_fixed_support=dict(clipped_frac=sup.get("clipped_frac"), clipped_energy_frac=sup.get("clipped_energy_frac"),
                                          coverage_limited=bool((sup.get("clipped_frac") or 0.0) > P["coverage_flag"])),
               fine_in_vs_out_of_sample=M["fine_in_vs_out_of_sample"], events_per_packet=legs["R2"]["events_per_packet"],
               t_wall=row.get("t_wall_total"))
    if gray:
        out["classification"] = "GRAY"
    return out


def readings(states, part_b=None, part_c=None):
    """states: [read_state(row)]. -> C1..C5 with the A/B/C/D table."""
    P = PREREG
    by = {}
    for s in states:
        by.setdefault(s["ion"], {}).setdefault(s["axis"], {})[s["label"]] = s
    per_ion = {}
    for ion, axes in by.items():
        r = dict(axes={}, C1=None, C2=None, C3=None)
        ok_all, small, n_states, undefined = True, 0, 0, 0
        for ax, a in P["axes"].items():
            st = axes.get(ax, {})
            labels = [S.label(ax, v) for v in a["grid"]]
            present = [st[l] for l in labels if l in st]
            readable = [s for s in present if not s["gray"]]
            interior = st.get(S.label(ax, a["interior"]))
            transfer_holds = bool(readable) and all(s["transfers"] for s in readable) and len(readable) == len(labels)
            interp_ok = None if not interior or interior["gray"] else interior["interpolates"]
            r["axes"][ax] = dict(states={s["label"]: s["classification"] for s in present}, n_expected=len(labels), n_present=len(present),
                                 n_gray=len(present) - len(readable), transfer_holds=transfer_holds, interpolation_passes=interp_ok,
                                 needed=(not transfer_holds), table=[dict(label=s["label"], cls=s["classification"], k_rec=s["k_rec"],
                                                                          transfer_band=s["transfer"]["band"] if s["transfer"] else None,
                                                                          rows_never_trained=s["transfer"]["rows_never_trained_frac"] if s["transfer"] else None,
                                                                          clipped=s["outside_fixed_support"]["clipped_frac"],
                                                                          rec16=s["recomputed"].get(16, {}).get("band")) for s in present])
            for s in readable:
                n_states += 1
                if not s["exists_at_ng_t"]:
                    undefined += 1
                if s["k_rec"] is not None and s["k_rec"] <= P["c1_small_ng"]:
                    small += 1
            ok_all &= all(s["exists_at_ng_t"] for s in readable)
        # C1 existence
        r["C1"] = "RED" if undefined >= 2 else "YELLOW" if (undefined == 1 or (n_states and small / n_states < P["c1_small_frac"])) else "GREEN"
        # C2 transfer
        failing = [ax for ax, v in r["axes"].items() if not v["transfer_holds"]]
        if not failing:
            r["C2"] = "GREEN"
        elif all(r["axes"][ax]["interpolation_passes"] for ax in failing):
            r["C2"] = "YELLOW"
        else:
            r["C2"] = "RED"
        # C3 sufficiency over the controlled axes
        needed = [ax for ax in P["controlled_axes"] if r["axes"].get(ax, {}).get("needed")]
        interp_fail = [ax for ax in P["controlled_axes"] if r["axes"].get(ax, {}).get("interpolation_passes") is False]
        r["needed_axes"] = needed; r["interpolation_fails"] = interp_fail
        r["C3"] = "RED" if len(interp_fail) >= 2 else "GREEN" if len(needed) <= 2 else "YELLOW"
        per_ion[ion] = r
    out = dict(prereg=PREREG, per_ion=per_ion)
    dec = [i for i in P["decisive"] if i in per_ion]
    for c in ("C1", "C2", "C3"):
        vals = [per_ion[i][c] for i in dec]
        out[c] = "GRAY" if len(dec) < 2 else "RED" if "RED" in vals else "GREEN" if all(v == "GREEN" for v in vals) else "YELLOW"
    # C4 composability, C5 the realistic mixture
    out["C4"] = out["C5"] = "NOT RUN"
    if part_b:
        out["part_b"] = part_b
        out["C4"] = "GREEN" if part_b["k_mix"] is not None and part_b["k_direct"] is not None and part_b["k_mix"] <= 2 * part_b["k_direct"] \
            else "YELLOW" if part_b["k_direct"] is not None else "RED"
    if part_c:
        out["part_c"] = part_c
        k = part_c["k_rec"]
        out["C5"] = "GREEN" if k is not None and k <= 8 else "YELLOW" if k is not None else "RED"
    c2 = out["C2"]
    out["step4"] = ("transfers across states" if c2 == "GREEN" else
                    "exists everywhere and is tabulable in a small state vector" if (out["C1"] == "GREEN" and c2 == "YELLOW" and out["C3"] in ("GREEN", "YELLOW")) else
                    "reframed by the PI" if out["C1"] == "RED" else "the PI decides")
    return out


def read_blend(row, kind):
    """Part (b): k_mix / k_direct; part (c): k_rec and eps*."""
    M = A.metrics(row); gray = A.gray_checks(row, M)
    fam = {}
    for t, m in M["legs"].items():
        if "band" in m and "ng" in m:
            fam.setdefault(t.split("_ng")[0], {})[int(m["ng"])] = m
    out = dict(gray=gray, live_bands=M["live_bands"],
               table={f: {n: dict(band=m["band"]["max"], colour=m["colour"]["max"], event=m["event"]) for n, m in sorted(v.items())} for f, v in fam.items()})
    if kind == "b":
        out["k_mix"] = next((n for n in sorted(fam.get("Amix", {})) if passes(fam["Amix"][n])), None)
        out["k_direct"] = next((n for n in sorted(fam.get("Adirect", {})) if passes(fam["Adirect"][n])), None)
    else:
        out["k_rec"] = next((n for n in sorted(fam.get("Arec", {})) if passes(fam["Arec"][n])), None)
        es = A.eps_star(M, row) if any("eps" in m for m in M["legs"].values()) else None
        out["eps_star"] = None if es is None else dict(eps=es["eps"], max_dm=es["max_dm"])
    return out


def main():
    states, records = [], {}
    for p in sorted(HERE.glob("gate3_*_*.json")):
        if p.name.startswith(("gate3_partb", "gate3_partc", "gate3_verdict", "gate3_support")):
            continue
        row = json.loads(p.read_text()); check_prereg(row); records[p.name] = row
        states.append(read_state(row))
    pb = pc = None
    if (HERE / "gate3_partb_blend3.json").exists():
        pb = read_blend(json.loads((HERE / "gate3_partb_blend3.json").read_text()), "b")
    if (HERE / "gate3_partc_p1blend.json").exists():
        pc = read_blend(json.loads((HERE / "gate3_partc_p1blend.json").read_text()), "c")
    out = readings(states, pb, pc)
    out["states"] = states
    (HERE / "gate3_verdict.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for ion, r in out["per_ion"].items():
        print(f"\n{ion}: C1 {r['C1']}  C2 {r['C2']}  C3 {r['C3']}  needed {r['needed_axes']}  interp fails {r['interpolation_fails']}")
        for ax, v in r["axes"].items():
            print(f"  {ax}: " + "  ".join(f"{t['label']}:{t['cls']}(rec16 {t['rec16']:.3f}" + (f", fix {t['transfer_band']:.3f}, untrained {t['rows_never_trained']:.4f}" if t['transfer_band'] is not None else "") + ")"
                                       for t in v["table"]) + f"  | transfer {'holds' if v['transfer_holds'] else 'FAILS'}, interp {v['interpolation_passes']}")
    print(f"\nC1 {out['C1']}  C2 {out['C2']}  C3 {out['C3']}  C4 {out['C4']}  C5 {out['C5']}  -> step 4: {out['step4']}")
    return out


if __name__ == "__main__":
    main()
