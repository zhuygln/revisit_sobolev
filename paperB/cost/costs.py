#!/usr/bin/env python3
"""What the effective operator costs, from the existing records (the
referee's practical-utility point, the PI's 2026-10-06 decision): no new
transport. Three layers, kept apart:

  offline  -- generate the reference macroatom events and build R: the
              build reference's wall time (build seeds x packets), its event
              count;
  stored   -- the operator: the N_g x N_g matrix against the discrete exit
              tables (bytes, exit lines);
  online   -- transport with the operator against transport with the
              macroatom reference, same seeds and packets (wall time,
              events per packet), with the scalar closure for scale.

Per ion G1's final record (the state of Papers B's gates); the 13-ion blend
from G3 part (c); the radiation-field-driven macroatom's timing from the
R2M robustness check. And the limitation stated with the numbers: the
operator is built from the reference's own events, so at a new state or
mixture it does not yet remove the reference calculation.

    .venv/bin/python paperB/cost/costs.py          # -> paperB/cost/costs.json; --markdown
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
IONS = ("57LaII", "58CeII", "60NdII")
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II", "blend13": "13-ion blend"}


def g1_record(ion):
    cands = sorted((ROOT / "paperB/gate1").glob(f"gate1_{ion}*.json"), key=lambda p: json.loads(p.read_text())["n"])
    return (f"paperB/gate1/{cands[-1].name}", json.loads(cands[-1].read_text())) if cands else (None, None)


def layers(row, op_tag, build_tag="R2build", ref_tag="R2", scalar_tag="E0.10"):
    L, K = row["legs"], row["kernels"]
    b, r, o = L[build_tag], L[ref_tag], L[op_tag]
    k = K[op_tag]
    out = dict(
        n=row["n"], seeds=len(r["seeds"]), build_seeds=len(b["seeds"]),
        offline=dict(build_wall_s=b["t_wall"], build_packets=row["n"] * len(b["seeds"]), build_events_per_packet=b["events_per_packet"],
                     kernel_events=k["n_events"]),
        stored=dict(ng=k["ng"], matrix_entries=k["ng"] ** 2, matrix_bytes=8 * k["ng"] ** 2, exit_lines=k["n_exit_samples"],
                    serialized_bytes=k["serialized_bytes"], serialized_kb=k["serialized_bytes"] / 1024.0),
        online=dict(ref_wall_s=r["t_wall"], op_wall_s=o["t_wall"], op_over_ref=o["t_wall"] / r["t_wall"],
                    ref_events_per_packet=r["events_per_packet"], op_events_per_packet=o["events_per_packet"],
                    ref_interactions=r["n_interactions"], op_interactions=o["n_interactions"]))
    if scalar_tag in L:
        s = L[scalar_tag]
        out["online"].update(scalar_tag=scalar_tag, scalar_wall_s=s["t_wall"], scalar_over_ref=s["t_wall"] / r["t_wall"],
                             scalar_events_per_packet=s["events_per_packet"])
    return out


def main():
    g1v = json.loads((ROOT / "paperB/gate1/gate1_verdict.json").read_text())
    aud = json.loads((ROOT / "paperB/audit/exit_tables.json").read_text())
    out = dict(note="the operator is built from the reference's own events: at a new state or mixture it does not remove the reference calculation",
               per_record={})
    for ion in IONS:
        rel, row = g1_record(ion)
        ng = g1v["per_ion"][ion]["ng_star"]
        d = layers(row, f"A2_ng{ng}")
        d.update(record=rel, op_tag=f"A2_ng{ng}", ng_star=ng,
                 fine=dict(ng=128, bytes_total=aud[ion]["K128"]["bytes_total"], bytes_matrix=aud[ion]["K128"]["bytes_matrix"],
                           bytes_tables=aud[ion]["K128"]["bytes_tables"], n_exit=aud[ion]["n_distinct_exit_lines"]))
        out["per_record"][ion] = d
    pc = ROOT / "paperB/gate3/gate3_partc_p1blend.json"
    if pc.exists():
        row = json.loads(pc.read_text()); g3v = json.loads((ROOT / "paperB/gate3/gate3_verdict.json").read_text())
        kr = g3v["part_c"]["k_rec"]
        d = layers(row, f"Arec_ng{kr}"); d.update(record="paperB/gate3/gate3_partc_p1blend.json", op_tag=f"Arec_ng{kr}", ng_star=kr, n_ions=len(row["ions"]))
        out["per_record"]["blend13"] = d
    r2m = ROOT / "paperB/r2m/r2m_60NdII.json"
    if r2m.exists():
        row = json.loads(r2m.read_text()); L = row["legs"]
        out["r2m_60NdII"] = dict(n=row["n"], R2_wall_s=L["R2"]["t_wall"], R2M_wall_s=L["R2M"]["t_wall"],
                                 R2_events_per_packet=L["R2"]["events_per_packet"], R2M_events_per_packet=L["R2M"]["events_per_packet"],
                                 A2_ng8_wall_s=L["A2_ng8"]["t_wall"], A2M_ng8_wall_s=L["A2M_ng8"]["t_wall"])
    (HERE / "costs.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for key, d in out["per_record"].items():
        o, s, f = d["online"], d["stored"], d["offline"]
        print(f"{NAME[key]}: build {f['build_wall_s']:.0f} s ({f['kernel_events']:,} events); stored N_g={s['ng']}: {s['serialized_kb']:.0f} kB "
              f"({s['exit_lines']:,} exit lines, matrix {s['matrix_bytes']} B); online op/ref {o['op_wall_s']:.0f}/{o['ref_wall_s']:.0f} s = {o['op_over_ref']:.2f}"
              + (f"; scalar {o['scalar_wall_s']:.0f} s" if "scalar_wall_s" in o else ""))
    return out


if __name__ == "__main__" and "--markdown" not in sys.argv:
    main()


def markdown(out=None):
    out = out or json.loads((HERE / "costs.json").read_text())
    L = ["| record | packets × seeds | offline: build run [s] | kernel events | stored: N_g | matrix [B] | exit lines | operator [kB] | online: reference [s] | operator [s] | operator / reference | scalar ε = 0.10 [s] |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for key, d in out["per_record"].items():
        o, s, f = d["online"], d["stored"], d["offline"]
        L.append(f"| {NAME[key]} | {d['n']:,} × {d['seeds']} | {f['build_wall_s']:.0f} | {f['kernel_events']:,} | {s['ng']} | {s['matrix_bytes']:,} | "
                 f"{s['exit_lines']:,} | {s['serialized_kb']:.0f} | {o['ref_wall_s']:.0f} | {o['op_wall_s']:.0f} | {o['op_over_ref']:.2f} | "
                 + (f"{o['scalar_wall_s']:.0f} |" if "scalar_wall_s" in o else "— |"))
    if "r2m_60NdII" in out:
        r = out["r2m_60NdII"]
        L += ["", f"Nd II under the radiation-field-driven macroatom (R2M, {r['n']:,} packets): reference {r['R2_wall_s']:.0f} s downward against "
                  f"{r['R2M_wall_s']:.0f} s radiation-field-driven ({r['R2_events_per_packet']:.1f} against {r['R2M_events_per_packet']:.1f} events per packet); "
                  f"the 8-group operators {r['A2_ng8_wall_s']:.0f} s and {r['A2M_ng8_wall_s']:.0f} s."]
    return "\n".join(L)


if __name__ == "__main__" and "--markdown" in sys.argv:
    print(markdown())
