#!/usr/bin/env python3
"""Collect the Paper B headline numbers from the committed gate records into
paperB/FROZEN.json. Every number in the manuscript is a LaTeX macro generated
from this file by docs/paperB/latex_tables.py; nothing in the prose is typed
by hand (the standing rule after Paper IV's fabricated-table incident).

    .venv/bin/python paperB/freeze.py            # regenerate FROZEN.json
    .venv/bin/python paperB/freeze.py --check    # regenerate in memory and compare

Sources: paperB/gate1/gate1_verdict.json (G1, F67), paperB/r2m/r2m_verdict.json
(the robustness check, F68), paperB/audit/exit_tables.json (the exit-table
audit, F68), paperB/gate2/gate2_verdict.json (G2, F69), paperB/gate3/gate3_verdict.json
(G3, F70), paperB/scalar/minimax.json (the joint minimax scalar, 2026-10-06) and
paperB/cost/costs.json (the cost layers, 2026-10-06). A pure function of
the committed JSONs: no transport is run.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "paperB/FROZEN.json"
IONS = ("57LaII", "58CeII", "60NdII")
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}
DECISIVE = ("58CeII", "60NdII")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


def load(rel):
    p = ROOT / rel
    return json.loads(p.read_text()), rel, sha(p)


def build():
    g1, g1_rel, g1_sha = load("paperB/gate1/gate1_verdict.json")
    r2m, r2m_rel, r2m_sha = load("paperB/r2m/r2m_verdict.json")
    aud, aud_rel, aud_sha = load("paperB/audit/exit_tables.json")
    g2, g2_rel, g2_sha = load("paperB/gate2/gate2_verdict.json")
    g3, g3_rel, g3_sha = load("paperB/gate3/gate3_verdict.json")
    mm, mm_rel, mm_sha = load("paperB/scalar/minimax.json")
    co, co_rel, co_sha = load("paperB/cost/costs.json")
    g2r, g2r_rel, g2r_sha = load("paperB/gate2r/g2r_verdict.json")

    h = dict(sources={g1_rel: g1_sha, r2m_rel: r2m_sha, aud_rel: aud_sha, g2_rel: g2_sha, g3_rel: g3_sha, mm_rel: mm_sha, co_rel: co_sha,
                      g2r_rel: g2r_sha},
             ions=list(IONS), ion_names=NAME, decisive=list(DECISIVE))

    # ---- G1 (F67) ----
    h["g1"] = dict(B1=g1["B1"], B2=g1["B2"], B3=g1["B3"], decision=g1["decision"],
                   n=g1["prereg"]["n"], dm_max=g1["prereg"]["dm_max"], dcolour_max=g1["prereg"]["dcolour_max"],
                   per_ion={})
    for ion in IONS:
        r = g1["per_ion"][ion]
        h["g1"]["per_ion"][ion] = dict(
            ng_star=r["ng_star"], live=r["live_bands"], e_eps=r["e_eps"], e_R=r["e_R"],
            eps_star=r["eps_star"]["eps"], eps_max_dm=r["eps_star"]["max_dm"], eps_sed=r["eps_star"]["sed"],
            events_per_packet=r["R2"]["events_per_packet"],
            table=[dict(ng=t["ng"], band_max=t["band_max"], band_mean=t["band_mean"], sed=t["sed"],
                        event=t["event"], colour_max=t["colour_max"], n_exit=t["n_exit_samples"],
                        table_kb=t["table_kb"]) for t in r["table"]],
            seed_std_max=max(r["seed_std_R2"].values()))

    # ---- the R2M robustness check (F68) ----
    h["r2m"] = dict(per_ion={})
    for ion in IONS:
        r = r2m["per_ion"][ion]
        h["r2m"]["per_ion"][ion] = dict(
            survives=r["survives"], shift_max=r["shift_max"], shift=r["shift_R2M_minus_R2"],
            ev_R2M=r["events_per_packet"]["R2M"], ev_R2=r["events_per_packet"]["R2"],
            rebuilt_band=r["legs"]["A2M_ng8"]["band"]["max"], rebuilt_colour=r["legs"]["A2M_ng8"]["colour"]["max"],
            rebuilt_sed=r["legs"]["A2M_ng8"]["sed"], rebuilt_event=r["legs"]["A2M_ng8"]["event_vs_K128M"],
            downward_band=r["legs"]["A2_ng8"]["band"]["max"], r2_band=r["legs"]["R2"]["band"]["max"],
            n_exit_rebuilt=r["legs"]["A2M_ng8"]["n_exit_samples"], n_exit_downward=r["legs"]["A2_ng8"]["n_exit_samples"])

    # ---- the exit-table audit (F68) ----
    h["audit"] = dict(per_ion={})
    for ion in IONS:
        a = aud[ion]
        h["audit"]["per_ion"][ion] = dict(
            n_events=a["n_events"], n_lines=a["n_lines"], n_opacity=a["n_opacity"],
            n_exit=a["n_distinct_exit_lines"], exact=a["exact_line_fraction"],
            conc90=a["concentration"]["0.9"], conc99=a["concentration"]["0.99"],
            kb_total=a["K128"]["bytes_total"] / 1024.0, kb_matrix=a["K128"]["bytes_matrix"] / 1024.0,
            kb_tables=a["K128"]["bytes_tables"] / 1024.0,
            discovery=[[d["events"], d["distinct"]] for d in a["discovery"]])

    # ---- G2 (F69) ----
    h["g2"] = dict(H1=g2["H1"], H2=g2["H2"], decision=g2["decision"], per_ion={})
    for ion in IONS:
        r = g2["per_ion"][ion]
        h["g2"]["per_ion"][ion] = dict(
            k_local=r["k_local"], k_global=r["k_global"], L_star=r["L_star"], f_star=r["f_star"],
            H1=r["H1"], H2=r["H2"], live=r["live_bands"],
            local=[dict(k=int(k), band=v["band_max"], colour=v["colour_max"], joint=max(v["band_max"], v["colour_max"]),
                        event=v["event"], n_params=v["n_params"], passes=v["passes"])
                   for k, v in sorted(r["local"].items(), key=lambda kv: int(kv[0]))],
            glob=[dict(k=int(k), band=v["band_max"], colour=v["colour_max"], joint=max(v["band_max"], v["colour_max"]),
                       event=v["event"], n_params=v["n_params"], passes=v["passes"])
                  for k, v in sorted(r["global_nmf"].items(), key=lambda kv: int(kv[0]))],
            trunc=[dict(f=float(f), band=v["band_max"], rho=v["rho_exit"], n_exit=v["n_exit"])
                   for f, v in sorted(r["truncation"].items(), key=lambda kv: float(kv[0]))],
            diagnostic=r["diagnostic"])
    # the headline contrast at the largest archetype count both families reach
    K = 32

    def at(ion, fam, k):
        return next(e for e in h["g2"]["per_ion"][ion][fam] if e["k"] == k)

    def r_diag(ion):
        d = g2["per_ion"][ion]["diagnostic"]
        return None if not d else bool(d["events_vs_observables"])

    def r_diag_k(ion):
        d = g2["per_ion"][ion]["diagnostic"]
        return None if not d else int(d["k"])

    h["g2"]["contrast"] = {ion: dict(
        k=K, event_global=at(ion, "glob", K)["event"], event_local=at(ion, "local", K)["event"],
        band_global=at(ion, "glob", K)["band"], band_local=at(ion, "local", K)["band"],
        colour_global=at(ion, "glob", K)["colour"], colour_local=at(ion, "local", K)["colour"],
        matched_k=r_diag_k(ion), events_vs_observables_at_matched=r_diag(ion),
        params_global=at(ion, "glob", K)["n_params"], params_local=at(ion, "local", K)["n_params"],
        event_ratio=at(ion, "local", K)["event"] / at(ion, "glob", K)["event"],
        band_ratio=at(ion, "glob", K)["band"] / at(ion, "local", K)["band"]) for ion in IONS}
    # ---- G3 (F70): state transfer, existence, interpolation, mixtures ----
    P3 = g3["prereg"]
    h["g3"] = dict(C1=g3["C1"], C2=g3["C2"], C3=g3["C3"], C4=g3["C4"], C5=g3["C5"], step4=g3["step4"],
                   ng_t=P3["ng_t"], n_axes=len(P3["axes"]), axes=list(P3["axes"]),
                   n_states_per_ion=sum(len(a["grid"]) for a in P3["axes"].values()), per_ion={})
    for ion in IONS:
        r = g3["per_ion"][ion]
        sts = [x for x in g3["states"] if x["ion"] == ion and x["axis"] != "ref"]
        fresh = []
        for x in sts:
            r16 = x["recomputed"].get(str(P3["ng_t"])) or x["recomputed"].get(P3["ng_t"])
            t = x["transfer"]; i = x["interpolation"]
            fresh.append(dict(axis=x["axis"], label=x["label"], cls=x["classification"], coord_value=x["coord_value"],
                              k_rec=x["k_rec"], exists=bool(x["exists_at_ng_t"]), gray=bool(x["gray"]),
                              fresh_band=r16["band"], fresh_colour=r16["colour"],
                              anchor_band=t["band"] if t else None, anchor_colour=t["colour"] if t else None,
                              rows_never_trained=t["rows_never_trained_frac"] if t else None,
                              interp_whole_band=i["whole"]["band"] if i else None,
                              interp_whole_colour=i["whole"]["colour"] if i else None,
                              interp_matrix_band=i["matrix_only"]["band"] if i else None,
                              interp_matrix_colour=i["matrix_only"]["colour"] if i else None, lam=i["lam"] if i else None,
                              fresh_joint=max(r16["band"], r16["colour"]),
                              anchor_joint=max(t["band"], t["colour"]) if t else None,
                              interp_whole_joint=max(i["whole"]["band"], i["whole"]["colour"]) if i else None,
                              interp_matrix_joint=max(i["matrix_only"]["band"], i["matrix_only"]["colour"]) if i else None))
        readable = [f for f in fresh if not f["gray"]]
        failing = [f for f in readable if not f["exists"]]
        h["g3"]["per_ion"][ion] = dict(
            C1=r["C1"], C2=r["C2"], C3=r["C3"], needed_axes=r["needed_axes"], interpolation_fails=r["interpolation_fails"],
            transfer_holds_axes=[ax for ax, a in r["axes"].items() if a["transfer_holds"]],
            transfer_fails_axes=[ax for ax, a in r["axes"].items() if not a["transfer_holds"]],
            interp_passes_axes=[ax for ax, a in r["axes"].items() if a["interpolation_passes"] is True],
            n_states=len(fresh), n_readable=len(readable), n_gray=len(fresh) - len(readable),
            n_fresh_fail=len(failing), n_fresh_fail_band=sum(1 for f in failing if f["fresh_band"] > P3["dm_max"]),
            worst_fresh_band=max(f["fresh_band"] for f in readable), worst_fresh_colour=max(f["fresh_colour"] for f in readable),
            worst_fail_colour_over=max((f["fresh_colour"] - P3["dcolour_max"] for f in failing), default=0.0),
            least_fail_colour_over=min((f["fresh_colour"] - P3["dcolour_max"] for f in failing), default=0.0),
            k_rec_max=max(f["k_rec"] for f in readable if f["k_rec"] is not None),
            states=fresh)
    # the coupled-trajectory boundary: Nd II at 3 d
    nd = next(f for f in h["g3"]["per_ion"]["60NdII"]["states"] if f["label"] == "P3d")
    h["g3"]["nd_trajectory"] = dict(label=nd["label"], anchor_band=nd["anchor_band"], anchor_colour=nd["anchor_colour"],
                                    interp_whole_band=nd["interp_whole_band"], interp_whole_colour=nd["interp_whole_colour"],
                                    interp_matrix_band=nd["interp_matrix_band"], interp_matrix_colour=nd["interp_matrix_colour"],
                                    fresh_band=nd["fresh_band"], fresh_colour=nd["fresh_colour"], k_rec=nd["k_rec"],
                                    matrix_only_passes=bool(nd["interp_matrix_band"] <= P3["dm_max"] and nd["interp_matrix_colour"] <= P3["dcolour_max"]))
    pb, pc = g3["part_b"], g3["part_c"]
    amix_valid = {int(n): m for n, m in pb["table"]["Amix"].items() if not m.get("invalid")}
    best_n = min(amix_valid, key=lambda n: amix_valid[n]["band"])
    kd = pb["k_direct"]
    h["g3"]["part_b"] = dict(k_mix=pb["k_mix"], k_direct=kd, invalid_legs=pb["invalid_legs"],
                             amix_best_n=best_n, amix_best_band=amix_valid[best_n]["band"], amix_best_colour=amix_valid[best_n]["colour"],
                             adirect_band=pb["table"]["Adirect"][str(kd)]["band"], adirect_colour=pb["table"]["Adirect"][str(kd)]["colour"])
    kr = pc["k_rec"]
    h["g3"]["part_c"] = dict(k_rec=kr, band=pc["table"]["Arec"][str(kr)]["band"], colour=pc["table"]["Arec"][str(kr)]["colour"],
                             eps_star=pc["eps_star"]["eps"], eps_err=pc["eps_star"]["max_dm"], n_ions=13)
    # ---- the fair scalar comparator (the referee's point, 2026-10-06): the joint minimax eps from the existing grids ----
    h["minimax"] = dict(rule=mm["rule"], per_record={})
    for key, r in mm["per_record"].items():
        h["minimax"]["per_record"][key] = dict(
            record=r["record"], live=r["live_bands"],
            prereg_eps=r["prereg"]["eps"], prereg_band=r["prereg"]["max_band"], prereg_colour=r["prereg"]["max_colour"],
            prereg_joint=r["prereg"]["joint"],
            band_eps=r["band_minimax"]["eps"], band_band=r["band_minimax"]["max_band"],
            mm_eps=r["joint_minimax"]["eps"], mm_band=r["joint_minimax"]["max_band"], mm_colour=r["joint_minimax"]["max_colour"],
            mm_joint=r["joint_minimax"]["joint"], mm_interior=r["joint_minimax"]["interior"])
    # ---- what the operator costs, from the records (offline / stored / online) ----
    h["cost"] = dict(note=co["note"], per_record={}, r2m_60NdII=co.get("r2m_60NdII"))
    for key, d in co["per_record"].items():
        h["cost"]["per_record"][key] = dict(
            record=d["record"], op_tag=d["op_tag"], ng=d["stored"]["ng"], n=d["n"], seeds=d["seeds"],
            build_s=d["offline"]["build_wall_s"], build_packets=d["offline"]["build_packets"], kernel_events=d["offline"]["kernel_events"],
            matrix_bytes=d["stored"]["matrix_bytes"], exit_lines=d["stored"]["exit_lines"], op_kb=d["stored"]["serialized_kb"],
            ref_s=d["online"]["ref_wall_s"], op_s=d["online"]["op_wall_s"], op_over_ref=d["online"]["op_over_ref"],
            scalar_s=d["online"].get("scalar_wall_s"), scalar_over_ref=d["online"].get("scalar_over_ref"),
            ref_ev_per_pkt=d["online"]["ref_events_per_packet"], op_ev_per_pkt=d["online"]["op_events_per_packet"])
    ratios = [d["op_over_ref"] for d in h["cost"]["per_record"].values()]
    h["cost"]["op_over_ref_min"] = min(ratios); h["cost"]["op_over_ref_max"] = max(ratios)
    # ---- G2R (F71): the frequency-adjacency ablation ----
    P = g2r["prereg"]
    h["g2r"] = dict(reading=g2r["reading"], decision=g2r["decision"], n_scrambled=P["n_scrambled"], n_orderings=P["n_scrambled"] + 1,
                    k_grid=list(P["k_grid"]), per_ion={})
    for ion in g2r["per_ion"]:
        r = g2r["per_ion"][ion]
        cells = [dict(k=int(k), e_physical=c["e_physical"], band_physical=c["band_physical"], colour_physical=c["colour_physical"],
                      rank=c["rank"], p=c["p_value"], scr_min=c["scrambled"]["min"], scr_median=c["scrambled"]["median"],
                      scr_max=c["scrambled"]["max"], scr_frac_pass=c["scrambled"]["frac_pass"], physical_passes=c["physical_passes"],
                      event_physical=c["event_physical"], event_scr_median=c["event_scrambled"]["median"], event_rank=c["event_rank"],
                      scrambled=[o["e_joint"] for o in c["orderings"] if o["perm_id"] != 0],
                      scrambled_event=[o["event"] for o in c["orderings"] if o["perm_id"] != 0])
                 for k, c in sorted(r["cells"].items(), key=lambda kv: int(kv[0]))]
        at = next(c for c in cells if c["k"] == r["k_local"])
        h["g2r"]["per_ion"][ion] = dict(reading=r["reading"], k_local=r["k_local"], rank_at_k_local=r["rank_at_k_local"],
                                        n_rank1_cells=r["n_rank1_cells"], n_cells=len(cells), gray=r["gray"], cells=cells,
                                        at_k_local=dict(e_physical=at["e_physical"], scr_min=at["scr_min"], scr_median=at["scr_median"],
                                                        scr_frac_pass=at["scr_frac_pass"], p=at["p"]),
                                        first_cells=[c["k"] for c in cells if c["rank"] == 1],
                                        not_first_cells=[c["k"] for c in cells if c["rank"] != 1])
    dec = [i for i in DECISIVE if i in h["g2r"]["per_ion"]]
    h["g2r"]["n_cells"] = sum(h["g2r"]["per_ion"][i]["n_cells"] for i in dec)
    h["g2r"]["n_first"] = sum(h["g2r"]["per_ion"][i]["n_rank1_cells"] for i in dec)
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    h = build()
    txt = json.dumps(h, indent=1, default=float, sort_keys=True) + "\n"
    if a.check:
        if not OUT.exists():
            print("FROZEN.json missing"); sys.exit(1)
        if OUT.read_text() != txt:
            print("FROZEN.json is not the regeneration of the committed records"); sys.exit(1)
        print(f"freeze check OK: {len(h['sources'])} sources, {len(IONS)} ions")
        return
    OUT.write_text(txt)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
