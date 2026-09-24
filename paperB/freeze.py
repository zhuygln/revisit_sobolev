#!/usr/bin/env python3
"""Collect the Paper B headline numbers from the committed gate records into
paperB/FROZEN.json. Every number in the manuscript is a LaTeX macro generated
from this file by docs/paperB/latex_tables.py; nothing in the prose is typed
by hand (the standing rule after Paper IV's fabricated-table incident).

    .venv/bin/python paperB/freeze.py            # regenerate FROZEN.json
    .venv/bin/python paperB/freeze.py --check    # regenerate in memory and compare

Sources: paperB/gate1/gate1_verdict.json (G1, F67), paperB/r2m/r2m_verdict.json
(the robustness check, F68), paperB/audit/exit_tables.json (the exit-table
audit, F68) and paperB/gate2/gate2_verdict.json (G2, F69). A pure function of
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

    h = dict(sources={g1_rel: g1_sha, r2m_rel: r2m_sha, aud_rel: aud_sha, g2_rel: g2_sha},
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
            local=[dict(k=int(k), band=v["band_max"], event=v["event"], n_params=v["n_params"])
                   for k, v in sorted(r["local"].items(), key=lambda kv: int(kv[0]))],
            glob=[dict(k=int(k), band=v["band_max"], event=v["event"], n_params=v["n_params"])
                  for k, v in sorted(r["global_nmf"].items(), key=lambda kv: int(kv[0]))],
            trunc=[dict(f=float(f), band=v["band_max"], rho=v["rho_exit"], n_exit=v["n_exit"])
                   for f, v in sorted(r["truncation"].items(), key=lambda kv: float(kv[0]))],
            diagnostic=r["diagnostic"])
    # the headline contrast at the largest archetype count both families reach
    K = 32

    def at(ion, fam, k):
        return next(e for e in h["g2"]["per_ion"][ion][fam] if e["k"] == k)

    h["g2"]["contrast"] = {ion: dict(
        k=K, event_global=at(ion, "glob", K)["event"], event_local=at(ion, "local", K)["event"],
        band_global=at(ion, "glob", K)["band"], band_local=at(ion, "local", K)["band"],
        params_global=at(ion, "glob", K)["n_params"], params_local=at(ion, "local", K)["n_params"],
        event_ratio=at(ion, "local", K)["event"] / at(ion, "glob", K)["event"],
        band_ratio=at(ion, "glob", K)["band"] / at(ion, "local", K)["band"]) for ion in IONS}
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
