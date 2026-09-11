#!/usr/bin/env python3
"""Collect the Paper IV headline numbers from the committed run records into
paper4/FROZEN.json (the manuscript's numbers are LaTeX macros generated from
it by docs/paper4/latex_tables.py; nothing in the prose is typed by hand).

    .venv/bin/python paper4/freeze.py            # regenerate FROZEN.json
    .venv/bin/python paper4/freeze.py --check    # regenerate in memory and compare

Every entry names its source file; the file's SHA-256 is recorded. The
record is a pure function of the committed JSONs (no transport is run).
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
P4 = ROOT / "paper4"
OUT = P4 / "FROZEN.json"
BANDS = "grizJHK"
sys.path.insert(0, str(ROOT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def load(rel):
    p = ROOT / rel
    return json.loads(p.read_text()), rel, sha(p)


def dm(legs, a, r, bands=BANDS):
    return {b: float(legs[a]["mags"][b] - legs[r]["mags"][b]) for b in bands}


def headline():
    h = {"sources": {}}

    def note(rel, s):
        h["sources"][rel] = s

    # ---- the published states, corrected legs (F58, F61, F62) ----
    states = {"P1_t1": "P1 1 d", "P1_t2": "P1 2 d", "P1_t3": "P1 3 d", "P1_t5": "P1 5 d",
              "P1r1_t2": "P1r1", "P1r2_t2": "P1r2", "P2_t3.4": "P2 3.4 d"}
    h["states"] = {}
    for key, name in states.items():
        d, rel, s = load(f"paper4/phase9_final/legs_{key}.json"); note(rel, s)
        L = d["legs"]
        h["states"][key] = dict(name=name, shell=d.get("shell"), t_d=d.get("t_d"), n_lines=d.get("n_lines"),
                                n_opacity=d.get("n_opacity"), S_band=d.get("S_band"),
                                B2_R2=dm(L, "B2", "R2"), Bbin2_R2=dm(L, "Bbin2", "R2"), A2_R2=dm(L, "A2", "R2"),
                                Bth_Rth=dm(L, "Bth", "Rth"), Bbinth_Rth=dm(L, "Bbinth", "Rth"), Rth_R2=dm(L, "Rth", "R2"),
                                B1_R1=dm(L, "B1", "R1"), R1_R2=dm(L, "R1", "R2"), C2_R2=dm(L, "C2", "R2"))
    for tag, rel in (("gate2_P1", "paper4/phase9_final/gate2_P1_fixed.json"), ("gate2_P2", "paper4/phase9_final/gate2_P2_fixed.json"),
                     ("gate2_robust", "paper4/phase9_final/gate2_robust_fixed.json")):
        try:
            d, rel_, s = load(rel); note(rel_, s); h[tag] = d.get("verdict", d)
        except FileNotFoundError:
            pass
    # summary statistics used in the prose
    p1 = [h["states"][k] for k in ("P1_t1", "P1_t2", "P1_t3", "P1_t5")]
    h["summary"] = dict(
        B2_R2_z_P1_min=min(x["B2_R2"]["z"] for x in p1), B2_R2_z_P1_max=max(x["B2_R2"]["z"] for x in p1),
        B2_R2_K_P1_min=min(x["B2_R2"]["K"] for x in p1), B2_R2_K_P1_max=max(x["B2_R2"]["K"] for x in p1),
        Bbin2_R2_P1_min=min(min(x["Bbin2_R2"].values()) for x in p1), Bbin2_R2_P1_max=max(max(x["Bbin2_R2"].values()) for x in p1),
        Bbinth_Rth_P1_absmax=max(max(abs(v) for v in x["Bbinth_Rth"].values()) for x in p1),
        Bth_Rth_P1_min=min(min(x["Bth_Rth"].values()) for x in p1), Bth_Rth_P1_max=max(max(x["Bth_Rth"].values()) for x in p1),
        A2_R2_absmax=max(max(abs(v) for v in x["A2_R2"].values()) for x in h["states"].values()),
        B2_R2_P2_g=h["states"]["P2_t3.4"]["B2_R2"]["g"], B1_R1_P1_t2_absmax_zK=max(abs(h["states"]["P1_t2"]["B1_R1"][b]) for b in "zK"),
        R1_R2_P1_t2_z=h["states"]["P1_t2"]["R1_R2"]["z"])

    # ---- grid invariance (F60) ----
    h["grid"] = {}
    for tag, rel in (("single_s27", "paper4/phase3_legs/neighbour_P1_t2_s27_fixed.json"),
                     ("single_s28", "paper4/phase8_shells/single_P1_t2_s28_fixed.json"),
                     ("single_s29", "paper4/phase3_legs/neighbour_P1_t2_s29_fixed.json"),
                     ("shells4", "paper4/phase8_shells/shells_P1_t2_s28-31_fixed.json"),
                     ("split3", "paper4/phase8_shells/shells_P1_t2_split3_fixed.json"),
                     ("shells6", "paper4/phase8_shells/shells_P1_t2_f6_fixed.json"),
                     ("shells12", "paper4/phase8_shells/shells_P1_t2_f12_fixed.json"),
                     ("shells23", "paper4/phase8_shells/shells_P1_t2_f24_fixed.json"),
                     ("P2_shells32", "paper4/phase8_shells/shells_P2_t3.4_s0-31_fixed.json")):
        d, rel_, s = load(rel); note(rel_, s); L = d["legs"]
        h["grid"][tag] = dict(n_shell=d.get("n_shell", 1), B2_R2=dm(L, "B2", "R2"), Bbin2_R2=dm(L, "Bbin2", "R2"))
    h["grid"]["split_bit_identical"] = all(abs(h["grid"]["shells4"]["B2_R2"][b] - h["grid"]["split3"]["B2_R2"][b]) < 1e-9 for b in BANDS)
    d, rel, s = load("paper4/phase6_convergence/converge_P1_t2_s28_fixed.json"); note(rel, s)
    zs = [r["B2"]["dm_vs_R2"]["z"] for ax in d["axes"].values() for r in ax if "B2" in r and "dm_vs_R2" in r["B2"]]
    h["convergence"] = dict(B2_R2_z_min=min(zs), B2_R2_z_max=max(zs), n_settings=len(zs))

    # ---- the macroatom robustness test (F63) ----
    d, rel, s = load("paper4/phase10_fontes/macro_P1_t2_s28.json"); note(rel, s); L = d["legs"]
    h["macro"] = dict(B2_R2=dm(L, "B2", "R2"), Bbin2_R2=dm(L, "Bbin2", "R2"), B2M_R2M=dm(L, "B2M", "R2M"),
                      Bbin2M_R2M=dm(L, "Bbin2M", "R2M"), R2M_R2=dm(L, "R2M", "R2"),
                      ev_R2=L["R2"]["events_per_packet"], ev_R2M=L["R2M"]["events_per_packet"])

    # ---- the Fontes snapshot (F64) ----
    h["fontes_snapshot"] = {}
    for tag, rel in (("gsi_reemit", "paper4/phase10_fontes/fontes_t4_z45_reemit.json"),
                     ("gsi_reflect", "paper4/phase10_fontes/fontes_t4_z45_reflect.json"),
                     ("gsi_reflect_lossless", "paper4/phase10_fontes/fontes_t4_z45_reflect_lossless.json"),
                     ("jplt_reemit", "paper4/phase10_fontes/fontes_t4_z45_jplt_reemit.json"),
                     ("jplt_reflect", "paper4/phase10_fontes/fontes_t4_z45_jplt_reflect.json")):
        d, rel_, s = load(rel); note(rel_, s); L = d["legs"]
        e = {}
        for a, r in (("Bth", "Rth"), ("Bbinth", "Rth"), ("B2", "R2"), ("Bbin2", "R2"), ("Rth", "R2")):
            if a in L and r in L:
                e[f"{a}_{r}"] = dm(L, a, r)
                e[f"{a}_{r}_rK"] = float((L[a]["mags"]["r"] - L[a]["mags"]["K"]) - (L[r]["mags"]["r"] - L[r]["mags"]["K"]))
        e["ev"] = {k: L[k]["events_per_packet"] for k in L}
        e["n_lines"] = d.get("n_lines"); e["dataset"] = d.get("dataset")
        h["fontes_snapshot"][tag] = e
    h["fontes_snapshot"]["jplt_bbin2_terminated"] = False

    # ---- the light curve (F65) ----
    d, rel, s = load("paper4/phase10_fontes/prod_record/merged/summary.json"); note(rel, s)
    note("paper4/phase10_fontes/prod_record/merged/run.json", sha(ROOT / "paper4/phase10_fontes/prod_record/merged/run.json"))   # the summary derives from it
    lc = {}
    for leg, o in d["legs"].items():
        lc[leg] = dict(t_peak_d=o["t_peak_d"], L_peak=o["L_peak"], E_rad=o["E_rad"], W_tot=o["W_tot"], E_init=o["E_init"],
                       E_inj=o["E_inj"], E_end=o["E_end"], closure=o["closure"], f_capped_max=o["f_capped_max"])
        if "dL_peak" in o:
            mad = {}
            for b in BANDS:
                col = np.array([x[b] for x in o["dm_vs_ref"]], float)
                mad[b] = float(np.nanmax(np.abs(col))) if np.isfinite(col).any() else None
            lc[leg].update(dL_peak=o["dL_peak"], dt_peak=o["dt_peak"], dE_rad=o["dE_rad"], max_abs_dcolour=o["max_abs_dcolour"],
                           max_abs_dm=mad)
    h["lightcurve"] = dict(legs=lc, readings={k: (v.get("outcome") if isinstance(v, dict) else v) for k, v in d["readings"].items()},
                           t_grid_d=[float(t) / 86400.0 for t in d["t_grid"][:1] + d["t_grid"][-1:]], n_slabs=len(d["t_grid"]) - 1,
                           config={k: d["config"][k] for k in ("n_shell", "n_init", "n_heat", "max_events", "temperature", "init")})
    # the line-binned fluorescence leg at the two event caps (the merged record carries the 3e5 run)
    cap = {}
    for tag, rel_run in (("cap1e5", "paper4/phase10_fontes/prod_record/Bbin2/run.json"), ("cap3e5", "paper4/phase10_fontes/prod_record/Bbin2_cap3e5/run.json")):
        r, rel, s = load(rel_run); note(rel, s); T = r["tallies"]["Bbin2"]
        E_in = T[0]["E_carried_in"] + sum(x["E_inj_new"] for x in T)
        cap[tag] = dict(max_events=r["config"]["max_events"], f_capped_first=T[0]["f_capped"], E_capped_frac=sum(x["E_capped"] for x in T) / E_in,
                        events_mean_first=T[0]["events_mean"], E_rad=sum(x["E_esc"] for x in T), W_tot=sum(x["W"] for x in T))
    old, rel, s = load("paper4/phase10_fontes/prod_record/merged_cap1e5/summary.json"); note(rel, s)
    Lo = np.array(old["legs"]["Bbin2"]["L_esc"]); Ln = np.array(d["legs"]["Bbin2"]["L_esc"])
    cap["E_rad_change"] = cap["cap3e5"]["E_rad"] / cap["cap1e5"]["E_rad"] - 1.0
    cap["L_ratio_max_dev"] = float(np.max(np.abs(Ln / Lo - 1.0)))
    cap["L_peak_change"] = d["legs"]["Bbin2"]["L_peak"] / old["legs"]["Bbin2"]["L_peak"] - 1.0
    cap["cap1e5_dL_peak"] = old["legs"]["Bbin2"]["dL_peak"]; cap["cap1e5_max_abs_dcolour"] = old["legs"]["Bbin2"]["max_abs_dcolour"]
    h["lightcurve"]["bbin2_cap"] = cap
    try:
        d, rel, s = load("paper4/phase10_fontes/prod_record/merged_rad/summary.json"); note(rel, s)
        h["lightcurve"]["rad_variant"] = {leg: dict(t_peak_d=o["t_peak_d"], L_peak=o["L_peak"], dL_peak=o.get("dL_peak"), max_abs_dcolour=o.get("max_abs_dcolour"))
                                          for leg, o in d["legs"].items()}
    except FileNotFoundError:
        pass
    try:
        d, rel, s = load("paper4/phase10_fontes/prod_record/free_R2/run.json"); note(rel, s)
        t = d["tallies"]["R2"]; E0 = t[0]["E_carried_in"]; Ein = sum(x["E_inj_new"] for x in t); Eesc = sum(x["E_esc"] for x in t)
        h["lightcurve"]["free_streaming"] = dict(W=sum(x["W"] for x in t), closure=(E0 + Ein - Eesc - sum(x["W"] for x in t) - t[-1]["E_carried_out"]) / (E0 + Ein))
    except FileNotFoundError:
        pass
    # ---- the P1-xkn light curve (F66), when present ----
    p = ROOT / "paper4/phase10_fontes/prod_record/p1xkn/summary.json"
    if p.exists():
        d, rel, s = load("paper4/phase10_fontes/prod_record/p1xkn/summary.json"); note(rel, s)
        r, rel, s = load("paper4/phase10_fontes/prod_record/p1xkn/run.json"); note(rel, s)
        legs = {}
        for leg, o in d["legs"].items():
            e = dict(t_peak_d=o["t_peak_d"], L_peak=o["L_peak"], E_rad=o["E_rad"], W_tot=o["W_tot"], E_inj=o["E_inj"],
                     dL_peak=o.get("dL_peak"), max_abs_dcolour=o.get("max_abs_dcolour"), f_capped_max=o["f_capped_max"],
                     nan_band_cells=o["nan_band_cells"])
            if "dm_vs_ref" in o:
                dmv = {b: np.array([x[b] for x in o["dm_vs_ref"]], float) for b in BANDS}
                e["median_dm"] = {b: (float(np.nanmedian(v)) if np.isfinite(v).any() else None) for b, v in dmv.items()}
                for a, b in (("i", "K"), ("z", "K"), ("J", "K")):
                    v = dmv[a] - dmv[b]
                    e[f"median_d{a}{b}"] = float(np.nanmedian(v)); e[f"max_abs_d{a}{b}"] = float(np.nanmax(np.abs(v)))
            T = r["tallies"][leg]
            e["events_mean_first"] = T[0]["events_mean"]; e["events_mean_last"] = T[-1]["events_mean"]
            e["W_over_esc_first"] = T[0]["W"] / T[0]["E_esc"]; e["E_core_returned"] = sum(x["E_core"] for x in T)
            legs[leg] = e
        T = r["tallies"]["R2"]
        h["lightcurve_p1xkn"] = dict(legs=legs, readings=d["readings"], n_slabs=len(d["t_grid"]) - 1, config=r["config"],
                                     zone=dict(x_ph_first=T[0]["x_ph"], T_ph_first=T[0]["T_ph"], n_zone_first=T[0]["n_zone"],
                                               x_ph_last=T[-1]["x_ph"], T_ph_last=T[-1]["T_ph"], n_zone_last=T[-1]["n_zone"],
                                               core_frac_first=T[0]["core_frac"], v_max_c=r["v_max_c"]))
        ctrl = {}
        for tag in ("conv5", "conv10", "free"):
            q = f"paper4/phase10_fontes/prod_record/p1xkn_{tag}/summary.json"
            if (ROOT / q).exists():
                c, rel, s = load(q); note(rel, s); o = c["legs"]["R2"]
                ctrl[tag] = dict(E_rad=o["E_rad"], W_tot=o["W_tot"], E_end=o["E_end"], closure=o["closure"], n_slabs=len(c["t_grid"]) - 1)
                if tag == "free":
                    ratio = np.array(o["L_esc"]) / np.array(o["L_xkn"])
                    ctrl[tag]["L_ratio_first"] = float(ratio[0]); ctrl[tag]["L_ratio_min_after"] = float(ratio[1:].min()); ctrl[tag]["L_ratio_max_after"] = float(ratio[1:].max())
        h["lightcurve_p1xkn"]["controls"] = ctrl
    return h


def build():
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    return dict(paper="IV", git_head=head, headline=headline())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
    fz = build()
    if a.check:
        old = json.loads(OUT.read_text())
        same = json.dumps(old["headline"], sort_keys=True, default=float) == json.dumps(fz["headline"], sort_keys=True, default=float)
        print("FROZEN.json headline", "matches" if same else "DIFFERS from", "the records")
        sys.exit(0 if same else 1)
    OUT.write_text(json.dumps(fz, indent=1, default=float))
    print(f"wrote {OUT} at {fz['git_head'][:10]} with {len(fz['headline']['sources'])} sources")


if __name__ == "__main__":
    main()
