#!/usr/bin/env python3
"""G3U -- the paired-seed uncertainty of the near-threshold G3 readings
(paperB/prl_gate.md "G3U"). Two parts:

  affected()   the mechanically defined set: every decision leg of the
               decisive ions' G3 states and of the blend records whose max
               band OR max colour error lies within the window
               [0.07, 0.13] mag of the 0.10 threshold -- derived from the
               frozen records, never typed;
  read()       per rerun record and affected leg, the paired per-seed
               difference closure minus reference on the frozen live
               bands, the 12-seed point estimate of the joint decision
               statistic E_joint = max(max_b |dm_b|, max_c |dcolour_c|), its
               paired bootstrap (seed indices resampled, B = 10,000), the
               68 % and 95 % intervals, P(E_joint <= 0.10), and whether the
               reading is DECIDED (the 95 % interval excludes the threshold)
               or WITHIN NOISE. G3's frozen letters are not re-read.

    .venv/bin/python paperB/gate3u/analyse.py --affected     # the set, from the frozen G3 records
    .venv/bin/python paperB/gate3u/analyse.py                # g3u_*.json -> g3u_verdict.json; --markdown
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
A = _ilu.module_from_spec(_spec); _spec.loader.exec_module(A)
_spec3 = _ilu.spec_from_file_location("paperB_gate3_analyse", ROOT / "paperB/gate3/analyse.py")
G3 = _ilu.module_from_spec(_spec3); _spec3.loader.exec_module(G3)
from sobolev.photometry import COLORS                                  # noqa: E402

GATE3 = ROOT / "paperB/gate3"
PREREG = dict(window=[0.07, 0.13], dm_max=0.10, dcolour_max=0.10, seeds=list(range(1, 13)), build_seeds=[101, 102, 103],
              n_boot=10_000, boot_seed=0, ci68=[16.0, 84.0], ci95=[2.5, 97.5], operator_match=1e-10, determinism=1e-6,
              decisive=["58CeII", "60NdII"], ng_t=16,
              state_legs=["Afix_ng16", "Arec_ng16", "Arec_ng32", "Aint_ng16", "AintM_ng16"],
              blend_families=["Adirect", "Amix", "Arec"])
NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}


# ---- the frozen G3 records, exactly as gate3/analyse.py selects them ----
def frozen_records():
    """{name: row} -- the highest-packet-count record per (axis, label, ion), plus the two blend records."""
    best = {}
    for p in sorted(GATE3.glob("gate3_*_*.json")):
        if p.name.startswith(("gate3_partb", "gate3_partc", "gate3_verdict", "gate3_support")):
            continue
        row = json.loads(p.read_text())
        key = (row["axis"], row["label"], row["ion"])
        if key not in best or row["n"] > best[key][1]["n"]:
            best[key] = (p.name, row)
    out = {name: row for _, (name, row) in sorted(best.items())}
    for name in ("gate3_partb_blend3.json", "gate3_partc_p1blend.json"):
        if (GATE3 / name).exists():
            out[name] = json.loads((GATE3 / name).read_text())
    return out


def candidate_legs(row):
    """The decision legs of a record: for a state the transfer, fresh R16 and R32 and the two interpolants that exist;
    for a blend every valid leg of the families searched for a K*."""
    if "part" in row:
        M = A.metrics(row)
        return [t for t, m in M["legs"].items() if "band" in m and "ng" in m and t.split("_ng")[0] in PREREG["blend_families"]
                and not G3.leg_invalid(m)]
    return [t for t in PREREG["state_legs"] if t in row["legs"]]


def affected(records=None, window=None):
    """[(record name, leg, band max, colour max)] -- the near-threshold decision legs of the decisive ions and the blends."""
    lo, hi = window or PREREG["window"]
    records = records or frozen_records()
    out = []
    for name, row in records.items():
        if "part" not in row and row["ion"] not in PREREG["decisive"]:
            continue
        if "part" not in row and row["axis"] == "ref":
            continue
        M = A.metrics(row)
        for tag in candidate_legs(row):
            m = M["legs"][tag]
            b, c = m["band"]["max"], m["colour"]["max"]
            if lo <= b <= hi or lo <= c <= hi:
                out.append(dict(record=name, leg=tag, band=b, colour=c, ion=row.get("ion"), axis=row.get("axis"), label=row.get("label"),
                                part=row.get("part"), n=row["n"]))
    return out


# ---- the paired statistic ----
def paired_diffs(leg, ref, live):
    """Per seed: the band differences closure - reference on the live bands and the colour differences on the live pairs."""
    S = min(len(leg["mags_per_seed"]), len(ref["mags_per_seed"]))
    bands = np.array([[leg["mags_per_seed"][s][b] - ref["mags_per_seed"][s][b] for b in live] for s in range(S)])
    pairs = [(a, b) for a, b in COLORS if a in live and b in live]
    cols = np.array([[(leg["mags_per_seed"][s][a] - leg["mags_per_seed"][s][b]) - (ref["mags_per_seed"][s][a] - ref["mags_per_seed"][s][b])
                      for a, b in pairs] for s in range(S)])
    return bands, cols, pairs


def e_joint_of(bands, cols):
    """max over bands and colours of |mean over seeds|."""
    eb = np.abs(bands.mean(axis=0)).max() if bands.size else 0.0
    ec = np.abs(cols.mean(axis=0)).max() if cols.size else 0.0
    return float(max(eb, ec))


def bootstrap(bands, cols, n_boot=None, seed=None):
    """The paired bootstrap of E_joint over the seed index (every replicate recomputes the max over bands and colours)."""
    P = PREREG; n_boot = n_boot or P["n_boot"]; seed = P["boot_seed"] if seed is None else seed
    S = bands.shape[0]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, S, (n_boot, S))
    mb = bands[idx].mean(axis=1)                       # (n_boot, n_bands)
    eb = np.abs(mb).max(axis=1) if bands.shape[1] else np.zeros(n_boot)
    if cols.shape[1]:
        mc = cols[idx].mean(axis=1); ec = np.abs(mc).max(axis=1)
    else:
        ec = np.zeros(n_boot)
    e = np.maximum(eb, ec)
    return dict(n_boot=int(n_boot), seed=int(seed), ci68=[float(x) for x in np.percentile(e, P["ci68"])],
                ci95=[float(x) for x in np.percentile(e, P["ci95"])], p_pass=float(np.mean(e <= P["dm_max"])),
                mean=float(e.mean()), sd=float(e.std(ddof=1)))


def read_case(row, frozen, leg_tag):
    """One rerun record and one affected leg against its frozen record."""
    P = PREREG
    live = A.live_bands(frozen, "R2")[0]                  # the state's frozen live-band set
    leg, ref = row["legs"][leg_tag], row["legs"]["R2"]
    bands, cols, pairs = paired_diffs(leg, ref, live)
    S = bands.shape[0]
    e12 = e_joint_of(bands, cols)
    e3 = e_joint_of(bands[:3], cols[:3])
    bs = bootstrap(bands, cols)
    gray = []
    if S < len(P["seeds"]):
        gray.append(f"{S} of {len(P['seeds'])} seeds")
    # determinism: the per-seed scatter over the first three seeds must be the frozen record's
    for t in (leg_tag, "R2"):
        if t in frozen["legs"] and len(row["legs"][t]["mags_per_seed"]) >= 3:
            sd = {b: float(np.std([m[b] for m in row["legs"][t]["mags_per_seed"][:3]], ddof=1)) for b in live}
            dev = max(abs(sd[b] - frozen["legs"][t]["mags_seed_std"][b]) for b in live)
            if dev > P["determinism"]:
                gray.append(f"{t}: seeds 1-3 scatter differs from the frozen record by {dev:.1e} mag")
    om = row.get("operator_match", {}).get(leg_tag)
    if om is None or om > P["operator_match"]:
        gray.append(f"{leg_tag}: operator differs from the frozen record by {om} (or not checked)")
    fm = A.metrics(frozen)["legs"][leg_tag]
    lo, hi = bs["ci95"]
    status = "GRAY" if gray else "DECIDED_PASS" if hi < P["dm_max"] else "DECIDED_FAIL" if lo > P["dm_max"] else "WITHIN_NOISE"
    return dict(leg=leg_tag, live_bands=live, colours=[f"{a}-{b}" for a, b in pairs], n_seeds=S, gray=gray,
                frozen=dict(band=fm["band"]["max"], colour=fm["colour"]["max"], e_joint=float(max(fm["band"]["max"], fm["colour"]["max"])),
                            passes=G3.passes(fm)),
                e_joint_3seeds=e3, e_joint=e12, band_max=float(np.abs(bands.mean(axis=0)).max()) if bands.size else None,
                colour_max=float(np.abs(cols.mean(axis=0)).max()) if cols.size else None,
                per_band_mean={b: float(v) for b, v in zip(live, bands.mean(axis=0))},
                per_band_sd={b: float(v) for b, v in zip(live, bands.std(axis=0, ddof=1))} if S > 1 else None,
                bootstrap=bs, status=status, operator_match=om)


def read_record(row, frozen):
    return [read_case(row, frozen, t) for t in row["affected_legs"]]


def main():
    if "--affected" in sys.argv:
        aff = affected()
        for a in aff:
            print(f"{a['record']:40s} {a['leg']:12s} band {a['band']:.4f} colour {a['colour']:.4f}")
        print(f"{len(aff)} affected legs in {len({a['record'] for a in aff})} records")
        return aff
    frozen = frozen_records()
    out = dict(prereg=PREREG, affected=affected(frozen), records={})
    for p in sorted(HERE.glob("g3u_*.json")):
        if p.name == "g3u_verdict.json":
            continue
        row = json.loads(p.read_text())
        fr = frozen[row["frozen_record"]]
        out["records"][p.name] = dict(frozen_record=row["frozen_record"], ion=row.get("ion"), label=row.get("label"), part=row.get("part"),
                                      seeds=row["seeds"], n=row["n"], cases=read_record(row, fr))
    cases = [c for r in out["records"].values() for c in r["cases"]]
    out["summary"] = dict(n_cases=len(cases), n_affected=len(out["affected"]),
                          status={s: sum(1 for c in cases if c["status"] == s) for s in ("DECIDED_PASS", "DECIDED_FAIL", "WITHIN_NOISE", "GRAY")},
                          flips=[dict(record=n, leg=c["leg"], frozen_passes=c["frozen"]["passes"], e_joint=c["e_joint"])
                                 for n, r in out["records"].items() for c in r["cases"] if (c["e_joint"] <= PREREG["dm_max"]) != c["frozen"]["passes"]])
    (HERE / "g3u_verdict.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for n, r in out["records"].items():
        for c in r["cases"]:
            b = c["bootstrap"]
            print(f"{n:34s} {c['leg']:11s} frozen {c['frozen']['e_joint']:.3f} -> {c['n_seeds']}-seed {c['e_joint']:.3f} "
                  f"[{b['ci95'][0]:.3f}, {b['ci95'][1]:.3f}] P(pass) {b['p_pass']:.2f}  {c['status']}" + (f"  gray: {c['gray']}" if c["gray"] else ""))
    print(f"\n{out['summary']}")
    return out


if __name__ == "__main__" and "--markdown" not in sys.argv:
    main()


def markdown(out=None):
    out = out or json.loads((HERE / "g3u_verdict.json").read_text())
    L = ["| record | leg | frozen (3 seeds) band / colour | E_joint, 12 seeds | 68 % | 95 % | P(E_joint ≤ 0.10) | status |",
         "|---|---|---|---|---|---|---|---|"]
    for n, r in out["records"].items():
        for c in r["cases"]:
            b = c["bootstrap"]
            L.append(f"| {n.replace('g3u_', '').replace('.json', '')} | {c['leg']} | {c['frozen']['band']:.3f} / {c['frozen']['colour']:.3f} | "
                     f"{c['e_joint']:.3f} | [{b['ci68'][0]:.3f}, {b['ci68'][1]:.3f}] | [{b['ci95'][0]:.3f}, {b['ci95'][1]:.3f}] | {b['p_pass']:.2f} | "
                     f"{c['status'].lower().replace('_', ' ')} |")
    s = out["summary"]
    L += ["", f"{s['n_cases']} cases of {s['n_affected']} affected legs: " + ", ".join(f"{k.lower().replace('_', ' ')} {v}" for k, v in s["status"].items())
          + (f"; point estimate on the other side of the threshold from the frozen reading: {len(s['flips'])}" if s["flips"] else "; no point estimate crosses the threshold") + "."]
    return "\n".join(L)


if __name__ == "__main__" and "--markdown" in sys.argv:
    print(markdown())
