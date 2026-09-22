#!/usr/bin/env python3
"""Paper B gate G1: the metrics, epsilon*, the readings B1-B3 with the gray
conditions checked first, exactly as preregistered in paperB/prl_gate.md.

    .venv/bin/python paperB/gate1/analyse.py                # all gate1_*.json present
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper4/phase3_legs", ROOT / "paper3/phase12_grid"):
    sys.path.insert(0, str(p))
from sobolev.photometry import COLORS                      # noqa: E402
import verdict as V                                         # noqa: E402

# ---- the preregistration, as numbers (prl_gate.md is the text of record) ----
PREREG = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, n=300_000, seeds=[1, 2, 3],
              ng_grid=[2, 4, 8, 16, 32], ng_fine=128, eps_grid=[round(0.05 * k, 2) for k in range(21)],
              dm_max=0.10, dcolour_max=0.10, eps_substantial=0.20, eps_ratio_green=3.0, eps_ratio_red=1.5,
              gray_min_live=2, gray_seed_std=0.05, gray_identity=1e-10, gray_kernel_energy=1e-12,
              gray_trapped_frac=0.01, gray_fallback_frac=0.01, green_ng=8, max_ng=32)
IONS = ("57LaII", "58CeII", "60NdII")


def check_prereg(row, strict=True):
    """The record must be the preregistered experiment."""
    bad = []
    for k in ("state", "shell", "n", "ng_grid", "ng_fine", "eps_grid"):
        if row.get(k) != PREREG[k]:
            bad.append(f"{k}: {row.get(k)!r} != {PREREG[k]!r}")
    if row.get("seeds") != PREREG["seeds"]:
        bad.append(f"seeds: {row.get('seeds')} != {PREREG['seeds']}")
    if bad and strict:
        raise ValueError("record is not the preregistered experiment: " + "; ".join(bad))
    return bad


# ---- metrics ----
def block_expand(R_coarse, edges_coarse, edges_fine):
    """R_coarse (n_c x n_c) expanded onto the fine edges: every fine group
    inherits its coarse group's row, and each coarse column's mass is spread
    over the fine groups it contains in proportion to their width in log nu
    (the finest statement the coarse matrix makes)."""
    ec, ef = np.asarray(edges_coarse), np.asarray(edges_fine)
    cf = np.sqrt(ef[1:] * ef[:-1])
    gi = np.clip(np.searchsorted(ec, cf, side="right") - 1, 0, ec.size - 2)
    n_f = cf.size
    out = np.zeros((n_f, n_f))
    width = np.log(ef[1:] / ef[:-1])
    for j in range(ec.size - 1):
        cols = np.flatnonzero(gi == j)
        share = width[cols] / width[cols].sum() if cols.size else None
        for i in range(n_f):
            if cols.size:
                out[i, cols] = R_coarse[gi[i], j] * share
    return out


def m_event(kern_coarse, kern_fine):
    """Usage-weighted L1 between the coarse energy matrix (block-expanded)
    and the fine one, both from the same events."""
    Rf = np.asarray(kern_fine["R"]); cf = np.asarray(kern_fine["counts"], float)
    Rc = block_expand(np.asarray(kern_coarse["R"]), kern_coarse["edges"], kern_fine["edges"])
    rows = cf > 0
    d = np.abs(Rc[rows] - Rf[rows]).sum(axis=1)
    return float((d * cf[rows]).sum() / cf[rows].sum())


def m_sed(leg, ref):
    a, b = np.asarray(leg["L_nu"], float), np.asarray(ref["L_nu"], float)
    return float(np.abs(a - b).sum() / b.sum())


def m_band(leg, ref, live):
    dm = [leg["mags"][b] - ref["mags"][b] for b in live]
    return dict(max=float(np.max(np.abs(dm))), mean=float(np.mean(np.abs(dm))), per_band={b: float(d) for b, d in zip(live, dm)})


def m_colour(leg, ref, live):
    out = {}
    for a, b in COLORS:
        if a in live and b in live:
            out[f"{a}-{b}"] = float((leg["mags"][a] - leg["mags"][b]) - (ref["mags"][a] - ref["mags"][b]))
    return dict(max=float(max(abs(v) for v in out.values())) if out else float("nan"), per_colour=out)


def seed_noise(row, live):
    """R2's seed scatter, and the noise on a difference of two legs (sqrt 2 x)."""
    s = row["legs"]["R2"]["mags_seed_std"]
    return {b: float(s[b]) for b in live}


def metrics(row):
    ref = row["legs"]["R2"]
    live = V.live_bands(row, "R2")
    noise = seed_noise(row, live)
    fine = row["kernels"][f"K{row['ng_fine']}"]
    out = dict(live_bands=live, seed_std_R2=noise, legs={})
    n_int_ref = max(ref.get("n_interactions", 0), 1)
    for tag, leg in row["legs"].items():
        m = dict(mode=leg["mode"], identity=abs(leg["energy"]["identity_residual"]),
                 trapped_frac=leg["n_trapped"] / (row["n"] * len(row["seeds"])),
                 fallback_frac=leg.get("n_coherent_fallback", 0) / max(leg.get("n_interactions", 0), 1),
                 events_per_packet=leg["events_per_packet"], t_wall=leg["t_wall"])
        if tag != "R2":
            m["band"] = m_band(leg, ref, live); m["sed"] = m_sed(leg, ref); m["colour"] = m_colour(leg, ref, live)
        if tag in row["kernels"]:
            k = row["kernels"][tag]
            m.update(ng=k["ng"], n_params=k["ng"] ** 2, table_kb=k.get("table_kb"), kernel_energy=k["validate_energy"],
                     empty_rows=k["empty_rows"], event=m_event(k, fine) if tag != f"K{row['ng_fine']}" else 0.0)
        if "eps" in leg:
            m["eps"] = leg["eps"]
        out["legs"][tag] = m
    return out


def eps_star(M, row):
    """The grid value minimising mean |dm| over the live bands; flags an
    edge minimum whose neighbour is within the seed noise."""
    grid = sorted((m["eps"], tag) for tag, m in M["legs"].items() if "eps" in m)
    curve = [(e, M["legs"][t]["band"]["mean"], M["legs"][t]["band"]["max"]) for e, t in grid]
    k = int(np.argmin([c[1] for c in curve]))
    e_star, tag = grid[k]
    noise = float(np.sqrt(2.0) * np.mean(list(M["seed_std_R2"].values()))) if M["seed_std_R2"] else float("nan")
    edge = k in (0, len(curve) - 1)
    neighbour = curve[1][1] if k == 0 else curve[-2][1] if k == len(curve) - 1 else None
    edge_gray = bool(edge and neighbour is not None and abs(neighbour - curve[k][1]) <= noise)
    return dict(eps=e_star, tag=tag, mean_dm=curve[k][1], max_dm=curve[k][2], sed=M["legs"][tag]["sed"],
                colour_max=M["legs"][tag]["colour"]["max"], interior=not edge, edge_gray=edge_gray,
                curve=[dict(eps=e, mean_dm=a, max_dm=b) for e, a, b in curve], noise_on_dm=noise)


def gray_checks(row, M):
    P = PREREG; fired = []
    if len(M["live_bands"]) < P["gray_min_live"]:
        fired.append(f"1: {len(M['live_bands'])} live bands")
    bad = [b for b, s in M["seed_std_R2"].items() if s > P["gray_seed_std"]]
    if bad:
        fired.append(f"2: R2 seed scatter > {P['gray_seed_std']} in {bad}")
    for tag, m in M["legs"].items():
        if m["identity"] > P["gray_identity"]:
            fired.append(f"3: {tag} identity residual {m['identity']:.1e}")
        if m.get("kernel_energy", 0.0) > P["gray_kernel_energy"]:
            fired.append(f"3: {tag} kernel energy {m['kernel_energy']:.1e}")
        if m["trapped_frac"] > P["gray_trapped_frac"]:
            fired.append(f"4: {tag} chain-capped {m['trapped_frac']:.3f}")
        if m["fallback_frac"] > P["gray_fallback_frac"]:
            fired.append(f"4: {tag} empty-row fallback {m['fallback_frac']:.3f}")
    return fired


def readings(records):
    """B1-B3 over the ions; gray first. `records` = {ion: row}."""
    P = PREREG; per_ion = {}
    for ion, row in records.items():
        M = metrics(row); es = eps_star(M, row); gray = gray_checks(row, M)
        if es["edge_gray"]:
            gray.append(f"5: eps* at the grid edge {es['eps']} with the neighbour within the noise")
        ngs = sorted(m["ng"] for t, m in M["legs"].items() if "ng" in m and "band" in m)
        by_ng = {m["ng"]: m for t, m in M["legs"].items() if "ng" in m and "band" in m}
        ng_star = next((n for n in ngs if by_ng[n]["band"]["max"] <= P["dm_max"] and by_ng[n]["colour"]["max"] <= P["dcolour_max"]), None)
        e_R = by_ng[ng_star]["band"]["max"] if ng_star is not None else by_ng[max(ngs)]["band"]["max"]
        noise = es["noise_on_dm"]
        mono = {}
        for key, get in (("band", lambda m: m["band"]["max"]), ("sed", lambda m: m["sed"]), ("event", lambda m: m["event"])):
            vals = [get(by_ng[n]) for n in ngs]
            tol = noise if key == "band" else 0.0
            rises = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
            mono[key] = dict(values=vals, monotone=all(r <= tol for r in rises), max_rise=float(max(rises)) if rises else 0.0)
        per_ion[ion] = dict(gray=gray, live_bands=M["live_bands"], seed_std_R2=M["seed_std_R2"], ng_star=ng_star,
                            e_eps=es["max_dm"], e_R=e_R, eps_star=es, monotone=mono,
                            table=[dict(ng=n, n_params=n * n, table_kb=by_ng[n].get("table_kb"), band_max=by_ng[n]["band"]["max"],
                                        band_mean=by_ng[n]["band"]["mean"], sed=by_ng[n]["sed"], event=by_ng[n]["event"],
                                        colour_max=by_ng[n]["colour"]["max"], events_per_packet=by_ng[n]["events_per_packet"]) for n in ngs],
                            R2=dict(events_per_packet=M["legs"]["R2"]["events_per_packet"], t_wall=M["legs"]["R2"]["t_wall"]))
    live_ions = [i for i, r in per_ion.items() if not r["gray"]]
    out = dict(prereg=PREREG, per_ion=per_ion, ions_read=live_ions, ions_gray=[i for i in per_ion if per_ion[i]["gray"]])
    if len(live_ions) < 2:
        out.update(B1="GRAY", B2="GRAY", B3="GRAY", decision="GRAY", reason=f"only {len(live_ions)} ion(s) readable")
        return out
    # B1
    ns = [per_ion[i]["ng_star"] for i in live_ions]
    defined = [n for n in ns if n is not None]
    small = sum(1 for n in defined if n <= P["green_ng"])
    if small >= 2 and len(defined) == len(ns):
        B1 = "GREEN"
    elif len(defined) >= 2:
        B1 = "YELLOW"
    else:
        B1 = "RED"
    # B2
    verdicts = []
    for i in live_ions:
        r = per_ion[i]
        if r["e_eps"] >= P["eps_substantial"] and r["e_eps"] >= P["eps_ratio_green"] * r["e_R"]:
            verdicts.append("GREEN")
        elif r["e_eps"] <= P["eps_ratio_red"] * r["e_R"]:
            verdicts.append("RED")
        else:
            verdicts.append("YELLOW")
        r["B2"] = verdicts[-1]
    B2 = "RED" if verdicts.count("RED") >= 2 or (verdicts.count("RED") == 1 and len(verdicts) == 2) else "GREEN" if verdicts.count("GREEN") >= 2 else "YELLOW"
    # B3
    b3 = "GREEN"
    for i in live_ions:
        mo = per_ion[i]["monotone"]
        if not mo["band"]["monotone"] and mo["band"]["max_rise"] > per_ion[i]["eps_star"]["noise_on_dm"]:
            b3 = "RED"; break
        if not all(mo[k]["monotone"] for k in mo):
            b3 = "YELLOW"
    decision = "CONTINUE" if (B1 == "GREEN" and B2 == "GREEN" and b3 != "RED") else "STOP" if (B1 == "RED" or B2 == "RED") else "YELLOW"
    out.update(B1=B1, B2=B2, B3=b3, decision=decision)
    return out


def main():
    records = {}
    for ion in IONS:
        p = HERE / f"gate1_{ion}.json"
        if p.exists():
            row = json.loads(p.read_text()); check_prereg(row); records[ion] = row
    if not records:
        print("no records"); return 1
    out = readings(records)
    (HERE / "gate1_verdict.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for ion, r in out["per_ion"].items():
        print(f"\n{ion}: live {r['live_bands']}, R2 seed std {{{', '.join(f'{b} {s:.3f}' for b, s in r['seed_std_R2'].items())}}}"
              + (f"  GRAY: {r['gray']}" if r["gray"] else ""))
        print("   N_g  params  band_max  band_mean   sed    event  colour_max  ev/pkt")
        for t in r["table"]:
            print(f"  {t['ng']:4d}  {t['n_params']:6d}   {t['band_max']:.3f}     {t['band_mean']:.3f}   {t['sed']:.3f}  {t['event']:.3f}    {t['colour_max']:.3f}   {t['events_per_packet']:.1f}")
        es = r["eps_star"]
        print(f"  eps* = {es['eps']:.2f}: band_max {es['max_dm']:.3f} mean {es['mean_dm']:.3f} sed {es['sed']:.3f} colour {es['colour_max']:.3f}"
              f" (interior {es['interior']}); N_g* = {r['ng_star']}; B2 {r.get('B2', '-')}")
    print(f"\nB1 {out['B1']}  B2 {out['B2']}  B3 {out['B3']}  ->  {out['decision']}" + (f"  ({out['reason']})" if "reason" in out else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
