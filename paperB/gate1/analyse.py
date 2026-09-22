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
PREREG = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, n=300_000, seeds=[1, 2, 3], build_seeds=[101, 102, 103],
              ng_grid=[2, 4, 8, 16, 32], ng_fine=128, eps_grid=[round(0.05 * k, 2) for k in range(21)],
              dm_max=0.10, dcolour_max=0.10, eps_substantial=0.20, eps_ratio_green=3.0, eps_ratio_red=1.5,
              live_frac=0.01, live_seed_std=0.05,
              gray_min_live=2, gray_dropped_frac=0.01, gray_identity=1e-10, gray_kernel_energy=1e-12,
              gray_trapped_frac=0.01, gray_fallback_frac=0.01, green_ng=8, max_ng=32, b3_lo=8, b3_hi=32)
IONS = ("57LaII", "58CeII", "60NdII")


def check_prereg(row, strict=True):
    """The record must be the preregistered experiment."""
    bad = []
    for k in ("state", "shell", "ng_grid", "ng_fine", "eps_grid", "build_seeds"):
        if row.get(k) != PREREG[k]:
            bad.append(f"{k}: {row.get(k)!r} != {PREREG[k]!r}")
    if row.get("n", 0) < PREREG["n"]:                     # the packet count may only be raised (a gray remedy), never lowered
        bad.append(f"n: {row.get('n')} < {PREREG['n']}")
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
    """The expected total-variation loss of the outgoing redistribution
    distribution for an energy-weighted incoming interaction:
    m = sum_i W_i d_i / sum_i W_i, d_i = 1/2 sum_j |R^N_ij - R^128_ij|,
    W_i the energy absorbed in fine row i; the coarse matrix block-expanded
    onto the fine edges, the fine matrix from INDEPENDENT events (the
    evaluation seeds). In [0, 1]."""
    Rf = np.asarray(kern_fine["R"]); W = np.asarray(kern_fine["E_in"], float)
    Rc = block_expand(np.asarray(kern_coarse["R"]), kern_coarse["edges"], kern_fine["edges"])
    rows = W > 0
    d = 0.5 * np.abs(Rc[rows] - Rf[rows]).sum(axis=1)
    return float((d * W[rows]).sum() / W[rows].sum())


def m_sed(leg, ref, dnu):
    """Integrated: sum_b |L_nu - L_nu^R2| dnu_b / sum_b L_nu^R2 dnu_b."""
    a, b = np.asarray(leg["L_nu"], float), np.asarray(ref["L_nu"], float)
    return float((np.abs(a - b) * dnu).sum() / (b * dnu).sum())


def m_band(leg, ref, live):
    dm = [leg["mags"][b] - ref["mags"][b] for b in live]
    if not dm:                                             # no live band: the ion is gray (condition 1); nothing to read
        return dict(max=float("nan"), mean=float("nan"), per_band={})
    return dict(max=float(np.max(np.abs(dm))), mean=float(np.mean(np.abs(dm))), per_band={b: float(d) for b, d in zip(live, dm)})


def m_colour(leg, ref, live):
    out = {}
    for a, b in COLORS:
        if a in live and b in live:
            out[f"{a}-{b}"] = float((leg["mags"][a] - leg["mags"][b]) - (ref["mags"][a] - ref["mags"][b]))
    return dict(max=float(max(abs(v) for v in out.values())) if out else float("nan"), per_colour=out)


def live_bands(row, leg="R2"):
    """G1's rule: finite magnitude, >= 1 % of the window luminosity, and an
    acceptable Monte Carlo precision (seed scatter <= 0.05 mag). The 40 Mpc
    detectability of Paper IV is reported separately, not applied."""
    frac = V.band_fractions(row, leg); o = row["legs"][leg]
    live, dropped = [], []
    for b in "grizJHK":
        if not np.isfinite(o["mags"].get(b, np.nan)) or frac[b] < PREREG["live_frac"]:
            continue
        if o["mags_seed_std"].get(b, np.nan) > PREREG["live_seed_std"]:
            dropped.append(b)
        else:
            live.append(b)
    detectable = [b for b in live if o["mags"][b] <= V.MAG_LIMIT[b]]
    return live, dropped, detectable


def seed_noise(row, live):
    """R2's seed scatter, and the noise on a difference of two legs (sqrt 2 x)."""
    s = row["legs"]["R2"]["mags_seed_std"]
    return {b: float(s[b]) for b in live}


def metrics(row):
    ref = row["legs"]["R2"]
    live, dropped, detectable = live_bands(row, "R2")
    noise = seed_noise(row, live)
    fine = row["kernels"][f"K{row['ng_fine']}"]                         # from the evaluation seeds: independent of every kernel
    fine_build = row["kernels"].get(f"K{row['ng_fine']}build")
    edges = V.nu_edges(*row["lam_window"], row["n_spec"]); dnu = np.diff(edges)
    out = dict(live_bands=live, dropped_for_precision=dropped, detectable_40mpc=detectable, seed_std_R2=noise, legs={},
               fine_in_vs_out_of_sample=m_event(fine_build, fine) if fine_build else None)
    for tag, leg in row["legs"].items():
        m = dict(mode=leg["mode"], identity=abs(leg["energy"]["identity_residual"]),
                 trapped_frac=leg["n_trapped"] / (row["n"] * len(leg.get("seeds", row["seeds"]))),
                 fallback_frac=leg.get("n_coherent_fallback", 0) / max(leg.get("n_interactions", 0), 1),
                 events_per_packet=leg["events_per_packet"], t_wall=leg["t_wall"], seeds=leg.get("seeds", row["seeds"]))
        if tag not in ("R2", "R2build"):
            m["band"] = m_band(leg, ref, live); m["sed"] = m_sed(leg, ref, dnu); m["colour"] = m_colour(leg, ref, live)
        if tag in row["kernels"]:
            k = row["kernels"][tag]
            m.update(ng=k["ng"], n_matrix_dof=k["ng"] ** 2, n_exit_samples=k.get("n_exit_samples"), table_kb=k.get("table_kb"),
                     kernel_energy=k["validate_energy"], empty_rows=k["empty_rows"], kernel_source=k["source"],
                     event=m_event(k, fine) if not tag.startswith(f"K{row['ng_fine']}") else 0.0)
        if "eps" in leg:
            m["eps"] = leg["eps"]
        out["legs"][tag] = m
    return out


def eps_star(M, row):
    """The grid value minimising mean |dm| over the live bands; flags an
    edge minimum whose neighbour is within the seed noise."""
    grid = sorted((m["eps"], tag) for tag, m in M["legs"].items() if "eps" in m)
    curve = [(e, M["legs"][t]["band"]["mean"], M["legs"][t]["band"]["max"]) for e, t in grid]
    k = int(np.nanargmin([c[1] for c in curve])) if any(np.isfinite(c[1]) for c in curve) else 0
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
    if M["dropped_for_precision"]:
        fired.append(f"2: band(s) {M['dropped_for_precision']} carry >= 1 % of L_bol but exceed the precision rule (seed scatter > {P['live_seed_std']}): raise the packet count")
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
        conv = {}
        for key, get in (("band", lambda m: m["band"]["max"]), ("sed", lambda m: m["sed"]), ("event", lambda m: m["event"])):
            vals = {n: get(by_ng[n]) for n in ngs}
            tol = noise if key == "band" else 0.0
            v = [vals[n] for n in ngs]; rises = [v[i + 1] - v[i] for i in range(len(v) - 1)]
            lo, hi = P["b3_lo"], P["b3_hi"]
            conv[key] = dict(values=vals, local_non_monotone=any(r > tol for r in rises), max_rise=float(max(rises)) if rises else 0.0,
                             degrades_lo_hi=bool(lo in vals and hi in vals and vals[hi] > vals[lo] + tol),
                             converges_lo_hi=bool(lo in vals and hi in vals and vals[hi] <= vals[lo] + tol))
        per_ion[ion] = dict(gray=gray, live_bands=M["live_bands"], dropped_for_precision=M["dropped_for_precision"],
                            detectable_40mpc=M["detectable_40mpc"], seed_std_R2=M["seed_std_R2"], ng_star=ng_star,
                            e_eps=es["max_dm"], e_R=e_R, eps_star=es, convergence=conv,
                            fine_in_vs_out_of_sample=M["fine_in_vs_out_of_sample"],
                            table=[dict(ng=n, n_matrix_dof=n * n, n_exit_samples=by_ng[n].get("n_exit_samples"), table_kb=by_ng[n].get("table_kb"),
                                        band_max=by_ng[n]["band"]["max"], band_mean=by_ng[n]["band"]["mean"], sed=by_ng[n]["sed"],
                                        event=by_ng[n]["event"], colour_max=by_ng[n]["colour"]["max"], fallback_frac=by_ng[n]["fallback_frac"],
                                        events_per_packet=by_ng[n]["events_per_packet"]) for n in ngs],
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
    # B3 (revised 2026-09-22): the event-level metric is the structural diagnostic;
    # Red only for a representation pathology: m_event not converging from 8 to 32
    # groups, or a persistent worsening of BOTH observables from 8 to 32 beyond the
    # noise; an innocent local wiggle in an observable is Yellow.
    b3 = "GREEN"
    for i in live_ions:
        cv = per_ion[i]["convergence"]
        if cv["event"]["degrades_lo_hi"] or (cv["band"]["degrades_lo_hi"] and cv["sed"]["degrades_lo_hi"]):
            b3 = "RED"; break
        if cv["band"]["local_non_monotone"] or cv["sed"]["local_non_monotone"] or not cv["event"]["converges_lo_hi"]:
            b3 = "YELLOW"
    decision = "CONTINUE" if (B1 == "GREEN" and B2 == "GREEN" and b3 != "RED") else "STOP" if (B1 == "RED" or B2 == "RED") else "YELLOW"
    out.update(B1=B1, B2=B2, B3=b3, decision=decision)
    return out


def main():
    records = {}
    for ion in IONS:
        # the record at the highest packet count wins: a gray remedy reruns the whole grid at more packets
        cands = sorted(HERE.glob(f"gate1_{ion}*.json"), key=lambda q: json.loads(q.read_text()).get("n", 0))
        if cands:
            row = json.loads(cands[-1].read_text()); check_prereg(row); records[ion] = row
            print(f"{ion}: {cands[-1].name} (n = {row['n']:,})")
    if not records:
        print("no records"); return 1
    out = readings(records)
    (HERE / "gate1_verdict.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for ion, r in out["per_ion"].items():
        print(f"\n{ion}: live {r['live_bands']} (dropped {r['dropped_for_precision']}, detectable at 40 Mpc {r['detectable_40mpc']}), "
              f"R2 seed std {{{', '.join(f'{b} {s:.3f}' for b, s in r['seed_std_R2'].items())}}}, fine matrix in- vs out-of-sample {r['fine_in_vs_out_of_sample']:.4f}"
              + (f"  GRAY: {r['gray']}" if r["gray"] else ""))
        print("   N_g  N_g^2  exit_samples  kB    band_max  band_mean   sed    event  colour_max  fallback  ev/pkt")
        for t in r["table"]:
            print(f"  {t['ng']:4d}  {t['n_matrix_dof']:5d}  {t['n_exit_samples']:11d}  {t['table_kb']:6.1f}  {t['band_max']:.3f}     {t['band_mean']:.3f}   {t['sed']:.3f}  {t['event']:.4f}   {t['colour_max']:.3f}     {t['fallback_frac']:.4f}   {t['events_per_packet']:.1f}")
        es = r["eps_star"]
        print(f"  eps* = {es['eps']:.2f}: band_max {es['max_dm']:.3f} mean {es['mean_dm']:.3f} sed {es['sed']:.3f} colour {es['colour_max']:.3f}"
              f" (interior {es['interior']}); N_g* = {r['ng_star']}; B2 {r.get('B2', '-')}")
    print(f"\nB1 {out['B1']}  B2 {out['B2']}  B3 {out['B3']}  ->  {out['decision']}" + (f"  ({out['reason']})" if "reason" in out else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
