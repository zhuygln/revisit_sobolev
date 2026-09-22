#!/usr/bin/env python3
"""Paper B, the exit-table audit (representation only; no transport is
scored, no physics changes). What `RedistributionKernel` stores per exit
and what that storage is bounded by.

Established before this script was written: `from_branching_mc` already
aggregates the exits by unique frequency (np.unique(nu_out) with the
energy weights summed), so a stored entry is one distinct exit line, never
a repeated event. This audit measures, per ion on the G1 build events:
the number of distinct exit lines and whether every one is an exact line
rest frequency; the discovery curve (distinct lines against events
drawn); the concentration of the exit energy across lines (how many carry
50 / 90 / 99 / 99.9 %); what a per-group truncation at fraction f would
keep (the numbers G2's H2 will transport); and the serialized size split
into matrix and tables, checked against G1's recorded K128build bytes.

    .venv/bin/python paperB/audit/exit_tables.py --ion 57LaII      # then 58CeII, 60NdII (Nd at 1e6: ~10 min, ~3 GB)
    .venv/bin/python paperB/audit/exit_tables.py --figure           # docs/figures/paperB/exit_table_audit
"""
import argparse
import io
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paperB/gate1", ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper4/phase2_energy",
          ROOT / "paper3/phase11_observables", ROOT / "paper3/phase12_grid"):
    sys.path.insert(0, str(p))

from redistribution import RedistributionKernel                    # noqa: E402

FRACTIONS = (0.1, 0.2, 0.5, 0.9, 0.99, 0.999)
G1_N = {"57LaII": 300_000, "58CeII": 300_000, "60NdII": 1_000_000}   # the packet count of each ion's final G1 record


# ---- the pure functions (tested on synthetic events) ----
def distinct_exits(nu_out, w_e):
    """(unique exit frequencies, energy per unique exit) from the events."""
    vals, inv = np.unique(np.asarray(nu_out, float), return_inverse=True)
    return vals, np.bincount(inv, weights=np.asarray(w_e, float), minlength=vals.size)


def discovery_curve(nu_out, rng, n_points=14, n_min=1e4):
    """Distinct exit lines against the number of events drawn, in a random
    order (seeded): whether the count saturates."""
    nu_out = np.asarray(nu_out, float)
    order = rng.permutation(nu_out.size)
    ms = np.unique(np.geomspace(min(n_min, nu_out.size), nu_out.size, n_points).astype(int))
    seen = np.zeros(ms.size, int)
    # cumulative distinct count via first-occurrence positions
    _, first = np.unique(nu_out[order], return_index=True)
    first = np.sort(first)
    for i, m in enumerate(ms):
        seen[i] = int(np.searchsorted(first, m, side="left"))
    return [dict(events=int(m), distinct=int(s)) for m, s in zip(ms, seen)]


def concentration(e_line, fractions=FRACTIONS):
    """How many lines (sorted by energy) carry each fraction of the total."""
    e = np.sort(np.asarray(e_line, float))[::-1]
    c = np.cumsum(e) / e.sum()
    return {f"{f:g}": int(np.searchsorted(c, f, side="left") + 1) for f in fractions}


def kept_per_group(vals, e_line, edges, fractions=FRACTIONS):
    """What a per-output-group truncation at fraction f keeps: for each f the
    total number of lines kept and the maximum over groups; every populated
    group keeps at least one line."""
    g = np.clip(np.searchsorted(edges, vals, side="right") - 1, 0, edges.size - 2)
    out = {f"{f:g}": dict(kept=0, max_in_group=0) for f in fractions}
    populated = 0
    for j in np.unique(g):
        m = g == j
        e = np.sort(e_line[m])[::-1]; c = np.cumsum(e) / e.sum(); populated += 1
        for f in fractions:
            k = int(np.searchsorted(c, f, side="left") + 1)
            out[f"{f:g}"]["kept"] += k; out[f"{f:g}"]["max_in_group"] = max(out[f"{f:g}"]["max_in_group"], k)
    return out, int(populated)


def size_split(kern):
    """(total bytes as saved, bytes without the discrete tables)."""
    b = io.BytesIO(); kern.save(b); total = b.getbuffer().nbytes
    bare = RedistributionKernel(kern.edges, kern.R, kern.N_cum, kern.q_dep, kern.sub_cum, kern.counts, kern.metadata)
    b2 = io.BytesIO(); bare.save(b2)
    return int(total), int(b2.getbuffer().nbytes)


# ---- the run ----
def audit(ion, n=None, build_seeds=(101, 102, 103), ng_fine=128, verbose=True):
    import run_gate1 as G
    import legs as L
    from forest_mc import run_mc
    import sobolev.photometry as phot
    n = G1_N[ion] if n is None else int(n)
    t0 = time.time()
    st, zone, atom, n_ion = G.build(ion)
    lo, hi = (float(x) for x in phot.nu_edges(*L.LAM_WIN, 1))
    ev = []
    for s in build_seeds:
        r = run_mc(atom, zone["r_core"], zone["r_out"], zone["t_exp"], lo, hi, n, "sobolev_dmacro", seed=s, t_core=zone["t_core"],
                   relativity="worldline", max_steps=L.MAX_STEPS, chain_max=L.CHAIN_MAX, chain_overflow="absorb", packets="energy",
                   launch_weight="energy", collect_events=True)
        ev.append(r["events"]); del r
    nu_in = np.concatenate([e[0] for e in ev]); nu_out = np.concatenate([e[1] for e in ev])
    w_in = np.concatenate([e[2] for e in ev]); w_out = np.concatenate([e[3] for e in ev]); del ev
    vals, e_line = distinct_exits(nu_out, w_out * nu_out)
    nu0 = np.unique(atom.nu0_all)
    exact = float(np.isin(vals, nu0).mean())
    k_lo, k_hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    kern = RedistributionKernel.from_branching_mc(nu_in, nu_out, w_in, ng_fine, nu_lo=k_lo, nu_hi=k_hi, w_out=w_out)
    assert kern.disc_vals.size == vals.size and np.array_equal(np.sort(kern.disc_vals), vals)
    total, bare = size_split(kern)
    # the G1 record's K128build, built from the same events, must be the same size
    rec = None
    for cand in sorted(G.HERE.glob(f"gate1_{ion}*.json")):
        row = json.loads(cand.read_text())
        if row["n"] == n and row["build_seeds"] == list(build_seeds):
            rec = row["kernels"].get(f"K{ng_fine}build")
    kept, populated = kept_per_group(vals, e_line, kern.edges)
    # the N_g* kernel's split (G1's verdict)
    ng_star = None
    vp = G.HERE / "gate1_verdict.json"
    if vp.exists():
        ng_star = json.loads(vp.read_text())["per_ion"].get(ion, {}).get("ng_star")
    star = None
    if ng_star:
        ks = RedistributionKernel.from_branching_mc(nu_in, nu_out, w_in, int(ng_star), nu_lo=k_lo, nu_hi=k_hi, w_out=w_out)
        t_s, b_s = size_split(ks)
        star = dict(ng=int(ng_star), bytes_total=t_s, bytes_matrix=b_s, bytes_tables=t_s - b_s, n_exit=int(ks.disc_vals.size))
    out = dict(ion=ion, n=n, build_seeds=list(build_seeds), n_events=int(nu_out.size), n_packets_total=n * len(build_seeds),
               n_lines=int(atom.n_lines_total), n_opacity=int(atom.n_opacity), n_distinct_exit_lines=int(vals.size),
               exact_line_fraction=exact, distinct_over_lines=vals.size / atom.n_lines_total,
               exit_energy_in_opacity_lines=float(e_line[np.isin(vals, atom.op_nu)].sum() / e_line.sum()),
               discovery=discovery_curve(nu_out, np.random.default_rng(0)),
               concentration=concentration(e_line), kept_per_group=kept, populated_output_groups=populated,
               K128=dict(bytes_total=total, bytes_matrix=bare, bytes_tables=total - bare, n_exit=int(vals.size),
                         bytes_per_exit=(total - bare) / vals.size, g1_record_bytes=rec["serialized_bytes"] if rec else None,
                         g1_record_n_exit=rec["n_exit_samples"] if rec else None),
               K_ng_star=star, t_wall=time.time() - t0)
    if verbose:
        print(f"{ion}: {out['n_events']} events -> {vals.size} distinct exit lines ({exact:.3f} exact rest frequencies) of {atom.n_lines_total} lines; "
              f"99 % of the exit energy in {out['concentration']['0.99']} lines; K128 {total/1024:.0f} kB of which tables {(total-bare)/1024:.0f} kB"
              + (f" (G1 record {rec['serialized_bytes']/1024:.0f} kB)" if rec else "") + f"; {out['t_wall']:.0f} s", flush=True)
    return out


def figure(records, out_dir=ROOT / "docs/figures/paperB"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COL = {"57LaII": "#0072B2", "58CeII": "#E69F00", "60NdII": "#D55E00"}
    NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ion, r in records.items():
        c = COL[ion]
        axes[0].plot([d["events"] for d in r["discovery"]], [d["distinct"] for d in r["discovery"]], "o-", color=c, label=f"{NAME[ion]} ({r['n_lines']:,} lines)")
        fr = sorted(float(k) for k in r["concentration"]); axes[1].plot([r["concentration"][f"{f:g}"] for f in fr], fr, "o-", color=c, label=NAME[ion])
        axes[2].bar([f"{NAME[ion]}\nmatrix"], [r["K128"]["bytes_matrix"] / 1024], color=c, alpha=0.5)
        axes[2].bar([f"{NAME[ion]}\ntables"], [r["K128"]["bytes_tables"] / 1024], color=c)
    axes[0].set_xscale("log"); axes[0].set_yscale("log"); axes[0].set_xlabel("events drawn"); axes[0].set_ylabel("distinct exit lines"); axes[0].set_title("discovery of exit lines", fontsize=9); axes[0].legend(fontsize=7)
    axes[1].set_xscale("log"); axes[1].set_xlabel("exit lines, largest energy first"); axes[1].set_ylabel("fraction of the exit energy"); axes[1].set_title("energy concentration", fontsize=9); axes[1].legend(fontsize=7)
    axes[2].set_yscale("log"); axes[2].set_ylabel("kB"); axes[2].set_title("the 128-group kernel as saved", fontsize=9)
    fig.tight_layout()
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}), ("png", {"Software": None})):
        p = out_dir / f"exit_table_audit.{ext}"; fig.savefig(p, metadata=meta, dpi=150, bbox_inches="tight"); written.append(p)
    plt.close(fig)
    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", choices=sorted(G1_N))
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--figure", action="store_true")
    a = ap.parse_args()
    path = HERE / "exit_tables.json"
    records = json.loads(path.read_text()) if path.exists() else {}
    if a.ion:
        records[a.ion] = audit(a.ion, a.n)
        path.write_text(json.dumps(records, indent=1, default=float) + "\n")
        print("wrote", path)
    if a.figure:
        for p in figure(records):
            print("wrote", p)


if __name__ == "__main__" and "--markdown" not in sys.argv:
    main()


def markdown(records=None):
    """The report's table, machine-generated from exit_tables.json."""
    records = records or json.loads((HERE / "exit_tables.json").read_text())
    NAME = {"57LaII": "La II", "58CeII": "Ce II", "60NdII": "Nd II"}
    L = ["| ion | packets × seeds | events | lines in the list | opacity lines | distinct exit lines | exact rest frequencies | exit energy in opacity lines | lines carrying 50 / 90 / 99 / 99.9 % | kept per group at f = 0.1 / 0.2 / 0.5 / 0.9 / 0.99 / 0.999 | K128 kB: total = matrix + tables | G1 record | K at N_g*: kB (tables) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for ion, a in records.items():
        c = a["concentration"]; k = a["kept_per_group"]; K = a["K128"]; s = a["K_ng_star"]
        L.append(f"| {NAME[ion]} | {a['n']:,} × {len(a['build_seeds'])} | {a['n_events']:,} | {a['n_lines']:,} | {a['n_opacity']:,} | {a['n_distinct_exit_lines']:,} | {a['exact_line_fraction']:.3f} | "
                 f"{a['exit_energy_in_opacity_lines']:.3f} | {c['0.5']:,} / {c['0.9']:,} / {c['0.99']:,} / {c['0.999']:,} | "
                 " / ".join(f"{k[f]['kept']:,}" for f in ("0.1", "0.2", "0.5", "0.9", "0.99", "0.999") if f in k) + " | "
                 f"{K['bytes_total']/1024:.0f} = {K['bytes_matrix']/1024:.0f} + {K['bytes_tables']/1024:.0f} | {K['g1_record_bytes']/1024:.0f} kB, {K['g1_record_n_exit']:,} | "
                 + (f"N_g* = {s['ng']}: {s['bytes_total']/1024:.0f} ({s['bytes_tables']/1024:.0f}) |" if s else "— |"))
    L += ["", "| ion | distinct exit lines after 10⁴ / 10⁵ / 10⁶ / all events |", "|---|---|"]
    for ion, a in records.items():
        d = a["discovery"]
        def at(n):
            return next((x["distinct"] for x in d if x["events"] >= n), d[-1]["distinct"])
        L.append(f"| {NAME[ion]} | {at(1e4):,} / {at(1e5):,} / {at(1e6):,} / {d[-1]['distinct']:,} (of {d[-1]['events']:,} events) |")
    return "\n".join(L)


if __name__ == "__main__" and "--markdown" in sys.argv:
    print(markdown())
