"""Paper IV Phase 8: multi-shell R2 vs B2 on a benchmark state -- does spatial
integration across shells on opposite sides of the saturation crossover
cancel the single-zone closure error?

One ZonedAtom over the transported shells (`--shells a-b` of the committed
state, or `--fine N` sub-shells between the photosphere and the edge built
from the model's profile), the legs R2 / B2 / Bbin2 (A2 needs per-shell
kernels and is deferred), the energy ledger, per-shell tallies, the escaped
energy per band split by the last shell a packet interacted in (the
band-forming weights f_s(b)), and the comparison block:

    dm_s(b)      the single-zone closure error on each transported shell
                 (from the committed neighbour / zone files where they
                 exist, else run here)
    dm_multi(b)  the multi-shell closure error
    dm_mix(b)    sum_s f_s(b) dm_s(b), the linear mixture
    cancel(b)    1 - dm_multi(b) / dm_mix(b)  (~0: shells only average;
                 ~1: the error cancels across the crossover)

Gate 4, pre-declared (paper4/README.md): per live band, |dm_multi| >= 0.5
|dm_mix| with signs kept everywhere -> survives; 0.2-0.5 somewhere ->
partial cancellation; < 0.2 or a sign flip in a band the single zone had at
> 1 mag -> cancels; gray if the identity, the crossing bookkeeping or the
identical-shell tests are out of bounds.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase11_observables",
          ROOT / "paper3/phase12_grid", ROOT / "paper4/phase2_energy", ROOT / "paper4/phase1_benchmarks"):
    sys.path.insert(0, str(p))

from sobolev.constants import C                                   # noqa: E402
from sobolev import photometry as phot                            # noqa: E402
from sobolev.ejecta import EjectaState                            # noqa: E402
from sobolev.zoned_atom import ZonedAtom                          # noqa: E402
from sobolev.forest_stats import band_saturation                  # noqa: E402
from sobolev.energy_packets import energy_accounting              # noqa: E402
from forest_mc import run_mc                                      # noqa: E402
from observables import observe, LAM_WIN, N_SPEC, BAND3800        # noqa: E402
from grid import photometer, rss_mb                               # noqa: E402
import legs as L                                                  # noqa: E402

DAY = 86400.0
SHELL_LEGS = ("R2", "B2", "Bbin2")


def band_weights_by_shell(res_list, edges, n_shell):
    """Escaped energy per band split by the shell of the packet's last
    interaction (-1: never interacted), summed over seeds, normalised per band."""
    from sobolev.photometry import BANDS_PHOT
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    out = {}
    for b, (lo_a, hi_a) in BANDS_PHOT.items():
        lo, hi = C / (hi_a * 1e-8), C / (lo_a * 1e-8)
        tot = np.zeros(n_shell + 1)
        for r in res_list:
            esc = (r["fate"] == 1) & (r["nu_out_all"] >= lo) & (r["nu_out_all"] < hi)
            e = (r["w"] * r["nu_out_all"])[esc]
            s = r["last_shell"][esc] + 1                       # -1 -> 0 (no interaction)
            np.add.at(tot, s, e)
        out[b] = (tot / tot.sum()).tolist() if tot.sum() > 0 else [np.nan] * (n_shell + 1)
    return out


def run_shells(state, shells, n, legs=SHELL_LEGS, seeds=L.SEEDS, stages=("II",), tau_min=1e-3,
               emis_cut=1e-6, relativity="worldline", budget_s=None, verbose=True):
    t0 = time.time()
    zone = state.transport_zone(shells)
    atom = ZonedAtom.from_state(state, shells, stages=stages, tau_min=tau_min, emis_cut=emis_cut)
    t_build = time.time() - t0
    nb = (C / (BAND3800[1] * 1e-8), C / (BAND3800[0] * 1e-8))
    per_shell = []
    for k, s in enumerate(shells):
        v = atom.shell_view(k)
        per_shell.append(dict(shell=int(s), v_c=float(state.v_edges[s] / C), rho=float(state.rho[s]),
                              T=float(state.T_gas[s]), n_opacity=int(v.n_opacity),
                              **(band_saturation(v, *nb) if v.n_opacity >= 2 else {})))
    lo, hi = (float(x) for x in phot.nu_edges(*LAM_WIN, 1))
    l_core = phot.planck_luminosity(lo, hi, zone["r_core"], zone["t_core"])
    edges = phot.nu_edges(*LAM_WIN, N_SPEC)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    row = dict(state=state.meta.get("name"), t_d=state.t / DAY, shells=list(map(int, shells)), n_shell=len(shells),
               zone={k: v for k, v in zone.items()}, n=n, seeds=list(seeds), stages=list(stages), tau_min=tau_min,
               emis_cut=emis_cut, relativity=relativity, n_lines=int(atom.n_lines_total),
               n_opacity_union=int(atom.n_opacity), n_opacity_shell=atom.n_opacity_shell.tolist(),
               macro_entries=int(atom.macro.n_entries), macro_dep_entries=int(atom.macro.n_dep_entries),
               t_build=t_build, rss_mb_atom=rss_mb(), per_shell=per_shell, L_core_window=l_core,
               lam_window=list(LAM_WIN), n_spec=N_SPEC, legs={}, timing={}, git=L.git_sha())
    if verbose:
        print(f"{row['state']} t={row['t_d']:g} d shells {shells[0]}-{shells[-1]}: union {atom.n_opacity} lines "
              f"({', '.join(str(x) for x in atom.n_opacity_shell)}), macro {atom.macro.n_entries} entries "
              f"({atom.macro.n_dep_entries} shell-dependent), build {t_build:.0f}s, rss {rss_mb():.0f} MB", flush=True)
    for tag in legs:
        spec = L.LEGS[tag]
        tl = time.time()
        res = [run_mc(atom, zone["r_core"], zone["r_out"], zone["t_exp"], lo, hi, n, spec["mode"], seed=s,
                      t_core=zone["t_core"], relativity=relativity, max_steps=L.MAX_STEPS, packets="energy",
                      launch_weight="energy", wall_s=budget_s, reprocess=spec.get("reprocess")) for s in seeds]
        o = photometer(observe(res, l_core, spec["scale"]), edges, nu_c, phot.D_40MPC)
        per_seed = [photometer(observe([r], l_core, spec["scale"]), edges, nu_c, phot.D_40MPC)["mags"] for r in res]
        o["mags_seed_std"] = {b: float(np.std([m[b] for m in per_seed], ddof=1)) if len(res) > 1 else np.nan for b in o["mags"]}
        acc = [energy_accounting(r) for r in res]
        o["energy"] = {k: float(np.mean([a[k] for a in acc])) for k in acc[0] if k != "packets"}
        o["ledger"] = L.ledger(res)
        o["events_per_packet"] = float(np.mean([r["n_events"].mean() for r in res]))
        o["n_events_shell"] = np.sum([r["n_events_shell"] for r in res], axis=0).tolist()
        o["n_kpackets_shell"] = np.sum([r["n_kpackets_shell"] for r in res], axis=0).tolist()
        o["n_dead_end_shell"] = np.sum([r["n_dead_end_shell"] for r in res], axis=0).tolist()
        o["e_thermal_shell"] = (np.sum([r["e_thermal_shell"] for r in res], axis=0) / np.sum([r["accounting"]["E_inj"] for r in res])).tolist()
        o["n_cross_in"] = np.sum([r["n_cross_in"] for r in res], axis=0).tolist()
        o["n_cross_out"] = np.sum([r["n_cross_out"] for r in res], axis=0).tolist()
        o["band_weights_by_shell"] = band_weights_by_shell(res, edges, len(shells))
        o["steps"] = int(np.mean([r["steps"] for r in res]))
        o["mode"] = spec["mode"]; o["scale"] = spec["scale"]
        o["t_wall"] = row["timing"][tag] = time.time() - tl
        row["legs"][tag] = o
        if verbose:
            e = o["energy"]
            print(f"  {tag:6s} {spec['mode']:17s} {o['t_wall']:7.1f}s ev/pkt={o['events_per_packet']:6.1f} steps={o['steps']} "
                  f"esc={e['esc_frac']:.3f} core={e['core_frac']:.3f} work={e['work_frac']:+.4f} "
                  f"events/shell={o['n_events_shell']} g={o['mags'].get('g', np.nan):.2f} z={o['mags'].get('z', np.nan):.2f} "
                  f"K={o['mags'].get('K', np.nan):.2f}", flush=True)
    if "R2" in row["legs"]:
        ref = row["legs"]["R2"]
        for tag, o in row["legs"].items():
            o["dm_vs_R2"] = phot.delta_mag(o["mags"], ref["mags"])
            o["dcolor_vs_R2"] = {k: o["colors"][k] - ref["colors"][k] for k in ref["colors"]}
            o["dm_bol_vs_R2"] = phot.bol_delta_mag(o["L_bol"], ref["L_bol"])
    row["t_wall"] = time.time() - t0
    row["rss_mb"] = rss_mb()
    return row


def compare(row, single_files, weights_leg="R2"):
    """The cancellation block from the single-zone leg files of the transported
    shells: {shell: path}. dm_s(b) is each shell's single-zone B2 - R2; the
    mixture dm_mix(b) = sum_s f_s(b) dm_s(b) uses the band-forming weights
    f_s(b) (escaped energy per band by the shell of last interaction) of
    `weights_leg` in the multi-shell run (R2: where the reference's light in
    the band forms; the B2-weighted mixture is reported beside it), and
    cancel(b) = 1 - dm_multi / dm_mix. Also the photospheric shell's own
    single-zone error dm_ph(b), the number the single-zone verdicts quote."""
    shells = row["shells"]
    dm_s = {}
    for s, f in single_files.items():
        d = json.loads(Path(f).read_text())
        dm_s[int(s)] = d["legs"]["B2"]["dm_vs_R2"]
    out = {}
    if "B2" not in row["legs"] or not dm_s or weights_leg not in row["legs"]:
        return out
    fw = row["legs"][weights_leg]["band_weights_by_shell"]
    fw_b = row["legs"]["B2"]["band_weights_by_shell"]
    ph = shells[0]
    for b in row["legs"]["B2"]["dm_vs_R2"]:
        have = [k for k, s in enumerate(shells) if s in dm_s]
        w = np.array(fw[b][1:]); wb = np.array(fw_b[b][1:])     # drop the never-interacted column
        if not have or not np.isfinite(w).all():
            continue
        dm = np.array([dm_s[shells[k]][b] for k in have])
        def _mix(ww):
            ww = ww[have]
            return float(np.sum(ww * dm) / ww.sum()) if ww.sum() > 0 else np.nan
        mix, mix_b = _mix(w), _mix(wb)
        multi = row["legs"]["B2"]["dm_vs_R2"][b]
        out[b] = dict(dm_multi=multi, dm_mix=mix, dm_mix_B2w=mix_b, weights_leg=weights_leg,
                      dm_ph=float(dm_s[ph][b]) if ph in dm_s else np.nan,
                      dm_shells={str(shells[k]): float(dm_s[shells[k]][b]) for k in have},
                      weights={str(shells[k]): float(w[k] / w[have].sum()) for k in have},
                      weights_B2={str(shells[k]): float(wb[k] / wb[have].sum()) for k in have},
                      cancel=(1.0 - multi / mix) if (np.isfinite(mix) and mix != 0) else np.nan,
                      ratio=(multi / mix) if (np.isfinite(mix) and mix != 0) else np.nan,
                      ratio_B2w=(multi / mix_b) if (np.isfinite(mix_b) and mix_b != 0) else np.nan)
    return out


def gate4(row, seed_sigma=0.1):
    """The pre-declared Gate 4 reading of a comparison block, per live band of
    the multi-shell R2 (verdict.live_bands): ratio = dm_multi / dm_mix with
    dm_mix the R2-band-forming-weighted mixture of the shells' single-zone
    errors.  Survives: |ratio| >= 0.5 with the sign kept in every live band.
    Partial: some live band in 0.2-0.5.  Cancels: every live band < 0.2, or
    a sign flip in a band whose single-zone error exceeded 1 mag.  Gray:
    the energy identity, the crossing bookkeeping or a missing single-zone
    file; a band whose |dm_mix| is within `seed_sigma` of zero is not read
    (the ratio is noise there)."""
    sys.path.insert(0, str(HERE.parent / "phase3_legs"))
    from verdict import live_bands
    cmp_ = row.get("comparison") or {}
    if not cmp_:
        return dict(outcome="GRAY", reason="no comparison block")
    for k in ("R2", "B2"):
        o = row["legs"][k]
        if abs(o["ledger"].get("identity", 0.0)) > 1e-9 or o["energy"].get("dep_cm", 0.0) != 0.0:
            return dict(outcome="GRAY", reason=f"{k}: energy identity {o['ledger'].get('identity')} / dep_cm")
    live = [b for b in live_bands(row) if b in cmp_]
    read = {}
    for b in live:
        c = cmp_[b]
        if not np.isfinite(c["dm_mix"]) or abs(c["dm_mix"]) < seed_sigma:
            continue
        flip = np.sign(c["dm_multi"]) != np.sign(c["dm_mix"])
        big = max(abs(v) for v in c["dm_shells"].values()) > 1.0
        read[b] = dict(ratio=c["ratio"], flip=bool(flip), big_single=bool(big), dm_multi=c["dm_multi"], dm_mix=c["dm_mix"])
    if not read:
        return dict(outcome="GRAY", reason="no live band with |dm_mix| above the seed noise", live=live)
    if any(r["flip"] and r["big_single"] for r in read.values()) or all(abs(r["ratio"]) < 0.2 for r in read.values()):
        out = "CANCELS"
    elif all(abs(r["ratio"]) >= 0.5 and not r["flip"] for r in read.values()):
        out = "SURVIVES"
    else:
        out = "PARTIAL"
    return dict(outcome=out, live=live, bands=read)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state")
    ap.add_argument("--shells", default=None, help="a-b (inclusive) of the state's shells")
    ap.add_argument("--fine", type=int, default=None, help="N sub-shells between the photosphere and the edge")
    ap.add_argument("--legs", default=",".join(SHELL_LEGS))
    ap.add_argument("--n", type=int, default=300000)
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--stages", default="II")
    ap.add_argument("--tau-min", type=float, default=1e-3)
    ap.add_argument("--emis-cut", type=float, default=1e-6)
    ap.add_argument("--relativity", default="worldline")
    ap.add_argument("--budget", type=float, default=None)
    ap.add_argument("--single", default="", help="shell=file,... single-zone leg files for the comparison")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    st = EjectaState.from_json(a.state)
    if a.fine:
        import build
        st = build.refine_above_photosphere(st, a.fine)
        s0 = st.meta["photospheric_shell"]
        shells = list(range(s0, st.n_shell))
    elif a.shells:
        lo_, hi_ = (int(x) for x in a.shells.split("-"))
        shells = list(range(lo_, hi_ + 1))
    else:
        shells = list(range(st.meta.get("photospheric_shell", 0), st.n_shell))
    seeds = tuple(int(s) for s in a.seeds.split(","))
    row = run_shells(st, shells, a.n, tuple(a.legs.split(",")), seeds, tuple(a.stages.split(",")), a.tau_min,
                     None if a.emis_cut <= 0 else a.emis_cut, a.relativity or None, a.budget)
    row["fine"] = a.fine
    if a.single:
        files = dict(kv.split("=") for kv in a.single.split(","))
        row["comparison"] = compare(row, files)
        row["gate4"] = gate4(row)
        for b, c in row["comparison"].items():
            print(f"  {b}: dm_multi={c['dm_multi']:+.2f} dm_mix={c['dm_mix']:+.2f} ratio={c['ratio']:.2f} cancel={c['cancel']:+.2f} "
                  f"shells={ {k: round(v, 2) for k, v in c['dm_shells'].items()} } weights={ {k: round(v, 2) for k, v in c['weights'].items()} }")
    out = a.out or HERE / f"shells_{row['state']}_t{row['t_d']:g}_s{shells[0]}-{shells[-1]}{'_f%d' % a.fine if a.fine else ''}.json"
    Path(out).write_text(json.dumps(row, indent=1, default=float) + "\n")
    print(f"wrote {out} in {row['t_wall']:.0f}s")


if __name__ == "__main__":
    main()
