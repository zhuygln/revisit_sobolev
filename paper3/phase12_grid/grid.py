"""Paper III, phase 12: one (M_ej, v_ej, X_lan) model of the physically powered
grid -- the F41 scorecard on a heating-powered source.

The source is `sobolev.source.SourceModel`: (M_ej, v_ej) set L(t), T_eff(t) and
R_ph(t) through radioactive heating, Barnes+2016 thermalization and one-zone
diffusion at a FIXED grey kappa; X_lan sets only the line opacity of the
four-ion La/Ce/Pr/Nd blend (equal split). The transport is called exactly as
§4.37 called it -- launch from r_core = R_ph at t_core = T_eff, worldline
relativity, the five legs of `observables.LEGS` -- so nothing in the closure
comparison changed, only what illuminates the shell and how the shell's density
is set. `observables.epoch` is left untouched (the F40/F41 path stays
bit-identical); `run_epoch` here is its sibling with a state dict injected, a
wall-time guard, real-filter photometry alongside the top-hats, and per-run
timing.

Two things the F41 harness did not need, because F41's forests had S ~ 50
(§4.39 measures why both are needed at S ~ 10^4-10^5):

* the radiation energy is conserved -- `observe(core="conserving")`
  renormalizes every leg's escaped spectrum to the core's window luminosity
  (`photometry._scale`): the core re-emits what returns to it and the gas
  re-radiates what the packets deposit, with the escaped spectrum's shape.
  With the F41 absorbing core, 75% of the central model's packets return at
  3 d and the reference deposits 45% of the energy in level excitation, so a
  grouped closure that interacts less reads as a 1.9 mag *bolometric* error
  that is the harness's bookkeeping, not the closure's; `dm_bol_absorbing`,
  `f_return` and `f_dep` keep those numbers as diagnostics;
* re-absorption chains that have not left their emitting line after
  `CHAIN_MAX` draws are thermalized (`chain_overflow="absorb"`) instead of
  aborting the run; the trapped count is recorded per leg. `max_steps` is
  raised to 10^6 (the step loop is vectorized over packets, so the cap is on
  the slowest packet's history, not on work).

Wall guard: a probe of the reference leg (`N_PROBE` packets) times the
transport per packet -- at S ~ 1e5 a packet interacts several hundred times
and costs ~10 ms, a thousand times the S ~ 50 cost. The per-packet cost is
sublinear in n (the vectorized step loop has a per-step overhead: 500
packets cost 50-100 ms each, 5000 cost 2-9 ms, 125 000 cost 1.8 ms), so a
small probe overestimates and the projection is conservative;
if the epoch at `n` would exceed `budget_s` (`COST_RUNS` reference-run
equivalents -- the closure legs are cheap), every run of the epoch uses the
reduced `n_used` (floor `N_FLOOR`) and the row is marked `reduced_n`; if even
the floor would exceed `OVER_FACTOR x budget_s` the epoch is marked
`over_budget` with the projected cost and no transport is run; a single run
(probe included) that outlives `OVER_FACTOR x budget_s` is abandoned and the
epoch marked `wall`. A RuntimeError from `run_mc` marks the epoch `max_steps`
or `chain`. The JSON is rewritten
after every epoch so a killed run keeps its rows.

Reruns (§4.41): `--chain-max` changes the chain cap, `--n-override` skips the
probe and runs every leg at that n (so a cell can be reproduced at its stored
`n_used` with the stored seeds), `--t-scale` perturbs the source's launch
temperature at fixed L (`SourceModel(t_scale=...)`; the JSON gets a `_T<s>`
suffix) and `--t-scale-gas` scales T_gas with it. A cell redone at a larger
budget is `--epochs t --budget 5400 --out <fresh file>` (see `run_grid.py
--redo`); the row records the `budget_s` and `chain_max` it ran with.

Usage: python grid.py --mass 0.01 --v 0.1 --xlan 0.01 [--n 300000] [--budget 1500]
       python grid.py --mass 0.01 --v 0.1 --xlan 0.01 --dry-run   (atoms only)
"""
import argparse, json, resource, subprocess, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper3" / "phase11_observables"))
import observables as obs                      # noqa: E402  (adds the paper2/paper3 paths)
from observables import (BAND3800, BLEND, LAM_WIN, LEGS, N_SPEC, SEEDS,   # noqa: E402
                         build_atom, observe, band_saturation, run_mc,
                         RedistributionKernel, phot, C)
from sobolev.source import SourceModel, DAY   # noqa: E402

EPOCHS = (0.5, 1.0, 2.0, 3.0, 5.0, 7.0)
N_DEFAULT = 300_000
N_FLOOR = 20_000
BUDGET_S = 1500.0
OVER_FACTOR = 3.0        # an epoch may cost up to OVER_FACTOR x budget_s at N_FLOOR
N_PROBE = 5_000
MAX_STEPS = 1_000_000
CHAIN_MAX = 2_000
CORE = "conserving"
COST_RUNS = 5.0          # epoch cost in reference-run units: 3 seeds of the reference
                         # plus the four closure legs, which cost 2-20% of it each
RELATIVITY = "worldline"
NG = 32
PASSBANDS = phot.load_passbands()


SOURCE_KEYS = ("L", "Qdot", "f_th", "tau_d", "T_eff", "T_eff_grey", "R_ph", "fth_clamped",
               "v_ph", "v_ph_floored", "tau_grey", "t_scale", "t_scale_gas")


def model_name(m, v, x, t_scale=1.0):
    s = f"model_M{m:g}_v{v:g}_X{x:g}"
    return s if t_scale == 1.0 else s + f"_T{t_scale:g}"


def build_state_atom(st, x_lan):
    """The four-ion blend at a source state: (ForestAtom, n_ion dict)."""
    return build_atom("blend", st, x_lan)


def git_sha():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                              text=True, cwd=ROOT).stdout.strip()
    except Exception:
        return None


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def photometer(o, edges, nu_c, dist):
    """Add real-filter magnitudes to an `observe()` dict; keep the top-hats."""
    lnu = np.asarray(o["L_nu"], float)
    o["mags_tophat"] = o["mags"]
    o["colors_tophat"] = o["colors"]
    o["mags"] = phot.magnitudes(nu_c, lnu, PASSBANDS, dist, edges)
    o["colors"] = phot.colors(o["mags"])
    return o


def run_epoch(st, x_lan, n, budget_s=BUDGET_S, ng=NG, relativity=RELATIVITY, dry_run=False,
              chain_max=CHAIN_MAX, n_override=None, atom=None):
    """One epoch on a source state `st` (from `SourceModel.state`).

    `n_override` skips the probe and runs every leg at that n; `atom` injects a
    prebuilt (ForestAtom, n_ion) instead of the blend (tests, reruns).
    """
    t_d = st["t_exp"] / DAY
    t0 = time.time()
    atom, n_ion = atom if atom is not None else build_state_atom(st, x_lan)
    row = {"t_d": t_d, "n_ion": n_ion, "T_gas": st["T_gas"], "t_core": st["t_core"],
           "rho": st["rho"], "r_core": st["r_core"], "r_out": st["r_out"],
           "v_core": st["r_core"] / (C * st["t_exp"]), "v_out": st["r_out"] / (C * st["t_exp"]),
           "x_lan": x_lan, "chain_max": chain_max, "budget_s": budget_s,
           "source": {k: st[k] for k in SOURCE_KEYS if k in st},
           "n_opacity": int(atom.n_opacity),
           "tau_max": float(atom.op_tau.max()) if atom.n_opacity else 0.0,
           "rss_mb_atom": rss_mb()}
    nb = (C / (BAND3800[1] * 1e-8), C / (BAND3800[0] * 1e-8))
    if atom.n_opacity >= 2:
        row.update({f"band_{k}": v for k, v in band_saturation(atom, *nb).items()})
    if atom.n_opacity < 10:
        row.update(status="skipped", skipped=f"{atom.n_opacity} opacity lines", t_wall=time.time() - t0)
        return row
    if dry_run:
        row.update(status="dry_run", t_wall=time.time() - t0)
        return row

    lo, hi = (float(x) for x in phot.nu_edges(*LAM_WIN, 1))
    l_core = phot.planck_luminosity(lo, hi, st["r_core"], st["t_core"])
    row["L_core_window"] = l_core
    edges = phot.nu_edges(*LAM_WIN, N_SPEC)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    dist = phot.D_40MPC

    def mc(mode, seed, n_run, wall_s=None, **kw):
        return run_mc(atom, st["r_core"], st["r_out"], st["t_exp"], lo, hi, n_run, mode,
                      seed=seed, t_core=st["t_core"], relativity=relativity,
                      max_steps=MAX_STEPS, chain_max=chain_max, chain_overflow="absorb",
                      wall_s=wall_s, **kw)

    def failed(e):
        m = str(e)
        return "chain" if "chain" in m else "wall" if "wall" in m else "max_steps"

    # --- wall guard: probe the reference leg -------------------------------
    if n_override is not None:
        n_used = int(n_override)
        row.update(n_used=n_used, n_override=True, status="ok" if n_used >= n else "reduced_n")
    else:
        n_probe = min(n, N_PROBE)
        tp = time.time()
        try:
            pr = mc("sobolev_branch", SEEDS[0], n_probe, wall_s=budget_s)
        except RuntimeError as e:
            row.update(status=failed(e), error=str(e)[:200], t_wall=time.time() - t0)
            return row
        per_packet = (time.time() - tp) / n_probe
        n_fit = int(budget_s / (COST_RUNS * per_packet))
        n_used = n if n_fit >= n else max(N_FLOOR, n_fit)
        row.update(n_used=n_used, probe_s_per_packet=per_packet, probe_n=n_probe,
                   probe_f_return=phot.return_fraction(pr), probe_f_dep=phot.deposited_fraction(pr),
                   probe_n_trapped=int(pr["n_trapped"]),
                   projected_s=COST_RUNS * per_packet * n_used,
                   status="ok" if n_used == n else "reduced_n")
        if row["projected_s"] > OVER_FACTOR * budget_s:
            row.update(status="over_budget", t_wall=time.time() - t0)
            return row

    timing = {}
    try:
        ref_res, ei, eo = [], [], []
        tl = time.time()
        for s in SEEDS:
            r = mc("sobolev_branch", s, n_used, collect_events=True, wall_s=OVER_FACTOR * budget_s)
            ref_res.append(r); e = r["events"]; ei.append(e[0]); eo.append(e[1])
        timing["reference"] = time.time() - tl
        ref = photometer(observe(ref_res, l_core, CORE), edges, nu_c, dist)

        k_lo, k_hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
        nu_in, nu_out = np.concatenate(ei), np.concatenate(eo)
        kern = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size),
                                                      ng, nu_lo=k_lo, nu_hi=k_hi)
        legs = {}
        for tag, mode in LEGS:
            tl = time.time()
            kw = {"kernel": kern} if mode.endswith("_group") else {}
            o = photometer(observe([mc(mode, s, n_used, wall_s=OVER_FACTOR * budget_s, **kw)
                                    for s in SEEDS], l_core, CORE),
                           edges, nu_c, dist)
            o["dm"] = phot.delta_mag(o["mags"], ref["mags"])
            o["dcolor"] = {k: o["colors"][k] - ref["colors"][k] for k in ref["colors"]}
            o["dm_tophat"] = phot.delta_mag(o["mags_tophat"], ref["mags_tophat"])
            o["dcolor_tophat"] = {k: o["colors_tophat"][k] - ref["colors_tophat"][k]
                                  for k in ref["colors_tophat"]}
            o["dm_bol"] = phot.bol_delta_mag(o["L_bol"], ref["L_bol"])
            o["dm_bol_absorbing"] = phot.bol_delta_mag(o["L_bol_absorbing"], ref["L_bol_absorbing"])
            o["dF_b3800"] = o["b3800"] / ref["b3800"] - 1 if ref["b3800"] > 1e-3 else None
            o["t_wall"] = timing[tag] = time.time() - tl
            legs[tag] = o
    except RuntimeError as e:
        row.update(status=failed(e), error=str(e)[:200], t_wall=time.time() - t0)
        return row
    row.update(ref=ref, legs=legs, timing=timing, t_wall=time.time() - t0, rss_mb=rss_mb())
    return row


HEADER_DEFAULTS = {"t_scale": 1.0, "t_scale_gas": False}   # absent from pre-§4.41 JSONs


def run_model(m_ej, v_ej, x_lan, n=N_DEFAULT, epochs=EPOCHS, kappa=1.0, budget_s=BUDGET_S,
              out=None, dry_run=False, verbose=True, t_scale=1.0, t_scale_gas=False,
              chain_max=CHAIN_MAX, n_override=None):
    src = SourceModel(m_ej, v_ej, kappa=kappa, t_scale=t_scale, t_scale_gas=t_scale_gas)
    out = Path(out) if out else HERE / "grid" / f"{model_name(m_ej, v_ej, x_lan, t_scale)}.json"
    out.parent.mkdir(exist_ok=True)
    d = {"m_ej_msun": m_ej, "v_ej_c": v_ej, "x_lan": x_lan, "kappa_src": kappa,
         "v_ph_frac": src.v_ph_frac, "tau_d_s": src.tau_d, "t_peak_s": src.t_peak,
         "fth_clamped": src.fth_clamped, "composition": "blend", "ions": list(BLEND),
         "n": n, "seeds": list(SEEDS), "epochs": list(epochs), "relativity": RELATIVITY,
         "ng": NG, "lam_window": list(LAM_WIN), "n_spec": N_SPEC, "distance_cm": phot.D_40MPC,
         "filter_set": "real", "filters": {b: p.name for b, p in PASSBANDS.items()},
         "budget_s": budget_s, "over_factor": OVER_FACTOR, "n_floor": N_FLOOR,
         "core": CORE, "max_steps": MAX_STEPS, "chain_max": chain_max,
         "t_scale": t_scale, "t_scale_gas": t_scale_gas, "n_override": n_override,
         "git": git_sha(), "dry_run": dry_run, "rows": []}
    if verbose:
        print(f"{out.name}: M={m_ej} v={v_ej} X={x_lan}  tau_d={src.tau_d/DAY:.2f} d  "
              f"{n} packets x {len(SEEDS)} seeds  {'DRY RUN' if dry_run else ''}", flush=True)
    # resume: keep completed transport rows of an earlier run with the same
    # settings; guard-decided rows (over_budget, wall, ...) are re-decided
    done = {}
    if out.exists() and not dry_run:
        try:
            prev = json.loads(out.read_text())
            same = all(prev.get(k, HEADER_DEFAULTS.get(k)) == d[k]
                       for k in ("n", "budget_s", "core", "chain_max", "max_steps",
                                 "over_factor", "n_floor", "t_scale", "t_scale_gas", "n_override"))
            if same:
                done = {r["t_d"]: r for r in prev["rows"] if r.get("status") in ("ok", "reduced_n")}
        except Exception:
            done = {}
    t0 = time.time()
    for t_d in epochs:
        st = src.state(t_d)
        if t_d in done:
            r = dict(done[t_d], resumed=True)
        else:
            r = run_epoch(st, x_lan, n, budget_s=budget_s, dry_run=dry_run,
                          chain_max=chain_max, n_override=n_override)
        d["rows"].append(r)
        d["t_wall"] = time.time() - t0
        out.write_text(json.dumps(d, indent=1))
        if verbose:
            s = r["status"]
            line = (f"  t={t_d:4.1f} d  T_eff={st['T_eff']:6.0f} K  L={st['L']:.2e}  "
                    f"n_ion(Ce)={r['n_ion'].get('58CeII', 0):.2e}  n_op={r['n_opacity']:7d}  "
                    f"tau_max={r['tau_max']:.2e}  S={r.get('band_S_band', float('nan')):8.1f}  "
                    f"[{s}, {r['t_wall']:.0f}s, {r['rss_mb_atom']:.0f} MB]")
            if s in ("ok", "reduced_n"):
                cb = r["legs"]["C_both"]
                line += (f"  n_used={r['n_used']}  f_ret={r['ref']['f_return']:.2f}  "
                         f"f_dep={r['ref']['f_dep']:.2f}  "
                         f"dm_bol={cb['dm_bol']:+.3f} ({cb['dm_bol_absorbing']:+.2f} abs)  "
                         f"d(g-r)={cb['dcolor']['g-r']:+.3f}  d(i-J)={cb['dcolor']['i-J']:+.3f}")
            elif s == "over_budget":
                line += f"  projected {r['projected_s']/3600:.1f} h at n={N_FLOOR}"
            print(line, flush=True)
    d["complete"] = True
    out.write_text(json.dumps(d, indent=1))
    if verbose:
        print(f"wrote {out}  [{time.time()-t0:.0f}s]", flush=True)
    return d


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass", type=float, required=True)
    ap.add_argument("--v", type=float, required=True)
    ap.add_argument("--xlan", type=float, required=True)
    ap.add_argument("--n", type=float, default=N_DEFAULT)
    ap.add_argument("--epochs", default=",".join(str(e) for e in EPOCHS))
    ap.add_argument("--kappa", type=float, default=1.0)
    ap.add_argument("--budget", type=float, default=BUDGET_S)
    ap.add_argument("--out", default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--t-scale", type=float, default=1.0)
    ap.add_argument("--t-scale-gas", action="store_true")
    ap.add_argument("--chain-max", type=int, default=CHAIN_MAX)
    ap.add_argument("--n-override", type=int, default=None)
    a = ap.parse_args()
    run_model(a.mass, a.v, a.xlan, int(a.n), tuple(float(e) for e in a.epochs.split(",")),
              a.kappa, a.budget, a.out, a.dry_run, t_scale=a.t_scale, t_scale_gas=a.t_scale_gas,
              chain_max=a.chain_max, n_override=a.n_override)
