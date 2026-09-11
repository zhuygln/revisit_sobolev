"""Paper IV, Phase 10b: the Fontes et al. (2020) simplified problem as a
TIME-DEPENDENT light curve -- resolved Sobolev vs expansion vs line-binned
opacity, under complete thermal redistribution (eps = 1) and under
energy-conserving downward-macroatom fluorescence. Same ejecta, same atomic
data, same time grid, same heating; only the closure and the redistribution
change.

The transport is `run_mc` in time slabs (`t_stop` / `resume` / `launch_energy`,
see its docstring): a log grid of slabs from t0 to t1; per slab the state is
rebuilt at the slab's start (rho ∝ t^-3 from the Appendix C profile, T from
the chosen rule, LTE Saha), the zoned atom built from it, and every leg's
packet population transported from the previous checkpoint to the slab's
end with the slab's heating injected in the volume in proportion to mass
(uniform specific heating) at a uniform time within the slab. The trapped
radiation field at t0 (a T^4 V per shell, seven times the 4-16 d heating)
is the initial population. The innermost cell edge is a lossless mirror
(the centre of a sphere); cell 0's mass (2.5e-5 of the total) is dropped
and reported.

Temperature rules: `prescribed` = the Appendix C initial profile scaled
T ∝ t^-1 (adiabatic, radiation-dominated; no heating enters T; identical
for every leg); `radiation` = T_rad per shell from the packets present at
the slab boundary, (E_cm / a V)^(1/4), per leg (a one-slab lag; shells with
too few packets fall back to the prescribed T and are counted). Neither is
SuperNu's LTE energy equation; each is stated.

Every leg is guarded: `max_events` per packet per slab books a packet as
"capped" (its energy is reported as E_capped, the fraction f_capped is
printed next to every number), `wall_slab` per slab raises and the driver
records the failure and keeps the last checkpoint (a gray outcome for that
leg). `run.json` is rewritten atomically after every leg; `--resume`
continues each leg from its checkpoint with the stored per-slab seed, so a
resumed slab reproduces the killed one bit for bit.
"""
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3", ROOT / "paper3/phase11_observables",
          ROOT / "paper3/phase12_grid", ROOT / "paper4/phase2_energy", ROOT / "paper4/phase1_benchmarks",
          ROOT / "paper4/phase7_ionization", ROOT / "paper4/phase8_shells", str(HERE)):
    sys.path.insert(0, str(p))

from sobolev.constants import C, H                                # noqa: E402
from sobolev import photometry as phot                            # noqa: E402
from sobolev import timeslab as ts                                # noqa: E402
from sobolev.zoned_atom import ZonedAtom                          # noqa: E402
from forest_mc import run_mc                                      # noqa: E402
from observables import LAM_WIN, N_SPEC                           # noqa: E402
from grid import PASSBANDS, rss_mb                                # noqa: E402
import legs as L                                                  # noqa: E402
import fontes                                                     # noqa: E402
import build as p1build                                           # noqa: E402
from sobolev import xkn as xk                                     # noqa: E402

_XKN = {}


def xkn_model(cfg):
    """The published xkn secular component at the configured opacity (cached)."""
    key = (cfg.get("kappa"), cfg.get("n_modes", 1000))
    if key not in _XKN:
        _XKN[key] = xk.XknSecular(kappa=cfg["kappa"], n_modes=cfg.get("n_modes", 1000))
    return _XKN[key]

DAY = 86400.0
LEG_ID = {"R2": 1, "B2": 2, "Bbin2": 3, "Rth": 4, "Bth": 5, "Bbinth": 6}
REF_OF = {"R2": "R2", "B2": "R2", "Bbin2": "R2", "Rth": "Rth", "Bth": "Rth", "Bbinth": "Rth"}


def time_grid(t0_d, t1_d, n_slabs):
    return t0_d * DAY * (t1_d / t0_d) ** (np.arange(n_slabs + 1) / n_slabs)


def config_hash(cfg):
    return hashlib.sha256(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()[:16]


def write_json_atomic(path, obj):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=float))
    os.replace(tmp, path)


def seed_for(base, leg, k):
    return int(np.random.SeedSequence([int(base), LEG_ID[leg], int(k)]).generate_state(1)[0])


def build_state(cfg, t_s, T_override=None):
    """The model's state at the exact epoch t_s: the Fontes Appendix C ejecta,
    or the P1 composition on the published xkn secular structure (Phase 10c,
    `paper4/phase1_benchmarks/build.py::build_p1_xkn`)."""
    if cfg.get("model", "fontes") == "fontes":
        return fontes.build_fontes(t_d=t_s / DAY, n_shell=cfg["n_shell"], saha=True, t_s=t_s, T_override=T_override)
    if T_override is not None:
        raise NotImplementedError("the xkn model prescribes its temperatures (eq. 50); no override")
    return p1build.build_p1_xkn(t_s, kappa=cfg["kappa"], n_outer=cfg["n_outer"], x_lo=cfg["x_lo"], model=xkn_model(cfg))


def zone_of(cfg, st):
    """The transported shells of a state: fixed (`transport`, or 1..end) for
    Fontes; the shells outside the xkn photosphere for p1xkn (grows inward)."""
    if cfg.get("model", "fontes") == "p1xkn":
        return list(range(int(st.meta["photospheric_shell"]), st.n_shell))
    if cfg.get("transport"):
        lo_, hi_ = (int(x) for x in cfg["transport"].split("-")); return list(range(lo_, hi_ + 1))
    return list(range(1, st.n_shell))


def slab_sources(cfg, st, shells, t_a, t_b):
    """Energies injected in [t_a, t_b]: (E_boundary, E_heat_per_shell, weights,
    core_frac, t_core). Fontes: the heating law uniform in mass, nothing at
    the boundary. p1xkn: the xkn thick luminosity at the photosphere with
    T_ph, and each thin shell's own deposited heating (xkn eq. 47, 59)."""
    m_sh = st.shell_mass()[shells]
    if cfg.get("model", "fontes") == "fontes":
        E_heat, per = ts.heating_energy(m_sh, t_a, t_b)
        return 0.0, per, None, 0.0, float(st.T_gas[shells[0]])
    model = xkn_model(cfg)
    x_mid = 0.5 * (st.v_edges[1:] + st.v_edges[:-1])[shells] / (st.meta["v_max_c"] * C)
    E_thick = ts.thick_energy(model, t_a, t_b)
    E_heat, per = ts.thin_heating_energy(model, m_sh, np.minimum(x_mid, 1.0 - 1e-9), t_a, t_b)
    tot = E_thick + E_heat
    return E_thick, per, per, (E_thick / tot if tot > 0 else 0.0), float(st.meta["T_ph"])


def build_atom(cfg, st, shells):
    return ZonedAtom.from_state(st, shells, stages=tuple(cfg.get("stages", ("II", "III"))), tau_min=cfg["tau_min"],
                                emis_cut=cfg["emis_cut"], f_min=cfg["f_min"], dataset=cfg["dataset"])


def run_slab(cfg, st, atom, shells, leg, k, pop, sources, n_heat_k, t_a, t_b, out_dir, e_next):
    spec = L.LEGS[leg]
    E_bnd, E_heat_per, weights, core_frac, t_core = sources
    E_inject = float(E_bnd + E_heat_per.sum())
    lo, hi = (float(x) for x in phot.nu_edges(*cfg["lam_transport"], 1))
    seed_k = seed_for(cfg["seed"], leg, k)
    t0 = time.time()
    res = run_mc(atom, atom.r_edges[0], atom.r_edges[-1], t_a, lo, hi, int(n_heat_k), spec["mode"], seed=seed_k,
                 t_core=t_core, relativity="worldline", max_steps=cfg["max_steps"], packets="energy",
                 launch_weight="energy", launch="volume", launch_energy=E_inject, t_stop=t_b,
                 launch_core_frac=core_frac, launch_shell_weights=weights,
                 resume=pop.as_dict(), core=cfg["core"], core_max_passes=cfg["max_passes"],
                 max_events=cfg["max_events"], wall_s=cfg["wall_slab"], macro_kw=spec.get("macro_kw"))
    a = res["accounting"]
    obs = ts.escapes(res)
    np.savez(out_dir / f"esc_{leg}_{k:03d}.npz", **obs)
    pop_next = ts.carried_population(res, e_next)
    tally = dict(k=k, t_a=t_a, t_b=t_b, seed=seed_k, n_new=int(res["n_new"]), n_carried_in=int(pop.n),
                 E_carried_in=a["E_carried_in"], E_inj_new=a["E_inj_new"], E_inj_cm=a["E_inj_cm"], E_heat=E_inject,
                 E_boundary_in=float(E_bnd), E_heat_zone=float(E_heat_per.sum()), core_frac=float(core_frac), t_core=float(t_core),
                 n_zone=len(shells), s0=int(shells[0]), x_ph=float(st.meta.get("x_ph", np.nan)), T_ph=float(st.meta.get("T_ph", np.nan)),
                 L_thick=float(st.meta.get("L_thick", np.nan)), L_thin_xkn=float(st.meta.get("L_thin_xkn", np.nan)),
                 E_esc=a["E_esc"], E_core=a["E_core"], E_abs=a["E_abs"], E_capped=a["E_capped"], E_dep_lab=a["E_dep_lab"],
                 W=a["W"], E_carried_out=a["E_carried_out"], identity=a["identity_residual"],
                 n_paused=int(res["n_paused"]), n_capped=int(res["n_capped"]), n_escaped=int(res["n_escaped"]),
                 n_core=int(res["n_core"]), events_mean=float(res["n_events"].mean()) if res["n_events"].size else 0.0,
                 events_max=int(res["n_events"].max()) if res["n_events"].size else 0, steps=int(res["steps"]),
                 core_passes=int(res["n_core_passes_total"]), t_wall=time.time() - t0, rss_mb=rss_mb(),
                 f_capped=(a["E_capped"] / (a["E_carried_in"] + a["E_inj_new"])) if (a["E_carried_in"] + a["E_inj_new"]) > 0 else 0.0,
                 status="ok")
    return tally, pop_next


def analyse(out_dir, cfg=None, phot_bins=10, n_min_band=100):
    """Light curves from the escape records. Bolometric quantities on the
    slab grid; band photometry on `phot_bins` coarser time bins (groups of
    consecutive slabs) so that the red spectrum's optical bands hold enough
    packets; a band cell with fewer than `n_min_band` packets inside its
    passband is NaN and counted."""
    run = json.loads((out_dir / "run.json").read_text())
    cfg = run["config"]
    t_grid = np.array(run["t_grid"])
    n_slab = t_grid.size - 1
    groups = [g for g in np.array_split(np.arange(n_slab), min(phot_bins, n_slab)) if g.size]
    t_phot = np.array([t_grid[g[0]] for g in groups] + [t_grid[-1]])
    edges = phot.nu_edges(*LAM_WIN, N_SPEC); nu_c = np.sqrt(edges[1:] * edges[:-1])
    lo_w, hi_w = float(edges[0]), float(edges[-1])
    band_w = {b: (pb.bin_weights(edges) > 0) for b, pb in PASSBANDS.items()}
    legs = [l for l in cfg["legs"] if run["done"].get(l)]
    summary = dict(config=cfg, t_grid=t_grid.tolist(), t_phot=t_phot.tolist(), phot_bins=phot_bins, n_min_band=n_min_band,
                   legs={}, readings={})
    for leg in legs:
        ks = set(run["done"][leg])
        L_esc = np.full(n_slab, np.nan); L_win = np.full(n_slab, np.nan)
        L_obs = np.zeros(n_slab); n_esc = np.zeros(n_slab, int)
        mags = []; nan_bands = 0
        for g in groups:
            nu_g, e_g = [], []
            for k in g:
                f = out_dir / f"esc_{leg}_{k:03d}.npz"
                if k not in ks or not f.exists():
                    continue
                d = np.load(f)
                dt = t_grid[k + 1] - t_grid[k]
                L_esc[k] = d["e"].sum() / dt
                inw = (d["nu"] >= lo_w) & (d["nu"] < hi_w)
                L_win[k] = d["e"][inw].sum() / dt
                n_esc[k] = d["e"].size
                h, _ = np.histogram(d["t_obs"], bins=t_grid, weights=d["e"]); L_obs += h / np.diff(t_grid)
                nu_g.append(d["nu"][inw]); e_g.append(d["e"][inw])
            dt_g = t_grid[g[-1] + 1] - t_grid[g[0]]
            if nu_g:
                nu_g = np.concatenate(nu_g); e_g = np.concatenate(e_g)
                lnu = ts.lnu_absolute(nu_g, e_g, edges, dt_g)
                m = phot.magnitudes(nu_c, lnu, PASSBANDS, phot.D_40MPC, edges)
                cnt, _ = np.histogram(nu_g, bins=edges)
                for b in PASSBANDS:
                    if (cnt * band_w[b]).sum() < n_min_band:
                        m[b] = np.nan; nan_bands += 1
            else:
                m = {b: np.nan for b in phot.BANDS_PHOT}; nan_bands += len(phot.BANDS_PHOT)
            mags.append(m)
        tallies = run["tallies"][leg]
        W_tot = sum(tl["W"] for tl in tallies); E_rad = float(np.nansum(L_esc * np.diff(t_grid)))
        E_capped = sum(tl["E_capped"] for tl in tallies); E_core = sum(tl["E_core"] for tl in tallies)
        E_abs = sum(tl["E_abs"] for tl in tallies); E_inj = sum(tl["E_inj_new"] for tl in tallies)
        E_init = tallies[0]["E_carried_in"] if tallies else 0.0
        E_end = tallies[-1]["E_carried_out"] if tallies else 0.0
        closure = (E_init + E_inj - (E_rad + E_core + E_abs + E_capped + W_tot + E_end)) / max(E_init + E_inj, 1e-300)
        # peak by a 3-point parabola in (ln t, L)
        tm = np.sqrt(t_grid[1:] * t_grid[:-1]); ok = np.isfinite(L_esc)
        if ok.sum() >= 3:
            i = int(np.nanargmax(np.where(ok, L_esc, -np.inf)))
            i = min(max(i, 1), L_esc.size - 2)
            x = np.log(tm[i - 1:i + 2]); y = L_esc[i - 1:i + 2]
            if np.all(np.isfinite(y)) and y[1] >= y[0] and y[1] >= y[2] and (y[0] - 2 * y[1] + y[2]) < 0:
                a2, b2, c2 = np.polyfit(x, y, 2)
                xp = -b2 / (2 * a2); t_peak = float(np.exp(xp)); L_peak = float(a2 * xp * xp + b2 * xp + c2)
            else:
                t_peak, L_peak = float(tm[i]), float(L_esc[i])
        else:
            t_peak = L_peak = np.nan
        f_cap_max = max((tl["f_capped"] for tl in tallies), default=0.0)
        failed = [tl["k"] for tl in tallies if tl.get("status") != "ok"]
        summary["legs"][leg] = dict(L_esc=L_esc.tolist(), L_window=L_win.tolist(), L_obs=L_obs.tolist(), n_esc=n_esc.tolist(),
                                    mags=mags, t_peak_d=t_peak / DAY, L_peak=L_peak, E_rad=E_rad, W_tot=W_tot,
                                    E_init=E_init, E_inj=E_inj, E_core=E_core, E_abs=E_abs, E_capped=E_capped, E_end=E_end,
                                    closure=closure, f_capped_max=f_cap_max, failed_slabs=failed, nan_band_cells=nan_bands,
                                    t_obs_complete_until_d=float(t_grid[-1] * (1 - run.get("v_max_c", fontes.FONTES["v_max_c"])) / DAY),
                                    L_xkn=[(tl.get("L_thick", np.nan) + tl.get("L_thin_xkn", np.nan)) for tl in tallies],
                                    E_boundary_in=sum(tl.get("E_boundary_in", 0.0) for tl in tallies))
    # residuals vs the resolved leg of each redistribution
    for leg in legs:
        ref = REF_OF[leg]
        if ref == leg or ref not in summary["legs"]:
            continue
        A, R = summary["legs"][leg], summary["legs"][ref]
        dm = []; dcol = []
        for ma, mr in zip(A["mags"], R["mags"]):
            dm.append({b: (ma[b] - mr[b]) if (np.isfinite(ma[b]) and np.isfinite(mr[b])) else np.nan for b in ma})
            dcol.append({f"{b1}-{b2}": ((ma[b1] - ma[b2]) - (mr[b1] - mr[b2]))
                         if all(np.isfinite(v) for v in (ma[b1], ma[b2], mr[b1], mr[b2])) else np.nan for b1, b2 in phot.COLORS})
        allc = np.array([v for d in dcol for v in d.values()], float)
        A["dm_vs_ref"] = dm; A["dcolour_vs_ref"] = dcol
        A["max_abs_dcolour"] = float(np.nanmax(np.abs(allc))) if np.isfinite(allc).any() else np.nan
        A["dL_peak"] = (A["L_peak"] / R["L_peak"] - 1.0) if R["L_peak"] else np.nan
        A["dt_peak"] = (A["t_peak_d"] / R["t_peak_d"] - 1.0) if R["t_peak_d"] else np.nan
        A["dE_rad"] = (A["E_rad"] / R["E_rad"] - 1.0) if R["E_rad"] else np.nan
    # readings with gray outcomes
    def reading(trio, label):
        r, e, b = trio
        have = [l for l in trio if l in summary["legs"]]
        if len(have) < 3:
            return dict(outcome="GRAY", reason=f"legs missing: {[l for l in trio if l not in have]}")
        gray = [l for l in have if summary["legs"][l]["f_capped_max"] > 0.01 or summary["legs"][l]["failed_slabs"]]
        if gray:
            return dict(outcome="GRAY", reason=f"capped > 1 % or failed slabs in {gray}")
        Lp = {l: summary["legs"][l]["L_peak"] for l in trio}
        if not all(np.isfinite(v) for v in Lp.values()):
            return dict(outcome="GRAY", reason="no peak")
        spread = (max(Lp.values()) - min(Lp.values())) / Lp[r]
        # Poisson noise on L_peak from the escaped-packet count at the peak slab
        n_pk = {l: max(summary["legs"][l]["n_esc"]) for l in trio}
        noise = max(1.0 / np.sqrt(max(n, 1)) for n in n_pk.values())
        ordered = Lp[r] >= Lp[e] >= Lp[b]
        out = dict(spread=spread, ordered=bool(ordered), noise=noise, L_peak=Lp,
                   t_peak_d={l: summary["legs"][l]["t_peak_d"] for l in trio})
        if 3 * noise > spread and spread < 0.08:
            out["outcome"] = "GRAY"; out["reason"] = "peak differences within the Poisson noise"
        elif spread <= 0.08 and ordered:
            out["outcome"] = "GREEN"; out["reason"] = "peaks within 8 % and Sobolev > expansion > line-binned"
        elif spread <= 0.08:
            out["outcome"] = "GREEN-unordered"; out["reason"] = "peaks within 8 %, ordering differs"
        else:
            out["outcome"] = "RED"; out["reason"] = f"peak spread {spread:.2f} exceeds 8 %"
        return out
    summary["readings"]["R1_eps1_peaks"] = reading(("Rth", "Bth", "Bbinth"), "eps=1")
    summary["readings"]["R2_fluorescence_peaks"] = reading(("R2", "B2", "Bbin2"), "fluorescence")
    summary["readings"]["R3_max_colour_residual"] = {l: summary["legs"][l].get("max_abs_dcolour") for l in legs if l in summary["legs"] and "max_abs_dcolour" in summary["legs"][l]}
    write_json_atomic(out_dir / "summary.json", summary)
    try:
        _figure(out_dir, summary)
    except Exception as e:                                        # the figure is a convenience
        summary["figure_error"] = str(e)
    return summary


def _figure(out_dir, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = np.array(summary["t_grid"]); tm = np.sqrt(t[1:] * t[:-1]) / DAY
    tp = np.array(summary["t_phot"]); tpm = np.sqrt(tp[1:] * tp[:-1]) / DAY
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for leg, o in summary["legs"].items():
        ax[0].plot(tm, o["L_esc"], label=leg)
        m = [d["z"] for d in o["mags"]]
        ax[1].plot(tpm, m, "o-", label=leg)
        if "dcolour_vs_ref" in o:
            ax[2].plot(tpm, [d["i-J"] if "i-J" in d else np.nan for d in o["dcolour_vs_ref"]], "o-", label=leg)
    ax[0].set_yscale("log"); ax[0].set_xlabel("t [d]"); ax[0].set_ylabel("L_esc [erg/s]"); ax[0].legend(fontsize=7)
    ax[1].set_xlabel("t [d]"); ax[1].set_ylabel("z [AB mag at 40 Mpc]"); ax[1].invert_yaxis(); ax[1].legend(fontsize=7)
    ax[2].set_xlabel("t [d]"); ax[2].set_ylabel("Δ(i−J) vs resolved [mag]"); ax[2].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out_dir / "lightcurve.png", dpi=120); plt.close(fig)


def merge_runs(out_dir, sources):
    """Combine per-leg run directories (one process per leg) into `out_dir`:
    the esc files are copied, run.json is the union of done/tallies/T_rad
    (the configs must share the grid; per-leg n_scale and caps are kept as
    a table for the record)."""
    import shutil
    out_dir.mkdir(parents=True, exist_ok=True)
    merged = None
    for src in sources:
        src = Path(src); r = json.loads((src / "run.json").read_text())
        if merged is None:
            merged = dict(r); merged["legs_config"] = {}; merged["config"] = dict(r["config"]); merged["config"]["legs"] = []
            merged["done"] = {}; merged["tallies"] = {}; merged["T_rad"] = {}
        assert np.allclose(merged["t_grid"], r["t_grid"]), f"{src}: a different time grid"
        for leg in r["done"]:
            if not r["done"][leg]:
                continue
            merged["config"]["legs"].append(leg)
            merged["legs_config"][leg] = dict(n_scale=r["config"]["n_scale"].get(leg, 1.0), max_events=r["config"]["max_events"],
                                              temperature=r["config"]["temperature"], source=str(src))
            merged["done"][leg] = r["done"][leg]; merged["tallies"][leg] = r["tallies"][leg]; merged["T_rad"][leg] = r["T_rad"].get(leg, {})
            for k in r["done"][leg]:
                f = src / f"esc_{leg}_{k:03d}.npz"
                if f.exists():
                    shutil.copy(f, out_dir / f.name)
    write_json_atomic(out_dir / "run.json", merged)
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--merge", default=None, help="comma list of per-leg run directories to combine into --out before --analyse")
    ap.add_argument("--t0", type=float, default=4.0); ap.add_argument("--t1", type=float, default=16.0)
    ap.add_argument("--n-slabs", type=int, default=40)
    ap.add_argument("--n-shell", type=int, default=64)
    ap.add_argument("--transport", default=None, help="a-b shells (default 1 to the edge)")
    ap.add_argument("--core", default="reflect", help="absorb, reemit (lab-frame re-emission), reemit_cm (comoving-frame re-emission, boundary work booked), reflect"); ap.add_argument("--max-passes", type=int, default=10 ** 9)
    ap.add_argument("--legs", default="R2,B2,Bbin2,Rth,Bth,Bbinth")
    ap.add_argument("--n-scale", default="", help="per-leg packet multipliers, e.g. Bbinth=0.03,Bbin2=0.3")
    ap.add_argument("--n-init", type=int, default=200000); ap.add_argument("--n-heat", type=int, default=100000)
    ap.add_argument("--init", default="radiation", help="radiation | none")
    ap.add_argument("--temperature", default="prescribed", help="prescribed | radiation")
    ap.add_argument("--max-events", type=int, default=20000); ap.add_argument("--max-steps", type=int, default=300000)
    ap.add_argument("--wall-slab", type=float, default=7200.0)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--dataset", default=None)
    ap.add_argument("--tau-min", type=float, default=1e-3); ap.add_argument("--emis-cut", type=float, default=1e-6)
    ap.add_argument("--f-min", type=float, default=1e-3)
    ap.add_argument("--lam-transport", default="1000,128000")
    ap.add_argument("--resume", action="store_true"); ap.add_argument("--analyse", action="store_true")
    ap.add_argument("--phot-bins", type=int, default=10); ap.add_argument("--n-min-band", type=int, default=100)
    ap.add_argument("--max-slabs", type=int, default=None, help="stop after this many slabs (pilots)")
    ap.add_argument("--model", default="fontes", help="fontes | p1xkn (the P1 composition on the xkn secular structure)")
    ap.add_argument("--kappa", type=float, default=22.3, help="p1xkn: the grey opacity of the xkn photosphere (Tanaka Y_e = 0.2: 22.3)")
    ap.add_argument("--n-outer", type=int, default=24); ap.add_argument("--x-lo", type=float, default=0.75)
    ap.add_argument("--stages", default=None, help="ion stages (default: II,III for fontes, II for p1xkn)")
    a = ap.parse_args()
    out_dir = Path(a.out); out_dir.mkdir(parents=True, exist_ok=True)
    if a.merge:
        merge_runs(out_dir, a.merge.split(","))
        print(f"merged {a.merge} into {out_dir}")
    if a.analyse:
        s = analyse(out_dir, phot_bins=a.phot_bins, n_min_band=a.n_min_band)
        for k, v in s["readings"].items():
            print(k, v)
        for leg, o in s["legs"].items():
            print(f"{leg:6s} t_peak {o['t_peak_d']:.2f} d  L_peak {o['L_peak']:.3e}  E_rad {o['E_rad']:.3e}  W {o['W_tot']:.3e}  closure {o['closure']:+.2e}  f_capped_max {o['f_capped_max']:.3f}"
                  + (f"  dL_peak {o['dL_peak']:+.3f} dt_peak {o['dt_peak']:+.3f} max|dcolour| {o['max_abs_dcolour']:.2f}" if "dL_peak" in o else ""))
        return
    legs = a.legs.split(",")
    n_scale = {l: 1.0 for l in legs}
    for kv in filter(None, a.n_scale.split(",")):
        k_, v_ = kv.split("="); n_scale[k_] = float(v_)
    stages = tuple(a.stages.split(",")) if a.stages else (("II", "III") if a.model == "fontes" else ("II",))
    cfg = dict(model=a.model, kappa=a.kappa, n_outer=a.n_outer, x_lo=a.x_lo, stages=list(stages),
               t0=a.t0, t1=a.t1, n_slabs=a.n_slabs, n_shell=a.n_shell, transport=a.transport, core=a.core, max_passes=a.max_passes,
               legs=legs, n_scale=n_scale, n_init=a.n_init, n_heat=a.n_heat, init=a.init, temperature=a.temperature,
               max_events=a.max_events, max_steps=a.max_steps, wall_slab=a.wall_slab, seed=a.seed, dataset=a.dataset,
               tau_min=a.tau_min, emis_cut=None if a.emis_cut <= 0 else a.emis_cut, f_min=None if a.f_min <= 0 else a.f_min,
               lam_transport=[float(x) for x in a.lam_transport.split(",")], git=L.git_sha())
    t_grid = time_grid(a.t0, a.t1, a.n_slabs)
    run_path = out_dir / "run.json"
    if a.resume and run_path.exists():
        run = json.loads(run_path.read_text())
        if run["config_hash"] != config_hash({k: v for k, v in cfg.items() if k != "git"}):
            raise SystemExit("resume: the configuration differs from run.json")
    else:
        run = dict(config=cfg, config_hash=config_hash({k: v for k, v in cfg.items() if k != "git"}), t_grid=t_grid.tolist(),
                   done={l: [] for l in legs}, tallies={l: [] for l in legs}, T_rad={l: {} for l in legs}, notes=[])
    st0 = build_state(cfg, float(t_grid[0]))
    shells = zone_of(cfg, st0)
    m_all = st0.shell_mass(); m_sh = m_all[shells]
    run["shells"] = shells; run["mass_dropped_frac"] = float(1.0 - m_sh.sum() / m_all.sum())
    run["v_max_c"] = float(st0.meta.get("v_max_c", fontes.FONTES["v_max_c"]))
    # the per-slab injected energies for the packet allocation (the p1xkn zone
    # grows with time; the allocation uses the t0 zone's heating plus the
    # boundary luminosity, the exact per-slab values are recomputed in the loop)
    E_heat = np.array([sum(x for x in (slab_sources(cfg, st0, shells, t_grid[k], t_grid[k + 1])[0],
                                        slab_sources(cfg, st0, shells, t_grid[k], t_grid[k + 1])[1].sum()))
                       for k in range(a.n_slabs)])
    run["E_heat_per_slab"] = E_heat.tolist()
    lo, hi = (float(x) for x in phot.nu_edges(*cfg["lam_transport"], 1))
    n_last = a.n_slabs if a.max_slabs is None else min(a.n_slabs, a.max_slabs)
    for k in range(n_last):
        t_a, t_b = float(t_grid[k]), float(t_grid[k + 1])
        pending = [l for l in legs if k not in run["done"][l]]
        if not pending:
            continue
        st = build_state(cfg, t_a) if cfg["temperature"] == "prescribed" or k == 0 else None
        if st is not None:
            shells = zone_of(cfg, st)
        atom_shared = build_atom(cfg, st, shells) if st is not None and cfg["temperature"] == "prescribed" else None
        sources = slab_sources(cfg, st if st is not None else st0, shells, t_a, t_b)
        if cfg.get("model", "fontes") == "p1xkn":
            st_b = build_state(cfg, t_b); shells_next = zone_of(cfg, st_b)
        else:
            shells_next = shells
        e_next = st0.v_edges[shells_next[0]:shells_next[-1] + 2] * t_b
        for leg in pending:
            if cfg["temperature"] == "radiation" and k > 0:
                T_prev = run["T_rad"][leg].get(str(k - 1))
                T_over = None
                if T_prev is not None:
                    T_over = np.array(build_state(cfg, t_a).T_gas)
                    Tp = np.array(T_prev, float); ok = np.isfinite(Tp)
                    T_over[np.array(shells)[ok]] = Tp[ok]
                st_leg = build_state(cfg, t_a, T_override=T_over); atom = build_atom(cfg, st_leg, shells)
            else:
                st_leg = st if st is not None else build_state(cfg, t_a)
                atom = atom_shared if atom_shared is not None else build_atom(cfg, st_leg, shells)
            pop_path = out_dir / f"pop_{leg}.npz"
            if k == 0:
                if cfg["init"] == "radiation":
                    rng0 = np.random.default_rng(seed_for(cfg["seed"], leg, 10 ** 6))
                    pop, E0, _ = ts.initial_radiation(st_leg, int(round(cfg["n_init"] * n_scale[leg])), rng0, lo, hi, shells)
                else:
                    pop = ts.Population.empty()
            else:
                if not pop_path.exists():
                    run["notes"].append(f"{leg} slab {k}: no checkpoint (an earlier slab failed); leg stopped")
                    print(f"  {leg} slab {k}: no checkpoint, leg stopped (gray)", flush=True)
                    write_json_atomic(run_path, run)
                    continue
                pop = ts.Population.from_npz(pop_path)
            n_heat_k = int(round(cfg["n_heat"] * n_scale[leg] * E_heat[k] / E_heat.sum()))
            try:
                tally, pop_next = run_slab(cfg, st_leg, atom, shells, leg, k, pop, sources, n_heat_k, t_a, t_b, out_dir, e_next)
            except RuntimeError as e:
                run["tallies"][leg].append(dict(k=k, status=f"failed: {str(e)[:120]}", t_a=t_a, t_b=t_b, f_capped=0.0, W=0.0, E_capped=0.0,
                                                E_core=0.0, E_abs=0.0, E_inj_new=0.0, E_carried_in=0.0, E_carried_out=0.0))
                run["notes"].append(f"{leg} slab {k} failed: {e}")
                print(f"  {leg} slab {k}: FAILED {str(e)[:100]} -- the leg stops here (gray)", flush=True)
                write_json_atomic(run_path, run)
                break                                   # the next slab would resume a stale checkpoint
            if pop_path.exists():
                os.replace(pop_path, out_dir / f"pop_{leg}.prev.npz")
            pop_next.to_npz(pop_path)
            if cfg["temperature"] == "radiation":
                T_rad, n_sh_pop = ts.t_rad_from_population(pop_next, build_state(cfg, t_b), shells)
                run["T_rad"][leg][str(k)] = [None if not np.isfinite(x) else float(x) for x in T_rad]
                tally["n_T_fallback"] = int(np.sum(~np.isfinite(T_rad)))
            run["tallies"][leg].append(tally); run["done"][leg].append(k)
            write_json_atomic(run_path, run)
            print(f"slab {k:3d} {t_a/DAY:6.2f}-{t_b/DAY:6.2f} d {leg:6s} n_in {pop.n:7d} n_new {n_heat_k:6d} esc {tally['E_esc']:.3e} "
                  f"W {tally['W']:.3e} carried {tally['E_carried_out']:.3e} ev {tally['events_mean']:.0f}/{tally['events_max']} "
                  f"capped {tally['f_capped']:.4f} {tally['t_wall']:.0f}s rss {tally['rss_mb']:.0f}MB", flush=True)
    print("done; analyse with --analyse")


if __name__ == "__main__":
    main()
