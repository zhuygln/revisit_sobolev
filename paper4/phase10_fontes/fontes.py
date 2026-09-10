"""Paper IV, Phase 10 (the pivot): the Fontes et al. (2020, MNRAS 493, 4143)
simplified problem as a snapshot, under complete thermal redistribution and
under energy-conserving fluorescence, with the three opacity treatments.

Their Appendix C: rho(r, t) = rho0 (t/t0)^-3 (1 - r^2/(v_max t)^2)^3,
T(r, t) = T0 (t/t0)^-1 (1 - r^2/(v_max t)^2), v_max = 0.25 c, t0 = 4 d,
M = 1.4e-2 Msun -> rho0 = 2.51e-15 g cm^-3, T0 = 5.70e3 K, 100 % Nd,
bound-bound only, lines with f < 1e-3 dropped, 64 uniform cells in v from
0 to v_max, heating eps = eps0 (t/t0)^-1.3 uniform in mass; direct Sobolev
vs expansion vs line-binned, all with eps = 1 (complete thermal
redistribution). Their result: the three light-curve peaks within 8 %,
expansion brighter than line-binned by <= 5 % (Section 4.2.1).

What this snapshot keeps and what it changes (each stated in the output):
  * the profile, mass, composition, v_max, t0 and the f cut as written;
  * T(r) from their INITIAL profile at the chosen epoch (their code evolves
    T; ours is a fixed-state transport), LTE Saha ionization of Nd at that
    T and rho (Nd II and III carry the GSI lines; Nd I and IV have no data
    here), Boltzmann populations;
  * packets injected in the volume in proportion to mass (uniform specific
    heating) with the local Planck spectrum, in one snapshot; the observable
    is the escaping spectrum and its band magnitudes and the escaped energy
    fraction, not a light curve;
  * the innermost cell edge is a tiny re-emitting core (v_max/64);
  * classical (non-worldline) transport for every leg so that the thermal
    legs, which the zoned transport carries classically, and the macroatom
    legs are compared on the same footing.
Legs: Rth / Bth / Bbinth (eps = 1) and R2 / B2 / Bbin2 (downward
macroatom) -- the same three treatments under the two redistributions.
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
          ROOT / "paper3/phase12_grid", ROOT / "paper4/phase2_energy", ROOT / "paper4/phase1_benchmarks",
          ROOT / "paper4/phase7_ionization", ROOT / "paper4/phase8_shells"):
    sys.path.insert(0, str(p))

from sobolev.constants import C                                   # noqa: E402
from sobolev.source import MSUN                                   # noqa: E402
from sobolev import ejecta as ej                                  # noqa: E402
from sobolev import photometry as phot                            # noqa: E402
from sobolev.zoned_atom import ZonedAtom                          # noqa: E402
from sobolev.energy_packets import energy_accounting              # noqa: E402
from forest_mc import run_mc                                      # noqa: E402
from observables import observe, LAM_WIN, N_SPEC                  # noqa: E402
from grid import photometer, rss_mb                               # noqa: E402
import legs as L                                                  # noqa: E402
import saha_states                                                # noqa: E402
from shells import band_weights_by_shell                          # noqa: E402

DAY = 86400.0
FONTES = dict(name="Fontes2020C", m_msun=1.4e-2, v_max_c=0.25, t0_d=4.0, rho0=2.51e-15, T0=5700.0,
              eps0=8.2e8, eps_index=-1.3, composition="Nd", f_min=1e-3, n_cell=64,
              source="Fontes, Fryer, Hungerford, Wollaeger & Korobkin 2020, MNRAS 493, 4143, Appendix C")
LEGS = ("Rth", "Bth", "Bbinth", "R2", "B2", "Bbin2")


def build_fontes(t_d=4.0, n_shell=64, saha=True):
    """The Appendix C ejecta at epoch t_d as an EjectaState (uniform shells in
    v from 0 to v_max; the profile evaluated as cell volume averages for rho
    and at cell centres for T, as SuperNu does), pure Nd, LTE Saha."""
    t = t_d * DAY
    v_max = FONTES["v_max_c"] * C
    scale_rho = (t / (FONTES["t0_d"] * DAY)) ** -3
    scale_T = (t / (FONTES["t0_d"] * DAY)) ** -1
    v_edges = np.linspace(0.0, v_max, n_shell + 1)
    rho = np.empty(n_shell); T = np.empty(n_shell)
    for i in range(n_shell):
        vv = np.linspace(v_edges[i], v_edges[i + 1], 401)
        f = (1.0 - (vv / v_max) ** 2) ** 3
        rho[i] = FONTES["rho0"] * scale_rho * np.trapezoid(f * vv ** 2, vv) / np.trapezoid(vv ** 2, vv)
        vc = 0.5 * (v_edges[i] + v_edges[i + 1])
        T[i] = FONTES["T0"] * scale_T * (1.0 - (vc / v_max) ** 2)
    rho = np.maximum(rho, 1e-30); T = np.maximum(T, 300.0)
    X = {"Nd": np.ones(n_shell)}
    f_ion = {"Nd II": np.ones(n_shell), "Nd III": np.zeros(n_shell)}
    st = ej.EjectaState(t=t, v_edges=v_edges, rho=rho, T_gas=T, T_rad=T.copy(), X=X, f_ion=f_ion, Y_e=None,
                        meta=dict(FONTES, t_d=t_d, n_shell=n_shell, photospheric_shell=0, mass_msun=None))
    st.meta["mass_msun"] = st.mass() / MSUN
    if saha:
        st = saha_states.ionize(st)
    st.check()
    return st


def run_fontes(state, shells, n, legs=LEGS, seeds=L.SEEDS, tau_min=1e-3, emis_cut=1e-6, f_min=1e-3,
               budget_s=None, verbose=True, dataset=None, core="reemit"):
    t0 = time.time()
    zone = state.transport_zone(shells)
    m = state.shell_mass()
    core_frac = float(m[:shells[0]].sum() / m.sum())      # the interior's share of a uniform heating rate
    atom = ZonedAtom.from_state(state, shells, stages=("II", "III"), tau_min=tau_min, emis_cut=emis_cut, f_min=f_min,
                                dataset=dataset)
    t_build = time.time() - t0
    lo, hi = (float(x) for x in phot.nu_edges(*LAM_WIN, 1))
    edges = phot.nu_edges(*LAM_WIN, N_SPEC); nu_c = np.sqrt(edges[1:] * edges[:-1])
    l_core = phot.planck_luminosity(lo, hi, zone["r_core"], zone["t_core"])
    row = dict(state=state.meta["name"], t_d=state.t / DAY, shells=list(map(int, shells)), n_shell=len(shells),
               zone={k: v for k, v in zone.items()}, n=n, seeds=list(seeds), tau_min=tau_min, emis_cut=emis_cut,
               f_min=f_min, dataset=dataset or "gsi", relativity=None, launch="volume", launch_core_frac=core_frac, core=core,
               n_lines=int(atom.n_lines_total),
               n_opacity_union=int(atom.n_opacity), n_opacity_shell=atom.n_opacity_shell.tolist(),
               ions=atom.ions, t_build=t_build, rss_mb_atom=rss_mb(), lam_window=list(LAM_WIN), n_spec=N_SPEC,
               f_ion_NdIII=[float(state.f_ion["Nd III"][s]) for s in shells], T=[float(state.T_gas[s]) for s in shells],
               legs={}, timing={}, git=L.git_sha())
    if verbose:
        print(f"{row['state']} t={row['t_d']:g} d, {len(shells)} shells: {atom.ions}, {atom.n_lines_total} lines (f > {f_min}), "
              f"union {atom.n_opacity}, build {t_build:.0f}s, rss {rss_mb():.0f} MB", flush=True)
    for tag in legs:
        spec = L.LEGS[tag]
        tl = time.time()
        res = [run_mc(atom, zone["r_core"], zone["r_out"], zone["t_exp"], lo, hi, n, spec["mode"], seed=s,
                      t_core=float(state.T_gas[shells[0]]), relativity=None, max_steps=L.MAX_STEPS, packets="energy",
                      launch_weight="energy", launch="volume", launch_core_frac=core_frac, core=core, wall_s=budget_s) for s in seeds]
        # "absorbing" = no synthetic renormalisation: the escaping spectrum itself; legs share E_inj
        o = photometer(observe(res, l_core, "absorbing"), edges, nu_c, phot.D_40MPC)
        per_seed = [photometer(observe([r], l_core, "absorbing"), edges, nu_c, phot.D_40MPC)["mags"] for r in res]
        o["mags_seed_std"] = {b: float(np.std([m[b] for m in per_seed], ddof=1)) if len(res) > 1 else np.nan for b in o["mags"]}
        acc = [energy_accounting(r) for r in res]
        o["energy"] = {k: float(np.mean([a[k] for a in acc])) for k in acc[0] if k != "packets"}
        o["ledger"] = L.ledger(res)
        o["E_esc_frac"] = float(np.sum([r["accounting"]["E_esc"] for r in res]) / np.sum([r["accounting"]["E_inj"] for r in res]))
        o["events_per_packet"] = float(np.mean([r["n_events"].mean() for r in res]))
        o["n_events_shell"] = np.sum([r["n_events_shell"] for r in res], axis=0).tolist()
        o["n_core_passes"] = int(np.sum([r["n_core_passes_total"] for r in res]))
        o["band_weights_by_shell"] = band_weights_by_shell(res, edges, len(shells))
        o["mode"] = spec["mode"]; o["t_wall"] = time.time() - tl
        row["legs"][tag] = o
        if verbose:
            print(f"  {tag:7s} {spec['mode']:20s} {o['t_wall']:6.1f}s ev/pkt={o['events_per_packet']:6.2f} "
                  f"esc={o['E_esc_frac']:.4f} core_passes={o['n_core_passes']} "
                  + " ".join(f"{b}={o['mags'][b]:.2f}" for b in ("g", "z", "K")), flush=True)
    ref = {}
    for tag in legs:
        r_tag = "Rth" if tag.endswith("th") else "R2"
        if r_tag in row["legs"] and tag != r_tag:
            row["legs"][tag]["dm_vs_ref"] = phot.delta_mag(row["legs"][tag]["mags"], row["legs"][r_tag]["mags"])
            row["legs"][tag]["dE_esc_vs_ref"] = row["legs"][tag]["E_esc_frac"] / row["legs"][r_tag]["E_esc_frac"] - 1.0
            row["legs"][tag]["ref"] = r_tag
    row["t_wall"] = time.time() - t0; row["rss_mb"] = rss_mb()
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t", type=float, default=4.0, help="epoch in days (4.0 = t0; 6.3 = their peak)")
    ap.add_argument("--n-shell", type=int, default=64)
    ap.add_argument("--transport", default=None, help="a-b shells to transport (default: 1 to the edge)")
    ap.add_argument("--legs", default=",".join(LEGS))
    ap.add_argument("--n", type=int, default=300000)
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--tau-min", type=float, default=1e-3)
    ap.add_argument("--emis-cut", type=float, default=1e-6)
    ap.add_argument("--f-min", type=float, default=1e-3)
    ap.add_argument("--budget", type=float, default=None)
    ap.add_argument("--no-saha", action="store_true")
    ap.add_argument("--dataset", default=None, help="None = GSI cache; 'jplt' = the Japan-Lithuania line list")
    ap.add_argument("--core", default="reemit", help="inner boundary: reemit (thermalising interior) or reflect (mirror)")
    ap.add_argument("--state-out", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    st = build_fontes(a.t, a.n_shell, saha=not a.no_saha)
    if a.state_out:
        st.to_json(a.state_out)
    if a.transport:
        lo_, hi_ = (int(x) for x in a.transport.split("-")); shells = list(range(lo_, hi_ + 1))
    else:
        shells = list(range(1, st.n_shell))
    row = run_fontes(st, shells, a.n, tuple(a.legs.split(",")), tuple(int(s) for s in a.seeds.split(",")), a.tau_min,
                     None if a.emis_cut <= 0 else a.emis_cut, None if a.f_min <= 0 else a.f_min, a.budget,
                     dataset=a.dataset, core=a.core)
    out = a.out or (HERE / f"fontes_t{a.t:g}_n{a.n_shell}{'' if a.dataset is None else '_' + a.dataset}_{a.core}.json")
    Path(out).write_text(json.dumps(row, indent=1, default=float))
    print(f"wrote {out} in {row['t_wall']:.0f}s")


if __name__ == "__main__":
    main()
