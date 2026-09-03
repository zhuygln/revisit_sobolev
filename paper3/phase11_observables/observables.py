"""Paper III §10: what the cancellation boundary does to an observer.

F40 put a homologous Ce II ejecta history through the closure's sign change --
+64.6% too bright at 0.5 d, zero at 1.17 d, -28.4% too opaque at 1.5 d. Every
number in it is a transport residual in one 155 A diagnostic band. Nobody
observes a 155 A band ratio.

    Does the apparent agreement at the crossing epoch coexist with
    significant colour and spectral errors before and after it?

That is the question this driver answers, in magnitudes. `run_mc` already
returns the escaping packet list; every driver in the repo, `trajectory.py`
included, reduced it to a band ratio and discarded the spectrum.
`sobolev/photometry.py` turns it into absolute L_nu, L_bol and AB magnitudes.

Three departures from `trajectory.py`, all deliberate:

  * **Fixed launch window** (`LAM_WIN`), identical at every epoch and for every
    ion. `trajectory.py:82` sets it from each atom's own opacity extent -- for
    Ce II that is 1128 A to 36.8 um, so most packets land in a far-IR tail no
    filter samples and no two epochs are directly comparable. The expansion bin
    grid is unaffected: `forest_mc.py:471` clamps it to cover the full opacity
    range regardless of the launch window.
  * **Cooling core** (`--core cool`), t_core = T_gas(t). `trajectory.py` holds
    6000 K at every epoch, which freezes the injected continuum's shape; run
    `--core fixed` as the control that isolates the opacity contribution.
  * **The crossing is computed and stored.** F40's t = 1.17 d and S = 47.5 were
    interpolated by hand and persisted nowhere.

What this is not: a kilonova light curve. The source is an imposed blackbody
core, there is no radioactive heating and no energy equation, so L_bol(t)'s
shape is the core's. Only DIFFERENCES between treatments on the same source and
the same ejecta are claimed. Bandpasses are top-hats (see photometry.py).

Usage:
  python observables.py [--ion 58CeII|57LaII|blend] [--core cool|fixed] [--n N]
                        [--v f40|kn] [--relativity worldline]
                        [--mass M_ej_Msun] [--xlan X] [--epochs 2,3,4,...]
                        [--label TAG]

  # §4.37 primary: F40's history in magnitudes
  python observables.py --ion 58CeII --core cool --n 1000000

  # the physically normalized lanthanide-poor kilonova of §4.37.8
  python observables.py --ion 57LaII --core cool --v kn --relativity worldline \
      --mass 0.01 --xlan 1e-3 --epochs 2,3,4,6,8,12 --label blue --n 250000

Reduce with `summary.py`, which masks bands the reference barely emits in.
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paper3",
          ROOT / "paper3/phase0_reference", ROOT / "paper3/phase10_kilonova"):
    sys.path.insert(0, str(p))
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import R_CORE, R_OUT, T_EXP, nu_of
from redistribution import RedistributionKernel
from reference import SEEDS
from trajectory import A_OF, DATA, EPOCHS, ION_FRAC, RHO_1D, T_GAS_1D, X_LAN, state

from sobolev.constants import C
from sobolev.forest_stats import band_saturation
from sobolev.normalization import from_conditions
from sobolev import photometry as phot

# The launch window: fixed for every epoch, ion and leg. Its BLUE edge matters
# physically, not just for bookkeeping -- most lanthanide opacity sits in the
# UV, and it is UV photons fluorescing redward that refill the optical bands, so
# truncating it changes the band-3800 residual (a 2000 A cut moves F40's 0.5 d
# value from +64.6% to +54.8%). 1000 A is below every GSI II line, so nothing is
# lost; the RED edge at 3 um only avoids spending packets on a far-IR tail that
# no filter samples and that carries almost no Planck energy below 7000 K.
DAY_S = 86400.0
LAM_WIN = (1000.0, 30000.0)
N_SPEC = 200                    # spectral bins across the window
BAND3800 = (3800.0, 3955.0)     # the F40 diagnostic band, for continuity
BLEND = ("57LaII", "58CeII", "59PrII", "60NdII")

LEGS = (("A_redist", "sobolev_group"), ("B_opacity", "expansion_branch"),
        ("C_both", "expansion_group"), ("C_binned", "binned_group"))

# Homologous velocities of the shell's inner and outer edges, in units of c.
# `trajectory.py` inherits Paper I's shell, 1000-3000 km/s -- kilonova DENSITIES
# and temperatures at supernova velocities. Sobolev tau does not care (tau ~ n t
# is velocity-free in homologous flow), but the wavelength interval a packet
# sweeps does: dlam/lam = dv/c is 0.0067 here against ~0.2 for real ejecta, so a
# packet meets ~30x fewer lines than it would in a kilonova. Since that interval
# is precisely what a grouped closure coarse-grains, the F40 geometry is a
# LOWER bound on the closure error. `--v` runs the kilonova case.
V_F40 = (R_CORE / T_EXP / C, R_OUT / T_EXP / C)   # 0.00334 - 0.0100 c
V_KN = (0.05, 0.20)                              # kilonova ejecta

# `trajectory.py` documents RHO_1D = 2e-17 as chosen so Ce II crosses inside the
# observable window. At kilonova velocities that is not an ejecta mass: a
# uniform sphere out to 0.2c holds only 5.9e-6 Msun at that density. But the
# right comparison is n_ion, not rho, and there the tuned value turns out to be
# almost exactly a LANTHANIDE-POOR component's:
#
#   ejecta                                   rho(1 d)   n_ion(Ce II, 1 d)
#   F40 tuned                                2.0e-17    8.6e3
#   blue  M=0.01 Msun v=0.3c X_lan=1e-3      1.0e-14    4.4e4    (5x F40)
#   red   M=0.03 Msun v=0.2c X_lan=0.1       1.0e-13    4.4e7    (5100x F40)
#
# So the boundary lies in the blue component's range and ~3.7 decades below the
# red component's. `--mass` derives rho(1 d) from (M_ej, v_max) instead of
# taking it as given, so that claim is run rather than asserted.
MSUN = 1.989e33


def rho_1d_from_mass(m_ej_msun, v_max_c):
    """rho at 1 d for a uniform sphere of mass m_ej out to v_max."""
    r = v_max_c * C * DAY_S
    return m_ej_msun * MSUN / ((4.0 / 3.0) * np.pi * r**3)


def build_atom(ion, st, x_lan=X_LAN):
    """The epoch's atom under the ASTROPHYSICAL normalization, plus its n_ion.

    A blend splits x_lan equally between its ions, so the total lanthanide mass
    fraction matches the single-ion runs and only the composition differs.
    """
    if ion == "blend":
        x = x_lan / len(BLEND)
        specs, n_ions = [], {}
        for s in BLEND:
            lev, tr = DATA / f"{s}_levels_calib.txt", DATA / f"{s}_transitions_calib.txt"
            n, _ = from_conditions(lev, tr, st["T_gas"], st["t_exp"], st["rho"],
                                   x, ION_FRAC, A_OF[s])
            specs.append((lev, tr, n)); n_ions[s] = float(n)
        return ForestAtom.from_gsi_blend(specs, st["T_gas"], st["t_exp"],
                                         tau_min=1e-3), n_ions
    lev, tr = DATA / f"{ion}_levels_calib.txt", DATA / f"{ion}_transitions_calib.txt"
    n, _ = from_conditions(lev, tr, st["T_gas"], st["t_exp"], st["rho"],
                           x_lan, ION_FRAC, A_OF[ion])
    return ForestAtom.from_gsi(lev, tr, st["T_gas"], n, st["t_exp"],
                               tau_min=1e-3), float(n)


def observe(res_list, l_core, core="absorbing"):
    """Seed-averaged L_nu, L_bol, magnitudes and colours for one leg.

    `core` selects the inner boundary the spectrum is normalized to
    (`photometry._scale`); F40/F41 used the absorbing core, the §4.39 grid the
    equilibrium one. `f_return` and `n_trapped` are diagnostics of the run.
    """
    edges = phot.nu_edges(*LAM_WIN, N_SPEC)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    lnu = np.mean([phot.emergent_lnu(r, edges, l_core, core) for r in res_list], axis=0)
    lbol = float(np.mean([phot.bolometric(r, l_core, core) for r in res_list]))
    mags = phot.magnitudes(nu_c, lnu)
    return {"L_nu": lnu.tolist(), "L_bol": lbol, "mags": mags,
            "colors": phot.colors(mags), "core": core,
            "f_return": float(np.mean([phot.return_fraction(r) for r in res_list])),
            "f_dep": float(np.mean([phot.deposited_fraction(r) for r in res_list])),
            "n_trapped": int(sum(r.get("n_trapped", 0) for r in res_list)),
            "L_bol_absorbing": float(np.mean([phot.bolometric(r, l_core) for r in res_list])),
            "b3800": float(np.mean([band_ratio(r, *nu_of(*BAND3800),
                                               weight="energy")[0] for r in res_list]))}


def epoch(ion, t_d, n, core_law, ng=32, vel=V_F40, relativity=None,
          rho_1d=RHO_1D, x_lan=X_LAN):
    st = state(t_d, core_law)
    st["r_core"] = vel[0] * C * st["t_exp"]
    st["r_out"] = vel[1] * C * st["t_exp"]
    st["rho"] = rho_1d * t_d ** -3
    atom, n_ion = build_atom(ion, st, x_lan)
    if atom.n_opacity < 10:
        return {"t_d": t_d, "skipped": f"{atom.n_opacity} opacity lines"}
    lo, hi = phot.nu_edges(*LAM_WIN, 1)          # the fixed launch window, Hz
    lo, hi = float(lo), float(hi)
    l_core = phot.planck_luminosity(lo, hi, st["r_core"], st["t_core"])

    ref_res, ei, eo = [], [], []
    for s in SEEDS:
        r = run_mc(atom, st["r_core"], st["r_out"], st["t_exp"], lo, hi, n,
                   "sobolev_branch", seed=s, t_core=st["t_core"],
                   relativity=relativity, collect_events=True)
        ref_res.append(r); e = r["events"]; ei.append(e[0]); eo.append(e[1])
    ref = observe(ref_res, l_core)

    # the kernel keeps trajectory.py's per-epoch grid on the atom's own
    # opacity extent -- it is the closure's group structure, not a reporting
    # grid, and each epoch is an independent snapshot
    k_lo, k_hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    nu_in, nu_out = np.concatenate(ei), np.concatenate(eo)
    kern = RedistributionKernel.from_branching_mc(nu_in, nu_out, np.ones(nu_in.size),
                                                  ng, nu_lo=k_lo, nu_hi=k_hi)

    legs = {}
    for tag, mode in LEGS:
        rows = []
        for s in SEEDS:
            kw = {"kernel": kern} if mode.endswith("_group") else {}
            rows.append(run_mc(atom, st["r_core"], st["r_out"], st["t_exp"], lo, hi,
                               n, mode, seed=s, t_core=st["t_core"],
                               relativity=relativity, **kw))
        o = observe(rows, l_core)
        o["dm"] = phot.delta_mag(o["mags"], ref["mags"])
        o["dcolor"] = {k: o["colors"][k] - ref["colors"][k] for k in ref["colors"]}
        o["dm_bol"] = phot.bol_delta_mag(o["L_bol"], ref["L_bol"])
        o["dF_b3800"] = o["b3800"] / ref["b3800"] - 1 if ref["b3800"] > 1e-3 else None
        legs[tag] = o

    nb = (C / (BAND3800[1] * 1e-8), C / (BAND3800[0] * 1e-8))
    out = {"t_d": t_d, "n_ion": n_ion, "T_gas": st["T_gas"], "t_core": st["t_core"],
           "v_core": vel[0], "v_out": vel[1], "x_lan": x_lan,
           "rho": st["rho"], "r_core": st["r_core"], "L_core_window": l_core,
           "n_opacity": int(atom.n_opacity), "tau_max": float(atom.op_tau.max()),
           "ref": ref, "legs": legs}
    out.update({f"band_{k}": v for k, v in band_saturation(atom, *nb).items()})
    return out


def crossing_epoch(rows, tag, key="dF_b3800"):
    """Epoch and band saturation at which a leg's signed error crosses zero.

    F40 quoted t = 1.17 d and S = 47.5 from a hand interpolation that no code
    performs and no file stores. Linear in t (the sampling variable), log in S
    (which spans four decades). Detects both directions, as
    `paper3/synthetic/boundary.py:72` does -- the first version of that helper
    tested one direction only and silently missed half the crossings.
    """
    r = [x for x in rows if "skipped" not in x
         and x["legs"][tag].get(key) is not None]
    r.sort(key=lambda x: x["t_d"])
    for a, b in zip(r, r[1:]):
        va, vb = a["legs"][tag][key], b["legs"][tag][key]
        if va * vb >= 0:
            continue
        f = -va / (vb - va)
        t = a["t_d"] + f * (b["t_d"] - a["t_d"])
        Sa, Sb = a["band_S_band"], b["band_S_band"]
        S = (float(np.exp(np.log(Sa) + f * (np.log(Sb) - np.log(Sa))))
             if Sa > 0 and Sb > 0 else None)
        return {"t_d": float(t), "S_band": S,
                "direction": "neg_to_pos" if va < 0 else "pos_to_neg"}
    return None


def main(ion, n, core_law, out_name, vel=V_F40, relativity=None,
         rho_1d=RHO_1D, x_lan=X_LAN, epochs=EPOCHS):
    print(f"{ion}: rho(1d) = {rho_1d:.1e} g/cm3, X_lan = {x_lan}, "
          f"T_gas(1d) = {T_GAS_1D:.0f} K, core = {core_law}, "
          f"v = {vel[0]:.4g}-{vel[1]:.4g} c, relativity = {relativity}, "
          f"window {LAM_WIN[0]:.0f}-{LAM_WIN[1]:.0f} A, {n} packets x {len(SEEDS)} seeds",
          flush=True)
    hdr = (f"{'t/d':>5s}{'T_gas':>7s}{'T_core':>7s}{'S_band':>9s}{'dF3800':>9s}"
           f"{'dm_bol':>8s}{'dm_g':>8s}{'dm_r':>8s}{'dm_J':>8s}"
           f"{'d(g-r)':>8s}{'d(i-J)':>8s}")
    print(hdr, flush=True)
    rows, t0 = [], time.time()
    for t_d in epochs:
        r = epoch(ion, t_d, n, core_law, vel=vel, relativity=relativity,
                  rho_1d=rho_1d, x_lan=x_lan)
        rows.append(r)
        if "skipped" in r:
            print(f"{t_d:5.2f}  skipped -- {r['skipped']}", flush=True); continue
        L = r["legs"]["C_both"]
        def g(d, k, w=8, p=3):
            v = d.get(k)
            return " " * (w - 4) + "--  " if v is None or not np.isfinite(v) \
                else f"{v:+{w}.{p}f}"
        print(f"{t_d:5.2f}{r['T_gas']:7.0f}{r['t_core']:7.0f}{r['band_S_band']:9.1f}"
              f"{g({'d': None if L['dF_b3800'] is None else 100*L['dF_b3800']}, 'd', 8, 1)}%"
              f"{g(L,'dm_bol')}{g(L['dm'],'g')}"
              f"{g(L['dm'],'r')}{g(L['dm'],'J')}{g(L['dcolor'],'g-r')}"
              f"{g(L['dcolor'],'i-J')}", flush=True)
    cross = {tag: crossing_epoch(rows, tag) for tag, _ in LEGS}
    for tag, c in cross.items():
        if c:
            print(f"  {tag}: band-3800 error crosses zero at t = {c['t_d']:.2f} d, "
                  f"S = {c['S_band']:.1f} ({c['direction']})", flush=True)
    (HERE / out_name).write_text(json.dumps(
        {"ion": ion, "n": n, "core_law": core_law, "seeds": list(SEEDS),
         "v_core": vel[0], "v_out": vel[1], "relativity": relativity,
         "rho_1d": rho_1d, "x_lan": x_lan,
         "lam_window": LAM_WIN, "n_spec": N_SPEC,
         "T_gas_1d": T_GAS_1D, "epochs": list(epochs), "bands": phot.BANDS_PHOT,
         "distance_cm": phot.D_40MPC, "crossing": cross, "rows": rows}, indent=1))
    print(f"wrote {out_name}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="58CeII")
    ap.add_argument("--core", default="cool", choices=("cool", "fixed"))
    ap.add_argument("--n", type=float, default=1e6)
    ap.add_argument("--v", default="f40", choices=("f40", "kn"),
                    help="shell velocities: f40 = Paper I's 0.003-0.01c, "
                         "kn = kilonova 0.05-0.2c")
    ap.add_argument("--relativity", default=None, choices=(None, "worldline"))
    ap.add_argument("--mass", type=float, default=None,
                    help="ejecta mass in Msun; derives rho(1 d) from (M, v_max) "
                         "instead of using F40's tuned RHO_1D")
    ap.add_argument("--xlan", type=float, default=X_LAN)
    ap.add_argument("--label", default=None, help="suffix for the output file")
    ap.add_argument("--epochs", default=None,
                    help="comma-separated epochs in days (default: F40's)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    vel = V_F40 if a.v == "f40" else V_KN
    rho = RHO_1D if a.mass is None else rho_1d_from_mass(a.mass, vel[1])
    tag = (f"observables_{a.ion}_{a.core}"
           + ("" if a.v == "f40" else f"_{a.v}")
           + ("" if a.relativity is None else f"_{a.relativity}")
           + ("" if a.label is None else f"_{a.label}"))
    eps = EPOCHS if a.epochs is None else tuple(float(x) for x in a.epochs.split(","))
    main(a.ion, int(a.n), a.core, a.out or tag + ".json", vel, a.relativity,
         rho, a.xlan, eps)
