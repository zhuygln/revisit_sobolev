"""Breadth sweep: windows x epochs x ion mixes.

Tests whether F6 (strength floor + width wing), F7 (non-monotonic density)
and F9 (Sobolev-proper vs expansion separation) are general or artifacts of
the single 3850-3950 A / day-1 / La II condition mapped so far.

Axes
----
windows : 4, auto-selected by optical-depth richness over the combined
          La II + Ce II line list, spanning optical -> NIR.
epochs  : t = 0.5, 1, 3 days at FIXED ejecta mass, so rho ~ t^-3 and the
          radii r = v t scale with epoch. tau_S ~ n_l t ~ t^-2, giving a ~36x
          strength range across the epoch axis. The realized tau_max is
          recorded per condition -- it, not t, is the map's true abscissa.
ions    : La II | La II + Ce II | La II + Ce II + Ce III.

Population control (research_requirements section 17) is preserved by giving
each SPECIES its own element group with ionization disabled. Ce III cannot
share Z=58 with Ce II under that scheme, so it is carried in the Z=59 slot
with Ce III atomic data and the Ce mass number: with ionization disabled Z is
only a label, and this keeps the II/III ratio fixed by hand (as a composition
choice) instead of handing it to Saha.

Per condition: 2 SEDONA runs (resolved, expansion) + both analytic legs.
The deterministic solver is skipped -- Delta is a same-code differential and
the solver is the slow leg.
"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.atomic_data import load_gsi
from sobolev.constants import C, SIGMA_CLASSICAL
from sobolev.populations import boltzmann_fractions_from_levels, statistical_weight
from sobolev.sobolev_leg import expansion_damp, sobolev_attenuation

from sobolev.sedona import sedona_cmd, sedona_home, sedona_timeout

SEDONA_HOME = sedona_home()
M_P = 1.67262192e-24
CM1_TO_EV = 1.239841984e-4
SB = 5.670374e-5

T_SHELL = 3000.0
T_CORE = 6000.0
V_CORE, V_MAX = 1.0e8, 3.0e8
N_ZONES = 100
V_D = 1.0e7  # 100 km/s
EPOCHS_DAY = [0.5, 1.0, 3.0]
TAU_REF = 5.0  # tau_max at t=1 day in the reference window, La II only

SPECIES = {  # label -> (Z slot, mass number, levels file, transitions file)
    "LaII": (57, 139, "57LaII_levels_calib.txt", "57LaII_transitions_calib.txt"),
    "CeII": (58, 140, "58CeII_levels_calib.txt", "58CeII_transitions_calib.txt"),
    "CeIII": (59, 140, "58CeIII_levels_calib.txt", "58CeIII_transitions_calib.txt"),
}
MIXES = {
    "La": ["LaII"],
    "LaCe": ["LaII", "CeII"],
    "LaCeCe3": ["LaII", "CeII", "CeIII"],
}

print("loading atomic data ...", flush=True)
DATA = {}
for name, (z, a, lev_f, tr_f) in SPECIES.items():
    levels = load_gsi(ROOT / "data" / lev_f)
    lines = load_gsi(ROOT / "data" / tr_f)
    DATA[name] = dict(z=z, a=a, levels=levels, lines=lines)
    print(f"  {name}: {len(levels)} levels, {len(lines)} lines", flush=True)


def window_lines(name, lo, hi, t_exp):
    """Window-selected lines with Boltzmann populations and tau per unit rho."""
    d = DATA[name]
    lam = d["lines"]["WV_Transition"].to_numpy()
    win = d["lines"][(lam >= lo) & (lam < hi)].reset_index(drop=True)
    if len(win) == 0:
        return None
    frac = boltzmann_fractions_from_levels(d["levels"], T_SHELL)
    pop = frac[win["Lower"].to_numpy()]
    g_l = statistical_weight(win["J_Lower"].to_numpy())
    f_lu = 10 ** win["Log(gf)"].to_numpy() / g_l
    # tau per unit TOTAL mass density, at mass fraction x (applied by caller)
    tau_per_rho = (
        SIGMA_CLASSICAL * f_lu * pop * win["WV_Transition"].to_numpy() * 1e-8
        * t_exp / (d["a"] * M_P)
    )
    return dict(win=win, pop=pop, f_lu=f_lu, tau_per_rho=tau_per_rho, **d)


# ---------------- window selection ----------------
print("scanning windows ...", flush=True)
t_ref = 86400.0
cands = []
for lo in range(3500, 10000, 100):
    tot, strong = 0, 0.0
    for nm in ("LaII", "CeII"):
        w = window_lines(nm, lo, lo + 100, t_ref)
        if w is None:
            continue
        tot += len(w["win"])
        strong = max(strong, w["tau_per_rho"].max())
    if tot > 20:
        cands.append((lo, tot, strong))
# spread the picks across wavelength: best-by-strength within each band
BANDS = [(3500, 4500), (4500, 6000), (6000, 8000), (8000, 10000)]
WINDOWS = []
for blo, bhi in BANDS:
    inb = [c for c in cands if blo <= c[0] < bhi]
    if inb:
        best = max(inb, key=lambda c: c[2])
        WINDOWS.append((float(best[0]), float(best[0] + 100)))
print("windows:", WINDOWS, flush=True)

PARAM = """sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "{atom}"
grid_type    = "grid_1D_sphere"
model_file   = "{mod}"
hydro_module = "homologous"
transport_nu_grid  = {{{nu_lo:.5e}, {nu_hi:.5e}, 4.17e-5, 1}}
spectrum_nu_grid   = {{{nu_lo:.5e}, {nu_hi:.5e}, 1.0e-4, 1}}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1
texp             = {texp:.6e}
tstep_time_start = {texp:.6e}
core_n_emit      = 1e6
core_radius      = {rcore:.6e}
core_temperature = {tcore:.1f}
core_luminosity  = {lum:.6e}
opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = {bb}
opacity_line_expansion      = {exp}
opacity_epsilon             = 1
line_velocity_width         = {vd:.3e}
output_write_radiation = 0
"""


def write_atom(path, specs, lo, hi, t_exp):
    with h5py.File(path, "w") as f:
        for s in specs:
            n_lev = len(s["levels"])
            g = f.create_group(str(s["z"]))
            g.attrs["n_ions"] = np.int64(2)
            g.attrs["n_levels"] = np.int64(n_lev + 1)
            g.attrs["n_lines"] = np.int64(len(s["win"]))
            g.create_dataset("ion_chi", data=np.array([9.9999e4, 9.9999e4]))
            g.create_dataset("ion_ground", data=np.array([0, n_lev], dtype=np.int64))
            g.create_dataset(
                "level_E",
                data=np.concatenate(
                    [s["levels"]["Energy"].to_numpy() * CM1_TO_EV, [0.0]]
                ),
            )
            g.create_dataset(
                "level_g",
                data=np.concatenate(
                    [statistical_weight(s["levels"]["J"].to_numpy()).astype(np.int64), [1]]
                ),
            )
            g.create_dataset("level_i", data=np.array([0] * n_lev + [1], dtype=np.int64))
            g.create_dataset("line_A", data=s["win"]["A"].to_numpy())
            g.create_dataset("line_l", data=s["win"]["Lower"].to_numpy().astype(np.int64))
            g.create_dataset("line_u", data=s["win"]["Upper"].to_numpy().astype(np.int64))


def band_from_spectrum(path, lam_lo, lam_hi, red_lo, red_hi, r_core):
    s = np.loadtxt(path, comments="#")
    nu, lum = s[:, 0], s[:, 1]
    lam = C / nu * 1e8
    red = (lam > red_lo) & (lam < red_hi) & (lum > 0)
    if red.sum() < 3 or np.mean(lum[red]) <= 0:
        return None
    ratio = lum / np.mean(lum[red])
    m = (lam > lam_lo) & (lam < lam_hi)
    o = np.argsort(lam[m])
    return np.trapezoid(ratio[m][o], lam[m][o]) / (lam_hi - lam_lo)


results = []
t0_all = time.time()
for mix_name, members in MIXES.items():
    for (lo, hi) in WINDOWS:
        for t_day in EPOCHS_DAY:
            t_exp = t_day * 86400.0
            r_core, r_out = V_CORE * t_exp, V_MAX * t_exp
            lum = 4 * np.pi * r_core**2 * SB * T_CORE**4
            x_frac = 1.0 / len(members)

            specs = []
            for nm in members:
                w = window_lines(nm, lo, hi, t_exp)
                if w is not None:
                    specs.append(w)
            if not specs:
                continue

            # fixed ejecta mass: rho set so tau_max = TAU_REF for La II alone
            # in the reference window at 1 day, then scaled as t^-3.
            ref = window_lines("LaII", 3850.0, 3950.0, 86400.0)
            rho_1day = TAU_REF / (ref["tau_per_rho"].max() * 1.0)
            rho = rho_1day * (1.0 / t_day) ** 3
            taus_all = np.concatenate(
                [s["tau_per_rho"] * rho * x_frac for s in specs]
            )
            tau_max = float(taus_all.max())
            n_lines = int(sum(len(s["win"]) for s in specs))

            tag = f"{mix_name}_w{int(lo)}_t{t_day:g}"
            atom = HERE / f"atom_{mix_name}_w{int(lo)}_t{t_day:g}.hdf5"
            write_atom(atom, specs, lo, hi, t_exp)
            mod = HERE / f"model_{tag}.mod"
            v_edges = np.linspace(V_CORE, V_MAX, N_ZONES + 1)[1:]
            with open(mod, "w") as fh:
                fh.write("1D_sphere standard\n")
                fh.write(f"{N_ZONES}\t{r_core:.6e}\t{t_exp:.6e} {len(specs)} \n")
                fh.write(" ".join(f"{s['z']}.{s['a']}" for s in specs) + "\n")
                for v in v_edges:
                    fh.write(
                        f"{v*t_exp:.6e} {v:.6e} {rho:.6e} {T_SHELL:.6e} "
                        + " ".join(f"{x_frac:.4f}" for _ in specs) + "\n"
                    )

            lam_blue, lam_red_edge = lo * 0.985, hi * 1.010
            nu_lo, nu_hi = C / (lam_red_edge * 1e-8), C / (lam_blue * 1e-8)
            band_lo, band_hi = lo * 0.9885, hi + 1.0
            red_lo, red_hi = hi * 1.0015, hi * 1.008

            fluxes = {}
            for mode, bb, ex in [("bb", 1, 0), ("exp", 0, 1)]:
                run = HERE / f"run_{tag}_{mode}"
                run.mkdir(exist_ok=True)
                (run / "param.lua").write_text(
                    PARAM.format(
                        atom=f"../{atom.name}", mod=f"../{mod.name}",
                        nu_lo=nu_lo, nu_hi=nu_hi, texp=t_exp, rcore=r_core,
                        tcore=T_CORE, lum=lum, bb=bb, exp=ex, vd=V_D,
                    )
                )
                r = subprocess.run(
                    sedona_cmd(), cwd=run, capture_output=True, text=True,
                    env={**os.environ, "SEDONA_HOME": SEDONA_HOME}, timeout=sedona_timeout(2000),
                )
                if r.returncode != 0:
                    fluxes[mode] = None
                    print(f"FAIL {tag} {mode} rc={r.returncode}", flush=True)
                    continue
                fluxes[mode] = band_from_spectrum(
                    run / "spectrum_1.dat", band_lo, band_hi, red_lo, red_hi, r_core
                )

            # analytic legs on the same band
            nu_grid = np.geomspace(nu_lo, nu_hi, 3000)
            lam_grid = C / nu_grid * 1e8
            lines_an = []
            for s in specs:
                for lam_k, f_k, p_k in zip(
                    s["win"]["WV_Transition"], s["f_lu"], s["pop"]
                ):
                    lines_an.append(
                        (C / (lam_k * 1e-8), f_k, p_k * x_frac / (s["a"] * M_P))
                    )
            sob = sobolev_attenuation(nu_grid, lines_an, r_core, r_out, t_exp, rho)
            exa = sobolev_attenuation(
                nu_grid, lines_an, r_core, r_out, t_exp, rho, damp=expansion_damp
            )
            mm = (lam_grid > band_lo) & (lam_grid < band_hi)
            oo = np.argsort(lam_grid[mm])
            f_sob = float(
                np.trapezoid(sob[mm][oo], lam_grid[mm][oo]) / (band_hi - band_lo)
            )
            f_exa = float(
                np.trapezoid(exa[mm][oo], lam_grid[mm][oo]) / (band_hi - band_lo)
            )

            row = dict(
                mix=mix_name, window=lo, t_day=t_day, tau_max=tau_max,
                n_lines=n_lines, rho=rho, f_res=fluxes.get("bb"),
                f_exp=fluxes.get("exp"), f_sob=f_sob, f_exp_ana=f_exa,
            )
            if row["f_res"] and row["f_exp"]:
                row["d_exp"] = (row["f_exp"] - row["f_res"]) / row["f_res"]
                row["d_sob"] = (f_sob - row["f_res"]) / row["f_res"]
                print(
                    f"{tag:22s} lines={n_lines:5d} tau_max={tau_max:9.2f}  "
                    f"res={row['f_res']:.4f}  D_sob={row['d_sob']:+7.1%}  "
                    f"D_exp={row['d_exp']:+7.1%}", flush=True
                )
            results.append(row)
            (HERE / "breadth_results.json").write_text(json.dumps(results, indent=1))

print(f"done in {(time.time()-t0_all)/60:.1f} min; {len(results)} conditions",
      flush=True)
