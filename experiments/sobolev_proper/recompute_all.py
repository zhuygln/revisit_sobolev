"""Recompute every Sobolev-vs-expansion number through the shared band-ratio
helper, and report what moved.

The previous values compared a correctly-normalized analytic leg against a
SEDONA band flux normalized by RAW LUMINOSITY, leaving the Planck slope in
the answer. That produced a spurious ~5-7% "v_D-independent Sobolev floor".
Delta_expansion, being a same-code differential, should barely move -- that
invariance is the check that this helper changed only what it should.

No new transport runs: every spectrum is already on disk.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.sobolev_leg import sobolev_attenuation, expansion_damp
from sobolev.spectra import band_average, band_ratio

FOREST = ROOT / "experiments/laII_forest"
T_EXP, R_CORE, R_OUT, T_CORE = 86400.0, 8.64e12, 2.592e13, 6000.0
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)

d = np.load(FOREST / "forest_lines.npz")
N_ION = float(d["n_ion"])
LINES = [(C / (l * 1e-8), f, p) for l, f, p in zip(d["lam"], d["f_lu"], d["pop"])]
nu = np.geomspace(7.50e14, 7.95e14, 4000)
lam_an = C / nu * 1e8


def analytic(scale, damp=None):
    att = sobolev_attenuation(nu, LINES, R_CORE, R_OUT, T_EXP, N_ION * scale,
                              damp=damp)
    return band_average(lam_an, att, BAND)


def sedona(run):
    return band_ratio(FOREST / run / "spectrum_1.dat", BAND, MARGIN,
                      R_CORE, T_CORE)


rows = []
old = json.loads((FOREST / "sweep_results.json").read_text())

print("=== tau_max x v_D grid ===")
print("tau_max  v_D | D_Sob new    old |  D_exp new    old")
for tau in (0.5, 5.0, 50.0):
    f_sob = analytic(tau / 5.0)
    f_exp_ana = analytic(tau / 5.0, damp=expansion_damp)
    for vd in (10, 30, 100, 300):
        fb, fe = sedona(f"sweep_tau{tau:g}_vd{vd:g}_bb"), sedona(f"sweep_tau{tau:g}_vd{vd:g}_exp")
        o = [x for x in old if x["tau_max"] == tau and x["v_d_kms"] == vd][0]
        d_sob, d_exp = (f_sob - fb) / fb, (fe - fb) / fb
        d_sob_old = (f_sob - o["bb"]) / o["bb"]
        rows.append(dict(tau_max=tau, v_d=vd, f_res=fb, f_exp=fe, f_sob=f_sob,
                         f_exp_ana=f_exp_ana, d_sob=d_sob, d_exp=d_exp))
        print(f"{tau:6g} {vd:4d} | {d_sob:+8.1%} {d_sob_old:+8.1%} |"
              f" {d_exp:+8.1%} {o['delta_sob']:+8.1%}")

print("\n=== thermal-width frontier (tau_max = 5) ===")
print(" v_D | F_resolved  D_Sobolev   D_expansion")
frontier = []
f_sob5 = analytic(1.0)
for vd, rb, re in ((1, "tsweep_T3000_vd1_bb", "tsweep_T3000_vd1_exp"),
                   (3, "tsweep_T3000_vd3_bb", "tsweep_T3000_vd3_exp"),
                   (10, "sweep_tau5_vd10_bb", "sweep_tau5_vd10_exp"),
                   (30, "sweep_tau5_vd30_bb", "sweep_tau5_vd30_exp"),
                   (100, "sweep_tau5_vd100_bb", "sweep_tau5_vd100_exp"),
                   (300, "sweep_tau5_vd300_bb", "sweep_tau5_vd300_exp")):
    fb, fe = sedona(rb), sedona(re)
    frontier.append(dict(v_d=vd, f_res=fb, f_exp=fe, f_sob=f_sob5,
                         d_sob=(f_sob5 - fb) / fb, d_exp=(fe - fb) / fb))
    print(f"{vd:4d} | {fb:10.4f} {(f_sob5-fb)/fb:+10.1%} {(fe-fb)/fb:+12.1%}")

# headline forest condition
fb = sedona("run_bb"); fe = sedona("run_exp")
print(f"\n=== La II forest headline (v_D = 100 km/s, tau_max = 5) ===")
print(f"  SEDONA resolved  {fb:.4f}")
print(f"  Sobolev proper   {f_sob5:.4f}   Delta = {(f_sob5-fb)/fb:+.1%}")
print(f"  expansion        {fe:.4f}   Delta = {(fe-fb)/fb:+.1%}")

(Path(__file__).parent / "separation_results.json").write_text(json.dumps(
    dict(grid=rows, frontier=frontier,
         headline=dict(res=fb, sob=f_sob5, exp=fe)), indent=1))
print("\nwrote separation_results.json")
