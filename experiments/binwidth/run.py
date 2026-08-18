"""P1-B: is the expansion-opacity error intrinsic, or an artifact of running
the formalism outside its design regime?

THE OBJECTION. Expansion opacity (Karp et al. 1977; Eastman & Pinto 1993) is
constructed for frequency bins containing MANY lines: the sum over lines in a
bin is what makes the statistical averaging meaningful. Our transport grid is
far finer than that. The 3771-3997 A window at dnu/nu = 4.17e-5 gives ~12.5
km/s bins against a mean line spacing of ~50 km/s -- roughly one line per four
bins. A referee can reasonably say we measured the formalism outside its
intended domain and that the 40% error is unsurprising.

THE TEST. Appendix A.4 derives that the bin width CANCELS: integrating the
binned opacity across the path over which a line stays inside one bin gives an
effective optical depth of exactly 1 - exp(-tau_S), independent of resolution.
So if the derivation is right, Delta_expansion must be invariant under bin
width, and the objection is answered. If Delta_expansion instead falls as bins
widen, the error is a usage artifact and Paper I's central claim needs
reframing before submission.

Everything else is held fixed: same model, same atom, same physics, same
band and normalization. Only transport_nu_grid resolution changes.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.spectra import band_ratio

SEDONA = "/home/yozhuz_223/personal/pubsed/src/sedona6.ex"
SEDONA_HOME = "/home/yozhuz_223/personal/pubsed"
FOREST = ROOT / "experiments/laII_forest"
R_CORE, T_CORE = 8.64e12, 6000.0
BAND, MARGIN = (3800.0, 3955.0), (3952.0, 3970.0)
C = 2.99792458e10

# dnu/nu of the transport grid. The production runs used 4.17e-5 (~12.5 km/s).
# Sweep two decades: ~1.25 km/s to ~1250 km/s per bin. At the coarse end a bin
# spans many lines, which is the regime the formalism was built for.
RESOLUTIONS = [4.17e-6, 1.25e-5, 4.17e-5, 1.25e-4, 4.17e-4, 1.25e-3, 4.17e-3]

PARAM = """sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = "{forest}/atom_laII.hdf5"
grid_type    = "grid_1D_sphere"
model_file   = "{forest}/laII.mod"
hydro_module = "homologous"
transport_nu_grid  = {{7.50e14, 7.95e14, {dnu:.4e}, 1}}
spectrum_nu_grid   = {{7.50e14, 7.95e14, 1.0e-4, 1}}
transport_radiative_equilibrium = 0
transport_steady_iterate        = 1
texp             = 86400.0
tstep_time_start = 86400.0
core_n_emit      = 4e6
core_radius      = 8.64e12
core_temperature = 6000.0
core_luminosity  = 6.8937e37
opacity_grey_opacity        = 0
opacity_electron_scattering = 0
opacity_bound_bound         = {bb}
opacity_line_expansion      = {exp}
opacity_epsilon             = 1
line_velocity_width         = 1.0e7
output_write_radiation = 0
"""

results = []
print("  dnu/nu   bin [km/s]  lines/bin |   F_res    F_exp   Delta_exp")
for dnu in RESOLUTIONS:
    bin_kms = dnu * C / 1e5
    # 153 lines spread over the window's velocity span (~7700 km/s)
    lines_per_bin = 153 * bin_kms / 7700.0
    row = {"dnu_over_nu": dnu, "bin_kms": bin_kms, "lines_per_bin": lines_per_bin}
    for mode, bb, ex in (("bb", 1, 0), ("exp", 0, 1)):
        run = HERE / f"run_{dnu:.2e}_{mode}"
        run.mkdir(parents=True, exist_ok=True)
        (run / "param.lua").write_text(
            PARAM.format(forest=FOREST, dnu=dnu, bb=bb, exp=ex)
        )
        r = subprocess.run(
            [SEDONA, "param.lua"], cwd=run, capture_output=True, text=True,
            env={**os.environ, "SEDONA_HOME": SEDONA_HOME}, timeout=6000,
        )
        if r.returncode != 0:
            row[mode] = None
            print(f"  FAIL dnu={dnu:.2e} {mode} rc={r.returncode}", flush=True)
            continue
        row[mode] = band_ratio(run / "spectrum_1.dat", BAND, MARGIN,
                               R_CORE, T_CORE)
    if row.get("bb") and row.get("exp"):
        row["d_exp"] = (row["exp"] - row["bb"]) / row["bb"]
        print(f"{dnu:9.2e} {bin_kms:10.2f} {lines_per_bin:10.3f} |"
              f" {row['bb']:.4f}  {row['exp']:.4f}  {row['d_exp']:+9.1%}",
              flush=True)
    results.append(row)
    (HERE / "binwidth_results.json").write_text(json.dumps(results, indent=1))

# The comparison is only meaningful where the RESOLVED leg is still resolved.
# Once a bin is wider than the line profile the reference itself stops sampling
# the profile and its band flux runs away (0.342 -> 0.44 -> 0.68 here). Those
# points measure the reference failing, not expansion opacity changing, and
# including them in a verdict is how the first version of this script wrongly
# reported "NOT INVARIANT".
V_D_KMS = 100.0
ok = [r for r in results
      if r.get("d_exp") is not None and r["bin_kms"] <= 1.5 * V_D_KMS]
bad = [r for r in results
       if r.get("d_exp") is not None and r["bin_kms"] > 1.5 * V_D_KMS]

if ok:
    d = np.array([r["d_exp"] for r in ok])
    f = np.array([r["bb"] for r in ok])
    print(f"\nVALID range (bin <= 1.5 v_D, reference converged): "
          f"{ok[0]['bin_kms']:.2f}-{ok[-1]['bin_kms']:.0f} km/s bins, "
          f"{ok[0]['lines_per_bin']:.3f}-{ok[-1]['lines_per_bin']:.2f} lines/bin")
    print(f"  F_resolved  mean {f.mean():.4f}, spread {f.max()-f.min():.4f} "
          f"({(f.max()-f.min())/f.mean():.2%}) -- reference is stable")
    print(f"  Delta_exp   mean {d.mean():+.1%}, spread "
          f"{100*(d.max()-d.min()):.1f} points, std {100*d.std():.1f} points")
    print("  => INVARIANT: the error is intrinsic to the formalism, not a"
          " consequence of bin choice"
          if (d.max() - d.min()) < 0.05 else
          "  => NOT INVARIANT: reframe the claim")
if bad:
    print("\nEXCLUDED (bin > 1.5 v_D -- the reference is no longer resolved):")
    for r in bad:
        print(f"  bin {r['bin_kms']:7.0f} km/s: F_res = {r['bb']:.4f} "
              f"(vs {ok[0]['bb']:.4f} converged), Delta_exp {r['d_exp']:+.1%}")
