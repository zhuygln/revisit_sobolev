"""E2: thermal-width convergence -- SEDONA's resolved+thermal vs the MC's
Sobolev+thermal, as the Doppler width shrinks.

Inside the MC the Sobolev legs are delta resonances and the closure depends
on bin width, not v_D, so Delta_closure(v_D) is constant by construction. What
can be tested is the SEDONA side: Paper I established resolved -> Sobolev as
v_D -> 0, so the 3.7% by which the MC's window-confined Sobolev+thermal leg
sits above SEDONA's RE N=1 resolved run at 100 km/s should shrink at 10 and
1 km/s, while the expansion side (MC 0.905 vs SEDONA 0.900) should not move.
Prediction stated before the 10 km/s runs: SEDONA bb band 0.835 -> ~0.866,
Delta_SEDONA +7.8% -> ~+4.5%.

Also the closure's own bin-width systematic under re-emission: the MC
expansion legs at dnu/nu = 4.17e-6 / 4.17e-5 / 4.17e-4 (1.25 / 12.5 / 125 km/s
bins; SEDONA's grids at 10 / 100 km/s widths are the first two).
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C
from sobolev.formal_transfer import planck_bnu
from forest_mc import ForestAtom, band_ratio, run_mc
from run_forest import (BAND, BLUEWING, FOREST, RED, WINDOW, LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, nu_of)

LF = ROOT / "experiments/laII_forest"
T_CORE = 6000.0


def sed_bands(path, scale):
    s = np.loadtxt(path, comments="#"); nu, l = s[:, 0], s[:, 1]; g = l > 0; nu, l = nu[g], l[g]
    lam = C / nu * 1e8; r = l / (4 * np.pi**2 * R_CORE**2 * planck_bnu(nu, T_CORE)) / scale
    out = {}
    for name, (lo, hi) in (("bluewing", BLUEWING), ("forest", FOREST), ("band", BAND), ("red", RED)):
        m = (lam > lo) & (lam < hi); out[name] = float(r[m].mean())
    return out


def main(n=2_000_000, seeds=(1, 2, 3)):
    # the RE-off red-margin scale, so that RE-off margins read 1
    s0 = np.loadtxt(LF / "run_bb/spectrum_1.dat", comments="#"); nu0, l0 = s0[:, 0], s0[:, 1]; g = l0 > 0
    lam0 = C / nu0[g] * 1e8; r0 = l0[g] / (4 * np.pi**2 * R_CORE**2 * planck_bnu(nu0[g], T_CORE))
    scale = r0[(lam0 > 3952) & (lam0 < 3970)].mean()

    sedona = {}
    for vd, tag in ((100, ""), (10, "_vd10"), (1, "_vd1")):
        for mode in ("bb", "exp"):
            vals = []
            for seed in (1, 2, 3):
                p = LF / f"re_prod_N1_{mode}{tag}_s{seed}" / "spectrum_1.dat"
                if p.exists():
                    vals.append(sed_bands(p, scale))
            if vals:
                sedona[f"{mode}_vd{vd}"] = {k: float(np.mean([v[k] for v in vals])) for k in vals[0]} | {"n_seeds": len(vals)}

    d = np.load(LF / "forest_lines.npz"); n_ion = float(d["n_ion"])
    atom = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3, opacity_window=nu_of(*FOREST))
    lo, hi = nu_of(*WINDOW)
    mc = {}
    for label, mode, dnu in (("sob_thermal", "sobolev_thermal", 4.17e-5),
                             ("exp_thermal_bin1.25", "expansion_thermal", 4.17e-6),
                             ("exp_thermal_bin12.5", "expansion_thermal", 4.17e-5),
                             ("exp_thermal_bin125", "expansion_thermal", 4.17e-4),
                             ("exp_absorb_bin1.25", "expansion_absorb", 4.17e-6),
                             ("exp_absorb_bin12.5", "expansion_absorb", 4.17e-5),
                             ("exp_absorb_bin125", "expansion_absorb", 4.17e-4)):
        rows = []
        for seed in seeds:
            res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=seed,
                         emit_window=nu_of(*FOREST), dnu_over_nu=dnu)
            rows.append({k: band_ratio(res, *nu_of(*b), weight="energy")[0]
                         for k, b in (("bluewing", BLUEWING), ("forest", FOREST), ("band", BAND), ("red", RED))})
        mc[label] = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]} | {"band_std": float(np.std([r["band"] for r in rows], ddof=1))}

    print("SEDONA RE N=1 (energy ratio on the RE-off margin scale) vs MC (energy-weighted), Paper I's atom")
    print(f"{'leg':22s} {'bluewing':>9s} {'forest':>8s} {'band':>8s} {'red':>7s}")
    for k, v in sedona.items():
        print(f"SEDONA {k:15s} {v['bluewing']:9.3f} {v['forest']:8.3f} {v['band']:8.4f} {v['red']:7.3f}   ({v['n_seeds']} seeds)")
    for k, v in mc.items():
        print(f"MC     {k:15s} {v['bluewing']:9.3f} {v['forest']:8.3f} {v['band']:8.4f} {v['red']:7.3f}")
    print("\nDelta (exp vs resolved/Sobolev) in 3800-3955 A:")
    for vd in (100, 10, 1):
        if f"bb_vd{vd}" in sedona and f"exp_vd{vd}" in sedona:
            print(f"  SEDONA v_D={vd:3d} km/s   {100*(sedona[f'exp_vd{vd}']['band']/sedona[f'bb_vd{vd}']['band']-1):+6.2f}%")
    print(f"  MC (v_D-independent) {100*(mc['exp_thermal_bin12.5']['band']/mc['sob_thermal']['band']-1):+6.2f}%   "
          f"[bin 1.25: {100*(mc['exp_thermal_bin1.25']['band']/mc['sob_thermal']['band']-1):+.2f}%, bin 125: {100*(mc['exp_thermal_bin125']['band']/mc['sob_thermal']['band']-1):+.2f}%]")
    for vd in (100, 10, 1):
        if f"bb_vd{vd}" in sedona:
            print(f"  MC Sob+thermal vs SEDONA resolved+thermal at {vd} km/s: {100*(mc['sob_thermal']['band']/sedona[f'bb_vd{vd}']['band']-1):+.2f}% (band)")
    (HERE / "e2_thermal_width.json").write_text(json.dumps(dict(sedona=sedona, mc=mc), indent=1))
    print("wrote e2_thermal_width.json")


if __name__ == "__main__":
    main()
