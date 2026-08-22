"""Paper II Phase 1: La II under five treatments, one code, two experiments.

PART A -- calibration on Paper I's window atom. Opacity only from the lines
in 3850-3950 A (the atom Paper I's SEDONA runs and analytic legs carried),
flat launch over SEDONA's transport window 3770-3997 A, and for the thermal
legs re-emission confined to that window -- which is what SEDONA's radiative
equilibrium does, because its emissivity CDF lives on the transport grid.
Must reproduce: sobolev_absorb 0.351, expansion_absorb 0.485 (Paper I's
analytic legs), and RE N=1 resolved 0.836 / expansion 0.900 (+7.7%).
The branch leg on this atom is the first physics result: the window lines'
upper levels decay through 71 channels spanning 3105-6333 A, only 16% of
them back into the window.

PART B -- the full atom. All 949 opacity lines, 1142-17,697 A; launch
photon-number-weighted Planck at T_core = 6000 K over that range, so UV pumps
carry their physical weight; thermal re-emission over the whole LTE
emissivity (mostly infrared at 3000 K). The three legs, plus the two
attenuation controls. This is the experiment Paper I could not do.

Band quantity: escaped / launched per unit frequency in 3800-3955 A -- with a
flat launch the transmitted fraction, with a Planck launch the emergent-to-
incident ratio. Also 3955-3995 A (immediately redward; line-free only in
the window atom) and the destination of the band's own launched photons.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C
from forest_mc import ForestAtom, band_ratio, run_mc, spectrum

R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE = 8.64e12, 2.592e13, 86400.0, 3000.0, 6000.0
BAND = (3800.0, 3955.0); RED = (3955.0, 3978.0)
BLUEWING = (3805.0, 3850.0)                 # the forest's 3000 km/s blueshifted absorption
WINDOW = (3770.0, 3997.0); FOREST = (3850.0, 3950.0)
# SEDONA RE N=1 (seed 1), raw L/L_cont on the RE-off red-margin scale, for comparison:
SEDONA_RE = {"bb": dict(bluewing=0.418, forest=0.990, band=0.835, red=1.108),
             "exp": dict(bluewing=0.632, forest=0.998, band=0.900, red=1.083),
             "off_bb": dict(bluewing=0.166, forest=0.360, band=0.346, red=1.002),
             "off_exp": dict(bluewing=0.446, forest=0.470, band=0.498, red=1.000)}
LEV, TR = ROOT / "data/57LaII_levels_calib.txt", ROOT / "data/57LaII_transitions_calib.txt"


def nu_of(lam_lo, lam_hi):
    return C / (lam_hi * 1e-8), C / (lam_lo * 1e-8)


def destination(res, band):
    lo, hi = nu_of(*band)
    inb = (res["nu_launch"] >= lo) & (res["nu_launch"] < hi)
    out = res["nu_out_all"][inb]; esc = np.isfinite(out)
    if not esc.any():
        return dict(escaped=0.0)
    lam = C / out[esc] * 1e8
    return dict(escaped=float(esc.mean()), in_band=float(np.mean((lam >= band[0]) & (lam < band[1]))),
                redward=float(np.mean(lam >= band[1])), blueward=float(np.mean(lam < band[0])),
                beyond_1um=float(np.mean(lam > 10000)))


def run_part(name, atom, legs, nu_min, nu_max, n, seed, t_core, edges):
    print(f"\n=== {name}: {atom.n_opacity} opacity lines, launch {C/nu_max*1e8:.0f}-{C/nu_min*1e8:.0f} A"
          f"{' (Planck %d K)' % t_core if t_core else ' (flat)'}, {n:.0e} packets/leg ===", flush=True)
    out, spectra = {}, {}
    for tag, mode, win in legs:
        t0 = time.time()
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, nu_min, nu_max, n, mode, seed=seed,
                     emit_window=nu_of(*win) if win else None, t_core=t_core)
        b, be = band_ratio(res, *nu_of(*BAND)); r, re_ = band_ratio(res, *nu_of(*RED))
        bw, _ = band_ratio(res, *nu_of(*BLUEWING)); fo, _ = band_ratio(res, *nu_of(*FOREST))
        b_e, _ = band_ratio(res, *nu_of(*BAND), weight="energy")
        acc = res["accounting"]
        dest = destination(res, BAND)
        out[tag] = dict(mode=mode, emit_window=win, band=b, band_err=be, red=r, red_err=re_,
                        bluewing=bw, forest=fo, band_energy_weighted=b_e,
                        accounting=dict(identity_residual=acc["identity_residual"],
                                        E_dep_cm_over_E_inj=acc["E_dep_cm"] / acc["E_inj"],
                                        W_over_E_interacting=acc["W"] / max(acc["E_interacting"], 1e-300),
                                        E_esc_over_E_inj=acc["E_esc"] / acc["E_inj"],
                                        E_core_over_E_inj=acc["E_core"] / acc["E_inj"]),
                        destination_of_band_launches=dest, n_escaped=res["n_escaped"], n_core=res["n_core"],
                        n_absorbed=res["n_absorbed"], n_interactions=res["n_interactions"],
                        steps=res["steps"], seconds=time.time() - t0)
        sp, _ = spectrum(res, edges); spectra[tag] = sp
        d = dest
        print(f"  {tag:28s} bluewing {bw:.3f} forest {fo:.3f} BAND {b:.4f}+-{be:.4f} (E-wt {b_e:.4f}) red {r:.3f} | "
              f"E: dep_cm {acc['E_dep_cm']/acc['E_inj']:+.4f} W {acc['W']/max(acc['E_interacting'],1e-300):+.2e} id {acc['identity_residual']:.1e} | band launches: esc {d.get('escaped',0):.3f} "
              f"in-band {d.get('in_band',0):.2f} red {d.get('redward',0):.2f} blue {d.get('blueward',0):.2f} >1um {d.get('beyond_1um',0):.2f}"
              f"  [{res['steps']} steps, {time.time()-t0:.0f}s]", flush=True)
    return out, spectra


def main(n, seed, tag):
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz"); n_ion = float(d["n_ion"])
    results = dict(n_packets=n, seed=seed, geometry=dict(r_core=R_CORE, r_out=R_OUT, t_exp=T_EXP, T=T_SHELL, n_ion=n_ion, T_core=T_CORE))

    # ---- Part A: window atom, window launch
    atom_w = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3, opacity_window=nu_of(*FOREST))
    lo, hi = nu_of(*WINDOW)
    edges_w = np.geomspace(lo, hi, 400)
    legs_a = [("sobolev_absorb", "sobolev_absorb", None), ("expansion_absorb", "expansion_absorb", None),
              ("sobolev_thermal_window", "sobolev_thermal", FOREST), ("expansion_thermal_window", "expansion_thermal", FOREST),
              ("sobolev_thermal_fullemis", "sobolev_thermal", None),
              ("sobolev_branch", "sobolev_branch", None)]
    A, spA = run_part("PART A  window atom (Paper I's)", atom_w, legs_a, lo, hi, n, seed, None, edges_w)
    print("  SEDONA for comparison (raw L/L_cont, RE-off margin scale):")
    for k, v in SEDONA_RE.items():
        print(f"    {k:10s} bluewing {v['bluewing']:.3f} forest {v['forest']:.3f} BAND {v['band']:.3f} red {v['red']:.3f}")
    g = lambda k: A[k]["band"]
    results["part_A"] = dict(atom=dict(n_opacity=atom_w.n_opacity), legs=A, differentials={
        "exp_vs_sob_absorb  [Paper I analytic: +38.1%]": g("expansion_absorb") / g("sobolev_absorb") - 1,
        "exp_vs_sob_thermal_window  [SEDONA RE N=1: 0.900/0.835 = +7.8%]": g("expansion_thermal_window") / g("sobolev_thermal_window") - 1,
        "branch_vs_sob_thermal_window": g("sobolev_branch") / g("sobolev_thermal_window") - 1,
        "branch_vs_sob_absorb": g("sobolev_branch") / g("sobolev_absorb") - 1,
        "exp_thermal_window_vs_branch": g("expansion_thermal_window") / g("sobolev_branch") - 1,
    })
    print("  differentials:"); [print(f"    {k:52s} {100*v:+7.2f}%") for k, v in results["part_A"]["differentials"].items()]

    # ---- Part B: full atom, Planck launch
    atom_f = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    lo_f, hi_f = atom_f.op_nu.min() * 0.995, atom_f.op_nu.max() * 1.005
    edges_f = np.geomspace(lo_f, hi_f, 3000)
    legs_b = [("sobolev_absorb", "sobolev_absorb", None), ("expansion_absorb", "expansion_absorb", None),
              ("sobolev_thermal", "sobolev_thermal", None), ("expansion_thermal", "expansion_thermal", None),
              ("sobolev_branch", "sobolev_branch", None)]
    B, spB = run_part("PART B  full atom, Planck launch", atom_f, legs_b, lo_f, hi_f, n, seed, T_CORE, edges_f)
    g = lambda k: B[k]["band"]
    results["part_B"] = dict(atom=dict(n_opacity=atom_f.n_opacity, launch_A=[C/hi_f*1e8, C/lo_f*1e8]), legs=B, differentials={
        "exp_vs_sob_absorb": g("expansion_absorb") / g("sobolev_absorb") - 1,
        "exp_vs_sob_thermal": g("expansion_thermal") / g("sobolev_thermal") - 1,
        "branch_vs_sob_thermal  [what fluorescence does]": g("sobolev_branch") / g("sobolev_thermal") - 1,
        "branch_vs_sob_absorb  [fill-in by fluorescence]": g("sobolev_branch") / g("sobolev_absorb") - 1,
        "exp_thermal_vs_branch  [closure code vs physics]": g("expansion_thermal") / g("sobolev_branch") - 1,
    })
    print("  differentials:"); [print(f"    {k:52s} {100*v:+7.2f}%") for k, v in results["part_B"]["differentials"].items()]

    (HERE / f"forest_results{tag}.json").write_text(json.dumps(results, indent=1))
    np.savez(HERE / f"forest_spectra{tag}.npz", edges_A=edges_w, edges_B=edges_f,
             **{f"A_{k}": v for k, v in spA.items()}, **{f"B_{k}": v for k, v in spB.items()})
    print(f"\nwrote forest_results{tag}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e5); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--tag", default="")
    a = ap.parse_args(); main(int(a.n), a.seed, a.tag)
