"""E7: which cascades refill the optical band, and from where?

For escaped packets with lambda_out in 3800-3955 A (sobolev_branch, full
La II, Planck 6000 K launch): the pathway key is (first absorbing line ->
last emitting line), ranked by escaped energy; for each top key the pump
wavelength, the pumped upper level (index and configuration), the exit
line's wavelength and tau, the mean number of events and re-absorptions, and
the share of the band's escaped energy. Plus a launch-band breakdown of the
band's escaped energy -- the test of the "UV pumps" attribution: the 6000 K
Planck photon budget in 1142-2500 A is only 0.32% of the total, less than the
band's refill, so the pumps are expected to be 2500-4500 A.

A "cascade" here is one emission per absorption event (the exit lower level's
excitation returns to the pool); multi-event pathways are re-absorptions by
other lines elsewhere.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from sobolev.constants import C, H
from forest_mc import ForestAtom, run_mc
from run_forest import LEV, TR, R_CORE, R_OUT, T_EXP, T_SHELL, T_CORE, nu_of, BAND

LAUNCH_BANDS = [(1142, 2500), (2500, 3300), (3300, 3800), (3800, 3955), (3955, 4500), (4500, 6000), (6000, 17697)]


def main(n, seed, mode):
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz"); n_ion = float(d["n_ion"])
    atom = ForestAtom.from_gsi(LEV, TR, T_SHELL, n_ion, T_EXP, tau_min=1e-3)
    lo, hi = atom.op_nu.min() * 0.995, atom.op_nu.max() * 1.005
    res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=seed, t_core=T_CORE)
    blo, bhi = nu_of(*BAND)
    esc = (res["fate"] == 1) & (res["nu_out_all"] >= blo) & (res["nu_out_all"] < bhi)
    E_out = H * res["nu_out_all"][esc]; E_tot = E_out.sum()
    lam_in = C / res["nu_launch"][esc] * 1e8
    print(f"{mode}, {n:.0e} packets: {esc.sum()} escaped in 3800-3955 A; band escaped energy / launched-in-band energy = "
          f"{E_tot / (H * res['nu_launch'][(res['nu_launch'] >= blo) & (res['nu_launch'] < bhi)]).sum():.3f}")
    out = {"mode": mode, "n": n, "seed": seed, "n_escaped_in_band": int(esc.sum())}

    # --- launch-band breakdown of the band's escaped energy
    print("\nwhere the band's escaped energy was LAUNCHED:")
    lb = {}
    for lo_, hi_ in LAUNCH_BANDS:
        m = (lam_in >= lo_) & (lam_in < hi_); lb[f"{lo_}-{hi_}"] = float(E_out[m].sum() / E_tot)
        print(f"  {lo_:5d}-{hi_:5d} A  {100*E_out[m].sum()/E_tot:6.2f}%")
    out["launch_band_share"] = lb
    direct = (res["n_events"][esc] == 0)
    print(f"  (transmitted without any interaction: {100*E_out[direct].sum()/E_tot:.2f}%)")
    out["direct_share"] = float(E_out[direct].sum() / E_tot)

    # --- pathways: (first_line -> last_line)
    fl, ll = res["first_line"][esc], res["last_line"][esc]
    inter = fl >= 0
    keys = fl[inter].astype(np.int64) * 100000 + ll[inter]
    uniq, inv = np.unique(keys, return_inverse=True)
    e_key = np.bincount(inv, weights=E_out[inter]); order = np.argsort(-e_key)
    lev = atom.levels; conf = lev["Configuration"].to_numpy() if "Configuration" in lev else None
    lam_all = C / atom.nu0_all * 1e8
    print(f"\ntop pathways (first absorbing line -> last emitting line), by share of the band's escaped energy "
          f"(interacting packets carry {100*E_out[inter].sum()/E_tot:.1f}% of it):")
    print(f"  {'share':>6s}  {'pump A':>8s} {'tau_pump':>8s}  {'upper level':24s} {'exit A':>8s} {'tau_exit':>8s}  {'events':>6s} {'reabs':>5s}")
    rows = []
    for j in order[:25]:
        k = uniq[j]; f_, l_ = int(k // 100000), int(k % 100000); sel = inv == j
        u = int(atom.upper_all[f_])
        rows.append(dict(first_line=f_, last_line=l_, share=float(e_key[j] / E_tot), pump_A=float(lam_all[f_]),
                         tau_pump=float(atom.tau_all[f_]), upper=u, upper_conf=(str(conf[u]) if conf is not None else ""),
                         upper_E_cm=float(atom.level_energy_cm[u]), exit_A=float(lam_all[l_]), tau_exit=float(atom.tau_all[l_]),
                         mean_events=float(res["n_events"][esc][inter][sel].mean()), mean_reabs=float(res["n_reabs"][esc][inter][sel].mean())))
        r = rows[-1]
        print(f"  {100*r['share']:6.2f}% {r['pump_A']:8.1f} {r['tau_pump']:8.2f}  {r['upper']:4d} {r['upper_conf'][:18]:18s} {r['exit_A']:8.1f} {r['tau_exit']:8.2f}  {r['mean_events']:6.2f} {r['mean_reabs']:5.2f}")
    out["pathways"] = rows
    # how concentrated? cumulative share of the top N
    cum = np.cumsum(e_key[order]) / E_tot
    for N in (5, 10, 25, 50, 100):
        if N <= cum.size:
            print(f"  top {N:3d} pathways carry {100*cum[N-1]:.1f}% of the band's escaped energy ({uniq.size} pathways in all)")
    out["cumulative_share"] = {str(N): float(cum[N - 1]) for N in (5, 10, 25, 50, 100) if N <= cum.size}
    # pump wavelength distribution among interacting band photons (first absorbing line)
    lam_pump = lam_all[fl[inter]]
    print("\npump-line wavelength of the band's interacting escaped energy:")
    for lo_, hi_ in LAUNCH_BANDS:
        m = (lam_pump >= lo_) & (lam_pump < hi_)
        print(f"  {lo_:5d}-{hi_:5d} A  {100*E_out[inter][m].sum()/E_out[inter].sum():6.2f}%")
    (HERE / f"e7_pathways_{mode}.json").write_text(json.dumps(out, indent=1))
    print(f"wrote e7_pathways_{mode}.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=2e6); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mode", default="sobolev_branch")
    a = ap.parse_args(); main(int(a.n), a.seed, a.mode)
