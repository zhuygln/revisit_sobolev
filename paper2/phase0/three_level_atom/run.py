"""P2-0B: does the branching Monte Carlo reproduce known branching ratios?

Paper II's question is whether Paper I's expansion-opacity bias survives
radiative redistribution. Answering it needs a code in which a photon absorbed
in one line can leave in another. The P2-0A audit established that the public
SEDONA cannot do this at all, so the instrument had to be built
(`branching_mc.py`); this script is its acceptance test.

Three measurements, each against a number known in advance:

  1. the fluorescent yield equals A32 / (A31 + A32);
  2. the interaction probability equals 1 - exp(-tau_S);
  3. the pure-absorption spectrum equals Paper I's analytic Sobolev leg.

(3) is the one that matters. `sobolev.sobolev_leg.sobolev_attenuation` was
written for Paper I, is independently tested, and computes the same emergent
spectrum by integrating over impact parameter rather than by sampling rays.
Agreement pins the geometry, the resonance-plane solve and tau_S together.

Nothing here is a claim about kilonova spectra. It is an instrument
calibration, and the atom is synthetic.

Run:  python paper2/phase0/three_level_atom/run.py
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from sobolev.constants import C
from sobolev.sobolev_leg import sobolev_attenuation

from branching_mc import run_mc, three_level_atom

# Paper I's La II geometry: v_core = 1000 km/s, shell out to 5 r_core.
R_CORE = 8.64e12
R_OUT = 5.0 * R_CORE
T_EXP = 86400.0

LAM_13, LAM_32 = 4000e-8, 6000e-8          # cm
NU_13, NU_32 = C / LAM_13, C / LAM_32
F_OSC = 0.1

N_PACKETS = 200000


def pump_band():
    """Frequencies whose 1->3 resonance lies inside the shell for every ray."""
    ct = C * T_EXP
    return (NU_13 / (1.0 - 1.2 * R_CORE / ct),
            NU_13 / (1.0 - 0.8 * np.sqrt(R_OUT**2 - R_CORE**2) / ct))


def n_ground_for_tau(tau):
    probe = three_level_atom(NU_13, NU_32, F_OSC, 1.0, 1.0, 1.0)
    return tau / probe.tau(T_EXP)[0]


def measure_branching():
    print("\n1. BRANCHING RATIO -- measured fluorescent yield vs A32/(A31+A32)")
    print("   A31   A32 |  predicted   measured   (n_int)   deviation")
    rows = []
    n_g = n_ground_for_tau(3.0)
    for a31, a32 in [(4.0, 1.0), (2.0, 1.0), (1.0, 1.0), (1.0, 3.0), (1.0, 0.0)]:
        atom = three_level_atom(NU_13, NU_32, F_OSC, n_g, a31, a32)
        out = run_mc(atom, R_CORE, R_OUT, T_EXP, *pump_band(),
                     n_packets=N_PACKETS, seed=11, interaction="branch")
        n_int = int(out["n_interacted"])
        p_pred = a32 / (a31 + a32)
        p_meas = out["first_branch"][1] / n_int
        err = np.sqrt(max(p_pred * (1 - p_pred), 1e-12) / n_int)
        sigma = (p_meas - p_pred) / err if err > 0 else 0.0
        print(f"  {a31:4.1f}  {a32:4.1f} |   {p_pred:.4f}     {p_meas:.4f}   "
              f"({n_int:6d})   {sigma:+6.2f} sigma")
        rows.append({"a31": a31, "a32": a32, "p_predicted": p_pred,
                     "p_measured": p_meas, "n_interacted": n_int,
                     "sigma": float(sigma)})
    return rows


def measure_interaction_probability():
    print("\n2. INTERACTION PROBABILITY -- measured vs 1 - exp(-tau_S)")
    print("     tau |  predicted   measured   deviation")
    rows = []
    for tau in [0.1, 0.3, 1.0, 3.0, 10.0]:
        atom = three_level_atom(NU_13, NU_32, F_OSC, n_ground_for_tau(tau),
                                1.0, 0.0)
        out = run_mc(atom, R_CORE, R_OUT, T_EXP, *pump_band(),
                     n_packets=N_PACKETS, seed=12, interaction="absorb")
        assert out["n_offered"] == out["n_packets"], "band is not fully covered"
        p_pred = -np.expm1(-tau)
        p_meas = out["n_absorbed"] / out["n_packets"]
        err = np.sqrt(p_pred * (1 - p_pred) / out["n_packets"])
        sigma = (p_meas - p_pred) / err
        print(f"  {tau:6.1f} |   {p_pred:.4f}     {p_meas:.4f}   {sigma:+6.2f} sigma")
        rows.append({"tau": tau, "p_predicted": p_pred, "p_measured": p_meas,
                     "sigma": float(sigma)})
    return rows


def measure_spectrum_against_paper1(tau=1.5, n_bins=40, n_packets=600000):
    """The decisive check: pure absorption vs Paper I's analytic Sobolev leg.

    Compared as BIN AVERAGES on both sides. The trough edges are near-vertical
    -- the resonance plane leaves the shell over a fraction of a percent in
    frequency -- so comparing a histogram bin against the analytic value at the
    bin centre reports a large spurious disagreement in exactly those two bins
    while everything else agrees to <1%.
    """
    print(f"\n3. SPECTRUM vs Paper I analytic Sobolev leg (tau = {tau})")
    n_g = n_ground_for_tau(tau)
    atom = three_level_atom(NU_13, NU_32, F_OSC, n_g, 1.0, 0.0)

    nu_min, nu_max = 0.995 * NU_13, 1.025 * NU_13
    out = run_mc(atom, R_CORE, R_OUT, T_EXP, nu_min, nu_max,
                 n_packets=n_packets, seed=13, interaction="absorb")

    edges = np.linspace(nu_min, nu_max, n_bins + 1)
    h_in, _ = np.histogram(out["nu_launch"], bins=edges)
    h_out, _ = np.histogram(out["nu_out"], bins=edges)
    mc = h_out / h_in

    sub_f = (np.arange(48) + 0.5) / 48
    sub = edges[:-1, None] + sub_f[None, :] * np.diff(edges)[:, None]
    an = sobolev_attenuation(sub.ravel(), [(NU_13, F_OSC, 1.0)],
                             R_CORE, R_OUT, T_EXP, n_g, n_p=400)
    an = an.reshape(sub.shape).mean(axis=1)

    err = np.sqrt(np.maximum(mc * (1 - mc), 1e-12) / h_in)
    sigma = (mc - an) / np.maximum(err, 1e-4)
    print(f"   trough depth: MC {mc.min():.4f}  analytic {an.min():.4f}")
    print(f"   worst bin {abs(sigma).max():.2f} sigma, "
          f"RMS {np.sqrt((sigma**2).mean()):.2f} sigma "
          f"over {n_bins} bins at {n_packets//n_bins} packets/bin")
    print(f"   max |MC - analytic| = {np.abs(mc - an).max():.4f}")
    return {"tau": tau, "trough_mc": float(mc.min()),
            "trough_analytic": float(an.min()),
            "worst_sigma": float(abs(sigma).max()),
            "rms_sigma": float(np.sqrt((sigma**2).mean())),
            "max_abs_diff": float(np.abs(mc - an).max()),
            "nu_mid": (0.5 * (edges[:-1] + edges[1:])).tolist(),
            "mc": mc.tolist(), "analytic": an.tolist()}


def make_figure(spec, results_dir):
    """Emergent spectrum with and without the fluorescent channel."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_g = n_ground_for_tau(3.0)
    nu_min, nu_max = 0.9 * NU_32, 1.05 * NU_13
    edges = np.linspace(nu_min, nu_max, 260)
    mid = 0.5 * (edges[:-1] + edges[1:])
    lam = 1e8 * C / mid                       # Angstrom

    fig, ax = plt.subplots(2, 1, figsize=(7.2, 7.0),
                           gridspec_kw={"height_ratios": [1.3, 1]})

    for a32, colour, label in [(0.0, "0.35", r"$A_{32}=0$ (resonant only)"),
                               (1.0, "tab:red", r"$A_{32}=A_{31}$ (50% fluorescence)")]:
        atom = three_level_atom(NU_13, NU_32, F_OSC, n_g, 1.0, a32)
        out = run_mc(atom, R_CORE, R_OUT, T_EXP, nu_min, nu_max,
                     n_packets=400000, seed=14, interaction="branch")
        h_in, _ = np.histogram(out["nu_launch"], bins=edges)
        h_out, _ = np.histogram(out["nu_out"], bins=edges)
        ax[0].step(lam, h_out / np.maximum(h_in, 1), where="mid",
                   color=colour, lw=1.2, label=label)

    ax[0].axvline(1e8 * LAM_13, color="tab:blue", ls=":", lw=1)
    ax[0].axvline(1e8 * LAM_32, color="tab:orange", ls=":", lw=1)
    ax[0].text(1e8 * LAM_13, 2.35, r" $\lambda_{13}$ pump", color="tab:blue",
               fontsize=9, va="top")
    ax[0].text(1e8 * LAM_32, 2.35, r" $\lambda_{32}$ fluorescent",
               color="tab:orange", fontsize=9, va="top")
    ax[0].set_ylabel("emergent / launched")
    # Headroom for the P-Cygni emission peak: clipping it would hide the fact
    # that resonant scattering redistributes rather than removes photons.
    ax[0].set_ylim(0, 2.4)
    ax[0].legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.45, 0.92))
    ax[0].set_title("Three-level branching atom: photons removed at "
                    r"$\lambda_{13}$ reappear at $\lambda_{32}$", fontsize=10)

    ax[1].plot(1e8 * C / np.array(spec["nu_mid"]), spec["mc"], "o", ms=3.5,
               color="tab:red", label="branching MC (pure absorption)")
    ax[1].plot(1e8 * C / np.array(spec["nu_mid"]), spec["analytic"], "-",
               color="k", lw=1.2, label="Paper I analytic Sobolev leg")
    ax[1].set_xlabel(r"wavelength [$\AA$]")
    ax[1].set_ylabel("transmitted fraction")
    ax[1].legend(fontsize=9)
    ax[1].set_title(rf"Validation: $\tau_S={spec['tau']}$, worst bin "
                    rf"{spec['worst_sigma']:.1f}$\sigma$", fontsize=10)

    fig.tight_layout()
    path = results_dir / "fig_p2_three_level.png"
    fig.savefig(path, dpi=150)
    print(f"\nfigure -> {path}")
    return path


def main():
    print("P2-0B: three-level branching atom")
    print(f"  r_core {R_CORE:.3e} cm, r_out {R_OUT:.3e} cm, t_exp {T_EXP:.0f} s")
    print(f"  v_core {R_CORE/T_EXP/1e5:.0f} km/s, v_out {R_OUT/T_EXP/1e5:.0f} km/s")
    print(f"  lambda_13 {1e8*LAM_13:.0f} A, lambda_32 {1e8*LAM_32:.0f} A, "
          f"f = {F_OSC}")

    branching = measure_branching()
    interaction = measure_interaction_probability()
    spectrum = measure_spectrum_against_paper1()

    worst_branch = max(abs(r["sigma"]) for r in branching)
    worst_interact = max(abs(r["sigma"]) for r in interaction)
    ok = (worst_branch < 4 and worst_interact < 4
          and spectrum["worst_sigma"] < 5)
    print("\n" + "=" * 66)
    print(f"branching   worst {worst_branch:.2f} sigma")
    print(f"interaction worst {worst_interact:.2f} sigma")
    print(f"spectrum    worst {spectrum['worst_sigma']:.2f} sigma")
    print("VERDICT:", "instrument validated -- ready for Phase 1"
          if ok else "FAILED -- do not use for Phase 1")
    print("=" * 66)

    results = {"branching": branching, "interaction": interaction,
               "spectrum": spectrum, "validated": bool(ok),
               "geometry": {"r_core": R_CORE, "r_out": R_OUT, "t_exp": T_EXP,
                            "lambda_13_A": 1e8 * LAM_13,
                            "lambda_32_A": 1e8 * LAM_32, "f_osc": F_OSC}}
    (HERE / "results.json").write_text(json.dumps(results, indent=1))
    print(f"results -> {HERE / 'results.json'}")
    make_figure(spectrum, HERE)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
