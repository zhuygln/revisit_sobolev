"""Validation of Paper II's branching Monte Carlo (P2-0B).

The MC is a new instrument, so the burden is to show it reproduces things
already known by other means before it is used to measure anything new.

Four independent handles, weakest to strongest:

  1. photon-number conservation -- catches packets vanishing;
  2. interaction probability = 1 - exp(-tau_S) -- catches the sampler;
  3. branching fractions = A_uj / sum A_uj -- the quantity P2-0B exists for;
  4. the pure-absorption spectrum against `sobolev.sobolev_leg` -- the real
     test. That analytic function was written independently, is used in
     Paper I, and is covered by its own tests, so agreement pins the geometry,
     the impact-parameter weighting, the resonance-plane solve AND tau_S all
     at once. A transport bug that survives (1)-(3) will not survive this.

Every test is seeded. A tolerance here is a statement about Monte Carlo noise
at the stated packet count, not a fudge factor: each is set from the binomial
standard error of that specific measurement.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper2/phase0/three_level_atom"))

from sobolev.constants import C
from sobolev.sobolev_leg import sobolev_attenuation

from branching_mc import Atom, run_mc, three_level_atom

# Paper I's La II geometry, so the regime is the one already characterised.
R_CORE = 8.64e12
R_OUT = 5.0 * R_CORE
T_EXP = 86400.0

LAM_13 = 4000e-8          # pump line, cm
LAM_32 = 6000e-8          # fluorescent escape channel, cm
NU_13 = C / LAM_13
NU_32 = C / LAM_32


def _atom(a31=1.0, a32=1.0, n_ground=1.0e5, f_osc=0.1):
    return three_level_atom(NU_13, NU_32, f_osc, n_ground, a31, a32)


def _wide_band():
    """A launch band straddling both lines, with room on either side."""
    return 0.9 * NU_32, 1.1 * NU_13


def _pump_band():
    """A launch band in which EVERY packet is offered the 1->3 resonance.

    The resonance plane sits at z = c t (1 - nu_0/nu), so only a narrow slice
    of frequencies resonates inside the shell at all: here v/c runs from
    0.0033 at the core to 0.0167 at r_out, so the window is ~1.4% wide. Over
    `_wide_band` that is 3% of the packets, which is why the first version of
    the interaction-probability test could not measure anything.

    A ray at impact parameter p enters the shell at z_lo = sqrt(r_c^2 - p^2)
    and leaves at z_hi = sqrt(r_out^2 - p^2). Requiring r_core < z_res <
    sqrt(r_out^2 - r_core^2) therefore guarantees a crossing for EVERY ray,
    whatever its p -- so the offered fraction is exactly 1 and the measured
    interacted fraction is a clean estimate of 1 - exp(-tau).
    """
    ct = C * T_EXP
    z_lo = 1.2 * R_CORE
    z_hi = 0.8 * np.sqrt(R_OUT**2 - R_CORE**2)
    return NU_13 / (1.0 - z_lo / ct), NU_13 / (1.0 - z_hi / ct)


# --------------------------------------------------------------------------
# 1. Conservation
# --------------------------------------------------------------------------


def test_photon_number_is_conserved_between_escape_and_core():
    """Branching scatters packets; it must not create or destroy them.

    The core is the only sink, so N_launched = N_escaped + N_core_absorbed
    exactly -- an integer identity, not a statistical one.
    """
    out = run_mc(_atom(), R_CORE, R_OUT, T_EXP, *_wide_band(),
                 n_packets=20000, seed=1, interaction="branch")
    assert out["nu_out"].size + out["n_core"] == out["n_packets"]
    # And the test would be vacuous if nothing ever came back to the core.
    assert out["n_core"] > 0


def test_absorb_mode_destroys_exactly_the_interacting_packets():
    out = run_mc(_atom(), R_CORE, R_OUT, T_EXP, *_wide_band(),
                 n_packets=20000, seed=2, interaction="absorb")
    assert (out["nu_out"].size + out["n_core"] + out["n_absorbed"]
            == out["n_packets"])
    assert out["n_absorbed"] > 0


# --------------------------------------------------------------------------
# 2. Interaction probability
# --------------------------------------------------------------------------


@pytest.mark.parametrize("tau_target", [0.1, 1.0, 5.0])
def test_interaction_probability_matches_one_minus_exp_minus_tau(tau_target):
    """Among packets actually offered the resonance, the interacted fraction
    must be 1 - exp(-tau_S).

    Run in 'absorb' mode so that "interacted" is unambiguous: exactly one
    interaction per packet, no re-emission to confuse the tally.
    """
    # Choose n_ground to hit the requested tau, using the code's own tau so
    # that this tests the sampler and not the prefactor (which test 4 pins).
    probe = _atom(n_ground=1.0)
    n_ground = tau_target / probe.tau(T_EXP)[0]
    atom = _atom(n_ground=n_ground)
    assert atom.tau(T_EXP)[0] == pytest.approx(tau_target, rel=1e-12)

    n = 40000
    out = run_mc(atom, R_CORE, R_OUT, T_EXP, *_pump_band(),
                 n_packets=n, seed=3, interaction="absorb")

    # The band is constructed so every ray crosses the plane; if that ever
    # stops holding, the denominator below is silently wrong.
    offered = out["n_offered"]
    assert offered == n, f"{n - offered} packets missed the resonance"
    frac = out["n_absorbed"] / offered
    expected = 1.0 - np.exp(-tau_target)
    err = np.sqrt(expected * (1 - expected) / offered)
    assert abs(frac - expected) < 4 * err


# --------------------------------------------------------------------------
# 3. Branching -- what P2-0B is for
# --------------------------------------------------------------------------


@pytest.mark.parametrize("a31,a32", [(4.0, 1.0), (1.0, 1.0), (1.0, 3.0)])
def test_first_interaction_branching_matches_analytic_ratios(a31, a32):
    """P(3->2) measured = A32 / (A31 + A32), to binomial noise."""
    atom = _atom(a31=a31, a32=a32)
    out = run_mc(atom, R_CORE, R_OUT, T_EXP, *_pump_band(),
                 n_packets=40000, seed=4, interaction="branch")

    counts = out["first_branch"]
    n_int = counts.sum()
    assert n_int == out["n_interacted"]
    assert n_int > 3000, "too few interactions to measure a ratio"

    p32 = a32 / (a31 + a32)
    measured = counts[1] / n_int
    err = np.sqrt(p32 * (1 - p32) / n_int)
    assert abs(measured - p32) < 4 * err


def test_pure_resonant_atom_emits_nothing_at_the_fluorescent_line():
    """A32 = 0 must switch the fluorescent channel off completely.

    Not a statistical statement: zero packets may emerge anywhere near
    lambda_32. This is the degenerate limit that would catch a branching table
    that ignores its rates.
    """
    atom = _atom(a31=1.0, a32=0.0)
    out = run_mc(atom, R_CORE, R_OUT, T_EXP, *_wide_band(),
                 n_packets=20000, seed=5, interaction="branch")

    assert out["first_branch"][1] == 0
    assert out["first_branch"][0] == out["n_interacted"] > 0

    # No emergent packet within +-5% of nu_32 beyond those simply launched
    # there and never touched -- i.e. the count must equal the launch count in
    # that window, since nothing can be added to or removed from it.
    lo, hi = 0.95 * NU_32, 1.05 * NU_32
    n_out = np.sum((out["nu_out"] > lo) & (out["nu_out"] < hi))
    n_in = np.sum((out["nu_launch"] > lo) & (out["nu_launch"] < hi))
    assert n_out == n_in


def test_fluorescent_packets_appear_only_redward_when_branching_is_on():
    """With A32 > 0 the nu_32 window must GAIN packets, and by the amount the
    branching ratio predicts (up to those that fall back into the core)."""
    a31, a32 = 1.0, 1.0
    atom = _atom(a31=a31, a32=a32)
    out = run_mc(atom, R_CORE, R_OUT, T_EXP, *_pump_band(),
                 n_packets=40000, seed=6, interaction="branch")

    lo, hi = 0.95 * NU_32, 1.05 * NU_32
    n_out = np.sum((out["nu_out"] > lo) & (out["nu_out"] < hi))
    n_in = np.sum((out["nu_launch"] > lo) & (out["nu_launch"] < hi))
    assert n_in == 0, "pump band must not overlap the fluorescent window"
    gained = n_out - n_in
    predicted = out["first_branch"][1]

    assert gained > 0
    # Fluorescent packets escape freely (level 2 is unpopulated, so nu_32
    # carries no opacity) except those re-emitted inward into the core.
    assert gained <= predicted
    assert gained > 0.5 * predicted


# --------------------------------------------------------------------------
# 4. The real test: transport against Paper I's analytic Sobolev leg
# --------------------------------------------------------------------------


def _analytic_binned(edges, atom, n_ground, n_sub=48):
    """Bin-AVERAGED analytic transmission, which is what the MC measures.

    A histogram bin reports the mean transmission over the bin, and packets
    are launched uniformly in nu, so the like-for-like analytic quantity is
    the uniform average of T(nu) across the bin -- not T evaluated at the bin
    centre.

    This is not a detail. The trough has near-vertical edges (the resonance
    plane leaves the shell over a fraction of a percent in frequency), and at
    those edges the midpoint value and the bin average differ by tens of
    percent while agreeing to <1% everywhere else. The first version of this
    test compared midpoint-to-average and reported a 27-sigma "disagreement"
    in exactly the two edge bins -- the same failure mode as the red-margin
    straddle that once manufactured a 7% error in Paper I's breadth sweep.
    """
    lo, hi = edges[:-1], edges[1:]
    frac = (np.arange(n_sub) + 0.5) / n_sub
    sub = lo[:, None] + frac[None, :] * (hi - lo)[:, None]
    t = sobolev_attenuation(
        sub.ravel(), [(NU_13, atom.lines[0]["f_osc"], 1.0)],
        R_CORE, R_OUT, T_EXP, n_ground, n_p=400,
    )
    return t.reshape(sub.shape).mean(axis=1)


@pytest.mark.parametrize("tau_target", [0.3, 1.5])
def test_absorption_spectrum_matches_analytic_sobolev_attenuation(tau_target):
    """Pure absorption must reproduce `sobolev_attenuation` frequency by
    frequency.

    That function integrates exp(-sum tau) over impact parameter across the
    core disk, counting only resonance planes between the core surface and the
    outer boundary. The MC arrives at the same quantity by an entirely
    different route -- sampling rays and resonance crossings -- so agreement
    is a real check on the geometry, not a tautology.
    """
    probe = _atom(n_ground=1.0)
    n_ground = tau_target / probe.tau(T_EXP)[0]
    atom = _atom(n_ground=n_ground, a32=0.0)

    nu_min, nu_max = 0.995 * NU_13, 1.025 * NU_13
    n = 400000
    out = run_mc(atom, R_CORE, R_OUT, T_EXP, nu_min, nu_max,
                 n_packets=n, seed=7, interaction="absorb")

    edges = np.linspace(nu_min, nu_max, 41)
    h_in, _ = np.histogram(out["nu_launch"], bins=edges)
    h_out, _ = np.histogram(out["nu_out"], bins=edges)
    mid = 0.5 * (edges[:-1] + edges[1:])
    transmitted = h_out / h_in

    analytic = _analytic_binned(edges, atom, n_ground)

    # The trough must actually be there, or this compares two flat lines.
    assert analytic.min() < 0.85

    err = np.sqrt(np.maximum(transmitted * (1 - transmitted), 1e-12) / h_in)
    resid = np.abs(transmitted - analytic) / np.maximum(err, 1e-4)
    assert resid.max() < 5.0, (
        f"worst bin {resid.argmax()}: MC {transmitted[resid.argmax()]:.4f} "
        f"vs analytic {analytic[resid.argmax()]:.4f}"
    )


def test_analytic_comparison_would_fail_on_a_wrong_tau():
    """A control on the control.

    Paper I twice had a 'probe' that passed because it sat outside the region
    it meant to test. So: perturb tau by 20% and confirm the comparison above
    actually rejects it. If this test fails, the previous one proves nothing.
    """
    probe = _atom(n_ground=1.0)
    n_ground = 1.5 / probe.tau(T_EXP)[0]
    atom = _atom(n_ground=n_ground, a32=0.0)

    nu_min, nu_max = 0.995 * NU_13, 1.025 * NU_13
    n = 400000
    out = run_mc(atom, R_CORE, R_OUT, T_EXP, nu_min, nu_max,
                 n_packets=n, seed=7, interaction="absorb")

    edges = np.linspace(nu_min, nu_max, 41)
    h_in, _ = np.histogram(out["nu_launch"], bins=edges)
    h_out, _ = np.histogram(out["nu_out"], bins=edges)
    mid = 0.5 * (edges[:-1] + edges[1:])
    transmitted = h_out / h_in

    wrong = _analytic_binned(edges, atom, 1.2 * n_ground)
    err = np.sqrt(np.maximum(transmitted * (1 - transmitted), 1e-12) / h_in)
    resid = np.abs(transmitted - wrong) / np.maximum(err, 1e-4)
    assert resid.max() > 5.0


# --------------------------------------------------------------------------
# Guards on the atom itself
# --------------------------------------------------------------------------


def test_unpopulated_line_carries_no_opacity_but_remains_a_branch():
    atom = _atom()
    assert atom.tau(T_EXP)[1] == 0.0
    assert 1 in atom.branches[3][0]


def test_branch_sampler_respects_cumulative_boundaries():
    atom = Atom([
        {"nu0": NU_13, "lower": 1, "upper": 3, "f_osc": 0.1,
         "n_lower": 1e5, "A": 3.0},
        {"nu0": NU_32, "lower": 2, "upper": 3, "f_osc": 0.0,
         "n_lower": 0.0, "A": 1.0},
    ])
    assert atom.branch(3, 0.0) == 0
    assert atom.branch(3, 0.74) == 0
    assert atom.branch(3, 0.76) == 1
    assert atom.branch(3, 0.999) == 1


def test_fluorescent_line_must_be_redder():
    with pytest.raises(ValueError):
        three_level_atom(NU_32, NU_13, 0.1, 1e5, 1.0, 1.0)
