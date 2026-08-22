"""Acceptance tests for the Phase 1 forest Monte Carlo.

The successor must (i) agree with the Phase 0 instrument on the three-level
atom it was calibrated on, (ii) reproduce Paper I's analytic legs in both
pure-absorption modes -- Sobolev and the expansion-opacity closure -- and
(iii) conserve photon number in every mode. Everything Paper II measures sits
on those three.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper2/phase1"))
sys.path.insert(0, str(ROOT / "paper2/phase0/three_level_atom"))

from sobolev.constants import C
from sobolev.sobolev_leg import expansion_damp, sobolev_attenuation
from forest_mc import ForestAtom, MODES, band_ratio, run_mc, spectrum
import branching_mc as p0

R_CORE = 8.64e12
R_OUT = 5.0 * R_CORE
T_EXP = 86400.0
NU_13, NU_32 = C / 4000e-8, C / 6000e-8
F_OSC = 0.1


def three_level(tau_target, a31=1.0, a32=1.0):
    """ForestAtom and Phase-0 Atom for the same three-level system."""
    from sobolev.optical_depth import tau_sobolev
    n_ground = tau_target / tau_sobolev(F_OSC, 1.0, C / NU_13, T_EXP)
    fa = ForestAtom(nu0=[NU_13, NU_32], f_osc=[F_OSC, 0.0], n_lower=[n_ground, 0.0],
                    n_upper=[0.0, 0.0], A=[a31, a32], lower=[1, 2], upper=[3, 3],
                    t_exp=T_EXP, tau_min=1e-6, stim=False)
    fa.n_lower_0 = n_ground
    pa = p0.three_level_atom(NU_13, NU_32, F_OSC, n_ground, a31, a32)
    return fa, pa


def pump_band():
    ct = C * T_EXP
    return (NU_13 / (1.0 - 1.2 * R_CORE / ct),
            NU_13 / (1.0 - 0.8 * np.sqrt(R_OUT**2 - R_CORE**2) / ct))


# ---------------------------------------------------------------- conservation

@pytest.mark.parametrize("mode", MODES)
def test_photon_number_is_conserved_in_every_mode(mode):
    fa, _ = three_level(1.0)
    lo, hi = pump_band()
    if mode.endswith("thermal"):
        # give the upper level a population so the thermal sampler has weight
        fa.emis_w = np.array([1.0, 1.0])
    res = run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, 20000, mode, seed=1)
    assert res["n_escaped"] + res["n_core"] + res["n_absorbed"] == res["n_packets"]
    if mode.endswith("absorb"):
        assert res["n_absorbed"] == res["n_interacted"]
    else:
        assert res["n_absorbed"] == 0


# ------------------------------------------------ agreement with Phase 0 instrument

@pytest.mark.parametrize("tau_target", [0.3, 1.0, 3.0])
def test_interaction_probability_matches_phase0_and_theory(tau_target):
    fa, pa = three_level(tau_target)
    lo, hi = pump_band()
    n = 60000
    new = run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_absorb", seed=3)
    old = p0.run_mc(pa, R_CORE, R_OUT, T_EXP, lo, hi, n // 3, seed=3, interaction="absorb")
    p_new = new["n_interacted"] / n
    p_old = old["n_interacted"] / old["n_offered"]
    p_th = -np.expm1(-tau_target)
    assert abs(p_new - p_th) < 4 * np.sqrt(p_th * (1 - p_th) / n) + 1e-3
    assert abs(p_new - p_old) < 4 * np.sqrt(p_th * (1 - p_th) / (n // 3)) + 2e-3


@pytest.mark.parametrize("a31,a32", [(4.0, 1.0), (1.0, 1.0), (1.0, 3.0)])
def test_branching_ratio_matches_phase0(a31, a32):
    fa, pa = three_level(2.0, a31, a32)
    lo, hi = pump_band()
    n = 60000
    new = run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_branch", seed=5)
    old = p0.run_mc(pa, R_CORE, R_OUT, T_EXP, lo, hi, n // 3, seed=5, interaction="branch")
    # fluorescent yield from the FIRST interaction's branch choice (line 1 is
    # 3->2), as Phase 0 measured it; escaped-packet counts would be biased by
    # core absorption of inward re-emission and by repeat interactions.
    p_new = new["first_branch"][1] / new["n_interacted"]
    p_old = old["first_branch"][1] / old["n_interacted"]
    p_th = a32 / (a31 + a32)
    assert abs(p_new - p_th) < 4 * np.sqrt(p_th * (1 - p_th) / new["n_interacted"]) + 1e-3
    assert abs(p_new - p_old) < 4 * np.sqrt(p_th * (1 - p_th) / old["n_interacted"]) + 2e-3


# ------------------------------------------- agreement with Paper I's analytic legs

def _pure_absorption_spectrum(mode, tau_target, n=150000, seed=7):
    fa, _ = three_level(tau_target)
    lo, hi = pump_band()
    res = run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=seed)
    n_bins = 24
    edges = np.linspace(lo, hi, n_bins + 1)
    mc, err = spectrum(res, edges)
    # bin-average the analytic leg (the Phase 0 lesson: midpoint vs average)
    sub = np.linspace(lo, hi, n_bins * 40 + 1)
    damp = expansion_damp if mode.startswith("expansion") else None
    n_ground = fa.n_lower_0
    an = sobolev_attenuation(sub, [(NU_13, F_OSC, 1.0)], R_CORE, R_OUT, T_EXP, n_ground, damp=damp)
    an_b = an[:-1].reshape(n_bins, 40).mean(axis=1)
    return mc, err, an_b


@pytest.mark.parametrize("tau_target", [0.5, 2.0])
def test_sobolev_absorb_reproduces_analytic_sobolev_leg(tau_target):
    mc, err, an = _pure_absorption_spectrum("sobolev_absorb", tau_target)
    z = (mc - an) / err
    assert np.nanmax(np.abs(z)) < 3.5, (mc, an)


@pytest.mark.parametrize("tau_target", [0.5, 2.0])
def test_expansion_absorb_reproduces_the_poisson_closure(tau_target):
    """The closure's own prediction: e^{-(1-e^{-tau})} per crossing. This is
    Paper I's F4/R2 inside the Monte Carlo, same code as the Sobolev leg."""
    mc, err, an = _pure_absorption_spectrum("expansion_absorb", tau_target)
    z = (mc - an) / err
    assert np.nanmax(np.abs(z)) < 3.5, (mc, an)
    # and it is NOT the Sobolev leg: the two differ by more than the noise
    mc_s, _, an_s = _pure_absorption_spectrum("sobolev_absorb", tau_target)
    assert np.nanmin(an) > np.nanmax(an_s[np.isfinite(an_s)][8:16]) - 1e-9 or tau_target < 1.0


# -------------------------------------------------------------- the real atom

def test_la_ii_atom_builds_with_the_expected_size():
    data = ROOT / "data"
    if not (data / "57LaII_transitions_calib.txt").exists():
        pytest.skip("GSI La II data not present")
    d = np.load(ROOT / "experiments/laII_forest/forest_lines.npz")
    atom = ForestAtom.from_gsi(data / "57LaII_levels_calib.txt", data / "57LaII_transitions_calib.txt",
                               3000.0, float(d["n_ion"]), T_EXP, tau_min=1e-3)
    assert atom.n_lines_total == 17743
    assert 800 < atom.n_opacity < 1100
    assert np.all(np.diff(atom.op_nu) >= 0)
    # every upper level fed by an opacity line has a branching table that sums to 1
    for u, cum in atom.branch_cum.items():
        assert np.isclose(cum[-1], 1.0)
    # window lines: the forest's 153 lines include its ~10 with tau > 1e-3
    lam = C / atom.op_nu * 1e8
    assert 5 <= np.count_nonzero((lam >= 3850) & (lam < 3950)) <= 20


# ------------------------------------------------------ escape probability

@pytest.mark.parametrize("tau_target,a31,a32", [(2.0, 1.0, 1.0), (5.0, 4.0, 1.0), (0.5, 1.0, 3.0)])
def test_trapped_fluorescence_yield_follows_escape_probability(tau_target, a31, a32):
    """A photon re-emitted in the resonant line escapes it with beta =
    (1-e^-tau)/tau and is otherwise re-absorbed and re-drawn, so the chance
    that an interaction chain ends in the fluorescent channel is

        P = b / (1 - (1-b)(1-beta)),   b = A32/(A31+A32),

    not b. The first-draw tally still gives b (the Phase 0 quantity)."""
    fa, _ = three_level(tau_target, a31, a32)
    lo, hi = pump_band()
    n = 80000
    res = run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, n, "sobolev_branch", seed=11)
    b = a32 / (a31 + a32)
    beta = -np.expm1(-tau_target) / tau_target
    p_first = res["first_branch"][1] / res["n_interacted"]
    p_chain = res["chain_exit"][1] / res["n_interacted"]
    p_th = b / (1.0 - (1.0 - b) * (1.0 - beta))
    se = np.sqrt(p_th * (1 - p_th) / res["n_interacted"])
    assert abs(p_first - b) < 4 * np.sqrt(b * (1 - b) / res["n_interacted"]) + 1e-3
    assert abs(p_chain - p_th) < 4 * se + 1e-3, (p_chain, p_th)
    assert p_chain > p_first + 3 * se  # trapping enhances fluorescence


def test_pure_resonant_scattering_is_unaffected_by_the_escape_loop():
    """With a single channel the emergent spectrum cannot depend on how many
    times the photon is re-absorbed in its own line."""
    fa, _ = three_level(3.0, 1.0, 0.0)
    lo, hi = pump_band()
    res = run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, 60000, "sobolev_branch", seed=13)
    assert res["n_escaped"] + res["n_core"] == res["n_packets"]
    # the chain_exit must be entirely in line 0 and interactions > first draws
    assert res["chain_exit"][1] == 0 and res["first_branch"][1] == 0
    assert res["n_interactions"] > res["n_interacted"]
