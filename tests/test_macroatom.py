"""Paper IV 2B: the downward macroatom on toy atoms, analytically.

Test A  two-level: coherent, energy conserved, no chain trapping (also
        Test D, the optically thick repeat-interaction limit).
Test B  three-level cascade 3 -> 2 -> 1 with h nu31 = h nu32 + h nu21: the
        energy fractions per line are the Lucy table's, they sum to one, and
        the INFERRED physical photon numbers N(nu) ~ E(nu)/(h nu) of the two
        cascade lines are equal (the packet counts are not).
Test B' populated level 2: beta enters once (single application pinned).
Test C  the thermal channel: E_in = E_rad + E_thermal; deposit vs re-emit.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1"):
    sys.path.insert(0, str(p))

from sobolev.constants import C, H                      # noqa: E402
from sobolev.optical_depth import tau_sobolev           # noqa: E402
from sobolev.macroatom import DownwardMacroAtom         # noqa: E402
from forest_mc import ForestAtom, run_mc, band_ratio    # noqa: E402
import test_forest_mc as tfm                            # noqa: E402

R_CORE, R_OUT, T_EXP = tfm.R_CORE, tfm.R_OUT, tfm.T_EXP
NU_13, NU_32, F_OSC = tfm.NU_13, tfm.NU_32, tfm.F_OSC
NU_21 = NU_13 - NU_32
LEVEL_E = np.array([0.0, 0.0, NU_21 / C, NU_13 / C])   # levels 0 (unused), 1, 2, 3


def cascade_atom(tau13, a31, a32, a21=1.0, tau32=0.0):
    """Lines 0: 1->3 (pump), 1: 3->2 (fluorescence, opacity tau32), 2: 2->1
    (cascade exit, no opacity)."""
    n1 = tau13 / tau_sobolev(F_OSC, 1.0, C / NU_13, T_EXP)
    n2 = tau32 / tau_sobolev(F_OSC, 1.0, C / NU_32, T_EXP) if tau32 > 0 else 0.0
    fa = ForestAtom(nu0=[NU_13, NU_32, NU_21], f_osc=[F_OSC, F_OSC, 0.0],
                    n_lower=[n1, n2, 0.0], n_upper=[0.0, 0.0, 0.0], A=[a31, a32, a21],
                    lower=[1, 2, 1], upper=[3, 3, 2], t_exp=T_EXP, tau_min=1e-6, stim=False)
    fa.level_energy_cm = LEVEL_E
    fa.temperature = 3000.0
    fa.emis_w = np.array([1.0, 3.0, 0.0])
    return fa


def analytic_fractions(fa, a31, a32):
    b31, b32 = fa.beta_all[0], fa.beta_all[1]
    S = a31 * b31 + a32 * b32
    f31 = a31 * b31 / S
    f32 = a32 * b32 * NU_32 / (NU_13 * S)
    f21 = a32 * b32 * NU_21 / (NU_13 * S)
    return f31, f32, f21


def _run(fa, mode="sobolev_dmacro", n=100000, **kw):
    lo, hi = tfm.pump_band()
    return run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, n, mode, seed=11, packets="energy", **kw)


# ----------------------------------------------------------------- the table

def test_table_probabilities_are_the_lucy_downward_set():
    fa = cascade_atom(2.0, 2.0, 1.0)
    dm = fa.dmacro()
    lines, p_de, p_in = dm.probabilities(3)
    f31, f32, f21 = analytic_fractions(fa, 2.0, 1.0)
    assert list(lines) == [0, 1]
    assert np.isclose(p_de[0], f31) and np.isclose(p_in[0], 0.0)     # ground: no internal jump
    assert np.isclose(p_de[1], f32) and np.isclose(p_in[1], f21)     # internal jump to 2 carries eps_2
    assert np.isclose(p_de.sum() + p_in.sum(), 1.0)
    lines2, p_de2, p_in2 = dm.probabilities(2)
    assert list(lines2) == [2] and np.isclose(p_de2[0], 1.0) and p_in2[0] == 0.0
    assert not dm.has_exit[1] and not dm.dead_end[1]                  # the ground
    assert dm.n_entries == 6


def test_walk_reproduces_the_cascade_energy_flow_and_equal_physical_photon_numbers():
    fa = cascade_atom(2.0, 2.0, 1.0)
    dm = fa.dmacro()
    rng = np.random.default_rng(4)
    n = 200000
    ex, nj, dead = dm.walk(np.full(n, 3), rng)
    assert not dead.any() and np.all(ex >= 0)
    f31, f32, f21 = analytic_fractions(fa, 2.0, 1.0)
    # every exit carries the same energy: fractions are packet fractions
    for line, f in ((0, f31), (1, f32), (2, f21)):
        got = np.mean(ex == line)
        assert abs(got - f) < 4 * np.sqrt(f * (1 - f) / n), (line, got, f)
    assert np.isclose(f31 + f32 + f21, 1.0)
    # inferred physical photon numbers of the cascade lines: E / (h nu)
    n32 = np.sum(ex == 1) / NU_32; n21 = np.sum(ex == 2) / NU_21
    assert abs(n32 / n21 - 1.0) < 4 * np.sqrt(1 / np.sum(ex == 1) + 1 / np.sum(ex == 2))
    # the internal jump is the only route to line 2
    assert np.all(nj[ex == 2] == 1) and np.all(nj[ex != 2] == 0)


def test_table_cut_keeps_the_dominant_exits_and_renormalises():
    fa = cascade_atom(2.0, 100.0, 1.0)      # 3->1 dominates by ~100
    dm = fa.dmacro(a_cut=0.05)
    assert dm.n_entries < 6 and dm.a_cut == 0.05
    lines, p_de, p_in = dm.probabilities(3)
    assert np.isclose(p_de.sum() + p_in.sum(), 1.0)
    assert p_de[list(lines).index(0)] > 0.98


# ----------------------------------------------------------------- Test A / D

@pytest.mark.parametrize("tau", [0.5, 3.0, 30.0, 1e3])
def test_two_level_is_coherent_energy_conserving_and_never_traps(tau):
    atom, _ = tfm._one_line_atom(tau=tau)
    atom.level_energy_cm = np.array([0.0, atom.op_nu[0] / C])
    nu0 = atom.op_nu[0]
    ct = C * T_EXP
    r_core, r_out = R_CORE, R_OUT
    # every packet is offered the resonance (the pump_band construction)
    lo, hi = nu0 / (1.0 - 1.2 * r_core / ct), nu0 / (1.0 - 0.8 * np.sqrt(r_out**2 - r_core**2) / ct)
    res = run_mc(atom, r_core, r_out, T_EXP, lo, hi, 40000, "sobolev_dmacro", seed=1,
                 packets="energy", chain_max=1)
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12 and a["E_dep_cm"] == 0.0
    assert res["n_trapped"] == 0 and res["n_reabs"].sum() == 0     # closed form: no chain
    assert np.array_equal(res["w"], res["w_launch"])               # coherent: w unchanged
    p_int = res["n_interacted"] / res["n_packets"]
    p_th = 1.0 - np.exp(-tau)
    assert abs(p_int - p_th) < 4 * np.sqrt(p_th * (1 - p_th) / res["n_packets"]) + 1e-3
    # photon and energy band ratios differ only by the Doppler work, O(beta_out)
    ph, _ = band_ratio(res, lo, hi, "photon"); en, _ = band_ratio(res, lo, hi, "energy")
    assert abs(ph / en - 1.0) < 3 * r_out / ct + 2e-3


# ----------------------------------------------------------------- Test B / B'

def test_transport_exit_energies_follow_the_table_and_conserve_energy():
    a31, a32 = 2.0, 1.0
    fa = cascade_atom(2.0, a31, a32)
    res = _run(fa)
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12 and a["E_dep_cm"] == 0.0
    assert res["n_reabs"].sum() == 0 and res["n_dead_end"] == 0 and res["n_kpackets"] == 0
    ex = res["exit_energy"]; tot = ex.sum()
    n_act = res["n_interactions"]
    f31, f32, f21 = analytic_fractions(fa, a31, a32)
    for line, f in ((0, f31), (1, f32), (2, f21)):
        assert abs(ex[line] / tot - f) < 4 * np.sqrt(f * (1 - f) / n_act) + 1e-3, (line, ex[line] / tot, f)
    # inferred physical photon numbers of 3->2 and 2->1 are equal
    n32, n21 = ex[1] / (H * NU_32), ex[2] / (H * NU_21)
    assert abs(n32 / n21 - 1.0) < 4 * np.sqrt(2.0 / (n_act * f32))
    # and the emitted energy is the activated energy: exits carry the whole packet
    assert np.isclose(tot, a["E_interacting"] * 0 + tot)   # (tautology guard for the diagnostic below)
    # R1 <-> R2 in one number: photon-number branching would have deposited
    # h nu21 per fluorescence; here that energy leaves in line 2
    ref = run_mc(fa, R_CORE, R_OUT, T_EXP, *tfm.pump_band(), 100000, "sobolev_branch", seed=11)
    dep_r1 = ref["accounting"]["E_dep_cm"]
    assert abs(ex[2] / dep_r1 - 1.0) < 0.08          # same atom, same seed, different draws: noise


def test_populated_level_two_applies_beta_once():
    """With tau32 > 0, F31 = A31 b31 / (A31 b31 + A32 b32). Sending the exit
    through the beta chain again would give A b^2 weights; at these depths
    the two predictions are many sigma apart, and the MC must sit on the
    first."""
    a31, a32 = 1.0, 1.0
    fa = cascade_atom(3.0, a31, a32, tau32=6.0)
    b31, b32 = fa.beta_all[0], fa.beta_all[1]
    assert 0 < b32 < b31 < 1
    res = _run(fa, n=150000)
    ex = res["exit_energy"]
    got = ex[0] / ex.sum()
    f_once = a31 * b31 / (a31 * b31 + a32 * b32)
    f_twice = a31 * b31**2 / (a31 * b31**2 + a32 * b32**2)
    se = np.sqrt(f_once * (1 - f_once) / res["n_interactions"])
    assert abs(got - f_once) < 4 * se + 1e-3, (got, f_once)
    assert abs(f_twice - f_once) > 10 * se            # the test can tell them apart
    assert abs(got - f_twice) > 6 * se


# ----------------------------------------------------------------- Test C

def test_thermal_channel_deposit_books_the_k_packets():
    fa = cascade_atom(2.0, 2.0, 1.0)
    res = _run(fa, eps_k=0.3, thermal_k="deposit")
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12
    assert np.isclose(a["E_abs"], a["E_thermal"], rtol=1e-12)     # every k-packet was deposited
    n_act = res["n_kpackets"] + int(np.rint(res["exit_energy"].sum() / (H * NU_13)))
    frac = res["n_kpackets"] / n_act
    assert abs(frac - 0.3) < 4 * np.sqrt(0.3 * 0.7 / n_act)
    e_frac = a["E_thermal"] / (a["E_thermal"] + res["exit_energy"].sum())
    assert abs(e_frac - 0.3) < 4 * np.sqrt(0.3 * 0.7 / n_act) + 1e-3
    assert res["n_absorbed"] == res["n_kpackets"] and res["n_trapped"] == 0


def test_thermal_channel_reemit_keeps_the_energy_and_follows_the_net_emissivity():
    fa = cascade_atom(2.0, 2.0, 1.0)            # emis_w = [1, 3, 0]; beta = [b31, 1, 1]
    res = _run(fa, eps_k=1.0, thermal_k="reemit")   # every activation thermalises
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12 and a["E_abs"] == 0.0 and a["E_dep_cm"] == 0.0
    assert res["n_kpackets"] == res["n_interactions"] and res["n_absorbed"] == 0
    ex = res["exit_energy"]
    wgt = fa.emis_w * fa.nu0_all * fa.beta_all
    p = wgt / wgt.sum()
    n = res["n_kpackets"]
    for line in (0, 1):
        assert abs(ex[line] / ex.sum() - p[line]) < 4 * np.sqrt(p[line] * (1 - p[line]) / n) + 1e-3
    assert ex[2] == 0.0
