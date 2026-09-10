"""Paper IV test 2A.0: the packet energy contract in an expanding flow,
before any fluorescence test.

A packet is (E_cm, nu_cm, nu_lab, mu, r, t). E_lab is conserved in free
flight; E_cm is conserved at every atomic interaction; the two differ by the
local Doppler factor. `sobolev/energy_packets.py` states the contract and
`run_mc(packets="energy")` implements it; this file checks both against the
analytic transformations, classical and worldline.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1"):
    sys.path.insert(0, str(p))

from sobolev.constants import C, H                     # noqa: E402
from sobolev import energy_packets as ep               # noqa: E402
from forest_mc import run_mc                           # noqa: E402
import test_forest_mc as tfm                           # noqa: E402

T_EXP = 86400.0
CT = C * T_EXP


@pytest.mark.parametrize("relativity", [None, "worldline"])
def test_doppler_factor_is_the_analytic_transformation(relativity):
    rng = np.random.default_rng(0)
    r = rng.uniform(0.02, 0.1, 200) * CT          # beta up to 0.1
    mu = rng.uniform(-1.0, 1.0, 200)
    d = ep.doppler_factor(r, mu, CT, relativity)
    beta = r / CT
    if relativity is None:
        assert np.allclose(d, 1.0 / (1.0 - beta * mu), rtol=1e-14)
    else:
        gamma = 1.0 / np.sqrt(1.0 - beta**2)
        assert np.allclose(d, 1.0 / (gamma * (1.0 - beta * mu)), rtol=1e-14)
        # worldline differs from classical at second order only
        d_cl = ep.doppler_factor(r, mu, CT, None)
        assert np.all(np.abs(d / d_cl - 1.0) < 1.5 * beta**2)
    # comoving <-> lab round trip
    nu = rng.uniform(4e14, 1e15, 200)
    nu_cm = ep.comoving_frequency(nu, r, mu, CT, relativity)
    assert np.allclose(nu_cm * d, nu, rtol=1e-13)
    assert np.allclose(ep.lab_frequency(nu_cm, r, mu, CT, relativity), nu, rtol=1e-13)


@pytest.mark.parametrize("relativity", [None, "worldline"])
def test_coherent_scattering_conserves_e_cm_and_moves_e_lab_by_the_doppler_ratio(relativity):
    """The contract's statement for one scattering at (r, mu_in -> mu_out)."""
    rng = np.random.default_rng(1)
    r = rng.uniform(0.02, 0.1, 500) * CT
    mu_in, mu_out = rng.uniform(-1, 1, 500), rng.uniform(-1, 1, 500)
    w = rng.uniform(0.5, 2.0, 500); nu_in = rng.uniform(4e14, 1e15, 500)
    e_lab_in, e_cm_in = ep.energies(w, nu_in, r, mu_in, CT, relativity)
    e_lab_out = ep.scattered_lab_energy(e_lab_in, r, mu_in, mu_out, CT, relativity)
    nu_out = ep.lab_frequency(ep.comoving_frequency(nu_in, r, mu_in, CT, relativity), r, mu_out, CT, relativity)
    _, e_cm_out = ep.energies(e_lab_out / (H * nu_out), nu_out, r, mu_out, CT, relativity)
    assert np.allclose(e_cm_out, e_cm_in, rtol=1e-12)
    d_in = ep.doppler_factor(r, mu_in, CT, relativity); d_out = ep.doppler_factor(r, mu_out, CT, relativity)
    assert np.allclose(e_lab_out / e_lab_in, d_out / d_in, rtol=1e-12)


def test_reweight_is_exact_for_coherent_scattering():
    w = np.array([0.3, 1.0, 7.5]); nu = np.array([4e14, 7.5e14, 1e15])
    assert np.array_equal(ep.reweight(w, nu, nu), w)          # x / x == 1.0 exactly
    assert np.allclose(ep.reweight(w, nu, nu / 2) * (nu / 2), w * nu, rtol=1e-15)


@pytest.mark.parametrize("relativity", [None, "worldline"])
def test_free_flight_conserves_lab_energy(relativity):
    """No resonance is reachable (launch band below the only line, and the
    comoving frequency only falls): every packet keeps nu_lab and w, the
    Doppler work is zero, and the identity closes with no deposit."""
    atom, _ = tfm._one_line_atom(tau=5.0)
    nu0 = atom.op_nu[0]
    res = run_mc(atom, CT / 30.0, CT / 10.0, T_EXP, nu0 * 0.90, nu0 * 0.98, 20000,
                 "sobolev_absorb", seed=3, relativity=relativity, packets="energy")
    assert res["n_interactions"] == 0
    assert np.array_equal(res["nu_final"], res["nu_launch"])
    assert np.array_equal(res["w"], res["w_launch"])
    a = res["accounting"]
    assert a["W"] == 0.0 and a["E_dep_cm"] == 0.0 and a["E_dep_lab"] == 0.0
    assert abs(a["identity_residual"]) < 1e-12
    assert np.isclose(a["E_esc"] + a["E_core"], a["E_inj"], rtol=1e-12)


@pytest.mark.parametrize("relativity", [None, "worldline"])
def test_coherent_scattering_in_transport_leaves_w_and_e_cm_and_books_the_work(relativity):
    """sobolev_tla with eps = 0 is pure coherent scattering: under energy
    packets w must not change (nu_rest == nu_abs,cm), the comoving deposit
    is exactly zero, and E_dep_lab is the Doppler work -- the lab energy each
    packet lost, summed over the run, O(beta) of the interacting energy."""
    fa, _ = tfm.three_level(3.0, 1.0, 0.0)
    fa.emis_w = np.array([1.0, 1.0])
    lo, hi = tfm.pump_band()
    r_core, r_out = (tfm.R_CORE, tfm.R_OUT) if relativity is None else (CT / 30.0, CT / 10.0)
    res = run_mc(fa, r_core, r_out, T_EXP, lo, hi, 30000, "sobolev_tla", seed=6, eps=0.0,
                 relativity=relativity, packets="energy")
    assert res["n_interactions"] > 0
    assert np.array_equal(res["w"], res["w_launch"])
    assert res["e_dep_cm"].sum() == 0.0
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12
    # the work equals the lab energy change of every packet that died
    done = res["fate"] > 0
    lab_loss = float(np.sum(res["w_launch"][done] * H * res["nu_launch"][done]
                            - res["w"][done] * H * res["nu_final"][done]))
    assert np.isclose(a["E_dep_lab"], lab_loss, rtol=1e-9, atol=1e-12 * a["E_inj"])
    assert a["W"] == a["E_dep_lab"]
    beta_out = r_out / CT
    assert 0 < abs(a["W"]) / a["E_interacting"] < 3 * beta_out
