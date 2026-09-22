"""Chapter 4: the comoving frequency sweep, the resonance radius, the toy
Sobolev depth, and the two quantities that both carry the name Sobolev."""
import numpy as np

from rtedu import C, SIGMA_CLASSICAL
from rtedu.sobolev import (comoving_frequency, resonance_position, tau_sobolev_toy, interaction_probability,
                           escape_probability, resonance_crossings, first_interaction, toy_line_list)


def test_comoving_frequency_falls_outward_and_resonates_where_it_should():
    t = 2.0 * 86400.0; nu = 5e14
    r = np.linspace(0, 0.1 * C * t, 5)
    nc = comoving_frequency(nu, r, t)
    assert np.all(np.diff(nc) < 0) and np.isclose(nc[-1], 0.9 * nu)
    nu_line = 0.95 * nu
    r_res = resonance_position(nu, nu_line, t)
    assert np.isclose(comoving_frequency(nu, r_res, t), nu_line)
    assert np.isnan(resonance_position(nu, 1.05 * nu, t))          # a bluer line is never reached


def test_toy_sobolev_depth_has_the_production_formula_shape():
    """tau_S = (pi e^2 / m_e c) f n_l lambda_0 t: linear in every factor."""
    tau = tau_sobolev_toy(1e6, 0.1, 5e-5, 1e5)
    assert np.isclose(tau, SIGMA_CLASSICAL * 0.1 * 1e6 * 5e-5 * 1e5)
    assert np.isclose(tau_sobolev_toy(2e6, 0.1, 5e-5, 1e5), 2 * tau)


def test_sobolev_interaction_probability():
    """P_int = 1 - e^{-tau}: direct counting of the coin at N = 10^5 agrees
    within 4 binomial sigma for tau = 0.1, 1, 3."""
    rng = np.random.default_rng(8); n = 100_000
    for tau in (0.1, 1.0, 3.0):
        p = float(interaction_probability(tau))
        assert np.isclose(p, 1 - np.exp(-tau))
        hits = np.mean(rng.random(n) < p)
        assert abs(hits - p) < 4 * np.sqrt(p * (1 - p) / n)


def test_escape_probability_limits():
    """beta = (1 - e^{-tau}) / tau is NOT the interaction probability:
    beta -> 1 as tau -> 0, beta -> 1/tau as tau -> infinity, and
    P_int = tau * beta exactly."""
    assert escape_probability(0.0) == 1.0 and abs(escape_probability(1e-6) - 1.0) < 1e-6
    assert abs(escape_probability(100.0) * 100.0 - 1.0) < 1e-12
    tau = np.array([0.01, 0.5, 2.0, 20.0])
    assert np.allclose(interaction_probability(tau), tau * escape_probability(tau))
    assert np.all(escape_probability(tau[:-1]) > interaction_probability(tau[:-1]) / tau[:-1] - 1e-15)


def test_resonance_crossings_in_order_and_only_reachable_lines():
    t = 3.0 * 86400.0; nu = 6e14; r_out = 0.2 * C * t
    nu_lines = nu * np.array([1.02, 0.99, 0.95, 0.90, 0.75])       # the first is bluer, the last beyond r_out
    tau_lines = np.array([5.0, 1.0, 0.5, 2.0, 1.0])
    cr = resonance_crossings(nu, nu_lines, tau_lines, 0.0, r_out, t, rng=None)
    assert [c["line"] for c in cr] == [1, 2, 3]
    assert np.all(np.diff([c["r_res"] for c in cr]) > 0) and not any(c["interacted"] for c in cr)
    assert np.allclose([c["p_int"] for c in cr], 1 - np.exp(-tau_lines[[1, 2, 3]]))
    # with a coin the walk stops at the first interaction; over many packets
    # the probability of escaping all three is prod(e^{-tau})
    rng = np.random.default_rng(9); n = 20_000
    esc = sum(first_interaction(nu, nu_lines, tau_lines, 0.0, r_out, t, rng) is None for _ in range(n)) / n
    p = float(np.exp(-tau_lines[[1, 2, 3]].sum()))
    assert abs(esc - p) < 4 * np.sqrt(p * (1 - p) / n)


def test_toy_line_list_is_teachable():
    """Depths between 0.01 and 2 (both regimes of beta), eight of twelve lines
    within reach of the edge, and an escape probability a 2x10^4-packet run
    can measure."""
    L = toy_line_list()
    assert L["tau"].min() > 0.005 and L["tau"].max() < 2.0 and (L["tau"] > 1).sum() >= 2
    reach = L["frac"] > 1.0 - L["v_max_c"]
    assert reach.sum() == 8 and 0.02 < np.exp(-L["tau"][reach].sum()) < 0.2
