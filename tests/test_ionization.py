"""Paper IV Phase 7: the Saha module -- analytic limits, charge neutrality,
and the declared policies."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sobolev import ionization as ion                # noqa: E402
from sobolev.constants import K_B                    # noqa: E402


def test_saha_ratio_is_the_textbook_formula():
    T, n_e = 5000.0, 1e9
    r = ion.saha_ratio(T, n_e, 5.5, 1.0, 1.0)
    expect = 2.0 * ion.SAHA_CONST * T ** 1.5 * np.exp(-5.5 * ion.EV / (K_B * T)) / n_e
    assert np.isclose(r, expect, rtol=1e-12)
    # hydrogen-like sanity: 13.6 eV at 1e4 K, n_e = 1e14: ratio of order 1
    assert 0.01 < ion.saha_ratio(1e4, 1e14, 13.598, 1.0, 2.0) < 100


def test_stage_fractions_sum_to_one_and_have_the_right_limits():
    f = ion.stage_fractions(4000.0, 1e8, (5.5, 11.0, 22.0), (1.0, 1.0, 1.0))
    assert f.shape == (4,) and np.isclose(f.sum(), 1.0)
    cold = ion.stage_fractions(1500.0, 1e12, (5.5, 11.0, 22.0), (1.0, 1.0, 1.0))
    hot = ion.stage_fractions(30000.0, 1e6, (5.5, 11.0, 22.0), (1.0, 1.0, 1.0))
    assert cold[0] > 0.99 and hot[3] > 0.99
    # monotone: raising T never lowers the mean charge
    zs = [np.sum(np.arange(4) * ion.stage_fractions(T, 1e8, (5.5, 11.0, 22.0), (1.0, 1.0, 1.0))) for T in (2000, 3000, 5000, 8000, 12000)]
    assert np.all(np.diff(zs) >= 0)


def test_lanthanides_are_singly_ionised_at_the_benchmark_conditions():
    """At 3400 K and rho = 3e-15 the plan's II-only assumption is a good one;
    at 2300 K the neutral stage appears -- the 'scale' policy matters late."""
    X = {"Ce": 7e-3, "Nd": 1.4e-2, "bulk": 1.0 - 2.1e-2}
    part = {"Ce": (30.0, 20.0), "Nd": (60.0, 40.0)}
    n_e, fr = ion.solve_ionization(3400.0, 3.3e-15, X, part)
    assert 1e5 < n_e < 1e12
    assert fr["Ce"][1] > 0.95 and fr["Nd"][1] > 0.95 and fr["Ce"][2] < 0.05
    n_e2, fr2 = ion.solve_ionization(2300.0, 7e-16, X, part)
    assert fr2["Ce"][0] > 0.005                        # a neutral fraction appears (0.97 %)
    n_e3, fr3 = ion.solve_ionization(2300.0, 7e-16, X, part, z1_policy="drop")
    assert fr3["Ce"][0] == 0.0 and abs(sum(fr3["Ce"]) - 1) < 1e-12
    # hotter: III takes over
    _, fr4 = ion.solve_ionization(9000.0, 3.3e-15, X, part)
    assert fr4["Ce"][2] > fr4["Ce"][1]


def test_charge_neutrality_holds_for_the_solution():
    X = {"La": 4e-3, "bulk": 0.996}
    part = {"La": (10.0, 5.0)}
    T, rho = 4000.0, 1e-14
    n_e, fr = ion.solve_ionization(T, rho, X, part)
    n_la = rho * 4e-3 / (138.905 * ion.M_U); n_b = rho * 0.996 / (ion.BulkSpecies().A * ion.M_U)
    fb = ion.stage_fractions(T, n_e, ion.BulkSpecies().chi_ev, ion.BulkSpecies().u_ratio)
    q = n_la * np.sum(np.arange(4) * np.array(fr["La"])) + n_b * np.sum(np.arange(4) * fb)
    assert abs(q - n_e) / n_e < 1e-9
    assert set(ion.CHI_EV) >= {"La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"}
    with pytest.raises(ValueError):
        ion.u_ratios_for("Ce", 4000.0, 1.0, 1.0, z1_policy="?")
