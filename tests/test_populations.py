"""Boltzmann populations: analytic two-level checks and limits.

Wrong populations shift every tau by the same factor in both treatments, so
these bugs cancel in E_Sob but corrupt absolute optical depths -- same failure
mode as a wrong prefactor, and equally invisible without direct tests.
"""

from pathlib import Path

import numpy as np

from sobolev.atomic_data import load_gsi
from sobolev.constants import C, H, K_B
from sobolev.populations import (
    boltzmann_fractions,
    boltzmann_fractions_from_levels,
    partition_function,
    statistical_weight,
)

LEVELS_EXCERPT = Path(__file__).parent / "data" / "57LaII_levels_calib_excerpt.txt"


def test_two_level_analytic():
    # g = (1, 3), E = (0, 1000 cm^-1): fractions follow the closed form.
    g = np.array([1.0, 3.0])
    e = np.array([0.0, 1000.0])
    t = 5000.0
    boltz = 3.0 * np.exp(-H * C * 1000.0 / (K_B * t))
    expected = np.array([1.0, boltz]) / (1.0 + boltz)
    assert np.allclose(boltzmann_fractions(g, e, t), expected, rtol=1e-12)
    assert np.isclose(partition_function(g, e, t), 1.0 + boltz, rtol=1e-12)


def test_low_temperature_limit_is_ground_state():
    g = np.array([5.0, 7.0, 9.0])
    e = np.array([0.0, 1000.0, 2000.0])
    frac = boltzmann_fractions(g, e, 10.0)
    assert np.isclose(frac[0], 1.0)
    assert frac[1:].max() < 1e-10


def test_high_temperature_limit_is_statistical_weights():
    g = np.array([5.0, 7.0, 9.0])
    e = np.array([0.0, 1000.0, 2000.0])
    frac = boltzmann_fractions(g, e, 1.0e9)
    assert np.allclose(frac, g / g.sum(), rtol=1e-4)


def test_la_ii_levels_from_file():
    levels = load_gsi(LEVELS_EXCERPT)
    assert levels.shape == (20, 10)
    frac = boltzmann_fractions_from_levels(levels, 5000.0)
    assert np.isclose(frac.sum(), 1.0, rtol=1e-12)
    # The rigorous Boltzmann invariant: fractions PER STATISTICAL WEIGHT
    # decrease monotonically with energy (the file is energy-sorted). Raw
    # fractions need not: at 5000 K the J=3 level at 1016 cm^-1 out-populates
    # the J=2 ground state because g=7 beats the mild Boltzmann factor.
    g = statistical_weight(levels["J"].to_numpy())
    assert np.all(np.diff(frac / g) < 0)
    assert g[0] == 5.0
