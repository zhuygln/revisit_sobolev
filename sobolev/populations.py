"""Boltzmann level populations (Stage B of babystep_plan.md section 16).

    n_l / n_ion = g_l exp(-E_l / kT) / Z(T),   Z(T) = sum_i g_i exp(-E_i / kT)

Ionization balance (Saha, Stage C) is deliberately NOT here yet: fractions are
relative to the total population of one ionization stage, so uncertainty in the
transfer comparison cannot hide in the ionization state.

Energies are taken in cm^-1, matching the GSI level files; statistical weights
are g = 2J + 1.
"""

import numpy as np

from .constants import C, H, K_B

# hc in erg cm: converts an energy in cm^-1 to erg.
HC = H * C


def statistical_weight(j):
    """g = 2J + 1."""
    return 2.0 * np.asarray(j, dtype=float) + 1.0


def partition_function(g, energy_cm, temperature):
    """Z(T) = sum_i g_i exp(-E_i hc / kT) over the supplied level list.

    The sum runs over whatever levels are passed in; a truncated level list
    truncates Z. For GSI files all levels below the ionization threshold are
    included, which is the standard choice.
    """
    g = np.asarray(g, dtype=float)
    e_erg = HC * np.asarray(energy_cm, dtype=float)
    return np.sum(g * np.exp(-e_erg / (K_B * temperature)))


def boltzmann_fractions(g, energy_cm, temperature):
    """n_i / n_ion for every level in the list, normalized by Z(T) of that list."""
    g = np.asarray(g, dtype=float)
    e_erg = HC * np.asarray(energy_cm, dtype=float)
    weights = g * np.exp(-e_erg / (K_B * temperature))
    return weights / np.sum(weights)


def boltzmann_fractions_from_levels(levels_df, temperature):
    """Level fractions for a GSI levels DataFrame (as returned by load_gsi).

    Uses the J and Energy columns; returns an array aligned with the DataFrame
    rows (which the GSI files sort by energy).
    """
    return boltzmann_fractions(
        statistical_weight(levels_df["J"].to_numpy()),
        levels_df["Energy"].to_numpy(),
        temperature,
    )
