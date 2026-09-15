"""Redistribution models: what comes out of a line interaction.

Chapter 5: the scalar epsilon (two-level-atom) model,
    R^(eps) = (1 - eps) R_coherent + eps R_thermal,
coherent = re-emit in the same line, thermal = draw a line from the atom's
thermal emissivity. The matrix models of Part III are added in E3.
"""
import numpy as np


class EpsilonRedistribution:
    """After an interaction in line `k`: with probability 1 - eps the packet
    is re-emitted in the same line (coherent scattering); with probability
    eps it thermalises and is re-emitted in a line drawn from the thermal
    emissivity `emis` (energy-weighted, normalised)."""

    def __init__(self, eps, emis):
        self.eps = float(eps)
        self.emis = np.asarray(emis, float) / np.sum(emis)
        self.cum = np.cumsum(self.emis)

    def __call__(self, rng, k):
        if rng.random() < self.eps:
            return int(np.searchsorted(self.cum, rng.random()))
        return int(k)

    @property
    def n_parameters(self):
        return 1
