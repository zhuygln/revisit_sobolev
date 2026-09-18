"""A slab with an opacity that may depend on frequency (chapters 1 and 3).

Beer-Lambert: dI/ds = -chi I  =>  I(s) = I_0 e^{-chi s} = I_0 e^{-tau}.
The Monte Carlo transmission of `packets.propagate_slab` reproduces it.
"""
import numpy as np

from .packets import propagate_slab


class Slab:
    """A plane-parallel slab of geometric depth `depth` and opacity `kappa`,
    a number or a callable of frequency (per unit length)."""

    def __init__(self, depth, kappa):
        self.depth = float(depth)
        self._kappa = kappa

    def kappa(self, nu=None):
        return float(self._kappa(nu)) if callable(self._kappa) else float(self._kappa)

    def optical_depth(self, nu=None, s=None):
        """tau along the normal to depth s (the whole slab by default)."""
        s = self.depth if s is None else s
        return self.kappa(nu) * np.asarray(s, float)

    def transmission(self, nu=None, s=None):
        """The analytic Beer-Lambert transmission e^{-tau}."""
        return np.exp(-self.optical_depth(nu, s))

    def mc_transmission(self, rng, n, nu=None, albedo=0.0):
        """The Monte Carlo transmitted fraction of n packets (dict of
        packets.propagate_slab); equals e^{-tau} within the binomial noise
        when albedo = 0."""
        return propagate_slab(rng, self.optical_depth(nu), n, albedo=albedo)


class ThreeGroupSlab:
    """Chapter 3: three frequency groups B, V, IR through one slab, with
    kappa_B > kappa_V > kappa_IR. Frequency is a label: a packet keeps its
    group for ever (no redistribution before chapter 5); the groups differ
    only in how far a packet travels before it interacts."""

    GROUPS = ("B", "V", "IR")

    def __init__(self, depth, kappa_B, kappa_V, kappa_IR):
        assert kappa_B > kappa_V > kappa_IR > 0.0
        self.depth = float(depth)
        self.kappa = dict(B=float(kappa_B), V=float(kappa_V), IR=float(kappa_IR))

    def optical_depth(self, group):
        return self.kappa[group] * self.depth

    def transmission(self, group):
        return float(np.exp(-self.optical_depth(group)))

    def run(self, rng, n_per_group, albedo=0.0):
        """Equal numbers of packets per group; returns per-group dicts with
        tau, the analytic transmission, and the Monte Carlo result."""
        out = {}
        for g in self.GROUPS:
            mc = propagate_slab(rng, self.optical_depth(g), n_per_group, albedo=albedo)
            out[g] = dict(tau=self.optical_depth(g), transmission_analytic=self.transmission(g), **mc)
        return out
