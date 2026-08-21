"""Impact-parameter quadrature shared by the analytic and resolved legs.

Why this exists. The attenuation of a ray at impact parameter p is a STEP
function of p for every Sobolev resonance -- the ray crosses line k iff
p > sqrt(r_core^2 - z_k^2) -- so any p-quadrature carries an O(1/n) step
error whose sign depends on where each step falls between the nodes. Two legs
sampled with different nodes therefore carry different O(1/n) errors, and
their DIFFERENCE inherits the mismatch rather than the physics. That is how
Delta_Sobolev changed sign between 12 and 96 solver rays in the high-beta
pilot while the analytic leg sat on 200 rays throughout.

The cure is not a cleverer rule but the SAME rule in both legs: with matched
nodes the step error is common and cancels in the difference, and what
remains is O(v_D / Delta v_shell / n) from the resolved profile smoothing the
step. Midpoint is the right choice for a step integrand -- Gauss-Legendre's
high order buys nothing against a discontinuity and its clustered nodes make
the step error erratic -- and it is also what the existing legs already use,
so the default paths stay bit-identical.

`core_rays_midpoint` reproduces, to the last ulp, the expression that
`sobolev_attenuation` has used since the beginning; tests pin that.
"""

from dataclasses import dataclass

import numpy as np


def core_rays_midpoint(r_core, n_p):
    """The legacy core-disk midpoint rule, bit-for-bit.

    This is exactly the expression that lived inline in `sobolev_attenuation`;
    it is kept as its own function so the identity can be asserted rather
    than assumed.
    """
    return np.linspace(0.0, r_core, n_p, endpoint=False) + r_core / (2 * n_p)


@dataclass(frozen=True)
class RaySet:
    """Impact parameters `p`, quadrature weights `w` (already including the
    p dp measure), and a mask `is_core` for rays that start on the core."""

    p: np.ndarray
    w: np.ndarray
    is_core: np.ndarray

    @property
    def n_core(self):
        return int(self.is_core.sum())

    def average(self, values):
        """Core-disk average of `values` (shape (n_core, ...) or (n_p, ...)),
        i.e. sum(w v)/sum(w) over the core rays. This is the quantity both
        attenuation legs return."""
        v = np.asarray(values)
        if v.shape[0] == self.p.size:
            v = v[self.is_core]
        w = self.w[self.is_core]
        return np.tensordot(w, v, axes=(0, 0)) / w.sum()

    @staticmethod
    def midpoint(r_core, r_out, n_core, n_env=0):
        """Midpoint nodes on [0, r_core] (core) and [r_core, r_out] (envelope).

        The core nodes are `core_rays_midpoint(r_core, n_core)`; the weights
        are p dp, which for the core integrate to r_core^2/2 exactly.
        """
        p_core = core_rays_midpoint(r_core, n_core)
        dp_core = r_core / n_core
        if n_env > 0:
            p_env = np.linspace(r_core, r_out, n_env, endpoint=False)
            p_env = p_env + 0.5 * (p_env[1] - p_env[0] if n_env > 1 else r_out - r_core)
            dp_env = (r_out - r_core) / n_env
        else:
            p_env = np.empty(0)
            dp_env = 0.0
        p = np.concatenate([p_core, p_env])
        w = np.concatenate([p_core * dp_core, p_env * dp_env])
        is_core = np.concatenate([np.ones(n_core, bool), np.zeros(p_env.size, bool)])
        return RaySet(p, w, is_core)

    @staticmethod
    def gauss_legendre(r_core, r_out, n_core, n_env=0):
        """Gauss-Legendre nodes on the same two intervals. Provided for the
        convergence comparison; NOT recommended for matched legs (see the
        module docstring)."""
        x, wx = np.polynomial.legendre.leggauss(n_core)
        p_core = 0.5 * r_core * (x + 1.0)
        w_core = 0.5 * r_core * wx * p_core
        if n_env > 0:
            x, wx = np.polynomial.legendre.leggauss(n_env)
            p_env = r_core + 0.5 * (r_out - r_core) * (x + 1.0)
            w_env = 0.5 * (r_out - r_core) * wx * p_env
        else:
            p_env, w_env = np.empty(0), np.empty(0)
        p = np.concatenate([p_core, p_env])
        w = np.concatenate([w_core, w_env])
        is_core = np.concatenate([np.ones(n_core, bool), np.zeros(p_env.size, bool)])
        return RaySet(p, w, is_core)
