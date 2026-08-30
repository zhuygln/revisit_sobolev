"""The P9 composition rule: RedistributionKernel.mix().

The plan's required tests for a mixture operator are the limit X_s -> 0 (the
mixture must collapse onto the surviving species) and energy conservation of
the mixed rows. Both are exact statements, not tolerances, because convex
mixing preserves the conservation object identically.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paper3"))
from redistribution import RedistributionKernel


def _kernel(seed, n_groups=8, lo=1e14, hi=1e15, n=20000):
    """A kernel from synthetic events with a species-specific exit structure."""
    rng = np.random.default_rng(seed)
    nu_in = np.exp(rng.uniform(np.log(lo * 1.01), np.log(hi * 0.99), n))
    # each "species" exits onto its own discrete line set
    lines = np.exp(np.linspace(np.log(lo * 1.02), np.log(hi * 0.98), 30 + seed))
    nu_out = rng.choice(lines, size=n)
    return RedistributionKernel.from_branching_mc(
        nu_in, nu_out, np.ones(n), n_groups, nu_lo=lo, nu_hi=hi)


def test_single_species_mix_is_the_identity():
    k = _kernel(1)
    w = np.ones((k.n_groups, 1))
    m = RedistributionKernel.mix([k], w)
    assert np.allclose(m.R, k.R)
    assert np.allclose(m.N_cum, k.N_cum)
    assert np.allclose(m.q_dep, k.q_dep)
    assert np.array_equal(m.disc_off, k.disc_off)
    assert np.allclose(m.disc_vals, k.disc_vals)


def test_mixture_limit_recovers_the_surviving_species():
    """X_b -> 0: the rule must reduce to species a exactly, not approximately."""
    a, b = _kernel(1), _kernel(2)
    w = np.zeros((a.n_groups, 2)); w[:, 0] = 1.0
    m = RedistributionKernel.mix([a, b], w)
    assert np.allclose(m.R, a.R)
    assert np.allclose(m.N_cum, a.N_cum)
    assert np.allclose(m.q_dep, a.q_dep)
    assert np.allclose(m.disc_vals, a.disc_vals)
    assert np.array_equal(m.disc_off, a.disc_off)


def test_mixed_rows_conserve_energy_exactly():
    a, b = _kernel(1), _kernel(2)
    rng = np.random.default_rng(0)
    x = rng.uniform(size=a.n_groups)
    w = np.stack([x, 1 - x], axis=1)
    m = RedistributionKernel.mix([a, b], w)
    assert m.validate_energy() < 1e-12


def test_mixed_rows_are_between_the_species_rows():
    """A convex combination cannot leave the interval its parents span."""
    a, b = _kernel(1), _kernel(2)
    w = np.full((a.n_groups, 2), 0.5)
    m = RedistributionKernel.mix([a, b], w)
    live = ~a.empty_rows & ~b.empty_rows
    lo = np.minimum(a.R[live], b.R[live]); hi = np.maximum(a.R[live], b.R[live])
    assert np.all(m.R[live] >= lo - 1e-12) and np.all(m.R[live] <= hi + 1e-12)


def test_mixed_kernel_samples_inside_the_grid():
    a, b = _kernel(1), _kernel(2)
    w = np.full((a.n_groups, 2), 0.5)
    m = RedistributionKernel.mix([a, b], w)
    rng = np.random.default_rng(3)
    nu = np.exp(rng.uniform(np.log(1.1e14), np.log(9e14), 5000))
    out = m.sample_nu_out(nu, rng)
    ok = np.isfinite(out)
    assert ok.mean() > 0.9
    assert out[ok].min() >= m.edges[0] and out[ok].max() <= m.edges[-1]
    # exits are drawn from the union of the two species' discrete line sets
    union = np.union1d(a.disc_vals, b.disc_vals)
    assert np.isin(out[ok], union).all()


def test_edges_must_match():
    a = _kernel(1)
    b = _kernel(2, lo=2e14)
    with pytest.raises(ValueError, match="identical edges"):
        RedistributionKernel.mix([a, b], np.full((a.n_groups, 2), 0.5))


def test_weight_shape_is_checked():
    a, b = _kernel(1), _kernel(2)
    with pytest.raises(ValueError, match="n_groups"):
        RedistributionKernel.mix([a, b], np.full((a.n_groups, 3), 1 / 3))


def test_unpopulated_species_rows_do_not_dilute():
    """If species b never populated row i, the mixture's row i is species a's,
    whatever weight b nominally carries -- otherwise a blend would inherit
    b's empty rows as spurious deposition."""
    a, b = _kernel(1), _kernel(2)
    b.counts = b.counts.copy(); b.counts[0] = 0.0
    b.empty_rows = b.counts <= 0
    w = np.full((a.n_groups, 2), 0.5)
    m = RedistributionKernel.mix([a, b], w)
    assert np.allclose(m.R[0], a.R[0])
