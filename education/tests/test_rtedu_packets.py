"""Chapter 2: the interaction depth of a packet is exponentially distributed.
The tests state the statistical prediction they check (PI amendment 7)."""
import numpy as np

from rtedu.packets import sample_tau, propagate_slab, walk_2d


def test_tau_sampling_is_exponential():
    """With a fixed seed, the Kolmogorov statistic of 10^5 draws against
    F(tau) = 1 - e^{-tau} must lie below the 1 % critical value
    1.63 / sqrt(N); the sample mean and variance (both 1 for the unit
    exponential) must agree within 4 standard errors."""
    n = 100_000
    tau = sample_tau(np.random.default_rng(2), n)
    assert np.all(np.isfinite(tau)) and tau.min() >= 0.0
    tau_sorted = np.sort(tau)
    F = -np.expm1(-tau_sorted)
    i = np.arange(1, n + 1)
    D = max(np.max(i / n - F), np.max(F - (i - 1) / n))
    assert D < 1.63 / np.sqrt(n), f"KS statistic {D:.4f}"
    assert abs(tau.mean() - 1.0) < 4.0 / np.sqrt(n)            # std error of the mean is 1/sqrt(N)
    assert abs(tau.var() - 1.0) < 4.0 * np.sqrt(8.0 / n)       # var of the sample variance ~ (mu4 - sigma^4)/N = 8/N


def test_slab_transmission_without_scattering_is_a_single_draw():
    """With albedo 0 the transmitted fraction is P(tau_draw > tau_slab) = e^{-tau}
    within 4 binomial sigma, and no transmitted packet interacted."""
    rng = np.random.default_rng(3)
    for tau in (0.3, 1.0, 3.0):
        n = 40_000
        out = propagate_slab(rng, tau, n, albedo=0.0)
        p = np.exp(-tau); sig = np.sqrt(p * (1 - p) / n)
        assert abs(out["transmitted"] - p) < 4 * sig
        assert out["mean_interactions_transmitted"] == 0.0
        assert abs(out["transmitted"] + out["absorbed"] + out["reflected"] - 1.0) < 1e-12


def test_elastic_scattering_lengthens_the_path_and_reflects():
    """Pure elastic scattering (albedo 1) absorbs nothing: every packet is
    transmitted or reflected, transmitted packets have interacted, and their
    path exceeds the slab depth."""
    out = propagate_slab(np.random.default_rng(4), 2.0, 20_000, albedo=1.0)
    assert out["absorbed"] == 0.0 and out["still_alive"] == 0
    assert out["reflected"] > 0.2 and out["transmitted"] > 0.1
    assert out["mean_interactions_transmitted"] > 1.0 and out["mean_path_transmitted"] > 2.0


def test_walk_2d_ends_outside_the_disc():
    p = walk_2d(np.random.default_rng(5), radius=4.0)
    assert np.hypot(*p[-1]) >= 4.0 and np.all(np.hypot(p[:-1, 0], p[:-1, 1]) < 4.0)
