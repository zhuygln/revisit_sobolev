"""Chapter 6: the toy atom's branching and cascades; chapter 5's epsilon model."""
import numpy as np

from rtedu.atom import five_level_atom, H_ERG_S
from rtedu.redistribution import EpsilonRedistribution


def test_branching_probabilities_sum_to_one():
    a = five_level_atom()
    for level in range(1, a.n_levels):
        k, p = a.branching(level)
        assert k.size > 0 and abs(p.sum() - 1.0) < 1e-12 and np.all(p > 0)
        assert np.allclose(p, a.A[k] / a.A[k].sum())
    k0, p0 = a.branching(0)
    assert k0.size == 0 and p0.size == 0


def test_cascade_reaches_the_ground_and_conserves_energy():
    """Every cascade from level 4 ends at the ground and the photon
    energies sum to h nu(4 -> 0) exactly: fluorescence conserves energy by
    construction, whatever route it takes."""
    a = five_level_atom(); rng = np.random.default_rng(1)
    E4 = H_ERG_S * a.nu[0]                     # line 0 is 4 -> 0
    for _ in range(200):
        c = a.cascade(rng, 4)
        assert a.lower[c[-1]] == 0 and np.all(a.upper[c[1:]] == a.lower[c[:-1]])
        assert abs(a.photon_energies(c).sum() - E4) < 1e-12 * E4


def test_populations_and_line_list_scale_as_they_should():
    a = five_level_atom()
    n = a.populations(4000.0, 100.0)
    assert abs(n.sum() - 100.0) < 1e-9 and np.all(np.diff(n) < 0)
    tau1 = a.line_list(4000.0, 100.0, 1e5); tau2 = a.line_list(4000.0, 200.0, 2e5)
    assert np.allclose(tau2, 4 * tau1)
    hot = a.populations(20000.0, 100.0)
    assert hot[4] > n[4]                        # a hotter atom populates its upper levels more


def test_epsilon_model_limits():
    """eps = 0 always re-emits in the absorbing line; eps = 1 draws from the
    thermal emissivity; in between the thermal fraction is eps within the
    binomial noise."""
    a = five_level_atom(); emis = a.thermal_emissivity(4000.0, 100.0)
    rng = np.random.default_rng(2)
    r0 = EpsilonRedistribution(0.0, emis); assert all(r0(rng, 3) == 3 for _ in range(100))
    r1 = EpsilonRedistribution(1.0, emis); draws = np.array([r1(rng, 3) for _ in range(20_000)])
    freq = np.bincount(draws, minlength=a.n_lines) / draws.size
    assert np.all(np.abs(freq - emis) < 4 * np.sqrt(emis * (1 - emis) / draws.size) + 1e-9)
    r5 = EpsilonRedistribution(0.5, emis); same = np.mean([r5(rng, 3) == 3 for _ in range(20_000)])
    assert abs(same - (0.5 + 0.5 * emis[3])) < 4 * np.sqrt(0.25 / 20_000)
