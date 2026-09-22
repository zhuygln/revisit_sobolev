"""Chapter 1: Beer-Lambert, analytically and by Monte Carlo."""
import numpy as np

from rtedu.slab import Slab, ThreeGroupSlab


def test_beer_lambert():
    """I(s)/I_0 = e^{-kappa s} at every depth, and the Monte Carlo
    transmission agrees within 4 binomial sigma."""
    slab = Slab(depth=2.0, kappa=1.5)
    s = np.linspace(0, 2, 11)
    assert np.allclose(slab.transmission(s=s), np.exp(-1.5 * s))
    assert np.isclose(slab.optical_depth(), 3.0)
    n = 50_000
    mc = slab.mc_transmission(np.random.default_rng(1), n)
    p = np.exp(-3.0); sig = np.sqrt(p * (1 - p) / n)
    assert abs(mc["transmitted"] - p) < 4 * sig


def test_frequency_dependent_opacity_is_evaluated_per_frequency():
    slab = Slab(depth=1.0, kappa=lambda nu: 2.0 / nu)
    assert np.isclose(slab.transmission(nu=2.0), np.exp(-1.0)) and np.isclose(slab.transmission(nu=4.0), np.exp(-0.5))


def test_three_groups_order_by_opacity_and_keep_their_frequency():
    """kappa_B > kappa_V > kappa_IR: transmission increases B -> V -> IR;
    with elastic scattering the bluer group interacts more and walks further;
    every packet ends in the group it started in (there is no
    redistribution before chapter 5), which the per-group bookkeeping makes
    true by construction."""
    tg = ThreeGroupSlab(depth=1.0, kappa_B=3.0, kappa_V=1.0, kappa_IR=0.3)
    res = tg.run(np.random.default_rng(6), 20_000, albedo=0.0)
    assert res["B"]["transmitted"] < res["V"]["transmitted"] < res["IR"]["transmitted"]
    for g in tg.GROUPS:
        p = res[g]["transmission_analytic"]; sig = np.sqrt(p * (1 - p) / 20_000)
        assert abs(res[g]["transmitted"] - p) < 4 * sig
    sc = tg.run(np.random.default_rng(7), 20_000, albedo=1.0)
    assert sc["B"]["mean_interactions_transmitted"] > sc["V"]["mean_interactions_transmitted"] > sc["IR"]["mean_interactions_transmitted"]
    assert sc["B"]["mean_path_transmitted"] > sc["IR"]["mean_path_transmitted"]
