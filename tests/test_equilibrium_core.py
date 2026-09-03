"""The two harness additions of Paper III §4.39: an inner boundary in
radiative equilibrium (a normalization, not a transport change) and the
thermalization of re-absorption chains that never escape their line.

Both must be inert by default -- every earlier result used the absorbing core
and `chain_overflow="raise"` -- and both must keep the energy identity closed.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paper2/phase1"))
from forest_mc import run_mc                      # noqa: E402
from sobolev import photometry as phot             # noqa: E402
from sobolev.constants import H                    # noqa: E402

sys.path.insert(0, str(ROOT / "tests"))
from test_forest_mc import three_level, pump_band, R_CORE, R_OUT, T_EXP   # noqa: E402


def _res(tau=3.0, a31=1.0, a32=0.0, seed=3, **kw):
    fa, _ = three_level(tau, a31, a32)
    lo, hi = pump_band()
    return fa, lo, hi, run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, 20000, "sobolev_branch",
                              seed=seed, **kw)


def test_equilibrium_core_recycles_the_returned_energy():
    """Pure resonant scattering off a tau = 3 line sends a good fraction of the
    packets back to the core. The absorbing core loses that energy; the
    equilibrium core returns it, so the emergent luminosity is the injected
    window luminosity less the adiabatic deposition -- for the same packets."""
    fa, lo, hi, res = _res()
    a = res["accounting"]
    l_core = 1.0e40
    f_ret = phot.return_fraction(res)
    assert 0.02 < f_ret < 0.8, f_ret
    l_abs = phot.bolometric(res, l_core)
    l_eq = phot.bolometric(res, l_core, core="equilibrium")
    assert l_abs == pytest.approx(l_core * a["E_esc"] / a["E_inj"])
    assert l_eq == pytest.approx(l_core * a["E_esc"] / (a["E_esc"] + a["E_dep_lab"]), rel=1e-12)
    assert l_eq > l_abs
    assert l_eq == pytest.approx(l_core, rel=0.05)        # deposition is O(v/c)
    edges = phot.nu_edges(1000.0, 30000.0, 200)
    for core in ("absorbing", "equilibrium"):
        lnu = phot.emergent_lnu(res, edges, l_core, core)
        assert np.sum(lnu * np.diff(edges)) == pytest.approx(phot.bolometric(res, l_core, core), rel=1e-12)
    # the spectral SHAPE is the same object under both normalizations
    r = phot.emergent_lnu(res, edges, l_core, "equilibrium") / np.maximum(phot.emergent_lnu(res, edges, l_core), 1e-300)
    live = phot.emergent_lnu(res, edges, l_core) > 0
    assert np.allclose(r[live], r[live][0])


def test_equilibrium_core_equals_absorbing_when_nothing_returns():
    fa, lo, hi, res = _res(tau=1e-4)
    assert phot.return_fraction(res) == pytest.approx(0.0, abs=2e-3)
    l_core = 3.0e39
    assert phot.bolometric(res, l_core, "equilibrium") == pytest.approx(phot.bolometric(res, l_core), rel=5e-3)


def test_unknown_core_rejected():
    fa, lo, hi, res = _res()
    with pytest.raises(ValueError):
        phot.bolometric(res, 1.0, core="mirror")


def test_chain_overflow_default_raises_and_absorb_thermalizes():
    """A tau = 400 resonance line (beta = 1/400) chains ~400 draws per
    interaction; with chain_max = 20 the default aborts and "absorb" books
    the trapped packets as fate 3 with the identity still closed to roundoff.
    """
    with pytest.raises(RuntimeError, match="chain"):
        _res(tau=400.0, chain_max=20)
    fa, lo, hi, res = _res(tau=400.0, chain_max=20, chain_overflow="absorb")
    a = res["accounting"]
    assert res["n_trapped"] > 0
    assert res["n_absorbed"] == res["n_trapped"]
    assert a["E_abs"] > 0
    assert abs(a["identity_residual"]) < 1e-12
    assert res["n_escaped"] + res["n_core"] + res["n_absorbed"] == res["n_packets"]
    # the thermalized packets' energy is exactly what the equilibrium core recycles
    l_core = 1.0e40
    assert phot.bolometric(res, l_core, "equilibrium") == pytest.approx(
        l_core * a["E_esc"] / (a["E_esc"] + a["E_dep_lab"]), rel=1e-12)


def test_chain_overflow_is_inert_when_chains_terminate():
    """With a generous chain_max nothing is trapped and the run is bit-identical
    to the default."""
    fa, lo, hi, r0 = _res(tau=3.0)
    fa, lo, hi, r1 = _res(tau=3.0, chain_overflow="absorb", chain_max=10000)
    assert r1["n_trapped"] == 0
    assert np.array_equal(r0["nu_out"], r1["nu_out"])
    assert np.array_equal(r0["fate"], r1["fate"])


def test_conserving_core_emits_the_window_luminosity_for_every_leg():
    """Under "conserving" the escaped energy IS l_core: a grey rescaling of the
    escaped spectrum, so colours are untouched and only leakage out of the
    histogram's range can change the in-window luminosity."""
    fa, lo, hi, res = _res(tau=3.0)
    l_core = 2.0e40
    assert phot.bolometric(res, l_core, "conserving") == pytest.approx(l_core, rel=1e-12)
    edges = phot.nu_edges(1000.0, 30000.0, 200)
    lnu_c = phot.emergent_lnu(res, edges, l_core, "conserving")
    lnu_a = phot.emergent_lnu(res, edges, l_core, "absorbing")
    assert np.sum(lnu_c * np.diff(edges)) == pytest.approx(l_core, rel=1e-12)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    live = lnu_a > 0
    ma = phot.magnitudes(nu_c, np.where(live, lnu_a, 1e-300), phot.BANDS_PHOT)
    mc = phot.magnitudes(nu_c, np.where(live, lnu_c, 1e-300), phot.BANDS_PHOT)
    ca, cc = phot.colors(ma), phot.colors(mc)
    for k in ca:
        if np.isfinite(ca[k]):
            assert cc[k] == pytest.approx(ca[k], abs=1e-9)
    assert phot.deposited_fraction(res) + phot.return_fraction(res) + res["accounting"]["E_esc"] / res["accounting"]["E_inj"] == pytest.approx(1.0, abs=1e-12)


def test_wall_limit_raises_and_is_inert_when_unset():
    with pytest.raises(RuntimeError, match="wall"):
        _res(tau=3.0, wall_s=0.0)
    fa, lo, hi, r0 = _res(tau=3.0)
    fa, lo, hi, r1 = _res(tau=3.0, wall_s=3600.0)
    assert np.array_equal(r0["nu_out"], r1["nu_out"])
