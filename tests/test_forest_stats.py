"""Forest statistics: the phase-diagram axes.

These pin behaviour that later experiments will read as physics, so the limits
matter more than the typical values -- a weak forest must give E/S = 1, a
coherent kernel must give zero redistribution range, and crowding must be
invariant to where in wavelength the forest sits.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paper2/phase1"))
from forest_mc import ForestAtom

from sobolev.constants import C
from sobolev.forest_stats import (SE_sums, band_saturation, crowding,
                                  forest_summary, redistribution_range,
                                  saturation_stats, spacing_stats)

T_EXP = 86400.0


def _atom(nu0, tau, t_exp=T_EXP):
    """A ForestAtom with prescribed optical depths, one upper level per line."""
    from sobolev.optical_depth import tau_sobolev
    nu0 = np.asarray(nu0, float); tau = np.asarray(tau, float)
    f = np.full(nu0.size, 0.01)
    n_l = np.array([t / tau_sobolev(fi, 1.0, C / n, t_exp)
                    for t, fi, n in zip(tau, f, nu0)])
    n = nu0.size
    return ForestAtom(nu0=nu0, f_osc=f, n_lower=n_l, n_upper=np.zeros(n),
                      A=np.full(n, 1e8), lower=np.zeros(n, int),
                      upper=np.arange(1, n + 1), t_exp=t_exp,
                      tau_min=1e-6, stim=False)


def test_saturation_census():
    a = _atom([5e14, 6e14, 7e14], [0.05, 0.5, 8.0])
    s = saturation_stats(a)
    assert s["n_opacity"] == 3
    assert s["n_tau_gt1"] == 1
    assert s["n_tau_gt01"] == 2
    assert s["tau_max"] == pytest.approx(8.0, rel=1e-6)


def test_weak_forest_has_E_over_S_one():
    """E = sum(1-e^-tau) -> sum tau as tau -> 0: no saturation deficit."""
    a = _atom([5e14, 6e14, 7e14], [1e-4, 2e-4, 3e-4])
    s = SE_sums(a)
    assert s["E_over_S"] == pytest.approx(1.0, abs=1e-3)
    assert s["deficit"] == pytest.approx(0.0, abs=1e-6)


def test_saturated_forest_has_E_over_S_small():
    a = _atom([5e14, 6e14, 7e14], [20.0, 20.0, 20.0])
    s = SE_sums(a)
    assert s["E_over_S"] == pytest.approx(3.0 / 60.0, rel=1e-3)
    assert s["deficit"] == pytest.approx(57.0, rel=1e-3)


def test_crowding_is_scale_free():
    """The same forest shifted in wavelength has the same crowding."""
    tau = [5.0] * 4
    a = _atom(np.array([4.0e14, 4.4e14, 4.8e14, 5.2e14]), tau)
    b = _atom(np.array([4.0e14, 4.4e14, 4.8e14, 5.2e14]) * 3.0, tau)
    assert crowding(a)["n_sat_per_lnlam"] == pytest.approx(
        crowding(b)["n_sat_per_lnlam"], rel=1e-9)


def test_crowding_counts_only_saturated_lines():
    a = _atom([5e14, 6e14, 7e14, 8e14], [0.01, 0.01, 5.0, 5.0])
    c = crowding(a)
    assert c["n_sat_per_lnlam"] == pytest.approx(c["n_per_lnlam"] / 2, rel=1e-9)


def test_spacing_is_positive_and_ordered():
    a = _atom([5e14, 6e14, 7e14], [1.0, 1.0, 1.0])
    s = spacing_stats(a, v_doppler=1e7)
    assert s["dv_min"] > 0 and s["dv_median"] >= s["dv_min"]
    assert s["overlap_median"] > 0


def test_coherent_redistribution_has_zero_range():
    nu = np.linspace(4e14, 9e14, 5000)
    r = redistribution_range(nu, nu.copy())
    assert r["mean_abs_dlnlam"] == pytest.approx(0.0, abs=1e-15)
    assert r["median_abs_dlnlam"] == pytest.approx(0.0, abs=1e-15)


def test_same_group_fraction_is_one_for_coherent_scattering():
    nu = np.linspace(4e14, 9e14, 5000)
    edges = np.geomspace(3.9e14, 9.1e14, 17)
    r = redistribution_range(nu, nu.copy(), edges=edges)
    assert r["same_group_frac"] == pytest.approx(1.0)
    assert r["mean_abs_dgroup"] == pytest.approx(0.0)


def test_redward_redistribution_has_positive_signed_range():
    """nu_out < nu_in is a redward move: ln(nu_in/nu_out) > 0."""
    rng = np.random.default_rng(0)
    nu_in = rng.uniform(6e14, 9e14, 20000)
    nu_out = nu_in * 0.8
    r = redistribution_range(nu_in, nu_out)
    assert r["mean_dlnlam"] == pytest.approx(-np.log(0.8), rel=1e-6)
    assert r["mean_abs_dlnlam"] == pytest.approx(-np.log(0.8), rel=1e-6)


def test_blueward_redistribution_is_signed_negative():
    rng = np.random.default_rng(1)
    nu_in = rng.uniform(4e14, 6e14, 20000)
    r = redistribution_range(nu_in, nu_in * 1.25)
    assert r["mean_dlnlam"] < 0
    assert r["mean_abs_dlnlam"] > 0


def test_wider_redistribution_crosses_more_groups():
    rng = np.random.default_rng(2)
    nu_in = rng.uniform(5e14, 8e14, 40000)
    edges = np.geomspace(3e14, 1.2e15, 33)
    near = redistribution_range(nu_in, nu_in * 0.99, edges=edges)
    far = redistribution_range(nu_in, nu_in * 0.6, edges=edges)
    assert far["mean_abs_dgroup"] > near["mean_abs_dgroup"]
    assert far["same_group_frac"] < near["same_group_frac"]


def test_forest_summary_is_one_row():
    a = _atom([5e14, 6e14, 7e14], [0.5, 5.0, 50.0])
    nu = np.array([5e14, 6e14, 7e14] * 100, float)
    row = forest_summary(a, v_doppler=1e7, events=(nu, nu * 0.9))
    for k in ("n_opacity", "tau_max", "E_over_S", "n_sat_per_lnlam",
              "dv_median", "mean_abs_dlnlam"):
        assert k in row and np.isfinite(row[k])


def test_empty_events_do_not_crash():
    r = redistribution_range(np.array([np.nan]), np.array([np.nan]))
    assert r["n_events"] == 0


def test_band_saturation_restricts_to_the_band():
    a = _atom([5e14, 6e14, 7e14, 8e14], [5.0, 5.0, 0.01, 0.01])
    b = band_saturation(a, 4.5e14, 6.5e14)
    assert b["n_band"] == 2 and b["n_sat_band"] == 2
    b2 = band_saturation(a, 6.5e14, 8.5e14)
    assert b2["n_band"] == 2 and b2["n_sat_band"] == 0


def test_band_saturation_S_matches_the_lines_inside():
    a = _atom([5e14, 6e14, 7e14], [2.0, 3.0, 100.0])
    b = band_saturation(a, 4.5e14, 6.5e14)
    assert b["S_band"] == pytest.approx(5.0, rel=1e-6)


def test_empty_band_is_reported_not_crashed():
    a = _atom([5e14, 6e14], [1.0, 1.0])
    b = band_saturation(a, 1e15, 2e15)
    assert b["n_band"] == 0 and b["n_sat_band"] == 0
