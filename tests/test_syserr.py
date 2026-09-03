"""`syserr.py`: the masked one-mode statistic, its nulls, and the allowance numbers."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
P12 = ROOT / "paper3" / "phase12_grid"
sys.path.insert(0, str(P12))
import syserr as E


def _rank1(m=27, n=38, seed=1):
    rng = np.random.default_rng(seed)
    return np.outer(rng.normal(1, 0.3, m), rng.normal(0, 1, n))


def test_rank1_masked_recovers_rank1_under_random_mask():
    X = _rank1()
    rng = np.random.default_rng(2)
    X[rng.random(X.shape) < 0.4] = np.nan
    f1, u, v = E.rank1_masked(X)
    assert f1 > 1 - 1e-8
    M = np.isfinite(X)
    assert np.allclose(np.outer(u, v)[M], X[M], atol=1e-6)
    assert abs(np.linalg.norm(v) - 1) < 1e-12


def test_sign_scrambled_rank1_is_far_below_one():
    X = _rank1()
    null = E.sign_scramble_null(X, n_draws=50, seed=0)
    assert null["p95"] < 0.5
    assert null["median"] > 0.0


def test_iid_noise_matches_marchenko_pastur_scale():
    rng = np.random.default_rng(3)
    m, n = 27, 17
    X = rng.normal(size=(m, n))
    f1, _, _ = E.rank1_masked(X)
    mp = E.mp_scale(m, n)
    assert 0.5 * mp < f1 < 2 * mp
    f_svd, ncol = E.svd_filled(X, min_n=1)
    assert ncol == n and abs(f_svd - f1) < 1e-9   # no gaps: the masked fit is the SVD


def test_threshold_fractions_on_hand_vector():
    X = np.array([[0.2, -0.6, 1.5, np.nan], [np.nan, 0.9, -1.2, 0.1]])
    f = E.threshold_fractions(X)
    assert f["0.5"] == [4, 6] and f["1"] == [2, 6]


def test_sign_pattern_counts_coepochal_pairs_only():
    keys = [("g", 1.0), ("K", 1.0), ("g", 2.0), ("K", 3.0)]
    X = np.array([[-1.0, 2.0, -1.0, 2.0],      # g<0, K>0 at 1 d; no K at 2 d, no g at 3 d
                  [1.0, 2.0, -1.0, -1.0],      # g>0 at 1 d
                  [-1.0, np.nan, -1.0, 2.0]])  # K masked at 1 d
    assert E.sign_pattern(X, keys) == [1, 2]


def test_chi2_equiv_scales_as_inverse_sigma_squared():
    X = np.array([[1.0, -2.0, np.nan], [0.5, 0.5, 0.5]])
    c1, c5 = E.chi2_equiv(X, 1.0), E.chi2_equiv(X, 0.5)
    assert np.allclose(np.array(c5), 4 * np.array(c1))
    assert np.isclose(c1[0], 2.5) and np.isclose(c1[1], 0.25)


def test_matrix_places_entries_by_key():
    sens = {"points": {"(0.01, 0.1, 0.01)": {"legs": {"L": {"status": "ok", "keys": [["g", 1.0], ["K", 1.0]],
                                                              "d_rt": [-0.5, 1.0], "fit": [-0.4, 0.2]}}},
                       "(0.003, 0.1, 0.01)": {"legs": {"L": {"status": "ok", "keys": [["g", 1.0]],
                                                               "d_rt": [-0.7], "fit": [0.0]}}}}}
    X, pts, keys = E.matrix(sens, "L")
    assert pts[0] == "(0.003, 0.1, 0.01)" and keys == [("g", 1.0), ("K", 1.0)]
    assert X[1, 0] == -0.5 and X[1, 1] == 1.0 and X[0, 0] == -0.7 and np.isnan(X[0, 1])
    R, _, _ = E.matrix(sens, "L", field="residual")
    assert np.isclose(R[1, 1], 0.8)


@pytest.mark.skipif(not (P12 / "syserr.json").exists(), reason="derived JSON not present")
def test_committed_one_mode_fraction():
    d = json.loads((P12 / "syserr.json").read_text())
    cb, ar = d["legs"]["C_both"], d["legs"]["A_redist"]
    assert abs(cb["one_mode"]["f1"] - 0.80) < 0.02
    assert abs(ar["one_mode"]["f1"] - 0.31) < 0.02
    assert cb["one_mode"]["f1"] > cb["null_sign_scramble"]["max"]
    assert cb["sign_pattern_gK"][0] == cb["sign_pattern_gK"][1]
