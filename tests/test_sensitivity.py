"""`sensitivity.py` on a synthetic grid with a known answer."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paper3" / "phase12_grid"))
import sensitivity as S
from run_grid import M_GRID, V_GRID, X_GRID


def synthetic_vecs(a_true, noise=0.003, orthogonal=False, seed=0):
    """A smooth m(M, v, X) on the grid, d_RT = sum a_theta d_theta + noise."""
    rng = np.random.default_rng(seed)
    epochs = (0.5, 1.0, 2.0, 3.0, 5.0, 7.0)
    keys = [(b, t) for t in epochs for b in S.BANDS]
    # per-observable sensitivities (mag per e-fold), different per band/epoch
    cM = {k: -0.4 + 0.05 * i for i, k in enumerate(keys)}
    cv = {k: 0.3 * np.sin(i) for i, k in enumerate(keys)}
    cX = {k: 0.8 * (1 if k[0] in "griz" else -0.5) * (1 + 0.1 * k[1]) for k in keys}

    def m(M, v, X):
        return {k: 20.0 + cM[k] * np.log(M / 0.01) + cv[k] * np.log(v / 0.1) + cX[k] * np.log(X / 0.01)
                for k in keys}

    vecs = {}
    for M in M_GRID:
        for v in V_GRID:
            for X in X_GRID:
                mm = m(M, v, X)
                d_rt = {}
                for leg in S.LEGS:
                    if leg == "A_redist":
                        d_rt[leg] = {k: rng.normal(0, noise) for k in keys}
                    elif orthogonal:
                        # a pattern with zero weighted projection on all three directions
                        d_rt[leg] = {k: 0.3 * (1 if i % 2 else -1) * (1 + 0.01 * i) for i, k in enumerate(keys)}
                    else:
                        d_rt[leg] = {k: a_true[0] * cM[k] + a_true[1] * cv[k] + a_true[2] * cX[k]
                                     + rng.normal(0, noise) for k in keys}
                vecs[(M, v, X)] = (keys, mm, d_rt, {k: True for k in keys})
    return vecs


def test_recovers_projection_on_linear_grid():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    p = S.analyse_point(vecs, (0.01, 0.1, 0.01), "C_both", noise_floor=0.0)
    assert p["status"] == "ok"
    assert p["a"] == pytest.approx([0.1, -0.2, 0.8], abs=0.02)
    assert p["R"] < 0.1 and p["cls"] == "C-A"
    assert p["derivative_info"][0]["scheme"] == "central"


def test_orthogonal_error_is_distinct():
    vecs = synthetic_vecs((0, 0, 0), orthogonal=True)
    p = S.analyse_point(vecs, (0.01, 0.1, 0.01), "C_both", noise_floor=0.0)
    assert p["status"] == "ok"
    assert p["cls"] == "C-B" and p["R"] > 0.5
    assert p["chi2_RT_N"] > 4


def test_small_error_is_cc():
    vecs = synthetic_vecs((0.001, 0.0, 0.002))
    p = S.analyse_point(vecs, (0.01, 0.1, 0.01), "C_both", noise_floor=0.0)
    assert p["cls"] == "C-C"


def test_edge_point_uses_one_sided_and_flags_it():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    p = S.analyse_point(vecs, (0.003, 0.05, 0.001), "C_both", noise_floor=0.0)
    assert p["status"] == "ok"
    assert all(i["one_sided"] for i in p["derivative_info"])
    assert p["a"] == pytest.approx([0.1, -0.2, 0.8], abs=0.02)   # the synthetic m is linear in ln theta


def test_mask_intersection_and_secants():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    keys, m0, d, mask = vecs[(0.01, 0.1, 0.01)]
    mask[keys[0]] = False
    p = S.analyse_point(vecs, (0.01, 0.1, 0.01), "C_both", noise_floor=0.0)
    assert p["N"] == len(keys) - 1
    s = S.secant_disagreement(vecs, (0.01, 0.1, 0.01), "X")
    assert s["disagreement"] == pytest.approx(0.0, abs=1e-9)


def test_classify_thresholds():
    assert S.classify({"chi2_RT_N": 3.0, "R": 0.9, "signif": [0, 0, 0]}) == "C-C"
    assert S.classify({"chi2_RT_N": 30.0, "R": 0.2, "signif": [1, 5, 0]}) == "C-A"
    assert S.classify({"chi2_RT_N": 30.0, "R": 0.5, "signif": [1, 5, 0]}) == "C-B"
    assert S.classify({"chi2_RT_N": 30.0, "R": 0.2, "signif": [1, 2, 0]}) == "C-B"
