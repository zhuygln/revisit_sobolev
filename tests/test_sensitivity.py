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
                vecs[(M, v, X)] = (keys, mm, d_rt, {k: True for k in keys},
                                   {"T_eff": {t: 5000.0 for t in epochs}, "floored_rows": {}})
    return vecs


def test_recovers_projection_on_linear_grid():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    p = S.analyse_point(vecs, (0.01, 0.1, 0.01), "C_both", noise_floor=0.0)
    assert p["status"] == "ok"
    assert p["a"] == pytest.approx([0.1, -0.2, 0.8], abs=0.02)
    assert p["R"] < 0.1 and p["cls"] == "C-A"
    assert p["derivative_info"][0]["scheme"] == "central"
    assert not p["extrapolated"] and max(p["a_over_dln"]) < 1.0


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
    keys, m0, d, mask = vecs[(0.01, 0.1, 0.01)][:4]
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


# --- §4.41 generalization: floor mask, core convention, nuisance columns ---------

def _with_error(vecs, point, err):
    """Overwrite C_both's error at `point` with err[(band, t)]."""
    keys, m0, d, mask, info = vecs[point]
    d["C_both"] = {k: err[k] for k in keys}


def test_floored_mask_from_flag_and_from_v_core():
    model = {"v_ej_c": 0.1}
    assert S.row_floored({"v_core": 0.05, "source": {}}, model)
    assert not S.row_floored({"v_core": 0.09, "source": {}}, model)
    assert S.row_floored({"v_core": 0.09, "source": {"v_ph_floored": True}}, model)
    assert not S.row_floored({"v_core": 0.05, "source": {"v_ph_floored": False}}, model)


def test_floored_rows_change_N():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    point = (0.01, 0.1, 0.01)
    p0 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0)
    keys, m0, d, mask, info = vecs[point]
    for k in keys:
        if k[1] == 7.0:
            mask[k] = False          # what vectors(floored="exclude") does to a floored epoch
    p1 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0)
    assert p0["N"] - p1["N"] == len(S.BANDS)


def test_core_shift_is_grey_and_matches_dm_bol():
    o = {"L_bol": 2.0e40, "L_bol_absorbing": 1.0e40}
    assert S.core_shift(o, "conserving") == 0.0
    assert S.core_shift(o, "absorbing") == pytest.approx(2.5 * np.log10(2.0))
    with pytest.raises(ValueError):
        S.core_shift(o, "equilibrium")


def test_grey_per_epoch_error_absorbed_under_T1_not_T0():
    vecs = synthetic_vecs((0, 0, 0), orthogonal=True)
    point = (0.01, 0.1, 0.01)
    keys = vecs[point][0]
    _with_error(vecs, point, {k: 0.5 * np.sin(k[1]) for k in keys})    # grey per epoch
    p0 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0, nuisance=())
    p1 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0, nuisance=("L_t",))
    assert p0["cls"] == "C-B" and p0["R"] > 0.5
    assert p1["R"] < 1e-6 and p1["chi2_res_dof"] < 1e-6
    assert p1["cls"] != "C-A"            # nothing physical is significant
    assert p1["columns"][3:] == [f"L_{t:g}" for t in (0.5, 1, 2, 3, 5, 7)]
    assert p1["R_nuisance_only"] < 1e-6


def test_planck_temperature_error_absorbed_under_T2_not_T1():
    vecs = synthetic_vecs((0, 0, 0), orthogonal=True)
    point = (0.01, 0.1, 0.01)
    keys = vecs[point][0]
    dT = S.tbb_derivative(5000.0)
    _with_error(vecs, point, {k: 0.3 * dT[k[0]] for k in keys})      # a 35 % hotter photosphere
    p1 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0, nuisance=("L_t",))
    p2 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0, nuisance=("L_t", "T_bb"))
    assert p1["R"] > 0.1                 # the synthetic X direction is itself a blue/NIR colour
    assert p2["R"] < 1e-6
    assert p2["a_nuisance"]["T_bb"] == pytest.approx(0.3, abs=1e-6)
    assert p2["cls"] != "C-A"


def test_tbb_derivative_is_fixed_L_colour_direction():
    d = S.tbb_derivative(5000.0)
    assert d["g"] < d["r"] < d["i"] < d["z"] < d["J"] < d["H"] < d["K"]   # hotter = bluer
    assert d["g"] < 0 < d["K"]


def test_orthogonal_error_stays_CB_under_T3_and_nuisance_significance_never_makes_CA():
    vecs = synthetic_vecs((0, 0, 0), orthogonal=True)
    point = (0.01, 0.1, 0.01)
    p3 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0, nuisance=("L_t", "T_bb", "2c"))
    assert p3["status"] == "ok" and "2c" in p3["nuisance_used"]
    assert p3["cls"] == "C-B" and p3["R"] > 0.3
    # a fake point where only a nuisance column is significant
    fake = dict(p3, R=0.1, chi2_RT_N=50.0, signif=[0.1, 0.2, 0.3] + [10.0] * (p3["p"] - 3))
    assert S.classify(fake, n_phys=3) == "C-B"
    assert S.classify(dict(fake, signif=[5.0, 0, 0] + [10.0] * (p3["p"] - 3)), n_phys=3) == "C-A"


def test_2c_skipped_at_blue_point_and_lin_measure():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    p = S.analyse_point(vecs, (0.01, 0.1, 0.001), "C_both", noise_floor=0.0, nuisance=("2c",))
    assert "2c" not in p["nuisance_used"] and "skipped" in p["nuisance_notes"]["2c"]
    p = S.analyse_point(vecs, (0.01, 0.1, 0.01), "C_both", noise_floor=0.0, nuisance=("2c",))
    assert "2c" in p["nuisance_used"] and p["lin_2c"] >= 0


def test_dof_and_underdetermined():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    point = (0.01, 0.1, 0.01)
    keys, m0, d, mask, info = vecs[point]
    for k in keys:
        mask[k] = k[1] in (1.0, 3.0) and k[0] in "griz" and k != ("z", 3.0)     # 7 observables
    p0 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0)
    assert p0["N"] == 7 and p0["dof"] == 4 and not p0["underdetermined"] and p0["low_N"]
    p1 = S.analyse_point(vecs, point, "C_both", noise_floor=0.0, nuisance=("L_t",))
    # the synthetic X column is grey within griz at one epoch, so it is collinear
    # with the two L_t columns: rank 4 of 5, dof 3 -> underdetermined
    assert p1["p"] == 5 and p1["rank"] == 4 and p1["dof"] == 3 and p1["cls"] == "underdetermined"


def test_dead_column_has_zero_significance():
    d = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    D = np.array([d, np.zeros(5)]).T
    p = S.project(d, D, np.full(5, 0.1))
    assert p["signif"][1] == 0.0 and p["rank"] == 1 and p["dof"] == 4


GRID = ROOT / "paper3" / "phase12_grid" / "grid"


@pytest.mark.skipif(not (GRID / "model_M0.01_v0.1_X0.01.json").exists(), reason="grid JSON not present")
def test_absorbing_core_rederivation_matches_stored_dm_bol():
    """The absorbing-core per-band shift is the grey factor L_bol/L_bol_absorbing
    (per leg, minus the reference's); it must reproduce the stored
    `dm_bol_absorbing` on every committed row of the central model."""
    import json
    model = json.loads((GRID / "model_M0.01_v0.1_X0.01.json").read_text())
    n = 0
    for r in model["rows"]:
        if r.get("status") not in ("ok", "reduced_n"):
            continue
        s_ref = S.core_shift(r["ref"], "absorbing")
        for leg, o in r["legs"].items():
            assert o["dm_bol"] + S.core_shift(o, "absorbing") - s_ref == pytest.approx(o["dm_bol_absorbing"], abs=1e-9)
            n += 1
    assert n > 0
    keys, m_cons, d_cons, _, _ = S.vectors(model, core="conserving")
    keys, m_abs, d_abs, _, _ = S.vectors(model, core="absorbing")
    for k in keys:
        if np.isfinite(m_cons[k]):
            assert m_abs[k] >= m_cons[k] - 1e-9          # absorbing core is never brighter
    # colours are convention-invariant: d_abs - d_cons is grey within an epoch
    for leg in S.LEGS:
        for t in model["epochs"]:
            sh = [d_abs[leg][(b, t)] - d_cons[leg][(b, t)] for b in S.BANDS if np.isfinite(d_cons[leg][(b, t)])]
            if sh:
                assert np.ptp(sh) < 1e-9
