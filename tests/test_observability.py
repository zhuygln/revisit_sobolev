"""Phase 3A `observe.py`: the noise model, the scenario masks and chi2_RT,obs."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paper3" / "phase13_observability"))
sys.path.insert(0, str(ROOT / "tests"))
import observe as O            # noqa: E402
import sensitivity as S        # noqa: E402
from test_sensitivity import synthetic_vecs   # noqa: E402


def test_sigma_at_depth_and_bright_limit():
    assert O.sigma_of(23.5, 23.5, 0.03) == pytest.approx(np.sqrt(0.03 ** 2 + 0.2172 ** 2), abs=1e-4)
    assert O.sigma_of(15.0, 23.5, 0.03) == pytest.approx(0.03, abs=1e-3)      # SNR ~ 1e4: the floor
    assert O.sigma_of(24.5, 23.5, 0.03) > O.sigma_of(23.5, 23.5, 0.03)      # fainter = noisier


def test_scenario_masks():
    epochs = (0.5, 1.0, 2.0, 3.0, 5.0, 7.0)
    keys = [(b, t) for t in epochs for b in S.BANDS]
    m_ref = {k: 20.0 for k in keys}
    m_ref[("g", 0.5)] = 24.0                    # below every depth
    mask, sig = O.scenario_obs("dense", keys, m_ref)
    assert len(mask) == 41 and ("g", 0.5) not in mask
    mask, sig = O.scenario_obs("optical", keys, m_ref)
    assert all(k[0] in "griz" for k in mask) and len(mask) == 23
    mask, sig = O.scenario_obs("sparse", keys, m_ref)
    assert {k[1] for k in mask if k[0] in "griz"} == {1.0, 3.0, 7.0}
    assert {k[1] for k in mask if k[0] in "JHK"} == {2.0, 5.0}
    assert all(sig[k] == pytest.approx(O.sigma_of(20.0, O.SCENARIOS["sparse"]["depth"][k[0]],
                                                  O.SCENARIOS["sparse"]["sys"][k[0]])) for k in mask)


def test_chi2_obs_is_sum_of_d_over_sigma_and_gate():
    vecs = synthetic_vecs((0, 0, 0), orthogonal=True)
    point = (0.01, 0.1, 0.01)
    r = O.analyse(vecs, point, "C_both", "dense", "T0", 0.0)
    assert r["status"] == "ok" and r["N_obs"] == 42
    d, s = np.array(r["d_rt"]), np.array(r["sigma"])
    assert r["chi2_RT_obs"] == pytest.approx(np.sum((d / s) ** 2))
    assert r["chi2_RT_N"] == pytest.approx(r["chi2_RT_obs"] / r["N_obs"])
    assert r["eligible"] and r["detectable"] and r["cls"] == "C-B" and r["survives"]
    nir = np.array([k[0] in "JHK" for k in map(tuple, r["keys"])])
    assert r["nir_share"] == pytest.approx(np.sum((d[nir] / s[nir]) ** 2) / r["chi2_RT_obs"])
    r = O.analyse(vecs, point, "C_both", "optical", "T0", 0.0)
    assert r["N_nir"] == 0 and r["nir_share"] == 0.0


def test_absorbable_error_does_not_survive():
    vecs = synthetic_vecs((0.1, -0.2, 0.8))
    r = O.analyse(vecs, (0.01, 0.1, 0.01), "C_both", "dense", "T0", 0.0)
    assert r["detectable"] and r["cls"] == "C-A" and not r["survives"]
