"""Paper B G3U (paired-seed uncertainty): the affected set is derived from
the frozen G3 records by the preregistered window, the paired bootstrap
recovers a known interval on synthetic seeds, and a case is gray when the
operator does not match the frozen record or the seed scatter differs."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


U = _load("paperB_analyse_g3u", "paperB/gate3u/analyse.py")
BANDS = "grizJHK"


def test_affected_set_is_the_near_threshold_decision_legs_of_the_frozen_records():
    if not (ROOT / "paperB/gate3/gate3_partc_p1blend.json").exists():
        pytest.skip("no G3 records")
    aff = U.affected()
    assert len(aff) == 33 and len({a["record"] for a in aff}) == 18
    for a in aff:
        assert 0.07 <= a["band"] <= 0.13 or 0.07 <= a["colour"] <= 0.13
        assert a["ion"] in ("58CeII", "60NdII") or a["part"] in ("b", "c")
    legs = {(a["record"], a["leg"]) for a in aff}
    assert ("gate3_P_P3d_60NdII.json", "Aint_ng16") in legs          # the Nd trajectory boundary
    assert ("gate3_J_J2500_58CeII.json", "Afix_ng16") in legs        # the colour-only J-axis failure
    assert ("gate3_partc_p1blend.json", "Arec_ng4") in legs          # C5's k_rec
    assert not any(a["ion"] == "57LaII" for a in aff)
    assert U.affected(window=(0.0999, 0.1001)) == []


def _rows(n_seeds=12, shift=0.06, sd=0.02, seed=0, live=("g", "r", "i", "z", "J", "H", "K")):
    """Synthetic closure and reference legs whose paired difference is shift + noise(sd) in every band."""
    rng = np.random.default_rng(seed)
    ref = [{b: 20.0 + 0.5 * i + rng.normal(0, 0.3) for i, b in enumerate(BANDS)} for _ in range(n_seeds)]
    clo = [{b: m[b] + shift + rng.normal(0, sd) for b in BANDS} for m in ref]
    mk = lambda ps: dict(mags_per_seed=ps, mags={b: float(np.mean([p[b] for p in ps])) for b in BANDS},
                         mags_seed_std={b: float(np.std([p[b] for p in ps[:3]], ddof=1)) for b in BANDS})
    return mk(clo), mk(ref)


def test_paired_bootstrap_recovers_the_interval_of_a_known_shift():
    clo, ref = _rows(shift=0.06, sd=0.02)
    bands, cols, pairs = U.paired_diffs(clo, ref, list(BANDS))
    assert bands.shape == (12, 7) and cols.shape[1] == len(pairs) == 5
    e = U.e_joint_of(bands, cols)
    assert 0.04 < e < 0.09
    bs = U.bootstrap(bands, cols, n_boot=4000, seed=1)
    assert bs["ci95"][0] <= e <= bs["ci95"][1] and bs["ci68"][0] >= bs["ci95"][0] and bs["ci68"][1] <= bs["ci95"][1]
    assert bs["ci95"][1] < 0.10 and bs["p_pass"] > 0.95             # decided: well below the threshold
    # a shift just under the threshold with more scatter is within noise; the max over bands and colours
    # is biased upward (the winner's curse the interval is meant to carry), so the point sits above the shift
    clo2, ref2 = _rows(shift=0.09, sd=0.05, seed=3)
    b2, c2, _ = U.paired_diffs(clo2, ref2, list(BANDS))
    e2 = U.e_joint_of(b2, c2)
    bs2 = U.bootstrap(b2, c2, n_boot=4000, seed=1)
    assert e2 > 0.09 and bs2["ci95"][0] < 0.10 < bs2["ci95"][1] and 0.05 < bs2["p_pass"] < 0.95
    # reproducible
    assert U.bootstrap(bands, cols, n_boot=500, seed=7) == U.bootstrap(bands, cols, n_boot=500, seed=7)


def _frozen_and_rerun(match=0.0, scatter_dev=0.0):
    clo, ref = _rows(shift=0.06, sd=0.02)
    frozen = dict(n=3000, seeds=[1, 2, 3], lam_window=[1000.0, 30000.0], n_spec=200, ng_fine=16,
                  legs={"R2": dict(ref, mode="sobolev_dmacro", energy=dict(identity_residual=0.0), n_trapped=0, n_coherent_fallback=0,
                                   n_interactions=10, events_per_packet=1.0, t_wall=1.0, L_nu=[1.0] * 200, L_bol=1.0, colors={}),
                        "Arec_ng16": dict(clo, mode="sobolev_group", energy=dict(identity_residual=0.0), n_trapped=0, n_coherent_fallback=0,
                                          n_interactions=10, events_per_packet=1.0, t_wall=1.0, L_nu=[1.0] * 200, L_bol=1.0, colors={})},
                  kernels={"Arec_ng16": dict(R=[[1.0]], ng=16, validate_energy=0.0, empty_rows=0, source="R2build", n_exit_samples=1, table_kb=1.0),
                           "K16": dict(R=[[1.0]], edges=[1.0, 2.0], E_in=[1.0], ng=16, validate_energy=0.0, empty_rows=0, source="R2", n_exit_samples=1, table_kb=1.0)})
    frozen["legs"]["R2"]["mags_seed_std"] = {b: v + scatter_dev for b, v in frozen["legs"]["R2"]["mags_seed_std"].items()}
    rerun = dict(frozen, affected_legs=["Arec_ng16"], operator_match={"Arec_ng16": match})
    return rerun, frozen


def test_read_case_is_gray_on_operator_mismatch_or_scatter_mismatch(monkeypatch):
    # the frozen live bands come from G1's rule; here every band is live
    monkeypatch.setattr(U.A, "live_bands", lambda row, leg="R2": (list(BANDS), [], list(BANDS)))
    monkeypatch.setattr(U.A, "metrics", lambda row: dict(legs={"Arec_ng16": dict(band=dict(max=0.06), colour=dict(max=0.05), identity=0.0,
                                                                                  kernel_energy=0.0, trapped_frac=0.0, fallback_frac=0.0)}))
    rerun, frozen = _frozen_and_rerun()
    c = U.read_case(rerun, frozen, "Arec_ng16")
    assert c["gray"] == [] and c["status"] == "DECIDED_PASS" and c["n_seeds"] == 12 and c["frozen"]["passes"] is True
    rerun, frozen = _frozen_and_rerun(match=1e-6)
    assert "operator" in " ".join(U.read_case(rerun, frozen, "Arec_ng16")["gray"])
    rerun, frozen = _frozen_and_rerun(scatter_dev=1e-3)
    assert "scatter" in " ".join(U.read_case(rerun, frozen, "Arec_ng16")["gray"])
