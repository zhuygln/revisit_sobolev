"""Paper B R2M robustness: the readings on synthetic records (survives / does
not / gray), the prereg check, and a smoke run of the runner."""
import importlib.util
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


R = _load("paperB_analyse_r2m", "paperB/r2m/analyse_r2m.py")
BANDS = "grizJHK"


def _leg(mags, mode, seed_std=0.01, ev=5.0):
    return dict(mode=mode, mags={b: float(m) for b, m in zip(BANDS, mags)}, mags_seed_std={b: seed_std for b in BANDS},
                L_nu=np.ones(200).tolist(), L_bol=1e40, colors={}, energy=dict(identity_residual=1e-16), n_trapped=0,
                n_coherent_fallback=0, n_interactions=1000, events_per_packet=ev, t_wall=1.0,
                ledger=dict(E_esc=0.5, E_core=0.4, W=0.1))


def _kernel(ng, source):
    edges = np.geomspace(1e14, 1e15, ng + 1)
    return dict(source=source, ng=ng, n_events=100 * ng, empty_rows=0, validate_energy=0.0, edges=edges.tolist(),
                R=np.eye(ng).tolist(), counts=[100.0] * ng, E_in=[100.0] * ng, n_exit_samples=10 * ng,
                serialized_bytes=100 * ng, table_kb=0.1 * ng)


def _record(shift=1.5, dm_a2m=0.05, dm_a2=0.8, seed_std_m=0.01):
    ref = [20.0] * 7; refm = [m - shift for m in ref]
    legs = {"R2build": _leg(ref, "sobolev_dmacro"), "R2": _leg(ref, "sobolev_dmacro"),
            "R2Mbuild": _leg(refm, "sobolev_macro", seed_std_m), "R2M": _leg(refm, "sobolev_macro", seed_std_m, ev=15.0),
            "A2_ng8": _leg([m + dm_a2 for m in refm], "sobolev_group"), "A2M_ng8": _leg([m + dm_a2m for m in refm], "sobolev_group")}
    legs["R2build"]["seeds"] = legs["R2Mbuild"]["seeds"] = [101, 102, 103]
    kernels = {"A2_ng8": _kernel(8, "R2build"), "A2M_ng8": _kernel(8, "R2Mbuild"), "K128M": _kernel(128, "R2M"), "K128Mbuild": _kernel(128, "R2Mbuild")}
    return dict(n=300_000, seeds=[1, 2, 3], build_seeds=[101, 102, 103], ng_check=8, ng_fine=128, macro_W=0.5,
                state="paper4/phase1_benchmarks/P1_t2.json", shell=28, lam_window=[1000.0, 30000.0], n_spec=200, legs=legs, kernels=kernels)


@pytest.fixture(autouse=True)
def _equal_band_fractions(monkeypatch):
    monkeypatch.setattr(R.A.V, "band_fractions", lambda row, leg="R2": {b: 1.0 / 7 for b in BANDS})
    monkeypatch.setattr(R.A.V, "nu_edges", lambda lo, hi, n: np.geomspace(1e14, 3e15, n + 1))


def test_survives_when_the_rebuilt_operator_is_within_threshold_and_the_shift_is_reported():
    out = R.readings({"58CeII": _record(shift=1.5, dm_a2m=0.05, dm_a2=0.8)})
    r = out["per_ion"]["58CeII"]
    assert r["survives"] == "YES" and not r["gray"] and out["survives"]["58CeII"] == "YES"
    assert abs(r["shift_max"] - 1.5) < 1e-12 and all(abs(v + 1.5) < 1e-12 for v in r["shift_R2M_minus_R2"].values())
    assert abs(r["legs"]["A2_ng8"]["band"]["max"] - 0.8) < 1e-12          # the downward-trained operator's reference dependence
    assert r["legs"]["A2M_ng8"]["event_vs_K128M"] >= 0.0 and r["fine_in_vs_out_of_sample"] == 0.0


def test_does_not_survive_when_the_rebuilt_operator_misses():
    out = R.readings({"60NdII": _record(dm_a2m=0.25)})
    assert out["per_ion"]["60NdII"]["survives"] == "NO"


def test_gray_on_the_precision_rule_and_on_the_identity():
    out = R.readings({"60NdII": _record(seed_std_m=0.08)})
    r = out["per_ion"]["60NdII"]
    assert r["survives"] == "GRAY" and any(f.startswith("2:") for f in r["gray"]) and r["dropped_for_precision"] == list(BANDS)
    row = _record(); row["legs"]["A2M_ng8"]["energy"]["identity_residual"] = 1e-6
    r = R.readings({"58CeII": row})["per_ion"]["58CeII"]
    assert r["survives"] == "GRAY" and any(f.startswith("3: A2M_ng8") for f in r["gray"])


def test_prereg_check_refuses_a_different_experiment_and_accepts_a_raised_n():
    row = _record(); assert R.check_prereg(row) == []
    row["n"] = 1_000_000; assert R.check_prereg(row) == []
    row["macro_W"] = 1.0
    with pytest.raises(ValueError):
        R.check_prereg(row)
    row = _record(); row["n"] = 100_000
    assert R.check_prereg(row, strict=False)


from sobolev import atomic_cache as ac                          # noqa: E402


@pytest.mark.skipif(not (ac.CACHE_DIR / "57LaII.npz").exists(), reason="La II cache not built")
def test_runner_smoke_la(tmp_path):
    G = _load("paperB_run_r2m", "paperB/r2m/run_r2m.py")
    row = G.run("57LaII", n=2000, seeds=(1,), build_seeds=(101,), out=tmp_path / "r.json", verbose=False)
    assert set(row["legs"]) == {"R2build", "R2", "R2Mbuild", "R2M", "A2_ng8", "A2M_ng8"}
    assert set(row["kernels"]) == {"A2_ng8", "A2M_ng8", "K128M", "K128Mbuild"}
    assert row["kernels"]["A2M_ng8"]["source"] == "R2Mbuild" and row["kernels"]["K128M"]["source"] == "R2M"
    assert row["macro_W"] == 0.5 and row["legs"]["R2M"]["mode"] == "sobolev_macro" and row["legs"]["R2M"]["seeds"] == [1]
    for tag, leg in row["legs"].items():
        assert abs(leg["energy"]["identity_residual"]) < 1e-10, tag
    for tag, k in row["kernels"].items():
        assert k["validate_energy"] < 1e-12, tag
    assert (tmp_path / "r.json").exists()
