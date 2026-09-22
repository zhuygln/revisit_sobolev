"""Paper B gate G1: the preregistered readings on synthetic records, eps*,
the event-level metric, and a smoke run of the runner."""
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


A = _load("paperB_analyse", "paperB/gate1/analyse.py")
BANDS = "grizJHK"


def _leg(mags, mode="sobolev_group", ng=None, eps=None, L_nu=None, seed_std=0.01):
    d = dict(mode=mode, mags={b: float(m) for b, m in zip(BANDS, mags)}, mags_seed_std={b: seed_std for b in BANDS},
             L_nu=(np.ones(200) if L_nu is None else L_nu).tolist(), L_bol=1e40, colors={},
             energy=dict(identity_residual=1e-16), n_trapped=0, n_coherent_fallback=0, n_interactions=1000,
             events_per_packet=5.0, t_wall=1.0)
    if eps is not None:
        d["eps"] = eps
    return d


def _kernel(ng, edges_lo=1e14, edges_hi=1e15, R=None, counts=None, E_in=None, source="R2build"):
    edges = np.geomspace(edges_lo, edges_hi, ng + 1)
    R = np.eye(ng) if R is None else np.asarray(R)
    counts = np.full(ng, 100.0) if counts is None else np.asarray(counts, float)
    E_in = counts if E_in is None else np.asarray(E_in, float)
    return dict(source=source, ng=ng, n_events=int(counts.sum()), empty_rows=int((counts <= 0).sum()), validate_energy=0.0,
                edges=edges.tolist(), R=R.tolist(), counts=counts.tolist(), E_in=E_in.tolist(), n_exit_samples=ng * 10,
                serialized_bytes=ng * 100, table_kb=ng * 0.1)


def _record(dm_by_ng, eps_curve, live=("g", "r", "i", "z", "J", "H", "K"), seed_std=0.01):
    """A synthetic record: R2 at mags 20 in every band (all live at 40 Mpc: brighter
    than the limits, each band 1/7 of L_bol), the A2 legs offset by dm_by_ng[ng]
    in every live band, the eps legs by eps_curve[eps]."""
    ref = [20.0] * 7
    legs = {"R2": _leg(ref, "sobolev_dmacro", seed_std=seed_std), "R2build": _leg(ref, "sobolev_dmacro", seed_std=seed_std)}
    legs["R2build"]["seeds"] = [101, 102, 103]
    for ng, dm in dm_by_ng.items():
        legs[f"A2_ng{ng}"] = _leg([m + dm for m in ref], ng=ng)
    for e, dm in eps_curve.items():
        legs[f"E{e:.2f}"] = _leg([m + dm for m in ref], "sobolev_tla", eps=e)
    kernels = {f"A2_ng{ng}": _kernel(ng) for ng in dm_by_ng}; kernels["K128"] = _kernel(128, source="R2"); kernels["K128build"] = _kernel(128)
    row = dict(n=300_000, seeds=[1, 2, 3], build_seeds=[101, 102, 103], ng_grid=sorted(dm_by_ng), ng_fine=128, eps_grid=sorted(eps_curve),
               state="paper4/phase1_benchmarks/P1_t2.json", shell=28, lam_window=[1000.0, 30000.0], n_spec=200, legs=legs, kernels=kernels)
    # the live-band rule needs band fractions: give every band the same L_nu share via `legs[..]["mags"]`
    # (verdict.band_fractions reads a per-band flux; emulate with equal fractions)
    return row


@pytest.fixture(autouse=True)
def _equal_band_fractions(monkeypatch):
    monkeypatch.setattr(A.V, "band_fractions", lambda row, leg="R2": {b: 1.0 / 7 for b in BANDS})
    monkeypatch.setattr(A.V, "nu_edges", lambda lo, hi, n: np.geomspace(1e14, 3e15, n + 1))


GRID = {round(0.05 * k, 2): 0.6 - 0.02 * k for k in range(21)}       # minimum 0.20 at eps = 1: an edge minimum


def test_green_when_two_ions_reach_the_threshold_by_eight_groups_and_eps_is_far():
    dm = {2: 0.5, 4: 0.3, 8: 0.08, 16: 0.05, 32: 0.03}
    eps_curve = {round(0.05 * k, 2): 0.9 - 0.02 * k if k < 10 else 0.7 + 0.02 * (k - 10) for k in range(21)}   # interior min 0.70
    recs = {ion: _record(dm, eps_curve) for ion in A.IONS}
    out = A.readings(recs)
    assert out["B1"] == "GREEN" and out["B2"] == "GREEN" and out["B3"] == "GREEN" and out["decision"] == "CONTINUE"
    assert all(r["ng_star"] == 8 for r in out["per_ion"].values()) and not any(r["gray"] for r in out["per_ion"].values())


def test_red_when_no_ion_reaches_the_threshold_or_eps_is_nearly_as_good():
    dm = {2: 0.5, 4: 0.4, 8: 0.3, 16: 0.25, 32: 0.2}
    eps_curve = {round(0.05 * k, 2): 0.9 - 0.02 * k if k < 10 else 0.7 + 0.02 * (k - 10) for k in range(21)}
    out = A.readings({ion: _record(dm, eps_curve) for ion in A.IONS})
    assert out["B1"] == "RED" and out["decision"] == "STOP"
    dm2 = {2: 0.5, 4: 0.3, 8: 0.08, 16: 0.05, 32: 0.03}
    eps_close = {round(0.05 * k, 2): 0.3 - 0.02 * k if k < 10 else 0.10 + 0.002 * (k - 10) for k in range(21)}   # min 0.10 interior
    out2 = A.readings({ion: _record(dm2, eps_close) for ion in A.IONS})
    assert out2["B1"] == "GREEN" and out2["B2"] == "RED" and out2["decision"] == "STOP"


def test_yellow_and_gray_paths():
    dm = {2: 0.5, 4: 0.3, 8: 0.15, 16: 0.08, 32: 0.03}                # N_g* = 16
    eps_curve = {round(0.05 * k, 2): 0.9 - 0.02 * k if k < 10 else 0.7 + 0.02 * (k - 10) for k in range(21)}
    out = A.readings({ion: _record(dm, eps_curve) for ion in A.IONS})
    assert out["B1"] == "YELLOW" and out["decision"] == "YELLOW"
    # gray: one ion with a luminous band dropped for precision (seed scatter above 0.05) is excluded; two remain -> still read
    recs = {ion: _record({2: 0.5, 4: 0.3, 8: 0.08, 16: 0.05, 32: 0.03}, eps_curve, seed_std=0.01) for ion in A.IONS}
    recs["60NdII"] = _record({2: 0.5, 4: 0.3, 8: 0.08, 16: 0.05, 32: 0.03}, eps_curve, seed_std=0.08)
    out = A.readings(recs)
    assert out["ions_gray"] == ["60NdII"] and out["B1"] == "GREEN"
    # two gray ions -> the gate is gray
    recs["58CeII"] = _record({2: 0.5, 4: 0.3, 8: 0.08, 16: 0.05, 32: 0.03}, eps_curve, seed_std=0.08)
    assert A.readings(recs)["decision"] == "GRAY"
    # a local wiggle in the band error is Yellow, not Red (revised B3)
    wiggle = {2: 0.5, 4: 0.3, 8: 0.04, 16: 0.06, 32: 0.02}
    out = A.readings({ion: _record(wiggle, eps_curve) for ion in A.IONS})
    assert out["B1"] == "GREEN" and out["B3"] == "YELLOW" and out["decision"] == "CONTINUE"
    # a persistent worsening of both observables from 8 to 32 is Red
    worse = {2: 0.5, 4: 0.3, 8: 0.05, 16: 0.10, 32: 0.30}
    recs = {ion: _record(worse, eps_curve) for ion in A.IONS}
    for r in recs.values():                                            # make the SED worse at 32 than at 8 too
        r["legs"]["A2_ng32"]["L_nu"] = (np.ones(200) * 1.5).tolist(); r["legs"]["A2_ng8"]["L_nu"] = (np.ones(200) * 1.01).tolist()
    out = A.readings(recs)
    assert out["B3"] == "RED" and out["decision"] != "CONTINUE"


def test_eps_star_takes_the_grid_minimum_and_flags_an_edge():
    row = _record({2: 0.5, 4: 0.08, 8: 0.05, 16: 0.04, 32: 0.03}, GRID)
    M = A.metrics(row); es = A.eps_star(M, row)
    assert es["eps"] == 1.0 and not es["interior"] and abs(es["mean_dm"] - 0.20) < 1e-12
    assert es["edge_gray"] is False                    # the neighbour differs by 0.02 > sqrt(2) * 0.01
    flat = dict(GRID); flat[0.95] = 0.205
    es2 = A.eps_star(A.metrics(_record({2: 0.5}, flat)), None)
    assert es2["edge_gray"] is True


def test_event_metric_is_zero_at_full_resolution_energy_weighted_and_bounded():
    fine = _kernel(8, R=np.eye(8))
    assert A.m_event(fine, fine) == 0.0
    coarse = _kernel(2, R=np.eye(2))
    assert 0.0 < A.m_event(coarse, fine) <= 1.0
    # a coarse matrix that is exactly the block average of a fine one, with uniform within-block rows, expands to the fine one
    Rf = np.zeros((4, 4)); Rf[:2, :2] = 0.5; Rf[2:, 2:] = 0.5
    fine4 = _kernel(4, R=Rf); coarse2 = _kernel(2, R=np.eye(2))
    assert A.m_event(coarse2, fine4) < 1e-12
    # energy weighting: a row with no absorbed energy does not count, one that carries all of it dominates
    Rf2 = np.eye(4); Rf2[0] = [0.0, 1.0, 0.0, 0.0]                  # fine row 0 disagrees with the coarse identity
    fine_w = _kernel(4, R=Rf2, E_in=[1.0, 0.0, 0.0, 0.0]); fine_u = _kernel(4, R=Rf2, E_in=[0.0, 1.0, 1.0, 1.0])
    c4 = _kernel(4, R=np.eye(4))
    assert A.m_event(c4, fine_w) == 1.0 and A.m_event(c4, fine_u) == 0.0


def test_check_prereg_refuses_a_different_experiment():
    row = _record({2: 0.5, 4: 0.3, 8: 0.08, 16: 0.05, 32: 0.03}, GRID)
    assert A.check_prereg(row) == []
    row["n"] = 1000
    with pytest.raises(ValueError):
        A.check_prereg(row)


@pytest.mark.skipif(not (ROOT / "data/cache/57LaII.npz").exists(), reason="La II cache not built")
def test_runner_smoke_on_la_ii(tmp_path):
    R = _load("paperB_run", "paperB/gate1/run_gate1.py")
    row = R.run("57LaII", n=1500, seeds=(1,), ng_grid=(4,), eps_grid=(0.0, 1.0), ng_fine=16, out=tmp_path / "s.json", verbose=False, build_seeds=(101,))
    assert set(row["legs"]) == {"R2build", "R2", "E0.00", "E1.00", "A2_ng4"} and set(row["kernels"]) == {"A2_ng4", "K16", "K16build"}
    for tag, leg in row["legs"].items():
        assert abs(leg["energy"]["identity_residual"]) < 1e-10, tag
    assert row["legs"]["R2build"]["seeds"] == [101] and row["legs"]["R2"]["seeds"] == [1]
    assert row["kernels"]["A2_ng4"]["source"] == "R2build" and row["kernels"]["K16"]["source"] == "R2"
    for k in row["kernels"].values():
        assert k["validate_energy"] < 1e-12 and k["serialized_bytes"] > 0 and len(k["E_in"]) == k["ng"]
    assert row["legs"]["E1.00"]["eps"] == 1.0 and row["ng_grid"] == [4] and row["eps_grid"] == [0.0, 1.0] and row["build_seeds"] == [101]
    assert "n_coherent_fallback" in row["legs"]["A2_ng4"]
    M = A.metrics(row)
    assert M["fine_in_vs_out_of_sample"] >= 0.0 and 0.0 <= M["legs"]["A2_ng4"]["event"] <= 1.0
