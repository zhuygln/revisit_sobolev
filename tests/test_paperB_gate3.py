"""Paper B gate G3: the state domain, the frozen support, whole-operator
interpolation (energy exact, union exit support, endpoint limits), the
support diagnostic, the A/B/C/D reading and C1-C3 on synthetic states, and
a smoke run of the runner through an interior interpolation."""
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper3", ROOT / "paperB/gate3"):
    sys.path.insert(0, str(p))

from redistribution import RedistributionKernel                       # noqa: E402
import operators3 as O3                                                 # noqa: E402
import states as S                                                      # noqa: E402


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


G3 = _load("paperB_analyse3", "paperB/gate3/analyse.py")


def _events(n=40_000, n_lines=60, seed=0, shift=1.0):
    rng = np.random.default_rng(seed)
    lines = np.sort(rng.uniform(2e14, 8e14, n_lines))
    nu_in = rng.uniform(2.1e14, 7.9e14, n)
    centre = nu_in * shift * rng.choice([1.0, 0.85, 0.7], n, p=[0.5, 0.3, 0.2])
    idx = np.clip(np.searchsorted(lines, centre) + rng.integers(-2, 3, n), 0, n_lines - 1)
    nu_out = lines[idx]; w_in = rng.uniform(0.5, 1.5, n)
    return dict(nu_in=nu_in, nu_out=nu_out, w_in=w_in, w_out=w_in * nu_in / nu_out), lines


def _kern(ev, ng=16, lo=1.9e14, hi=8.1e14):
    return RedistributionKernel.from_branching_mc(ev["nu_in"], ev["nu_out"], ev["w_in"], ng, nu_lo=lo, nu_hi=hi, w_out=ev["w_out"])


def test_state_domain_and_frozen_support_are_as_preregistered():
    assert [ax for ax in S.AXES] == ["T", "D", "J", "P"] and len(S.all_states()) == 15
    for ax, a in S.AXES.items():
        assert a["interior"] in a["grid"] and a["ref"] in a["bracket"] and a["interior"] not in a["bracket"]
    sup = json.loads((ROOT / "paperB/gate3/gate3_support.json").read_text())
    for ion in S.IONS:
        lo, hi = sup["per_ion"][ion]["nu_lo"], sup["per_ion"][ion]["nu_hi"]
        for key, r in sup["per_state"].items():
            if key.endswith(ion):
                assert lo <= r["nu_lo"] * 0.995 + 1e-3 and hi >= r["nu_hi"] * 1.005 - 1e-3   # the union covers every state
    assert sup["blend3"]["nu_lo"] == min(v["nu_lo"] for v in sup["per_ion"].values())
    # the coordinates are logs of the moved quantity
    assert abs(S.coord_of("T", 3000.0) - np.log(3000.0)) < 1e-12 and abs(S.coord_of("J", 7000.0) - np.log(7000.0)) < 1e-12
    assert abs(S.coord_of("P", 3) - np.log(3 * 86400.0)) < 1e-12


def test_whole_operator_interpolation_is_energy_exact_on_the_union_support_and_recovers_the_endpoints():
    ea, _ = _events(seed=0, shift=1.0); eb, _ = _events(seed=1, shift=0.9)
    ka, kb = _kern(ea), _kern(eb)
    for lam, ref in ((0.0, ka), (1.0, kb)):
        k = O3.interpolate_whole(ka, kb, lam, coord="logT", endpoints=["ref", "T2500"])
        assert k.validate_energy() < 1e-12
        live = ~ref.empty_rows
        assert np.allclose(k.R[live], ref.R[live])
        # at an endpoint the union tables reproduce that endpoint's distribution (absent lines carry zero weight)
        for j in range(k.n_groups):
            va, wa = O3._group_weights(ref, j, "energy"); vu, wu = O3._group_weights(k, j, "energy")
            if va.size:
                m = np.isin(vu, va)
                assert np.allclose(wu[m], wa) and np.allclose(wu[~m], 0.0)
    k = O3.interpolate_whole(ka, kb, 0.4)
    tr = k.metadata["transform"]
    assert tr["kind"] == "interp_whole" and tr["n_exit_union"] >= max(tr["n_exit_a"], tr["n_exit_b"])
    assert tr["rows_both"] + tr["rows_one_endpoint"] == int((~k.empty_rows).sum()) and k.validate_energy() < 1e-12
    out = k.sample_nu_out(ea["nu_in"][:3000], np.random.default_rng(0), rows="energy")
    assert np.all(np.isin(out[np.isfinite(out)], k.disc_vals))
    km = O3.interpolate_matrix(ka, kb, 0.4, "a")
    assert np.allclose(km.R, k.R) and km.disc_vals is ka.disc_vals and km.metadata["transform"]["kind"] == "interp_matrix"
    with pytest.raises(ValueError):
        O3.interpolate_whole(ka, _kern(eb, lo=1e14), 0.5)


def test_support_diagnostic_counts_frequencies_outside_the_frozen_edges():
    ev, _ = _events()
    k = _kern(ev, lo=3e14, hi=7e14)                                       # a support narrower than the events
    k2 = O3.support_diag()(k, ev)
    tr = k2.metadata["transform"]
    expect = float(((ev["nu_in"] < 3e14) | (ev["nu_in"] > 7e14)).mean())
    assert tr["kind"] == "support" and abs(tr["clipped_frac"] - expect) < 1e-12 and 0 < tr["clipped_frac"] < 1
    assert abs(tr["below_frac"] + tr["above_frac"] - tr["clipped_frac"]) < 1e-12


# ---- readings on synthetic state records ----
BANDS = "grizJHK"


def _leg(mags, mode="sobolev_group", seed_std=0.01, fallback=0.0):
    return dict(mode=mode, mags={b: float(m) for b, m in zip(BANDS, mags)}, mags_seed_std={b: seed_std for b in BANDS},
                L_nu=np.ones(200).tolist(), L_bol=1e40, colors={}, energy=dict(identity_residual=1e-16), n_trapped=0,
                n_coherent_fallback=int(fallback * 1000), n_interactions=1000, events_per_packet=5.0, t_wall=1.0)


def _kernel(ng, source, transform=None):
    d = dict(source=source, ng=ng, n_events=100 * ng, empty_rows=0, validate_energy=0.0, edges=np.geomspace(1e14, 1e15, ng + 1).tolist(),
             R=np.eye(ng).tolist(), counts=[100.0] * ng, E_in=[100.0] * ng, n_exit_samples=10 * ng, serialized_bytes=1000, table_kb=1.0)
    if transform:
        d["transform"] = transform
    return d


def _state(axis, value, ion, rec=None, fix=None, whole=None, clipped=0.0, untrained=0.0):
    rec = rec or {2: 0.3, 4: 0.2, 8: 0.05, 16: 0.03, 32: 0.02}
    ref = [20.0] * 7
    legs = {"R2": _leg(ref, "sobolev_dmacro"), "R2build": _leg(ref, "sobolev_dmacro")}; legs["R2build"]["seeds"] = [101, 102, 103]
    kernels = {"K128": _kernel(128, "R2", dict(kind="support", clipped_frac=clipped, clipped_energy_frac=clipped)), "K128build": _kernel(128, "R2build")}
    for n, dm in rec.items():
        legs[f"Arec_ng{n}"] = _leg([m + dm for m in ref]); kernels[f"Arec_ng{n}"] = _kernel(n, "R2build")
    interp = None
    if axis != "ref":
        legs["Afix_ng16"] = _leg([m + fix for m in ref], fallback=untrained); kernels["Afix_ng16"] = _kernel(16, "anchor:ref")
        if whole is not None:
            legs["Aint_ng16"] = _leg([m + whole for m in ref]); legs["AintM_ng16"] = _leg([m + whole + 0.05 for m in ref])
            kernels["Aint_ng16"] = _kernel(16, "interp_whole", dict(kind="interp_whole", rows_both=10, rows_one_endpoint=2, n_exit_a=100, n_exit_b=90, n_exit_union=120))
            kernels["AintM_ng16"] = _kernel(16, "interp_matrix")
            interp = dict(lam=0.5, coord=S.AXES[axis]["coord"], endpoints=["ref", "x"], coord_values=[0, 1, 0.5])
    return dict(gate="G3", axis=axis, value=value, label=S.label(axis, value), ion=ion, coord=None if axis == "ref" else S.AXES[axis]["coord"],
                coord_value=None, n=G3.PREREG["n"][ion], seeds=[1, 2, 3], build_seeds=[101, 102, 103], ng_grid=[2, 4, 8, 16, 32], ng_t=16,
                ng_fine=128, lam_window=[1000.0, 30000.0], n_spec=200, legs=legs, kernels=kernels, interpolation=interp)


@pytest.fixture(autouse=True)
def _equal_band_fractions(monkeypatch):
    monkeypatch.setattr(G3.A.V, "band_fractions", lambda row, leg="R2": {b: 1.0 / 7 for b in BANDS})
    monkeypatch.setattr(G3.A.V, "nu_edges", lambda lo, hi, n: np.geomspace(1e14, 3e15, n + 1))


def _domain(ion, fix_by_axis, whole_by_axis=None, rec=None):
    """Every state of the domain for one ion; fix_by_axis[ax] = the transfer error on that axis."""
    out = [G3.read_state(_state("ref", None, ion, rec=rec))]
    for ax, a in S.AXES.items():
        for v in a["grid"]:
            w = (whole_by_axis or {}).get(ax) if v == a["interior"] else None
            out.append(G3.read_state(_state(ax, v, ion, rec=rec, fix=fix_by_axis[ax], whole=w)))
    return out


def test_classification_a_b_c_d_and_gray_isolation():
    ion = "58CeII"
    s = G3.read_state(_state("T", 3000.0, ion, fix=0.03, whole=0.02)); assert s["classification"] == "A" and s["transfers"]
    s = G3.read_state(_state("T", 3000.0, ion, fix=0.3, whole=0.02)); assert s["classification"] == "B"
    s = G3.read_state(_state("T", 3000.0, ion, fix=0.3, whole=0.3)); assert s["classification"] == "C" and s["exists_at_ng_t"]
    s = G3.read_state(_state("T", 3000.0, ion, rec={n: 0.5 for n in (2, 4, 8, 16, 32)}, fix=0.3, whole=0.3)); assert s["classification"] == "D"
    # the transfer leg's fallback is the mechanism, never gray; it is reported as rows never trained
    s = G3.read_state(_state("T", 4000.0, ion, fix=0.3, untrained=0.2))
    assert not s["gray"] and s["classification"] == "C" and abs(s["transfer"]["rows_never_trained_frac"] - 0.2) < 1e-12
    # a clipped fraction above 5 % flags coverage-limited, and is never gray either
    s = G3.read_state(_state("T", 4000.0, ion, fix=0.03, clipped=0.1))
    assert s["outside_fixed_support"]["coverage_limited"] and not s["gray"] and s["classification"] == "A"
    # G1's conditions still gray the recomputed legs and R2
    row = _state("T", 4000.0, ion, fix=0.03); row["legs"]["Arec_ng16"]["energy"]["identity_residual"] = 1e-6
    s = G3.read_state(row); assert s["gray"] and s["classification"] == "GRAY"


def test_c1_c2_c3_readings_over_the_domain():
    ok = {ax: 0.03 for ax in S.AXES}
    states = _domain("58CeII", ok) + _domain("60NdII", ok)
    out = G3.readings(states)
    assert out["C1"] == "GREEN" and out["C2"] == "GREEN" and out["C3"] == "GREEN" and out["step4"] == "transfers across states"
    # transfer fails on T and J for Ce, interpolation passes there -> C2 Yellow, tabulable; two axes needed -> C3 Green
    fail = dict(ok, T=0.3, J=0.3)
    states = _domain("58CeII", fail, whole_by_axis={"T": 0.02, "J": 0.02}) + _domain("60NdII", ok)
    out = G3.readings(states)
    ce = out["per_ion"]["58CeII"]
    assert ce["C2"] == "YELLOW" and ce["needed_axes"] == ["T", "J"] and ce["C3"] == "GREEN"
    assert out["C2"] == "YELLOW" and out["step4"].startswith("exists everywhere")
    # transfer fails on an axis where interpolation also fails -> C2 Red
    states = _domain("58CeII", fail, whole_by_axis={"T": 0.3, "J": 0.02}) + _domain("60NdII", ok)
    assert G3.readings(states)["per_ion"]["58CeII"]["C2"] == "RED"
    # all three controlled axes needed but interpolable -> C3 Yellow; two interpolation failures -> C3 Red
    fail3 = dict(ok, T=0.3, D=0.3, J=0.3)
    assert G3.readings(_domain("58CeII", fail3, {"T": 0.02, "D": 0.02, "J": 0.02}) + _domain("60NdII", ok))["per_ion"]["58CeII"]["C3"] == "YELLOW"
    assert G3.readings(_domain("58CeII", fail3, {"T": 0.3, "D": 0.3, "J": 0.02}) + _domain("60NdII", ok))["per_ion"]["58CeII"]["C3"] == "RED"
    # existence failing at two states -> C1 Red for that ion
    bad = _domain("60NdII", ok, rec={n: 0.5 for n in (2, 4, 8, 16, 32)})
    assert G3.readings(_domain("58CeII", ok) + bad)["per_ion"]["60NdII"]["C1"] == "RED"


from sobolev import atomic_cache as ac                                  # noqa: E402


@pytest.mark.skipif(not (ac.CACHE_DIR / "57LaII.npz").exists(), reason="La II cache not built")
def test_runner_smoke_through_an_interior_interpolation(tmp_path, monkeypatch):
    monkeypatch.setenv("G3_KDIR", str(tmp_path / "k")); monkeypatch.setenv("G3_OUT", str(tmp_path))
    R = _load("paperB_run_gate3", "paperB/gate3/run_gate3.py")
    kw = dict(n=2000, seeds=(1,), build_seeds=(101,), ng_grid=(4, 16), verbose=False)
    ref = R.run_state("ref", None, "57LaII", **kw)
    assert "Afix_ng16" not in ref["legs"] and (tmp_path / "k" / "ref_57LaII_rec_ng16.npz").exists()
    with pytest.raises(FileNotFoundError):
        R.run_state("T", 3000.0, "57LaII", **kw)                        # the endpoint T2500 has not run
    end = R.run_state("T", 2500.0, "57LaII", **kw)
    assert end["kernels"]["Afix_ng16"]["injected"]["kind"] == "anchor" and end["kernels"]["Afix_ng16"]["E_in"] is None
    mid = R.run_state("T", 3000.0, "57LaII", **kw)
    assert set(mid["legs"]) >= {"R2build", "R2", "Arec_ng4", "Arec_ng16", "Afix_ng16", "Aint_ng16", "AintM_ng16"}
    assert mid["interpolation"]["endpoints"] == ["T2500", "ref"] and 0.0 < mid["interpolation"]["lam"] < 1.0
    for tag, k in mid["kernels"].items():
        assert k["validate_energy"] < 1e-12, tag
    for tag, leg in mid["legs"].items():
        assert abs(leg["energy"]["identity_residual"]) < 1e-10, tag
    sup = mid["kernels"]["K128"]["transform"]; assert sup["kind"] == "support" and sup["clipped_frac"] == 0.0
    assert np.allclose(mid["kernels"]["Aint_ng16"]["edges"], mid["kernels"]["Arec_ng16"]["edges"])   # one frozen support
    s = G3.read_state(dict(mid, ng_t=16)); assert s["classification"] in "ABCD"
