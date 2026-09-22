"""Paper B gate G2: the derived-operator utilities, the three families on
synthetic events, the readings H1/H2 on synthetic records, and a smoke run."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper3", ROOT / "paperB/gate2"):
    sys.path.insert(0, str(p))

from redistribution import RedistributionKernel                      # noqa: E402
import operators as OP                                                 # noqa: E402


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


G2 = _load("paperB_analyse2", "paperB/gate2/analyse.py")
BANDS = "grizJHK"


def _events(n=60_000, n_lines=80, seed=0):
    """Exits on discrete lines; the exit distribution depends smoothly on the
    absorbed frequency (a local law), energy conserved per event."""
    rng = np.random.default_rng(seed)
    lines = np.sort(rng.uniform(2e14, 8e14, n_lines))
    nu_in = rng.uniform(2.1e14, 7.9e14, n)
    # exit line drawn near the absorbed frequency (locality) with a red tail
    centre = nu_in * rng.choice([1.0, 0.85, 0.7], n, p=[0.5, 0.3, 0.2])
    idx = np.clip(np.searchsorted(lines, centre) + rng.integers(-2, 3, n), 0, n_lines - 1)
    nu_out = lines[idx]
    w_in = rng.uniform(0.5, 1.5, n); w_out = w_in * nu_in / nu_out
    return dict(nu_in=nu_in, nu_out=nu_out, w_in=w_in, w_out=w_out), lines


def _fine(ev, ng=16):
    return RedistributionKernel.from_branching_mc(ev["nu_in"], ev["nu_out"], ev["w_in"], ng, nu_lo=2e14, nu_hi=8e14, w_out=ev["w_out"])


def test_with_matrix_keeps_energy_exact_and_empty_rows_empty():
    ev, _ = _events(); k = _fine(ev)
    R2 = np.full_like(k.R, 1.0 / k.n_groups)
    k2 = k.with_matrix(R2, dict(transform=dict(kind="test")))
    assert k2.validate_energy() < 1e-12 and np.array_equal(k2.empty_rows, k.empty_rows)
    assert np.all(k2.R[k.empty_rows] == 0) and k2.metadata["transform"]["kind"] == "test"
    assert k2.disc_vals is k.disc_vals and np.allclose(k2.E_cum[~k.empty_rows][:, -1], 1.0)
    with pytest.raises(ValueError):
        k.with_matrix(np.zeros((3, 3)))


def test_truncate_exits_identity_at_one_subset_below_and_renormalised():
    ev, lines = _events(); k = _fine(ev)
    k1, n1 = k.truncate_exits(1.0)
    assert n1 == k.disc_vals.size and np.array_equal(k1.disc_vals, k.disc_vals) and np.allclose(k1.disc_cum_E, k.disc_cum_E)
    k9, n9 = k.truncate_exits(0.9)
    assert 0 < n9 < n1 and np.all(np.isin(k9.disc_vals, k.disc_vals)) and np.array_equal(k9.R, k.R)
    for j in range(k9.n_groups):
        a, b = k9.disc_off[j], k9.disc_off[j + 1]
        if b > a:
            assert abs(k9.disc_cum[b - 1] - 1.0) < 1e-12 and abs(k9.disc_cum_E[b - 1] - 1.0) < 1e-12
            assert np.all(np.diff(k9.disc_vals[a:b]) > 0)                  # frequency order kept
    # every populated group keeps at least one line, and the kept lines carry >= 90 % of the group's energy
    for j in range(k.n_groups):
        a, b = k.disc_off[j], k.disc_off[j + 1]
        if b > a:
            assert k9.disc_off[j + 1] > k9.disc_off[j]
            w = np.diff(np.hstack([0.0, k.disc_cum_E[a:b]]))
            kept = np.isin(k.disc_vals[a:b], k9.disc_vals[k9.disc_off[j]:k9.disc_off[j + 1]])
            assert w[kept].sum() >= 0.9 - 1e-9
    # transport samples only kept lines
    out = k9.sample_nu_out(ev["nu_in"][:5000], np.random.default_rng(0), rows="energy")
    assert np.all(np.isin(out[np.isfinite(out)], k9.disc_vals))


def test_nmf_rank_one_is_an_outer_product_and_error_falls_with_rank():
    ev, _ = _events(); k = _fine(ev)
    live = ~k.empty_rows; V = k.R[live]
    W1, H1, c1 = OP.nmf(V, 1); W4, H4, c4 = OP.nmf(V, 4)
    assert np.linalg.matrix_rank(W1 @ H1) == 1 and c4["rel_frobenius"] < c1["rel_frobenius"]
    assert c1["converged"] and c4["converged"] and c4["iters"] <= c4["max_iters"]
    g = OP.nmf_rank(4)(k, ev)
    tr = g.metadata["transform"]
    assert tr["kind"] == "nmf" and tr["archetypes"] == 4 and tr["n_params"] == 2 * k.n_groups * 4 and g.validate_energy() < 1e-12
    assert np.allclose(g.R[live].sum(axis=1), V.sum(axis=1)) and 0.0 <= tr["in_sample_tv"] <= 1.0
    assert tr["n_rows_zeroed"] == 0 and tr["zeroed_energy_share"] == 0.0
    # more archetypes reproduce the events better
    tv = [OP.nmf_rank(kk)(k, ev).metadata["transform"]["in_sample_tv"] for kk in (1, 2, 8)]
    assert tv[0] >= tv[1] >= tv[2]


def test_a_row_the_factorisation_zeroes_becomes_empty_and_energy_still_closes():
    """The G2 defect found on the La II control: rank 1 and 2 sent one live
    row to zero, the row sum no longer matched q_dep and validate_energy
    returned 1.0 (gray condition 3). Such a row is now declared empty, so
    transport falls back to coherent scattering there, and the identity
    holds exactly."""
    ev, _ = _events(); k = _fine(ev)
    live = np.flatnonzero(~k.empty_rows)
    R = k.R.copy()
    g = k.with_matrix(R)
    assert g.validate_energy() < 1e-12
    counts = k.counts.copy(); counts[live[0]] = 0.0
    g2 = k.with_matrix(R, counts=counts)
    assert g2.empty_rows[live[0]] and np.all(g2.R[live[0]] == 0) and g2.validate_energy() < 1e-12
    out = g2.sample_nu_out(np.full(200, np.sqrt(g2.edges[live[0]] * g2.edges[live[0] + 1])), np.random.default_rng(0), rows="energy")
    assert np.all(np.isnan(out))                       # the caller scatters these coherently


def test_local_on_fine_tables_is_sampling_equivalent_to_the_coarse_kernel():
    ev, lines = _events(n=120_000); k = _fine(ev, 16)
    coarse = RedistributionKernel.from_branching_mc(ev["nu_in"], ev["nu_out"], ev["w_in"], 4, nu_lo=2e14, nu_hi=8e14, w_out=ev["w_out"])
    loc = OP.local(4)(k, ev)
    assert loc.metadata["transform"]["archetypes"] == 4 and loc.metadata["transform"]["n_params"] == 16 and loc.validate_energy() < 1e-12
    # the fine rows within a coarse block are identical
    assert np.allclose(loc.R[0], loc.R[3]) and not np.allclose(loc.R[0], loc.R[4])
    # the same absorbed frequencies sampled through both: the exit-line histograms agree to the binomial noise
    nu_abs = np.random.default_rng(5).uniform(2.1e14, 7.9e14, 200_000)
    a = coarse.sample_nu_out(nu_abs, np.random.default_rng(1), rows="energy"); b = loc.sample_nu_out(nu_abs, np.random.default_rng(2), rows="energy")
    ha = np.array([np.sum(a == l) for l in lines], float); hb = np.array([np.sum(b == l) for l in lines], float)
    sig = np.sqrt(ha + hb + 1.0)
    assert np.max(np.abs(ha - hb) / sig) < 5.0 and np.isfinite(a).mean() > 0.99


# ---- readings on synthetic records ----
def _leg(mags, mode="sobolev_group", seed_std=0.01):
    return dict(mode=mode, mags={b: float(m) for b, m in zip(BANDS, mags)}, mags_seed_std={b: seed_std for b in BANDS},
                L_nu=np.ones(200).tolist(), L_bol=1e40, colors={}, energy=dict(identity_residual=1e-16), n_trapped=0,
                n_coherent_fallback=0, n_interactions=1000, events_per_packet=5.0, t_wall=1.0)


def _kernel(ng, source, R=None, transform=None, n_exit=1000):
    edges = np.geomspace(1e14, 1e15, ng + 1)
    d = dict(source=source, ng=ng, n_events=100 * ng, empty_rows=0, validate_energy=0.0, edges=edges.tolist(),
             R=(np.eye(ng) if R is None else np.asarray(R)).tolist(), counts=[100.0] * ng, E_in=[100.0] * ng,
             n_exit_samples=n_exit, serialized_bytes=1000, table_kb=1.0)
    if transform:
        d["transform"] = transform
    return d


def _g1(dm_by_ng, ev_by_ng):
    ref = [20.0] * 7
    legs = {"R2": _leg(ref, "sobolev_dmacro"), "R2build": _leg(ref, "sobolev_dmacro")}; legs["R2build"]["seeds"] = [101, 102, 103]
    kernels = {"K128": _kernel(128, "R2"), "K128build": _kernel(128, "R2build")}
    for ng, dm in dm_by_ng.items():
        legs[f"A2_ng{ng}"] = _leg([m + dm for m in ref])
        kernels[f"A2_ng{ng}"] = _kernel(ng, "R2build", R=np.eye(ng))      # block-expanded against the identity K128: a fixed TV distance per N
    row = dict(n=300_000, seeds=[1, 2, 3], build_seeds=[101, 102, 103], ng_grid=sorted(dm_by_ng), ng_fine=128, eps_grid=[0.0],
               state="paper4/phase1_benchmarks/P1_t2.json", shell=28, lam_window=[1000.0, 30000.0], n_spec=200, legs=legs, kernels=kernels)
    return row


def _g2(ion, dm_global, dm_trunc, rho_trunc, ev_global, ref_off=0.0, ctrl_off=0.0, conv=1e-6):
    ref = [20.0 + ref_off] * 7
    legs = {"R2": _leg(ref, "sobolev_dmacro"), "R2build": _leg(ref, "sobolev_dmacro"), "A2_ng8": _leg([m + 0.03 for m in ref]),
            "L128_ng8": _leg([m + 0.03 + ctrl_off for m in ref])}
    legs["R2build"]["seeds"] = [101, 102, 103]
    kernels = {"K128": _kernel(128, "R2"), "K128build": _kernel(128, "R2build"), "A2_ng8": _kernel(8, "R2build"),
               "L128_ng8": _kernel(128, "R2build", transform=dict(kind="local", n_coarse=8, archetypes=8, n_params=64, in_sample_tv=0.3))}
    for k, dm in dm_global.items():
        legs[f"G_k{k}"] = _leg([m + dm for m in ref])
        # a matrix whose TV distance from the identity K128 is ev_global[k]: identity mixed with uniform
        e = ev_global[k]; R = (1 - e) * np.eye(128) + e * np.full((128, 128), 1 / 128)
        kernels[f"G_k{k}"] = _kernel(128, "R2build", R=R, transform=dict(kind="nmf", k=k, archetypes=k, n_params=256 * k, in_sample_tv=e,
                                                                          rel_frobenius=0.5, rel_change_tail=conv, iters=600, seed=0))
    for f, dm in dm_trunc.items():
        legs[f"T_f{f:g}"] = _leg([m + dm for m in ref])
        kernels[f"T_f{f:g}"] = _kernel(128, "R2build", n_exit=int(1000 * rho_trunc[f]),
                                       transform=dict(kind="truncate", f=f, n_exit_kept=int(1000 * rho_trunc[f]), n_exit_total=1000, archetypes=128, n_params=128 ** 2))
    return dict(ion=ion, n=G2.PREREG["n"][ion], seeds=[1, 2, 3], build_seeds=[101, 102, 103], k_grid=[1, 2, 4, 8, 16, 32], f_grid=[0.1, 0.2, 0.5, 0.9, 0.99, 0.999],
                ng_control=8, ng_fine=128, state="paper4/phase1_benchmarks/P1_t2.json", shell=28, lam_window=[1000.0, 30000.0], n_spec=200,
                legs=legs, kernels=kernels)


@pytest.fixture(autouse=True)
def _equal_band_fractions(monkeypatch):
    monkeypatch.setattr(G2.A.V, "band_fractions", lambda row, leg="R2": {b: 1.0 / 7 for b in BANDS})
    monkeypatch.setattr(G2.A.V, "nu_edges", lambda lo, hi, n: np.geomspace(1e14, 3e15, n + 1))


LOCAL_DM = {2: 0.5, 4: 0.3, 8: 0.2, 16: 0.06, 32: 0.04}                   # K*_local = 16 (a Ce-like ion)
LOCAL_EV = {2: 0.6, 4: 0.5, 8: 0.4, 16: 0.3, 32: 0.2}
ND_DM = {2: 0.05, 4: 0.04, 8: 0.03, 16: 0.02, 32: 0.02}                    # K*_local = 2 (an Nd-like ion)
TRUNC_DM = {0.1: 0.5, 0.2: 0.08, 0.5: 0.05, 0.9: 0.04, 0.99: 0.03, 0.999: 0.03}
TRUNC_RHO = {0.1: 0.02, 0.2: 0.04, 0.5: 0.12, 0.9: 0.36, 0.99: 0.67, 0.999: 0.87}   # Ce II's audited curve


def _case(ion, dm_g, ev_g, dm_t=None, rho_t=None, **kw):
    return _g2(ion, dm_g, dm_t or TRUNC_DM, rho_t or TRUNC_RHO, ev_g, **kw)


GLOBAL_16 = {1: 0.9, 2: 0.6, 4: 0.5, 8: 0.4, 16: 0.08, 32: 0.05}           # K*_global = 16
EV_GLOBAL_BETTER = {1: 0.5, 2: 0.4, 4: 0.3, 8: 0.2, 16: 0.1, 32: 0.05}     # the global family fits the events better


def test_h1_green_on_archetype_counts_alone():
    """Green needs only K*_global >= K*_local on both decisive ions; the
    event-level fit is a diagnostic (the PI's amendment of 2026-09-22)."""
    g1 = {i: _g1(LOCAL_DM, LOCAL_EV) for i in G2.IONS}; g1["60NdII"] = _g1(ND_DM, LOCAL_EV)
    g2 = {i: _case(i, GLOBAL_16, EV_GLOBAL_BETTER) for i in G2.IONS}
    out = G2.readings(g2, g1)
    ce = out["per_ion"]["58CeII"]
    assert ce["k_local"] == 16 and ce["k_global"] == 16 and ce["H1"] == "GREEN"
    assert ce["diagnostic"]["events_vs_observables"] is True                    # global fits events better, local the light
    nd = out["per_ion"]["60NdII"]
    assert nd["k_local"] == 2 and nd["k_global"] == 16 and nd["H1"] == "GREEN"
    assert out["H1"] == "GREEN" and out["decision"] == "WRITE" and not out["ions_gray"]


def test_h1_green_even_when_the_global_family_fits_the_events_worse():
    """The amendment's substance: the diagnostic may be False and Green stands."""
    g1 = {i: _g1(LOCAL_DM, LOCAL_EV) for i in G2.IONS}; g1["60NdII"] = _g1(ND_DM, LOCAL_EV)
    ev_worse = {k: 0.9 for k in GLOBAL_16}                                      # global fits the events WORSE at every k
    g2 = {i: _case(i, GLOBAL_16, ev_worse) for i in G2.IONS}
    out = G2.readings(g2, g1)
    ce = out["per_ion"]["58CeII"]
    assert ce["diagnostic"]["global_fits_events_better"] is False and ce["diagnostic"]["events_vs_observables"] is False
    assert ce["H1"] == "GREEN" and out["H1"] == "GREEN" and out["decision"] == "WRITE"


def test_h1_green_when_the_global_family_never_passes():
    g1 = {i: _g1(LOCAL_DM, LOCAL_EV) for i in G2.IONS}; g1["60NdII"] = _g1(ND_DM, LOCAL_EV)
    never = {k: 0.9 for k in GLOBAL_16}
    g2 = {i: _case(i, never, EV_GLOBAL_BETTER) for i in G2.IONS}
    out = G2.readings(g2, g1)
    assert out["per_ion"]["58CeII"]["k_global"] is None and out["H1"] == "GREEN"


def test_h1_red_when_a_low_rank_law_passes_with_four_times_fewer_archetypes():
    g1 = {i: _g1(LOCAL_DM, LOCAL_EV) for i in G2.IONS}
    dm_g = {1: 0.9, 2: 0.6, 4: 0.05, 8: 0.04, 16: 0.03, 32: 0.02}              # K*_global = 4 <= 16 / 4
    g2 = {i: _case(i, dm_g, {k: 0.1 for k in dm_g}) for i in G2.IONS}
    out = G2.readings(g2, g1)
    assert out["per_ion"]["58CeII"]["H1"] == "RED" and out["H1"] == "RED" and out["decision"] == "REFRAME"


def test_h1_yellow_when_the_single_archetype_null_passes_on_the_nd_like_ion():
    g1 = {i: _g1(ND_DM, LOCAL_EV) for i in G2.IONS}                             # K*_local = 2: Red is out of reach
    dm_g = {1: 0.05, 2: 0.04, 4: 0.03, 8: 0.02, 16: 0.02, 32: 0.02}             # K*_global = 1 < 2
    g2 = {i: _case(i, dm_g, {k: 0.1 for k in dm_g}) for i in G2.IONS}
    out = G2.readings(g2, g1)
    nd = out["per_ion"]["60NdII"]
    assert nd["k_global"] == 1 and nd["H1"] == "YELLOW" and out["H1"] == "YELLOW" and out["decision"] == "PI"


def test_h2_ladder_green_yellow_red_on_the_retained_fraction():
    g1 = {i: _g1(LOCAL_DM, LOCAL_EV) for i in G2.IONS}
    # Green: the smallest passing retained fraction is 0.02 <= 0.10
    out = G2.readings({i: _case(i, GLOBAL_16, EV_GLOBAL_BETTER) for i in G2.IONS}, g1)
    ce = out["per_ion"]["58CeII"]
    assert ce["f_star"] == 0.2 and abs(ce["L_star"] - 0.04) < 1e-12 and ce["H2"] == "GREEN" and out["H2"] == "GREEN"
    # Yellow: nothing passes below 0.12
    dm_t = {0.1: 0.5, 0.2: 0.4, 0.5: 0.08, 0.9: 0.04, 0.99: 0.03, 0.999: 0.03}
    out = G2.readings({i: _case(i, GLOBAL_16, EV_GLOBAL_BETTER, dm_t=dm_t) for i in G2.IONS}, g1)
    assert abs(out["per_ion"]["58CeII"]["L_star"] - 0.12) < 1e-12 and out["H2"] == "YELLOW"
    # Red: the first passing fraction is 0.67 > 0.50
    dm_t = {0.1: 0.5, 0.2: 0.4, 0.5: 0.3, 0.9: 0.2, 0.99: 0.05, 0.999: 0.03}
    out = G2.readings({i: _case(i, GLOBAL_16, EV_GLOBAL_BETTER, dm_t=dm_t) for i in G2.IONS}, g1)
    assert abs(out["per_ion"]["58CeII"]["L_star"] - 0.67) < 1e-12 and out["H2"] == "RED"
    # Red also when no truncation passes at all
    dm_t = {f: 0.9 for f in TRUNC_DM}
    out = G2.readings({i: _case(i, GLOBAL_16, EV_GLOBAL_BETTER, dm_t=dm_t) for i in G2.IONS}, g1)
    assert out["per_ion"]["58CeII"]["L_star"] is None and out["H2"] == "RED"


def test_la_is_a_control_and_does_not_decide_either_reading():
    g1 = {i: _g1(LOCAL_DM, LOCAL_EV) for i in G2.IONS}
    g2 = {i: _case(i, GLOBAL_16, EV_GLOBAL_BETTER) for i in G2.IONS}
    dm_g_bad = {1: 0.02, 2: 0.02, 4: 0.02, 8: 0.02, 16: 0.02, 32: 0.02}         # La alone would read Red
    g2["57LaII"] = _case("57LaII", dm_g_bad, {k: 0.1 for k in dm_g_bad}, dm_t={f: 0.9 for f in TRUNC_DM})
    out = G2.readings(g2, g1)
    assert out["per_ion"]["57LaII"]["H1"] == "RED" and out["per_ion"]["57LaII"]["H2"] == "RED"
    assert out["decisive"] == ["58CeII", "60NdII"] and out["H1"] == "GREEN" and out["H2"] == "GREEN"


def test_gray_on_reference_mismatch_control_and_nmf_convergence():
    g1 = {i: _g1(LOCAL_DM, LOCAL_EV) for i in G2.IONS}
    g2 = {i: _case(i, GLOBAL_16, EV_GLOBAL_BETTER) for i in G2.IONS}
    g2["58CeII"] = _case("58CeII", GLOBAL_16, EV_GLOBAL_BETTER, ref_off=1e-3)
    out = G2.readings(g2, g1)
    assert any(f.startswith("6:") for f in out["per_ion"]["58CeII"]["gray"]) and out["decision"] == "GRAY"
    g2["58CeII"] = _case("58CeII", GLOBAL_16, EV_GLOBAL_BETTER, ctrl_off=0.2)
    out = G2.readings(g2, g1)
    assert any(f.startswith("7:") for f in out["per_ion"]["58CeII"]["gray"])
    g2["58CeII"] = _case("58CeII", GLOBAL_16, EV_GLOBAL_BETTER, conv=1e-2)
    out = G2.readings(g2, g1)
    assert any(f.startswith("8:") for f in out["per_ion"]["58CeII"]["gray"])
    assert all(r["H1"] == "GRAY" and r["H2"] == "GRAY" for i, r in out["per_ion"].items() if r["gray"])


def test_prereg_check():
    row = _case("58CeII", GLOBAL_16, EV_GLOBAL_BETTER); assert G2.check_prereg(row) == []
    row["n"] = 100_000
    with pytest.raises(ValueError):
        G2.check_prereg(row)


from sobolev import atomic_cache as ac                                  # noqa: E402


@pytest.mark.skipif(not (ac.CACHE_DIR / "57LaII.npz").exists(), reason="La II cache not built")
def test_runner_smoke_la(tmp_path):
    R = _load("paperB_run_gate2", "paperB/gate2/run_gate2.py")
    row = R.run("57LaII", n=2000, seeds=(1,), build_seeds=(101,), k_grid=(1,), f_grid=(0.9,), out=tmp_path / "r.json", verbose=False)
    assert set(row["legs"]) == {"R2build", "R2", "A2_ng8", "L128_ng8", "G_k1", "T_f0.9"}
    assert set(row["kernels"]) == {"A2_ng8", "L128_ng8", "G_k1", "T_f0.9", "K128", "K128build"}
    for tag, leg in row["legs"].items():
        assert abs(leg["energy"]["identity_residual"]) < 1e-10, tag
    for tag, k in row["kernels"].items():
        assert k["validate_energy"] < 1e-12, tag
    assert row["kernels"]["G_k1"]["transform"]["kind"] == "nmf" and row["kernels"]["T_f0.9"]["transform"]["n_exit_kept"] < row["kernels"]["K128build"]["n_exit_samples"]
    assert row["kernels"]["L128_ng8"]["ng"] == 128 and row["kernels"]["L128_ng8"]["transform"]["archetypes"] == 8
