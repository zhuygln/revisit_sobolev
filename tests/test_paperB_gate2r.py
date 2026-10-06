"""Paper B G2R (the frequency-adjacency ablation): the block coarse-graining
on the fine matrix reproduces G2's local family exactly on the frozen
records, the permuted operator keeps the block structure and the energy
rows, the orderings are reproducible from their seeds, and the readings
follow the preregistered ladder on synthetic cells."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper3", ROOT / "paperB/gate2", ROOT / "paperB/gate2r"):
    sys.path.insert(0, str(p))

from redistribution import RedistributionKernel                      # noqa: E402
import operators as OP                                                 # noqa: E402
import operators_perm as OPP                                           # noqa: E402


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


G2R = _load("paperB_analyse_g2r", "paperB/gate2r/analyse.py")
G2_RECORDS = {ion: ROOT / f"paperB/gate2/gate2_{ion}.json" for ion in ("58CeII", "60NdII")}


@pytest.mark.parametrize("ion", sorted(G2_RECORDS))
def test_physical_ordering_reproduces_the_frozen_g2_local_family(ion):
    p = G2_RECORDS[ion]
    if not p.exists():
        pytest.skip("no G2 record")
    row = json.loads(p.read_text())
    fine = row["kernels"]["K128build"]
    R, E = np.asarray(fine["R"]), np.asarray(fine["E_in"])
    for k in (2, 4, 8, 16, 32):
        frozen = np.asarray(row["kernels"][f"L128_ng{k}"]["R"])
        got = OPP.block_coarsen(R, E, k, None)
        got[np.asarray(fine["counts"]) <= 0] = 0.0                # with_matrix zeroes the rows never populated
        assert np.allclose(got, frozen, atol=1e-12, rtol=0), (ion, k, np.abs(got - frozen).max())


def _events(n=40_000, n_lines=60, seed=0):
    rng = np.random.default_rng(seed)
    lines = np.sort(rng.uniform(2e14, 8e14, n_lines))
    nu_in = rng.uniform(2.1e14, 7.9e14, n)
    centre = nu_in * rng.choice([1.0, 0.85, 0.7], n, p=[0.5, 0.3, 0.2])
    idx = np.clip(np.searchsorted(lines, centre) + rng.integers(-2, 3, n), 0, n_lines - 1)
    nu_out = lines[idx]
    w_in = rng.uniform(0.5, 1.5, n); w_out = w_in * nu_in / nu_out
    return dict(nu_in=nu_in, nu_out=nu_out, w_in=w_in, w_out=w_out)


def test_identity_transform_equals_local_and_a_permutation_keeps_rows_and_blocks():
    ev = _events()
    k = RedistributionKernel.from_branching_mc(ev["nu_in"], ev["nu_out"], ev["w_in"], 16, nu_lo=2e14, nu_hi=8e14, w_out=ev["w_out"])
    ident = OPP.permuted_local(4, np.arange(16), 0)(k, ev)
    loc = OP.local(4)(k, ev)
    assert np.allclose(ident.R, loc.R, atol=1e-12) and ident.metadata["transform"]["physical"] is True
    assert ident.metadata["transform"]["n_params"] == 16 and ident.metadata["transform"]["archetypes"] == 4
    perm = np.random.default_rng(1001).permutation(16)
    scr = OPP.permuted_local(4, perm, 1, 1001)(k, ev)
    assert not np.allclose(scr.R, loc.R)
    # energy rows: every live row of a block carries the block's row sum, and the tables are untouched
    live = ~k.empty_rows
    assert scr.validate_energy() <= loc.validate_energy() + 1e-12
    assert scr.disc_vals is k.disc_vals and np.array_equal(scr.empty_rows, k.empty_rows)
    inv = np.argsort(perm); Rp = scr.R[np.ix_(perm, perm)]
    for I in range(4):
        rows = Rp[I * 4:(I + 1) * 4]
        rows = rows[live[perm][I * 4:(I + 1) * 4]]
        if rows.shape[0] > 1:
            assert np.allclose(rows, rows[0]), "the rows of one block are identical"
    assert scr.metadata["transform"]["perm_seed"] == 1001 and scr.metadata["transform"]["physical"] is False
    with pytest.raises(ValueError):
        OPP.block_coarsen(k.R, OP.energy_in(k, ev), 4, np.arange(15))
    with pytest.raises(ValueError):
        OPP.block_coarsen(k.R, OP.energy_in(k, ev), 5, None)


def test_orderings_are_reproducible_and_distinct():
    a = OPP.orderings(31, 128); b = OPP.orderings(31, 128)
    assert len(a) == 32 and a[0][2] is None and np.array_equal(a[0][1], np.arange(128))
    assert all(np.array_equal(x[1], y[1]) and x[2] == y[2] for x, y in zip(a, b))
    assert a[1][2] == 1001 and a[31][2] == 1031
    assert len({tuple(p) for _, p, _ in a}) == 32
    assert G2R.PREREG["perm_seeds"] == [s for _, _, s in a]


def _cell(k, e_phys, e_scr, ev_phys=0.2, ev_scr=None):
    ev_scr = ev_scr if ev_scr is not None else [0.25] * len(e_scr)
    orders = [dict(perm_id=0, e_joint=e_phys, event=ev_phys, passes=e_phys <= 0.1)]
    orders += [dict(perm_id=i + 1, e_joint=e, event=v, passes=e <= 0.1) for i, (e, v) in enumerate(zip(e_scr, ev_scr))]
    scr = np.array(e_scr)
    return dict(k=k, gray=[], rank=int(1 + np.sum(scr <= e_phys)), p_value=None, e_physical=e_phys,
                scrambled=dict(n=len(e_scr), min=float(scr.min()), median=float(np.median(scr)), max=float(scr.max()), frac_pass=0.0),
                event_physical=ev_phys, event_scrambled=dict(min=0.25, median=0.25, max=0.25), event_rank=1, orderings=orders)


def test_readings_follow_the_preregistered_ladder():
    scr = list(np.linspace(0.2, 0.8, 31))
    win = {k: _cell(k, 0.05, scr) for k in (2, 4, 8, 16, 32)}
    assert G2R.read_ion(win, 16)["reading"] == "GREEN"
    # first at K*_local but first in only three cells: Yellow
    three = dict(win); three[2] = _cell(2, 0.5, scr); three[4] = _cell(4, 0.5, scr)
    assert G2R.read_ion(three, 16)["reading"] == "YELLOW"
    # outside the top quartile at K*_local: Red, whatever the other cells
    red = dict(win); red[16] = _cell(16, 0.5, scr)
    assert G2R.read_ion(red, 16)["rank_at_k_local"] > 8 and G2R.read_ion(red, 16)["reading"] == "RED"
    # a tie counts against the physical ordering
    tie = dict(win); tie[16] = _cell(16, 0.2, scr)
    assert G2R.read_ion(tie, 16)["rank_at_k_local"] == 2 and G2R.read_ion(tie, 16)["reading"] == "YELLOW"
    # a missing block count is gray
    assert G2R.read_ion({k: win[k] for k in (2, 4, 8, 16)}, 16)["reading"] == "GRAY"
    out = G2R.readings({"58CeII": win, "60NdII": win}, {"58CeII": 16, "60NdII": 2})
    assert out["reading"] == "GREEN" and out["decision"] == "ADJACENCY"
    out = G2R.readings({"58CeII": red, "60NdII": win}, {"58CeII": 16, "60NdII": 2})
    assert out["reading"] == "RED" and out["decision"] == "REFRAME"
    out = G2R.readings({"58CeII": win}, {"58CeII": 16, "60NdII": 2})
    assert out["reading"] == "GRAY"


def test_check_prereg_refuses_a_foreign_record():
    row = dict(state="paper4/phase1_benchmarks/P1_t2.json", shell=28, seeds=[1, 2, 3], build_seeds=[101, 102, 103], ng_fine=128,
               n_scrambled=31, perm_seed0=1000, perm_seeds=G2R.PREREG["perm_seeds"], k=16, n=300_000, ion="58CeII")
    assert G2R.check_prereg(row) == []
    bad = dict(row, perm_seeds=[None] + [7 + m for m in range(1, 32)])
    with pytest.raises(ValueError):
        G2R.check_prereg(bad)
    assert G2R.check_prereg(dict(row, k=3), strict=False)
