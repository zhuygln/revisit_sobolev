"""Paper IV WP1: the compact per-ion cache reproduces `load_gsi` and
`ForestAtom.from_gsi` exactly, and streams from the zip archive."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper2/phase1"):
    sys.path.insert(0, str(p))

from sobolev import atomic_cache as ac                     # noqa: E402
from sobolev.atomic_data import load_gsi                   # noqa: E402
from sobolev.constants import C                            # noqa: E402
from sobolev.populations import parse_j, statistical_weight  # noqa: E402
from forest_mc import ForestAtom, run_mc                   # noqa: E402

DATA = ROOT / "data"
needs_la = pytest.mark.skipif(not (DATA / "57LaII_transitions_calib.txt").exists(),
                              reason="GSI La II data not present")
needs_zip = pytest.mark.skipif(not ac.TRANSITIONS_ZIP.exists(), reason="GSI archives not present")


def test_excerpt_round_trips_through_the_cache(tmp_path):
    """The committed 20-row excerpts: every cached array equals load_gsi's."""
    exc = ROOT / "tests/data"
    (tmp_path / "x").mkdir()
    for kind in ("levels", "transitions"):
        (tmp_path / "x" / f"57LaII_{kind}_calib.txt").write_bytes((exc / f"57LaII_{kind}_calib_excerpt.txt").read_bytes())
    path = ac.build_cache("57LaII", tmp_path / "x", tmp_path / "x", out_dir=tmp_path / "cache")
    d = ac.load_cached("57LaII", cache_dir=tmp_path / "cache", build=False)
    tr = load_gsi(tmp_path / "x/57LaII_transitions_calib.txt"); lev = load_gsi(tmp_path / "x/57LaII_levels_calib.txt")
    assert d["n_lines"] == len(tr) and d["n_levels"] == len(lev)
    assert np.array_equal(d["lower"], tr["Lower"].to_numpy()) and np.array_equal(d["upper"], tr["Upper"].to_numpy())
    assert np.array_equal(d["nu0"], C / (tr["WV_Transition"].to_numpy() * 1e-8))
    g_l = statistical_weight(tr["J_Lower"].to_numpy())
    assert np.array_equal(d["f_lu"], 10 ** tr["Log(gf)"].to_numpy() / g_l)
    assert np.array_equal(d["A"], tr["A"].to_numpy())
    assert np.array_equal(d["E_lev"], lev["Energy"].to_numpy(dtype=float))
    assert np.array_equal(d["g_lev"], statistical_weight(parse_j(lev["J"].to_numpy())).astype(np.float32))
    assert path.exists() and ac.available(tmp_path / "cache") == ["57LaII"]


@needs_la
def test_la_ii_atom_from_cache_is_bit_identical_to_from_gsi(tmp_path):
    n_ion = float(np.load(ROOT / "experiments/laII_forest/forest_lines.npz")["n_ion"])
    ac.build_cache("57LaII", DATA, DATA, out_dir=tmp_path)
    a = ForestAtom.from_gsi(DATA / "57LaII_levels_calib.txt", DATA / "57LaII_transitions_calib.txt",
                            3000.0, n_ion, 86400.0, tau_min=1e-3)
    b = ForestAtom.from_cached([("57LaII", n_ion)], 3000.0, 86400.0, tau_min=1e-3, cache_dir=tmp_path)
    for k in ("nu0_all", "tau_all", "beta_all", "emis_w", "A_all", "op_idx", "op_tau", "op_nu",
              "br_perm", "br_cum_A", "br_cum_Ab", "level_energy_cm"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k
    lo, hi = a.op_nu.min() * 0.995, a.op_nu.max() * 1.005
    ra = run_mc(a, 8.64e12, 2.592e13, 86400.0, lo, hi, 5000, "sobolev_branch", seed=5, t_core=6000.0)
    rb = run_mc(b, 8.64e12, 2.592e13, 86400.0, lo, hi, 5000, "sobolev_branch", seed=5, t_core=6000.0)
    assert np.array_equal(ra["nu_out_all"], rb["nu_out_all"], equal_nan=True)
    assert np.array_equal(ra["fate"], rb["fate"])


@needs_zip
def test_small_ion_streams_from_the_archive(tmp_path):
    """Yb III (0.3 MB in the archive) straight from the zip members."""
    path = ac.build_cache("70YbIII", ac.LEVELS_ZIP, ac.TRANSITIONS_ZIP, out_dir=tmp_path)
    d = ac.load_cached("70YbIII", cache_dir=tmp_path, build=False)
    assert d["n_lines"] > 100 and d["n_levels"] > 10
    assert np.all(d["upper"] < d["n_levels"]) and np.all(d["lower"] < d["n_levels"])
    assert np.all(d["E_lev"][d["upper"]] > d["E_lev"][d["lower"]])
    assert np.all(d["nu0"] > 0) and np.all(d["A"] >= 0) and np.all(d["g_lev"] >= 1)
    assert path.stat().st_size < 50 * d["n_lines"] + 20000
