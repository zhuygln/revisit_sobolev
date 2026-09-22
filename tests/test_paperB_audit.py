"""Paper B exit-table audit: the pure functions on synthetic events."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


X = _load("paperB_audit", "paperB/audit/exit_tables.py")


def _events(n=50_000, n_lines=60, seed=0):
    rng = np.random.default_rng(seed)
    lines = np.sort(rng.uniform(2e14, 8e14, n_lines))
    p = 1.0 / np.arange(1, n_lines + 1) ** 1.5; p /= p.sum()          # Zipf: a few lines carry most exits
    idx = rng.choice(n_lines, n, p=p)
    return lines, lines[idx], rng.uniform(0.5, 1.5, n)


def test_distinct_exits_aggregates_exactly_by_frequency():
    lines, nu_out, w = _events()
    vals, e = X.distinct_exits(nu_out, w)
    assert vals.size <= lines.size and np.all(np.isin(vals, lines)) and abs(e.sum() - w.sum()) < 1e-9 * w.sum()
    for v, ev in zip(vals[:5], e[:5]):
        assert abs(ev - w[nu_out == v].sum()) < 1e-9 * w.sum()


def test_discovery_curve_is_monotone_and_saturates_at_the_line_count():
    lines, nu_out, _ = _events()
    curve = X.discovery_curve(nu_out, np.random.default_rng(1), n_min=100)
    d = [c["distinct"] for c in curve]
    assert all(b >= a for a, b in zip(d, d[1:])) and d[-1] == np.unique(nu_out).size <= lines.size
    assert curve[-1]["events"] == nu_out.size


def test_concentration_and_per_group_truncation():
    lines, nu_out, w = _events()
    vals, e = X.distinct_exits(nu_out, w)
    c = X.concentration(e)
    assert 1 <= c["0.5"] <= c["0.9"] <= c["0.99"] <= c["0.999"] <= vals.size and c["0.5"] < vals.size
    edges = np.geomspace(2e14, 8e14, 9)
    kept, populated = X.kept_per_group(vals, e, edges, fractions=(0.5, 0.99, 1.0))
    assert kept["1"]["kept"] == vals.size and kept["0.5"]["kept"] <= kept["0.99"]["kept"] <= vals.size
    assert 1 <= populated <= 8 and kept["0.5"]["max_in_group"] >= 1


def test_size_split_isolates_the_tables():
    lines, nu_out, w = _events()
    nu_in = nu_out[::-1]
    k = X.RedistributionKernel.from_branching_mc(nu_in, nu_out, w, 8, w_out=w)
    total, bare = X.size_split(k)
    n_exit = k.disc_vals.size
    assert total > bare and abs((total - bare) - 3 * 8 * n_exit - 8 * (k.disc_off.size)) < 2000   # three float64 arrays + offsets + headers
