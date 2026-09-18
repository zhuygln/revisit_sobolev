"""The transactional results pipeline (PI amendment 4)."""
import json

import numpy as np

from rtedu import results


def test_results_json_roundtrip(tmp_path):
    g = tmp_path / "generated"; r = tmp_path / "results.json"
    results.record("ch01", {"p": 0.123456789, "arr": np.array([1.0, 2.5]), "n": 3, "flag": True, "name": "x"}, generated=g)
    results.record("ch02", {"q": 1e-7}, generated=g)
    merged = results.merge(generated=g, results=r)
    assert merged["ch01"]["p"] == 0.123457 and merged["ch01"]["arr"] == [1.0, 2.5] and merged["ch01"]["n"] == 3
    assert merged["ch01"]["flag"] is True and merged["ch01"]["name"] == "x"
    assert results.load("ch02", results=r)["q"] == 1e-7
    assert list(json.loads(r.read_text()).keys()) == ["ch01", "ch02"]
    # re-recording ch01 and merging leaves no stale key
    results.record("ch01", {"p": 0.5}, generated=g)
    assert "arr" not in results.merge(generated=g, results=r)["ch01"]
