"""Format regression for the GSI transition file parser.

Runs against a 20-row excerpt of the real La II file checked into tests/data/,
so the test works without the full (gitignored) dataset.
"""

from pathlib import Path

import numpy as np

from sobolev.atomic_data import load_gsi, nearest_neighbour_velocity_spacing

EXCERPT = Path(__file__).parent / "data" / "57LaII_transitions_calib_excerpt.txt"


def test_load_gsi_excerpt():
    df = load_gsi(EXCERPT)
    assert df.shape == (20, 19)
    assert list(df.columns[:2]) == ["Lower", "E_Lower"]
    assert "WV_Transition" in df.columns and "Log(gf)" in df.columns
    # First data row of the real file: the La II ground-state transition.
    assert np.isclose(df["E_Upper"].iloc[0], 14147.98)
    assert (df["Type"] == "E1").all()
    assert set(df["Method_Lower"]) <= {"uncalib", "shifted", "xmatch"}


def test_spacing_on_excerpt():
    df = load_gsi(EXCERPT)
    dv = nearest_neighbour_velocity_spacing(df["WV_Transition"].to_numpy())
    assert len(dv) == len(df) - 1
    assert (dv >= 0).all()
