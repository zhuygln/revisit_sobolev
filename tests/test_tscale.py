"""§4.43: the MC temperature direction removes the grey term exactly."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "paper3" / "phase12_grid"))
import tscale as T  # noqa: E402
import sensitivity as S  # noqa: E402


def _row(t, mags, L):
    return {"t_d": t, "status": "ok", "n_used": 1000, "T_gas": 5000.0,
            "source": {"T_eff": 5000.0, "L": 1e40, "R_ph": 1e14, "v_ph": 0.1},
            "ref": {"mags": mags, "L_bol": L}}


def test_grey_scaling_gives_zero_direction():
    base = {b: 20.0 + i for i, b in enumerate(S.BANDS)}
    lo = {"rows": [_row(1.0, base, 1e40)]}
    # a pure grey brightening by 2 mag: L x 6.31, every band 2 mag brighter
    hi = {"rows": [_row(1.0, {b: m - 2.0 for b, m in base.items()}, 1e40 * 10 ** 0.8)]}
    d, checks = T.mc_direction({T.SCALES[0]: lo, T.SCALES[1]: hi})
    assert all(abs(v) < 1e-9 for v in d.values())
    assert checks[1.0]["same_L_R_v"] and abs(checks[1.0]["grey_mag"] - 2.0) < 1e-9


def test_colour_change_is_recovered_per_unit_lnT():
    base = {b: 20.0 for b in S.BANDS}
    lo = {"rows": [_row(2.0, base, 1e40)]}
    shifted = dict(base); shifted["g"] -= 0.5; shifted["K"] += 0.25
    hi = {"rows": [_row(2.0, shifted, 1e40)]}
    d, _ = T.mc_direction({T.SCALES[0]: lo, T.SCALES[1]: hi})
    dl = np.log(T.SCALES[1] / T.SCALES[0])
    assert abs(d[("g", 2.0)] + 0.5 / dl) < 1e-9 and abs(d[("K", 2.0)] - 0.25 / dl) < 1e-9
    assert d[("r", 2.0)] == 0.0
