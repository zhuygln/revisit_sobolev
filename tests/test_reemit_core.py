"""Paper IV 2A.3: the re-emitting core validates the exact normalisation.

An opaque core that relaunches every returning packet with the launch
spectrum makes each relaunch an independent draw, so the multi-pass
emergent spectrum is the single-pass one times 1/(1 - f_return):
`photometry._scale(core="equilibrium")`, which `energy_balance.scale_exact`
names. The explicit relaunch (`run_mc(core="reemit")`) must agree with it
per band; that is the whole point of the mode, which never runs in
production.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1"):
    sys.path.insert(0, str(p))

from sobolev.constants import C                        # noqa: E402
from sobolev import photometry                         # noqa: E402
from sobolev import energy_balance as eb               # noqa: E402
from forest_mc import run_mc                           # noqa: E402
import test_forest_mc as tfm                           # noqa: E402
from test_macroatom import cascade_atom                # noqa: E402

T_EXP = tfm.T_EXP
T_CORE = 6000.0


@pytest.mark.parametrize("relativity", [None, "worldline"])
def test_relaunch_reproduces_the_geometric_series_per_band(relativity):
    fa = cascade_atom(2.0, 1.0, 1.0)
    ct = C * T_EXP
    # a thin shell so that a fair fraction of the scattered packets returns
    r_core = tfm.R_CORE if relativity is None else ct / 30.0
    r_out = 1.6 * r_core
    lo = tfm.NU_13 / (1.0 - 1.2 * r_core / ct)
    hi = tfm.NU_13 / (1.0 - 0.8 * np.sqrt(r_out**2 - r_core**2) / ct)
    kw = dict(seed=21, packets="energy", t_core=T_CORE, launch_weight="energy",
              relativity=relativity)
    n = 300000
    single = run_mc(fa, r_core, r_out, T_EXP, lo, hi, n, "sobolev_dmacro", **kw)
    multi = run_mc(fa, r_core, r_out, T_EXP, lo, hi, n, "sobolev_dmacro", core="reemit",
                   core_max_passes=400, **kw)
    a1, a2 = single["accounting"], multi["accounting"]
    f_ret = a1["E_core"] / a1["E_inj"]
    assert 0.1 < f_ret < 0.9, f_ret
    assert a2["E_core"] / a2["E_inj"] < 1e-3                 # (almost) nothing left uncycled
    assert abs(a2["identity_residual"]) < 1e-12
    assert multi["n_core_passes_total"] > 0.5 * n * f_ret / (1 - f_ret)
    L = 1.0
    # bands: the pump window and the two fluorescence lines
    edges = np.array([tfm.NU_32 * 0.9 * (tfm.NU_13 - tfm.NU_32) / tfm.NU_32, tfm.NU_32 * 0.9,
                      tfm.NU_32 * 1.1, lo, hi])
    edges = np.sort(edges)
    l_single = eb.band_luminosities(single, edges, L, "equilibrium")
    l_multi = eb.band_luminosities(multi, edges, L, "absorbing")
    counts, _ = np.histogram(multi["nu_out"], edges)
    for k in range(edges.size - 1):
        if counts[k] < 2000:
            continue
        rel = 4.0 / np.sqrt(counts[k]) + 0.01          # 4 sigma plus the plan's 1 % criterion
        assert abs(l_multi[k] / l_single[k] - 1.0) < rel, (k, l_multi[k] / l_single[k], rel)
    # the exact scale is the equilibrium one; what Paper III's grey "conserving"
    # scale would still add on top of it is exactly the adiabatic work -- a
    # physical loss, not a renormalisation -- and it is O(beta)
    assert eb.scale_exact(single, L) == photometry._scale(single, L, "equilibrium")
    assert np.isclose(eb.renorm_ratio(multi), a2["W"] / a2["E_esc"], rtol=1e-9)
    assert abs(a2["W"]) / a2["E_inj"] < 3 * r_out / ct


def test_pass_cap_absorbs_and_keeps_the_identity():
    fa = cascade_atom(2.0, 1.0, 1.0)
    lo, hi = tfm.pump_band()
    res = run_mc(fa, tfm.R_CORE, tfm.R_OUT, T_EXP, lo, hi, 20000, "sobolev_dmacro", seed=3,
                 packets="energy", t_core=T_CORE, launch_weight="energy", core="reemit",
                 core_max_passes=1)
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12
    assert a["E_core"] > 0                                   # the cap bit
    assert np.max(res["n_core_passes"]) == 1


def test_reemit_requires_energy_packets_and_a_planck_core():
    fa = cascade_atom(2.0, 1.0, 1.0)
    lo, hi = tfm.pump_band()
    with pytest.raises(ValueError):
        run_mc(fa, tfm.R_CORE, tfm.R_OUT, T_EXP, lo, hi, 100, "sobolev_branch", core="reemit", t_core=T_CORE)
    with pytest.raises(ValueError):
        run_mc(fa, tfm.R_CORE, tfm.R_OUT, T_EXP, lo, hi, 100, "sobolev_dmacro", core="reemit",
               packets="energy")
