"""core="reemit_cm": the re-emitting inner boundary thermalises and re-emits
in its own comoving frame (isotropic there, comoving Planck at t_core,
aberrated to the lab) and books the lab-energy change as boundary work, so
the energy identity still closes to roundoff; the lab-frame variant
("reemit") keeps the lab energy and books nothing. Packets are launched in
the volume so that half of them reach the boundary."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1"):
    sys.path.insert(0, str(p))
from sobolev.constants import C                                  # noqa: E402
from forest_mc import run_mc                                     # noqa: E402
from test_macroatom import T_EXP                                 # noqa: E402
import test_forest_mc as tfm                                     # noqa: E402
import test_zoned_run_mc as tz                                   # noqa: E402


def _run(core, seed=3):
    fa, a = tz.forest()
    r_core = 0.6 * tfm.R_OUT
    r_edges = np.array([r_core, 0.75 * tfm.R_OUT, tfm.R_OUT])
    z2 = tz.zoned(a, 2, r_edges, scale=[1.0, 1.0]); z2.rho = np.array([1.0, 1.0])
    lo, hi = tfm.pump_band()
    return run_mc(z2, r_core, tfm.R_OUT, T_EXP, lo, hi, 20000, "sobolev_dmacro", seed=seed, packets="energy",
                  t_core=6000.0, launch_weight="energy", relativity="worldline", core=core, core_max_passes=50,
                  launch="volume")


def test_reemit_cm_closes_the_identity_and_does_boundary_work():
    a = _run("reemit"); b = _run("reemit_cm")
    for res in (a, b):
        assert abs(res["accounting"]["identity_residual"]) < 1e-10
        assert res["n_core_passes_total"] > 1000                # the boundary is exercised
    beta = 0.6 * tfm.R_OUT / (C * T_EXP)
    assert 0.005 < beta < 0.3
    # the boundary recedes from the interior and advances on the packets that
    # come back to it: in its frame they arrive blueshifted, leave isotropic
    # and are forward-beamed to the lab, so the radiation gains lab energy
    # (negative work) and more of it escapes than with the lab-frame variant
    assert b["accounting"]["W"] < a["accounting"]["W"]
    assert b["accounting"]["E_esc"] > a["accounting"]["E_esc"]
    # per pass the gain is O(beta): bounded by (1 + beta)^2 per pass on average
    gain = b["accounting"]["E_esc"] / a["accounting"]["E_esc"]
    passes = b["n_core_passes_total"] / b["n_packets"]
    assert 1.0 < gain < (1.0 + beta) ** (2.0 * passes + 2.0)
