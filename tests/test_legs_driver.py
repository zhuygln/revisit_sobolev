"""Paper IV: the legs driver runs the ladder end to end on a small zone and
returns the per-leg photometry, energy accounting and differences."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper4/phase2_energy", ROOT / "paper4/phase1_benchmarks"):
    sys.path.insert(0, str(p))

from sobolev import atomic_cache as ac          # noqa: E402
import legs                                     # noqa: E402
import build                                    # noqa: E402

needs = pytest.mark.skipif(not all((ac.CACHE_DIR / f"{i}.npz").exists() for i in ("57LaII", "58CeII")),
                           reason="La II / Ce II caches not built")


@needs
def test_ladder_runs_on_a_two_element_cut_of_p2():
    st = build.build_p2(n_shell=8)
    keep = {"La", "Ce"}
    x_keep = sum(st.X[s][0] for s in keep)
    st.X = {s: v for s, v in st.X.items() if s in keep or s == "bulk"}
    st.X["bulk"] = 1.0 - sum(st.X[s] for s in keep)
    st.f_ion = {k: v for k, v in st.f_ion.items() if k.split()[0] in keep}
    st.check()
    zone = st.local_zone(0)
    atom, n_ion = legs.atom_for_zone(st, 0)
    assert set(n_ion) == {"57LaII", "58CeII"} and atom.n_opacity > 10
    row = legs.run_legs(zone, atom, 3000, legs=("R1", "R1E", "R2", "B2", "A2"), seeds=(1,),
                        ng=8, relativity=None, verbose=False)
    L = row["legs"]
    assert set(L) == {"R1", "R1E", "R2", "B2", "A2"}
    for tag, o in L.items():
        assert abs(o["energy"]["identity_residual"]) < 1e-10, tag
        assert np.isfinite(o["mags"]["g"]) and np.isfinite(o["mags"]["K"])
        assert "dm_vs_R2" in o and "dm_vs_R1" in o
    assert L["R2"]["energy"]["dep_cm_frac"] == 0.0 and L["B2"]["energy"]["dep_cm_frac"] == 0.0
    assert L["R1"]["energy"]["packets"] == "photon" and L["R2"]["energy"]["packets"] == "energy"
    assert all(abs(v) < 1e-12 for v in L["R2"]["dm_vs_R2"].values())
    assert L["R2"]["scale"] == "equilibrium" and L["R1"]["scale"] == "conserving"
    assert row["timing"]["R2"] > 0 and row["n_opacity"] == atom.n_opacity
