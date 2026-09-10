"""Paper IV Gate 1: the committed benchmark states reproduce the pinned
published inputs (values ledgered in data/README.md, 'Paper IV')."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "paper4/phase1_benchmarks"):
    sys.path.insert(0, str(p))

from sobolev.constants import C                 # noqa: E402
from sobolev.source import MSUN                 # noqa: E402
from sobolev import ejecta as ej                # noqa: E402
import build                                    # noqa: E402

DAY = 86400.0
P1_DIR = ROOT / "paper4/phase1_benchmarks"
needs_states = pytest.mark.skipif(not (P1_DIR / "P2_t3.4.json").exists(), reason="benchmark states not built")


def test_p1_reproduces_the_xkn_secular_inputs():
    st = build.build_p1(2.0)
    assert abs(st.mass() / MSUN - 2.64e-2) / 2.64e-2 < 1e-6
    m = st.shell_mass()
    v_rms = np.sqrt(np.sum(m * st.v_mid ** 2) / m.sum())
    assert abs(v_rms / C - 0.06) / 0.06 < 1e-2          # 32 shells of the (1-x^2)^3 profile
    assert st.Y_e == 0.20 and st.meta["s_kb"] == 10.0 and abs(st.meta["tau_ms"] - 17) < 1
    # the density prescription: rho / rho_centre = (1 - x^2)^3 at shell centres, to the discretisation
    x = st.v_mid / (st.meta["v_max_c"] * C)
    pred = (1 - x ** 2) ** 3
    assert np.allclose(st.rho / st.rho[0], pred / pred[0], rtol=0.05, atol=1e-3)
    assert st.check(m_target=2.64e-2 * MSUN)["mass_rel"] < 1e-6
    from sobolev.abundances import Z_OF, lanthanide_fraction
    x_lan = sum(st.X[s][0] for s in st.X if s in Z_OF and Z_OF[s] <= 70)
    assert abs(x_lan - 0.11) < 1e-12                       # Tanaka et al. 2020, Table 1, Ye = 0.20
    assert abs(sum(v[0] for v in st.X.values()) - 1.0) < 1e-12
    assert "Sr" in st.X and "Pt" in st.X and "bulk" not in st.X     # the full Ye-0.21a element list
    assert set(st.f_ion) == {f"{s} II" for s in st.X if s in Z_OF and Z_OF[s] <= 70}
    # the lanthanide PATTERN is Ye-0.21a's own (Dy-rich), not solar
    assert st.X["Dy"][0] > st.X["Ce"][0] > st.X["La"][0]
    assert st.meta["photospheric_shell"] >= 0 and np.all(np.diff(st.T_gas) <= 1e-9)


def test_p2_is_the_gillanders_2026_model_exactly():
    st = build.build_p2()
    assert st.t == 3.4 * DAY
    assert st.v_edges[0] == 0.15 * C and st.v_edges[-1] == 0.35 * C
    assert np.all(st.T_gas == 3200.0) and np.all(st.T_rad == 3200.0)
    rho_in = 4e-15 * (2.0 / 3.4) ** 3 * (0.15 * C / 1.4e9) ** -3
    assert np.isclose(st.meta["rho_at_v_in"], rho_in, rtol=1e-12)
    # shell densities follow v^-3 (volume averages; check at the shell centres)
    slope = np.polyfit(np.log(st.v_mid), np.log(st.rho), 1)[0]
    assert abs(slope + 3.0) < 0.01
    from sobolev.abundances import Z_OF, full_composition
    x_ln = sum(st.X[s][0] for s in st.X if s in Z_OF and Z_OF[s] <= 70)
    assert abs(x_ln - 2.5e-3) < 1e-15
    assert abs(sum(v[0] for v in st.X.values()) - 1.0) < 1e-12
    base = full_composition("gillanders2022_Ye-0.29a")
    assert abs(sum(v for k, v in base.items() if k in Z_OF and Z_OF[k] <= 70) - 4.99e-2) < 2e-4   # Table 2
    assert np.isclose(st.X["Ce"][0] / st.X["Nd"][0], base["Ce"] / base["Nd"])              # pattern kept
    assert np.isclose(st.X["Sn"][0], base["Sn"] * (1 - 2.5e-3) / (1 - sum(v for k, v in base.items() if k in Z_OF and Z_OF[k] <= 70)))
    st.check()
    z = st.local_zone(0)
    assert z["r_core"] == 0.15 * C * st.t and z["r_out"] == 0.35 * C * st.t and z["T_gas"] == 3200.0
    # the mass of the line-forming region is an output: 4 pi rho0 t0^3 v0^3 ln(v_out/v_in)
    m_pred = 4 * np.pi * 4e-15 * (2 * DAY) ** 3 * 1.4e9 ** 3 * np.log(0.35 / 0.15)
    assert abs(st.mass() / m_pred - 1) < 1e-3


def test_dataset_table_is_the_one_ledgered():
    """The committed CSVs carry the dataset's sha256 and reproduce Table 2's X_LN."""
    import hashlib
    from sobolev.abundances import DATA, full_composition, lanthanide_fraction
    for name, x in (("gillanders2022_Ye-0.29a", 4.99e-2), ("gillanders2022_Ye-0.21a", 2.99e-1)):
        txt = (DATA / f"{name}.csv").read_text()
        assert "9f06304bd86b8c2756793e5eabaa0ba0f9585aa2f1520f54b74b4e977f819bd1" in txt   # table sha256
        assert "10.17034/404fbfbe-5f47-42ff-a7d0-12e7c447ebff" in txt
        assert abs(lanthanide_fraction(full_composition(name, renormalise=False)) - x) < 2e-4
    src = ROOT / "data/gillanders2022/composition_profiles_complete.ascii"
    if src.exists():
        assert hashlib.sha256(src.read_bytes()).hexdigest() == "9f06304bd86b8c2756793e5eabaa0ba0f9585aa2f1520f54b74b4e977f819bd1"


def test_robustness_states():
    r1 = build.build_p1(2.0, variant="P1r1"); r2 = build.build_p1(2.0, variant="P1r2")
    assert abs(r1.meta["x_lan_actual"] - 0.2989) < 5e-4 and abs(r2.meta["x_lan_actual"] - 0.11) < 1e-12
    assert "bulk" in r2.X and "bulk" not in r1.X
    r1.check(m_target=2.64e-2 * MSUN); r2.check(m_target=2.64e-2 * MSUN)


@needs_states
def test_committed_states_match_the_builders():
    for t_d in build.EPOCHS_P1:
        on_disk = ej.EjectaState.from_json(P1_DIR / f"P1_t{t_d:g}.json")
        fresh = build.build_p1(t_d)
        assert np.allclose(on_disk.rho, fresh.rho, rtol=1e-12) and np.allclose(on_disk.T_gas, fresh.T_gas, rtol=1e-12)
        on_disk.check(m_target=2.64e-2 * MSUN)
    on_disk = ej.EjectaState.from_json(P1_DIR / "P2_t3.4.json")
    fresh = build.build_p2()
    assert np.allclose(on_disk.rho, fresh.rho, rtol=1e-12)
    on_disk.check()
    d = json.loads((P1_DIR / "P2_t3.4.json").read_text())
    assert d["meta"]["source"].startswith("Gillanders et al. 2026")
