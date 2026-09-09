"""Paper IV Phase 1: the two benchmark ejecta states, pinned to published
models (values ledgered in data/README.md, "Paper IV").

P1  the secular component of the xkn radiative-transfer comparison
    (Ricigliano et al. 2024, MNRAS 529, 647, sec. 5.2): M = 2.64e-2 Msun,
    v_rms = 0.06c, Y_e = 0.20, s = 10 k_B/baryon, tau_exp ~ 17 ms, density
    rho ~ (1 - (v/v_max)^2)^3 (eq. 25). A SECULAR component, not
    dynamical-ejecta-like. L(t) and the grey photosphere from
    sobolev.source.SourceModel with the model's mass and v_max and a grey
    kappa stated below; T(v) the Eddington grey run. Lanthanide fraction
    and pattern PROVISIONAL (see the ledger).
P2  the Gillanders et al. 2026 3.4-d AT2017gfo model (MNRAS 548, stag748,
    Table 3): v_in = 0.15c, v_out = 0.35c, rho = rho0 (t0/t)^3 (v/v0)^-3
    with rho0 = 4e-15 g/cm^3, t0 = 2 d, v0 = 14000 km/s, T = 3200 K,
    X_LN = 2.5e-3 (the Ye-0.29a profile with lanthanides reduced 20x).
    The single zone is the published inner-boundary state. Pattern
    PROVISIONAL (solar_r) until the 2022 supplementary list is transcribed.

Gate 1 (pre-declared): check() passes; the committed states reproduce the
pinned inputs (tests/test_benchmarks.py).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

from sobolev.constants import C                                   # noqa: E402
from sobolev.source import MSUN, SourceModel                      # noqa: E402
from sobolev import abundances as ab                              # noqa: E402
from sobolev import ejecta as ej                                  # noqa: E402

DAY = 86400.0
N_SHELL = 32

P1 = dict(name="P1", m_msun=2.64e-2, v_rms_c=0.06, Y_e=0.20, s_kb=10.0, tau_ms=17.0,
          x_lan=0.10, pattern="solar_r", kappa=10.0,
          source="Ricigliano et al. 2024, MNRAS 529, 647, sec. 5.2, eq. 25",
          provisional=["x_lan (0.10, not tabulated by the source)", "pattern (solar_r)",
                       "kappa (10 cm^2/g grey, lanthanide-rich)"])
P2 = dict(name="P2", t_d=3.4, v_in_c=0.15, v_out_c=0.35, rho0=4e-15, t0_d=2.0, v0_kms=14000.0,
          gamma=-3.0, T=3200.0, x_lan=2.5e-3, pattern="solar_r",
          source="Gillanders et al. 2026, MNRAS 548, stag748, Table 3 and sec. 4",
          provisional=["pattern (solar_r; Ye-0.29a's own lanthanide list is in the 2022 supplement)"])
EPOCHS_P1 = (1.0, 2.0, 3.0, 5.0)


def composition(x_lan, pattern, n):
    w = ab.pattern(pattern)
    X = {s: np.full(n, x_lan * v) for s, v in w.items() if v > 0}
    X["bulk"] = np.full(n, 1.0 - x_lan)
    f_ion = {f"{s} II": np.ones(n) for s in X if s != "bulk"}       # Phase 1: fixed II
    return X, f_ion


def build_p1(t_d, n_shell=N_SHELL):
    m = P1["m_msun"] * MSUN
    t = t_d * DAY
    v_max = ej.xkn_vmax_from_vrms(P1["v_rms_c"] * C)
    v_edges, rho = ej.xkn_profile(m, t, v_max, n_shell, v_min_frac=0.02)
    n = rho.size
    X, f_ion = composition(P1["x_lan"], P1["pattern"], n)
    st = ej.EjectaState(t=t, v_edges=v_edges, rho=rho, T_gas=np.ones(n), T_rad=np.ones(n),
                        X=X, f_ion=f_ion, Y_e=P1["Y_e"],
                        meta=dict(P1, v_max_c=float(v_max / C), t_d=t_d, n_shell=n_shell))
    src = SourceModel(P1["m_msun"], float(v_max / C), kappa=P1["kappa"])
    L = float(src.luminosity(t))
    T, s, t_eff = ej.grey_temperature(st, P1["kappa"], L)
    st.T_gas = T; st.T_rad = T
    st.meta.update(L=L, T_eff=t_eff, photospheric_shell=s, v_ph_c=float(st.v_edges[s] / C),
                   rho_ph=float(st.rho[s]), T_ph=float(T[s]))
    st.check(m_target=m)
    return st


def build_p2(n_shell=N_SHELL):
    t = P2["t_d"] * DAY
    v_in, v_out = P2["v_in_c"] * C, P2["v_out_c"] * C
    v0 = P2["v0_kms"] * 1e5
    rho_of_v = lambda v: P2["rho0"] * (P2["t0_d"] * DAY / t) ** 3 * (v / v0) ** P2["gamma"]
    v_edges = np.geomspace(v_in, v_out, n_shell + 1)
    # exact shell densities of the published profile (volume averages), no
    # renormalisation: the mass is an output, not an input, of this model
    rho = np.empty(n_shell)
    for i in range(n_shell):
        vv = np.linspace(v_edges[i], v_edges[i + 1], 401)
        rho[i] = np.trapezoid(rho_of_v(vv) * vv ** 2, vv) / np.trapezoid(vv ** 2, vv)
    X, f_ion = composition(P2["x_lan"], P2["pattern"], n_shell)
    st = ej.EjectaState(t=t, v_edges=v_edges, rho=rho, T_gas=np.full(n_shell, P2["T"]),
                        T_rad=np.full(n_shell, P2["T"]), X=X, f_ion=f_ion, Y_e=0.29,
                        meta=dict(P2, n_shell=n_shell, rho_at_v_in=float(rho_of_v(v_in)),
                                  photospheric_shell=0, mass_msun=None))
    st.meta["mass_msun"] = st.mass() / MSUN
    st.check()
    return st


def main(out_dir=HERE, plot=True):
    out_dir = Path(out_dir)
    written = []
    for t_d in EPOCHS_P1:
        st = build_p1(t_d)
        p = out_dir / f"P1_t{t_d:g}.json"; st.to_json(p); written.append(p)
        if plot:
            st.plot(out_dir / f"gate1_P1_t{t_d:g}.png", kappa=P1["kappa"], atomic_mass=ab.ATOMIC_MASS)
        print(f"P1 t={t_d:g} d: M={st.mass() / MSUN:.4e} Msun, v_max={st.meta['v_max_c']:.4f}c, "
              f"L={st.meta['L']:.3e}, T_eff={st.meta['T_eff']:.0f} K, photosphere shell {st.meta['photospheric_shell']} "
              f"at v={st.meta['v_ph_c']:.4f}c, rho={st.meta['rho_ph']:.3e}, T={st.meta['T_ph']:.0f} K")
    st = build_p2()
    p = out_dir / "P2_t3.4.json"; st.to_json(p); written.append(p)
    if plot:
        st.plot(out_dir / "gate1_P2_t3.4.png", atomic_mass=ab.ATOMIC_MASS)
    print(f"P2 t=3.4 d: M(0.15-0.35c)={st.meta['mass_msun']:.4e} Msun, rho(v_in)={st.meta['rho_at_v_in']:.3e}, "
          f"T={P2['T']} K, X_LN={sum(st.X[s][0] for s in st.X if s != 'bulk'):.3e}")
    return written


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()
    main(plot=not a.no_plot)
