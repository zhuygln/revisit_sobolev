"""Paper B gate G3: the state domain. One coordinate moves at a time from the
reference state theta0 (P1, 2 d, shell 28; G1/G2's), on four axes:

    T  T_gas of the atom           (coordinate log T_gas)
    D  number density n_ion x s    (coordinate log n_ion, i.e. log rho at fixed composition)
    J  source temperature t_core   (coordinate log T_core; the shape of the incident J_nu,
                                    decoupled from T_gas -- Paper III assumed it irrelevant)
    P  the physical trajectory     (P1 at 1, 3, 5 d, shells 29, 27, 26; coordinate log t)

`build_state(axis, value, ion)` returns everything a run needs. The interior
grid point of each axis is bracketed by theta0 and one other state, so the
anchor is always an interpolation endpoint.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "paper2/phase1", ROOT / "paperB/gate1"):
    sys.path.insert(0, str(p))

from sobolev.ejecta import EjectaState                    # noqa: E402
from sobolev.abundances import ATOMIC_MASS                 # noqa: E402
from forest_mc import ForestAtom                            # noqa: E402
from run_gate1 import IONS, TAU_MIN                         # noqa: E402

BENCH = ROOT / "paper4/phase1_benchmarks"
REF_T_D = 2
SHELL_OF_DAY = {1: 29, 2: 28, 3: 27, 5: 26}
T_REF = 3401.0

AXES = {
    "T": dict(coord="logT", grid=[2500.0, 3000.0, 4000.0, 5000.0], interior=3000.0, bracket=(2500.0, T_REF), ref=T_REF),
    "D": dict(coord="logn", grid=[0.1, 0.3, 3.0, 10.0], interior=0.3, bracket=(0.1, 1.0), ref=1.0),
    "J": dict(coord="logTcore", grid=[2500.0, 5000.0, 7000.0], interior=5000.0, bracket=(T_REF, 7000.0), ref=T_REF),
    "P": dict(coord="logt", grid=[1, 3, 5], interior=3, bracket=(2, 5), ref=2),
}


def label(axis, value):
    if axis == "ref":
        return "ref"
    if axis == "D":
        return f"D{value:g}"
    if axis == "P":
        return f"P{int(value)}d"
    return f"{axis}{int(round(value))}"


def ref_zone_and_density(ion, day=REF_T_D):
    st = EjectaState.from_json(BENCH / f"P1_t{day}.json")
    sh = SHELL_OF_DAY[day]
    zone = st.local_zone(sh)
    el = IONS[ion]
    n_ion = float(st.n_ion(el, "II", sh, ATOMIC_MASS[el]))
    return zone, n_ion, sh


def build_state(axis, value, ion, tau_min=TAU_MIN):
    """-> dict(axis, value, label, coord, coord_value, zone, atom, n_ion, T_gas, t_core, day, shell)."""
    if axis == "ref":
        zone, n_ion, sh = ref_zone_and_density(ion)
        T = float(zone["T_gas"]); day = REF_T_D
    elif axis == "T":
        zone, n_ion, sh = ref_zone_and_density(ion)
        T = float(value); zone = dict(zone, T_gas=T); day = REF_T_D
    elif axis == "D":
        zone, n_ion, sh = ref_zone_and_density(ion)
        n_ion = n_ion * float(value); T = float(zone["T_gas"]); day = REF_T_D
    elif axis == "J":
        zone, n_ion, sh = ref_zone_and_density(ion)
        zone = dict(zone, t_core=float(value)); T = float(zone["T_gas"]); day = REF_T_D
    elif axis == "P":
        day = int(value)
        zone, n_ion, sh = ref_zone_and_density(ion, day)
        T = float(zone["T_gas"])
    else:
        raise ValueError(axis)
    atom = ForestAtom.from_cached([(ion, n_ion)], T, zone["t_exp"], tau_min=tau_min)
    coord = None if axis == "ref" else AXES[axis]["coord"]
    cval = dict(T=np.log(T), D=np.log(n_ion), J=np.log(float(zone["t_core"])), P=np.log(float(zone["t_exp"])))
    return dict(axis=axis, value=value, label=label(axis, value), coord=coord,
                coord_value=None if coord is None else float(cval[axis]),
                zone=zone, atom=atom, n_ion=n_ion, T_gas=T, t_core=float(zone["t_core"]), day=day, shell=sh)


def all_states():
    """[(axis, value)] over the domain: the reference and the 14 others."""
    out = [("ref", None)]
    for ax, a in AXES.items():
        out += [(ax, v) for v in a["grid"]]
    return out


def coord_of(axis, value, ion=None):
    """The frozen interpolation coordinate of a state (log of the moved quantity)."""
    if axis == "T":
        return float(np.log(value))
    if axis == "D":
        _, n_ref, _ = ref_zone_and_density(ion)
        return float(np.log(n_ref * value))
    if axis == "J":
        return float(np.log(value))
    if axis == "P":
        return float(np.log(value * 86400.0))
    raise ValueError(axis)
