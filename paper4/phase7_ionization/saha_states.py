"""Paper IV Phase 7: LTE Saha ionization on the benchmark states.

Reads a Phase 1 state, replaces its fixed-II ion fractions by the Saha
solution of sobolev/ionization.py in every shell (partition functions of
II and III from the GSI level lists in the compact cache; the neutral
stage by policy; the non-lanthanide bulk as one proxy species), records
n_e, f_I..f_IV per element and shell, and writes <state>_saha.json. The
transport then blends II and III (`legs.py --stages II,III`).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

from sobolev import ionization as ion                       # noqa: E402
from sobolev.abundances import Z_OF                         # noqa: E402
from sobolev.atomic_cache import load_cached                # noqa: E402
from sobolev.ejecta import EjectaState                      # noqa: E402

STAGES = ("I", "II", "III", "IV")


def partition_functions(elements, T):
    out = {}
    from sobolev.abundances import GSI_IONS
    for el in elements:
        u = []
        for st in ("II", "III"):
            name = f"{Z_OF[el]}{el}{st}"
            if name not in GSI_IONS:            # Pm III: no GSI data -> U_III = U_II, declared
                u.append(u[-1]); continue
            d = load_cached(name)
            u.append(ion.partition_function_gsi(d["g_lev"], d["E_lev"], T))
        out[el] = tuple(u)
    return out


def ionize(state, z1_policy="scale", bulk=ion.BulkSpecies()):
    # the lanthanides carry the GSI partition functions and the line opacity;
    # every other named element enters charge neutrality with the proxy
    # energies inside solve_ionization
    elements = [el for el in state.X if el in Z_OF and Z_OF[el] <= 70]
    n = state.n_shell
    f_ion = {f"{el} {st}": np.zeros(n) for el in elements for st in STAGES}
    n_e = np.zeros(n)
    for s in range(n):
        T = float(state.T_gas[s])
        part = partition_functions(elements, T)
        X = {el: float(state.X[el][s]) for el in state.X}
        ne, fr = ion.solve_ionization(T, float(state.rho[s]), X, part, bulk=bulk, z1_policy=z1_policy)
        n_e[s] = ne
        for el in elements:
            for j, st in enumerate(STAGES):
                f_ion[f"{el} {st}"][s] = fr[el][j]
    state.f_ion = f_ion
    state.n_e = n_e
    state.meta = dict(state.meta, ionization="saha_lte", z1_policy=z1_policy,
                      bulk=dict(A=bulk.A, chi_ev=list(bulk.chi_ev), u_ratio=list(bulk.u_ratio)))
    state.check()
    return state


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("states", nargs="+")
    ap.add_argument("--z1-policy", default="scale")
    a = ap.parse_args()
    for path in a.states:
        st = EjectaState.from_json(path)
        ionize(st, a.z1_policy)
        s = st.meta.get("photospheric_shell", 0)
        out = Path(path).with_name(Path(path).stem + "_saha.json")
        st.to_json(out)
        lan = [el for el in st.X if el in Z_OF and Z_OF[el] <= 70]
        f2 = np.mean([st.f_ion[f"{el} II"][s] for el in lan])
        f3 = np.mean([st.f_ion[f"{el} III"][s] for el in lan])
        f1 = np.mean([st.f_ion[f"{el} I"][s] for el in lan])
        print(f"{st.meta.get('name')} t={st.t / 86400:g} d shell {s}: T={st.T_gas[s]:.0f} K n_e={st.n_e[s]:.3e} "
              f"<f_I>={f1:.3f} <f_II>={f2:.3f} <f_III>={f3:.3f}  Ce: {[round(x, 3) for x in (st.f_ion[f'Ce {q}'][s] for q in STAGES)]} -> {out.name}")


if __name__ == "__main__":
    main()
