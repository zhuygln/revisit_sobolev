"""LTE (Saha) ionization of the lanthanides (Paper IV Phase 7).

Everything before Paper IV forced every lanthanide singly ionised
(`ION_FRAC = 1`; `sobolev/populations.py` says why Saha was kept out of the
transfer comparison). Phase 7 adds the ionization balance as a separate,
controlled step:

    n_{j+1} n_e / n_j = 2 U_{j+1}/U_j (2 pi m_e k T / h^2)^{3/2} exp(-chi_j / k T)

with the ionization energies chi_j from NIST ASD (ledgered in data/README.md)
and the partition functions U_II, U_III from the GSI level lists (the same
truncation the populations use). Two things the data do not carry, both
declared and both varied once in Phase 7's sensitivity table:

  * U_I -- GSI has no neutral-stage level lists. `z1_policy="scale"` takes
    U_I = U_II (the two stages' low-lying level structures are similar in
    the lanthanides); "drop" ignores the neutral stage (f_I = 0). At the
    temperatures of the benchmarks the neutral fraction is 1e-4 at 3400 K
    and a few per cent at 2300 K, so the policy matters late.
  * n_e -- set by the non-lanthanide bulk of the ejecta, represented by ONE
    proxy species (`BulkSpecies`: mass number, chi_I, chi_II, a flat
    partition-function ratio) at the state's bulk mass fraction. Charge
    neutrality n_e = sum_Z sum_j j n_{Z,j} is solved by bisection in log n_e
    over every species, lanthanides included.

Stages without atomic data (I, IV) carry no line opacity but stay in the
mass balance; the transport sees f_II and f_III.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .constants import H, K_B, M_E
from .populations import partition_function

EV = 1.602176634e-12          # erg
M_U = 1.66053906660e-24       # g
SAHA_CONST = (2.0 * np.pi * M_E * K_B / H ** 2) ** 1.5     # (2 pi m_e k / h^2)^{3/2}, cm^-3 K^-3/2

# NIST ASD v5.12 (retrieved 2026-09-09): successive ionization energies, eV
#   (I -> II, II -> III, III -> IV)
CHI_EV = {
    "La": (5.5769, 11.18496, 19.1773), "Ce": (5.5386, 10.956, 20.1974),
    "Pr": (5.4702, 10.631, 21.6237), "Nd": (5.52475, 10.783, 22.09),
    "Pm": (5.58187, 10.938, 22.44), "Sm": (5.643722, 11.078, 23.55),
    "Eu": (5.670385, 11.240, 24.84), "Gd": (6.14980, 12.076, 20.54),
    "Tb": (5.8638, 11.513, 21.82), "Dy": (5.939061, 11.647, 22.89),
    "Ho": (6.0215, 11.781, 22.79), "Er": (6.1077, 11.916, 22.70),
    "Tm": (6.184402, 12.065, 23.66), "Yb": (6.254160, 12.179185, 25.053),
}
STAGES = ("I", "II", "III", "IV")


@dataclass
class BulkSpecies:
    """The non-lanthanide ejecta as one proxy species for charge neutrality.
    Defaults: a second-peak r-process element (A ~ 100, chi_I 6.0 eV,
    chi_II 12.0 eV -- Zr/Sr/Te-like), flat partition-function ratios."""
    A: float = 100.0
    chi_ev: tuple = (6.0, 12.0, 25.0)
    u_ratio: tuple = (1.0, 1.0, 1.0)      # U_{j+1}/U_j for j = I, II, III


def saha_ratio(T, n_e, chi_ev, u_hi, u_lo):
    """n_{j+1} / n_j for one ionization step."""
    return 2.0 * (u_hi / u_lo) * SAHA_CONST * T ** 1.5 * np.exp(-chi_ev * EV / (K_B * T)) / n_e


def stage_fractions(T, n_e, chi_ev, u_ratios):
    """(f_I, ..., f_{n+1}) from n successive ratios. chi_ev: n energies;
    u_ratios: n values of U_{j+1}/U_j."""
    r = [saha_ratio(T, n_e, c, u, 1.0) for c, u in zip(chi_ev, u_ratios)]
    rel = np.concatenate([[1.0], np.cumprod(r)])
    return rel / rel.sum()


def partition_function_gsi(g_lev, E_lev, T):
    return float(partition_function(g_lev, E_lev, T))


def u_ratios_for(element, T, u_II, u_III, z1_policy="scale"):
    """U_{j+1}/U_j for I->II, II->III, III->IV of a lanthanide with the GSI
    partition functions of II and III; U_I by policy; U_IV = U_III
    (the IV stage carries no data and no opacity)."""
    if z1_policy == "scale":
        u_I = u_II
    elif z1_policy == "drop":
        u_I = np.inf            # the neutral stage is suppressed
    else:
        raise ValueError(f"z1_policy must be 'scale' or 'drop', got {z1_policy!r}")
    return (u_II / u_I if np.isfinite(u_I) else np.inf, u_III / u_II, 1.0)


def _fractions(T, n_e, chi, u_ratios):
    if not np.isfinite(u_ratios[0]):          # drop the neutral stage
        f = stage_fractions(T, n_e, chi[1:], u_ratios[1:])
        return np.concatenate([[0.0], f])
    return stage_fractions(T, n_e, chi, u_ratios)


def solve_ionization(T, rho, X, partition, bulk=BulkSpecies(), z1_policy="scale", atomic_mass=None,
                     n_e_bounds=(1e-2, 1e20)):
    """Charge-neutral LTE ionization of a zone.

    T : K; rho : g cm^-3; X : {element: mass fraction} with "bulk" for the
    non-lanthanides; partition : {element: (U_II, U_III)} at T;
    atomic_mass : {element: A}. Returns (n_e, {element: (f_I..f_IV)}).
    """
    from .abundances import ATOMIC_MASS
    am = ATOMIC_MASS if atomic_mass is None else atomic_mass
    species = []
    for el, x in X.items():
        if x <= 0:
            continue
        if el == "bulk":
            n = rho * x / (bulk.A * M_U)
            species.append((el, n, bulk.chi_ev, bulk.u_ratio))
        elif el in CHI_EV and el in partition:
            n = rho * x / (am[el] * M_U)
            u_II, u_III = partition[el]
            species.append((el, n, CHI_EV[el], u_ratios_for(el, T, u_II, u_III, z1_policy)))
        else:
            # a named non-lanthanide element (a full published composition):
            # its own mass, the bulk proxy's ionization energies -- declared
            n = rho * x / (am.get(el, bulk.A) * M_U)
            species.append((el, n, bulk.chi_ev, bulk.u_ratio))

    def charge(log_ne):
        n_e = 10.0 ** log_ne
        tot = 0.0
        for el, n, chi, ur in species:
            f = _fractions(T, n_e, chi, ur)
            tot += n * np.sum(np.arange(f.size) * f)
        return tot - n_e

    lo, hi = np.log10(n_e_bounds)
    if charge(lo) <= 0.0:
        # the gas is neutral to below the lower bound (a cold, thin outer
        # cell): the electron density is at the bound and every fraction is
        # evaluated there -- declared rather than raised
        log_ne = lo
    else:
        log_ne = brentq(charge, lo, hi, xtol=1e-12)
    n_e = 10.0 ** log_ne
    fr = {el: tuple(float(v) for v in _fractions(T, n_e, chi, ur)) for el, n, chi, ur in species}
    return n_e, fr
