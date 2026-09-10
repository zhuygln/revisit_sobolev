"""Lanthanide abundance patterns for the Paper IV benchmarks.

Paper III's composition was an equal-mass split of one lanthanide fraction
over La II, Ce II, Pr II and Nd II (`four_ion_split`, kept here so the
relabelled X_4Ln of that grid has a definition in code). Paper IV replaces
it by a pattern X_Z = w_Z X_lan over the available lanthanides La-Yb
(Pm and Lu carry no GSI data), with w_Z from a published source that is
recorded in data/README.md before use. Patterns are keyed by name; each
carries its provenance in `PATTERN_SOURCES`.

Pattern tables live in data/abundances/<name>.csv (columns Z, symbol,
w_mass) so that a benchmark's own nucleosynthesis yield can be dropped in
without code changes; `pattern()` normalises whatever it reads to
sum(w) = 1 over the elements the caller keeps.
"""
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1] / "data" / "abundances"

# standard atomic weights (IUPAC conventional values), amu, H-U plus free neutrons.
# r-process ejecta are neutron-rich, so an element's mean mass differs from the
# solar standard weight by a few per cent; the dataset gives element (isotope-
# summed) mass fractions, and n_ion ~ 1/A at that level is irrelevant here.
_SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se "
            "Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy "
            "Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U").split()
_MASSES = (1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990, 24.305, 26.982,
           28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845,
           58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906,
           91.224, 92.906, 95.95, 98.0, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60,
           126.90, 131.29, 132.91, 137.33, 138.905, 140.116, 140.908, 144.242, 145.0, 150.36, 151.964, 157.25,
           158.925, 162.500, 164.930, 167.259, 168.934, 173.045, 174.967, 178.49, 180.95, 183.84, 186.21, 190.23,
           192.22, 195.08, 196.97, 200.59, 204.38, 207.2, 208.98, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0,
           232.04, 231.04, 238.03)
ELEMENT_Z = {s: i + 1 for i, s in enumerate(_SYMBOLS)}
ELEMENT_Z["Neut"] = 0
ATOMIC_MASS = dict(zip(_SYMBOLS, _MASSES))
ATOMIC_MASS["Neut"] = 1.0087
Z_OF = {"La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
        "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71}
SYMBOL_OF = {z: s for s, z in Z_OF.items()}
LANTHANIDES = tuple(s for s, z in sorted(Z_OF.items(), key=lambda kv: kv[1]) if z <= 70)
GSI_ELEMENTS = tuple(s for s in LANTHANIDES if s != "Pm") + ("Pm",)   # Pm II exists, Pm III does not

# ions with GSI level/transition files in the archive (La-Yb II and III, Pm II only)
GSI_IONS = tuple(f"{Z_OF[s]}{s}{st}" for s in LANTHANIDES for st in ("II", "III")
                 if not (s == "Pm" and st == "III"))

PATTERN_SOURCES = {}     # name -> provenance string, filled by register_pattern / the CSV header


def four_ion_split(x_4ln):
    """Paper III's composition: X_4Ln split equally over La, Ce, Pr, Nd."""
    return {s: x_4ln / 4.0 for s in ("La", "Ce", "Pr", "Nd")}


def register_pattern(name, weights, source):
    """A pattern given in code (tests, placeholders). `weights`: {symbol: w}."""
    w = {k: float(v) for k, v in weights.items()}
    PATTERN_SOURCES[name] = source
    _REGISTRY[name] = w
    return w


_REGISTRY = {}


def pattern(name, elements=None):
    """{symbol: w_Z} with sum(w) = 1 over `elements` (default: every element
    the pattern lists that has GSI data, i.e. La-Yb). Reads
    data/abundances/<name>.csv unless the name was registered in code."""
    if name in _REGISTRY:
        raw = dict(_REGISTRY[name])
    else:
        path = DATA / f"{name}.csv"
        if not path.exists():
            raise FileNotFoundError(f"no pattern {name!r}: {path} missing and not registered")
        raw, src = {}, []
        for line in path.read_text().splitlines():
            if line.startswith("#"):
                src.append(line[1:].strip()); continue
            if not line.strip():
                continue
            z, sym, w = line.split(",")[:3]
            raw[sym.strip()] = float(w)
        PATTERN_SOURCES.setdefault(name, " ".join(src))
    keep = LANTHANIDES if elements is None else tuple(elements)
    w = {s: raw.get(s, 0.0) for s in keep}
    tot = sum(w.values())
    if tot <= 0:
        raise ValueError(f"pattern {name!r} has no weight on {keep}")
    return {s: v / tot for s, v in w.items()}


def mass_fractions(x_lan, name, elements=None):
    """X_Z = w_Z x_lan for the pattern; sum(X_Z) == x_lan exactly."""
    w = pattern(name, elements)
    return {s: x_lan * v for s, v in w.items()}


def full_composition(name, renormalise=True):
    """{symbol: mass fraction} for EVERY element of a data/abundances/<name>.csv
    (e.g. a Gillanders et al. 2022 profile), renormalised to sum to 1 (the
    published tables sum to 0.998-1.000)."""
    path = DATA / f"{name}.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    comp, src = {}, []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            src.append(line[1:].strip()); continue
        if not line.strip():
            continue
        z, sym, w = line.split(",")[:3]
        comp[sym.strip()] = float(w)
    PATTERN_SOURCES.setdefault(name, " ".join(src))
    if renormalise:
        tot = sum(comp.values())
        comp = {k: v / tot for k, v in comp.items()}
    return comp


def lanthanide_fraction(comp):
    return sum(v for k, v in comp.items() if k in Z_OF and Z_OF[k] <= 70)


def rescale_lanthanides(comp, x_ln, renormalise=True):
    """Scale the Z = 57-70 mass fractions of a full composition to sum to
    `x_ln` (Gillanders et al. 2026's 20x reduction of Ye-0.29a; Tanaka et al.
    2020's X_lan for a Y_e pattern); with `renormalise` the other elements
    are scaled to 1 - x_ln so the total stays 1 (stated in the ledger)."""
    x0 = lanthanide_fraction(comp)
    if x0 <= 0:
        raise ValueError("composition has no lanthanides to rescale")
    f_ln = x_ln / x0
    rest0 = 1.0 - x0
    f_rest = (1.0 - x_ln) / rest0 if (renormalise and rest0 > 0) else 1.0
    return {k: v * (f_ln if (k in Z_OF and Z_OF[k] <= 70) else f_rest) for k, v in comp.items()}


def number_to_mass(number_weights):
    """Convert number (abundance) weights to mass weights, {symbol: N} ->
    {symbol: N A} (unnormalised)."""
    return {s: n * ATOMIC_MASS[s] for s, n in number_weights.items()}


def as_array(x, elements=LANTHANIDES):
    return np.array([x.get(s, 0.0) for s in elements])
