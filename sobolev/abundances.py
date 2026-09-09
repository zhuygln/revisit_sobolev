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

# standard atomic weights (IUPAC 2013 conventional values), amu
ATOMIC_MASS = {"La": 138.905, "Ce": 140.116, "Pr": 140.908, "Nd": 144.242, "Pm": 145.0,
               "Sm": 150.36, "Eu": 151.964, "Gd": 157.25, "Tb": 158.925, "Dy": 162.500,
               "Ho": 164.930, "Er": 167.259, "Tm": 168.934, "Yb": 173.045, "Lu": 174.967}
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


def number_to_mass(number_weights):
    """Convert number (abundance) weights to mass weights, {symbol: N} ->
    {symbol: N A} (unnormalised)."""
    return {s: n * ATOMIC_MASS[s] for s, n in number_weights.items()}


def as_array(x, elements=LANTHANIDES):
    return np.array([x.get(s, 0.0) for s in elements])
