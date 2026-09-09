"""Compact per-ion cache of the GSI line lists (Paper IV WP1).

`atomic_data.load_gsi` reads a whole transition file into a pandas frame with
seven string columns: 565 bytes per line on Ce II (measured), 16 GB for the
27 La-Yb ions of the archive -- unbuildable on a 24 GB machine. The transport
needs six numbers per line and three per level. This module streams a file
(from the extracted text or straight from the zip member, never writing the
5.6 GB of text) in chunks, keeping only those, and caches them as an `.npz`
at 36 bytes per line:

    nu0    f8   rest frequency, Hz, from WV_Transition (Angstrom)
    f_lu   f8   absorption oscillator strength 10^log(gf) / g_lower
    A      f8   Einstein A, s^-1
    lower  i4   lower level index (into the level table)
    upper  i4   upper level index
    E_lev  f8   level energies, cm^-1
    g_lev  f4   level statistical weights 2J + 1

f_lu and A stay f8 so that an atom built from the cache is bit-identical to
one built by `from_gsi` (the test pins it); the memory of a blend is set by
the atom's own tables, not by these. Everything else `load_gsi` carries
(configurations, LS terms, calibration methods) is diagnostic and stays in
the text files for `load_gsi`.

The cache is built once per ion (`build_cache`) and read with
`load_cached`; `ForestAtom.from_cached` builds atoms and blends from it
with exactly the populations and optical depths `from_gsi` would give
(pinned by tests/test_atomic_cache.py to 1e-12 on La II).
"""
import io
import zipfile
from pathlib import Path

import numpy as np

from .constants import C
from .populations import parse_j, statistical_weight

DATA = Path(__file__).resolve().parents[1] / "data"
CACHE_DIR = DATA / "cache"
LEVELS_ZIP = DATA / "GSI_lanthanides_calibrated_levels.zip"
TRANSITIONS_ZIP = DATA / "GSI_lanthanides_calibrated_transitions.zip"
TR_COLS = ["Lower", "J_Lower", "Upper", "WV_Transition", "Log(gf)", "A"]
LEV_COLS = ["Index", "Energy", "J"]


def _header_index(fh):
    """(index of the header row, header names) of a GSI file: the row after
    the last dashed separator, exactly as `load_gsi` finds it."""
    last_sep, cand, cand_idx = None, None, None
    for i, raw in enumerate(fh):
        line = raw.decode() if isinstance(raw, bytes) else raw
        if line.lstrip().startswith("---"):
            last_sep = i
            continue
        if last_sep is not None and i == last_sep + 1:
            cand, cand_idx = line.split(), i        # the row after a separator
            continue
        tok = line.split()
        if cand is not None and tok and tok[0].lstrip("-").isdigit():
            break                                   # first data row: the last candidate is the header
    if cand is None:
        raise ValueError("no dashed separator / header row -- not a GSI file?")
    return cand_idx, cand


def _open(source, member):
    """A text file handle: `source` is a directory of extracted files, a
    single file path, or a zip archive holding `member`."""
    source = Path(source)
    if source.suffix == ".zip":
        zf = zipfile.ZipFile(source)
        names = [n for n in zf.namelist() if n.endswith(member) and "/._" not in n and not n.startswith("._")]
        if not names:
            raise FileNotFoundError(f"{member} not in {source}")
        return io.TextIOWrapper(zf.open(names[0]), encoding="utf-8")
    path = source / member if source.is_dir() else source
    return open(path, encoding="utf-8")


def _read_columns(source, member, usecols, dtypes, chunksize=500_000):
    import pandas as pd
    with _open(source, member) as fh:
        header_idx, names = _header_index(fh)
    missing = [c for c in usecols if c not in names]
    if missing:
        raise ValueError(f"{member}: columns {missing} not in header {names}")
    parts = {c: [] for c in usecols}
    with _open(source, member) as fh:
        reader = pd.read_csv(fh, sep=r"\s+", skiprows=header_idx + 1, names=names,
                             usecols=usecols, dtype=dtypes, chunksize=chunksize,
                             engine="c")
        for chunk in reader:
            for c in usecols:
                parts[c].append(chunk[c].to_numpy())
    return {c: (np.concatenate(v) if v else np.empty(0)) for c, v in parts.items()}


def build_cache(ion, levels_source=None, transitions_source=None, out_dir=CACHE_DIR):
    """Stream one ion's level and transition files into `out_dir/<ion>.npz`.

    ion : e.g. "58CeII". Sources default to the extracted text in data/ if
    present, else the zip archives. Returns the cache path."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    lev_member, tr_member = f"{ion}_levels_calib.txt", f"{ion}_transitions_calib.txt"
    if levels_source is None:
        levels_source = DATA if (DATA / lev_member).exists() else LEVELS_ZIP
    if transitions_source is None:
        transitions_source = DATA if (DATA / tr_member).exists() else TRANSITIONS_ZIP
    lev = _read_columns(levels_source, lev_member, LEV_COLS,
                        {"Index": np.int64, "Energy": np.float64, "J": str})
    idx = lev["Index"]
    if not np.array_equal(idx, np.arange(idx.size)):
        raise ValueError(f"{ion}: level Index is not 0..n-1 (got {idx[:5]}...)")
    g_lev = statistical_weight(np.array([parse_j(j) for j in lev["J"]], float)).astype(np.float32)
    tr = _read_columns(transitions_source, tr_member, TR_COLS,
                       {"Lower": np.int64, "J_Lower": str, "Upper": np.int64,
                        "WV_Transition": np.float64, "Log(gf)": np.float64, "A": np.float64})
    g_l = statistical_weight(np.array([parse_j(j) for j in tr["J_Lower"]], float))
    f_lu = 10 ** tr["Log(gf)"] / g_l          # the expression from_gsi uses, verbatim
    nu0 = C / (tr["WV_Transition"] * 1e-8)
    path = out_dir / f"{ion}.npz"
    np.savez(path, ion=ion, nu0=nu0, f_lu=f_lu, A=tr["A"].astype(np.float64),
             lower=tr["Lower"].astype(np.int32), upper=tr["Upper"].astype(np.int32),
             E_lev=lev["Energy"].astype(np.float64), g_lev=g_lev,
             n_lines=np.int64(nu0.size), n_levels=np.int64(idx.size))
    return path


def load_cached(ion, cache_dir=CACHE_DIR, build=True):
    """The cached arrays of one ion as a dict (builds the cache if absent
    and `build`)."""
    path = Path(cache_dir) / f"{ion}.npz"
    if not path.exists():
        if not build:
            raise FileNotFoundError(path)
        build_cache(ion, out_dir=cache_dir)
    d = np.load(path)
    return {k: d[k] for k in ("nu0", "f_lu", "A", "lower", "upper", "E_lev", "g_lev")} | \
        dict(ion=str(d["ion"]), n_lines=int(d["n_lines"]), n_levels=int(d["n_levels"]))


def available(cache_dir=CACHE_DIR):
    return sorted(p.stem for p in Path(cache_dir).glob("*.npz")) if Path(cache_dir).exists() else []
