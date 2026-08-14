"""GSI v2 lanthanide line list: loading, filtering, and line-crowding statistics.

Phase 0C only needs one ion (La II) and only asks one question: does line
crowding on thermal-width velocity scales actually occur in calibrated data?
See babystep_plan.md section 6.

Nothing here touches transport. Keep it that way -- the transport consequences
are a SEDONA question for a later phase.
"""

import numpy as np

from .constants import C


def load_gsi(path):
    """Read a GSI per-ion transition file into a DataFrame.

    Format (GSI Database for Kilonova Radiative Transfer, Floers et al.,
    Zenodo record 19335084): a free-text preamble, a column-documentation
    block bracketed by dashed separator lines, a whitespace-separated header
    row immediately after the last separator, then one row per transition.

    Columns (verbatim from the file): Lower, E_Lower, J_Lower, P_Lower,
    Config_Lower, LS_Lower, Method_Lower, Upper, E_Upper, J_Upper, P_Upper,
    Config_Upper, LS_Upper, Method_Upper, Type, E_Transition, WV_Transition,
    Log(gf), A. Method_* is uncalib / shifted / xmatch; WV_Transition is in
    Angstrom.
    """
    import pandas as pd

    with open(path) as fh:
        lines = fh.readlines()

    # The header row follows the last dashed separator line.
    separators = [i for i, line in enumerate(lines) if line.lstrip().startswith("---")]
    if not separators:
        raise ValueError(f"{path}: no dashed separator lines -- not a GSI transition file?")
    header_idx = separators[-1] + 1
    names = lines[header_idx].split()

    df = pd.read_csv(path, sep=r"\s+", skiprows=header_idx + 1, names=names)
    if df.shape[1] != len(names):
        raise ValueError(
            f"{path}: {df.shape[1]} data columns but {len(names)} header names"
        )
    return df


def nearest_neighbour_velocity_spacing(wavelengths):
    """Velocity spacing to the next line, in cm/s.

        dv_i = c (lambda_{i+1} - lambda_i) / lambda_i

    Returns an array one shorter than the input. Sorts defensively, because an
    unsorted list silently produces negative spacings.
    """
    lam = np.sort(np.asarray(wavelengths, dtype=float))
    return C * np.diff(lam) / lam[:-1]


def overlap_parameter(dv_spacing, v_doppler):
    """O_i = v_doppler / dv_i. O_i > 1 means the neighbour sits inside one Doppler width."""
    return v_doppler / np.asarray(dv_spacing, dtype=float)
