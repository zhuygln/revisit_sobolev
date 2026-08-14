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
    """Read a GSI v2 line file into a DataFrame.

    Expected columns (babystep_plan.md section 6):
        wavelength, log_gf, E_lower, E_upper, method_lower, method_upper

    Left unimplemented until you have the file in hand -- the real format
    should drive this function, not a guess made before downloading it.
    """
    raise NotImplementedError("Phase 0C, Session 3 -- inspect the GSI file first")


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
