"""Toy photometric bands over the five-level atom's ten lines, by wavelength
at escape (lab frame). Six top-hat bands; a magnitude is -2.5 log10 of the
band's escaped energy (arbitrary zero point, so only differences and
colours mean anything)."""
import numpy as np

from . import C

BANDS = {"g": (400.0, 500.0), "r": (500.0, 700.0), "i": (700.0, 1000.0), "J": (1000.0, 1300.0), "H": (1300.0, 1800.0), "K": (1800.0, 2600.0)}
ORDER = list(BANDS)


def band_of(nu):
    """Band name of each frequency (None if outside every band)."""
    lam_nm = 1e7 * C / np.asarray(nu, float)
    out = np.full(lam_nm.shape, None, dtype=object)
    for b, (lo, hi) in BANDS.items():
        out[(lam_nm >= lo) & (lam_nm < hi)] = b
    return out


def band_fluxes(nu, energy=None):
    """Escaped energy per band."""
    nu = np.asarray(nu, float); energy = np.ones(nu.size) if energy is None else np.asarray(energy, float)
    b = band_of(nu)
    return {name: float(energy[b == name].sum()) for name in ORDER}


def magnitudes(fluxes):
    return {b: (-2.5 * np.log10(f) if f > 0 else np.nan) for b, f in fluxes.items()}


def colours(mags, pairs=(("g", "r"), ("r", "i"), ("i", "J"), ("J", "K"))):
    return {f"{a}-{b}": mags[a] - mags[b] for a, b in pairs}
