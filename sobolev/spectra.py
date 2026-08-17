"""Reading SEDONA spectra and reducing them to a band-averaged ratio.

ONE convention, in one place, because getting it wrong has now produced three
separate spurious results:

  * a 7% inflation in the breadth sweep, from a normalization margin that
    straddled SEDONA's final (partial, unreliable) spectrum bin;
  * a spurious ~5-7% "v_D-independent Sobolev floor" in the optical-depth
    sweeps, from normalizing by raw luminosity instead of the continuum
    ratio, which leaves the Planck slope across the band in the answer;
  * and, earlier, the thermal-emission mismatch of Finding F8.

The common shape: SAME-CODE differentials (SEDONA resolved vs SEDONA
expansion) are immune, because both sides carry the same bias and it
cancels. CROSS-CODE comparisons -- an analytic or deterministic leg against
SEDONA -- are not, and are only as good as the shared normalization. Any new
comparison of that kind must go through this module.
"""

import numpy as np

from .constants import C
from .formal_transfer import planck_bnu


def load_spectrum(path):
    """Return (nu, L_nu) from a SEDONA spectrum file, positive entries only."""
    s = np.loadtxt(path, comments="#")
    nu, lum = s[:, 0], s[:, 1]
    good = lum > 0
    return nu[good], lum[good]


def band_ratio(
    path,
    band,
    red_margin,
    r_core,
    t_core,
    return_spectrum=False,
):
    """Band-averaged L/L_continuum for a SEDONA spectrum.

    The spectrum is divided by the ANALYTIC continuum
    L_cont = 4 pi^2 r_core^2 B_nu(T_core) before normalizing, which removes
    the Planck slope across the band. The remaining scalar -- SEDONA's
    lightbulb pours the whole core luminosity into the transport window, so
    absolute L_nu is high by 1/f_window -- is taken out using the mean ratio
    over a line-free red margin.

    Parameters
    ----------
    band : (lam_lo, lam_hi) in Angstrom, the absorbed region to average over.
    red_margin : (lam_lo, lam_hi) in Angstrom, line-free and, importantly,
        clear of the last spectrum bin: SEDONA's final bin is partial and its
        flux collapses, so a margin that touches it depresses the reference
        and inflates every band value.

    Returns the band-averaged ratio, or (ratio, lam, ratio_spectrum) if
    `return_spectrum`.
    """
    nu, lum = load_spectrum(path)
    lam = C / nu * 1e8
    cont = 4.0 * np.pi**2 * r_core**2 * planck_bnu(nu, t_core)

    red = (lam > red_margin[0]) & (lam < red_margin[1])
    if red.sum() < 3:
        raise ValueError(
            f"{path}: only {red.sum()} points in red margin {red_margin}"
        )
    scale = np.mean(lum[red] / cont[red])
    ratio = lum / scale / cont

    avg = band_average(lam, ratio, band)
    if return_spectrum:
        return avg, lam, ratio
    return avg


def band_average(lam, ratio, band):
    """Band-average an already-normalized ratio spectrum (e.g. an analytic
    leg), using the same convention as `band_ratio`.

    Divides by the ACTUAL span integrated rather than the nominal band width:
    grid points rarely land on the band edges, and dividing by the nominal
    width biases the result low by roughly one grid spacing over the band.
    Small (~2e-4 here) but it is free to avoid, and it makes a line-free
    spectrum return exactly 1.
    """
    lam = np.asarray(lam, dtype=float)
    m = (lam > band[0]) & (lam < band[1])
    order = np.argsort(lam[m])
    x = lam[m][order]
    if x.size < 2:
        raise ValueError(f"band {band} contains {x.size} points")
    return float(np.trapezoid(np.asarray(ratio)[m][order], x) / (x[-1] - x[0]))
