"""Synthetic photometry: turning a Monte Carlo packet list into magnitudes.

Every result in this project so far is a *transport residual* -- a band ratio
against a reference run. That is the right quantity for asking whether a closure
is accurate, and the wrong one for asking whether anybody would notice. This
module converts `forest_mc.run_mc` output into absolute L_nu, bolometric
luminosity and AB magnitudes, so the closure error can be quoted in the units an
observer works in.

Three deliberate limitations, all of which must be stated wherever the numbers
appear:

1. **The default bandpasses are top-hats**, not transmission curves. Since
   2026-09-02 `data/filters/` carries the DECam griz and 2MASS JHKs curves from
   the SVO Filter Profile Service, and `Passband` integrates a histogram L_nu
   through them exactly; `load_passbands()` returns them keyed like
   `BANDS_PHOT` so every downstream key is filter-set-agnostic. Results quote
   which set they used (`filter_set`). The quantity this module exists to
   produce is still a *difference between two transport treatments on
   identical bands*, in which the profile shape is second order; an absolute
   magnitude computed here is for scale only and is not photometry.
2. **The source is imposed.** `run_mc` illuminates a line-blanketed shell with a
   blackbody core whose temperature and radius are inputs. There is no
   radioactive heating and no energy equation, so the *shape* of L_bol(t) is
   the core's, not the ejecta's. Only differences between treatments on the
   same source carry information.
3. **The launch window is finite.** L_bol here means "escaping energy within the
   launch window", not the true bolometric luminosity. Keep the window fixed
   across every epoch and leg being compared, or the comparison is meaningless.

Conventions follow what the repo already uses: the continuum normalization
`4 pi^2 r_core^2 B_nu(T)` is `sobolev/spectra.py:64`'s, and `planck_bnu` is
imported from `formal_transfer` rather than redefined.
"""
from pathlib import Path

import numpy as np

from .constants import C, H
from .formal_transfer import planck_bnu

# AB zero point: 3631 Jy = 3631e-23 erg s^-1 cm^-2 Hz^-1
F_AB = 3631.0e-23

# Top-hat bandpasses, (lam_lo, lam_hi) in Angstrom, on standard effective
# wavelengths. The GSI II line lists carry opacity from ~1100 A past 2 um, so
# every one of these sits inside the forest rather than on a bare continuum.
BANDS_PHOT = {
    "g": (4000.0, 5500.0),
    "r": (5500.0, 6900.0),
    "i": (6900.0, 8200.0),
    "z": (8200.0, 9500.0),
    "J": (11000.0, 14000.0),
    "H": (15000.0, 18000.0),
    "K": (20000.0, 24000.0),
}

COLORS = (("g", "r"), ("r", "i"), ("i", "z"), ("i", "J"), ("J", "K"))

# 40 Mpc -- the AT2017gfo distance, used only to put absolute magnitudes on a
# familiar scale. Nothing in this module's conclusions depends on it.
D_40MPC = 40.0 * 3.085677581e24

# Real transmission curves: DECam griz + 2MASS JHKs (the AT2017gfo follow-up
# system), keyed like BANDS_PHOT so `COLORS` and every result key carry over.
FILTER_DIR = Path(__file__).resolve().parents[1] / "data" / "filters"
FILTER_FILES = {"g": "decam_g", "r": "decam_r", "i": "decam_i", "z": "decam_z",
                "J": "2mass_J", "H": "2mass_H", "K": "2mass_Ks"}


class Passband:
    """A transmission curve T(lambda) and its exact integral over a histogram.

    `emergent_lnu` is a histogram on frequency bin `edges`, so the photon-
    counting band average `int f_nu T dnu/nu / int T dnu/nu` is exactly
    `sum_b f_b W_b / sum_b W_b` with `W_b = int_bin T dnu/nu` -- bin-averaged
    weights, not centre sampling (DECam g spans ~19 of the 200 log bins). The
    weights are computed on a 0.2 A resampling of the curve; since
    |dnu/nu| = |dlambda/lambda| the integral is done in wavelength.
    """

    def __init__(self, lam, T, name=""):
        lam, T = np.asarray(lam, float), np.asarray(T, float)
        o = np.argsort(lam)
        self.lam, self.T, self.name = lam[o], np.clip(T[o], 0.0, None), name
        self.lam_lo, self.lam_hi = float(self.lam[0]), float(self.lam[-1])

    @classmethod
    def from_file(cls, path, name=""):
        a = np.loadtxt(path)
        return cls(a[:, 0], a[:, 1], name or Path(path).stem)

    @property
    def lam_eff(self):
        """Photon-weighted effective wavelength, int lam T dlam / int T dlam."""
        return float(np.trapezoid(self.lam * self.T, self.lam) / np.trapezoid(self.T, self.lam))

    def bin_weights(self, edges_hz, dlam=0.2):
        """W_b = int_bin T dnu/nu for each frequency bin; NaN-free, zero outside.

        Raises ValueError if the curve is not fully inside the grid, because a
        clipped filter is a different filter.
        """
        edges_hz = np.asarray(edges_hz, float)
        lam_edges = np.sort(C / edges_hz) * 1e8                # A, ascending
        if self.lam_lo < lam_edges[0] or self.lam_hi > lam_edges[-1]:
            raise ValueError(f"{self.name}: {self.lam_lo:.0f}-{self.lam_hi:.0f} A "
                             f"outside grid {lam_edges[0]:.0f}-{lam_edges[-1]:.0f} A")
        lam = np.arange(self.lam_lo, self.lam_hi + dlam, dlam)
        T = np.interp(lam, self.lam, self.T)
        k = np.digitize(lam, lam_edges) - 1                     # wavelength-bin index
        w_lam = np.bincount(k, weights=T / lam * dlam, minlength=len(lam_edges) - 1)
        return w_lam[::-1]                                      # to ascending-frequency order


def load_passbands(names=tuple(FILTER_FILES), directory=FILTER_DIR):
    """The real filter set, keyed like BANDS_PHOT (the 2MASS Ks curve is 'K')."""
    return {b: Passband.from_file(Path(directory) / f"{FILTER_FILES[b]}.dat", name=FILTER_FILES[b])
            for b in names}


def nu_edges(lam_lo, lam_hi, n_bins):
    """Ascending frequency bin edges spanning [lam_lo, lam_hi] Angstrom."""
    return np.geomspace(C / (lam_hi * 1e-8), C / (lam_lo * 1e-8), n_bins + 1)


def planck_luminosity(nu_lo, nu_hi, r_core, temperature, n=20001):
    """Core luminosity emergent in [nu_lo, nu_hi], erg/s.

    `4 pi^2 r_core^2 int B_nu dnu` -- the same normalization
    `sobolev/spectra.py:34` divides SEDONA spectra by, integrated over the
    window rather than evaluated per bin. Over the full spectrum this reduces
    to `4 pi r^2 sigma T^4`, which is what the test asserts.
    """
    nu = np.geomspace(nu_lo, nu_hi, n)
    return 4.0 * np.pi**2 * r_core**2 * float(np.trapezoid(planck_bnu(nu, temperature), nu))


def _scale(res, l_core_window, core):
    """Synthetic-erg to erg/s conversion for the two inner boundaries.

    "absorbing": the core swallows what returns to it, so `E_inj` synthetic
    ergs are `l_core_window` and the emergent luminosity is the escaped
    fraction of that (every result up to F41).

    "equilibrium": the core re-emits what it absorbs as the same blackbody it
    launches, and packets thermalized in the ejecta (fate 3) are returned to
    it. A re-emitted packet is statistically a fresh launch, so the multi-pass
    emergent spectrum is the single-pass escaped spectrum scaled by the
    geometric series 1/(1 - f_return): `E_inj - E_core - E_abs` synthetic
    ergs are `l_core_window`. Emergent luminosity then equals the core's
    window luminosity minus what the flow deposits, for every closure alike.

    "conserving": additionally, the energy the packets leave in the gas
    (`E_dep_lab`: the level-energy difference of every fluorescence, plus the
    O(v/c) work) is re-radiated with the escaped spectrum's own shape, so
    `E_esc` synthetic ergs are `l_core_window` and every leg emerges with the
    core's window luminosity. This is the inner boundary AND the radiative
    equilibrium a diffusion-powered source implies, in the crudest form that
    needs no gas-temperature iteration. Paper III §4.39 needs it because the
    photon-number-conserving branching of `run_mc` does not conserve energy:
    at S ~ 10^4 a packet interacts ~300 times before it escapes, the
    reference deposits 45% of the injected energy in the gas and returns 75%
    of the packets to the core, while a grouped closure that interacts less
    does neither -- so the absorbing-core "bolometric" error of 1.9 mag is
    the harness's energy bookkeeping, not a property of the closure. Under
    "conserving" the in-window Delta m_bol measures only leakage out of the
    launch window; colours are unchanged by any of the three (each is a grey
    rescaling).
    """
    a = res["accounting"]
    if core == "absorbing":
        return l_core_window / a["E_inj"]
    if core == "equilibrium":
        return l_core_window / (a["E_inj"] - a["E_core"] - a["E_abs"])
    if core == "conserving":
        return l_core_window / a["E_esc"]
    raise ValueError(f"core must be 'absorbing', 'equilibrium' or 'conserving', got {core!r}")


def emergent_lnu(res, edges, l_core_window, core="absorbing"):
    """Absolute emergent L_nu on `edges` (Hz), erg s^-1 Hz^-1.

    `forest_mc.spectrum()` returns escaped/launched, a ratio to the injected
    continuum; nothing in the repo produced an absolute MC spectrum before this.
    The scale comes from the injection, not the escape: `E_inj` synthetic ergs
    correspond to `l_core_window`, so `L_nu = hist(w h nu) * l/E_inj / dnu`.
    Integrating the result recovers `l_core_window * E_esc/E_inj` by
    construction, which the caller may use as L_bol. `core="equilibrium"`
    replaces E_inj by the energy that leaves the core for good -- see `_scale`.

    Packets escaping outside `edges` are dropped, so `edges` should span the
    launch window.
    """
    esc = res["fate"] == 1
    w_esc = res["w"][esc]
    nu_out = res["nu_out"]
    e_bin, _ = np.histogram(nu_out, edges, weights=w_esc * H * nu_out)
    return e_bin * _scale(res, l_core_window, core) / np.diff(edges)


def band_flux_nu(nu_c, l_nu, band, distance_cm=D_40MPC, edges=None):
    """Photon-weighted mean f_nu over a band, erg s^-1 cm^-2 Hz^-1.

    `int f_nu T dnu/nu / int T dnu/nu` -- the photon-counting average a CCD
    measures, which is what AB magnitudes are defined against. A tuple `band`
    is a top-hat (lam_lo, lam_hi) integrated by trapezoid over the bin centres
    (unchanged since F41); a `Passband` needs the histogram `edges` and uses
    the exact bin weights. Returns NaN if the band is not covered by the grid.
    """
    if isinstance(band, Passband):
        if edges is None:
            raise ValueError("a Passband needs the histogram bin edges")
        try:
            w = band.bin_weights(edges)
        except ValueError:
            return np.nan
        f_nu = np.asarray(l_nu, float) / (4.0 * np.pi * distance_cm**2)
        return float(np.sum(f_nu * w) / np.sum(w))
    nu_lo, nu_hi = C / (band[1] * 1e-8), C / (band[0] * 1e-8)
    m = (nu_c >= nu_lo) & (nu_c <= nu_hi)
    if m.sum() < 2:
        return np.nan
    f_nu = np.asarray(l_nu)[m] / (4.0 * np.pi * distance_cm**2)
    x = nu_c[m]
    num = float(np.trapezoid(f_nu / x, x))
    den = float(np.trapezoid(1.0 / x, x))
    return num / den


def ab_magnitude(nu_c, l_nu, band, distance_cm=D_40MPC, edges=None):
    """AB magnitude in one band (top-hat tuple or Passband). NaN if no flux."""
    f = band_flux_nu(nu_c, l_nu, band, distance_cm, edges)
    if not np.isfinite(f) or f <= 0.0:
        return np.nan
    return -2.5 * np.log10(f / F_AB)


def magnitudes(nu_c, l_nu, bands=None, distance_cm=D_40MPC, edges=None):
    """AB magnitude in every band. `nu_c` are bin centres matching `l_nu`;
    `edges` the histogram edges, needed when `bands` holds Passbands."""
    bands = BANDS_PHOT if bands is None else bands
    return {b: ab_magnitude(nu_c, l_nu, w, distance_cm, edges) for b, w in bands.items()}


def colors(mags, pairs=COLORS):
    """Colour indices from a magnitude dict; NaN propagates."""
    return {f"{a}-{b}": mags.get(a, np.nan) - mags.get(b, np.nan) for a, b in pairs}


def delta_mag(mags_approx, mags_ref):
    """Per-band magnitude error of an approximate treatment.

    Negative means the approximation is too BRIGHT. Distance and zero point
    cancel exactly, which is why this -- and not the absolute magnitude -- is
    the quantity the closure comparison is allowed to claim.
    """
    return {b: mags_approx.get(b, np.nan) - mags_ref.get(b, np.nan) for b in mags_ref}


def bolometric(res, l_core_window, core="absorbing"):
    """Escaping luminosity within the launch window, erg/s.

    Not the true bolometric luminosity: energy outside the launch window is
    never injected, and there is no heating source in the ejecta. With
    `core="equilibrium"` it is the core's window luminosity less what the flow
    absorbs adiabatically (`_scale`).
    """
    return res["accounting"]["E_esc"] * _scale(res, l_core_window, core)


def return_fraction(res):
    """Fraction of injected energy that comes back to the core or is
    thermalized in the ejecta -- what an equilibrium core recycles."""
    a = res["accounting"]
    return (a["E_core"] + a["E_abs"]) / a["E_inj"]


def deposited_fraction(res):
    """Fraction of injected energy left in the gas (lab frame): fluorescence
    level-energy differences plus the O(v/c) work -- what "conserving" adds."""
    a = res["accounting"]
    return a["E_dep_lab"] / a["E_inj"]


def bol_delta_mag(l_approx, l_ref):
    """Bolometric magnitude error, -2.5 log10(L_approx / L_ref)."""
    if l_ref <= 0 or l_approx <= 0:
        return np.nan
    return -2.5 * np.log10(l_approx / l_ref)
