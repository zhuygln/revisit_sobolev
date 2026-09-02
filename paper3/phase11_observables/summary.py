"""Reduce §10's scorecards to the tables §4.37 quotes, with the faint-band mask.

Two things this does that reading the JSON naively does not.

**It masks bands the reference barely emits in.** The ejecta cool, and a band
whose reference flux has fallen to a part in 10^5 of the bolometric is measured
from a handful of packets: at 12 d the blue-kilonova run puts 2.5e-6 of its
luminosity in g, and the "0.75 mag" error there is a few photons, not physics.
Any band carrying less than `FRAC_MIN` of the reference L_bol is dropped. The
mask is computed from the stored reference spectrum, so it needs no re-run.

**It reports the `sobolev_group` leg as an empirical noise floor.** That leg
approximates only the redistribution, which §4.35/§4.37 measure at <= 0.013 mag
wherever a run is well sampled; a run in which it reads 0.3 mag is telling you
its own resolution, and nothing smaller in it means anything.

**It can re-photometer the stored spectra through real filters.** Every row
stores the 200-bin L_nu of the reference and of each leg, so `--filters real`
recomputes mags, colours, dm and dcolor through DECam griz + 2MASS JHKs
(`photometry.load_passbands`) with no MC rerun; `--filters tophat` (default)
reads the stored top-hat values. Gate 1 of the 2026-09-02 plan is the
comparison of the two (`gate1.py`).

Usage: python summary.py [--frac 0.01] [--filters tophat|real] [files...]
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from sobolev import photometry as phot
from sobolev.constants import C

FRAC_MIN = 0.01          # a band must carry >= 1% of the reference L_bol
LEGS = ("A_redist", "B_opacity", "C_both", "C_binned")
COLS = ("g-r", "r-i", "i-z", "i-J", "J-K")
FILTER_SETS = ("tophat", "real")


def rephotometer(d, filter_set="real"):
    """Recompute every row's mags / colours / dm / dcolor from the stored L_nu.

    Mutates and returns `d`; sets `d["filter_set"]`. With "tophat" this
    reproduces the stored values (a consistency check), with "real" it uses
    the transmission curves. Bolometric quantities are untouched.
    """
    edges = phot.nu_edges(*d["lam_window"], d["n_spec"])
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    dist = d.get("distance_cm", phot.D_40MPC)
    bands = phot.load_passbands() if filter_set == "real" else phot.BANDS_PHOT
    for r in d["rows"]:
        if "skipped" in r:
            continue
        ref = r["ref"]
        ref["mags"] = phot.magnitudes(nu_c, np.asarray(ref["L_nu"], float), bands, dist, edges)
        ref["colors"] = phot.colors(ref["mags"])
        for L in r["legs"].values():
            L["mags"] = phot.magnitudes(nu_c, np.asarray(L["L_nu"], float), bands, dist, edges)
            L["colors"] = phot.colors(L["mags"])
            L["dm"] = phot.delta_mag(L["mags"], ref["mags"])
            L["dcolor"] = {c: L["colors"][c] - ref["colors"][c] for c in L["colors"]}
    d["filter_set"] = filter_set
    return d


def band_fractions(row, lam_window, n_spec):
    """Each band's share of the reference bolometric luminosity."""
    edges = phot.nu_edges(*lam_window, n_spec)
    nu_c = np.sqrt(edges[1:] * edges[:-1])
    lnu = np.asarray(row["ref"]["L_nu"], float)
    tot = float(np.sum(lnu * np.diff(edges)))
    out = {}
    for b, (lam_lo, lam_hi) in phot.BANDS_PHOT.items():
        m = (nu_c >= C / (lam_hi * 1e-8)) & (nu_c <= C / (lam_lo * 1e-8))
        out[b] = float(np.trapezoid(lnu[m], nu_c[m])) / tot if tot > 0 else 0.0
    return out


def live_bands(row, lam_window, n_spec, frac_min=FRAC_MIN):
    f = band_fractions(row, lam_window, n_spec)
    return {b for b, v in f.items() if v >= frac_min}, f


def worst(path, frac_min=FRAC_MIN, filter_set="tophat"):
    d = json.loads(Path(path).read_text())
    if filter_set != "tophat":
        rephotometer(d, filter_set)
    lw, ns = d["lam_window"], d["n_spec"]
    out = {}
    for leg in LEGS:
        wm = wc = 0.0; wmt = wct = None
        for r in d["rows"]:
            if "skipped" in r:
                continue
            live, _ = live_bands(r, lw, ns, frac_min)
            L = r["legs"][leg]
            for b in live:
                v = L["dm"].get(b)
                if v is not None and np.isfinite(v) and abs(v) > wm:
                    wm, wmt = abs(v), (r["t_d"], b)
            for c in COLS:
                a, b = c.split("-")
                if a not in live or b not in live:
                    continue
                v = L["dcolor"].get(c)
                if v is not None and np.isfinite(v) and abs(v) > wc:
                    wc, wct = abs(v), (r["t_d"], c)
        out[leg] = {"dm": wm, "dm_at": wmt, "dcolor": wc, "dcolor_at": wct}
    return d, out


def report(path, frac_min=FRAC_MIN, filter_set="tophat"):
    d, w = worst(path, frac_min, filter_set)
    print(f"\n===== {Path(path).name}  [{filter_set} filters] =====")
    print(f"  {d['ion']}  core={d['core_law']}  rel={d.get('relativity')}  "
          f"v={d.get('v_core', float('nan')):.4g}-{d.get('v_out', float('nan')):.4g}c  "
          f"n={d['n']}x{len(d['seeds'])}  bands >= {100*frac_min:g}% of L_bol")
    print(f"  {'leg':10s}{'worst |dm|':>12}{'at':>16}{'worst |dcol|':>14}{'at':>16}")
    for leg in LEGS:
        v = w[leg]
        print(f"  {leg:10s}{v['dm']:12.3f}{str(v['dm_at']):>16}"
              f"{v['dcolor']:14.3f}{str(v['dcolor_at']):>16}")
    print(f"  noise floor (A_redist) = {w['A_redist']['dm']:.3f} mag")

    lw, ns = d["lam_window"], d["n_spec"]
    for leg in ("C_both", "C_binned"):
        print(f"\n  --- {leg} ---")
        print(f"  {'t/d':>6}{'S':>9}{'dF3800':>10}{'dm_bol':>9}"
              + "".join(f"{b:>8}" for b in "grizJHK")
              + "".join(f"{c:>9}" for c in COLS) + "   live")
        for r in d["rows"]:
            if "skipped" in r:
                continue
            live, _ = live_bands(r, lw, ns, frac_min)
            L = r["legs"][leg]
            f = L.get("dF_b3800")
            fs = "     --   " if f is None else f"{100*f:+9.1f}%"
            def q(v, b=None, w=8):
                if b is not None and b not in live:
                    return " " * (w - 2) + ". "
                return " " * (w - 4) + "--  " if v is None or not np.isfinite(v) \
                    else f"{v:+{w}.3f}"
            print(f"  {r['t_d']:6.2f}{r['band_S_band']:9.1f}{fs}{q(L['dm_bol'])}"
                  + "".join(q(L["dm"].get(b), b) for b in "grizJHK")
                  + "".join(q(L["dcolor"].get(c), None if all(x in live for x in c.split('-')) else "x", 9)
                            for c in COLS)
                  + "   " + "".join(sorted(live, key="grizJHK".index)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--frac", type=float, default=FRAC_MIN)
    ap.add_argument("--filters", choices=FILTER_SETS, default="tophat")
    ap.add_argument("files", nargs="*")
    a = ap.parse_args()
    files = a.files or sorted(str(p) for p in HERE.glob("observables_*.json"))
    for f in files:
        report(f if "/" in f else str(HERE / f), a.frac, a.filters)
