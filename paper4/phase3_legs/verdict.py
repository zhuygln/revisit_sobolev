"""Paper IV Gate 2, pre-declared (paper4/README.md): the verdict from the
leg files of paper4/phase2_energy/legs.py.

Primary B2 - R2; controls A2 - R2 and C2 - R2; live bands only. A band at
an epoch is live when its magnitude is finite, the band carries >= 1 % of
the reference L_bol and it is brighter than the depth limit at 40 Mpc
(griz 23.5, JHK 21.5) -- the Paper III mask (sensitivity.py), applied to R2.

    Green   |dm_B2| >= 0.5 mag in >= 2 live bands at >= 2 epochs;
            |dm_A2| <= 0.1 mag everywhere; g brighter (dm_g < 0) and
            K fainter (dm_K > 0) kept. Also reports whether C2 is within
            0.1 mag of B2 (redistribution compresses, opacity is the culprit).
    Yellow  0.2 <= max|dm_B2| < 0.5 mag, controls in bounds.
    Red     max|dm_B2| < 0.2 mag AND ||R1 - R2|| > 1/2 ||R1 - B1||.
    Gray    anything else: a control out of bounds, B2 - R2 small with
            R1 - R2 also small, a sign pattern against the causal reading.

The ladder decomposition R1E - R1 (bookkeeping), R2 - R1E (transition
probabilities), B2 - R2 (opacity coarse-graining) is printed per band.
Convergence (Phase 6) and the adequacy trigger are checked by their own
drivers; a failure there turns any outcome Gray.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (ROOT, ROOT / "paper3/phase12_grid"):
    sys.path.insert(0, str(_p))

from sensitivity import MAG_LIMIT, FRAC_MIN, BANDS   # noqa: E402
from sobolev.photometry import BANDS_PHOT, nu_edges  # noqa: E402
from sobolev.constants import C                      # noqa: E402

THRESH = dict(green_mag=0.5, green_bands=2, green_epochs=2, control_max=0.1, yellow_min=0.2,
              red_max=0.2, compress_within=0.1)


def band_fractions(row, leg="R2"):
    """Fraction of the reference's window luminosity in each band (top-hat
    edges of sensitivity's mask, on the stored L_nu)."""
    lnu = np.asarray(row["legs"][leg]["L_nu"], float)
    edges = nu_edges(*row["lam_window"], row["n_spec"])
    nu_c = np.sqrt(edges[1:] * edges[:-1]); dnu = np.diff(edges)
    tot = float(np.sum(lnu * dnu))
    out = {}
    for b, (lo_a, hi_a) in BANDS_PHOT.items():
        m = (nu_c >= C / (hi_a * 1e-8)) & (nu_c < C / (lo_a * 1e-8))
        out[b] = float(np.sum(lnu[m] * dnu[m]) / tot) if tot > 0 else 0.0
    return out


def live_bands(row, leg="R2"):
    frac = band_fractions(row, leg)
    mags = row["legs"][leg]["mags"]
    return [b for b in BANDS if np.isfinite(mags.get(b, np.nan)) and frac[b] >= FRAC_MIN
            and mags[b] <= MAG_LIMIT[b]]


def norm(dm, bands):
    v = np.array([dm[b] for b in bands if np.isfinite(dm.get(b, np.nan))])
    return float(np.sqrt(np.sum(v ** 2))) if v.size else np.nan


def analyse(rows):
    per_epoch = []
    for row in rows:
        L = row["legs"]
        live = live_bands(row)
        d = dict(name=row.get("name"), t_d=row.get("t_d"), shell=row.get("shell"), live=live)
        for tag in ("B2", "A2", "C2", "Bbin2", "Cbin2", "R1", "R1E", "B1", "B1E"):
            if tag in L:
                d[f"dm_{tag}"] = {b: L[tag]["dm_vs_R2"][b] for b in live}
        if "B2" in L:
            d["max_B2"] = max((abs(v) for v in d["dm_B2"].values()), default=np.nan)
            d["n_B2_big"] = sum(abs(v) >= THRESH["green_mag"] for v in d["dm_B2"].values())
            d["sign_ok"] = (d["dm_B2"].get("g", 0.0) < 0) and (d["dm_B2"].get("K", 0.0) > 0) if {"g", "K"} <= set(live) else None
        if "A2" in L:
            d["max_A2"] = max((abs(v) for v in d["dm_A2"].values()), default=np.nan)
        if "C2" in L and "B2" in L:
            d["C2_minus_B2_max"] = max((abs(L["C2"]["dm_vs_R2"][b] - L["B2"]["dm_vs_R2"][b]) for b in live), default=np.nan)
        if all(t in L for t in ("R1", "B1", "R2")):
            d["norm_R1_R2"] = norm({b: L["R2"]["dm_vs_R1"][b] for b in live}, live)
            d["norm_R1_B1"] = norm({b: L["B1"]["dm_vs_R1"][b] for b in live}, live)
        # the ladder decomposition per band
        if all(t in L for t in ("R1", "R1E", "R2", "B2")):
            d["ladder"] = {b: dict(bookkeeping=L["R1E"]["dm_vs_R1"][b], probabilities=L["R2"]["dm_vs_R1"][b] - L["R1E"]["dm_vs_R1"][b],
                                   opacity=L["B2"]["dm_vs_R2"][b]) for b in live}
        per_epoch.append(d)
    return per_epoch


def verdict(per_epoch):
    have = [d for d in per_epoch if "max_B2" in d]
    if not have:
        return "gray", "no B2 - R2 available"
    max_b2 = max(d["max_B2"] for d in have)
    control_ok = all(d.get("max_A2", 0.0) <= THRESH["control_max"] for d in have)
    signs_ok = all(d["sign_ok"] for d in have if d.get("sign_ok") is not None)
    n_epochs_big = sum(d["n_B2_big"] >= THRESH["green_bands"] for d in have)
    red_ratio = [d["norm_R1_R2"] > 0.5 * d["norm_R1_B1"] for d in have if "norm_R1_R2" in d]
    if not control_ok:
        return "gray", "control A2 - R2 exceeds 0.1 mag in a live band"
    if n_epochs_big >= THRESH["green_epochs"] and signs_ok:
        return "green", f"|dm_B2| >= 0.5 mag in >= 2 live bands at {n_epochs_big} epochs, signs kept"
    if n_epochs_big >= 1 and len(have) < THRESH["green_epochs"] and signs_ok:
        return "green-provisional", f"criterion met at the {len(have)} epoch(s) available; needs {THRESH['green_epochs']}"
    if max_b2 >= THRESH["yellow_min"] and not signs_ok:
        return "gray", (f"max |dm_B2| = {max_b2:.2f} mag but the sign pattern (g brighter, K fainter) "
                        "is not kept: contradicts the causal reading, diagnose before interpreting")
    if THRESH["yellow_min"] <= max_b2 < THRESH["green_mag"]:
        return "yellow", f"max |dm_B2| = {max_b2:.2f} mag, signs kept"
    if max_b2 < THRESH["red_max"] and red_ratio and all(red_ratio):
        return "red", f"max |dm_B2| = {max_b2:.2f} mag and ||R1 - R2|| > 1/2 ||R1 - B1||"
    return "gray", f"max |dm_B2| = {max_b2:.2f} mag, signs_ok = {signs_ok}, red ratio = {red_ratio}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rows = [json.loads(Path(f).read_text()) for f in a.files]
    per = analyse(rows)
    for d in per:
        print(f"{d['name']} t={d['t_d']:g} d shell {d['shell']}: live {d['live']}")
        for tag in ("R1E", "R1", "B1", "B2", "Bbin2", "A2", "C2", "Cbin2"):
            if f"dm_{tag}" in d:
                print(f"   dm({tag} - R2): " + "  ".join(f"{b}={v:+.2f}" for b, v in d[f'dm_{tag}'].items()))
        if "ladder" in d:
            print("   ladder (bookkeeping / probabilities / opacity):")
            for b, l in d["ladder"].items():
                print(f"     {b}: {l['bookkeeping']:+.2f} / {l['probabilities']:+.2f} / {l['opacity']:+.2f}")
        for k in ("max_B2", "max_A2", "C2_minus_B2_max", "norm_R1_R2", "norm_R1_B1", "sign_ok"):
            if k in d:
                print(f"   {k} = {d[k]}")
    v, why = verdict(per)
    print(f"\nGATE 2: {v.upper()} -- {why}")
    if a.out:
        Path(a.out).write_text(json.dumps(dict(thresholds=THRESH, per_epoch=per, verdict=v, why=why),
                                          indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
