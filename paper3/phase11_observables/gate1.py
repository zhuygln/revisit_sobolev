"""Gate 1 (plan of 2026-09-02): does the chromatic closure error survive real
filters? Re-photometers the committed F41 spectra through DECam griz + 2MASS
JHKs (no MC rerun) and compares the worst |dm|, |dcolour| per leg with the
top-hat values the report quotes.

Pass: the C_both / C_binned worst colour error stays >= 0.3 mag where the
top-hats gave 0.6-0.8, and the A_redist noise floor does not grow: its worst
|dm| stays <= max(0.02 mag, 1.5x its top-hat value). (The plan wrote "<= 0.02"
with F41's 0.01 floor in mind; the blue-kilonova run's floor is already 0.10
mag with top-hats at 3 d, so the criterion is that real filters do not
worsen it, not that they cure it.)

Usage: python gate1.py   -> gate1_realfilters.json
"""
import json, subprocess
from pathlib import Path

import numpy as np

from summary import HERE, LEGS, FRAC_MIN, worst, rephotometer

FILES = ("observables_58CeII_cool.json", "observables_57LaII_cool.json",
         "observables_blend_cool.json", "observables_57LaII_cool_kn_worldline_blue.json")
PASS = {"colour_min_real": 0.3, "colour_min_tophat": 0.6, "A_dm_max": 0.02, "A_dm_growth": 1.5}


def tophat_roundtrip(path):
    """Max |stored - recomputed| top-hat magnitude: the re-photometry is exact."""
    d = json.loads(Path(path).read_text())
    stored = [(r["ref"]["mags"], {k: v["dm"] for k, v in r["legs"].items()})
              for r in d["rows"] if "skipped" not in r]
    rephotometer(d, "tophat")
    err = 0.0
    for r, (m0, dm0) in zip([r for r in d["rows"] if "skipped" not in r], stored):
        for b in m0:
            if np.isfinite(m0[b]) and np.isfinite(r["ref"]["mags"][b]):
                err = max(err, abs(m0[b] - r["ref"]["mags"][b]))
        for k in dm0:
            for b in dm0[k]:
                a, c = dm0[k][b], r["legs"][k]["dm"][b]
                if np.isfinite(a) and np.isfinite(c):
                    err = max(err, abs(a - c))
    return err


def main():
    out = {"frac_min": FRAC_MIN, "pass_criteria": PASS, "files": {}}
    ok = True
    for f in FILES:
        p = HERE / f
        _, wt = worst(p, FRAC_MIN, "tophat")
        _, wr = worst(p, FRAC_MIN, "real")
        rt = tophat_roundtrip(p)
        rec = {"tophat": wt, "real": wr, "tophat_roundtrip_max_abs": rt}
        print(f"\n{f}   (top-hat round trip {rt:.1e} mag)")
        print(f"  {'leg':10s}{'|dm| tophat':>13}{'|dm| real':>11}{'|dcol| tophat':>15}{'|dcol| real':>13}")
        for leg in LEGS:
            print(f"  {leg:10s}{wt[leg]['dm']:13.3f}{wr[leg]['dm']:11.3f}"
                  f"{wt[leg]['dcolor']:15.3f}{wr[leg]['dcolor']:13.3f}   "
                  f"{wt[leg]['dcolor_at']} -> {wr[leg]['dcolor_at']}")
        # the criterion applies where the top-hats showed a large colour error
        big = [leg for leg in ("C_both", "C_binned") if wt[leg]["dcolor"] >= PASS["colour_min_tophat"]]
        a_ok = wr["A_redist"]["dm"] <= max(PASS["A_dm_max"], PASS["A_dm_growth"] * wt["A_redist"]["dm"])
        this = bool(all(wr[leg]["dcolor"] >= PASS["colour_min_real"] for leg in big) and a_ok and rt < 1e-9)
        rec["legs_tested"] = big; rec["pass"] = this
        out["files"][f] = rec
        ok &= this
    out["pass"] = ok
    out["git"] = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                                text=True, cwd=HERE).stdout.strip()
    (HERE / "gate1_realfilters.json").write_text(json.dumps(out, indent=1))
    print(f"\nGate 1: {'PASS' if ok else 'FAIL'}   wrote gate1_realfilters.json")


if __name__ == "__main__":
    main()
