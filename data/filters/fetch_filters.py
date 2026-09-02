"""Regenerate the transmission curves in this directory from the SVO Filter
Profile Service (http://svo2.cab.inta-csic.es/theory/fps/).

    python fetch_filters.py            # writes decam_*.dat, 2mass_*.dat

The files are committed; this script exists so their origin is reproducible.
`tests/test_passbands.py` pins their SHA-256, so a re-fetch that changes a
curve (SVO does revise them) fails loudly rather than silently moving the
photometry.
"""
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
URL = "http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id={}"
FILTERS = {"decam_g": "CTIO/DECam.g", "decam_r": "CTIO/DECam.r",
           "decam_i": "CTIO/DECam.i", "decam_z": "CTIO/DECam.z",
           "2mass_J": "2MASS/2MASS.J", "2mass_H": "2MASS/2MASS.H",
           "2mass_Ks": "2MASS/2MASS.Ks"}

if __name__ == "__main__":
    for name, fid in FILTERS.items():
        with urllib.request.urlopen(URL.format(fid), timeout=60) as r:
            (HERE / f"{name}.dat").write_bytes(r.read())
        print(f"{name}.dat <- {fid}")
