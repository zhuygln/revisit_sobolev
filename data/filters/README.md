# Filter transmission curves

Retrieved 2026-09-02 from the SVO Filter Profile Service
(http://svo2.cab.inta-csic.es/theory/fps/, `getdata.php?format=ascii&id=<ID>`)
by `fetch_filters.py`. Two columns: wavelength (Å), transmission. DECam curves
are SVO's DetectorType = 1 (photon counter) system throughput; the 2MASS RSRs
are already per-photon. Both are integrated as photon-counting AB by
`sobolev.photometry.Passband`.

DECam griz + 2MASS JHKs is the AT2017gfo-like follow-up system; every curve
lies inside the 1000–30000 Å launch window used by Paper III's observables.
λ_eff below is the photon-weighted ∫λT dλ / ∫T dλ, which `Passband.lam_eff`
returns and `tests/test_passbands.py` checks to 1 %.

| file | SVO id | range (Å) | λ_eff (Å) | sha256[:16] |
|---|---|---|---|---|
| `decam_g.dat` | `CTIO/DECam.g` | 3815–5585 | 4827 | `35936fc027e8e885` |
| `decam_r.dat` | `CTIO/DECam.r` | 5400–7350 | 6432 | `15c0d56823563a98` |
| `decam_i.dat` | `CTIO/DECam.i` | 6750–8700 | 7827 | `f975476f14a3ef22` |
| `decam_z.dat` | `CTIO/DECam.z` | 8250–10150 | 9179 | `c5e6ab3c3806b0a0` |
| `2mass_J.dat` | `2MASS/2MASS.J` | 10620–14500 | 12411 | `9d0c86e1279f17ff` |
| `2mass_H.dat` | `2MASS/2MASS.H` | 12890–19140 | 16514 | `f6c6074ba5dae493` |
| `2mass_Ks.dat` | `2MASS/2MASS.Ks` | 19000–23990 | 21656 | `58969bfd998b9080` |
