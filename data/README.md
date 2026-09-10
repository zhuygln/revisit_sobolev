# data/

Raw atomic data lives here and is **not** tracked by git — the files are large and
externally versioned. Only this README is committed.

## Phase 0C

One ion from the **GSI v2 calibrated lanthanide database** (La–Yb II/III) is enough.
Start with **La II**.

Record here what you actually downloaded, so a result can be traced back to its input:

| File | Ion | Source / version | Date |
|---|---|---|---|
| `57LaII_transitions_calib.txt` | La II | GSI Database for Kilonova Radiative Transfer, Zenodo record [19335084](https://zenodo.org/records/19335084) (latest under concept DOI `10.5281/zenodo.15835360`; paper DOI `10.1103/jxqw-7ynk`), published 2026-03, CC-BY 4.0. Extracted from `GSI_lanthanides_calibrated_transitions.zip`; 17,743 E1 transitions, methods `xmatch`/`shifted`. | 2026-08-14 |
| `57LaII_levels_calib.txt` | La II | Same record, `GSI_lanthanides_calibrated_levels.zip`; 472 levels. | 2026-08-14 |
| `58CeII_transitions_calib.txt` / `58CeII_levels_calib.txt` | Ce II | Same record; 2,829 levels. Note: half-integer J written as fractions (`7/2`) — handled by `sobolev.populations.parse_j`. | 2026-08-15 |
| `58CeIII_transitions_calib.txt` / `58CeIII_levels_calib.txt` | Ce III | Same record; extracted for future use, not yet in experiments. | 2026-08-15 |
| `60NdII_transitions_calib.txt` / `60NdII_levels_calib.txt` | Nd II | Same record; 9,994 levels, 3,336,077 transitions (687 MB — the largest ion in the archive). Half-integer J as fractions, as Ce II. Extracted to unblock Paper III R5/P10; not yet in experiments. | 2026-08-29 |

A 20-row excerpt of the transitions file is committed at
`tests/data/57LaII_transitions_calib_excerpt.txt` for format regression tests.

## Re-download

The whole record was re-fetched on **2026-08-29** onto a second machine and
the per-ion files verified against the counts above (La II 472 levels /
17,743 transitions, Ce II 2,829 levels — all exact). Both source archives
(`GSI_lanthanides_calibrated_{levels,transitions}.zip`, 875 MB together) are
kept here so the remaining 23 lanthanide ions can be extracted without another
download; they are gitignored along with everything else in `data/`.

The archives are Mac-made: skip the `__MACOSX/` entries and extract the
per-ion `.txt` files **flat** into `data/`, which is where `load_gsi` looks.
Commands are in [../docs/sedona/SETUP.md](../docs/sedona/SETUP.md) §2.

## Paper IV (2026-09-09)

Two committed exceptions to the "nothing in data/ is tracked" rule:
`data/abundances/` (small pattern tables, provenance in each file's header)
and this ledger. `data/cache/` holds the compact per-ion `.npz` caches that
`sobolev/atomic_cache.py` builds from the archives (derived, gitignored).

| Item | Content | Source | Date |
|---|---|---|---|
| `abundances/solar_r.csv` | solar-system r-process residual pattern by mass, La–Lu: `w = X_solar × f_r` | Prantzos, Abia, Cristallo, Limongi & Chieffi 2020, MNRAS 491, 1832, **Table 4** (Lodders et al. 2009 mass fractions; the study's s/r/p fractions); transcribed from arXiv:1911.02545 by text extraction and checked row by row (`w = X × f_r` to 0.2 %) | 2026-09-09 |
| P1 anchors | secular component of the xkn radiative-transfer comparison: `M_sec = 2.64e-2 Msun`, `v_rms = 0.06c`, `Y_e = 0.20`, `s = 10 k_B/baryon`, `τ_exp = c/v_rms ≈ 17 ms`; density `ρ(t,x) = ρ0 (t0/t)^3 (1 − x^2)^3`, `x = v/v_max` (**eq. 25**); heating from Wanajo et al. 2014 tracers, thermalisation Barnes et al. 2016 | Ricigliano et al. 2024, MNRAS 529, 647, §5.1–5.2 (arXiv:2311.15709) | 2026-09-09 |
| `gillanders2022/additional_paper_resources.zip` (gitignored) and `abundances/gillanders2022_Ye-0.29a.csv`, `abundances/gillanders2022_Ye-0.21a.csv` (committed) | the complete composition profiles of Gillanders et al. 2022 Table 2 — every element with its mass fraction for the twelve Y_e bins (`composition_profiles_complete.ascii`, sha256 `9f06304b…819bd1`; zip sha256 `793b5ae2…8784bc`, 5 468 757 bytes); the two profiles used are extracted verbatim (sums 0.9998 / 0.9984 before renormalisation; X_LN(57–70) = 4.986×10⁻² / 2.989×10⁻¹, Table 2's 4.99×10⁻² / 2.99×10⁻¹) | Gillanders, Smartt, Sim, Bauswein & Goriely 2022, MNRAS 515, 631; QUB dataset DOI 10.17034/404fbfbe-5f47-42ff-a7d0-12e7c447ebff (downloaded by the PI in a browser: the server challenges scripted downloads) | 2026-09-10 |
| P1 composition | `Ye−0.21a` with the lanthanides (Z = 57–70) rescaled to `X_lan = 0.11` and the other elements to 0.89; the 0.11 is Tanaka, Kato, Gaigalas & Kawaguchi 2020, MNRAS 496, 1369, **Table 1** at Y_e = 0.20 (X(La) = 1.1×10⁻¹ from Wanajo et al. 2014 yields; per-element yields only in their Fig. 11). Robustness: `P1r1` = Ye−0.21a at its own X_LN = 0.30, `P1r2` = `solar_r` at 0.11 (the previous provisional pattern) | PI decision 2026-09-10 | 2026-09-10 |
| P2 composition | `Ye−0.29a` with Z = 57–70 rescaled to `X_LN = 2.5×10⁻³` exactly (the 2026 paper's quoted value; a 20× reduction to its rounding) and the other elements to 1 − X_LN | Gillanders et al. 2026, §4 | 2026-09-10 |
| standard atomic weights | `sobolev/abundances.py::ATOMIC_MASS`, IUPAC conventional values | — | 2026-09-09 |
