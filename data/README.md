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
