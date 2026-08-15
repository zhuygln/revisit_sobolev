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

A 20-row excerpt of the transitions file is committed at
`tests/data/57LaII_transitions_calib_excerpt.txt` for format regression tests.
