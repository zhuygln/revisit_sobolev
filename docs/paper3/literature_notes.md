# Paper III literature search (2026-09-03)

Rule: nothing is cited from memory. Every entry below was verified by its
abstract (arXiv or publisher page, fetched on the date given) before it went
into `references.bib`, and the journal metadata by Crossref on the same day.
Sentences quoted from a paper's body were read from the full text (PDF
extracted locally), not inferred.

## Query 1 — the σ_sys model-error allowance in kilonova inference

Search: "kilonova light curve fitting systematic error 1 mag model
uncertainty added in quadrature Bayesian inference" (2026-09-03).

| candidate | verified | verdict |
|---|---|---|
| Coughlin et al. 2018, MNRAS 480, 3871 (arXiv:1805.09371) | abstract + full text: "We assign model uncertainties of 1 mag added in quadrature with the statistical error in the measured photometry … designed to capture difficult-to-quantify systematic uncertainties, such as those in the electron fraction and heating rate" | **cite** — the origin of the 1 mag convention |
| Pang et al. 2023, Nat. Commun. 14, 8352 (arXiv:2205.08513) | abstract + full text: likelihood eq. (16) with σ_sys added in quadrature to every point; σ_sys = 0.5, 1, 2 mag compared, "0.5 mag is underestimating the systematic uncertainty … the 2 mag is overcompensating … we concluded that the σ_sys of 1 mag is a sensible choice" | **cite** — NMMA's allowance and its rationale |
| Jhawar et al. 2025, PRD 111, 043046 (arXiv:2410.21978) | abstract: time- and filter-dependent systematic uncertainties; "a systematic error below 1 mag between 1 to 5 days after the merger" for AT2017gfo | **cite** — the allowance becoming structured; our result says what structure a transport closure would give it |
| Hussenot-Desenonges et al. 2025 (arXiv:2505.21392) | abstract: "the systematic error margin hyperparameter σ_sys, which can be exploited as a metric for a model's goodness-of-fit" | **cite** — σ_sys read as a misfit statistic, which is what our χ²-equivalent compares to |
| Dietrich et al. 2020, Science 370, 1450 (arXiv:2002.11355) | abstract verified; the 1 mag statement is in the supplement, not retrieved | not cited (Pang 2023 covers NMMA) |
| Brethauer et al. 2024, ApJ 975, 213 | abstract: atomic-data and thermalization choices move X_lan estimates by an order of magnitude and masses by 25–40 %; "observable properties including color and decay rate prove highly model-dependent" | **cite** — the systematic budget our closure error sits inside |

## Query 2 — line-interaction closure errors propagated into multiband observables

Search: "kilonova radiative transfer expansion opacity vs line-by-line
comparison bias colours" (2026-09-03), plus the codes named in the plan.

| candidate | verified | verdict |
|---|---|---|
| Kawaguchi, Shibata & Tanaka 2018, ApJL 865, L21 (arXiv:1806.04088) | abstract: 2-D transport with multiple ejecta components reproduces AT2017gfo's light curves and photospheric velocity | **cite** — a binned-opacity code used for inference on AT2017gfo |
| Wollaeger et al. 2018, MNRAS 478, 3298 (arXiv:1705.07084) | abstract: SuperNu, LTE opacities for representative elements, morphology/composition study | **cite** — same role |
| Bulla 2019, MNRAS 489, 5037 (arXiv:1906.04205) | abstract: POSSIS, wavelength- and time-dependent opacities, parameter-space studies | **cite** — the model family behind NMMA's Bu2019lm |
| Vieira et al. 2023, ApJ 944, 123 (arXiv:2209.06951) | abstract: SPARK, Bayesian abundance retrieval on AT2017gfo spectra through TARDIS | **cite** — spectral inference through a Sobolev/macroatom code |
| Shingles et al. 2023; Collins et al. 2023, 2026 | already in Paper I's bib | cited — the line-by-line reference class |
| Tanaka et al. 2020; Kato et al. 2024; Fontes et al. 2020, 2026; Flörs et al. (gsi_atomic); Domoto et al. 2022; Banerjee et al. 2022, 2024; Gillanders et al. 2022; Morag et al. 2026 | already in Paper I's bib | cited as before |
| no 2025–26 paper found that measures the expansion-opacity closure's error on multiband photometry against a same-code line-by-line reference | — | the gap the paper states; the search was one query plus the plan's list, so the statement is worded as "we are not aware of", not "there is none" |

## Query 3 — the observations and the passbands

| candidate | verified | verdict |
|---|---|---|
| Tanvir et al. 2017, ApJL 848, L27 | abstract: NIR counterpart, "much slower evolution in the near-infrared Ks-band compared to the optical" | **cite** — the NIR is where the signature lives |
| Chornock et al. 2017, ApJL 848, L19 | abstract: NIR spectra 1.5–10.5 d, peaks near 1.07 and 1.55 µm, model with 0.04 M_⊙, 0.1c, X_lan 10⁻² | **cite** — the central grid point is that model's parameters |
| Watson et al. 2019, Nature 574, 497 | abstract: Sr identification | **cite** — spectral identifications rely on line treatment |
| Flaugher et al. 2015, AJ 150, 150 | abstract | **cite** — DECam |
| Skrutskie et al. 2006, AJ 131, 1163 | abstract | **cite** — 2MASS |
| Rodrigo, Solano & Bayo 2012 (IVOA report); Rodrigo et al. 2024, A&A 689, A93 | the SVO FPS citation page names these; abstracts verified | **cite** — the curves in `data/filters/` |

## Budget

16 new entries against the plan's "≤ 12": the four extra are the passband
references (three) and Hussenot-Desenonges 2025, which is the one paper that
reads σ_sys as a goodness-of-fit statistic. Every entry has a DOI or arXiv id.
Main-text citations are counted by `check_structure.py` against the 50 cap.

## Source-model references (added 2026-09-03, Methods only)

| key | verified | how |
|---|---|---|
| `korobkin2012` MNRAS 426, 1940 | title/authors/volume/pages via Crossref API | the heating-rate fitting form ε(t) ∝ [1/2 − arctan((t−t0)/σ)/π]^1.3 |
| `barnes2016` ApJ 829, 110 | Crossref | thermalization efficiency form, bilinear interpolation in (log M, v) |
| `arnett1982` ApJ 253, 785 | Crossref | one-zone diffusion luminosity |

Budget: 19 new entries against the plan's ≤12; the three above are Methods-only citations of standard ingredients and are not counted toward the ≤50 main-text limit.
