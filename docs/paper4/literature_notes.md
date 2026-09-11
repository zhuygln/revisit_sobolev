# Paper IV literature notes (2026-09-11)

Rule (inherited from `docs/paper3/literature_notes.md`): nothing is cited
from memory. `docs/paper4/references.bib` is `docs/paper3/references.bib`
(every entry verified there on 2026-09-03) plus the four entries below,
each verified on 2026-09-11 against its arXiv abstract page and the
publisher DOI before it went into the file. Sentences attributed to a
paper in the manuscript are paraphrases of its abstract or of a passage
read in the full text, never inferred.

## Entries added for Paper IV

| key | verified | used for |
|---|---|---|
| `lucy2005` — Lucy 2005, A&A 429, 19 (astro-ph/0409249), DOI 10.1051/0004-6361:20041656 | abstract + journal metadata | indivisible energy packets in time-dependent Monte Carlo transport (the light-curve method, Sect. 2.1 and 2.4) |
| `kasenbadnell2013` — Kasen, Badnell & Barnes 2013, ApJ 774, 25 (arXiv:1303.5788), DOI 10.1088/0004-637X/774/1/25 | abstract + journal metadata | the lanthanide line forest as the source of kilonova opacity (Introduction) |
| `ricigliano2024` — Ricigliano et al. 2024, MNRAS 529, 647 (arXiv:2311.15709), DOI 10.1093/mnras/stae572 | abstract + full text (v1): eq. 22–29 the thick-ejecta diffusion solution and the photosphere, eq. 47–50 and 59 the thin-ejecta luminosity and temperatures; Table 2 the secular-ejecta parameters; the public `xkn` code read for the conventions recorded in `sobolev/xkn.py` | the P1 state (Sect. 3) and the P1-xkn transport replacement test |
| `jplt2021` — Japan-Lithuania Opacity Database for Kilonova, version 2.1 (Kato et al.), http://dpc.nifs.ac.jp/DB/Opacity-Database/ | the database page and its version history; the Nd II/III data files retrieved 2026-09-10 (`sobolev/atomic_cache.py::build_cache_jplt`); the accompanying paper is `kato2024` (verified for Paper III) | the independent line list of Sect. 3 and 4.3 |

## Claims the manuscript attributes to the literature

| claim in the manuscript | source | how verified |
|---|---|---|
| The expansion-opacity formalism replaces the lines of a bin by an effective opacity ∝ Σ(1 − e^{−τ}) | `karp1977`, `eastman1993` | Paper III notes (full text of Eastman & Pinto 1993 §2) |
| Fontes et al. 2020 compared line-by-line, expansion and line-binned transport on a pure-Nd problem and found agreement to within several per cent at the peak, resolved Sobolev the brightest, expansion next, line-binned the faintest | `fontes2020` | full text: Appendix C (the problem's parameters are transcribed in `paper4/phase10_fontes/fontes.py`) and their light-curve comparison of the three treatments |
| The line-binned opacity is the one used in SuperNu-based kilonova models | `wollaeger2018`, `fontes2020` | Fontes et al. 2020 (full text: the line-binned opacity is the treatment used with SuperNu) |
| ARTIS treats fluorescence line by line with the Lucy macroatom; Shingles et al. argue it matters for the NIR | `shingles2023`, `collins2023`, `lucy2002`, `lucy2003` | Paper III notes (abstracts + full text of Shingles et al. 2023) |
| Collins et al. 2026 and Morag et al. 2026 identify line-opacity treatment as a leading uncontrolled choice | `collins2026`, `morag2026` | Paper III notes (abstracts) |
| P1's structure, the (1 − x²)³ profile, the photosphere at τ = 2/3, κ = 22.3 cm² g⁻¹ at Y_e = 0.2 | `ricigliano2024`, `tanaka2020` | full text (above); κ = 22.3 is the Tanaka et al. (2020) Y_e parametrisation as implemented in the public xkn code (`kappa_2_ye.py`); X_lan = 0.11 from Tanaka et al. 2020 Table 1 (recorded in the P1 state's metadata) |
| The Gillanders Ye-0.21a / Ye-0.29a abundance patterns and the 20× lanthanide reduction for the AT2017gfo-like state | `gillanders2022`, `collins2026` (Gillanders et al. 2026 in the state's metadata) | Paper III / Paper IV Phase 1 (`paper4/phase1_benchmarks/build.py` records the dataset DOI) |

Not cited from memory anywhere: numbers quoted from other papers appear
in the manuscript only as qualitative statements ("several per cent",
"the brightest"), never as digits, so that no literal from another paper
can be confused with a result of this one (`check_structure.py` enforces
the literal-number ban on the prose).
