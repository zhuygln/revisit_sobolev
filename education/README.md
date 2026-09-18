# From Photon Random Walks to Effective Atomic Fluorescence

*A Computational Introduction to Kilonova Radiative Transfer* — the education
subproject of `revisit_sobolev`: a Quarto book, executable from first
principles, that builds toward the planned effective-redistribution ("Paper
B") calculation. The PI's specification and every decision are verbatim in
[`plan_review.md`](plan_review.md).

**The rule: educational demos ≠ publication evidence.** The demos teach and
validate concepts; the production code (`sobolev/`,
`paper2/phase1/forest_mc.py`, `paper3/redistribution/kernel.py`) uses the
same ideas with real atomic data and controlled physics. The toy package
`src/rtedu/` never imports production code; only the bridge chapter (E4)
compares the two.

## Layout

| path | what |
|---|---|
| `chapters/NN_*.qmd` | the polished teaching text, one chapter per file, every chapter in the same eleven-section template |
| `notebooks/NN_*.ipynb` | Level-1 code: 20–50 transparent lines per idea, validated against `rtedu`; executed headlessly into `notebooks/executed/` (ignored) |
| `src/rtedu/` | Level-2 code: the reusable toy package (not installed; `PYTHONPATH=$(CURDIR)/src`) |
| `tests/` | the physics as tests (`make education-test`; kept out of the research suite) |
| `scripts/` | the video generators, `check_no_literals.py`, `check_determinism.py` |
| `figures/`, `videos/` | generated, committed (PNG + PDF; GIF + an HTML/JS fragment) |
| `data/generated/chNN.json` → `data/results.json` | every number the chapters quote; transactional |

## Build

    cd education
    make education          # tests -> notebooks -> videos -> render (_book/index.html)
    make education-fast     # the same without the videos
    make education-test     # the physics tests only
    make education-check    # tests + the no-literals and determinism checks

Quarto comes from `pip install quarto-cli` (the `education` extra of
`pyproject.toml`); videos are GIF (pillow) plus an HTML5/JS fragment
(matplotlib) because ffmpeg is absent; MP4 is written when an ffmpeg writer
exists.

## Reproducibility

Every displayed number comes from code. Notebooks seed their generator from
`rtedu.SEEDS`, write `data/generated/chNN.json`, and `rtedu.results.merge`
builds `data/results.json`; `make education-notebooks` deletes both first.
Chapters quote values only through inline `{python}` expressions;
`scripts/check_no_literals.py` enforces that inside `::: {.result}` regions
(pedagogical constants elsewhere are free). `results.json`, figures and GIFs
are byte-stable across builds; HTML animations are compared by their
embedded frames (`scripts/check_determinism.py`).

## Milestones

| stage | chapters | exit criterion (PI) | state |
|---|---|---|---|
| E1 | 0–4 | you can derive and implement how a packet finds its next line interaction: ν_com(r), the next resonance, its τ_S, the coin P_int = 1 − e^{−τ_S} (the chapter-4 exercise) | this PR |
| E2 | 5–8 | you can explain precisely what ε approximates and why a macroatom is more physical | |
| E3 | 9–12 | you can derive R_ij, identify the information it discards, and explain when a neural surrogate would be justified | |
| E4 | 13 + bridge | one complete toy RT simulation produces an SED and light curve with ε, R and the macroatom | |

E1 stops exactly at: *given a packet and a line forest, how the solver finds
the next Sobolev resonance and decides whether the packet interacts.* No
`ToyAtom`, ε, macroatom or R_ij code is in E1.
