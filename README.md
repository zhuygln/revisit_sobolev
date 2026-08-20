# revisit_sobolev

Kilonova radiative-transfer codes cannot resolve the millions of narrow
lanthanide absorption lines that shape a kilonova's light, so they all use a
shortcut. That shortcut is usually called "the Sobolev approximation." This
project shows it is really **two** approximations with error budgets differing
by roughly an order of magnitude, and measures both against resolved
calculations on calibrated atomic data.

**The result.** Against frequency-resolved transport as truth, on realistic
La II line forests:

| treatment | band-flux error |
|---|---|
| Sobolev approximation proper (`e^−τs` per line) | **≲0.5%** at τ_max = 5, ≲2% overall |
| expansion opacity (its usual implementation) | **+38–48%** |

The gap comes from a per-resonance substitution `τ → 1−e^−τ` that belongs to
expansion opacity alone — not to Sobolev's locality or isolation assumptions —
and it does not vanish as line widths shrink. The Sobolev approximation's own
residual, by contrast, switches off entirely at physical line widths: it is a
finite-region boundary effect, not a failure of either Sobolev assumption. Across 36 conditions spanning
four wavelength windows, three epochs and three ion mixtures, the separation
holds, with the realized maximum line optical depth as the controlling
variable at fixed line width and geometry.

## The paper

**[docs/paper/manuscript.pdf](docs/paper/manuscript.pdf)** — 19 pp in MNRAS
format, readable directly on GitHub. Written for a reader with no
radiative-transfer background: a primer builds every concept from scratch and
an appendix derives every formula.

Source is [manuscript.tex](docs/paper/manuscript.tex); rebuild with
`cd docs/paper && make` (needs pdflatex + bibtex).

## Findings

| # | Finding | Detail |
|---|---|---|
| F1 | Gradients symmetric about a resonance cancel the Sobolev error at leading order, so gradient-based validity diagnostics are conservative | [§4.2](docs/results_report.md) |
| F2 | Resolved transport cost 378× expansion opacity on a standard example | §4.4 |
| F3 | The O(v_bulk/c) offset is an artifact of frozen-snapshot integration, not relativity — superseded by F11 | §4.5, §4.14 |
| F4 | Expansion opacity attenuates by `e^−(1−e^−τ)` per crossing, not `e^−τ` | §4.5 |
| F5 | That error is per-resonance: it does not average away with line count, only as τ→0 | §4.6 |
| F6 | Δ_expansion = a strength-set floor plus a width-growing term; the floor survives to thermal widths | §4.8 |
| F7 | Δ is non-monotonic in forest density — maximal for strong *sparse* forests, suppressed at full blanketing | §4.10 |
| F8 | Solver-vs-Monte-Carlo offsets traced to the thermal-emission convention, not physics; matched, the codes agree to 0.1% on the La II forest — inside the measured 0.3% Monte Carlo noise | §4.11 |
| F9 | The strength floor is expansion opacity’s alone; Sobolev’s own residual is small and width-driven (see F12) | §4.12 |
| F10 | The separation is universal across windows, epochs and ion mixes; realized τ_max controls it at fixed v_D and geometry | §4.13 |
| F11 | Neglect of light-travel-time evolution is a distinct approximation: frozen τ/τ_S = (1−β)/γ vs worldline 1/γ, which has no O(β) term. SEDONA in steady-iterate mode solves the frozen problem, so transport treatment is a third convention cross-code comparisons must match | §4.14 |
| F12 | Overlap is inert in pure absorption (optical depths add exactly); the Sobolev residual is a finite-region boundary effect ∝ v_D/Δv_shell, negligible at thermal widths | §4.15 |
| F13 | The expansion-opacity error is bin-width invariant from 0.025 to 41 lines per bin — intrinsic to the formalism, not a usage artifact | §4.16 |
| F14 | Neither candidate reference code can answer Paper II alone: SEDONA has no line branching at all, TARDIS has no expansion opacity. A ~250-line branching Sobolev MC, validated to ≤2.1σ against the analytic Sobolev leg, supplies the same-code differential both lack | §4.17 |

Full write-up with figures and numbers:
**[docs/results_report.md](docs/results_report.md)**.

## Layout

```
sobolev/       the package: constants, line profiles, optical depths,
               Boltzmann populations, GSI parsing, the deterministic
               formal solver, and the analytic Sobolev leg
experiments/   one directory per experiment; generators and comparison
               scripts committed, SEDONA run outputs gitignored
notebooks/     Phase 0: single-line toy model, GSI line spacing
paper2/        Paper II Phase 0: the SEDONA fluorescence source audit, the
               TARDIS install record, and the branching Monte Carlo
tests/         68 tests pinning the physics of every module
docs/          results report, lab notebook, planning inputs, paper/
data/          raw atomic data (gitignored; provenance in data/README.md)
outputs/       working figures (gitignored; committed copies in docs/figures/)
```

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]" h5py
pytest                      # 68 passed
```

Atomic data is not committed — see [data/README.md](data/README.md) for the
Zenodo record and re-download instructions. SEDONA lives outside this repo;
its no-root WSL2 build recipe is in the lab notebook, §5.

Per-experiment reproduction commands are in
[docs/results_report.md](docs/results_report.md) §7.

## How this got built

The project ran as a staged feasibility study — a weekend toy model first,
then cross-code validation, then real atomic data, then parameter maps — with
a go/no-go decision at each stage. That history, including the dead ends and
the bugs that mattered, is in
**[docs/lab_notebook.md](docs/lab_notebook.md)**; the organized results are in
the report §4. The original plans are preserved unmodified as
[docs/babystep_plan.md](docs/babystep_plan.md) and
[docs/research_requirements.md](docs/research_requirements.md).

## Status and open items

Phase 0 through the validity maps is complete; the manuscript covers
everything above. Outstanding:

1. ~~The residual Sobolev error.~~ **Resolved (F12).** Most of it was a
   normalization artifact; the remainder is a finite-region boundary effect
   that vanishes at physical line widths. Overlap was excluded both
   analytically (optical depths add) and numerically.
2. **Realistic velocities.** All results use 1000–3000 km/s shells. The
   machinery now reaches 0.1–0.3c: both the solver and the analytic leg carry
   a co-evolving medium and worldline transport. Two corrections came out of
   building it — the dominant term is light-travel *dilution* ((1−β)² = 51% at
   β = 0.3), not the 4.6% relativistic one, and the CD/CP resonance surfaces
   stay irrelevant under worldline transport, where the locus is linear.
   A pilot finds Δ_Sobolev β-independent but Δ_expansion **growing** with β;
   it is synthetic-forest only and not yet a finding (lab notebook §9n).
3. **Scattering and fluorescence** — also the regime where the analytic
   Sobolev leg stops being exact, so a per-line Sobolev *transport* scheme
   would have to be built rather than computed. This is now Paper II, and its
   Phase 0 is done: the instrument exists and is calibrated (F14, `paper2/`).
   What it has not yet done is measure anything.
