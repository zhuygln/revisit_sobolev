# revisit_sobolev

Kilonova radiative-transfer codes cannot resolve the millions of narrow
lanthanide absorption lines that shape a kilonova's light, so they all use a
shortcut. That shortcut is usually called "the Sobolev approximation." This
project shows it is really **two** approximations with error budgets differing
by roughly an order of magnitude, and measures both against resolved
calculations on calibrated atomic data.

**The result.** Against frequency-resolved transport as truth, on realistic
La II line forests:

| treatment | band-flux error (vs a deterministic finite-profile reference on identical rays) |
|---|---|
| Sobolev approximation proper (`e^−τs` per line) | **≲1%** at v_D ≤ 30 km/s for τ_max ≤ 5, **≲0.3%** at v_D ≤ 10 km/s — a boundary effect ∝ v_D, vanishing at thermal widths |
| expansion opacity (its usual implementation) | **+38%** at thermal widths (the pure Poisson-vs-Bernoulli gap), **+45%** in the Monte Carlo implementation at 100 km/s |

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
| F15 | **The mechanism, made exact (referee revision).** Expansion opacity preserves the expected interaction count per crossing E = Σ(1−e^−τ) identically — Karp's mean free path — and applies Poisson survival e^−E to what is a Bernoulli product e^−S; F_exp/F_Sob = ⟨e^D⟩ weighted by transmission, D = S−E, to 1e-13. F4, F5, F7, F13 are corollaries | §4.18 |
| F16 | Δ_Sobolev against a deterministic reference on identical rays is +3.1% at v_D = 100 km/s in all four transport modes, +0.3% at 10, +0.03% at 1 km/s; the breadth median was a stale normalization (+3.65% → +0.26%); SEDONA resolved validates the reference to ±0.3% | §4.18 |
| F17 | Seed-matched SEDONA pairs correlate at +0.95–0.999; paired Δ_exp scatter 0.02–0.28% against 0.2–0.4% from quadrature. Headline over 10 seeds: +44.90% ± 0.04 | §4.18 |
| F18 | With radiative equilibrium on, the emergent-band differential falls from +44% to +7.7% (redistribution alone) and +5.0% (T converged); the band fills 0.34 → 0.85 and the removed flux reappears redward. Sign preserved, magnitude not | §4.18 |
| F19 | **Paper II Phase 1.** The Sobolev escape probability on re-emission, β = (1−e^−τ)/τ, is decisive for redistribution (invisible for resonant scattering); the SEDONA RE comparison caught its absence in one look | §4.19 |
| F20 | Paper I's RE fill-in is a transport-window artifact: with the whole ion's LTE emissivity the band refills only 0.348 → 0.355 and the closure differential is back to +38% | §4.19 |
| F21 | Fluorescence refills the La II band by +260% under a 6000 K continuum (Sobolev + A-branching, whole ion); an expansion-opacity code with thermal redistribution lands **−37.5%** below it — the closure is too bright in pure absorption and too faint with fluorescence | §4.19 |

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
paper2/        Paper II: Phase 0 (SEDONA audit, TARDIS record, three-level
               branching MC) and Phase 1 (whole-ion Sobolev/expansion MC
               with fluorescence, the La II measurement)
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

Phase 0 through the validity maps is complete, and the manuscript has been
revised for a major-revision report (August 2026): the thesis is now that
the two approximations carried under one name preserve different properties
of a line forest (F15), Δ_Sobolev is formed against a deterministic reference
(F16), every Monte Carlo number is a seed-matched pair (F17), and one
radiative-equilibrium check is included (F18). The response letter is
`docs/paper/response_to_referee.md`. Outstanding:

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
3. **Scattering and fluorescence** — Paper II. Phase 0 built and calibrated
   the instrument (F14); Phase 1 (`paper2/phase1/`) carries the whole La II
   ion through five treatments in one vectorized code and has measured:
   fluorescence, not thermal re-emission, refills the optical band, and an
   expansion-opacity code lands 37.5% below the physics once it does
   (F19–F21). Next: NLTE/T feedback, more ions, and a self-consistent
   incident spectrum.
