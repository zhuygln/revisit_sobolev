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

**[docs/paper/manuscript.pdf](docs/paper/manuscript.pdf)** — 21 pp in MNRAS
format, readable directly on GitHub. Written for a reader with no
radiative-transfer background: a primer builds every concept from scratch and
an appendix derives every formula.

Source is [manuscript.tex](docs/paper/manuscript.tex); rebuild with
`cd docs/paper && make` (needs pdflatex + bibtex — see
[docs/sedona/SETUP.md](docs/sedona/SETUP.md) §4; conda-forge's texlive-core
does not work).

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
| F18 | Paper I's radiative-equilibrium check is confined to SEDONA's transport window (re-emission is drawn from the emissivity on the grid), so its +7.7%/+5.0% differentials are bookkeeping within one 230 Å window, not an emergent-band result; retracted as such | §4.18 |
| F19 | **Escape-probability branching (Paper II).** A photon re-emitted in a resonance zone escapes it only with β = (1−e^−τ)/τ and is otherwise re-absorbed and re-drawn, so radiative branches compete as A·β, not A; the chain's exit distribution is A_uj β_uj/Σ A β exactly (tested). Phase 0 lacked it; the SEDONA RE comparison caught its absence in one look | §4.19 |
| F20 | **Fluorescent optical refill.** Full La II, 6000 K continuum: direct A-branching raises the 3800–3955 Å emergent flux from 0.183 to 0.660 ± 0.003; the pumps are 3300–4500 Å (the far-UV contributes nothing), through 5d6p upper levels exiting the strong ground-connected lines, 971 pathways with the top 10 carrying 25% | §4.19 |
| F21 | **Closure sign reversal.** Expansion opacity is too transparent in pure absorption (+88%) but, with complete thermal redistribution (ε = 1, re-emitting from its own κ_exp B_ν), underproduces the same band by 38% relative to direct branching; at ε = 0 it overproduces it by 34%. With line identity carried through the bin (`expansion_branch`, Poisson absorption + exact A·β exit kernel) the opacity-representation error under fluorescence is +21% (density-limited — F24) | §4.19, §4.20 |
| F22 | **No scalar ε reproduces La II fluorescence (outcome B).** Sweeping the two-level-atom thermalisation parameter on the same opacity as the physics leg, ε_best is 0.02–0.06 in the red and blue, 0.36 in the UV, and the optical 4500–6000 Å band is outside the TLA's reach at any ε; the redistribution matrices show why — branching is a blue → optical channel, thermalisation a blue → red one, and ε interpolates between them. The expansion leg's ε_best ≈ 0.3 in the 3800–3955 Å band is a compensation with its +21% opacity error and does not carry to other bands | §4.20 |
| F23 | **The closure verdict survives 0.1c; the frozen shortcut does not.** With worldline-consistent transport (per-packet clock, exact Doppler, epoch-diluted τ, aberrated re-emission) on a 0.05–0.15c shell, no scalar ε reproduces branching either (ε_best 0.20–0.49 by band, red unreachable) — but calibrated ε values shift by 0.15–0.27 vs the slow shell and the redistribution structure changes, so v_bulk is a validity-map axis; a frozen first-order snapshot overstates blue transmission 3.8× at 0.1c and is unusable for closure work | §4.21 |
| F24 | **Outcome C, and the closure's density limit.** ε_best is ion-dependent: La → Ce shifts of 0.11–0.24 per band with reachability flipping in opposite directions, and the La+Ce blend tracks the dominant (Ce) forest — composition sets the calibration. The branching-aware Poisson closure is density-limited: +21% on La II's 949-line forest, +113% on Ce II's 22,960-line one, where per-bin saturation clipping leaves the closure three orders of magnitude too transparent in a band the true opacity blacks out | §4.22 |
| F25 | **The middle method works (Paper III).** A group-to-group redistribution operator with discrete exit-frequency tables reproduces explicit A·β branching on the same opacity: La II to ≤1.6% in every band with a 4×4 matrix, Ce II to 2.9%/0.7% at 32/64 groups, bolometric to ≤0.5%, tables of tens–hundreds of kB. Within-group re-emission must be discrete: a continuous PDF double-counts self-absorption with an error that grows under refinement | §4.23 |
| F26 | **The kernel's state space is (T_gas, τ_scale, ion).** Source spectrum transfers freely (≤1.4%, 4000–8000 K); the epoch axis collapses exactly onto τ_scale (a τ-matched 1 d kernel equals each epoch's own kernel, while the fixed kernel fails at 13% where τ_max = 34); T_gas is the one genuine axis (3000→5000 K transfer errs 9.6%, recomputed ≤1.5%) | §4.24 |

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
paper3/        Paper III: the reduced redistribution closure (kernel,
               frozen reference, compression sweep, T and epoch transfer)
tests/         174 tests pinning the physics of every module
docs/          results report, lab notebook, planning inputs, paper/
data/          raw atomic data (gitignored; provenance in data/README.md)
outputs/       working figures (gitignored; committed copies in docs/figures/)
```

## Setup

Needs Python >= 3.10 and numpy >= 2, so a distro 3.8/3.9 will not do.

```bash
python3 -m venv .venv              # or: conda create -p .venv python=3.12
.venv/bin/python -m pip install -e ".[dev]" h5py
.venv/bin/python -m pytest         # 174 passed
```

Atomic data is not committed — see [data/README.md](data/README.md) for the
Zenodo record and re-download instructions. Tests that need it skip cleanly
when `data/` is empty, so a data-less clone reports passes and 5 skips.

SEDONA lives outside this repo and is needed only by `experiments/`
(everything in `paper2/` and `paper3/` is pure Python). Its no-root WSL2
build recipe is in the lab notebook §5, with the working `Makefile.wsl` and
the full machine-setup walkthrough — Python, data, SEDONA, LaTeX — in
**[docs/sedona/SETUP.md](docs/sedona/SETUP.md)**.

The SEDONA experiment drivers find the binary through `SEDONA_HOME`
(default `~/personal/pubsed`) or `SEDONA_EXE`.

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
   ion through eight treatments in one vectorized code: fluorescence, not
   thermal re-emission, refills the optical band, and an expansion-opacity
   code with complete thermalisation lands 38% below the physics once it
   does (F19–F21). Phase 2 answered the decisive question: **no scalar
   thermalisation parameter ε reproduces direct La II fluorescence** across
   the spectrum (F22) — ε_best differs by band and one band is unreachable —
   whereas a closure that keeps line identity through the bin gets within
   21%. The verdict survives realistic bulk velocity: at 0.05–0.15c with
   worldline-consistent transport no scalar ε works either, while the
   calibrated ε values shift by 0.15–0.27 — and a frozen-snapshot code at
   0.1c overstates blue transmission 3.8× (F23). Ce II makes it outcome C:
   ε_best is ion-dependent, the La+Ce blend tracks the dominant forest, and
   the branching-aware closure itself is density-limited (+21% on La II,
   +113% on Ce II, whose 22,960-line forest the Poisson opacity cannot
   black out) (F24). Paper II's identity is that result; manuscript text
   and NLTE come afterwards. Paper III (`paper3/`) opens the constructive
   question — how much redistribution information is actually required? —
   and its first gate is passed: a discrete-table R_ij reproduces both
   forests (F25), La II at 4–8 groups and Ce II at 32–64 — and the
   kernel's state space is just (T_gas, τ_scale, ion): the source spectrum
   transfers freely and the epoch axis collapses onto τ_scale (F26).
