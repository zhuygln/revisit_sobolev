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
| F26 | **The kernel's state space is (T_gas, τ_scale, ion) — for La II.** Source spectrum transfers freely (≤1.4%, 4000–8000 K); the epoch axis collapses exactly onto τ_scale (a τ-matched 1 d kernel equals each epoch's own kernel, while the fixed kernel fails at 13% where τ_max = 34); T_gas is the one genuine axis (3000→5000 K transfer errs 9.6%, recomputed ≤1.5%). The source-spectrum half does not generalize — see F28 | §4.24, §4.25 |
| F27 | **Gate 2: compression is generic (outcome A).** A discrete-table R_ij reproduces explicit branching for all three ions on the same opacity: La II 1.62% and Ce II 7.48% at N_g = 4, Ce reaching 2.89%/0.74% at 32/64 — and **Nd II 0.44% at N_g = 4**, the easiest of the three despite 3,336,077 transitions, because the fixed τ_max = 5 normalization leaves it only 4,496 opacity lines against Ce's 22,960. The blue→blue block separates the hard ion from the easy ones but does not order them | §4.25 |
| F28 | **The kernel's state space is itself ion-dependent.** F26's structural claims carry to Nd II — the epoch axis still collapses onto τ_scale exactly, T_gas is still the genuine axis (6.4–19.2% fixed, ≤0.51% recomputed). Its third claim does not: Nd's fixed-kernel error is monotone in T_src and flips sign across the training point (+12.44 → −4.88% over 4000–8000 K) where La's is flat at ~1%. Denser groups draw a more evenly weighted absorbing-line mix, so reweighting the continuum reweights the rows | §4.25 |
| F29 | **The composition rule works at the 5% level (P9).** An opacity-weighted mixture of per-ion kernels — weights from the blend's opacity alone, no blend run — reaches 4.27% worst-band on La+Ce at N_g = 64 where a blend-trained kernel gets 1.37%, beating the best single-ion control by 2.4×. The gain is located: `ce_only` fails at −10.3% in the optical, the blue → optical branching channel that La carries despite being 5% of the opacity, and the rule repairs it to −0.7%. Composition therefore leaves the state space at the 5% level — a per-ion library suffices — while explicit blend training is still needed below ~2% | §4.26 |
| F30 | **The opacity is the binding constraint, not the redistribution (P11).** With the identical R_ij, exact line opacity errs 0.92% (La II) / 2.21% (Ce II); grouping the opacity takes that to 14–18% / 91–127% whichever single-scalar rule is used — Στ too opaque on La (−14.5%), Σ(1−e^−τ) too transparent (+17.8%), both far too transparent on Ce. F15 as a design constraint: one scalar per bin cannot carry both the attenuation and the interaction count, and scattering needs both. The κ_grouped + R_ij target architecture fails on dense forests; a bin carrying both is the candidate A bin carrying **both** quantities (survival from S, line draw from p) was then tested: it works on La II (21.32% → **8.66%**, saturated band +21.3% → −0.7%) and fails on Ce II (112.86% → **139.27%**), where more opacity makes the band *brighter* because it is refilled by fluorescence faster than absorbed — redistribution-limited, not attenuation-limited. And `dual_group` is bit-identical to `binned_group`: a pure R_ij closure never draws a line in the bin, so it cannot use the second quantity at all | §4.27 |
| F31 | **One remembered line is not the missing state — and the failure splits by line spacing.** Carrying exactly one extra number per packet (the frequency last emitted at, crediting that line's τ to the next free path) buys a factor 2.2 on La II (14.49% → **6.50%**, saturated band −14.5% → −6.3%), the cheapest grouped-opacity repair found. On Ce II it moves 126.66% → 116.31% and 91.29% → 81.66% — real but no rescue. The events/packet counter shows memory removes a comparable fraction of the excess interactions in both forests (39% La, 34% Ce), so the Ce error is not driven by excess interactions: it is the fluorescent refill of §4.19–4.20, and no local interaction bookkeeping reaches it. Sparse forests need one remembered line; ~~dense forests need the resonance *sequence*~~ — that second clause is **retracted by F33**, which finds no benefit from depth beyond m = 1 on Ce | §4.28, §4.30 |
| F32 | **The redistribution operator is local in frequency, not low-rank — so "few modes" is the wrong explanation for its compressibility.** Effective dimension never saturates: participation ratio grows as N_g^0.64–0.75 across La/Ce/Nd with PR/N_g falling 0.5 → 0.2, and NMF rank-8 of a 25-row operator still misses row-L1 0.47 (23% total variation per row). Rank *anti*-correlates with compressibility — Ce has the lowest energy-operator dimension (PR 1.60) and is the hardest to compress. Decisive at matched parameter count: 16 numbers as a coarse 4×4 matrix give 1.62% where 912 numbers as a rank-16 factorization give 11.30%. Coarsening averages neighbouring groups; truncation projects onto modes; only the first works. This unifies Paper III — redistribution is smooth at the group scale so it coarse-grains (F25/F27), the opacity is a comb whose ordering decides a packet's fate so it does not (F30/F31) | §4.29 |
| F33 | **Resonance-sequence depth is not the missing information, and the density limit is not about density (retracts F31's second clause).** Making memory a depth rather than a switch: La converges by m = 4 and Nd by m = 8, each gaining ≤0.5 points beyond m = 1, and **Ce gains nothing** — 116.31 → 116.79% across a sixteen-fold increase in remembered history. Memory is a between-step correction; the dense-forest failure is within-step. Nd II, run through P11 for the first time, breaks the density reading outright: 4.7× La's opacity lines and **1/12 its error** (expansion + A·β 1.79%, a working closure). What orders the three ions is **band-local saturation** — saturated lines inside the failing band (Nd 1, La 4, Ce 24), or Στ there (11.4, 22.4, 89.6) — not total line count. Three points cannot fix the exponent; that is what the synthetic phase diagram is for | §4.30 |
| F34 | **Band-local saturation controls the grouped-closure failure; redistribution does not.** 96 synthetic forests with independently dialled crowding, saturation, spacing and redistribution range: Spearman ρ = +0.91 for Στ in the band and +0.86 for N_sat there, against **+0.25 and −0.31 for the two redistribution axes** — an independent confirmation of F33 from the opposite direction. The synthetic family collapses as ΔF = 0.162·N_sat^0.58 (scatter ×1.95), and the real atoms sit at ratio 0.71 / 0.40 / 1.25, two of three inside that scatter. Building the forests required first measuring the real τ distribution *inside* the failing band: mostly weak lines with a saturated tail, ln-spread 1.7–2.05 across all three ions. Band-to-forest geometry is eliminated as the residual cause (10× change → 16% effect); the emergent-cascade `ladder` forests match the real interpolation at matched N_sat (59–62% vs 62%) where the dialled ones overshoot (78–94%). A partial collapse: not yet a general law | §4.31 |

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
   and both its gates are passed. A discrete-table R_ij reproduces
   explicit branching on the same opacity (F25), and it does so for every
   ion tried: La II and Nd II at four groups, Ce II at 32–64, which is
   Gate 2's outcome A (F27). The epoch axis collapses exactly onto
   τ_scale and T_gas is the one genuine thermodynamic axis (F26), both
   ion-independent — but whether the source spectrum can be dropped is
   not: it transfers freely for La and fails at 12% for Nd (F28), so
   T_src must be checked per ion before a kernel is tabulated. Composition,
   by contrast, does leave the state space at the 5% level: an
   opacity-weighted mixture of per-ion kernels reaches 4.27% on the La+Ce
   blend with no blend training, beating the dominant ion alone by 2.4×
   (F29). P11 then put the two halves together and returned a negative
   result worth having: with the identical kernel, exact line opacity errs
   ≤2.2%, but grouping the opacity costs 14–18% on La II and 91–127% on
   Ce II by either single-scalar rule (F30). The redistribution half was
   never the hard half — one scalar per bin cannot carry both the
   attenuation and the interaction count, and scattering needs both. So
   κ_grouped + R_ij is not usable on dense forests. A bin carrying both
   quantities was then tried: it works on La II (21.3% → 8.7%) and fails
   on Ce II (112.9% → 139.3%), where the deep band is redistribution- not
   attenuation-limited — and a pure R_ij closure cannot use the second
   quantity at all, since it never draws a line in the bin. Restoring
   line identity at emission is the remaining lever, and it costs the
   thing grouping was meant to buy. E1 then asked why the redistribution
   half compresses at all, and the answer is not "few modes": the
   operator's effective dimension never saturates (PR ~ N_g^0.65) and
   rank truncation fails badly, while 16 numbers as a coarse matrix beat
   912 as a rank-16 factorization (F32). Redistribution compresses
   because it is *smooth* in frequency, not low-rank — and the opacity
   does not because a comb of resonances is not smooth at any group
   scale. That is one statement covering both halves.
