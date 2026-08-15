# revisit_sobolev

How large is the error introduced by the **Sobolev approximation** in r-process
(kilonova) ejecta, and where in kilonova parameter space does that error matter?

## Phase 0 status: COMPLETE — GO

All four weekend deliverables from [docs/babystep_plan.md](docs/babystep_plan.md) are done:

- **0A** — `tau_exact` reproduces the analytic Sobolev optical depth to machine
  precision (E_Sob = 2×10⁻¹⁶) with a monotone convergence study.
- **0B / Figure 1** — a tanh population gradient breaks the approximation on a clean
  ε² power law (ε = v_D/v_scale), reaching ~5% error at ε = 1. Finding: a gradient
  symmetric about the resonance cancels at leading order, so the resonance must sit
  off-centre — gradient-only validity diagnostics are conservative.
- **0C / Figure 2** — real GSI La II data (3000 Å–3 μm): 11.3% of the broader subset
  has a neighbour within 10 km/s (3.3% within 3 km/s); even the high-confidence
  subset has a 1.2% tail. The thermal-width crowding regime exists.
- **0D** — SEDONA (public `dnkasen/pubsed`) built serially-invoked under WSL2 (deps
  via a micromamba env + Lua 5.1.5 source build + OpenMPI, custom
  `src/makefiles/Makefile.wsl`); the 1-D `spherical_lightbulb` example propagates
  2×10⁵ particles in ~1 s and writes a spectrum. Both line treatments needed later
  exist in this release: `opacity_line_expansion` (Sobolev) and
  `opacity_bound_bound` (resolved), `defaults/sedona_defaults.lua:115,118`.

## Week 1 status (babystep_plan.md §18)

- **Voigt profile** — implemented via the Faddeeva function; normalization,
  Gaussian-limit, and Lorentzian-wing tests.
- **Boltzmann populations** — partition function and level fractions straight from
  the GSI levels files (parsed by the same `load_gsi`); analytic two-level and
  temperature-limit tests, plus the per-statistical-weight monotonicity invariant
  on real La II levels.
- **SEDONA Type Ia pair** (`param_d20_lte_exp` vs `param_d20_lte_bb`, run in
  `~/personal/pubsed/examples/supernova/TypeIa/spectrum/run_{exp,bb}/`): the two
  configs differ **only** in the line treatment — the built-in controlled
  experiment this project needs. Results (Figure 3, `outputs/`):
  - cost: **14 s (Sobolev) vs 88 min (resolved)** for the same 4 iterations —
    a ~378× ratio, empirically confirming the §11 frequency-resolution cost and
    the necessity of narrow wavelength windows;
  - physics: the resolved treatment blankets the UV harder and redistributes
    ~3× more flux into the NIR (8000–12000 Å) — qualitatively the effect that
    matters for lanthanide-rich kilonova ejecta (toy model; qualitative only).

## Week 2 status (babystep_plan.md §19)

- **Formal transfer solver** (`sobolev/formal_transfer.py`) — deterministic
  impact-parameter ray integration through a homologous shell around a blackbody
  core; validated against the analytic blackbody-sphere luminosity, the Sobolev
  trough `e^-τs`, two-line independence, and emission fill-in.
- **Minimal 1-line cross-validation** (`experiments/minimal_1line/`, Figure 4):
  fake Ly-α (SEDONA's 2-level atom), τ_S = 2, neutral 2000 K shell. Trough
  depths vs analytic 0.1353:
  - Python solver **0.1372**, SEDONA resolved bound-bound **0.1420** — the two
    independent codes agree across the whole profile;
  - SEDONA expansion opacity **0.4285** = the predicted single-line failure
    `exp(-(1-e^-τs))` = 0.4212. **Expansion opacity ≠ Sobolev line transfer for
    a single strong line** — the statistical treatment attenuates by the
    bin-averaged opacity. First quantitative discrepancy of the project.
- **Documented systematic**: the observer-frame formal solution and the
  comoving-frame Sobolev formula differ at O(v_bulk/c) (7% at 0.077c, matching
  `exp(-τ_S(1−z/ct))` exactly). Week 3 comparisons must place resonances at low
  bulk velocity or match frame conventions.

- **2 → 20 line ladder** (`experiments/line_ladder/`, Figure 5): custom SEDONA
  atoms (2 lines at 1500 km/s; 20-line forest at 750 km/s), τ_S = 0.5 per line.
  Forest-center depth: analytic Sobolev **0.2405**, Python solver **0.2407**,
  SEDONA resolved **0.2435** — locked within ~1%, sawtooth reproduced tooth for
  tooth. SEDONA expansion opacity **0.3257** vs its own analytic prediction
  0.3194: the per-line `(1−e^−τ)` substitution error is per-resonance and does
  **not** average away with line count — only as τ→0. Expansion opacity is a
  weak-line approximation, not a many-line one.

Week 2 is complete.

**Detailed documentation:** [docs/results_report.md](docs/results_report.md)
(methods, all experiments with figures, findings register F1–F6, caveats,
reproduction guide) and [docs/lab_notebook.md](docs/lab_notebook.md)
(chronological log including dead ends and fixes).

## Weeks 3–4 status (babystep_plan.md §20)

- **Real GSI La II forest** (`experiments/laII_forest/`, Figure 6): the
  3850–3950 Å window — 153 calibrated lines, 3 with τ_S > 1 — at T = 3000 K,
  day 1, v_D = 100 km/s. The SEDONA atom file is generated *from* the GSI data
  (all 472 levels, ionization disabled) so both codes carry byte-identical
  levels, populations, and line strengths. Band-averaged L/L_cont:
  Python solver **0.355**, SEDONA resolved **0.343** (3.6% apart);
  SEDONA expansion **0.497** — **Δ_Sob = +45% band flux**, the first
  validity-map datapoint. Real forests fail *worse* than the uniform ladder:
  a few strong lines dominate, and `(1−e^−τ)` caps each crossing at 1 while
  the true attenuation saturates to black.

- **First validity-map slice** (`sweep.py`, Figure 7): 24 SEDONA runs,
  τ_max ∈ {0.5, 5, 50} × v_D ∈ {10, 30, 100, 300} km/s. Δ_Sob decomposes into
  **(1)** a v_D-independent floor set by line strength — the per-crossing
  `(1−e^−τ)` substitution: ~3% at τ_max = 0.5, ~40% at 5, ~45% at 50 — and
  **(2)** a wing term growing with v_D (expansion opacity is width-blind while
  resolved profiles absorb continuum as they widen): +91% at (τ=50, 300 km/s).
  Practical reading: the error does **not** vanish as v_D → thermal — it
  converges to the strength floor.

- **T axis + thermal-width frontier** (`tsweep.py`, Figure 8): T is a weak
  axis in this window (+41→+47% over 2500–5000 K at pinned τ_max = 5 — the
  same strong lines dominate throughout). The frontier confirms the floor:
  **Δ_Sob ≈ 38% flat from 10 km/s down to 1 km/s**, essentially at the La
  thermal width (0.6 km/s). Cost scales ≈ 1/v_D for *both* modes (~35 min/run
  at 1 km/s) — fine transport grids, not the resolved treatment, drive it.

- **Multi-ion blend** (`experiments/multiion/`, Figure 9): La II + Ce II,
  2,529 lines (Ce II floods the window), strong-line spacing down to 8 km/s.
  The blend saturates the window to near-black in every treatment and
  **Δ_Sob drops to +14.6%** — the error is *non-monotonic* in forest density
  (F7): maximal for strong sparse forests, suppressed at full blanketing.
  The dangerous regime is the intermediate one, where spectral features form.

- **Resolved-legs gap closed** (F8): the solver-vs-SEDONA offset was neither
  profile wings (3×10⁻⁵ effect) nor resolution (converged), but **shell
  thermal emission** — the solver integrates S = B_ν(T_shell), SEDONA at
  fixed T deposits absorbed energy without re-emitting. Like-for-like the two
  independent codes agree to **~1%**, at the Monte Carlo noise floor, in both
  the single-ion and blended forests. Δ_Sob is unaffected (same-code
  differential).

Next: other windows/epochs; a true per-line Sobolev leg to separate
approximation error from expansion-opacity implementation error; fold F7–F8
into the manuscript.

---

This repo started as a **Phase 0** weekend feasibility check, not a paper: reproduce
the Sobolev limit in a toy problem, break it deliberately, and look at whether real
lanthanide data occupy the failure regime. Everything else is deferred.

## Layout

```
sobolev/     toy solver: constants, profiles, optical depths, atomic data
notebooks/   00 single line (Phase 0A+0B), 01 GSI line spacing (Phase 0C)
tests/       profile normalization + the Phase 0A Sobolev-limit gate
data/        raw atomic data, gitignored except its README
docs/        the two planning documents
outputs/     figures (gitignored)
```

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

All 7 tests pass. The La II data files are not committed — see `data/README.md` for
the Zenodo provenance; re-download and drop them in `data/` to run notebook 01.

Go/no-go after each phase. See [docs/babystep_plan.md](docs/babystep_plan.md) for the
full roadmap and [docs/research_requirements.md](docs/research_requirements.md) for the
data and code landscape.
