# Sobolev Validity in Kilonova Ejecta — Results Report

**Status:** Phase 0 through the first Weeks 3–4 validity-map slice, complete.
**Date:** 2026-08-14. **Repo:** `zhuygln/revisit_sobolev`.
Companion document: [lab_notebook.md](lab_notebook.md) (chronological log,
including dead ends and fixes). Planning documents:
[babystep_plan.md](babystep_plan.md), [research_requirements.md](research_requirements.md).

---

## 1. Research question

> How large is the error introduced by the Sobolev approximation in r-process
> (kilonova) ejecta, and where in parameter space does it matter?

The strategy (per the planning documents) is not to build a kilonova pipeline,
but to construct a small, fully-controlled harness in which the *same* model,
atomic data, and populations are pushed through resolved-profile and
Sobolev-class line treatments, and the difference is measured.

## 2. The three-way harness

Every experiment from Week 2 onward runs the same physical setup through three
independent calculations:

| Leg | Code | Method |
|---|---|---|
| Analytic | pencil & paper (+ numpy) | Sobolev optical depths, p-averaged over the core disk |
| Deterministic | `sobolev/formal_transfer.py` | formal solution of dI/ds = −αI + j along impact-parameter rays |
| Monte Carlo | SEDONA (public `dnkasen/pubsed`) | resolved bound-bound **and** expansion-opacity modes |

Agreement of the first two with SEDONA's *resolved* mode validates the
harness; the gap to SEDONA's *expansion-opacity* mode is the measured
Sobolev-class error, Δ_Sob.

### 2.1 Geometry and conventions

- 1-D homologous shell (v = r/t) around a blackbody core; rays with
  p < r_core start at the core with I = B_ν(T_core); envelope rays start at 0.
- Pure-absorption LTE: S = B_ν(T_shell), `opacity_epsilon = 1`, no electron
  scattering, no radiative equilibrium (T fixed at the model value).
- Shell temperatures are chosen so B_ν(T_shell) ≈ 0 in the band — the line
  acts on the core continuum only, which all three legs treat identically.
- Comoving frequency to first order: ν′ = ν(1 − z/ct); resonance surfaces are
  z = const planes for every ray (a special property of homologous flow).

### 2.2 Population control

SEDONA computes LTE (Saha–Boltzmann) populations internally. To guarantee
*identical* populations across legs (research_requirements.md §17), the
SEDONA atomic-data files are **generated from the same source data the Python
side uses**, with ionization disabled (χ = 10⁵ eV), so both codes run
Boltzmann over byte-identical level lists. Verified in the SEDONA gas state
(e.g. n_e/n_H = 10⁻¹⁰ in the minimal model).

## 3. Infrastructure

### 3.1 Python package (`sobolev/`)

| Module | Contents | Key validation |
|---|---|---|
| `constants.py` | CGS constants; σ_cl = πe²/m_ec **derived, not hardcoded** | an earlier hardcoded value was wrong by 18% — caught by re-deriving |
| `profiles.py` | Gaussian + Voigt (Faddeeva/`wofz`) | ∫φdν = 1 to 1e-6; Voigt→Gaussian at γ=0 to 1e-12; Lorentzian wing amplitude analytic |
| `optical_depth.py` | resolved τ_exact, Sobolev τ_S, tanh gradient sweep | τ_exact → τ_S to 2×10⁻¹⁶; ε² breakdown law |
| `populations.py` | partition function, Boltzmann fractions from GSI levels | two-level closed form; T→0/∞ limits; per-g monotonicity on real La II |
| `atomic_data.py` | GSI file parser (transitions *and* levels), spacing stats | format regression on committed 20-row excerpts |
| `formal_transfer.py` | the deterministic solver (§2) | blackbody-sphere luminosity exact; e^−τs trough; 2-line independence & blend; emission fill-in |

Test suite: **19 tests, all green** (`pytest`).

### 3.2 SEDONA build (WSL2, no root)

Public release `github.com/dnkasen/pubsed`, built via
`src/makefiles/Makefile.wsl` (kept in the clone, not this repo):
GSL + HDF5 + OpenMPI from a micromamba env at `~/personal/sedona-deps/env`,
Lua 5.1.5 built in its source tree, `OMPI_CXX=g++` to point the conda MPI
wrapper at the system compiler. `sedona.h` hardcodes `MPI_PARALLEL`, so MPI is
mandatory; single-rank runs need no `mpirun`. Examples need
`SEDONA_HOME=~/personal/pubsed` and — important — the **full environment**
(a sanitized env fails with rc=14 from the OpenMPI runtime).

### 3.3 Atomic data

GSI Database for Kilonova Radiative Transfer (Flörs et al.), Zenodo record
**19335084** (latest under concept DOI `10.5281/zenodo.15835360`; paper DOI
`10.1103/jxqw-7ynk`), CC-BY 4.0. La II: 472 levels, 17,743 E1 transitions,
calibration flags `xmatch`/`shifted`. Raw files gitignored under `data/`
(provenance table in `data/README.md`); 20-row excerpts committed under
`tests/data/` for format regression. Internal consistency check: f from the
`A` column vs f from `Log(gf)` agree to 1.5×10⁻⁴ across the study window.

## 4. Experiments and results

### 4.1 Phase 0A — Sobolev limit reproduced (Figure 1, left half)

Single line, constant population, homologous slab. The resolved integral
τ_exact converges to the analytic τ_S:

| grid points | E_Sob |
|---|---|
| 201 | 8.2×10⁻¹ |
| 2,001 | 1.0×10⁻⁴ |
| 20,001 | 2.2×10⁻¹² |
| 200,001 | **2.0×10⁻¹⁶** (machine precision) |

The convergence study is itself a test: an unconverged grid mimics a physical
Sobolev deviation.

### 4.2 Phase 0B — controlled breakdown (Figure 1)

![Figure 1](figures/fig1_sobolev_error_vs_gradient.png)

tanh population gradient of scale v_scale; control parameter
ε = v_D/v_scale. Measured E_Sob follows an exact ε² power law (log-slope
2.0000) from 10⁻⁴ to ~0.3, reaching **4.8% at ε = 1**.

**Finding F1 (symmetry cancellation).** With the resonance at the tanh
*center*, the error cancels identically — the Gaussian is even and the tanh
odd (first order), and the center is the inflection point where n″ = 0
(second order). The experiment only works with the resonance placed
off-center (u = 0.5 used throughout). Consequence: validity diagnostics built
on the gradient magnitude G = v_D|d ln n/dv| alone are *conservative* — the
true error depends on where the resonance sits on the gradient.

### 4.3 Phase 0C — line crowding in real data (Figure 2)

![Figure 2](figures/fig2_laII_line_spacing.png)

Nearest-neighbour velocity spacings of GSI La II lines, 3000 Å–3 μm window:

| subset | lines | P(Δv<1 km/s) | <3 | <10 |
|---|---|---|---|---|
| high confidence (both levels `xmatch`, log gf > −1) | 321 | 0.0% | 0.0% | 1.2% |
| broader (`xmatch`/`shifted`, log gf > −2) | 7,950 | 1.0% | 3.3% | 11.3% |

The thermal-width crowding regime **exists** in calibrated single-ion data;
multi-ion blends only add lines. (Crowding is necessary, not sufficient, for
transport consequences — the harness measures the consequences.)

### 4.4 Phase 0D / Week 1 — SEDONA runs (Figure 3)

![Figure 3](figures/fig3_sedona_typeIa_exp_vs_bb.png)

`spherical_lightbulb` smoke test: 2×10⁵ particles, 1.1 s. Type Ia toy pair
(`param_d20_lte_exp` vs `param_d20_lte_bb`, identical but for the line
treatment):

**Finding F2 (cost).** Same model, same 4 iterations: expansion opacity
**14 s**, resolved bound-bound **88 min** — a **378×** ratio, empirically
confirming the frequency-resolution cost analysis (babystep_plan.md §11) and
the necessity of narrow wavelength windows.

Physics: the resolved treatment blankets the UV harder and moves ~3× more
flux into 8000–12000 Å — the UV→NIR redistribution mechanism relevant to
lanthanide-rich ejecta (toy model; qualitative).

### 4.5 Week 2 — solver validation and the minimal model (Figure 4)

![Figure 4](figures/fig4_minimal_1line_threeway.png)

Minimal 1-element/1-ion/1-line model: SEDONA's shipped 2-level atom (fake
Ly-α, f = 0.6647), τ_S = 2, neutral 2000 K shell. Trough depths against
analytic e^−τs = 0.1353:

| leg | trough |
|---|---|
| Python formal solver | 0.1372 |
| SEDONA resolved bound-bound | 0.1420 |
| SEDONA expansion opacity | **0.4285** |

**Finding F3 (frame systematic).** The observer-frame formal solution and the
comoving-frame Sobolev formula differ at O(v_bulk/c): a resonance at
z/ct = 7.7% produced a trough matching exp(−τ_S(1−z/ct)) to four digits.
Documented in `formal_transfer.py`; tests place resonances at v/c < 1%.
Cross-code comparisons must control this (low-velocity resonances, or matched
frame conventions).

**Finding F4 (expansion opacity ≠ Sobolev line transfer).** The
expansion-opacity trough is not e^−τs but **exp(−(1−e^−τs)) = 0.4212**
(measured 0.4285): a photon crossing one resonance is attenuated by the
bin-averaged opacity, capped at one "effective τ" per crossing. For a single
strong line the two differ by 3.2× in transmitted flux. "Sobolev mode" in
transport codes usually means expansion opacity — the distinction matters.

### 4.6 Week 2 — the 2 → 20 line ladder (Figure 5)

![Figure 5](figures/fig5_line_ladder.png)

Custom GSI-format SEDONA atoms: 2 lines at 1500 km/s separation; 20-line
forest at 750 km/s; τ_S = 0.5 per line. Forest-center depth, five ways:

| leg | depth |
|---|---|
| analytic Sobolev, p-averaged | 0.2405 |
| Python formal solver | 0.2407 |
| SEDONA resolved | 0.2435 |
| analytic expansion prediction | 0.3194 |
| SEDONA expansion | 0.3257 |

The resolved trio locks to ~1% and reproduces the sawtooth tooth-for-tooth.

**Finding F5 (the error is per-resonance).** Expansion opacity follows its
own per-line (1−e^−τ) substitution: the error does **not** average away with
line count — only as τ → 0. Expansion opacity is a *weak-line* approximation,
not a many-line one.

Technical notes that will save future time: (i) the analytic staircase must
be p-averaged over the core disk — resonance planes at z < r_core are crossed
by only the outer rays, and plane-counting without the average is visibly too
shallow; (ii) SEDONA atom files require the group attributes
`n_ions/n_levels/n_lines` — the reader segfaults without them.

### 4.7 Weeks 3–4 — real La II forest (Figure 6)

![Figure 6](figures/fig6_laII_forest.png)

Window 3850–3950 Å (auto-selected for the best τ distribution): 153
calibrated La II lines; τ_S > 1: 3 (max 5.0), 0.1–1: 6. Conditions
T = 3000 K, day 1, shell 1000–3000 km/s, v_D = 100 km/s, n_ion = 2146 cm⁻³
(ρ = 5.0×10⁻¹⁹ g/cm³). Band-averaged L/L_cont over 3800–3955 Å:

| leg | band flux |
|---|---|
| Python formal solver | 0.3549 |
| SEDONA resolved | 0.3426 (3.6% from solver) |
| SEDONA expansion | 0.4965 (**Δ_Sob = +44.9%**) |

The resolved complex saturates to black across 3812–3860 Å; expansion opacity
never drops below ~0.2. Real forests fail *worse* than the uniform ladder
(+45% vs +35%) because they are strong-line dominated and (1−e^−τ) caps each
crossing.

### 4.8 Weeks 3–4 — first validity-map slice (Figure 7)

![Figure 7](figures/fig7_validity_slice.png)

24 SEDONA runs: τ_max ∈ {0.5, 5, 50} × v_D ∈ {10, 30, 100, 300} km/s, same
window and conditions. Δ_Sob = (F_exp − F_bb)/F_bb (same-code differential,
immune to cross-code systematics):

| τ_max \ v_D | 10 | 30 | 100 | 300 km/s |
|---|---|---|---|---|
| 0.5 | +3.0% | +2.3% | +2.9% | +2.8% |
| 5 | +37.3% | +40.0% | +43.3% | +51.2% |
| 50 | +44.6% | +47.4% | +58.8% | +91.1% |

**Finding F6 (two-component error).** Δ_Sob decomposes into:
1. a **v_D-independent floor set by line strength** — the per-crossing
   substitution: ~3% at τ_max = 0.5, ~40% at 5, ~45% at 50;
2. a **wing term growing with v_D** — expansion opacity is width-blind (flat
   in v_D at every strength) while resolved profiles absorb continuum
   between lines as they widen: at (τ=50, 300 km/s) the error doubles
   to +91%.

Practical reading: pushing the resolved calculation toward true thermal
widths does **not** rehabilitate expansion opacity — the error converges to
the strength floor, not to zero. Lanthanide-rich kilonova ejecta live in the
strong-line regime, where the floor dominates.

### 4.9 Weeks 3–4 — temperature axis and thermal-width frontier (Figure 8)

![Figure 8](figures/fig8_T_and_frontier.png)

**Temperature axis** (τ_max pinned at 5 by rescaling n_ion at each T, so the
measurement isolates population *redistribution* from the strength scale):

| T [K] | n_ion [cm⁻³] | lines τ>1 / 0.1–1 | Δ_Sob |
|---|---|---|---|
| 2500 | 1866 | 3 / 5 | +40.8% |
| 3000 | 2146 | 3 / 6 | +43.6% |
| 4000 | 2647 | 3 / 6 | +46.7% |
| 5000 | 2934 | 4 / 5 | +46.2% |

T is a **weak axis** for this window: the same few strong lines dominate at
every temperature, so redistribution barely moves the error (±3% over a
factor-2 range in T). Caveat: at 5000 K the shell's thermal emission at
3900 Å reaches ~30% of the core surface brightness, so absolute fluxes
include fill-in; Δ_Sob remains a same-code differential.

**Thermal-width frontier** (T = 3000 K, τ_max = 5):

| v_D [km/s] | Δ_Sob | resolved wall time |
|---|---|---|
| 10 | +37.3% | ~25 s |
| 3 | +38.1% | 719 s |
| 1 | +37.8% | 2109 s |

**Finding F6 confirmed to 1 km/s** — within MC noise of the La thermal width
(0.6 km/s): Δ_Sob is flat below 10 km/s at the strength-set floor of ~38%.
The wing term is gone; the per-crossing substitution error remains in full.
The cost curve scales ≈ 1/v_D (transport bins), reaching ~35 min per run at
1 km/s — and the *expansion* runs pay the same grid cost (2305 s at 1 km/s):
fine grids, not the resolved mode per se, drive the expense.

### 4.10 Weeks 3–4 — multi-ion blend: La II + Ce II (Figure 9)

![Figure 9](figures/fig9_multiion_forest.png)

Two-element forest in the same window (equal mass fractions; total density
set so the strongest line of either species has τ_S = 5; Ce III extracted but
deliberately excluded — two stages of one element would hand the II/III split
to Saha and break exact population control). Ce II floods the window:
**2,376 lines** vs La's 153 (2,529 total; 40 with τ > 0.1; minimum spacing
between strong lines **8 km/s** — genuine sub-Doppler cross-species
blending). Band-averaged L/L_cont:

| leg | band flux |
|---|---|
| Python solver | 0.2426 |
| SEDONA resolved | 0.2277 (6.5% from solver) |
| SEDONA expansion | 0.2610 (**Δ_Sob = +14.6%**) |

**Finding F7 (non-monotonic density dependence).** The blend saturates the
whole window to near-black in *every* treatment — with enough crossings,
even capped per-crossing contributions sum to blackness — so the relative
band-flux error *shrinks* to +14.6%. Δ_Sob is non-monotonic in forest
density: ~τ²-small for weak forests, maximal (+45%) for strong sparse
forests, declining again toward full line blanketing. The dangerous regime
for expansion opacity is the *intermediate* one — strong lines that are not
yet fully blanketed — which is exactly where spectral features form.

Caveats: the resolved-legs gap widened to 6.5% in the dense blend,
plausibly the profile-wing difference (SEDONA's Voigt with 2,376 lines'
cumulative Lorentzian tails vs the solver's pure Gaussian) — a physical
setup difference to be resolved by adding Voigt to the solver. Data note:
Ce II carries half-integer J as fraction strings ('7/2'), now handled by
`sobolev.populations.parse_j`.

## 5. Findings register

| # | Finding | Where |
|---|---|---|
| F1 | Symmetric gradients cancel the Sobolev error at leading order; G-based diagnostics are conservative | §4.2 |
| F2 | Resolved vs expansion cost: 378× on the Type Ia toy | §4.4 |
| F3 | O(v_bulk/c) frame systematic between observer-frame integration and comoving Sobolev τ | §4.5 |
| F4 | Expansion opacity attenuates by exp(−(1−e^−τ)) per crossing, not e^−τ | §4.5 |
| F5 | That error is per-resonance: it does not average away with line count, only as τ→0 | §4.6 |
| F6 | Δ_Sob = strength-set floor + v_D-growing wing term; error survives the v_D→0 limit | §4.8 |
| F7 | Δ_Sob is non-monotonic in forest density: maximal for strong sparse forests, suppressed at full blanketing | §4.10 |

## 6. Caveats and limitations

- **Scope of "Sobolev mode":** SEDONA's expansion-opacity implementation is
  the measured Sobolev-class treatment. A direct per-line Sobolev interaction
  scheme (as in TARDIS-style codes) would show different (likely smaller)
  errors on F4/F5; the harness can be extended to test one.
- **Pure absorption only:** ε = 1, no scattering, no fluorescence, no NLTE.
  These are deferred by design (babystep_plan.md §9) — differences measured
  here are attributable to line transfer alone.
- **Controlled geometry:** shell velocities 1000–3000 km/s keep the O(v/c)
  frame systematic (F3) below the measurement level. Realistic kilonova
  velocities (0.1–0.3c) require frame-consistent treatment before the
  comparison is repeated there.
- **v_D = 10–300 km/s, not thermal:** the real La thermal width is
  ~0.6 km/s. F6 argues the floor persists in that limit, but the wing term at
  thermal widths remains to be measured (the expensive frontier).
- **Single ion, single window, single (T, t):** one slice of the eventual
  map. T sweeps (populations), multi-ion overlap, and other windows are next.
- **MC noise:** ~1–2% per SEDONA band flux at 2×10⁶ core particles; the
  resolved-trio agreements (1–4%) are at or near this floor.

## 7. Reproduction

```bash
# environment
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]" h5py
pytest                    # 19 passed

# data (once): Zenodo 19335084 -> data/, see data/README.md
# SEDONA (once): see lab_notebook.md "SEDONA build" entry

# notebooks (Figures 1-2)
jupyter nbconvert --execute --to notebook --inplace notebooks/*.ipynb

# minimal 1-line (Figure 4)
cd experiments/minimal_1line && python make_model.py
(cd run_bb  && SEDONA_HOME=~/personal/pubsed ~/personal/pubsed/src/sedona6.ex param.lua)
(cd run_exp && SEDONA_HOME=~/personal/pubsed ~/personal/pubsed/src/sedona6.ex param.lua)
python compare.py

# ladder (Figure 5)
cd ../line_ladder && python make_atoms.py   # then the four run_* dirs as above
python compare.py

# La II forest + sweep (Figures 6-7)
cd ../laII_forest && python setup.py        # atom + model + line data
# run_bb / run_exp as above, then:
python compare.py
python sweep.py && python fig7.py
```

## 8. Next steps (in rough order of value)

1. **T sweep** (2500–5000 K): populations shift and new lines activate — the
   G-parameter axis of the map.
2. **Multi-ion overlap:** add Ce II/III from the same GSI database; the
   crowding statistics (Figure 2) say blending grows quickly.
3. **Thermal-width frontier:** v_D 10 → 0.6 km/s; measure both Δ_Sob's wing
   term and the cost curve.
4. **A true per-line Sobolev leg** in the harness, to separate
   "Sobolev approximation error" from "expansion-opacity implementation
   error" — the distinction F4 exposed.
5. **Frame-consistent comparison at realistic velocities** (0.1–0.3c),
   resolving F3 properly rather than avoiding it.
