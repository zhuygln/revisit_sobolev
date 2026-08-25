# Sobolev Validity in Kilonova Ejecta — Results Report

**Status:** Phase 0 through the validity maps and mechanism isolation, complete
— findings F1–F14, 13 figures, 68 tests. **Date:** 2026-08-18. **Repo:** `zhuygln/revisit_sobolev`.

**Manuscript:** [paper/manuscript.pdf](paper/manuscript.pdf) — the paper drawn
from this report, written for readers without a radiative-transfer background.
Companion: [lab_notebook.md](lab_notebook.md) (chronological log, including
dead ends and fixes). Planning inputs, preserved unmodified:
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
| `formal_transfer.py` | the deterministic solver (§2); Gaussian or Voigt profiles, optional truncation | blackbody-sphere luminosity exact; e^−τs trough; 2-line independence & blend; emission fill-in |
| `spectra.py` | **the one** band-ratio convention for SEDONA spectra: continuum-ratio normalization, margin clear of the final bin | null spectrum returns exactly 1; scale- and Planck-slope-invariant; known trough recovered |
| `sobolev_leg.py` | p-averaged per-line Sobolev attenuation; `damp` switches to the expansion cap | single line → exp(−τ_S) to 1e-12; damp variant → exp(−(1−e^−τ)); pop_frac scales τ; partial shadowing bracketed |

Test suite: **68 tests, all green** (`pytest`).

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

**Finding F3 (frozen-snapshot artifact — see §4.14).** A resonance at
z/ct = 7.7% produced a trough matching exp(−τ_S(1−z/ct)) to four digits
rather than exp(−τ_S). Read first as a frame ambiguity and then as the
leading relativistic correction; **both readings were wrong**. It is an
artifact of integrating a frozen snapshot of the ejecta. The physical law is
τ/τ_S = 1/γ, which has no first-order term. Superseded by F11.

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
| Python formal solver | 0.3583 (worldline default) |
| SEDONA resolved | 0.3426. Like-for-like — emission off (F8) *and* solver in frozen mode to match SEDONA's frozen transport (§4.14) — the gap is **−0.53%** |
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
| Python solver | 0.2456 (worldline default; 0.2426 under the frozen mode used originally) |
| SEDONA resolved | 0.2277 (7.8% from the worldline solver, 6.5% from the frozen one — see §4.14 on matching transport treatments) |
| SEDONA expansion | 0.2610 (**Δ_Sob = +14.6%**) |

**Finding F7 (non-monotonic density dependence).** The blend saturates the
whole window to near-black in *every* treatment — with enough crossings,
even capped per-crossing contributions sum to blackness — so the relative
band-flux error *shrinks* to +14.6%. Δ_Sob is non-monotonic in forest
density: ~τ²-small for weak forests, maximal (+45%) for strong sparse
forests, declining again toward full line blanketing. The dangerous regime
for expansion opacity is the *intermediate* one — strong lines that are not
yet fully blanketed — which is exactly where spectral features form.

Data note: Ce II carries half-integer J as fraction strings ('7/2'), now
handled by `sobolev.populations.parse_j`.

### 4.11 Resolution of the resolved-legs gap (Finding F8)

The solver-vs-SEDONA offset (+3.6% single-ion, +6.5% blend) was tracked to
its cause. Two hypotheses were tested and **disconfirmed**:

- *Profile wings.* Running the solver with Voigt profiles (γ = A/2π) and
  SEDONA's hard-coded ±5-width truncation changed the band flux by
  **3×10⁻⁵**. The damping parameter here is a ≈ 3×10⁻⁵, so Lorentzian
  wings never rise above the Gaussian where opacity matters. (Ce II's
  f(A) vs f(log gf) is self-consistent to 10⁻⁴, matching La II.)
- *Path resolution.* Band flux is flat at 0.2410–0.2411 across 4→64
  z-points per Doppler width, and 0.2414→0.2408 across 20→320 rays. Both
  quadratures are converged.

**Finding F8 (the cause): shell thermal emission.** The solver integrates
the LTE source term S = B_ν(T_shell) along every ray; SEDONA, run at fixed
temperature with radiative equilibrium disabled, deposits absorbed packet
energy without re-emitting it. Setting T_shell → 0 in the solver isolates
the term exactly:

| experiment | solver *with* emission | *without* | SEDONA resolved |
|---|---|---|---|
| La II (single-ion) | 0.3538 (+3.3%) | 0.3390 (**−1.0%**) | 0.3426 |
| La II + Ce II blend | 0.2410 (+6.5%) | 0.2249 (**−1.2%**) | 0.2277 |

Like-for-like, the two independent codes agree to **~1%** — better than the
3.6% previously quoted. (§4.14 sharpens this to **−0.53%** once the transport
treatment is also matched; an intermediate claim of 0.06% there was retracted.)
The term is
larger than the naive B_ν(3000)/B_ν(6000) = 2×10⁻³ estimate because the
emitting shell subtends 9× the core's projected area, and it grows with
saturation, which is why the blend showed a bigger offset than the sparse
single-ion forest.

Interpretation: the solver's treatment is the physically complete one
(Kirchhoff's law demands an LTE gas emit what it absorbs); SEDONA in this
configuration performs a pure-attenuation calculation. Absolute band fluxes
quoted from SEDONA here are therefore attenuation-only. **Δ_Sob is
unaffected**, being a differential between two SEDONA runs sharing the
identical convention.

### 4.12 Separating the two errors (Figure 10, Finding F9)

![Figure 10](figures/fig10_sobolev_vs_expansion.png)

Everything before this section measured the **expansion-opacity
implementation**. This section measures the **Sobolev approximation proper**
— per-line `exp(-τ_S)` with delta-function resonances — and reports both
against the same truth.

In this pure-absorption LTE setup the Sobolev prediction is exact analytics:
the p-averaged staircase of Appendix A.5, now promoted to
`sobolev/sobolev_leg.py` and shared with the ladder experiment (whose five
numbers are unchanged — the regression check). A Monte Carlo implementation
of per-line Sobolev interactions would add only noise. The same function
yields the expansion prediction via `damp = 1 − e^(−τ)`.

La II forest, band-averaged L/L_cont against SEDONA resolved = 0.3426:

| treatment | band flux | Δ |
|---|---|---|
| Sobolev proper | 0.3508 | **+2.1%** |
| expansion (analytic cap) | 0.4844 | +41.0% |
| expansion (SEDONA) | 0.4973 | +44.8% |

Across the τ_max × v_D sweep (Δ in %, vs SEDONA resolved):

| τ_max | v_D | Δ Sobolev | Δ expansion |
|---|---|---|---|
| 0.5 | 10 | +0.2 | +3.2 |
| 0.5 | 300 | +0.0 | +2.8 |
| 5 | 10 | −0.2 | +38.5 |
| 5 | 300 | +9.2 | +53.0 |
| 50 | 10 | +1.1 | +45.3 |
| 50 | 300 | +32.9 | +93.8 |

**Corrected.** An earlier version of this table reported a v_D-independent
Sobolev error of 5–8%. That was a normalization artifact — SEDONA band fluxes
normalized by raw luminosity rather than by the continuum ratio, leaving the
Planck slope in the answer, compared against a correctly-normalized analytic
leg. Δ_expansion barely moved, being a same-code differential. See §4.15.

**Finding F9: F6's two error components have different owners.**

1. The **strength-set floor is expansion-opacity's alone**. Sobolev proper
   shows *no* strength floor: at v_D = 10 km/s it is +5–7% at τ_max = 0.5,
   5, and 50 alike, while expansion climbs +3% → +37% → +45%. The floor is
   the per-crossing `1 − e^(−τ)` cap, not a failure of Sobolev's locality or
   isolation assumptions.
2. The **v_D-growing wing term is a genuine Sobolev failure**, shared by
   both. A delta-function resonance has no width by construction, so the
   Sobolev leg is *exactly* v_D-independent (0.3507 at every v_D); all of
   the resolved calculation's width dependence is un-modelled. At
   (τ_max = 50, 300 km/s) this alone reaches +32.9% — but F12 later showed
   this to be a finite-region boundary artifact, not a Sobolev failure.

So the Sobolev approximation itself is accurate to **≲2%** at
kilonova-relevant line widths (≤ 30 km/s) and to **≲0.5%** at τ_max = 5 down
to 1 km/s, while its standard expansion-opacity implementation is off by
**+38–48%** in the same regime.
The two are conflated in the literature under one name; they are not the
same approximation and they fail for different reasons.

Two honest caveats. (i) At τ_max = 0.5 SEDONA's expansion result (+3.0%) is
*closer* to truth than Sobolev proper (+7.1%) — but the analytic cap there
is +10.1%, so this is a compensating error, not accuracy: SEDONA's binning
smears absorption into inter-line gaps, darkening the spectrum by ~6%
relative to the pure per-crossing cap and accidentally offsetting the cap's
over-transmission. (ii) The residual +5–7% Sobolev floor at small v_D is
real and un-explained by strength; the natural candidate is line overlap
within a Doppler width (the isolation assumption), which the crowding
statistics of §4.3 say is present.

### 4.13 Breadth sweep: is the separation universal? (Figure 11, Finding F10)

![Figure 11](figures/fig11_breadth.png)

36 conditions — 4 windows (4300, 4900, 7000, 9100 Å, auto-selected by
optical-depth richness, *none* of them the reference window) × 3 epochs
(0.5, 1, 3 days at fixed ejecta mass, so ρ ∝ t⁻³ and τ_S ∝ t⁻²) × 3 ion
mixes (La II; +Ce II; +Ce III) — 72 SEDONA runs in 41 minutes. Line counts
range 71 → 2504, realized τ_max spans 0.003 → 34.

| | median | max | τ_max > 3 | τ_max < 0.5 |
|---|---|---|---|---|
| Δ expansion | +5.1% | **+61.3%** | +20.0% | +0.3% |
| Δ Sobolev | +3.6% | **+11.3%** | +6.6% | +0.5% |

Median ratio for τ_max > 1: **expansion errs 3.5× more than Sobolev**.

**Finding F10: the separation is universal, and realized τ_max is the
controlling variable among the axes swept.** (The sweep holds shell geometry
and v_D fixed; §4.12 shows Δ does depend on v_D, so the collapse is onto
τ_max *within a fixed (v_D, Δv_shell) slice* — not a claim that τ_max alone
determines the error everywhere.) All 36 conditions collapse onto a single trend in
τ_max regardless of window, epoch, or ion mix — wavelength and composition
enter only through the τ they produce. Both errors vanish together as
τ_max → 0 (the weak-line sanity check), and both grow with τ_max, but
expansion opacity reaches +61% where Sobolev proper never exceeds +12%.
F9's separation, established on one window, holds across the optical–NIR.

The τ_max-only scaling is a useful practical result in itself: a modeller
can estimate the error from the strongest line in a band without knowing
which ion or epoch produced it.

**Caveat on the earlier sign flip.** The first pass of this sweep reported
Δ_Sobolev ≈ −7% and appeared to contradict §4.12. That was a normalization
artifact, not physics — see the notebook entry 9f. The corrected values are
above; Δ_expansion was unaffected because it is a same-code differential.

### 4.14 Frozen snapshot vs photon worldline (Finding F11)

Prompted by external review of the Appendix A.6 derivation, which was wrong.

That derivation differentiated the Doppler factor `D = γ(1−b)` **at fixed t**
— pushing a photon through a frozen snapshot. A photon takes time to cross a
resonance, and in homologous flow β = r/(ct) depends on the age as well as
the position. Along the worldline `db/dt = (1−b)/t`, against the frozen
`db/dz = 1/(ct)`: a factor (1−b) difference, **first order in β**, the same
order as the effect. The time term does not vanish just because the crossing
is brief — what matters is the gradient.

**The two laws, confirmed to five decimals** by integrating the resolved
opacity along the path with no Sobolev assumption anywhere:

| β | frozen (t fixed) | worldline (t advances) | 1/γ |
|---|---|---|---|
| 0.05 | 0.94881 | 0.99875 | 0.99875 |
| 0.10 | 0.89549 | 0.99499 | 0.99499 |
| 0.20 | 0.78384 | 0.97980 | 0.97980 |
| 0.30 | 0.66776 | 0.95394 | 0.95394 |

So frozen gives (1−β)/γ and the physical law gives **1/γ = 1 − β²/2**, with
**no first-order term**. The correction at β = 0.1/0.2/0.3 is 0.5%/2.0%/4.6%,
not the 10%/22%/33% the frozen law implies.

**Which problem does SEDONA solve? — a retraction.** I first read the forest
numbers (solver vs SEDONA: −1.04% frozen-1st-order, −0.53% frozen-exact,
−0.06% worldline) as SEDONA independently confirming the worldline law, and
said so. **That was wrong.** At β ≤ 0.01 the two laws differ by only ~1–2%,
comparable to other systematics, and I over-read an ordering.

The v/c sweep has 23× more leverage and says the opposite. Comparing SEDONA's
measured τ_eff = −ln(F/F_cont) against both laws at the exact resonance
velocity b_res = (1−k)/(1+k), k = (1−β_label)²:

| β_res | SEDONA τ/τ_S | frozen (1−b)/γ | worldline |
|---|---|---|---|
| 0.105 | 0.923 | 0.890 | 1.111 |
| 0.220 | 0.788 | 0.762 | 1.250 |
| 0.342 | 0.613 | **0.618** | 1.429 |

RMS deviation for β_res ≳ 0.1: **0.024 from frozen, 0.552 from worldline.**

Cause confirmed in the SEDONA source: `transport_steady_iterate` sets
`use_hydro_ = 0` and reads no time-stepping parameters — it iterates the
radiation field on a fixed grid at fixed epoch. **It is a frozen-snapshot
calculation.**

**Consequence for the harness.** Cross-code comparison requires matching the
*transport treatment* as a third convention, alongside thermal emission (F8)
and spectral normalization (§4.13). Like-for-like against this SEDONA
configuration means running the solver in frozen mode. At the velocities used
throughout (β ≤ 0.01) the choice shifts τ by ≲2% and changes no reported
result; at kilonova velocities it would dominate.

**Finding F11.** Neglect of light-travel-time evolution is a distinct
approximation — separate from Sobolev localization, from the independent-line
assumption, and from expansion-opacity binning. A fourth member of the family
this project is mapping.

**Consequence, opposite to the hypothesis that prompted the work.** At the
shell velocities used throughout (β ≤ 0.01) the correct relativistic
correction is 5×10⁻⁵. It therefore explains **no part** of the +5–11% Sobolev
residual, which belongs to overlap or something else. The v/c hypothesis is
eliminated, and the controlled-overlap experiments become the decisive test.

**Open discrepancy.** The single-line minimal model moved the other way
(solver-vs-SEDONA −3.4% → −5.6%). Time-anchoring was tested and is not the
cause (ray-start vs centre-plane differ by 0.4%). The gap is most plausibly
SEDONA-side discretization of a single sharp line — recall §4.12 found its
expansion mode ran ~6% darker than the analytic cap through bin smearing —
which would average out across a 153-line forest. Unresolved; it does not
affect Δ values, which are same-code differentials.

### 4.15 The residual was an artifact; what remains is geometry (F12)

![Figure 12](figures/fig12_boundary.png)

Two reviewer arguments closed this out.

**Overlap cannot be the cause — an analytic proof.** With fixed populations
and pure absorption, opacities add, so τ_ν = Σᵢ∫α_ν,ᵢds = Σᵢτ_ν,ᵢ *exactly*
and I = I₀exp(−Στ_S,ᵢ) whether profiles overlap or not. Driving two identical
lines from 20 Doppler widths apart to exact coincidence reproduces the
Sobolev prediction to **six decimal places** at every separation. Overlap
becomes real only with scattering, fluorescence or NLTE — all off here by
design. **The isolation assumption cannot be probed by an attenuation-only
harness**, so correlating the residual against a crowding statistic would
have measured a confound (strength, density and edge-proximity are mutually
correlated).

**Most of the residual was a normalization artifact.** Computing the
Sobolev residual for the *existing* v_D = 1 and 3 km/s runs exposed that
`sweep.py` normalized SEDONA band fluxes by **raw luminosity** in the red
margin rather than by the **continuum ratio**, leaving the Planck slope in
the answer while the analytic leg was correctly normalized. This is the third
instance of the class (cf. F8, and the breadth red-edge bug). Corrected, at
τ_max = 5:

| v_D [km/s] | Δ Sobolev | Δ expansion |
|---|---|---|
| 1 | **−0.3%** | +39.1% |
| 3 | −0.3% | +39.2% |
| 10 | −0.2% | +38.5% |
| 30 | +0.5% | +41.3% |
| 100 | +2.5% | +45.0% |
| 300 | +9.2% | +53.0% |

Δ_expansion moved by ≤2 points everywhere — same-code differentials are
immune, which is the check that the fix changed only what it should.

**Finding F12: what remains is geometric.** Where a resonance falls within a
few Doppler widths of a boundary of the line-forming region, the resolved
profile is clipped while Sobolev applies a hard step:

  τ_resolved/τ_S = ½[erf((z_hi−z_res)/Δ) − erf((z_lo−z_res)/Δ)], Δ = v_D·t

— unity deep inside, exactly ½ at an edge, zero outside. Direct integration
matches to four decimals at every width, and the curves collapse when plotted
against d/v_D (Figure 12, left). The band-averaged effect is the fraction of
the velocity span within a few widths of an edge. With two edges that is
~2 v_D/Δv_shell: 1% at v_D = 10 km/s, 30% at 300 km/s. That is an upper
bound rather than a prediction — only the clipped part of each profile is
misassigned, and the measured |Δ_Sob| is 0.2% and 9.2%. What it captures is
the linear scaling, hence the two orders of magnitude the effect loses at
physical widths.

**This is an artifact of the setup, not of the approximation.** It scales with
the thermal width, ~0.6 km/s for lanthanides against ejecta spans of
10⁴ km/s — utterly negligible in reality. Our sensitivity to it at
v_D = 100–300 km/s comes from using artificially broad lines to keep frequency
grids affordable.

**Net:** at physical line widths the Sobolev approximation is accurate to
≲0.5%, while expansion opacity errs by ~39%. Essentially all of the error
commonly attributed to "Sobolev" belongs to the expansion-opacity
construction.

### 4.16 Is the expansion-opacity error a bin-width artifact? (F13)

![Figure 13](figures/fig13_binwidth.png)

The sharpest referee objection to Paper I: expansion opacity is derived for
bins containing *many* lines, and our grids are far finer — 12.5 km/s bins
against ~50 km/s mean La II spacing, about one line per four bins. Have we
measured the formalism outside its design regime?

Appendix A.4 says no analytically (the bin width cancels; each crossing
contributes 1−e^−τ at any resolution). Confirmed numerically over two decades
of transport resolution at fixed physics:

| line list | lines/bin | F_resolved spread | Δ_expansion |
|---|---|---|---|
| La II (153 lines) | 0.025–2.5 | 0.13% | **+43.7% ± 1.2** |
| + Ce II (2529 lines) | 4.1–41 | 0.5% | **+15.3% ± 0.7** |

The two lists give different errors (F7 — the blend is denser and partly
blanketed) but neither depends on bin width, **including at 41 lines per bin**,
squarely inside the regime the construction targets. The error is intrinsic.

**A methodological trap, recorded because it nearly reversed the conclusion.**
Extending the sweep to bins several times wider than the line profile makes
Δ_expansion appear to collapse and change sign (+8.6%, −32.7%). It has not:
past that point the *resolved* leg stops resolving, its band flux climbing
0.342 → 0.442 → 0.679, and the comparison measures the reference rather than
the approximation. My first script computed a verdict over all resolutions and
printed "NOT INVARIANT → reframe the claim". **The reference's own convergence
must be established before its disagreement with anything else means
anything** — and excluded points must be shown with their reason, not dropped
silently.

### 4.17 Paper II Phase 0 — a branching instrument, built and calibrated (F14)

Paper II's question is whether the expansion-opacity bias of F9/F13 survives
radiative redistribution. Answering it needs a code in which a photon absorbed
in one line can leave in another. Two audits decided which code that would be.

**P2-0A: the public SEDONA cannot do it at all.**
`paper2/phase0/sedona_source_audit/NOTES.md` records the search. Outside
`sandbox/` there is no fluorescence path, no macroatom, no downbranching. What
exists is `opacity_epsilon`, a single global scalar
(`opacity/GasState_opacities.cpp:438`), splitting extinction into an absorptive
part; an interaction is then either coherent scattering or a redraw from the
zone's total emissivity (`transport/scatter.cpp:12,303`). **Neither channel
knows which line absorbed the photon.** Adding real branching means carrying
upper-level identity through the interaction — information the expansion-opacity
path deliberately discards, which is the same loss that produces F4's
per-crossing cap. That is core packet-interaction work, not a configuration
switch. Verdict: absent, not deprecated.

**P2-0D reconnaissance: TARDIS has the physics but not the data.**
`paper2/phase0/tardis_install/NOTES.md`. TARDIS 2026.8.10.dev7 installs
(only via the repo lockfile — conda-forge has no package and the PyPI name is
a 2015 stub) and `line_interaction_type: scatter | downbranch | macroatom` is
a first-class configuration option. But every mode fails at atomic-data load,
before transport: TARDIS's bundled downloader 404s because the LFS objects of
`tardis-regression-data` are gone server-side, and the one file still
retrievable (`tardis-atomdata`, `database_version v0.9`) predates the pandas-HDF
format the current reader expects. The blocker is upstream data distribution,
so **nothing has yet been learned about the modes themselves.**

**The instrument was therefore built here** —
`paper2/phase0/three_level_atom/branching_mc.py`, ~250 lines, the same shape as
the formal solver that produced this project's most defensible results. It is
also the only option that preserves the standing rule: TARDIS has branching but
no expansion opacity, SEDONA has expansion opacity but no branching, so a
comparison between them would be cross-code. One code that can vary *both* is
what Phase 1 needs.

Physics: 1-D homologous shell, opaque lightbulb core, Sobolev point
interactions, frozen populations, first-order Doppler. The three-level atom is
ground 1, metastable sink 2, upper 3; the 1→3 line carries the only opacity and
level 3 branches to 1 (resonant) or 2 (fluorescent, redder). Because comoving
frequency decreases monotonically along any ray, a packet re-emitted at a line
centre can never reach that resonance again — enforced by a positive distance
tolerance, not by special-casing.

**F11 does not bite here.** The worldline correction is τ→τ/γ, i.e. O(β²); at
the 0.003–0.017 c of this geometry that is ≤3×10⁻⁴, below MC noise at any
practical packet count. F11 matters because it distinguishes 1/γ from a
spurious O(β) law; there is no O(β) term to get wrong.

**Calibration** (`run.py`, 2×10⁵ packets; `results.json`):

| check | predicted | measured | worst |
|---|---|---|---|
| fluorescent yield, A₃₂/(A₃₁+A₃₂) at 5 ratios | 0.00–0.75 | matches | 0.96σ |
| interaction probability, 1−e^(−τ_S) at τ = 0.1–10 | 0.095–1.000 | matches | 1.41σ |
| pure-absorption spectrum vs the analytic Sobolev leg | — | max abs. diff 0.0070 | 2.07σ |

The third is the one that matters. `sobolev.sobolev_leg.sobolev_attenuation`
was written for Paper I, is independently tested, and reaches the same emergent
spectrum by integrating over impact parameter rather than by sampling rays.
Agreement pins the geometry, the resonance-plane solve and τ_S together; the
first two checks alone would not catch a transport bug.

**One methodological catch, and it is the familiar one.** The first version of
that comparison reported a 27σ disagreement in exactly two bins. It was
comparing a *bin-averaged* MC against a *midpoint-evaluated* analytic. The
trough edges are near-vertical — the resonance plane leaves the shell over a
fraction of a percent in frequency — so at those two bins the midpoint value
and the bin average differ enormously while everything else agrees to <1%. This
is the same failure as the red-margin straddle that once manufactured a 7%
error in the breadth sweep (§4.13), and the same as the bin-width verdict of
§4.16. Averaging both sides over the bin removed it entirely. **Three times
now, a "disagreement" has turned out to be the two sides being asked different
questions.**

The instrument is calibrated; it has measured nothing new. Phase 1 is where
branching gets switched on against an expansion-opacity leg in the same code.

### 4.18 The referee revision: the mechanism made exact, and a deterministic reference (F15–F18)

A major-revision report (August 2026) held that the paper had labelled the
difference between a deterministic finite-profile attenuation and a
statistical closure as proof the closure was "wrong". It was right, and the
tools to test the reframing were already here.

**The mechanism (F15).** Per ray, Sobolev transfer exponentiates S = Σ τ_k and
expansion opacity exponentiates E = Σ(1−e^−τ_k). E is what α_exp integrates
to — the expected number of line interactions per crossing for a photon that
is counted but not removed, Karp's mean-free-path statistic — and the closure
preserves it identically at any bin width. Transmission is e^−S (a Bernoulli
product) against e^−E (Poisson with the same mean); they agree while every
line is weak and separate by the saturation deficit D = S−E per ray once
lines saturate. F_exp/F_Sob = ⟨e^D⟩_w with w ∝ e^−S, exactly (pinned to 1e-13
in `tests/test_rays.py`). On the La II forest:

| τ_max | F_Sob | F_exp (analytic) | Δ_exp | ⟨E⟩/⟨S⟩ |
|---|---|---|---|---|
| 0.1 | 0.950 | 0.952 | +0.2% | 0.97 |
| 1 | 0.665 | 0.720 | +8.3% | 0.74 |
| 5 | 0.350 | 0.484 | **+38.2%** | 0.37 |
| 25 | 0.253 | 0.374 | +47.5% | 0.12 |

`experiments/r1_interaction_count/ray_diagnostic.py` draws it along one ray
(46 resonances: same count to 1e-15, survivals 0.051 vs 0.130, first-interaction
masses 0.949 vs 0.870) and `mc_check.py` confirms every panel by Monte Carlo.
`experiments/r8_saturation/plot.py` shows all 54 conditions on the 1:1 line
against ln⟨e^D⟩_w; the unweighted deficit is a Jensen bound that fails once
rays saturate; τ_max is a proxy with real scatter (the blanketed blends below
the sparse forests — F7).

**A deterministic reference (F16).** Both legs on one `RaySet`
(`sobolev/rays.py`; midpoint, because attenuation is a step function of p and
matched nodes make the O(1/n) error common — unmatched rays had made Δ_Sob
change sign between 12 and 96 solver rays). The resolved leg for uniform n_l,
Gaussian profile, pure absorption is closed-form (`resolved_attenuation`, the
boundary erf bracket per line per ray, cost independent of v_D; agrees with the
brute-force solver to 8e-5). Results, La II forest:

| v_D (km/s) | Δ_Sob^det (τ=5) | Δ_Sob^det (τ=50) | Δ_exp^det (τ=5) | Δ_exp SEDONA pair |
|---|---|---|---|---|
| 300 | +9.6% | +32.7% | +51.3% | +53.4% |
| 100 | +3.1% | +8.9% | +42.3% | +44.9% |
| 30 | +0.9% | +2.5% | +39.3% | +40.4% |
| 10 | +0.30% | +0.82% | +38.4% | +38.9% |
| 3 | +0.09% | — | +38.1% | +38.6% |
| 1 | +0.03% | — | +38.1% | +39.1% |

Δ_Sob^det at v_D = 100 is +3.09% in all four transport modes (first, classical,
exact, worldline+dilution): the mode dependence cancels in the matched pair.
SEDONA resolved (seed mean) validates the reference to −0.5% on the headline,
±0.3% over the grid (max 0.7%), −0.04% median / 0.31% spread over the 36
breadth conditions. **The breadth Δ_Sob row was stale**: `breadth/recompute.py`
still normalized by raw luminosity (the bug 1e2ba21 fixed elsewhere); through
`band_ratio` it is median +0.34% / max +7.8% (SEDONA ref) and +0.26% / +7.6%
(deterministic ref), down from +3.65% / +11.3%, while Δ_exp is unchanged
(+5.09% median, +62% max). The line-free window returns 0.997 ± 0.004.

**Seed-matched pairs (F17).** `mc_noise/seeds.py`, 90 runs: headline over 10
matched seeds F_res 0.3422 ± 0.0003, F_exp 0.4958 ± 0.0004, paired Δ_exp
+44.90% ± 0.04 (sem), corr +0.95; every grid point 3 pairs, correlations
+0.92–0.999, paired scatter 0.02–0.28% vs 0.2–0.4% quadrature. The
single-seed scatter sits within 1.5× the Poisson expectation.

**Radiative equilibrium (F18).** `laII_forest/re_run.py`, 3 matched pairs per
variant, blue-margin normalization (re-emission contaminates the red margin):

| variant | F_res | F_exp | Δ_exp |
|---|---|---|---|
| pure absorption (same normalization) | 0.343 | 0.495 | +44.3% |
| RE, one iteration (re-emit at input T) | 0.836 | 0.900 | **+7.7% ± 0.02** |
| RE, T converged (median 3000 → ~7700 K) | 0.855 | 0.898 | **+5.0% ± 0.7** |
| τ_max = 0.05 control | 0.993 | 0.999 | +0.6% |

The removed flux reappears at 3955–3995 Å (1.08–1.10 of continuum). Most of
the closure's departure in attenuation does not survive re-emission into the
emergent band; its sign does.

**Single-line benchmark.** The published 0.1372 (solver) is the frozen-snapshot
law e^−τ_S(1−β)/γ over the trough window (0.1371) — the problem a steady-iterate
code solves — reproduced to 0.1% under every variation; 0.1353 is the β→0
limit, not the target. SEDONA's 0.1420 sat on a grid with <2 bins per Doppler
width and 2×10⁶ packets; 3.2×10⁷ packets on the same grid give 0.1386 ± 0.0002.
The convergence ladder (`minimal_1line/ladder.py`, 59 runs, fixed seeds)
gives 0.1383 ± 0.0001 at the anchor (3.2×10⁷ packets, 8 bins/width, 5 seeds),
flat to ±0.3% across packets 5×10⁵–3.2×10⁷, grids 4×10⁻⁴–2×10⁻⁵, spectrum
grids and 25–400 zones; the zero-opacity control reads 1.0079, so the +0.9%
residual is a continuum offset (comoving-frame core emission), and corrected
the trough is 0.1372 — the frozen target to 0.1%. The expansion mode sits at
0.4438 (0.440 corrected) vs Poisson 0.4212: +4.5%, the exp-leg bin systematic
on one line. Figure: `docs/figures/fig_ladder.png`; table generator
`minimal_1line/ladder_table.py`.

### 4.19 Paper II Phase 1 — the whole ion, five treatments, one code (F19–F21)

#### 4.19.1 Instrument

`paper2/phase1/forest_mc.py` — Phase 0's successor: packets advanced in
lockstep with numpy, the next resonance by one `searchsorted` on the sorted
line frequencies, the expansion closure as a continuous absorber whose
interaction point is found by inverting the F13 cumulative, and the
treatments in one code:

| leg | absorption | re-emission |
|---|---|---|
| sobolev_absorb / expansion_absorb | Sobolev point / Poisson closure | none (Paper I's controls) |
| sobolev_thermal | Sobolev | LTE line emissivity A n_u, photon-number weighted, optionally window-confined; Sobolev escape probability β on the emitting line |
| expansion_thermal | Poisson closure | the closure's own Kirchhoff emissivity κ_exp B_ν per bin (saturating at 1−e^−τ per strong line), uniform within the bin, no β (E0) |
| sobolev_branch | Sobolev | A-weighted downward channel of the upper level, re-drawn on re-absorption with probability 1−β (F19) |
| sobolev_tla / expansion_tla | as above | thermalisation parameter ε: thermal with probability ε, coherent scattering otherwise, every event (E4) |
| expansion_branch | Poisson closure | absorbing line sampled within the bin ∝ (1−e^−τ_k)/E_b, exit by the A·β kernel (E8) |

Atom: all of La II from the GSI files — 17,743 lines in the branching tables
and emissivity, **949 with τ_S > 10⁻³ as opacity, 1148–17,609 Å** at 3000 K
(19 of 472 levels populated); stimulated emission as SEDONA applies it.
Geometry, epoch, density: Paper I's headline forest (v = 1000–3000 km/s,
day 1, n_ion for τ_max = 5 in 3850–3950 Å). Packets are photons carrying
h ν; launch flat (calibration) or photon-number Planck at 6000 K (physics)
over 1142–17,697 Å; 3 seeds × 2×10⁶ packets per leg; ~10 s per leg.
Band quantity: escaped/launched, photon- or energy-weighted (both reported;
<0.5% apart in the 3800–3955 Å band).

#### 4.19.2 Validation

- **Paper I's legs** (Part A, window atom, flat launch): sobolev_absorb 0.3475
  ± 0.0005 vs analytic 0.351; expansion_absorb 0.4839 vs 0.485; Δ = +39.3% vs
  +38.1% (46 acceptance tests in `tests/test_forest_mc.py`, including both
  pure-absorption legs against the analytic Sobolev leg and the Poisson
  closure bin by bin).
- **SEDONA's radiative equilibrium** (Part A, window-confined thermal legs vs
  RE N=1, sub-bands bluewing / forest / band / red): resolved+thermal MC
  0.454 / 1.025 / 0.866 / 1.113 vs SEDONA 0.418 / 0.990 / 0.835 / 1.108;
  expansion+thermal MC 0.632 / 1.007 / 0.905 / 1.085 vs 0.632 / 0.998 /
  0.900 / 1.083. The expansion leg matches to <1%; the Sobolev leg sits 3.7%
  high at v_D = 100 km/s — the Sobolev-vs-resolved offset (Paper I's +3.1%)
  carried through re-emission.
- **E2 — thermal-width convergence** (`e2_thermal_width.py`, SEDONA RE N=1
  at 10 km/s, 3 seeds × 2 modes; prediction written before the runs: bb band
  0.835 → ~0.866, Δ_SEDONA +7.8% → ~+4.5%, expansion ~0.90):

  | | bluewing | forest | band | red | Δ_exp (band) |
  |---|---|---|---|---|---|
  | SEDONA resolved, v_D = 100 | 0.418 | 0.990 | 0.8349 | 1.108 | +7.83% |
  | SEDONA resolved, v_D = 10 | 0.459 | 1.002 | 0.8546 | 1.108 | +5.22% |
  | SEDONA resolved, v_D = 1 (1 seed) | 0.463 | 1.003 | 0.8563 | 1.109 | |
  | SEDONA expansion, 100 / 10 | 0.630 / 0.617 | 0.999 / 1.003 | 0.9002 / 0.8992 | 1.082 / 1.084 | |
  | MC Sobolev + thermal (v_D-free) | 0.454 | 1.023 | 0.8622 | 1.113 | +4.62% |
  | MC expansion + thermal, bins 1.25 / 12.5 / 125 km/s | 0.630 / 0.632 / 0.497 | 1.005 / 1.005 / 1.071 | 0.9013 / 0.9021 / 0.9044 | 1.085 / 1.086 / 1.103 | |

  The resolved side moved toward the MC as predicted (band offset +3.3% →
  +0.9% → +0.7% at 100 / 10 / 1 km/s, sub-bands ≤ 2.1%; 10 → 1 km/s moves
  the band by 0.2%, so the width dependence is converged at 10 km/s), the expansion side did not move, and the MC row
  is constant because its Sobolev legs are delta resonances and its closure
  depends on bin width, not v_D — the bin-width leg shows the closure's
  re-emission placement is stable at 1.25–12.5 km/s bins and degrades only at
  125 km/s (blue wing 0.63 → 0.50). Gate A passed. The remaining +0.7% is
  the Sobolev-vs-resolved offset under re-emission at the level Paper I's
  F16 measured in attenuation (+3.1% at 100 km/s → sub-percent below 10).
- **The escape probability (F19).** Without β on re-emission the MC's thermal
  flux came out too blue and too spread (blue wing 0.71 vs SEDONA 0.42); with
  it, the agreement above. The trapped fluorescence yield b/(1−(1−b)(1−β)),
  the geometric re-absorption count (1−β)/β, and the A·β exit kernel are all
  tests.
- **E0 — the closure's own emissivity.** The expansion legs had re-emitted from
  the Sobolev line emissivity A n_u and then applied β: switching β off made
  the SEDONA agreement worse, bin-uniform placement alone worse still; the
  difference was the emissivity. κ_exp B_ν per bin (saturating), bin-uniform,
  no β, reproduces SEDONA's expansion RE to <1% in every sub-band. The
  "β on" agreement had been coincidental.
- **E1 — energy.** E_inj = E_esc + E_core + E_abs + E_dep closes to roundoff
  in every mode; the comoving exchange equals the level-energy difference per
  branching chain (10⁻¹² on the toy atom, 2×10⁻⁵ on La II); the O(v/c)
  Doppler work term is reported separately. Thermal legs at fixed T are
  bookkept, not conserved — why SEDONA iterates T.

#### 4.19.3 Phase-1 result

Full atom, Planck 6000 K launch, band 3800–3955 Å, 3 seeds × 2×10⁶:

| leg | band | re-emission goes |
|---|---|---|
| Sobolev, absorb | 0.183 ± 0.002 | — |
| expansion, absorb | 0.344 (+88% vs Sobolev) | — |
| Sobolev + thermal | 0.257 | 78% redward, 9% beyond 1 μm |
| expansion + thermal (κ_exp B_ν) | 0.408 | 61% redward |
| **Sobolev + fluorescence** | **0.660 ± 0.003** | 54% redward, 11% blueward, 35% in-band |

Paper I's atom (window, flat launch): fluorescence refills the band by only
+12.8% (0.348 → 0.392); 51% of the band's absorbed photons leave redward
through lines the window never had; window-confined thermal re-emission
refills it to 0.866 — which is why Paper I's RE check is window-bound.

#### 4.19.4 Current interpretation

- ε = 1 (complete thermalisation, re-emitting from the closure's own κ_exp
  B_ν) fails: −38.2% ± 0.3 below direct branching in the optical band, after
  +88% above it in pure absorption. The closure's error changes sign once
  fluorescence matters.
- This is **provisional**: it establishes only that ε = 1 fails, not that no
  scalar ε can work. The ε sweep (§4.20) decides.
- The pumps are 3300–4500 Å, not the far-UV (E7): 0.00% of the band's escaped
  energy was launched in 1142–2500 Å; 971 pathways, top 10 = 25%, dominated by
  5d6p upper levels exiting the strong ground-connected lines after 2–8
  re-absorptions.
- Limits: frozen LTE populations at 3000 K, one ion, first-order Doppler,
  photon packets, an opaque core absorbing ~12% of inward re-emission (as
  SEDONA's does), a Planck continuum for a photosphere, one emission per
  absorption event (the exit lower level's excitation returns to the pool).

#### 4.19.5 Literature connection

Kasen et al. (2006) and the SN Ia work around SEDONA compared direct
iron-peak fluorescence with a scalar two-level-atom thermalisation parameter
and found complete thermalisation (ε = 1) redistributes too much redward while
an intermediate ε ≈ 0.3 approximates the fluorescence reasonably; Fontes et
al. (2020) and Morag (2026) bound what the expansion-opacity closure can and
cannot represent; TARDIS (macroatom), ARTIS (line-by-line fluorescence) and
SUMO (NLTE) are the codes that carry line identity. The lanthanide question is
new: open-4f-shell networks are far richer than iron-peak ones, and whether
one scalar ε transfers is exactly what E4 tests.

### 4.20 Paper II Phase 2 — does any scalar ε reproduce La II fluorescence? (F22)

`paper2/phase1/e4_eps_sweep.py` (E4/E5), `e6_redistribution.py` (E6),
`e4_fig.py`; full La II atom, Planck 6000 K photon launch over 1142–17,697 Å,
3 seeds × 2×10⁶ packets per leg, energy-weighted band fractions F_b =
escaped/launched energy. Four legs: Sobolev + direct A·β branching (the
physics), expansion + branching (E8), Sobolev + TLA(ε), expansion + TLA(ε),
ε ∈ {0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1}. On the Sobolev leg the opacity is
identical to the physics leg, so only the redistribution closure differs; on
the expansion leg both the opacity representation and the closure differ.
Band edges checked against strong lines (3304 Å τ = 1.2 sits 4 Å inside the
blue band; the 3800–3955 Å edges are the Paper I ones). Figure
`docs/figures/fig_p2_eps_sweep.png`; χ² in `fig_p2_eps_chi2.png`.

#### 4.20.1 The sweep (E4)

| band [Å] | Sobolev+branch | expansion+branch | Sob+TLA ε=0 → 1 | exp+TLA ε=0 → 1 |
|---|---|---|---|---|
| UV 1142–3300 | 0.9373 ± 0.0025 | 0.9493 (+1.3%) | 0.972 → 0.901 | 0.979 → 0.920 |
| blue 3300–4500 | 0.8284 ± 0.0013 | 0.8629 (+4.2%) | 0.902 → 0.532 | 0.928 → 0.628 |
| optical 4500–6000 | 0.9705 ± 0.0005 | 0.9796 (+0.9%) | 0.952 → 0.838 | 0.964 → 0.869 |
| red 6000–9000 | 0.9936 ± 0.0002 | 0.9954 (+0.2%) | 0.988 → 1.119 | 0.990 → 1.092 |
| NIR 9000–17,697 | 1.0055 ± 0.0001 | 1.0049 (−0.1%) | 1.000 → 1.028 | 1.000 → 1.025 |
| 3800–3955 | 0.6590 ± 0.0032 | 0.7995 (+21.3%) | 0.824 → 0.257 | 0.884 → 0.410 |

F_b(ε) is monotonic in every band (falling below 4500 Å, rising above
6000 Å: thermalisation moves energy toward the 3000 K emissivity peak), so
ε_best per band is a single interpolated root, with the 3-seed scatter
propagated through the interpolation:

| band | Sobolev+TLA ε_best | expansion+TLA ε_best |
|---|---|---|
| UV | 0.36 [0.32, 0.41] | 0.68 [0.63, 0.73] |
| blue | 0.055 [0.054, 0.057] | 0.24 [0.23, 0.25] |
| optical | **not reachable** (target 0.9705 > F(0) = 0.952) | **not reachable** (> 0.964) |
| red | 0.016 [0.014, 0.018] | 0.026 [0.021, 0.031] |
| NIR (null control) | 0.13 | 0.22 |
| 3800–3955 | 0.065 [0.063, 0.068] | 0.32 [0.31, 0.33] |

**Verdict (Gate B): outcome B.** No scalar ε reproduces direct La II
fluorescence. With the opacity held fixed (Sobolev leg) the best ε is
near-zero in the blue and red (0.02–0.06), 0.36 in the UV, and the optical
4500–6000 Å band lies outside the TLA's whole range: direct branching
delivers more optical energy than even pure coherent scattering (ε = 0),
because the cascade feeds the optical from the blue (E6 below), which no
local two-level redistribution can do. The expansion leg's ε_best ≈ 0.3 in
the 3800–3955 Å band — the iron-peak-looking number — is a compensation
between an opacity that is 21% too transparent under fluorescence (E8) and a
closure that thermalises too much; it does not carry to the other bands
(0.03 red, 0.68 UV).

#### 4.20.2 Spectral χ² (E5)

200 log bins over the launch range, σ² = σ_TLA² + σ_branch², against
Sobolev+branch. χ²/dof: Sobolev+TLA 44 (ε=0), 47, 110, 178, 303, 406, 505,
545 (ε=1); expansion+TLA 80, 60, 53 (ε=0.2), 58, 94, 152, 229, 273. Minima
at ε = 0 and ε ≈ 0.2, neither a fit; the residual sits in the blue (3300–4500
Å, χ²/dof 15–23 at the minimum) and the optical (15–21), not the UV or NIR.
The ranking, not the absolute value, is the result.

#### 4.20.3 Redistribution matrices (E6)

P(λ_out | λ_in), 60 log bins, rows normalised by launched energy
(`e6_redistribution.npz`, `fig_p2_redistribution.png`). Block sums for
launches in 3300–4500 Å (the pump band):

| leg | stays blue | → optical 4500–6000 | → red 6000–9000 | → NIR | escaped |
|---|---|---|---|---|---|
| Sobolev + branch | 0.786 | **0.088** | 0.005 | 0.000 | 0.886 |
| expansion + thermal (ε=1) | 0.604 | 0.070 | **0.110** | 0.016 | 0.800 |
| Sobolev + thermal (ε=1) | 0.501 | 0.085 | 0.141 | 0.018 | 0.746 |
| expansion + TLA ε=0.3 | 0.795 | 0.039 | 0.042 | 0.006 | 0.882 |
| Sobolev + TLA ε=0.3 | 0.639 | 0.076 | 0.085 | 0.009 | 0.809 |

Branching is a blue → optical channel (8.8% of blue launched energy, 0.5% to
the red); the thermal closure is a blue → red channel (11–14% to the red);
intermediate ε interpolates between the two and so never reproduces the
first without the second. The optical rows of the physics are fed by the
blue and lose 3% back to the blue — the diagonal-plus-cascade structure in
the left panel of the figure against the diagonal-plus-Planck-smear of the
closure. Mean own-bin energy share: branch 0.867, expansion+thermal 0.864,
Sobolev+thermal 0.815, expansion+TLA(0.3) 0.925.

#### 4.20.4 The branching-aware closure (E8)

`expansion_branch` — Poisson absorption in the bin, absorbing line sampled
∝ (1−e^−τ_k)/E_b, exit by the exact A·β kernel — lands +21.3% ± 0.7 above
Sobolev+branch in 3800–3955 Å and within +0.2–4.2% in the wide bands, at
dν/ν = 4.17×10⁻⁵ (12.5 km/s bins). That is the opacity-representation error
under fluorescence, to be read against +88% in pure absorption and against
the −38% (ε=1) / +34% (ε=0) redistribution-closure error: carrying line
identity through the bin recovers most of the physics with an O(lines)
table (line identity per bin + per-level A·β CDF). It is bin-width
dependent and does not converge to the physics as the bins shrink
(`e8_binwidth.json`, 3 seeds): 3800–3955 Å 0.821 ± 0.003 / 0.800 ± 0.005 /
0.769 ± 0.004 at 1.25 / 12.5 / 125 km/s bins, i.e. +25% / +21% / +17% above
Sobolev+branch, the wide bands within 0.1–1% of each other except the blue
(0.860 / 0.863 / 0.892). What remains at the finest bins is the Poisson-vs-
Bernoulli opacity error itself, which Paper I showed is bin-independent in
the saturated regime; the self-absorption overlap (a re-emitted photon
sweeping the rest of its own bin) is the caveat that goes with it.

### 4.21 Paper II Phase 2.75 — the closure verdict at 0.1c (E13, F23)

`paper2/phase1/e13_worldline.py`; `forest_mc.run_mc(relativity="worldline")`
carries each packet's own clock — exact Doppler D = γ(1−β_z), the linear
resonance locus z_res = Z0(y²−1)/2 + p²/2Z0 from Paper I's addendum,
homologously expanding boundaries (the core overtakes launches with
μ < β_core; the caught fraction equals β_core² exactly), every optical depth
diluted to its resonance's own epoch τ_S(t_exp)(t_exp/t_res)²/γ in the
interaction and escape probabilities alike, and comoving-isotropic
re-emission aberrated to the lab. `relativity=None` is the classical
first-order transport on a time-frozen shell — the E13 frozen control, i.e.
what a naive code does at high β. The three velocity scales stay separate:
v_D never enters (delta resonances, §4.20's E2), Δv_shell sets the forest
sweep, v_bulk the transport convention.

#### 4.21.1 Single-line control (E13.3)

MC vs `sobolev_attenuation`, both conventions, τ_S = 5, band-mean
transmission (rms per frequency bin ≤ 0.003 everywhere):

| β_out | worldline MC / analytic | frozen-first-order MC / analytic | convention gap mean / max |
|---|---|---|---|
| 0.01 | 0.7919 / 0.7919 | 0.7907 / 0.7907 | −0.001 / 0.085 |
| 0.05 | 0.6158 / 0.6158 | 0.6048 / 0.6045 | −0.011 / 0.417 |
| 0.10 | 0.5580 / 0.5581 | 0.5318 / 0.5319 | −0.026 / 0.783 |
| 0.20 | 0.4963 / 0.4959 | 0.4320 / 0.4320 | −0.064 / 0.993 |

The frozen snapshot manufactures an O(β) artifact where the physical
worldline correction is ~β²/2 — the trap the control was built to catch. Two
instrument bugs were caught before any science ran (the moving-core root
selection, and a roundoff re-selection of the just-used resonance that let
16% of packets skip the remaining forest under worldline transport; the
full-forest worldline−classical differential now matches the analytic leg to
0.1% absolute; `tests/test_forest_mc.py`, 4 new tests).

#### 4.21.2 Matched line strength (E13.2)

τ_S = σ f n λ t has no shell-velocity dependence, so keeping n_ion, T, t_exp
from the slow reference makes the fast shell's τ distribution identical by
construction (τ_max 8.4, N(τ>1) = 42, N(τ>0.1) = 206). What changes is the
sweep: Δv_shell = 0.1c crosses ~15× more resonances per unit ln ν.

#### 4.21.3 Full La II, slow (0.0033–0.01c) vs fast (0.05–0.15c) shells

Worldline transport on both; energy-weighted; 3 seeds × 2×10⁶. The slow
shell reproduces §4.20's classical values shifted by the +1% dilution term
(branch band 0.6528 vs 0.6590). Selected legs:

| leg | blue (slow → fast) | optical (slow → fast) | NIR (slow → fast) | 3800–3955 (slow → fast) |
|---|---|---|---|---|
| Sobolev + branch | 0.828 → 0.043 | 0.971 → 0.355 | 1.005 → 1.085 | 0.653 → 0.044 |
| Sobolev + TLA ε=0 | 0.901 → 0.088 | 0.952 → 0.436 | 1.000 → 0.995 | 0.817 → 0.089 |
| Sobolev + TLA ε=1 | 0.533 → 0.030 | 0.838 → 0.219 | 1.028 → 1.151 | 0.258 → 0.027 |
| expansion + TLA ε=0.3 | 0.809 → 0.046 | 0.938 → 0.341 | 1.007 → 1.057 | 0.675 → 0.040 |
| expansion + TLA ε=1 | 0.629 → 0.040 | 0.871 → 0.256 | 1.025 → 1.155 | 0.411 → 0.037 |
| **Sobolev + branch, frozen-first-order** | **0.163** | **0.396** | 1.064 | **0.168** |

#### 4.21.4 The three verdicts (E13.6)

- **Transport convention dominates the observable.** On the fast shell the
  frozen control transmits 3.8× more blue-band energy than worldline
  transport (0.163 vs 0.043) — larger than any branch-vs-TLA difference.
  Frozen high-β results must not be used to assess the closure; every
  comparison below is worldline-vs-worldline.
- **The qualitative verdict survives: outcome B holds at 0.1c.** On the fast
  shell ε_best still differs irreconcilably by band — expansion leg: 0.20
  (3800–3955), 0.24 (optical), 0.46 (blue), 0.47 (UV), 0.49 (NIR), red
  unreachable (every TLA overfills it; branching does not). No scalar ε
  reproduces the branching spectrum at realistic bulk velocity either.
- **But calibrated ε values do not transfer.** |ε_best(0.1c) − ε_best(slow)|
  reaches 0.15–0.27 on the expansion leg (grid systematic ~0.05), and two
  bands change reachability status (optical becomes reachable, red becomes
  unreachable). Per the decision rule (|Δε_best| ≳ 0.2), v_bulk is a
  validity-map axis, not a nuisance parameter.

#### 4.21.5 Redistribution at 0.1c (E13.5)

Row-block sums of P(λ_out | λ_in) for Sobolev+branch (launches in
3300–4500 Å): slow shell — stays blue 0.785, → optical 0.088, → red 0.005,
escaped 0.886; fast shell — stays blue 0.015, → optical 0.121, → red 0.106,
escaped 0.273. The 0.1c sweep pushes every interacting packet through the
whole forest: the blue → optical fluorescence channel survives (and grows),
but it now competes with a comparable blue → red channel and 73% in-shell
deposition. The redistribution structure changes substantially — the same
conclusion as the ε shifts, seen mechanistically
(`e13_matrix_{slow,fast}.npz`).

### 4.22 Paper II Phase 4 — Ce II and the mixture: outcome C, and the closure's density limit (E9–E10, F24)

`paper2/phase1/e9_ceII.py`, `e10_blend.py`; `ForestAtom.from_gsi_blend`
(level-offset concatenation, branching kept ion-internal — ions share only
the radiation field). Ce II normalized by the identical recipe as the La II
reference (n_ion for a strongest classical τ of 5 in 3850–3950 Å at
3000 K → n_ion = 11,641 cm⁻³, 2,376 window lines): **22,960 opacity lines**
(La II: 949), τ_max 6.7, N(τ>1) = 149, opacity reaching the far-IR. The
blend carries each ion at its own reference density (23,909 opacity lines,
N(τ>1) = 191). Slow shell, classical transport, Planck 6000 K photon
launch, energy-weighted, 3 seeds × 2×10⁶.

#### 4.22.1 Ce II key legs (3800–3955 Å band, energy-weighted)

| leg | UV | blue | optical | 3800–3955 |
|---|---|---|---|---|
| Sobolev, absorb | 0.349 | 0.005 | 0.350 | **0.0001** |
| expansion, absorb | 0.380 | 0.112 | 0.412 | **0.224** |
| Sobolev + branch (physics) | 0.500 | 0.495 | 0.785 | 0.556 |
| expansion + branch (E8) | 0.496 | 0.592 | 0.810 | **1.183 (+113%)** |
| Sobolev + TLA ε=0 → 1 | 0.815 → 0.353 | 0.588 → 0.069 | 0.886 → 0.446 | 0.531 → 0.072 |
| expansion + TLA ε=0 → 1 | 0.810 → 0.382 | 0.662 → 0.161 | 0.899 → 0.496 | 0.746 → 0.316 |

#### 4.22.2 Outcome C: ε_best is ion-dependent

| band | La II (Sob leg) | Ce II | blend | La II (exp leg) | Ce II | blend |
|---|---|---|---|---|---|---|
| UV | 0.36 | 0.48 | 0.48 | 0.68 | 0.43 | 0.45 |
| blue | 0.06 | 0.04 | 0.04 | 0.24 | 0.09 | 0.09 |
| optical | unreach | 0.11 | 0.09 | unreach | 0.16 | 0.16 |
| red | 0.02 | 0.12 | 0.13 | 0.03 | 0.22 | 0.21 |
| NIR | 0.13 | 0.11 | 0.10 | 0.22 | 0.13 | 0.12 |
| 3800–3955 | 0.07 | **unreach** | **unreach** | 0.32 | 0.17 | 0.15 |

ε_best shifts by 0.11–0.24 per band between the ions, and reachability
flips in opposite directions: La II's optical band is unreachable while
Ce II's is reachable; La II's 3800–3955 Å band is reachable while Ce II's is
not (branching delivers 0.556, pure scattering on the same opacity tops out
at 0.531 — the cascade feeds the band from outside the TLA's reach). Ce II's
χ²/dof minima are 108 (expansion, ε≈0.1) and 167 (Sobolev) — worse than
La II's 44–53. The blend tracks Ce (the dominant forest) in every band:
ε_best is set by composition, so a universal scalar ε does not exist even
per band (outcome C on top of outcome B).

#### 4.22.3 The branching-aware closure hits a density limit

On La II, `expansion_branch` (Poisson absorption + exact A·β exit kernel)
came within +21% of the physics. On Ce II it **overfills the band by +113%**
(1.183 vs 0.556; blend +107%). The cause is in the absorption controls: with
2,376 lines in the 100 Å window the true (Bernoulli) transmission is black
(10⁻⁴) while the Poisson closure's per-bin saturation clipping transmits
0.224 — an opacity error of three orders of magnitude that no exit-kernel
correction can repair, and re-emission on top of an opacity that black
should occur deep enough that little escapes; the closure instead re-emits
where its too-transparent opacity puts the interactions. The line-identity
recommendation (§4.20.4) is therefore **density-limited**: it holds where
the Poisson and Bernoulli counts stay close (La II's 949-line forest), and
fails where saturation clipping dominates (Ce II's 22,960). Paper I's F15
mechanism, now measured under fluorescence.

#### 4.22.4 The mixture's answers to E10's three questions

1. Blanketing does not suppress the redistribution error: expansion+TLA
   ε=1 sits −44.6% below the blend's branching (La-only: −38%).
2. It does not move it either: the blend's band pattern is Ce II's.
3. ε_best does change with forest density/composition — by tracking the
   dominant ion (table above), which is the outcome-C statement itself.

## 5. Findings register

| # | Finding | Where |
|---|---|---|
| F1 | Symmetric gradients cancel the Sobolev error at leading order; G-based diagnostics are conservative | §4.2 |
| F2 | Resolved vs expansion cost: 378× on the Type Ia toy | §4.4 |
| F3 | ~~O(v_bulk/c) frame systematic~~ → an artifact of frozen-snapshot integration; superseded by F11 | §4.5, §4.14 |
| F4 | Expansion opacity attenuates by exp(−(1−e^−τ)) per crossing, not e^−τ | §4.5 |
| F5 | That error is per-resonance: it does not average away with line count, only as τ→0 | §4.6 |
| F6 | Δ_Sob = strength-set floor + v_D-growing wing term; error survives the v_D→0 limit | §4.8 |
| F7 | Δ_Sob is non-monotonic in forest density: maximal for strong sparse forests, suppressed at full blanketing | §4.10 |
| F8 | The resolved-legs offset is shell thermal emission, not profile wings or resolution; like-for-like the codes agree to ~1% | §4.11 |
| F9 | The strength floor belongs to expansion opacity alone; the v_D wing term is the finite-region boundary effect of F12, not a Sobolev failure. At v_D ≤ 30 km/s: Sobolev ≲2%, expansion +38–48% | §4.12 |
| F10 | The separation is universal across windows, epochs and ion mixes; realized τ_max controls it **at fixed v_D and geometry** (§4.12 shows the v_D dependence), and expansion errs ~3.5× more than Sobolev | §4.13 |
| F11 | Neglect of light-travel-time evolution is a distinct approximation: frozen τ/τ_S = (1−β)/γ vs worldline 1/γ. The physical law has no O(β) term | §4.14 |
| F12 | Overlap is inert in pure absorption (optical depths add exactly); the Sobolev residual is a finite-region boundary effect ∝ v_D/Δv_shell, negligible at thermal widths | §4.15 |
| F13 | The expansion-opacity error is bin-width invariant from 0.025 to 41 lines per bin — intrinsic to the formalism, not a usage artifact | §4.16 |
| F14 | Neither candidate reference code can answer Paper II alone: SEDONA has no line branching whatsoever, TARDIS has no expansion opacity. A ~250-line branching Sobolev MC, validated to ≤2.1σ against the analytic Sobolev leg, supplies the same-code differential both lack | §4.17 |
| F15 | The closure preserves the expected interaction count exactly and applies Poisson survival to a Bernoulli product; F_exp/F_Sob = ⟨e^D⟩_w exactly | §4.18 |
| F16 | Δ_Sob against a deterministic reference on identical rays: +3.1% at 100 km/s (all transport modes), +0.3% at 10, +0.03% at 1 km/s; the breadth median was stale (+3.65% → +0.26%) | §4.18 |
| F17 | Seed-matched pairs correlate at +0.95–0.999; paired Δ_exp scatter 0.02–0.28%; headline +44.90% ± 0.04 over 10 seeds | §4.18 |
| F18 | The radiative-equilibrium check is window-confined (re-emission drawn from the emissivity on SEDONA's transport grid); its +7.7%/+5.0% are bookkeeping within one 230 Å window, not an emergent-band result — retracted as such | §4.18 |
| F19 | Escape-probability branching: radiative branches compete as A·β, β = (1−e^−τ)/τ; the chain's exit distribution is A_uj β_uj/Σ A β (tested); Phase 0 lacked it, the SEDONA RE comparison caught it | §4.19 |
| F20 | Fluorescent optical refill: full La II under a 6000 K continuum, 0.183 → 0.660 ± 0.003; pumped from 3300–4500 Å (the far-UV contributes nothing), 971 pathways, top 10 = 25% | §4.19 |
| F21 | Closure sign reversal: +88% in pure absorption, −38% with complete thermalisation (ε = 1) and +34% with pure scattering (ε = 0), relative to direct branching; +21% once line identity is carried through the bin (`expansion_branch`; density-limited — F24) | §4.19, §4.20 |
| F22 | No scalar ε reproduces La II fluorescence (outcome B): ε_best 0.02–0.06 (red, blue), 0.36 (UV) on the same opacity, the optical 4500–6000 Å band unreachable at any ε; branching is a blue → optical channel, thermalisation a blue → red one | §4.20 |
| F23 | At v_bulk ~ 0.1c with worldline-consistent transport, outcome B survives (ε_best 0.20–0.49 by band, red unreachable) but calibrated ε values shift by 0.15–0.27 and the redistribution structure changes; the frozen-snapshot convention overstates blue transmission 3.8× and must not be used at high β | §4.21 |
| F24 | Outcome C: ε_best is ion-dependent (La vs Ce shifts 0.11–0.24, reachability flips in opposite directions; the La+Ce blend tracks Ce) — and the branching-aware Poisson closure is density-limited: +21% on La II's 949-line forest but +113% on Ce II's 22,960-line forest, where saturation clipping leaves the closure 3 orders of magnitude too transparent in the band | §4.22 |

## 6. Caveats and limitations

- **Both treatments are now measured, but not in a scattering code.** F9
  measures the Sobolev approximation proper alongside the expansion-opacity
  implementation. The Sobolev leg is *analytically exact* here **only because
  the setup is pure absorption with no scattering** — that is what makes the
  p-averaged staircase equal to per-line Sobolev transport. A per-line
  Sobolev scheme inside a Monte Carlo code with scattering or fluorescence
  (TARDIS-style macroatoms) is a genuinely different object and remains
  untested; F9's headline should not be read as covering it.
- **Pure absorption only:** ε = 1, no scattering, no fluorescence, no NLTE.
  Deferred by design (babystep_plan.md §9) — differences measured here are
  attributable to line transfer alone.
- **Controlled geometry:** shell velocities 1000–3000 km/s keep the O(v/c)
  frame systematic (F3) below the measurement level. Realistic kilonova
  velocities (0.1–0.3c) require frame-consistent treatment before the
  comparison is repeated there.
- **v_D = 1–300 km/s:** the frontier reached 1 km/s, within a factor ~2 of
  the La thermal width (0.6 km/s), where F6's floor is flat. The last factor
  of two is untested and costs ~1/v_D in transport bins.
- **Coverage:** the detailed maps use one window (3850–3950 Å) at day 1; the
  breadth sweep (§4.13) extends to four windows spanning 4300–9100 Å, three
  epochs and three ion mixtures, and confirms the separation there. What is
  still absent is a full lanthanide mixture — two elements and three
  ionization stages is not the dozen-plus species of real r-process ejecta,
  and F7 says density matters.
- **MC noise:** ~1–2% per SEDONA band flux at 2×10⁶ core particles. The
  resolved-leg agreement is ~1% once the F8 emission convention is matched
  (§4.11) — i.e. at the noise floor.
- **Emission convention:** absolute SEDONA band fluxes quoted here are
  attenuation-only (F8). Δ values are same-code differentials and unaffected.

## 7. Reproduction

```bash
# environment
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]" h5py
pytest                    # 68 passed

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

# La II forest + sweeps (Figures 6-8)
cd ../laII_forest && python setup.py        # atom + model + line data
# run_bb / run_exp as above, then:
python compare.py
python sweep.py && python fig7.py           # tau x v_D slice
python tsweep.py && python fig8.py          # T axis + thermal-width frontier

# multi-ion blend (Figure 9)
cd ../multiion && python setup.py
# run_bb / run_exp as above; the solver leg is slow, precompute it:
python solve_py.py && python compare.py

# Sobolev-proper vs expansion separation (Figure 10) -- analytic, no runs
cd ../sobolev_proper && python compare.py

# breadth: windows x epochs x ions (Figure 11)
cd ../breadth && python sweep.py            # ~36 conditions, 72 SEDONA runs
python recompute.py && python fig11.py      # band fluxes + figure

# the paper
cd ../../docs/paper && make                 # pdflatex x2 + bibtex -> manuscript.pdf
```

Long jobs: launch in the background with `python -u` and an **absolute** path
to the venv interpreter. A session restart resets the shell's working
directory, which silently breaks relative `../../.venv/bin/python`
invocations (exit 127) and kills in-flight runs.

## 8. Next steps (in rough order of value)

Completed since the first draft of this report: the T sweep and thermal-width
frontier (§4.9), multi-ion overlap (§4.10), the per-line Sobolev leg (§4.12)
— which resolved the approximation-vs-implementation distinction that F4
exposed — and the breadth sweep (§4.13), which established that the
separation is universal.

1. ~~Explain the residual Sobolev error.~~ **Done (§4.15).** It was mostly a
   normalization artifact; what remains is a finite-region boundary effect,
   not a failure of either Sobolev assumption. Overlap was excluded
   analytically and numerically.
2. **Frame-consistent comparison at realistic velocities** (0.1–0.3c).
   *Tooling done, measurement in pilot — see lab notebook §9m–9n.* Three
   things changed the shape of this item:

   - The controlling term is **not** the 0.5–4.6% frame systematic. A photon
     meets material moving at β at t_res = t₀/(1−β), so with n ∝ t⁻³ the
     medium contributes (1−β)² against transport's 1/γ: **51% against 4.6% at
     β = 0.3**, an order of magnitude the other way. Confirmed against direct
     integration of the resolved opacity.
   - The **CD/CP geometry objection is retired** for the physical treatment.
     Under worldline transport the resonance locus is linear in z and
     τ/τ_S(t_res) = 1/γ exactly for every impact parameter; the two-root
     quadratic belongs to the frozen snapshot.
   - **SEDONA can supply the time-dependent side** — its steady *mode* cannot,
     but the code's time-dependent mode dilutes ρ as t⁻³, expands the grid and
     carries photons across timesteps. All 165 param files here force the
     frozen mode and silently discard their own `hydro_module`.

   Pilot result (synthetic forest, converged in both grids, *not* a numbered
   finding): Δ_Sobolev is β-independent to 0.05 points, so Sobolev proper does
   not degrade at realistic velocity, while Δ_expansion **grows** +17.1% →
   +21.7% from β = 0.01 to 0.3. Only ~2/3 of that growth is the τ_max shift —
   the collapse hypothesis is refuted at ratio 1.35 — and the residual
   mechanism is unidentified. Needs the La II forest and a τ_max range before
   it can be quoted.
3. **Scattering and fluorescence.** Beyond their intrinsic importance, this
   is where the analytic Sobolev leg stops being exact, so it is the regime
   in which a per-line Sobolev *transport* scheme must be built and tested
   rather than computed in closed form.
4. **NLTE populations**, the deferred Stage D of babystep_plan.md §16.
