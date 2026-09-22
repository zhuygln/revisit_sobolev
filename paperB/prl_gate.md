# Paper B — gate G1, preregistered

Written before any production run (2026-09-21; revised 2026-09-22 on the
PI's seven pre-run changes, `plan_review.md`, before any result existed).
The state, the grids, the live-band rule, the metrics and the thresholds
below may not change after the run. What may: the packet count upward if a Gray condition is noise
(the whole grid is rerun, never one leg), and bug fixes with the run
repeated. G2–G4 are named at the end and are detailed here, each as its
own section, only after G1 is read.

## The question G1 answers

Does the compression of Paper III (F25/F27: a small group-to-group
operator reproduces explicit fluorescence on the same Sobolev opacity)
survive indivisible energy packets and the energy-conserving downward
macroatom reference?

## State

- `paper4/phase1_benchmarks/P1_t2.json`, shell 28 (`EjectaState.local_zone`):
  t = 2 d, T_gas = T_rad = t_core = 3401.0 K, ρ = 3.298×10⁻¹⁵ g cm⁻³,
  r_core = 5.223×10¹⁴ cm, r_out = 5.952×10¹⁴ cm, core law `local_shell`.
- One ion per experiment, alone, at its P1 number density
  `state.n_ion(el, "II", 28, ATOMIC_MASS[el])`: `57LaII` 3.293×10⁴,
  `58CeII` 9.538×10⁴, `60NdII` 1.369×10⁵ cm⁻³ (singly ionised, f_ion = 1,
  as in Paper IV's II-only states).
- `ForestAtom.from_cached([(ion, n_ion)], 3401.0, 172800.0, tau_min=1e-3)`,
  the GSI calibrated line list, stimulated emission on.
- Transport `relativity="worldline"`, `launch_weight="energy"`,
  `chain_max = 2000`, `chain_overflow="absorb"`, window 1000–30000 Å,
  200 spectral bins, the Paper IV legs runner (`paper4/phase2_energy/legs.py::run_legs`).
- Packets: 3×10⁵ per seed; every leg `packets="energy"`. **Construction and
  evaluation are split:** the kernels are built from the reference's events
  on the *build seeds* (101, 102, 103), pooled into one kernel per ion and
  N_g; the reference photometry, the independent fine matrix and every
  closure's transport use the *evaluation seeds* (1, 2, 3). No kernel is
  scored against the events it was built from.

## Legs per ion

| tag | mode | what |
|---|---|---|
| `R2build` | `sobolev_dmacro` | the kernel-construction run on the build seeds; its events (absorbed frequency, exit frequency, energy weights) are pooled into every R_N kernel; not scored |
| `R2` | `sobolev_dmacro` | **the reference: the energy-conserving downward macroatom** (Lucy's downward-only macroatom with indivisible energy packets; not a radiation-field-coupled or NLTE macroatom), on the evaluation seeds; `collect_events=True` for the independent fine matrix |
| `A2_ng2 … A2_ng32` | `sobolev_group` | R_ij at N_g ∈ {2, 4, 8, 16, 32}: kernel from R2build's events, `RedistributionKernel.from_branching_mc(nu_in, nu_out, w_in, N_g, nu_lo=op_nu.min×0.995, nu_hi=op_nu.max×1.005, w_out=w_out)`, log-spaced groups, discrete within-group exit tables, `rows="energy"`; rows never populated scatter coherently |
| `K128` | kernel only | the 128-group energy matrix from **R2's** events (the evaluation seeds), not transported: the independent fine reference of the event-level metric |
| `K128build` | kernel only | the same from R2build's events: the in-sample fine matrix, recorded only as a diagnostic (its distance from K128 is the sampling noise of the fine matrix itself) |
| `E0.00 … E1.00` | `sobolev_tla` | ε ∈ {0, 0.05, 0.10, …, 1.00} (21 values): with probability ε a thermal re-emission from the energy-weighted emissivity, otherwise resonant scattering in the same line |

After a pass only (not in the gate): `R2M`, the radiation-field-driven
macroatom with upward transitions under an imposed W = ½ B_ν(T_zone)
(`sobolev_macro`), and an `A2_ng8` rebuilt from R2M's events, one run per
ion, as Paper IV's F63 did. The hierarchy is: G1 reference = downward
macroatom → post-pass robustness = radiation-field-driven macroatom. A
claim of "full macroatomic physics" would need stronger validation than
either and is not made.

## Metrics (per leg, seed-averaged, with the seed standard deviation)

- **m_event** (the information the binning discards): with the N_g-group
  energy matrix block-expanded onto the 128-group edges, R^(N), and the
  128-group energy matrix from the **independent** evaluation events,
  R^(128),

      d_i = ½ Σ_j |R^(N)_ij − R^(128)_ij|,     m_event = Σ_i W_i d_i / Σ_i W_i,

  with W_i = Σ_{e ∈ i} w_in,e ν_in,e the energy absorbed in fine row i: the
  expected total-variation loss of the outgoing redistribution distribution
  for an energy-weighted incoming interaction, in [0, 1]. The in-sample
  fine matrix's own distance from the independent one (K128build vs K128)
  is recorded as the sampling floor of this metric.
- **m_sed**: the integrated distance
  Σ_b |L_ν,b − L_ν,b^{R2}| Δν_b / Σ_b L_ν,b^{R2} Δν_b over the 200
  log-spaced bins of the window (L_ν is a spectral density; the bin widths
  are not equal).
- **m_band**: max |Δm| and mean |Δm| over the live bands, Δm = m − m^{R2}.
- **m_colour**: max |Δcolour| over the pairs of `sobolev.photometry.COLORS`
  whose two bands are both live.
- **m_energy**: the identity residual of every leg (|E_esc + E_core + E_abs
  + E_dep,lab − E_inj| / E_inj) and `RedistributionKernel.validate_energy()`
  of every kernel.
- **m_complexity**, stated honestly: N_g² is the **coarse transition-matrix
  degrees of freedom**, not the model's complexity — the kernel that works
  also carries the discrete within-group exit tables that fixed F25's
  self-absorption artefact. Recorded per kernel: N_g², the number of stored
  exit samples, the serialized size in kB (`RedistributionKernel.save`,
  matrix plus tables), events per packet and wall time. The central figure
  uses N_g² on its axis, labelled as such, with the serialized size as a
  companion panel; no claim of "millions of transitions compressed into
  N_g² numbers" is made from this gate.
- **Live bands** (G1's rule): finite magnitude, ≥ 1 % of the window
  luminosity, and an acceptable Monte Carlo precision (R2's seed scatter
  ≤ 0.05 mag in the band). The 40 Mpc detectability cut of Paper IV is an
  observational selection and is reported separately, not applied.
- **ε\***: the grid value minimising mean |Δm| over the live bands; reported
  with its max |Δm|, m_sed, m_colour, the full curve, and whether the
  minimum is interior to the grid.

## Readings (Gray checked first)

**B1 — the kill gate.** For each ion, N_g* is the smallest N in
{2, 4, 8, 16, 32} with max |Δm| ≤ 0.10 mag over the live bands **and**
max |Δcolour| ≤ 0.10 mag.

- Green: N_g* ≤ 8 for at least two ions and N_g* ≤ 32 for the third.
- Yellow: N_g* ≤ 32 for at least two ions, without the Green condition.
- Red: N_g* undefined (no N ≤ 32 reaches the threshold) for two or more
  ions, or only one ion reaches it at any N ≤ 32.

The reference of every reading is the energy-conserving downward
macroatom (`R2`); the kernels are built from `R2build`.

**B2 — ε\* is not nearly as good.** e_ε = ε*'s max |Δm| over the live
bands; e_R = the N_g*-matrix's max |Δm| (R_32's if N_g* is undefined).

- Green: e_ε ≥ 0.20 mag and e_ε ≥ 3 e_R.
- Red: e_ε ≤ 1.5 e_R.
- Yellow: otherwise.

**B3 — convergence, not strict monotonicity** (revised before the run:
transport is nonlinear enough that a broadband error need not fall
monotonically even when the operator converges). m_event is the structural
diagnostic; m_sed is strongly expected to converge; m_band may wiggle.

- Green: overall convergence in all three metrics, and no statistically
  resolved degradation from N_g = 8 to 32 in any of them (m(32) ≤ m(8) +
  noise, the noise being the seed scatter for m_band and zero for the
  others).
- Yellow: a statistically resolved local non-monotonicity in an observable
  metric (m_band or m_sed), or m_event not converging from 8 to 32, without
  the Red condition.
- Red: a representation pathology — m_event **degrades** from 8 to 32, or
  both observables (m_band by more than the noise and m_sed) are worse at 32
  than at 8, or the empty-row fallback exceeds the gray limit (which fires
  first). An innocent 0.01 mag wiggle does not kill the PRL path.

**Decision.** Continue toward the PRL iff B1 Green and B2 Green and B3 not
Red. Stop the PRL attempt (a methods result at most) iff B1 Red or B2 Red.
Anything else is Yellow: the PI decides with the numbers in hand, and
nothing in this file is edited to make it Green.

**Gray conditions** (per ion, checked before any reading):

1. fewer than two live bands on R2;
2. a band carrying ≥ 1 % of the luminosity was dropped from the live set
   by the precision rule (seed scatter > 0.05 mag; the 0.10 mag threshold
   would sit within twice the noise): the remedy is more packets on the
   whole grid, not a smaller live set;
3. any leg's identity residual > 10⁻¹⁰, or any kernel's
   `validate_energy()` > 10⁻¹²;
4. any leg failed, or was chain-capped (`n_trapped`) on > 1 % of its
   packets, or drew an empty kernel row (coherent fallback) on > 1 % of its
   interactions;
5. the ε* minimum at a grid edge (ε = 0 or 1) with the neighbouring value
   within the seed noise of it.

A Gray ion is reported with the condition that fired and excluded from B1's
and B2's counts; the gate then needs the Green conditions on the remaining
ions, and is itself Gray if fewer than two ions remain.

## Expected noise (stated, not tuned)

Paper IV's R2 seed scatter on the 13-ion P1 blend at 3 × 3×10⁵ was
≤ 0.02 mag in the live bands. The 0.10 mag threshold is five times that;
Gray condition 2 fires if a luminous band's single-ion scatter is 2.5
times worse.

## The central figure

Transport error against the coarse transition-matrix size: m_band (max),
m_sed and m_event against N_g² for the three ions, ε*'s m_band as a
horizontal line per ion, the seed noise as a band, and a companion panel
of the serialized kernel size (matrix plus exit tables) against N_g,
`docs/figures/paperB/gate1_error_vs_complexity.{pdf,png}`.

## Records

`paperB/gate1/gate1_<ion>.json` (the `run_legs` row plus the kernels'
validation, table size, empty-row usage, and the grids as run),
`paperB/gate1/gate1_verdict.json` (every reading with the values it was
read from and every gray check that ran). `analyse.py` refuses a record
whose state, seeds, build seeds, packet count, N grid or ε grid differ
from this file.

## G2–G4 (named; each detailed here only after the previous gate is read)

- **G2, locality vs rank**: the 128-group energy matrix replaced by its
  rank-k non-negative factorisation with 2·128·k ≈ N_g*² parameters,
  transported through `sobolev_group` and scored with G1's metrics; F32's
  photon-packet result redone with energy packets.
- **G3, universality**: (a) state transfer — P1 at 1, 3, 5 d (shells 29,
  27, 26) and a T_gas sweep at 2 d, the 2 d kernel fixed vs recomputed;
  (b) species composability — `RedistributionKernel.mix` with
  `paper3/phase5_mixture/mixture.py::composition_weights` predicting a
  La+Ce+Nd macroatom run at the P1 densities with no blend fit;
  (c) realistic mixture — the 13-ion P1 blend at N_g = 2 … 32 (A2 at 32
  is within 0.05 mag, F61).
- **G4, σ_R vs σ_comp**: on the composition patterns that exist in this
  repository (Ye-0.21a at X_lan 0.11 and 0.30, solar_r at 0.11, Ye-0.29a),
  σ_R = the spread of R_{N_g*} − macroatom and **σ_comp** = the spread of
  the macroatom across those patterns. This is a composition proxy, not
  the nuclear-model ensemble; the name σ_nuclear is reserved for the
  future experiment that runs an external nuclear-realisation ensemble
  through the same transport.
