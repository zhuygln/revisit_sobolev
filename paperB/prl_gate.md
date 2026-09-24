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

---

## R2M robustness (post-pass; fixed 2026-09-22 before the run; the G1 sections above are unchanged)

The check the post-pass paragraph of "Legs per ion" announces, run now on
the PI's instruction. **Not a gate**: it asks whether the compression G1
found is tied to the downward macroatom, and reports.

- State, atoms, transport, window, live-band rule, thresholds: G1's.
- Packets: 3×10⁵ per seed for all three ions. Nd II's G1 remedy (10⁶) is
  not repeated here because the 10⁶ downward run peaked at 17.4 GB of
  resident memory on the 24 GB machine (every collected leg's events stay
  in memory for the run) and the radiation-field-driven macroatom carries
  about 3.5 times the events per packet (F63). Remedy if Gray condition 2
  fires on R2M: rerun the whole R2M leg set at 10⁶ once the runner releases
  a source leg's events after its kernels are built (a memory change to
  `run_legs`, no change to any number), never a smaller live set.
- Legs per ion (build seeds 101–103; evaluation seeds 1–3; every leg
  energy packets):

  | tag | mode | what |
  |---|---|---|
  | `R2build`, `R2` | `sobolev_dmacro` | G1's reference pair, rerun on the same seeds so that R2M − R2 is measured seed for seed |
  | `A2_ng8` | `sobolev_group` | G1's downward-trained 8-group operator (kernel from `R2build`), scored here against R2M: the reference dependence of the operator |
  | `R2Mbuild` | `sobolev_macro`, W = ½, T_rad = T_zone | the radiation-field-driven macroatom on the build seeds, `collect_events` |
  | `R2M` | `sobolev_macro`, W = ½ | **the reference of this section**, evaluation seeds, `collect_events` |
  | `A2M_ng8` | `sobolev_group` | the 8-group operator rebuilt from `R2Mbuild`'s events, transported on the evaluation seeds |
  | `K128M`, `K128Mbuild` | kernel only | the 128-group energy matrices from `R2M` and `R2Mbuild` (the independent fine reference and the in-sample one) |

- Readings per ion, Gray first (G1's conditions 1–4 applied with R2M as
  the reference leg):
  1. the reference shift: R2M − R2 per band live on R2M, and its maximum
     (F63 gave 1.3–1.8 mag on the 13-ion blend);
  2. **survives** iff `A2M_ng8` vs R2M has max |Δm| ≤ 0.10 mag over R2M's
     live bands and max |Δcolour| ≤ 0.10 mag (G1's thresholds, unchanged);
  3. the reference dependence: `A2_ng8` (downward-trained) vs R2M, the same
     two numbers, reported;
  4. m_event of the R2M-trained 8-group kernel against `K128M`, of the
     downward-trained one against `K128M`, and the `K128Mbuild` vs `K128M`
     floor.
- The statement made from it: "the compression survives (or does not
  survive) the radiation-field-driven macroatom at N_g = 8 for ion X",
  with the reference shift alongside. No claim about full macroatomic
  physics.
- Records `paperB/r2m/r2m_<ion>.json`, `paperB/r2m/r2m_verdict.json`;
  `analyse_r2m.py` refuses a record whose state, seeds, packet count, N_g
  or W differ from this section.

---

## G2 — locality vs rank, and the observable-relevant dimension (preregistered 2026-09-22; not run until the PI approves this section)

**The two hypotheses** (the PI's, `plan_review.md` 2026-09-22):

- **H1.** The success of R_ij is caused by local smoothness in frequency,
  not by a globally low-rank microscopic redistribution law.
- **H2.** Transport observables occupy a much smaller effective information
  space than the underlying fluorescence-event distribution.

**Design principle.** Every G2 operator is a 128-group energy matrix on G1's
edges (`op_nu.min×0.995 … max×1.005`, log-spaced) carrying **the same
128-group discrete exit tables**, built from `R2build`'s events on the
build seeds; every operator is transported through `sobolev_group`
(`rows="energy"`) on the evaluation seeds and scored against G1's `R2` with
G1's live-band rule, metrics and thresholds. The families differ only in
what they do to the group-to-group matrix (or, for H2, to the tables), so
the comparison isolates that. The matched axis of H1 is the **number of
archetypal exit distributions**: a local operator at N_g says "the exit
distribution is one of N_g, chosen by which contiguous frequency block the
absorbed frequency falls in"; a global rank-k operator says "it is a free
non-negative mixture of k, with the weights free per fine input group". At
equal count the global family is the more expressive; the parameter counts
(N_g² against 2·128·k) are recorded alongside, never used as the axis.

**State, packets, seeds.** G1's; La II and Ce II at 3×10⁵, Nd II at 10⁶
(each ion's final G1 record); evaluation seeds 1–3, build seeds 101–103.
`R2build` and `R2` are rerun in the G2 run (deterministic: `R2`'s
magnitudes must equal G1's record's to 10⁻⁶ mag, else Gray) so the fine
matrices and the reference are the same objects as G1's. The runner
releases a source leg's events after its kernels are built (the memory
change of 2026-09-22).

**Legs per ion** (`paperB/gate2/run_gate2.py`, families in `operators.py`):

| tag | operator | archetypes | matrix parameters |
|---|---|---|---|
| `A2_ng8` | G1's 8-group kernel (coarse tables), rerun | 8 | 64 |
| `L128_ng8` | **local control**: the 8-group energy matrix expanded onto the 128 groups, each coarse column's mass spread by the fine exit-energy marginal within it, on the 128-group tables. Sampling-equivalent to `A2_ng8` by construction (both draw the exit line with probability ∝ its exit energy within the coarse output group); transported to verify that the local family of G1 (`A2_ng2 … A2_ng32`, read from the G1 records) is the local family on shared tables | 8 | 64 |
| `G_k1 … G_k32` | **global family**: rank-k non-negative factorisation W H of the populated rows of the 128-group energy matrix (Lee–Seung multiplicative updates, initialisation seed 0, iterated to convergence: the relative Frobenius error is checked every 100 iterations and the loop stops when its relative change over that block reaches 10⁻⁶, capped at 20,000 iterations), rows rescaled to their original sums (energy exact); a live row the factorisation sends to zero is declared **empty**, so transport applies the kernel's standing convention for a row carrying no information — coherent scattering, counted in the fallback fraction and limited by gray condition 4 — and the number of such rows and the energy share they carry are recorded; k ∈ {1, 2, 4, 8, 16, 32}. k = 1 is the fully non-local null: one exit distribution for every input | k | 2·128·k |
| `T_f0.1 … T_f0.999` | **truncation family**: the full 128×128 matrix with each output group's exit table cut, by the selection rule fixed under H2 below (ranked by accumulated training-event energy within the group, retained until their share first reaches f, at least one per populated group, both weight tables renormalised); f ∈ {0.1, 0.2, 0.5, 0.9, 0.99, 0.999}. The grid was set after the exit-table audit (§4.63) so that H2's Green is reachable: per-group truncation retains, of the distinct exit lines, Ce II 2.3 / 4.3 / 12 / 36 / 67 / 87 % and Nd II 0.4 / 0.8 / 2.4 / 11 / 39 / 77 % at these f | 128 | 128² |
| `K128`, `K128build` | kernel only: the independent fine matrix (evaluation seeds) and the in-sample one | | |

The local family's curve is G1's `A2_ng{N}` legs, N ∈ {2, 4, 8, 16, 32},
against G1's `R2` — the same reference, seed for seed.

**Metrics** per operator: G1's m_band (max and mean |Δm| over the live
bands), m_colour, m_sed, m_event against `K128`; plus the archetype count,
the matrix parameter count, the in-sample TV distance of the derived
matrix from the fine one it was derived from, the number of stored exit
lines, the serialized size; for the NMF the relative Frobenius
reconstruction error and its relative change over the last 100 iterations.

**Readings** (Gray first). Ce II and Nd II decide; La II is reported and
never counted (0.7 interactions per packet at its P1 density: every
operator is close).

Gray, per ion: G1's conditions 1–4 on every leg; (6) `R2` differs from
G1's record by more than 10⁻⁶ mag in a live band; (7) the control
`L128_ng8` differs from `A2_ng8` by more than 2 σ of the difference of two
legs (σ_diff = √2 × R2's seed scatter) in any live band — then the
equivalence argument fails and the local family must be rerun on the
128-group tables before H1 is read; (8) an NMF rank whose relative
reconstruction change over the last 100 iterations exceeds 10⁻⁴ is
excluded from the k* search and reported; (9) no N_g* in G1's record.

**H1 — locality, not global low rank.** Per ion, two numbers on the same
transport criterion (G1's: max |Δm| ≤ 0.10 mag over the live bands and
max |Δcolour| ≤ 0.10 mag against `R2`):

    K*_local  = min { N_g in {2,4,8,16,32} : the local operator passes }   (G1's N_g*)
    K*_global = min { k  in {1,2,4,8,16,32} : the rank-k operator passes }

both undefined if no member of the family passes. The comparison is of
**archetype counts, not parameter counts**: at equal count the global
family holds many more free numbers (2·128·k against N_g²), so needing no
fewer archetypes than the local family is evidence for locality.

- Green: K*_global ≥ K*_local for **both** Ce II and Nd II (undefined
  counts as ≥: a global family that never passes is the strongest form of
  this).
- Red: the global family passes with K*_global ≤ K*_local/4 on **either**
  Ce II or Nd II.
- Yellow: everything between.
- Overall: Red if either decisive ion is Red; Green if both are Green;
  Yellow otherwise; Gray if either is Gray.

At G1's numbers only Ce II (K*_local = 16) can fire Red, at K*_global ≤ 4;
Nd II's K*_local = 2 puts Red out of its reach, and its informative
outcomes are Green (K*_global ≥ 2) and Yellow (K*_global = 1, the fully
non-local null passing where two local blocks are needed).

**The event-level fit is a mechanistic diagnostic, not part of the pass/
fail logic** (the PI's amendment of 2026-09-22). Reported per operator
alongside the readings: m_event against the independent 128-group matrix
and the in-sample TV distance from the matrix it was derived from. The
outcome of interest — *the global model fits the microscopic events
better while the local model reproduces the transport better* — is
recorded as a named diagnostic (`events_vs_observables`, true when at the
matched count k = K*_local the global operator's m_event is the smaller
and its max |Δm| the larger). It supports the reading that
microscopic-distribution fidelity is not observable-relevant fidelity; it
is a stronger and more specific claim than the locality hypothesis and is
therefore never required for Green.

**H2 — the observable-relevant dimension of the exit spectrum.** Per ion,

    L* = min { retained fraction of the distinct exit lines : the truncated
               operator passes the same transport criterion }

where the retained fraction is (exit lines kept) / (exit lines in
`K128build`), evaluated on the f grid; undefined if no f < 1 passes.

- Green: L* ≤ 0.10 for **both** Ce II and Nd II.
- Yellow: Green fails but at least one of them has L* ≤ 0.25.
- Red: either Ce II or Nd II requires L* > 0.50 (undefined counts as > 0.50).
- Residual band (both ions between 0.25 and 0.50, neither above it):
  Yellow, reported with the values. The ladder is read Red first, then
  Green, then Yellow.
- La II is reported as a control and **does not decide H2**: G1 showed it
  too thin at its P1 partial density to discriminate closures.

**The selection rule is fixed here and is a function of the
kernel-building sample alone.** Within each output group of the
`R2build` 128-group kernel, the distinct exit lines are ranked by their
accumulated training-event energy weight in that group's exit table; the
highest-ranked are retained until their share of the group's exit energy
first reaches f (at least one line per populated group); both the photon
and the energy weight tables are renormalised over the retained lines.
Nothing about transport enters the choice, and no subset is reselected
after seeing which preserves the light: the only freedom is f, and its
grid is frozen above.

**Decision.** The PRL's mechanism claim ("locality, not rank") is written
iff H1 is Green; H1 Red reframes the claim (the PI decides what, with the
numbers); Yellow is the PI's call. H2 sets the wording of "compact": with
H1 and H2 both Green, both layers of the effective operator are compact —
the group-to-group matrix and the within-group exit spectrum — and the
"the matrix is small but the table is not" objection is answered.

**Amended 2026-09-22** on the PI's review of PR #10, before any G2 transport existed: the event-level fit left H1's pass/fail logic and became a diagnostic; H1's comparison stated as archetype counts K*_local and K*_global; H2 read on the retained fraction L* with a three-level ladder and Ce II / Nd II decisive; the selection rule written out as a function of the build sample alone. The f grid had been extended to 0.1 and 0.2 earlier the same day, after the audit.

**What may not change after the run:** the k and f grids, the control
rule, the thresholds (G1's), the archetype-count axis, the NMF settings
(seed 0, the convergence rule above, rescaled rows). What may: the packet count upward
on a Gray (the whole ion rerun), and bug fixes with the run repeated.

### Amendment of 2026-09-22 (second): two implementation defects, the run repeated

Declared here in full because it was made **after part of a G2 run existed**.

*What was defective.* (1) At k = 1 and k = 2 the factorisation sent one
live row to exactly zero; the row's energy sum no longer matched its
`q_dep`, so `validate_energy()` returned 1.0 and gray condition 3 fired
on those two legs. (2) The preregistered 600 iterations did not compute
the preregistered object: at k = 16 and k = 32 the relative Frobenius
error was still falling by ~3×10⁻³ per 100 iterations, so gray condition
8 fired on those two legs. Neither is a property of the physics; both are
failures to compute the rank-k factorisation, and the preregistration's
standing remedy for a bug is to fix it and repeat the run.

*What had been seen when the defects were found.* The **La II** record
was complete and was inspected: its local family (max |Δm| 0.036 / 0.021
/ 0.013 / 0.009 / 0.010 at N_g = 2 … 32), its global family (0.97 / 0.07
/ 0.185 / 0.033 / 0.063 / 0.055 at k = 1 … 32, two legs gray on the
energy identity and two on convergence), its truncation family, and the
derived K*_local = 2, K*_global = 2, L* = 0.50. **La II decides neither
H1 nor H2** — it is the preregistered control, too thin at its P1 partial
density to discriminate closures (G1, F67). The **Ce II** run was in
progress and its partial record was deleted unread; **Nd II** had not
started. No decisive ion's numbers were seen before this amendment, and
no threshold, grid or reading was changed by it.

*The changes.* The NMF is iterated to convergence (relative change ≤ 10⁻⁶
over a 100-iteration block, cap 20,000, seed 0 unchanged) instead of a
fixed 600; gray condition 8 stays as the backstop at 10⁻⁴. A live row the
factorisation zeroes is declared empty and scatters coherently, as above.

*The direction of the second change, stated so it can be judged.* Zeroing
such a row makes the low-k global operators slightly worse, which pushes
K*_global up and so leans **toward** H1 Green. The alternative — falling
back to the row's own un-factorised distribution — leans toward Red and
would leak information the rank-k model does not contain. On La II the
affected row carried 1.5×10⁻⁶ of the absorbed energy and the coherent
fallback reached 0.07 % of interactions, far inside gray condition 4's
1 % limit; both numbers are recorded per leg for every ion.

*The run.* All three ions were rerun from scratch after the fix; no
record produced before it survives.

### Amendment of 2026-09-22 (third): the control fired, its declared remedy applied

The `L128_ng8` control of gray condition 7 fired on **both** decisive ions
in the first complete run: 2.6 σ on Ce II and 2.1 σ on Nd II. The
preregistration's remedy for exactly this case is written above — "then
the equivalence argument fails and the local family must be rerun on the
128-group tables before H1 is read" — and that is what was done. **No
threshold was touched.**

What the control was actually resolving, recorded because it matters for
reading the diagnostic: the per-band differences between the coarse
8-group kernel and its 128-table expansion are 0.002–0.035 mag, random in
sign, and every one of them is far inside the 0.10 mag science threshold.
The σ values that produced 2.6 and 2.1 come from the tightest bands
(Ce II H: seed scatter 0.0027 mag; Nd II J: 0.0034 mag), so the test was
resolving differences some thirty times below the threshold the gate
cares about. The empty-row explanation was checked and excluded: the fine
rows the 128-group kernel leaves empty carry exactly zero absorbed energy
and the coherent fallback is 0.00000 of interactions in both legs. The
test is a max over six or seven bands compared against a per-band 2 σ
limit, which fires roughly a quarter of the time under the null; that is
a defect of the control's design, not evidence of a broken operator, and
it is left standing rather than loosened after the fact.

*What changes.* The local family is now the `L128_ng{2,4,8,16,32}` legs
transported in the G2 run itself, on the same 128-group exit tables and
the same row occupancy as the global and truncation families; K*_local is
read from them. G1's own `A2_ng{N}` legs stay in the record as
`local_coarse`, reported for comparison. The control comparison stays in
the record as a diagnostic (`control`: the σ, the largest per-band
difference in magnitudes, and the per-band values) and no longer gates the
reading, because the assumption it was guarding — that G1's coarse legs
may stand in for the local family — is no longer made. Every threshold,
grid, metric and ladder is unchanged, and all three ions were rerun.

**Records** `paperB/gate2/gate2_<ion>.json`, `gate2_verdict.json`,
`docs/figures/paperB/gate2_locality_vs_rank.{pdf,png}` (max |Δm| against
the archetype count for both families, m_event alongside, the truncation
curve against the fraction of exit lines kept, and every operator on the
(m_event, max |Δm|) plane). `analyse.py` refuses a record whose state,
seeds, packet count, grids or control differ from this section.

---

## G3 — is the transport-relevant structure predictively tabulable? (preregistered 2026-09-22; not run until the PI approves this section)

G1 showed a small frequency-group operator reproduces energy-conserving
fluorescence where no scalar ε can; G2 showed the structure it exploits is
frequency locality, not global low rank, and that better reproduction of
microscopic events does not imply better reproduction of transported
observables. G3 asks the question those two leave (the PI, 2026-09-22):

> Can the locally coarse transport-relevant structure itself be tabulated
> predictively?

Its three parts keep the announced structure — **state transfer**,
**species composability**, **realistic mixture** — with the PI's two
priorities built in: *fixed and recomputed kernels are tested separately*,
and *the state axes are separated* rather than moved together.

### The two tests, kept apart

For a state θ and the reference state θ₀ (G1's: P1, 2 d, shell 28):

    R(θ₀) → θ    "fixed":      the θ₀ kernel transported at θ   — does it TRANSFER?
    R(θ)  → θ    "recomputed": a kernel built at θ, transported at θ — does the low-dimensional
                                representation still EXIST there?

The failures mean different things and are never merged into one number.
If `R_8(θ)` succeeds while `R_8(θ₀)` fails at θ, the compression principle
survives and what is needed is a state-dependent table; that outcome does
not count against the paper.

### The state axes, separated (Paper III's assumption is not inherited)

Paper III treated the source spectrum as irrelevant; its own results had
Nd II sensitive to the incident spectrum where La II was not, so that
assumption is dropped and the axis is measured. Each axis moves **one**
coordinate from θ₀, every other coordinate held:

| axis | coordinate | grid (θ₀ in bold) | what it changes |
|---|---|---|---|
| T | `T_gas` of the atom | 2500, 3000, **3401**, 4000, 5000 K | the LTE populations, hence which lines carry opacity (42,695 → 202,630 for Nd II) and the branching β |
| D | number density, `n_ion × s` | ×0.1, ×0.3, **×1**, ×3, ×10 | the Sobolev τ scale (τ_max 120 → 12,017 for Nd II), hence escape probabilities and the re-absorption chain |
| J | source temperature `t_core` (the shape of the incident J_ν, decoupled from T_gas) | 2500, **3401**, 5000, 7000 K | which input groups are illuminated, hence which rows of R the light actually uses |
| P | the physical trajectory: P1 at 1, 3, 5 d, shells 29, 27, 26 | 1 d (T 4411, ρ 1.03×10⁻¹⁴), **2 d**, 3 d (2987, 1.97×10⁻¹⁵), 5 d (2341, 7.35×10⁻¹⁶) | all coordinates together, as the ejecta actually moves |

14 states besides θ₀. Ions: La II, Ce II, Nd II, each alone at that state's
own `n_ion`; **Ce II and Nd II decide, La II is reported as a control**, as
in G2. Packets: La II and Ce II 3×10⁵ per seed, Nd II 10⁶; evaluation
seeds 1–3, build seeds 101–103; every other transport setting G1's.

### One frozen support per ion (the PI's amendment, 2026-09-24)

Every kernel of G3 — the anchor, every recomputed operator, every
interpolated one — is laid on **one frequency support per ion, frozen
before the run**: the union over the whole G3 domain of each state's
opacity range with G1's 0.5 % margin, computed by
`paperB/gate3/support.py` and committed as `gate3_support.json`:
La II 1.209×10¹⁴–2.624×10¹⁵ Hz, Ce II 7.129×10¹¹–2.661×10¹⁵ Hz,
Nd II 3.858×10¹²–2.619×10¹⁵ Hz; part (b) uses the three-ion union. The
128 fine groups and every coarse grid are log-spaced on that support, so
every matrix in G3 shares its edges, and "outside fixed support" can only
mean a frequency beyond the union (recorded per state as the clipped
fraction; by construction it is zero). G1's N_g* values were measured on
each ion's own range; G3's are on the frozen support and are recorded as
such.

### Legs per state and ion

| tag | what |
|---|---|
| `R2build` | the downward-macroatom build run **at θ** (build seeds), events pooled; its N_g = 16 kernel on the frozen support is saved for the interpolation legs |
| `R2` | the reference **at θ** (evaluation seeds): the target of every comparison at θ |
| `Arec_ng{2,4,8,16,32}` | kernels built from `R2build(θ)`, transported at θ — the *existence* test; `Arec_ng16` is the fresh compact operator of the four-way reading, and the grid gives K*_rec(θ) |
| `Afix_ng16` | **the anchor**: the N_g = 16 kernel built from `R2build(θ₀)` on the frozen support, transported at θ — the *transfer* test. **N_g = 16 for all three ions**, one controlled representation size across the comparison (Ce II's G1 N_g* is 16; an R_8 failing on transfer could not be told from Ce II being under-resolved before transfer began; Nd II and La II are deliberately over-resolved) |
| `Aint_ng16`, `AintM_ng16` | at the interior point only: the whole-operator interpolation and its matrix-only diagnostic (below) |
| `K128`, `K128build` | the fine matrices at θ, for m_event; `K128` also records the state's clipped fraction against the frozen support |

Metrics are G1's throughout (max and mean |Δm| over the live bands,
max |Δcolour|, m_sed, m_event vs `K128`, the energy identity), against
**`R2` at the same state**, with G1's thresholds 0.10 mag and 0.10 mag.

**"Row never trained" against "outside fixed support", reported
separately.** A frequency absorbed at θ in a fine row the anchor (or an
endpoint) never populated is transported by the kernel's standing
convention for a row with no information — coherent scattering — and its
fraction of the leg's interactions is recorded as `rows_never_trained`. A
frequency beyond the frozen edges is clipped into an edge group and its
fraction is recorded as `clipped_frac` (the energy fraction alongside).
**For the transfer and interpolation legs neither is a gray condition**:
they are the mechanism of a transfer failure, not a defect of the run. A
state whose clipped fraction exceeds 5 % is flagged `coverage-limited`.

### The interpolation test — is the whole operator smoothly tabulable?

Where transfer fails on an axis, the question is whether a *table* over
that axis would work, which is what "tabulated predictively" means. The
effective closure is not the 16×16 transition matrix alone: it carries
the conditional exit-frequency distribution and the deposition channel,
and the feasibility study shows the participating forest and its exit
support change substantially across the domain. So the preregistered test
interpolates **the whole effective operator**, and never touches the
interior state's own macroatom events.

At the interior grid point of each axis (T = 3000 K, D = ×0.3,
J = 5000 K, P = 3 d), bracketed by θ₀ and one other grid state, two legs
per ion, both at N_g = 16 on the frozen support, built from the two
endpoints' saved `R2build` kernels alone:

- `Aint_ng16`, **whole-operator interpolation**, linear in the frozen
  coordinate: the energy and photon transition rows (a row populated at
  both endpoints is interpolated; at one endpoint only, that endpoint's
  row is used and the count recorded; at neither, it stays empty), the
  deposition channel q_dep, and per output group the conditional
  exit-frequency weights on the **union of the two endpoints' exit lines,
  a line absent at one endpoint zero-filled there**, renormalised. The
  energy identity holds exactly by construction; the union size and the
  two endpoint sizes are recorded.
- `AintM_ng16`, the **matrix-only diagnostic**: the same interpolated
  transition rows carrying θ₀'s exit tables unchanged. Its distance from
  `Aint_ng16` isolates what the conditional exit spectrum contributes.

**The coordinates are frozen here:** log T_gas on the T axis, log n_ion
(log ρ at fixed composition) on D, log T_core on J, and log t along the
trajectory P. The interpolation weight is λ = (c − c_a)/(c_b − c_a) in
that coordinate. Passing means the operator's dependence on that
coordinate is smooth enough to tabulate and interpolate; failing means it
is not, and that coordinate needs dense sampling or is not a good
coordinate.

### Species composability and the realistic mixture

- **Part (b), composability.** At θ₀, a La II + Ce II + Nd II atom at the
  P1 partial densities. Legs: `R2` (the blend's own macroatom, the
  target), `Amix_ng{2,…,32}` (`RedistributionKernel.mix` of the three
  single-ion kernels on identical edges, weights from
  `paper3/phase5_mixture/mixture.py::composition_weights`, **no blend
  fit**), `Adirect_ng{2,…,32}` (a kernel trained on the blend's own
  events: the upper bound the mixing rule is measured against). 10⁶
  packets.
- **Part (c), realistic mixture.** At θ₀, the full 13-ion P1 blend
  (`legs.py::atom_for_zone`, II only). Legs: `R2`, `Arec_ng{2,…,32}`, and
  the ε grid of G1 for the scalar comparison. 10⁶ packets. Paper IV's F61
  put a 32-group kernel within 0.05 mag of R2 on this blend at 3×10⁵.

### Readings (Gray first; Ce II and Nd II decide, La II reported)

**The four-way reading per state** (the PI's, retained exactly), read in
this order:

- **A** — `Afix_ng16` passes: the fixed anchor transfers → very strong
  universality.
- **B** — transfer fails but `Aint_ng16` passes (interior points): a
  small tabulated R(θ) → arguably the most useful outcome.
- **C** — interpolation fails (or is not available at an endpoint) but
  `Arec_ng16` passes: the compression survives, but the state manifold is
  more nonlinear.
- **D** — the fresh compact operator itself fails: the actual challenge to
  the closure claim.

Over the domain:

- **C1 — existence.** Green: `Arec_ng16` passes at every one of the 14
  states for both decisive ions, and K*_rec ≤ 8 at two thirds of them.
  Yellow: passes everywhere but not the ≤ 8 majority, or fails at exactly
  one state. Red: fails at two or more states for either decisive ion.
- **C2 — transfer.** Per axis, transfer *holds* if `Afix_ng16` passes at
  every state on that axis for both decisive ions. Green: holds on all
  four axes. Yellow: fails on at least one axis, and on every failing axis
  the whole-operator interpolation passes. Red: fails on an axis where the
  whole-operator interpolation also fails.
- **C3 — a compact sufficient state vector.** An axis is *needed* if
  transfer fails on it. Green: at most two of T, D, J are needed by both
  decisive ions — R_ij depends on a small set of coordinates, not on the
  whole radiation field. Yellow: all three are needed but each is
  interpolable. Red: the whole-operator interpolation fails on two or more
  axes. The axis table is reported whatever the outcome, and **the J axis
  is reported explicitly per ion**, since Paper III's assumption that it
  does not matter is what is being tested.
- **C4 — composability.** Green: the mixed kernel passes both thresholds
  at some N ≤ 32 and its K* is within a factor of two of the directly
  trained kernel's. Yellow: the mixed kernel fails but the direct kernel
  passes. Red: neither passes.
- **C5 — realistic mixture.** Green: K*_rec ≤ 8 on the 13-ion blend;
  Yellow ≤ 32; Red otherwise. Reported with ε*'s error on the same blend.

**Decision.** The manuscript's fourth step (generality) is written as
"transfers across states" iff C2 is Green; as "exists everywhere and is
tabulable in a small state vector" iff C1 Green, C2 Yellow and C3 Green or
Yellow — which, per the PI, is a positive result and not a failure; the
PRL claim is reframed by the PI if C1 is Red. **A neural surrogate stays
out** unless C3 goes Red, i.e. unless the dependence is genuinely too
nonlinear for ordinary interpolation.

**Gray conditions** (per state and ion): G1's 1–4, applied to `R2` and the
recomputed legs only — the transfer and interpolation legs' fallback and
clipping are diagnostics, as above. Every state is read against its own
live-band set (R2 at that state), and the set is recorded. A gray state is
excluded from the counts of C1 and C2 and named in the record.

**Amended 2026-09-24 on the PI's hold of PR #12, before the section
froze:** transfer and interpolation at N_g = 16 for all three ions (was
R_8); one frozen support per ion over the domain; the interpolation test
made whole-operator (transition rows, conditional exit-frequency weights
on the union of endpoint exit lines, the deposition channel), with the
matrix-only version kept as a diagnostic; the coordinates frozen as
log T_gas, log n_ion, log T_core, log t; "row never trained" and "outside
fixed support" reported separately. No axis added, no grid expanded.

**What may not change after the run:** the axes and their grids, the
frozen supports, the coordinates, the packet counts, the N grid, N_g = 16
for the transfer and interpolation legs, the thresholds (G1's), the two
tests kept apart, the readings above. What may: the packet count upward on
a gray (the whole state rerun), and bug fixes with the affected states
repeated.

**Records** `paperB/gate3/gate3_ref_<ion>.json`,
`gate3_<axis>_<state>_<ion>.json`, `gate3_partb_blend3.json`,
`gate3_partc_p1blend.json`, `gate3_verdict.json`, the frozen
`gate3_support.json`; figures `docs/figures/paperB/gate3_*.{pdf,png}`.
`analyse.py` refuses a record whose seeds, packet count, grids or N_g = 16
differ from this section.
