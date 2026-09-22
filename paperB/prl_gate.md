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
| `G_k1 … G_k32` | **global family**: rank-k non-negative factorisation W H of the populated rows of the 128-group energy matrix (Lee–Seung multiplicative updates, 600 iterations, initialisation seed 0), rows rescaled to their original sums (energy exact), empty rows kept empty; k ∈ {1, 2, 4, 8, 16, 32}. k = 1 is the fully non-local null: one exit distribution for every input | k | 2·128·k |
| `T_f0.5 … T_f0.999` | **truncation family**: the full 128×128 matrix with each output group's exit table cut to the highest-energy-weight lines carrying fraction f of that group's exit energy (at least one per populated group), both weight tables renormalised over the kept lines; f ∈ {0.5, 0.9, 0.99, 0.999} | 128 | 128² |
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

**H1.** Per ion, N_g* (G1's) and k* = the smallest k in the grid whose
`G_k` has max |Δm| ≤ 0.10 mag over the live bands and max |Δcolour| ≤ 0.10
mag; undefined if none does.

- Green: k* ≥ N_g* (undefined counts as ≥) **and**, at the matched count
  k = N_g*, the global operator's m_event is smaller than the local one's
  (it fits the events better and is still not better on the light — the
  "despite" the PI asked for).
- Red: k* ≤ N_g*/4 — a free mixture of four times fewer archetypes
  reproduces the observables: the structure R_ij exploits is rank, not
  locality.
- Yellow: otherwise (in particular the case k* = 1 < N_g*: the single
  non-local archetype already passes; at G1's numbers this is the only
  non-Green outcome open to Nd II, whose N_g* = 2 puts Red out of reach;
  Ce II with N_g* = 16 carries the sharp test, Red iff k* ≤ 4).
- Overall: Green iff Ce II and Nd II are both Green; Red iff either is Red;
  Yellow otherwise; Gray if either is Gray.

**H2.** Per ion, f* = the smallest f in the grid whose `T_f` passes the same
thresholds; ρ_exit = (exit lines kept at f*) / (exit lines in
`K128build`).

- Green: ρ_exit ≤ 0.10 for both Ce II and Nd II — at least 90 % of the
  stored exit lines are invisible to the observables at the 0.1 mag level.
  Recorded alongside as the event-level counterpart: G1's m_event at N_g*
  (0.16–0.42 at 32 groups, F67).
- Red: no f < 1 passes for Ce II or Nd II — the full table is part of the
  physics the light sees.
- Yellow: otherwise.

**Decision.** The PRL's mechanism claim ("locality, not rank") is written
iff H1 is Green; H1 Red reframes the claim (the PI decides what, with the
numbers); Yellow is the PI's call. H2 sets the wording of "compact": Green
— the operator's information content is the matrix plus the top ρ_exit of
the exit lines; Red — the exit table is physics, not representation.

**What may not change after the run:** the k and f grids, the control
rule, the thresholds (G1's), the archetype-count axis, the NMF settings
(600 iterations, seed 0, rescaled rows). What may: the packet count upward
on a Gray (the whole ion rerun), and bug fixes with the run repeated.

**Records** `paperB/gate2/gate2_<ion>.json`, `gate2_verdict.json`,
`docs/figures/paperB/gate2_locality_vs_rank.{pdf,png}` (max |Δm| against
the archetype count for both families, m_event alongside, the truncation
curve against the fraction of exit lines kept, and every operator on the
(m_event, max |Δm|) plane). `analyse.py` refuses a record whose state,
seeds, packet count, grids or control differ from this section.
