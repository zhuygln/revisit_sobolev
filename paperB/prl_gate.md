# Paper B — gate G1, preregistered

Written before any production run (2026-09-21). The state, the grids, the
live-band rule, the metrics and the thresholds below may not change after
the run. What may: the packet count upward if a Gray condition is noise
(the whole grid is rerun, never one leg), and bug fixes with the run
repeated. G2–G4 are named at the end and are detailed here, each as its
own section, only after G1 is read.

## The question G1 answers

Does the compression of Paper III (F25/F27: a small group-to-group
operator reproduces explicit fluorescence on the same Sobolev opacity)
survive indivisible energy packets and the energy-conserving macroatom?

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
- Packets: 3 seeds (1, 2, 3) × 3×10⁵ per leg. Every leg `packets="energy"`.

## Legs per ion

| tag | mode | what |
|---|---|---|
| `R2` | `sobolev_dmacro` | the reference: the downward energy-conserving macroatom; `collect_events=True` |
| `A2_ng2 … A2_ng32` | `sobolev_group` | R_ij at N_g ∈ {2, 4, 8, 16, 32}: kernel from R2's own events, `RedistributionKernel.from_branching_mc(nu_in, nu_out, w_in, N_g, nu_lo=op_nu.min×0.995, nu_hi=op_nu.max×1.005, w_out=w_out)`, log-spaced groups, discrete within-group exit tables, `rows="energy"`; rows never populated scatter coherently |
| `K128` | kernel only | the 128-group kernel from the same events, not transported: the fine reference of the event-level metric |
| `E0.00 … E1.00` | `sobolev_tla` | ε ∈ {0, 0.05, 0.10, …, 1.00} (21 values): with probability ε a thermal re-emission from the energy-weighted emissivity, otherwise resonant scattering in the same line |

After a pass only (not in the gate): `R2M` (`sobolev_macro`, upward jumps
under W = ½ B_ν(T_zone)) and an `A2_ng8` rebuilt from R2M's events, one run
per ion, as Paper IV's F63 did.

## Metrics (per leg, seed-averaged, with the seed standard deviation)

- **m_event** (the information the binning discards, on events alone):
  the N_g-group energy matrix block-expanded onto the 128-group edges,
  against the 128-group energy matrix from the same events; L1 per row,
  weighted by the row's event count, summed. Zero at N_g = 128 by
  construction.
- **m_sed**: Σ_bins |L_ν − L_ν^{R2}| / Σ_bins L_ν^{R2} over the 200 bins of
  the window.
- **m_band**: max |Δm| and mean |Δm| over the live bands, Δm = m − m^{R2}.
- **m_colour**: max |Δcolour| over the pairs of `sobolev.photometry.COLORS`
  whose two bands are both live.
- **m_energy**: the identity residual of every leg (|E_esc + E_core + E_abs
  + E_dep,lab − E_inj| / E_inj) and `RedistributionKernel.validate_energy()`
  of every kernel.
- **m_complexity**: N_g², the kernel table size in kB (`save` size), events
  per packet, wall time.
- **Live bands**: `paper4/phase3_legs/verdict.py::live_bands` applied to
  R2: finite magnitude, ≥ 1 % of L_bol, brighter than the limits griz 23.5 /
  JHK 21.5 at 40 Mpc.
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

**B2 — ε\* is not nearly as good.** e_ε = ε*'s max |Δm| over the live
bands; e_R = the N_g*-matrix's max |Δm| (R_32's if N_g* is undefined).

- Green: e_ε ≥ 0.20 mag and e_ε ≥ 3 e_R.
- Red: e_ε ≤ 1.5 e_R.
- Yellow: otherwise.

**B3 — the ordering.** For every ion, m_band (max), m_sed and m_event fall
monotonically with N_g within the seed noise.

- Green: monotone for all three metrics on all ions.
- Yellow: a metric non-monotone by more than its noise.
- Red: m_band **rises** with N_g for any ion by more than its noise (the
  self-absorption artefact F25 found and fixed with discrete tables).

**Decision.** Continue toward the PRL iff B1 Green and B2 Green and B3 not
Red. Stop the PRL attempt (a methods result at most) iff B1 Red or B2 Red.
Anything else is Yellow: the PI decides with the numbers in hand, and
nothing in this file is edited to make it Green.

**Gray conditions** (per ion, checked before any reading):

1. fewer than two live bands on R2;
2. R2's seed scatter (`mags_seed_std`) > 0.05 mag in any live band (the
   0.10 mag threshold would then sit within twice the noise);
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
Gray condition 2 fires if the single-ion scatter is 2.5 times worse.

## The central figure

Transport error against effective-model complexity: m_band (max), m_sed and
m_event against N_g² for the three ions, ε*'s m_band as a horizontal line
per ion, the seed noise as a band, `docs/figures/paperB/gate1_error_vs_complexity.{pdf,png}`.

## Records

`paperB/gate1/gate1_<ion>.json` (the `run_legs` row plus the kernels'
validation, table size, empty-row usage, and the grids as run),
`paperB/gate1/gate1_verdict.json` (every reading with the values it was
read from and every gray check that ran). `analyse.py` refuses a record
whose state, seeds, packet count, N grid or ε grid differ from this file.

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
- **G4, σ_R vs σ_nuclear**: on the proxy ensemble that exists in this
  repository (Ye-0.21a at X_lan 0.11 and 0.30, solar_r at 0.11, Ye-0.29a),
  σ_R = the spread of R_{N_g*} − macroatom and σ_nuclear = the spread of
  the macroatom across realisations; stated as a proxy until an external
  nuclear ensemble is run through the same transport.
