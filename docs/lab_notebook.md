# Lab Notebook

Chronological log of what was actually run, what broke, and how it was fixed.
All entries 2026-08-14 (one long session). Companion:
[results_report.md](results_report.md) for the organized findings.

Commit hashes refer to this repo's history; every entry below is pushed.

---

## 1. Repo bootstrap (`6b0337a`, `3032a47`)

- Empty folder → git repo → GitHub `zhuygln/revisit_sobolev` (public).
- **Incident: NTFS alternate-data-stream leak.** Files copied from Windows
  arrive in WSL with companion files (`*.md:Zone.Identifier`, `:PG$Secure`,
  `:Avecto.Zone.Identifier`) containing the download URL and referrer. A
  `git add -A` swept them into the first push, briefly publishing a private
  conversation URL. Fixed by history rewrite (`--force-with-lease`) and
  `.gitignore` patterns `*:Zone.Identifier`, `*:*.Zone.Identifier`,
  `*:PG$Secure`. **Rule: check for `*:*` stream files after any copy from
  Windows.**
- **Account split:** the `gh` CLI is authenticated as a work Enterprise
  Managed User and *cannot* touch personal repos; plain git over SSH
  authenticates as `zhuygln`. Repo-local identity `zhuygln@gmail.com`
  (global git config carries the work email).

## 2. Phase 0 scaffold (`3397d9f`)

- Package skeleton with two deliberate stubs (`tau_exact`, `load_gsi`);
  failing tests as the Phase 0A specification.
- **Bug caught: hardcoded constant.** First version of
  `SIGMA_CLASSICAL = πe²/(m_e c)` was typed as 3.14e-2 — **18% wrong**
  (correct: 2.654e-2). Caught by re-deriving; now computed from E_ESU, M_E, C
  at import. It would have cancelled in E_Sob ratios but corrupted every
  absolute τ.
- Deviations from the plan doc's §24 layout (documented in the commit):
  notebooks 00+01 merged per §25; `formal_transfer.py` deferred (nothing in
  Phase 0 solves for intensity); both optical depths in one module so the
  prefactor/populations cannot drift apart.

## 3. Phase 0A/0B — toy model (`4bec3fc`, `6d47f17`)

- `tau_exact` implemented; gate test green: constant population reproduces
  τ_S to 2×10⁻¹⁶ at 2×10⁵ grid points; monotone convergence
  (N = 201 → 8.2e-1, 2001 → 1.0e-4, 20001 → 2.2e-12).
- **Physics trap found before it bit:** the planned tanh-gradient experiment
  with the resonance at the tanh center measures *zero* error at leading
  order (even profile × odd gradient; center is also the n″ = 0 inflection).
  Sweep places the resonance at u = 0.5 → clean ε² law (slope 2.0000),
  4.8% at ε = 1.
- **Notebook tooling bug:** nbformat cells written as `source.split("\n")`
  concatenate into one line (list-of-lines needs trailing `\n`s). Fixed by
  passing the source as a single string. Symptom: SyntaxError on execution.

## 4. Phase 0C — GSI data (`45bc18d`)

- Zenodo: the plan doc says "GSI v2"; record 15835361 is v1 (2025-07). The
  concept DOI `10.5281/zenodo.15835360` resolves to **19335084**
  (2026-03, matches the published PRD paper) — used that. 873 MB zip →
  extracted La II only, zips deleted, provenance in `data/README.md`.
- File format is self-documenting (column table between dashed separators);
  `load_gsi` autodetects the header after the *last* separator — works
  unchanged for both transitions and levels files.
- Crowding numbers (3000 Å–3 μm): broader subset 11.3% of lines within
  10 km/s of a neighbour; high-confidence subset 1.2%.

## 5. Phase 0D — SEDONA build (`fbb7d20`)

The no-root build recipe (WSL2, Ubuntu 22.04):

1. `git clone dnkasen/pubsed ~/personal/pubsed`
2. **CRLF trap:** global `core.autocrlf=true` checked the clone out with
   CRLF → `/bin/bash^M` shebang failures. Fix:
   `git config core.autocrlf input && git reset --hard HEAD`.
3. Deps without sudo: micromamba static binary → env at
   `~/personal/sedona-deps/env` with gsl, hdf5, zlib; Lua 5.1.5 from lua.org,
   `make generic` (30 s).
4. **MPI is not optional:** `src/main/sedona.h` hardcodes
   `#define MPI_PARALLEL 1` (the `#ifdef` guards are always true). Serial
   g++ build fails on `mpi.h`. Fix: `micromamba install openmpi`, and since
   the conda wrapper wants conda's compiler, `CXX = env OMPI_CXX=g++ mpicxx`.
5. `Makefile.wsl` (in the pubsed clone) carries all of this;
   `cd src && bash install.sh wsl` → `src/sedona6.ex` (the copy lands in
   `src/`, not the repo root).
6. Runs need `SEDONA_HOME=~/personal/pubsed`. Lightbulb: 2×10⁵ particles,
   1.1 s, spectrum written. Both opacity switches confirmed at
   `defaults/sedona_defaults.lua:115,118`.

## 6. Week 1 (`7ee2f6e`, `9398ea9`, `88524d4`)

- Voigt via `scipy.special.wofz`; normalization tolerance must budget for the
  1/Δν² Lorentzian tail (integrate to 2000 Doppler widths, rtol 1e-3).
- Boltzmann module. **Wrong test assumption caught by the data:** asserted
  the ground state dominates at 5000 K; actually the J=3 level at 1016 cm⁻¹
  out-populates it (g=7 beats the mild Boltzmann factor). The rigorous
  invariant is monotonicity of fraction *per statistical weight*.
- Type Ia example pair: expansion 14 s vs resolved 88 min (**378×**), same 4
  iterations. Ran the resolved one to completion rather than stop-lossing at
  1 h because it was already on iteration 4 of 4.

## 7. Week 2 — solver + minimal model (`786532d`, `eab2245`)

- Formal solver: p-z rays, cumulative-τ formal solution, vectorized per ray.
- **Instructive test failure → Finding F3:** the single-line trough came out
  0.4290 instead of e^−τs = 0.3996 — exactly exp(−τ_S(1−z/ct)) with the
  resonance at 7.7% of c. Observer-frame integration vs comoving Sobolev
  differ at O(v_bulk/c). Documented in the module; tests moved to 10-day
  epochs so test resonances sit below 1% of c.
- Also: first τ test used n0 giving τ_S = 229 (uninformative black trough) —
  and it made the two-line test pass *trivially* via `atol` at zero flux.
  Both re-scaled to τ ~ 1 with an explicit non-triviality assertion.
- Minimal SEDONA model design: T_shell = 2000 K because at 5000 K, Saha
  fully ionizes hydrogen at n ~ 5 cm⁻³ (checked before running); at 2000 K
  the SEDONA gas state shows n_e/n_H = 1e-10.
- **Normalization gotcha:** SEDONA's lightbulb pours the *entire*
  `core_luminosity` into the transport frequency window — absolute L_ν is
  high by 1/f_window (~15× here). Normalize by the red-side line-free
  continuum for comparisons.
- Result: solver 0.1372 / SEDONA bb 0.1420 / SEDONA exp 0.4285 vs
  exp(−(1−e^−τ)) = 0.4212 → Finding F4.

## 8. Week 2 — line ladder (`d48d4c6`, `f120ea2`)

- Custom atom files mirroring `2level_atomdata.hdf5`.
  **Segfault → fix:** the per-element group *attributes*
  `n_ions/n_levels/n_lines` are load-bearing; without them the reader
  segfaults on all four runs.
- **Analytic subtlety:** plane-counting staircases are too shallow; the
  attenuation must be p-averaged over the core disk (planes at z < r_core
  are crossed only by outer rays). After the fix: analytic 0.2405 vs solver
  0.2407 vs SEDONA bb 0.2435.
- Expansion follows its own prediction (0.3257 vs 0.3194) → Finding F5.

## 9. Weeks 3–4 — La II forest + sweep (`49351c5`, `fe7d73c`)

- Window scan (100 Å windows, 3500–9000 Å, ranked by 8th-strongest τ/n):
  **3850–3950 Å** wins; 153 lines, τ distribution 3 > 1 / 6 in 0.1–1 at
  n_ion = 2146 cm⁻³.
- Population control: atom file generated from GSI (all 472 levels, window
  lines only, χ = 1e5 eV). f(A) vs f(log gf) agree to 1.5e-4.
- Solver extended: per-line lower-level populations,
  `(nu0, f, pop_frac)` entries.
- Comparison-script bugs (both mine): plateau window placed beyond the
  shell's velocity extent (no absorption there — looked like "no effect");
  and the 15× lightbulb normalization again. Both fixed → Figure 6:
  solver 0.3549 / bb 0.3426 / exp 0.4965 (**Δ_Sob = +44.9%**).
- Sweep (24 runs): **first attempt failed instantly, rc=14** — a sanitized
  `env={SEDONA_HOME, PATH}` breaks the OpenMPI runtime; pass the full
  environment. Second attempt: 24/24 clean, ~30 s each → Figure 7 and
  Finding F6 (strength floor + v_D wing term).

## 9b. Weeks 3–4 — T sweep + thermal-width frontier (2026-08-15)

- `tsweep.py`: 12 runs. T axis at pinned τ_max = 5 (n_ion rescaled per T from
  fresh Boltzmann fractions); frontier at v_D = 3 and 1 km/s (46k / 140k
  transport bins, `timeout=3000` needed — the 1 km/s pair runs ~35 min each).
- Launched with `python -u` this time: the earlier sweep's `| tail` pipe
  buffered all progress output until exit, which made a live monitor useless.
- Results: T is a weak axis in this window (+40.8 → +46.7% over 2500–5000 K,
  same 3–4 strong lines dominate throughout). Frontier: Δ_Sob flat at ~38%
  from 10 km/s down to 1 km/s — the strength floor survives essentially to
  the thermal width. Cost scales ≈ 1/v_D for BOTH modes (the expansion run
  at 1 km/s took 2305 s): the fine transport grid, not the resolved
  treatment itself, is what costs.
- Note for future sweeps: at T = 5000 K the shell emits non-negligibly at
  3900 Å (B ratio ~0.3 to the core), so absolute fluxes include thermal
  fill-in; same-code differentials stay clean, but solver QC at high T needs
  the emission term compared explicitly.

## 9c. Weeks 3–4 — multi-ion blend (2026-08-15)

- Re-downloaded the GSI calibrated zips (deleted earlier to save disk),
  extracted Ce II and Ce III, deleted zips again.
- **Bug found by real data:** Ce II (odd-electron) writes half-integer J as
  fraction strings (`7/2`); `statistical_weight` assumed numeric J and died.
  Fixed with `parse_j` (accepts both), regression test added. La II never
  triggered this — even-electron ions have integer J.
- Blend: 2,529 window lines (Ce II contributes 2,376), 40 strong, minimum
  strong-line spacing 8 km/s. SEDONA runs ~44 s each (line count barely
  matters; the grid does). The Python solver became the slow leg (~16× the
  single-ion cost) — precomputed in parallel via `solve_py.py` with a cache
  the comparison script picks up.
- Result: Δ_Sob = +14.6%, *smaller* than single-ion (+44.9%) → Finding F7
  (non-monotonic density dependence; full blanketing suppresses the
  relative error). Resolved-legs gap widened to 6.5% — suspected Voigt-wing
  vs Gaussian difference across 2,376 lines; add Voigt to the solver before
  chasing this further.

## 9d. Closing the resolved-legs gap (2026-08-15)

Three rounds, two dead hypotheses, one answer.

1. **Voigt wings — wrong.** Extended the solver with Voigt (γ = A/2π) and a
   `cutoff_widths` option mimicking SEDONA's ±5-width truncation (found in
   `AtomicSpecies_opacities.cpp`; its comment says 20 widths, the code says
   5 — and SEDONA's `line_profile` parameter is declared in the defaults but
   never read, so the profile is always Voigt). Effect on the blend band
   flux: **3×10⁻⁵**. a ≈ 3×10⁻⁵ here; wings are irrelevant at these
   strengths. Also checked Ce II's f(A) vs f(log gf): consistent to 10⁻⁴,
   same as La II.
2. **Path/ray resolution — wrong.** 0.2410→0.2411 over 4→64 z-points per
   Doppler width; 0.2414→0.2408 over 20→320 rays. Converged.
3. **Shell thermal emission — right.** The solver integrates
   S = B_ν(T_shell); SEDONA at fixed T with radiative equilibrium off
   deposits absorbed energy without re-emitting. T_shell→0 in the solver:
   blend 0.2410→0.2249 vs SEDONA 0.2277 (−1.2%); single-ion 0.3538→0.3390
   vs 0.3426 (−1.0%). Like-for-like the codes agree at the MC noise floor.

Lesson worth keeping: the offset scaled with saturation (3.6% sparse →
6.5% blend), which *looked* like a line-count/blending effect and sent me
after profile wings first. The sign was the real clue all along — the
solver was too BRIGHT, and an extra emission source explains brightness far
more directly than a missing opacity would.

Process note: a session restart killed a 30-min solver run (no output file,
no completion record) and reset the shell's cwd, which then broke a
relative `../../.venv/bin/python` invocation (exit 127). Use absolute paths
for long background jobs.

## 9e. The Sobolev-proper leg (2026-08-16, `233db88`, `3fdbe9a`)

The paper had a framing problem: everything measured so far was the
*expansion-opacity implementation*, while the title claimed the *Sobolev
approximation*. Fixing it turned out to cost almost no compute.

- **The enabling realization:** in this pure-absorption LTE setup the
  p-averaged analytic staircase **is** the per-line Sobolev prediction,
  exactly. A Monte Carlo implementation of per-line Sobolev interactions
  would add only noise. So the work was refactoring, not simulation:
  `sobolev_staircase` was promoted out of `experiments/line_ladder/compare.py`
  into `sobolev/sobolev_leg.py`, generalized to accept per-line populations
  and arbitrary geometry, with the ladder's five numbers
  (0.2405 / 0.2407 / 0.2435 / 0.3194 / 0.3257) as the regression check. They
  are unchanged. The same function yields the expansion prediction through
  its `damp` argument.
- **Dead end worth recording.** The first run reported Δ_Sobolev = +7.0%
  *exceeding* Δ_expansion = +2.8% at τ_max = 0.5 — impossible, since the cap
  `1−e^−τ < τ` must always transmit *more*. Cause: I was comparing my
  analytic Sobolev against *SEDONA's* expansion run, not against the analytic
  cap. Adding the analytic expansion column restored the correct ordering
  (0.8199 > 0.7967) and turned the anomaly into a real observation:
  **SEDONA's expansion implementation runs ~6% darker than the pure
  per-crossing cap**, because binning smears absorption into the gaps between
  lines. At τ_max = 0.5 that darkening accidentally cancels the cap's
  over-transmission, which is why SEDONA's expansion looks *better* than
  Sobolev there — compensating errors, not accuracy.
  Lesson: compare analytic against analytic before concluding anything about
  a formalism.
- **Result (F9):** the strength-set floor is expansion opacity's alone
  (Sobolev proper is +5–7% at τ_max = 0.5, 5 and 50 alike), while the
  v_D-growing term is a genuine Sobolev failure (the leg is exactly
  v_D-independent by construction, so all width dependence is unmodelled).
- **Open:** the residual +5–7% Sobolev floor at small v_D is not explained by
  line strength; the candidate is overlap within a Doppler width. Early
  breadth-sweep numbers complicate it further — in the 4300 Å window
  Δ_Sobolev is *negative* (−0.2% to −2.5%), so the residual is not a fixed
  positive offset and may be window-dependent. Flagged, not resolved.

## 9f. Breadth sweep, and a normalization bug caught by a control (2026-08-16)

36 conditions (4 windows × 3 epochs × 3 ion mixes), 72 SEDONA runs, 41 min.
Windows auto-selected by τ-richness: 4300, 4900, 7000, 9100 Å — deliberately
none of them the reference window.

**The bug.** The first pass reported Δ_Sobolev ≈ −7% everywhere, apparently
contradicting §4.12's +5–7%. What gave it away was a *control that should
have been trivial*: the 9100 Å window at τ_max = 0.00 — no absorbing lines
at all — reported a resolved band flux of **1.0746** instead of 1.000.

Cause: SEDONA's final spectrum bin is a partial bin and collapses
(L/L_cont drops 48.5 → 15.4 in that last bin). My red normalization margin
ran to `hi*1.008` while the transport grid ended at `hi*1.010`, so the margin
straddled the bad bin, depressing the reference and inflating every band flux
by ~7%.

Effect by quantity, and this is the reusable lesson:
- **Δ_expansion was safe** — both fluxes carry the same bias and it cancels.
  The line-free controls correctly showed ≈0%.
- **Δ_Sobolev was contaminated** — the analytic leg is correctly normalized
  to 1, so it was compared against an inflated SEDONA number and picked up a
  spurious −7%, which masqueraded as a physical sign flip.

Same-code differentials are robust; cross-code comparisons are only as good
as the shared normalization. That is F8's lesson recurring in a new guise,
and it is now two for two — worth treating as a standing rule.

Fixed by `recompute.py`, which re-derives band fluxes from the saved spectra
with the margin pulled inside the grid (`hi*1.002` → `hi*1.006`); no SEDONA
re-runs needed. Line-free controls now give 0.9943 ± 0.0038. Corrected
Δ_Sobolev is positive and τ-tracking (+0.0% → +11.3%), consistent with the
reference window.

**Practice to keep:** every sweep should include a line-free control
condition. This one was free — the 9100 Å window happens to have no strong
La II lines — and it caught a 7% systematic that would otherwise have been
written up as a finding.

Result: **F10** — the F9 separation is universal, and realized τ_max is the
sole controlling variable (windows, epochs and ion mixes collapse onto one
trend).

## 9g. Bibliography verification (2026-08-17)

The manuscript carried a blanket red TODO: every reference was drafted from
memory. Checked all of them against publisher records. The scoreboard:

- **Wrong**: `metzger2019` — Metzger's *Kilonovae* review is Living Reviews
  in Relativity **23**, 1 (**2020**), not 2019. Key renamed `metzger2020`.
- **Misleading key**: `rothkasen2014` is ApJS 217, 9 (**2015**); the key
  carried the arXiv year. Renamed `rothkasen2015`.
- **Correct but thin**: the other ten were right on author/journal/volume but
  lacked DOIs, arXiv IDs and end pages. All now complete (Karp 161--178,
  Eastman 731--751, Castor 111--127, Kerzendorf 387--404, Tanaka
  1369--1392). Sobolev 1960 gained its translator and the note that it
  renders a 1947 Russian original — the bare 1960 date misleads.
- **Uncited**: five entries were in the file but never cited. Rather than
  deleting them they were placed where the text already implied them.
- **Missing**: the Discussion named "TARDIS-class codes" and the Introduction
  asserted lanthanide opacity dominance, both with no citation at all. Added
  Kerzendorf & Sim (2014), Tanaka et al. (2020, 2025) and Barnes & Kasen
  (2013). The data section now also says *why* GSI was chosen over the
  alternatives (experimentally calibrated wavelengths), which the draft had
  never justified.

**Standing constraint worth recording:** `gsi_atomic` could not be verified
from here — APS returns 403 and ADS blocks automated retrieval — so the
journal volume and year were left as an explicit narrow TODO rather than
guessed. **Resolved 2026-08-18**: the author supplied the published record
(Flörs, da Silva, Marques, Sampaio & Martínez-Pinedo, *Phys. Rev. D* **113**,
063041, 2026), which is now in the bibliography with the full author list,
and the TODO is gone. Note the published title uses a colon rather than the
arXiv version's "I."; the PRD form is what the entry carries.

The constraint itself still holds: publisher pages generally cannot be
scraped, so any citation that matters needs a human check. The useful
pattern was to write down exactly which *field* was unverified rather than
flagging the whole entry — it made the outstanding work a two-minute lookup
instead of a re-verification.

All 15 entries now resolve, no orphans, manuscript compiles clean at 23 pages.

## 9h. Manuscript restructure and documentation sweep (2026-08-17)

**Manuscript.** Retitled from "How Large Is the Error of the Sobolev-Class
Line-Transfer Approximation" to *Two Approximations Under One Name* — the old
title asked a question the paper outgrew once F9 showed there are two
approximations, not one. Results reordered so the separation leads instead of
arriving fifth, and the thermal-emission convention subsection moved out of
Results into Methods, where it belongs: it is a measurement convention, not a
finding about ejecta. Every seam the reorder created needed rewriting — the
separation section had opened "Everything above measures the implementation,"
which became false the moment it moved to the front.

While checking the reorder I found the Discussion still promising a per-line
Sobolev leg as future work, months after building it, and Reproducibility
still claiming 19 tests. Both fixed.

**Docs sweep.** The surrounding documentation had drifted behind the work
twice in three days. The worst case was a genuine contradiction rather than
staleness: the report's next-steps section still carried the *retracted*
claim that Δ_Sobolev is negative (−0.2 to −2.5%) — the normalization artifact
of §9f — while §4.13, in the same document, explained that it was an
artifact. A reader arriving at the end would have been told something the
middle had already withdrawn.

Also fixed: neither the README nor the report linked the manuscript at all,
despite the paper being the project's main deliverable; the README had grown
into a chronological diary (Phase 0 / Week 1 / Week 2 / Weeks 3–4) with
Layout and Setup buried 150 lines down, still opening with the
single-approximation question F9 superseded, and still listing 7 tests when
there are 27.

**Practice worth adopting:** update the report, notebook and README in the
*same commit* as the result they describe. Both drift episodes came from
deferring documentation to "after the next experiment," and both produced
statements that were not merely out of date but wrong.

## 9i. The A.6 derivation was wrong (2026-08-18)

External review caught it, and it is the most instructive error of the
project so far.

**What I did wrong.** Appendix A.6 derived the relativistic Sobolev depth by
differentiating `D = γ(1−b)` **at fixed t** — a frozen snapshot. But β = r/(ct)
in homologous flow depends on the ejecta age, and a photon takes time to
cross a resonance, so along the worldline `db/dt = (1−b)/t` rather than
`db/dz = 1/(ct)`. The difference is a factor (1−b): **first order**, the same
order as the effect I was computing.

The trap: I reasoned that the crossing time is a negligible fraction of t
(~v_D/c ~ 3e-4) and therefore that time-dependence was ignorable. Wrong test.
What enters is the *gradient*, and the time term is proportional to b itself.

**Settled numerically before writing any new code**: integrating the resolved
opacity along the path, with and without the age advancing, reproduces
(1−β)/γ and 1/γ respectively to five decimals. My published formula equals
(1−β)/γ exactly for radial rays — right algebra, wrong problem.

**An over-claim I made and retracted within the hour.** Switching the solver
default to worldline moved its forest agreement with SEDONA from −1.04%
(frozen 1st order) through −0.53% (frozen exact) to −0.06% (worldline), and I
announced that as an independent time-dependent code confirming the worldline
law. **Wrong.** At β ≤ 0.01 the two laws differ by only ~1–2%, comparable to
other systematics; I read an ordering out of noise-level differences because
it agreed with the conclusion I wanted.

The v/c sweep, with 23× more leverage, says the opposite: SEDONA's τ_eff
tracks the *frozen* law (RMS 0.024) not the worldline law (RMS 0.552), and at
β_res = 0.34 gives 0.613 against frozen 0.618 and worldline 1.43. The cause
is in the SEDONA source: `transport_steady_iterate` sets `use_hydro_ = 0` and
reads no time-stepping parameters. In this configuration it is a
frozen-snapshot calculation.

So the worldline law stands on the first-principles integration alone, which
is solid; SEDONA neither confirms nor refutes it, because it solves the other
problem. The real finding is that **transport treatment is a third convention
that cross-code comparison must match**, alongside thermal emission (F8) and
spectral normalization (§4.13) — and it is the one with the largest lever arm
at high velocity.

**Consequence, opposite to the campaign premise.** The whole mechanism-
isolation campaign was motivated by the idea that an O(v/c) effect might
explain the +5–11% Sobolev residual. The physical law has *no* O(v/c) term,
so at β ≤ 0.01 the correction is 5e-5 and explains none of it. The hypothesis
is dead and overlap is now the only live candidate.

**Test-design trap, second occurrence.** Worldline anchoring makes the ejecta
age grow along the ray, so b = z/(c t_exp + z − z0) saturates near 0.26 for
the test shell; a nominal β = 0.3 probe then falls *outside* the shell,
absorbs nothing, and every comparison passes trivially. Same failure mode as
the earlier high-β probe. There is now an explicit test asserting the probes
absorb.

**Also checked and cleared:** time anchoring (ray-start vs centre-plane)
changes results by 0.4%, so it is not behind the single-line discrepancy
noted in report §4.14.

**Process note.** Three readings of F3 — frame ambiguity, leading relativistic
correction, frozen-snapshot artifact — and the second was committed, pushed
and written into the manuscript before review caught it. Then, in the same
session, an over-claim about SEDONA confirming the fix, also pushed, also
retracted.

Two distinct lessons. (1) A derivation should be checked against a
first-principles numerical integration *before* it goes into the paper; that
check took minutes and would have caught the frozen error. (2) A cross-code
agreement is only evidence if the effect exceeds the systematics — I inferred
a law from ~1% orderings when the discriminating experiment, already running
in the background, had 23× the leverage and pointed the other way. Wait for
the high-leverage measurement before announcing the conclusion.

## 9j. The residual was mine, not Sobolev's (2026-08-18)

Two reviewer arguments, both right, and the second exposed a bug I had
already fixed twice elsewhere.

**Overlap can't be the cause, provably.** In pure absorption with fixed
populations opacities add, so tau = sum tau_i exactly and the emergent
intensity is exp(-sum tau_S,i) whether profiles overlap or not. My own
two-line unit test had asserted this since Week 2. I had been carrying
"overlap / isolation failure" as the leading hypothesis in the manuscript,
the report and the README for days without noticing that the harness cannot
express it. Numerically confirmed: two identical lines from 20 Doppler widths
apart to exact coincidence reproduce Sobolev to six decimals.

**Third instance of the normalization bug.** Asked to compute the Sobolev
residual for the v_D = 1 and 3 km/s runs that already existed, I found that
`sweep.py` normalized SEDONA band fluxes by RAW LUMINOSITY in the red margin
instead of the continuum ratio -- leaving the Planck slope across the band in
the answer, while the analytic leg was correctly normalized. That single
mismatch manufactured the entire "v_D-independent 5-8% Sobolev floor" that
three documents had been trying to explain.

Corrected, Delta_Sobolev is -0.3% at 1 km/s and -0.2% at 10 km/s, rising to
+9.2% only at 300 km/s. Delta_expansion moved by <2 points anywhere, because
same-code differentials cancel the bias -- exactly as in F8 and the breadth
red-edge case.

**Fix, structural this time.** Added `sobolev/spectra.py` with ONE band-ratio
routine and six tests, including a null spectrum that must return exactly 1.
Every experiment now calls it. Three occurrences of one bug class is enough
evidence that ad-hoc per-experiment normalization was the real defect.

**What the residual actually is (F12).** A finite-region boundary effect:
within a few Doppler widths of an edge the resolved profile is clipped while
Sobolev applies a step, giving an erf law confirmed to four decimals, with
band-averaged size ~v_D/Delta v_shell. It is a property of imposing a finite
line-forming region, not of the approximation, and at real thermal widths
(0.6 km/s against 10^4 km/s spans) it is negligible.

**Standing rule, now earned three times over:** same-code differentials are
robust; every cross-code comparison must share one normalization convention,
and every sweep should carry a null control that is required to return unity.

## 9k. Bin-width control, and a verdict I nearly published (2026-08-18)

Ran the P1-B control: does Delta_expansion depend on transport bin width? If
it does, the central claim is a usage artifact and Paper I needs reframing.

**The script told me it did.** It printed "NOT INVARIANT -> the error depends
on usage; reframe the claim", computed over all seven resolutions.

It was wrong, and the tell was in a column I had put there for context rather
than as a check: F_resolved. Stable at 0.3416-0.3420 for five points, then
0.4415, then 0.6789. Once the bin exceeds the line profile the RESOLVED leg
stops resolving; those points measure my reference failing, not expansion
opacity changing. Restricted to where the reference is converged
(bin <= 1.5 v_D), Delta_expansion is +43.7% +- 1.2 points across two decades.

**Closing the objection properly.** The valid range only reaches ~2.5
lines/bin, because resolving profiles (bin << v_D) and filling bins with lines
(bin >> spacing) pull opposite ways. The fix is denser lines, not coarser
bins: the La+Ce blend has 2529 lines in the same window and reaches 41
lines/bin with the reference still converged. Delta = +15.3% +- 0.7 there,
flat. Objection closed in both directions.

**Two lessons.**
1. Establish the reference's own convergence before interpreting its
   disagreement with anything else. I have now been bitten by the reference
   being wrong three times (frozen transport, normalization, and this), and
   each time the symptom looked like a property of the thing under test.
2. An automated verdict is only as good as its validity criterion. The script
   now excludes unconverged points AND prints them with the reason, because a
   silent exclusion is indistinguishable from cherry-picking.

## 9l. Paper II Phase 0: both reference codes were wrong for the job (2026-08-18)

Paper II needs a code where a photon absorbed in one line can leave in
another. I went looking for one. Neither candidate works.

**SEDONA (P2-0A).** I expected to find fluorescence switched off and to switch
it on. It is not switched off — it does not exist. Outside `sandbox/` the only
hit for fluor|branch|macroatom|downbranch in the whole source tree is a comment
about nuclear decay. `opacity_epsilon` is one global scalar, and an interaction
is either coherent scattering or a redraw from the zone's thermal pool. Neither
channel knows which line absorbed the photon, and the expansion-opacity path
throws that identity away by construction — the same information loss that
produces F4. Implementing branching is core packet-interaction work.

**TARDIS (P2-0D recon).** The opposite problem. `downbranch` and `macroatom`
are configuration-level options and the code is there. It installs — but only
via the repo lockfile: conda-forge has no package at all, and `pip install
tardis-sn` pulls a 2015-era stub that dies calling `ez_setup.py`. Then all
three modes fail identically at atomic-data load, before any transport. Its
bundled downloader 404s because `tardis-regression-data`'s LFS objects are gone
server-side (I checked three files — repo-wide), and the one file still
retrievable via LFS anywhere, from `tardis-atomdata`, is `database_version
v0.9`: plain HDF5 datasets, not the pandas tables the current reader wants.
`pip install carsus` doesn't exist either.

So I learned nothing about TARDIS's modes. That is worth stating plainly rather
than filing the install as a success: the blocker is upstream data
distribution, and all I verified is that the package imports.

**The decision.** Build the instrument. Not just because both candidates were
blocked, but because of the shape of the problem: TARDIS has branching and no
expansion opacity, SEDONA has expansion opacity and no branching. Comparing
them would be cross-code — the exact thing that has burned this project three
times. One code that varies both is what Phase 1 actually needs.
`paper2/phase0/three_level_atom/branching_mc.py`, ~250 lines.

**Calibration.** Branching yields match A₃₂/(A₃₁+A₃₂) to 0.96σ across five
ratios; interaction probabilities match 1−e^(−τ_S) to 1.41σ over τ = 0.1–10.
Those were never in much doubt. The real check was the pure-absorption
spectrum against `sobolev_attenuation` — Paper I's analytic leg, written
independently, reaching the same answer by integrating over impact parameter
instead of sampling rays. 2.07σ worst bin, max absolute difference 0.0070.

**And it failed the first time, at 27σ, in exactly two bins.** I was comparing
a bin-averaged MC against a midpoint-evaluated analytic. The trough edges are
near-vertical — the resonance plane leaves the shell over a fraction of a
percent in frequency — so at the edges midpoint and average differ hugely while
every other bin agreed to under 1%. Averaging both sides over the bin fixed it
completely.

That is the third time (§9k, §4.13, here) that an alarming disagreement was the
two sides being asked slightly different questions. The tell is the same every
time: the discrepancy is confined to where the quantity varies fastest, and it
is large where it appears at all rather than spread thin. I should check that
before checking the physics — and it is now written into the test's docstring
so the next person doesn't rediscover it.

Also worth recording: my first draft of the interaction-probability test
launched packets over a band 97% of which had no resonance in the shell at all,
and reported "too few packets to measure". The fix was to derive the band from
the geometry — require r_core < z_res < sqrt(r_out² − r_core²) and every ray
crosses, whatever its impact parameter — and then assert that the offered
fraction is exactly 1. A denominator you assume is a denominator you get wrong.

## 9m. The term I called small was the big one (2026-08-20)

Started the realistic-velocity work (open item 2) by quantifying what the
solver's worldline mode does not model, rather than by building anything.
That ordering paid immediately.

**My first estimate was wrong, and wrong in the safe direction.** I had the
light-travel factor as t_res = t0 (1 + beta): photon flies to the resonance
radius beta c t0, taking beta t0. That treats the resonance as sitting still.
It does not -- the fluid element is moving outward at beta c while the photon
chases it, so

    t_res = t0 / (1 - beta),

1.43 t0 at beta = 0.3, not 1.30. With n ~ t^-3 and tau_S ~ n t, the medium
contributes (1-beta)^2 and the transport 1/gamma:

    tau / tau_S(t0) = (1 - beta)^2 / gamma

Confirmed to five decimals against direct integration of the resolved opacity
with the density evolving along the path. At beta = 0.3 the 1/gamma deficit is
4.6% and the dilution deficit is 51% -- **the effect this project has been
calling "the relativistic correction" is 11x smaller than the one it was not
modelling at all.** Undiluted worldline transport gives 1.363; the physical
answer is 0.467; a factor 2.9 apart, and on opposite sides of 1.

**The hole was in the ground truth, not just in the solver.** `_ground_truth_tau`
anchors its clock at the resonance, so it measures tau/tau_S(t_res) and is by
construction blind to how the medium got there. It validated 1/gamma perfectly
and could never have caught this. The fix was a second integrator anchored at
the emission epoch -- which is the epoch a spectrum is actually labelled by.
Lesson, and it is the same one as sections 9k and 4.13 in a new costume:
**check what your reference is normalized against, not just whether it agrees.**

**The 0.26 saturation was this bug's fingerprint all along.** Section 9i
recorded that b = z/(c t_exp + z - z0) saturates near 0.26 for a beta_out = 0.35
shell and filed it as a test-design quirk. It is exactly beta_out/(1+beta_out),
and it is the frozen outer boundary: the photon chases a wall left behind at
t_exp. With r_out co-moving it goes away. I had the number, wrote it down, and
read it as an inconvenience rather than as evidence.

**What changed.** `emergent_luminosity` gains `dilution=`, which reads
`n_l_of_r` as the initial profile and dilutes it off the Lagrangian radius
(homology makes beta the label, so the element now at (r,t) sat at r t_exp/t),
and lets the outer boundary recede as r_out(t) = beta_out c t. Measured against
the frozen medium at fixed resonance velocity it reproduces (1-beta)^3 to ~2%,
the residual being the finite core.

The default is `None`, not False. Worldline mode above beta_out = 0.02 now
raises rather than picking one, because a silent 50% error is precisely this
project's characteristic failure and I would rather break a caller than repeat
it. Passing False explicitly is still allowed and the transport-law tests do
that, since holding one medium fixed across all three modes is how they isolate
the transport treatment.

**Also settled: SEDONA can do this, and we never asked it to.** The manuscript's
future-work line read as though the worldline comparison needed a code that
does not exist. It does not: `transport_steady_iterate` gates off hydro
(`SedonaClass.cpp:277`), and with it unset `hydro_homologous::step` divides rho
by e^3, expands the grid and advances t_now, while photons in flight are
carried across timesteps -- the worldline problem at timestep resolution.
Doppler and aberration are already exact. All 165 param files in `experiments/`
set `transport_steady_iterate = 1`, so every one of them silently discarded its
own `hydro_module = "homologous"`. Manuscript amended in both places; running
it is Stage 5.

**The analytic leg, and a geometry objection that dissolved.** Both Paper I
legs come from `sobolev_attenuation`, which had no relativity at all -- a
first-order plane z_res = c t (1 - nu0/nu) and a beta-free tau_S -- so no Delta
could be measured at high beta. Relativizing it turned out to be easier than
feared, for a reason worth recording.

Under worldline transport, with ct(z) = z + Z0,

    D = gamma (1 - beta_z) = Z0 / sqrt(Z0^2 + 2 Z0 z - p^2),

so D = nu0/nu is LINEAR in z: z_res = (Z0/2)(y^2-1) + p^2/(2 Z0), one root, no
discriminant, no branch to choose. And since |dD/dz| = D^3/Z0 with
D(z_res+Z0) = gamma Z0, the impact parameter cancels identically --
tau/tau_S(t_res) = 1/gamma for *every* ray.

The two-root quadratic behind Jeffery's CD/CP surfaces belongs to the FROZEN
problem. The manuscript flagged the non-planar locus as a limitation and
`vc_control` dodged it with a 1e-4 core; in the mode that is physically correct
the question does not arise. The frozen branch still carries both roots, and a
test asserts the upper one never falls inside a shell below 0.35c -- which is
the actual reason the plane picture was adequate for Paper I, now written down
instead of assumed.

`damp` is still applied to tau per crossing, so the expansion leg remains the
same code path with one keyword flipped. The classical path is untouched and
guarded by array_equal against hard-coded outputs, not allclose: adding the
modes must not move a single ulp. 104 tests green.

**Still open:** the Delta measurement itself (Stage 6), and running SEDONA's
time-dependent mode (Stage 5).

## 9n. A pilot at 0.1-0.3c, and two tables that meant nothing (2026-08-20)

With the analytic leg relativized (9m) the Delta measurement is possible for
the first time. Ran a pilot: synthetic 20-line forest, uniform density, one
epoch, tau_max = 5. `experiments/highbeta/pilot.py`.

**Two confounded tables came first, and both looked fine.** Kept in the script
under `--confounded`, because the failure mode is the interesting part.

The first scaled the shell but not the line list. At beta = 0.3 each line's
resonance region covers 30x more velocity space than at 0.01, so the forest
goes from sparse to fully blanketed -- and F7 says Delta is non-monotonic in
forest density, maximal for sparse forests and suppressed at full blanketing.
The table showed Delta_expansion collapsing from +51% to +1.6% and I nearly
believed it. It is F7's axis, measured by accident.

The second fixed that by scaling the wavelength window with the shell span,
but left v_D at 100 km/s. Then v_D/delta_v_shell varies 30x across the sweep
and the low-beta rows carry F12's finite-region boundary effect. Delta_Sobolev
read +19% at beta = 0.01, against Paper I's <=2%, which is the tell I should
have caught immediately.

**Standing rule, extended.** A sweep in one variable has to hold fixed every
dimensionless ratio that any earlier finding identified as controlling. Here
that is four of them -- tau_max (F10), lines per crossing (F7),
v_D/delta_v_shell (F12), beta_core/beta_out -- and I had a finding for each,
already written down, and still built two sweeps that violated them. The
findings list is not just a record; it is the checklist for designing the next
experiment.

**The reference was unconverged, again.** Third time. `emergent_luminosity`
averages over n_impact//2 core rays while `sobolev_attenuation` averages over
n_p = 200, so the resolved leg is the coarse one and Delta drifts as its
p-sampling settles:

    n_impact     12       24       48       96
    D_Sobolev  -0.93%   +0.08%   +0.50%   +0.59%     <- changes SIGN
    D_expans. +15.32%  +16.49%  +16.97%  +17.08%     (beta = 0.01)

I had already reported "D_Sobolev flat at -0.9%" before running this. It is
+0.6%. Frequency resolution, by contrast, was converged at n_nu = 200 already
(stable to +-0.12 points out to 3200) -- because the structure scale across
the band is the shell span, not the Doppler width. I had assumed the opposite
and printed a diagnostic saying 5000 points were needed.

**What saved it: the difference converges even when the endpoints do not.**
The ray error is beta-independent, so it cancels:

    n_impact        12      24      48      96
    d(D_expansion) +4.40   +4.57   +4.51   +4.58
    d(D_Sobolev)   -0.04   +0.06   -0.00   +0.05

An 8x refinement moves each endpoint by ~1.8 points and the difference by
0.18. Worth remembering as a design principle rather than a lucky escape: if
the quantity of interest is a difference, converge the difference, and say so.

**The result, at n_impact = 96 (converged):**

    beta      D_Sobolev    D_expansion
    0.01        +0.59%       +17.08%
    0.30        +0.64%       +21.66%

  1. Delta_Sobolev is beta-independent to 0.05 points. The Sobolev
     approximation proper does not degrade at realistic ejecta velocity.
     Its value here (+0.6%) also sits where Paper I says it should for
     tau_max = 5.
  2. Delta_expansion GROWS, by +4.58 points, ~27% relative. Expansion opacity
     is worse at the velocities kilonovae actually have, not better.

**The collapse hypothesis is refuted.** The plan predicted Delta_expansion at
high beta would land on Paper I's Delta(tau_max) curve evaluated at
tau_max^corr = tau_max (1-beta)^2/gamma -- no new mechanism, just a shift
along a curve already measured. Control at converged ray count:

    beta = 0.01 at tau_max = 2.337 (= 5 R):  D_exp = +16.06%
    beta = 0.30 at tau_max = 5.000:          D_exp = +21.66%
    ratio 1.349, against a predicted 1.000

Only about two thirds of the growth is the tau_max shift. Same ratio at 12
rays (1.374), so this is not a convergence artifact.

**Mechanism unidentified, and my first guess is already weak.** I suspected
the R gradient: R = (1-beta)^2/gamma varies from 0.81 to 0.47 across the band
at beta_out = 0.3, so relativity stretches the tau distribution rather than
rescaling it. But using the band-mean R instead of the edge value only moves
the ratio to ~1.30. A second candidate, untested: lines per crossing was held
fixed in nu_0, while the nu <-> beta_res map is nonlinear at high beta, so the
realized forest density may still drift across the band. Test that before
writing anything down as an explanation.

**Status.** Pilot only -- synthetic forest, uniform density, one epoch, one
tau_max. Not the Stage 6 measurement and not a numbered finding. Before it
becomes one it needs the La II forest and a tau_max range, and the first thing
to fix is the p-sampling mismatch between the legs: matching them removes most
of the cost, which is currently 10-20 minutes per converged point.

## 9o. The referee was right, and the numbers already knew (2026-08-21)

Major revision. The referee's central objection -- that the paper labelled
the difference between a deterministic attenuation and a statistical closure
as proof the closure is "wrong" -- is correct, and testing it on the La II
forest took an afternoon because every tool needed was already in the repo.

**The reframe, verified.** Per ray the Sobolev leg exponentiates S = sum tau_k
and the expansion leg exponentiates E = sum (1 - e^-tau_k). E is exactly what
alpha_exp integrates to -- the expected number of line interactions per
crossing for a photon counted but not removed, Karp's mean-free-path
statistic -- and the closure preserves it identically. Transmission is e^-S
(Bernoulli product) against e^-E (Poisson with the same mean). The gap closes
as tau -> 0 (+0.2% at tau_max = 0.1) and is +38.2% at tau_max = 5, the paper's
number. F_exp/F_Sob = <e^D>_w with w ~ e^-S, D = S - E, holds to 1e-13. F5,
F6, F7, F13 are all corollaries. The referee handed over a better thesis, not
a weaker one, and it is what the numbers had been saying.

**Three things the audit found that the referee did not know.**

1. The breadth Delta_Sob row (+3.6% median, +11.3% max) was STALE:
   `breadth/recompute.py:52` still normalized by raw luminosity -- the bug
   1e2ba21 fixed everywhere else and which retracted a "5-8% Sobolev floor"
   of the same sign and size. The breadth directory never migrated to
   `band_ratio`. Delta_exp is a same-code differential and unaffected.
   Pending: the v2 recompute.
2. The single-line "discrepancy" was the frame convention, exactly. The
   solver in frozen mode gives 0.1372 under every variation; the frozen law
   (1-beta)/gamma over the trough window is 0.1371. So the published 0.1372
   was always the frozen-mode number and 0.1353 is the beta -> 0 limit, not
   the target. SEDONA's 0.1420 was on a grid with < 2 bins per Doppler width
   and 2e6 packets; 3.2e7 packets on the SAME grid give 0.1386 +- 0.0002, so
   most of the 3.5% was sampling. The ladder is running.
3. "Partially correlated" was unjustified as written -- production runs were
   seeded from time(NULL) minutes apart -- but seed-matched pairs correlate
   at +0.97 and the paired Delta_exp scatter is 0.11% against 0.40% from
   quadrature. Matched seeds turn the claim into a measurement.

**The p-sampling mismatch, and why Delta_Sob could change sign.** Attenuation
is a step function of impact parameter; the analytic leg sampled it on 200
midpoint rays and the solver on n_impact//2, and the two O(1/n) step errors
did not cancel. With a shared RaySet (midpoint, not Gauss -- matched nodes
make the step error common) Delta_Sob on the forest is flat to 0.014 points
from 100 to 1600 rays. And it is +3.09% in ALL FOUR transport modes (first,
classical, exact, worldline+dilution): the mode dependence cancels in the
matched pair, which is the whole argument for matching.

**The erf leg.** The resolved calculation for uniform n_l, Gaussian profile,
pure absorption, first-order Doppler is closed-form: the boundary/run.py erf
bracket per line per ray, no z loop, cost independent of v_D. It agrees with
the brute-force solver on identical rays to 8e-5 and makes the 1 km/s
frontier free. Against it: Delta_Sob = +3.1% at 100 km/s, +0.9% at 30,
+0.30% at 10, +0.03% at 1 -- the F12 law -- and Delta_exp converges to +38.1%
at thermal widths, the pure Poisson number. SEDONA resolved validates the
reference to +-0.3% across the 18-point grid (max 0.7%); the "0.1%
agreement" was a single-seed coincidence.

**SEDONA's expansion leg is +2.3% above the analytic Poisson prediction**
(0.4955 vs 0.4845) and carries a ~2% bin-width systematic at fixed reference.
So +38% is the pure count-vs-transmission number and +45% is what the code
delivers. Both are quoted; the referee should not have to find that.

**Stimulated emission** is not negligible in the red breadth windows: 5e-3 at
9100 A and 3000 K, 2.5% on a tau ~ 5 line. SEDONA applies it; the repo did
not. Now folded into pop everywhere.

**RE, first look.** With radiative equilibrium on and T converged (N = 15),
the resolved band fills from 0.343 to 0.856 and expansion to 0.896: a +4.5%
differential against +45% in pure absorption. N = 1 and the other seeds
pending.

**Final numbers (same day, later).** Breadth v2: Delta_Sob median +0.26%
(det ref) / +0.34% (SEDONA ref), max +7.6 / +7.8%, from +3.65 / +11.3%;
Delta_exp unchanged. Seeds: 90 runs, headline 10 pairs +44.90% +- 0.04,
corr +0.95. RE: +7.7% (one iteration) and +5.0% (T converged) from +44%;
band fills 0.34 -> 0.85; flux reappears at 3955-3995 A. Predictor figure: 54
conditions on the 1:1 line. Ladder: 59/59, every rung within +0.5-1.4% of
the frozen target; anchor 0.1383 +- 0.0001; the zero-opacity control reads
1.0079 -- so the +0.9% residual is SEDONA's continuum, not the line, and
corrected the trough is 0.1372 = the frozen law to 0.1%. I had written a
speculative sentence attributing the residual to the comoving core emission
before the control came in; the control said the same thing with a number.
Run the null before writing the explanation.

**Process.** A thing worth keeping: I read the referee's Comment 1 as an
attack on the headline, spent ten minutes deciding whether to rebut it, and
then ran the test. The test took less time than the deciding. Run the test
first.

## 9p. Paper II Phase 1: the atom has to be the whole atom (2026-08-21)

Built the vectorized successor to the Phase 0 instrument and pointed it at
La II. Three things happened on the way to the first number, each the kind
this notebook exists to record.

**The atom decides the experiment.** Sizing the ion first: at 3000 K only 19
of 472 levels are populated, so opacity is sparse (949 lines with tau > 1e-3,
1148-23,363 A) but fluorescence channels are not -- the 9 upper levels that
carry the 3850-3950 A forest decay through 71 lines, 16% of them back into the
window by A. So the window atom cannot be the instrument: fluorescence leaves
through lines the window doesn't have. And Paper I's "red margin is line-free"
was true of Paper I's atom only: the full ion has strong lines at 3949, 3988
and 3995 A, and with them the 3955-3978 A band is absorbed to 0.01. Internally
consistent then; worth knowing now.

**Window-confined re-emission is an artifact, and so was Paper I's RE
fill-in.** First smoke test: thermal re-emission confined to the SEDONA window
on a full-atom launch gave a band of 9.9 -- the whole ion's absorption dumped
into 230 A. SEDONA's RE only made sense because launch, opacity and
re-emission were all confined to its transport grid. Reproducing it requires
the window atom, the window launch and the window emissivity together; the
physical version (whole-ion LTE emissivity at 3000 K) is mostly infrared, and
the band refills from 0.348 to 0.355 rather than to 0.87. So Paper I's RE
"+5-8%" is an upper limit on fill-in, now scoped as such in the manuscript.

**The escape probability.** With the window atom the MC's thermal legs gave
the right totals but the wrong shape: too blue, too spread (blue wing 0.71 vs
SEDONA 0.42). A photon re-emitted inside a resonance zone has still to cross
the rest of its own line's profile and escapes only with Sobolev's
beta = (1 - e^-tau)/tau; otherwise it is re-absorbed by the same line and
re-drawn. For pure resonant scattering that changes nothing observable --
same place, same frequency, isotropic either way -- which is exactly why
Phase 0's calibration could not see it. For thermal and branching
redistribution it decides which line the photon finally leaves through:
trapped in strong lines until a draw lands on a weak one. With it in, the MC
tracks SEDONA's RE spectrum sub-band by sub-band; the expansion leg to <1%,
the Sobolev leg 3.7% high -- the known Sobolev-vs-resolved offset at 100 km/s
carried through. A calibration leg against a code that does the resolved
physics found the missing term in one comparison. The trapped fluorescence
yield b/(1-(1-b)(1-beta)) is now a test.

**The numbers.** Paper I's atom: fluorescence refills the band by only
+12.8% (0.348 -> 0.392); 51% of absorbed photons leave redward. Full ion,
Planck 6000 K launch: Sobolev absorb 0.183, closure absorb 0.344 (+88%),
Sobolev + thermal 0.257, closure + thermal 0.412, Sobolev + fluorescence
0.660. So fluorescence from UV pumps refills the optical band by +260%, and
an expansion-opacity code -- whose only redistribution is thermal -- lands
37.5% BELOW the physics. Too bright in pure absorption, too faint with
fluorescence. That sign flip is the result, and it was not the one I was
expecting: I had assumed redistribution would shrink the closure's error.
It changes its sign instead.

**Caveats I want visible:** frozen LTE at 3000 K, one ion, first-order
Doppler, photon packets not energy packets, an opaque core that eats ~12% of
re-emission (SEDONA loses the same), and a Planck incident continuum standing
in for a photosphere. NLTE and T feedback are Phase 2.

**Next day (2026-08-22), three corrections before any new physics.**

*E0 -- the closure's own emissivity.* Review caught that the expansion legs
applied the emitting line's beta after a thermal re-emission drawn from the
Sobolev line emissivity A n_u. The continuous absorber re-absorbs through its
own bins and SEDONA's expansion mode has no beta, so that was a double count
-- but switching beta OFF made the agreement with SEDONA worse (blue wing
0.77 vs 0.63), and bin-uniform placement alone worse still (0.82). The real
difference was the emissivity: Kirchhoff for the closure's OWN opacity is
kappa_exp B_nu per bin, which saturates at (1 - e^-tau) per strong line; the
A n_u weighting gave the strongest lines ~3x too much re-emission. With
kappa_exp B_nu, bin-uniform, no beta: 0.632 / 1.007 / 0.905 / 1.085 against
SEDONA's 0.632 / 0.998 / 0.900 / 1.083. The "beta on" agreement had been a
coincidence. Headline moved 0.412 -> 0.408, -37.5% -> -38.2%.

*E1 -- energy.* Packets stay photons; every packet now carries h nu, and the
identity E_inj = E_esc + E_core + E_abs + E_dep_lab closes to roundoff in
every mode. The deposit splits into the comoving exchange (which for branching
equals the level-energy difference per chain, 1e-12 on the toy atom, 2e-5 on
La II -- the GSI wavelength-vs-Ritz consistency) and an O(v/c) Doppler work
term. The thermal legs at fixed T are bookkept, not conserved; that is why
SEDONA iterates T.

*E7 -- I had the pumps wrong.* "Fluorescence from UV pumps" is withdrawn: the
6000 K Planck photon budget in 1142-2500 A (0.32% of launches) is smaller than
the band's refill, and the census says 0.00% of the band's escaped energy
came from there. It is pumped from 3300-4500 A -- 28% from 3300-3800, 24%
in-band, 37% from 3955-4500 through deeper lower levels -- by 5d6p upper
levels exiting the strong ground-connected lines after 2-8 re-absorptions.
971 pathways, top 10 = 25%, top 100 = 52%: broad statistics with identifiable
families on top. The opacity extent was also overstated: 17,609 A at
tau > 1e-3, not 23,363 (that was tau > 1e-4).

*E2 prediction, written before the SEDONA 10 km/s runs finished:* resolved RE
N=1 band 0.835 -> ~0.866, Delta_SEDONA +7.8% -> ~+4.5%, expansion side stays
~0.90. First two seeds at 10 km/s: 0.856, 0.853 (resolved), 0.900 (expansion).

## 9q. The ε sweep: the answer is "which band?" (2026-08-22)

*Configuration.* `paper2/phase1/e4_eps_sweep.py` at `de908e3`+ (TLA modes
from `cb07b1f`), full La II from `data/57LaII_levels_calib.txt`
(sha256 870fd713b9b5b822…) and `data/57LaII_transitions_calib.txt`
(4867afa0a3ecf676…), T = 3000 K, n_ion from `forest_lines.npz` (τ_max = 5 in
3850–3950 Å), 949 opacity lines, Planck 6000 K photon launch over
1142–17,697 Å, 2×10⁶ packets × seeds 1–3 per leg, dν/ν = 4.17e-5 for the
expansion legs, ε ∈ {0, .1, .2, .3, .5, .7, .9, 1} on `sobolev_tla` and
`expansion_tla` against `sobolev_branch` and `expansion_branch`; ~8 s per
leg, whole sweep 18 legs × 3 seeds in ~8 min. Output `e4_eps_sweep.json`,
`e4_spectra.npz`, figures via `e4_fig.py`. Band edges printed against strong
lines: 3304 Å (τ 1.2) is 4 Å inside the blue band — noted, not moved.

*What came out* (report §4.20). Monotone F_b(ε) everywhere; the physics
value is outside the TLA's whole range in the optical 4500–6000 Å
(0.9705 ± 0.0005 against a TLA maximum of 0.952 at ε = 0 on the same opacity),
and inside it elsewhere at ε_best that differ by a factor of twenty across
bands (0.016 red, 0.055 blue, 0.36 UV on the Sobolev leg). χ²/dof minima 44
(ε = 0) and 53 (ε = 0.2, expansion leg). I had expected outcome B and half
expected the expansion leg to land near the iron-peak ε ≈ 0.3 in the
3800–3955 Å band; it does (0.32), and that is the trap — the same leg wants
0.68 in the UV and 0.03 in the red, and its band value is the product of an
opacity 21% too transparent (E8) and a closure that over-thermalises.

*E6 said why in one table.* Blue-launched energy: branching sends 8.8% to
the optical and 0.5% to the red; thermalisation sends 7–8% to the optical
and 11–14% to the red. A scalar ε mixes the two kernels; it cannot have the
first channel without the second, and the optical band therefore can never
be refilled to the physics' level. That is the sentence Paper II is about.

*E8.* `expansion_branch` at +21% in the band, +0.2–4% in the wide bands:
line identity through the bin is most of the fix. Bin width: +25% / +21% /
+17% at 1.25 / 12.5 / 125 km/s bins — it does not converge to the physics
with finer bins, it drifts the other way; the floor is the Poisson opacity
error Paper I measured, now under fluorescence.

*Bookkeeping.* The first E4 waiter died on its own timeout (exit 1) a few
minutes before the sweep finished; nothing was lost, the log has the
verdict. E2 came in while I wrote this: SEDONA 10 km/s resolved 0.8546
(3 seeds; predicted ~0.866 — half the predicted move, same direction),
expansion 0.8992 (unchanged, as predicted), Δ_SEDONA +5.2% (predicted ~+4.5%);
MC Sobolev+thermal vs SEDONA resolved now +0.9% in the band, forest sub-band
+2.1%. Gate A passed on its own terms. The 1 km/s run (1 seed, 1 h 47 min
— twice the estimate) came in at 0.8563: +0.69% from the MC, 0.2% from the
10 km/s value, so the resolved side is converged in width by 10 km/s and the
MC's delta-resonance Sobolev leg is the v_D → 0 limit SEDONA approaches, as
Paper I said it would.

## 9r. E13: the verdict survives 0.1c; the frozen shortcut does not (2026-08-24)

*Build.* `relativity="worldline"` in `forest_mc.run_mc` at `e3b5e04`: per-
packet clock, exact Doppler, the addendum's linear locus, moving boundaries,
epoch-diluted tau in both the interaction and the escape probability,
comoving-isotropic re-emission with aberration. Classical path untouched
(46 prior tests bit-for-bit).

*Two bugs the controls caught before any science.* (1) The moving-core
quadratic took the s = 0 root on the launch surface, so the expanding core
never swallowed a packet; fixed, and the caught fraction now equals
beta_core^2 to 4 decimals. (2) The one that mattered: at a just-used
resonance the worldline nu_cm (recomputed from gamma(1 - z/ct), not from the
locus inversion) could round a hair above the line, searchsorted re-selected
it, and the <= 1 cm guard sent the packet TO THE BOUNDARY past the whole
remaining forest -- 16% of interactions gone, transmission up 22% in the
band. The single-line control was blind to it (nothing left to skip); it
surfaced only because the slow-shell worldline branch leg refused to
reproduce 0.659. Fix: a 1e-12 relative searchsorted nudge; the classical
expressions are algebraically aligned and need none. Full-forest
worldline-minus-classical differential now matches `sobolev_attenuation`
to 0.1% absolute (new regression test, 80-point grid -- 40 points
under-resolve the troughs and fail honestly).

*Controls* (report 4.21.1): MC = analytic to <= 0.003 rms in both
conventions at beta_out = 0.01-0.20; the frozen-vs-worldline gap is -0.026
mean / 0.78 max at 0.1c. tau distribution identical across shells by
construction (tau_S carries no shell velocity).

*Result* (4.21.3-5). Fast shell 0.05-0.15c, worldline, 3 seeds x 2e6,
~10-45 s/leg: branch blue 0.043, band3800 0.044, optical 0.355, NIR 1.085.
Frozen-first-order control on the same shell: blue 0.163 -- 3.8x the
worldline value, bigger than any closure-vs-physics difference. eps_best on
the expansion leg: 0.20-0.49 by band, red unreachable from above -> outcome
B holds at 0.1c; but d_eps vs slow = 0.15-0.27, and blue-launched energy
that stayed blue at 2000 km/s (0.785) almost never does at 0.1c (0.015,
with 73% deposited in-shell). So: qualitative closure verdict robust,
calibrated eps not transferable, v_bulk a validity-map axis, and no frozen
high-beta result is usable for closure work.

*Caveats.* eps grids are coarse on the fast shell ({0, .2, .3, .5, 1}
expansion; {0, .065, .3, 1} Sobolev) -- eps_best interpolation carries ~0.05
grid systematic, small against the 0.15-0.27 shifts. Thermal-bin emissivity
weights and the E8 exit kernel stay at t_exp under worldline transport
(stated in the docstring); expansion_branch was not run at high beta.

## 9s. E9/E10: Ce II breaks the closure twice (2026-08-24)

*Setup.* Ce II normalized by setup.py's recipe verbatim (`e9_ceII.ce_n_ion`):
n_ion = 11,641 cm^-3 for tau_max = 5 in 3850-3950 A at 3000 K -- 2,376
window lines against La II's handful, and 22,960 opacity lines overall
(tau > 1e-3), reaching the far-IR. One data quirk: the Ce II GSI files carry
J as '7/2'-style strings, so 2J+1 must go through
`sobolev.populations.statistical_weight` (the first smoke test died on it,
silently, because I had backgrounded the whole launch chain -- lesson
re-learned: smoke in the foreground, launch in the background). The blend
(`ForestAtom.from_gsi_blend`) concatenates ions with level offsets;
branching stays ion-internal by construction and the structural test checks
every upper level's downward lines belong to its own ion.

*E9.* eps_best^La != eps_best^Ce in every band (0.11-0.24 shifts), and the
reachability flips run OPPOSITE ways: La optical unreachable / Ce optical
reachable; La band3800 reachable (0.07) / Ce band3800 unreachable -- direct
branching 0.556 vs 0.531 for pure scattering at eps = 0 on the same
opacity. Outcome C on top of outcome B. chi2/dof minima 108-167, worse than
La's 44-53.

*The one I did not expect: E8's closure is density-limited.* On La II,
expansion_branch had looked like "most of the fix" (+21%). On Ce II it
overfills the band +113% -- because the band is BLACK under Sobolev pure
absorption (1e-4) while the Poisson closure transmits 0.224. Three orders
of magnitude of opacity error from saturation clipping at 24 lines per
angstrom; no exit kernel repairs that. F15's mechanism, now with
fluorescence on top. The Paper II recommendation must carry this limit.

*E10.* The blend answers the plan's three questions cleanly: blanketing
neither suppresses (-44.6% vs La-only -38%) nor relocates the
redistribution error, and eps_best tracks the dominant forest (blend ==
Ce to 0.02 in every band). Composition sets the calibration -- the
outcome-C statement, measured a second way.

*Costs.* Ce legs 5-25 s each (3 seeds x 2e6); whole E9 ~10 min, E10 ~9 min.
Results in e9_ceII.json, e10_blend.json.

## 9t. Paper III R1-R4: the 4x4 matrix that worked, after the artifact that taught me why (2026-08-29)

*Build.* `paper3/` per the received plan: `forest_mc` gains one mode
(`sobolev_group`) and an rng-inert event collector (Gate 0 = 0.0 sigma,
bit-for-bit against e4/e9); `RedistributionKernel` trains on the
chain-collapsed (nu_abs -> nu_exit) events -- 1.18M (La), 4.57M (Ce).
Design choices recorded in paper3/README.md: event-level kernel,
photon-count sampling rows, R^E + q_dep as the conservation object
(exact to 1e-16), q_core left to transport.

*The v1 artifact -- worth a subsection of the methods paper on its own.*
With within-group re-emission from a continuous sub-histogram, the closure
UNDERSHOT the 3800-3955 band by an amount that GREW with refinement (-5% at
8 groups -> -21% at 128) while bolometric stayed sub-percent. Non-
convergence under refinement smells like transport, not information loss,
and the interactions-per-packet counter nailed it: 0.216 -> 0.272 vs 0.196
in branch. A histogram-sampled frequency lands just above the true exit
line half the time, and the packet immediately re-sweeps the line whose
escape probability the kernel's training already resolved. Exits are
discrete line frequencies; sample them exactly (float64 events -- float32
storage would have broken the equality the at-resonance skip relies on)
and the artifact vanishes.

*R3/R4 with discrete tables.* La II: worst band 1.6% at Ng = 4 (!),
chi2/dof 0.2-0.3 across the whole 200-bin spectrum, bolometric -0.01%,
fresh seeds confirm. Ce II: monotone 7.5% -> 2.9% (32) -> 0.7% (64).
Gate 1 excellent (La) / strong-to-excellent (Ce). Why the asymmetry: La's
blue->blue block is 0.79 -- redistribution is nearly input-independent, the
global discrete emission distribution does the work; Ce moves real energy
across bands (blue->blue 0.25) and needs the input dependence resolved.
chi2/dof < 1 partly reflects matched launch seeds with the reference
(correlated noise); the fresh-seed check covers it.

*Answer to the plan's boxed question:* yes -- the middle method is worth
pursuing. Costs: kernels build in seconds; group legs run at branch-leg
speed; tables 14-500 kB dominated by the exit-line list. R5 blocked: no
Nd II GSI files in data/. Next per plan: temperature/epoch transferability
(P5-P6), low-rank structure (P8), the mixture composition rule (P9).

## 9u. P5/P6: three predictions, three confirmations (2026-08-29)

Predicted before the runs, from the kernel's structure (its rows depend on
the state only through {tau_j} and the within-group absorbing-line mix):
source-T transfer near-exact; gas-T genuine (Boltzmann factors move each
line's tau by its own factor); epoch collapsing onto tau_scale (geometry
never enters the kernel; tau propto n t = t^-2 homologously).

All three landed (report 4.24): T_src 4000-8000 K with the fixed 6000 K
kernel -- worst 1.4%, indistinguishable from recomputed; T_gas 5000 K with
the 3000 K kernel -- 9.6%, recomputed 0.2%; the epoch tau-collapse EXACT
within noise at every epoch including t = 0.5 d where tau_max = 34 and the
fixed kernel fails at 13%. The tau-matched control is the elegant one: a
1 d atom with n_ion scaled t^-2 has the identical tau set, its kernel
trains in the 1 d geometry, and it reproduces the 0.5 d transport as well
as the 0.5 d kernel itself does.

So the Phase-7 question answered itself early: R_ij(T_gas, tau_scale,
ion), with T_src free and t folded into tau_scale. A production table is
~4 LTE temperatures x a few tau_scales x ions at tens of kB each.
Costs: each sweep 4 configs x 6-9 legs, whole chain ~13 min.

## 10. Standing environment notes

- Everything SEDONA lives *outside* this repo: code `~/personal/pubsed`,
  deps `~/personal/sedona-deps`. Rebuild: `cd ~/personal/pubsed/src && bash
  install.sh wsl`.
- Figures are working outputs in `outputs/` (gitignored); the copies the
  report references live in `docs/figures/` (committed).
- SEDONA run directories under `experiments/**/run_*/` are gitignored;
  generators + param templates are committed, so every run is reproducible.
- MC noise at `core_n_emit = 2e6` is ~1–2% per band flux; raise it before
  chasing sub-percent effects.
