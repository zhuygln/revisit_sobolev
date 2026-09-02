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

## 9v. Second machine: what a clone does not carry (2026-08-29)

Stood the project up on a new WSL2 box. The science all reproduced --
Paper III Gate 0 bit-for-bit on both ions (La II 1,176,141 events, Ce II
4,571,017, every band 0.0 sigma), R3 at Ng = 4 giving the same 1.62% worst
band, and the Paper II phase-1 driver returning the -38% closure gap. Under
numpy 2.5.2 the reference JSONs differ from the committed ones only in the
last decimal place (summation order) with `n_events` identical, so the event
stream itself is unchanged.

What did *not* survive the move, in rising order of how long it took to see:

1. *Python.* The box's default is 3.9; the package needs >= 3.10 and
   numpy >= 2. Made `.venv` a 3.12 env so the notebook's absolute
   `.venv/bin/python` launch idiom keeps working.
2. *Data.* `data/` is gitignored, so a clone has none of it -- but two tests
   (`test_worldline_forest_differential_matches_analytic`,
   `test_blend_atom_structure_and_identity`) lacked the `pytest.skip` guard
   the other three have and *failed* rather than skipped. Guards added; a
   data-less clone now reports 169 passed / 5 skipped, no red.
3. *The SEDONA build recipe.* `Makefile.wsl` is written in §5 as living "in
   the pubsed clone" -- which is outside this repo, so it was simply gone.
   Rewrote it from §5's prose and this time committed it, at
   `docs/sedona/Makefile.wsl`, with the whole machine walkthrough in
   `docs/sedona/SETUP.md`. It builds clean first try. Both §5 gotchas were
   real and still are: `sedona.h:10` hardcodes `#define MPI_PARALLEL 1`, and
   conda's mpicxx needs `OMPI_CXX=g++` under it. Added `-Wl,-rpath` so the
   env's libs resolve without `LD_LIBRARY_PATH`.
4. *The one that would have bitten silently.* All ten SEDONA drivers
   hardcoded `/home/yozhuz_223/personal/pubsed` -- the first machine's home.
   Every `experiments/` script was unrunnable anywhere else. Now
   `SEDONA_HOME` (default `~/personal/pubsed`) with a `SEDONA_EXE` override.
5. *LaTeX.* conda-forge `texlive-core` cannot build `pdflatex.fmt`: it ships
   `mktexlsr` as a shell script where the TeX Live perl layer wants the
   `mktexlsr.pl` module, so `fmtutil` dies in `@INC`. TinyTeX is the fix.
   The manuscript rebuilds to the same 21 pages, text identical but for
   `\today` -- and the committed PDF is 21 pp, not the 19 the README claimed
   (the README's "68 tests" was equally stale; it is 174).

*A blocked item that was never blocked.* R5/P10 wanted Nd II and the notes
said "no Nd II GSI files in data/". Nd II is in the same Zenodo record as
La and Ce -- it is just the biggest ion in the archive (687 MB of
transitions, 79% of the download), which is presumably why it never came
down. Extracted and parsed: 9,994 levels, 3,336,077 transitions, 10 s and
1.9 GB peak RSS, half-integer J as fractions exactly like Ce II. So Gate 2
(is compression generic across three ions?) and the atomic-data robustness
test are both open.

## 9w. R5: Nd II compresses, and Gate 2 lands on outcome A (2026-08-29)

*The blocked item was never blocked* -- Nd II is in the same Zenodo record as
La and Ce (notebook 9v). Added `nd_n_ion()` to the Phase 0 reference (the
setup.py recipe verbatim, as `ce_n_ion` does for Ce), the `build_atom` branch
and `--ion ndII`; `compression.py` needed only its `choices` list widened,
everything else was already keyed off the ion name. The change is inert for
the existing ions: La II re-ran bit-for-bit afterwards, 1,176,141 events and
every band 0.0 sigma.

*Gate 0 is skipped for Nd, deliberately.* It checks this wrapper against the
Paper II result for the same ion, and there is no Paper II Nd run. Rather
than invent a baseline, the check prints "GATE 0 N/A" with the reason: the
wrapper is already validated -- it reproduced La II and Ce II bit-for-bit --
and Nd changes only which files `build_atom` reads.

*The surprise is the forest size.* Nd II has 3,336,077 transitions, 188x La
and 8x Ce, but at the fixed tau_max = 5 normalization it yields only 4,496
opacity lines against Ce's 22,960. 57,916 Nd lines land in the 3850-3950 A
window (579 per angstrom, against Ce's 24), so pinning the strongest line at
tau = 5 drives n_ion to 1273 cm^-3 -- an order of magnitude below Ce's
11,641 -- and most of the forest falls below the tau > 1e-3 cut. Held to the
plan's own convention, the densest ion in the database is not the deepest
problem: every Nd band sits monotonically between La and Ce.

| ion | bol | blue->blue | worst band error at N_g = 4 / 32 / 64 |
|---|---|---|---|
| La II | 0.9626 | 0.785 | 1.62% / 0.94% / 0.33% |
| Ce II | 0.8256 | 0.247 | 7.48% / 2.89% / 0.74% |
| Nd II | 0.9252 | 0.572 | 0.44% / 0.25% / 0.11% |

*Gate 2: outcome A.* All three ions compress, and two of them at four
groups. Nd is the easiest of the three despite the largest line list.

*What this does to the 9t mechanism.* The blue->blue block was the proposed
explanation for why La compresses at 4 and Ce needs 32-64: high blue->blue
means redistribution is nearly input-independent and the global exit
distribution does the work. Nd tests it and the test is only half passed.
Ce still stands out as the ion with both the lowest block (0.247) and the
only real compression error. But Nd sits at 0.572, below La's 0.785, and
compresses *better* -- so the block does not order La against Nd. The honest
reading is that La and Nd are both sitting on the MC noise floor at every
group count (La's own sequence is non-monotone: 0.33% at 64 but 1.05% at
128, which is noise, not information loss), and this run cannot resolve
them. The block separates the hard ion from the easy ones; it does not rank
the easy ones. Distinguishing them needs more packets, not more groups.

*Costs.* Nd reference 17s for 3 seeds x 2e6 (1,648,874 events); each group
leg 14-15s; atom build 27s including the 687 MB parse, 2.5 GB peak RSS.
Tables are 298-371 kB against La's 14-500 kB -- dominated by the exit-line
list, as before.

*P5/P6 on Nd, same session.* Both sweeps were La-hardcoded in the same three
places (LEV/TR, n_ion from forest_lines.npz, FIXED = kernel_laII_ng32).
Factored `ion_inputs(ion)` out of build_atom and threaded `--ion` through
tsweep.py and epoch.py; La keeps its original output filenames so the
committed tsweep_src.json / epoch.json are not orphaned, and the refactor is
provably identical for La (same paths, n_ion equal to the last bit).

Two of F26's three claims carry to Nd and one does not.

- *tau_scale collapse: holds, and cleanly.* tau_matched == own at every
  epoch (0.91/0.91, 0.25/0.25, 0.15/0.15, 0.14/0.06%), including t = 0.5 d
  where tau_max = 26.4 and the fixed 1 d kernel errs 19.1%. Same structural
  result as La. Geometry still never enters the kernel.
- *T_gas: still the genuine axis*, and stronger than La -- fixed 6.4% (2500),
  11.1% (4000), 19.2% (5000) against La's 3.9/7.6/9.6; recomputed <= 0.51%.
- *T_src: fails.* La's fixed-kernel error is flat at 0.75-1.41% across
  4000-8000 K -- its noise floor, since the recomputed kernel scores the same
  -- which is what F26 rests on. Nd's is +12.44 / +3.89 / +0.25 / -4.88% at
  4000/5000/6000/8000: monotone, and it CHANGES SIGN across the 6000 K
  training point. Noise does not do that. It is a state-transfer signature.

The mechanism is 9u's, read the other way. The rows depend on the radiation
field through the within-group absorbing-line mix. La's 949 opacity lines
mean each group is dominated by a few strong lines whatever the incident
spectrum, so the mix barely moves; Nd's groups draw from a denser, more
evenly weighted set, so reweighting the continuum reweights which lines
absorb and the rows move with it. So (T_gas, tau_scale, ion) is La's state
space, not every ion's -- for Nd, T_src is a fourth axis, and whether it can
be dropped is a per-ion question. Written up as F28, with F26 amended rather
than left to look general.

*One blemish worth repeating at higher N.* At T_src = 5000 the recomputed
("own") kernel errs -2.63% in band3800 where it is <= 0.25% at the other
three temperatures. That puts Nd's band3800 noise floor above La's and makes
the +3.89% fixed point at 5000 marginal on its own. The 4000 K point and the
sign flip are well clear of it, so the verdict stands, but the individual
numbers should not be quoted until the sweep is repeated with more packets.

*Not done here:* P10's atomic-data robustness test (GSI vs independent Nd
data) still needs a second data source.

## 9x. P9: composition is cheap at 5%, and a 5% ion owns a whole band (2026-08-29)

*The rule.* `RedistributionKernel.mix()`: R_mix[i] = sum_s w[i,s] R_s[i] with
w[i,s] species s's share of sum(1 - e^-tau) in group i, taken from the blend's
OPACITY and never from a blend run -- otherwise it is a fit, not a rule.
Convex mixing preserves sum_j R + q_dep = 1 identically, which is a pleasant
free property. Per-ion kernels are trained on the blend's grid (nu_lo/nu_hi
passed explicitly) so the rows are mixable at all.

*Two bugs the smoke test caught, both worth recording.*

1. I weighted each species' exit table by its own ABSOLUTE run size as well
   as by w, counting the composition twice. Symptom: the mixture came out
   worse than its own 95% component, which is impossible for a sane convex
   rule and is what made me look. The per-group absorption budget has to be
   species-independent.
2. My first reading of the smoke test was wrong and the production run
   refuted it. At N_g = 8 / 5e4 the mixed kernel's rowL1 beat ce_only's while
   its transport was worse, and I took that to localize the failure to the
   exit tables. At production the ordering does not hold at all: la_only has
   the LOWEST rowL1 (1.042) and by far the worst transport (45.4%). The
   `live` mask only scores rows both kernels populate, and La leaves most of
   the blend grid empty, so its rowL1 is an average over a small easy subset.
   rowL1 is not a proxy for closure error. Do not use it as one.

*The result (report 4.26).* N_g = 64, La+Ce at reference densities, 94.9% Ce
by opacity: explicit 1.37%, mixed 4.27%, ce_only 10.31%, la_only 45.42%.

The controls were the whole point -- on a 95% Ce blend, "just use Ce" is the
cheap answer unless the rule beats it. It does, by 2.4x, and the gain sits in
one band. ce_only's failure is almost entirely optical, -10.3%; that is the
blue -> optical branching channel from 4.20, and La carries it while being 5%
of the opacity. The rule repairs it to -0.7% knowing nothing but the opacity
weights. A minority species owning a redistribution channel is the thing I
would not have guessed, and it is why the composition rule is worth having.

*The price.* 3.1x the explicit kernel's error at N_g = 64, 2.0x at 32. Clears
Gate 1's "strong" bar (<5% every band) at 64, not at 32; never reaches
"excellent" (<2%). So composition leaves the state space at the 5% level -- a
per-ion library plus opacity weights, no blend training -- and stays in it
below 2%.

*Untested and obvious next:* one blend, one ratio, one temperature. The
interesting case is comparable La:Ce, where neither dominates and the mixing
is not a small correction to one ion.

## 9y. P11: the redistribution half was never the hard half (2026-08-30)

*Setup.* Three opacity treatments carrying the SAME R_ij (N_g = 32, trained
once from sobolev_branch), all scored against sobolev_branch on the same
atom. Two new transport modes in forest_mc: `binned_group` (bin optical
depth sum tau -- the same binning with NO Poisson substitution) and
`expansion_group` (sum(1 - e^-tau), the substitution). `expansion_bins`
gained a weight argument; everything else in the not-sobolev path was
already opacity-agnostic, and the `outcome == "group"` block did not care
which opacity fed it. `binned_absorb` exists as the correctness check:
binned-exact must reproduce line-by-line Sobolev in pure absorption, because
optical depths add (F12), and it does (0.1852 vs 0.1796 in band3800 at 3e5;
expansion gives 0.3451, the F4 error).

*One bug of mine, caught by the counter again.* The `outcome == "group"`
block `continue`s before the end-of-step `tau_r` redraw, which only the
non-Sobolev legs need -- it was written when sobolev_group was the only
group mode. The new legs therefore kept an already-spent optical depth and
re-interacted immediately: 0.944 events/packet against the reference 0.196.
Redrawing tau_r in the group path fixed it (expansion_group 0.186 vs
expansion_branch 0.190). Third time the interactions-per-packet diagnostic
has caught a transport bug (9t, and twice here); it belongs in every leg.

*The result.*

| worst band | La II (949 lines) | Ce II (22,960) |
|---|---|---|
| R_ij alone, exact line opacity | 0.92% | 2.21% |
| + bin resolution (sum tau) | 14.49% | 126.66% |
| + Poisson (sum 1-e^-tau) | 17.84% | 91.29% |
| Poisson + exact A*beta [F24] | 21.32% | 112.86% |

The last row lands on F24's +21% / +113% on the same two forests, so the
harness is validated against a result taken independently.

*What it means.* The redistribution operator is not the problem and never
was. Same R_ij, exact opacity: 0.92% and 2.21%. Group the opacity by either
rule and it is 14-18% / 91-127%. Everything Paper III has established about
redistribution survives untouched; the opacity representation is the binding
constraint and it binds an order of magnitude harder on the dense forest.

*Why both rules fail, and this is the satisfying part.* On La they bracket
the truth from opposite sides -- sum tau too OPAQUE (-14.5%), Poisson too
TRANSPARENT (+17.8%) -- and the events/packet counter says why sum tau
misbehaves: 0.317 vs 0.196, because a tau = 8 line contributes 8 to sum tau
and therefore ~8 interactions where the physics has one. That is F15 read as
a design constraint. Expansion preserves the interaction count E exactly and
gets survival wrong; exact-sum preserves the attenuation S exactly and gets
the count wrong. One scalar per bin cannot carry both, and a scattering
problem needs both. On Ce the bracket collapses and both are too transparent.

*So the target architecture fails on dense forests.* kappa_grouped + R_ij is
not usable for anything Ce-like, and P12 must not be built on it.

*The constructive move, tested the same day, and it half works.* A bin
carrying BOTH quantities -- survival from S = sum tau, the within-bin line
draw from the p = 1-e^-tau distribution (`dual_*` modes, one extra array per
bin). On La II it does exactly what F15 predicts: the saturated band goes
+21.3% -> -0.7% and the worst band 21.32% -> 8.66%. There, F24's density
limit really IS the Poisson survival substitution.

On Ce II it fails, and worse than what it replaces: 139.27% against
expansion's 112.86%. The DIRECTION is what taught me something. S >= E
always, so S survival is strictly more opaque -- yet the band gets BRIGHTER
(1.33 vs 1.18 against a reference of 0.56). More opacity means more
interactions (0.961 ev/pkt vs 0.864), every interaction is a fluorescent
redraw, and on a forest this dense the band refills from elsewhere faster
than it absorbs. That is 4.19-4.20's fluorescent refill working against the
closure. The deep band on a dense forest is REDISTRIBUTION-limited, not
attenuation-limited, so no survival law can fix it. I had assumed the La
result would generalize and said so before Ce ran; it does not.

*And the repair is unavailable to the architecture anyway.* `dual_group` came
out bit-identical to `binned_group` -- every band, both ions, not close but
identical -- because the group path never draws a line within the bin, so the
p-distribution has nothing to attach to. The two-quantity bin only helps a
closure that restores line identity at absorption, which is precisely what
R_ij is defined not to do. So even on La II, where the fix works, it is not
available to kappa_grouped + R_ij.

*Where that leaves P11.* Four opacity rules tried; on the dense forest they
span 91-139% and none is usable. The redistribution operator is 0.92% / 2.21%
throughout. The problem is not which scalar the bin carries, it is that a bin
has no way to skip the line a packet was just emitted from, and on a dense
forest fluorescent refill amplifies that rather than damping it. Restoring
line identity at emission is the only lever left, and it costs exactly the
thing grouping was supposed to buy.

## 9z. The minimal-memory test: right mechanism, wrong scale (2026-08-30)

*The ask.* F30 said the grouped opacity cannot skip the line a packet was
just emitted from. Minimal test of that diagnosis: carry ONE number per
packet, the frequency last emitted at, and credit that line's own tau to the
next free-path draw. `run_mc(line_memory=True)`. No atomic level inspected
after emission -- the credit is an opacity lookup by frequency, and exits
below tau_min carry no opacity and cost nothing. Six lines of transport.

*Result.* La II 14.49% -> 6.50% (saturated band -14.5% -> -6.3%). One float
per packet, factor 2.2, and it beats the two-quantity bin's group variant
which costs an extra array per bin. Ce II 126.66% -> 116.31%, and
91.29% -> 79.46% on the Poisson opacity. Real, right direction, no rescue.

*So the last emitted line is NOT the minimal missing state variable*, and the
events/packet counter says why in a way I found genuinely clarifying. Memory
removes a COMPARABLE FRACTION of the excess interactions in both forests --
La 0.317 -> 0.270 against a reference 0.196, so 39% of the excess; Ce
0.918 -> 0.865 against 0.762, so 34%. The mechanism is doing the same job in
both. But a third of the excess interactions nearly halves the La error and
barely dents the Ce one. Therefore the Ce error is not driven by excess
interactions at all. It is 4.19-4.20's fluorescent refill: at 24 lines per
angstrom the deep band fills from elsewhere faster than it absorbs, and no
correction to LOCAL interaction bookkeeping can reach that.

*The clean statement of the two regimes.* Sparse forest: a packet's fate is
set by the few resonances it meets, so one remembered line recovers most of
what binning destroyed. Dense forest: the band is refilled by redistribution
from outside, and the missing information is not a scalar correction but the
RESONANCE SEQUENCE -- which lines, in what order, with how much bin between
them. A single number per bin cannot encode that and one remembered
frequency does not restore it.

*A confirming detail.* Memory helps the exact-sum opacity and NOT the Poisson
one (La 14.49 -> 6.50 vs 17.84 -> 20.03). Exactly the expected signature:
with S survival the emitting line carries its full tau into the bin so
re-absorbing on it is severe; with E every line contributes at most 1, so
self-absorption was never dominant there. The mechanism behaves as diagnosed
even where the cure is not enough.

*What this sets up.* If one bin scalar fails and one remembered line does not
rescue it, the next thing to try is not another scalar patch -- it is keeping
a compact ORDERED list of resonances per group, compressing their properties
rather than their count, since the count is what carries the ordering.

## 9aa. Stage 0, and a correction I predicted wrong (2026-08-30)

*The units bug was real.* `line_memory` credits the emitting line's optical
depth back to the packet, but `tau_r` is measured in the GRID's units -- the
weight `expansion_bins` was built from. The first implementation credited
`op_tau` everywhere, over-crediting a saturated line on the Poisson grid by
tau/(1-e^-tau), a factor 8 at tau = 8. Now `op_p` on the Poisson grid,
`op_tau` on the exact-sum and dual grids. (I got the dual case wrong on the
first pass too -- "dual" carries SURVIVAL on the exact-sum grid and only its
line-selection on p, so its credit is tau. Caught before running.)

*Effect: small, and confined to the expansion legs.* La 20.03 -> 19.84%,
Ce 79.46 -> 81.66%. Binned and dual unchanged, as they must be. F31's table
in report 4.28 and the README row are corrected.

*What I predicted and got wrong.* I expected the fix to remove the anomaly
that memory HELPS the exact-sum opacity and HURTS the Poisson one. It does
not: La is still 14.49 -> 6.50% with S and 17.84 -> 19.84% with E. The real
explanation is the sign of the failure being corrected, not a units artefact:

- memory always makes a grouped opacity MORE TRANSPARENT (it credits away
  optical depth, so the packet travels further before interacting);
- La's Poisson leg is already too transparent (+17.8%, the survival
  substitution), so more transparency is worse;
- Ce's error is the opposite kind -- its band is OVER-FILLED by fluorescent
  refill driven by excess interactions -- so removing excess interactions
  removes refill and the error falls, 91.29 -> 81.66%.

Memory's sign is therefore set by which failure mode dominates, not by which
opacity rule is in use. That is a sharper statement than the one 9z made, and
it only became visible once the units were right.

## 9ab. Stage 0: the phase-diagram infrastructure (2026-08-30)

*`sobolev/forest_stats.py`.* Five statistics existed only as copy-pasted
one-liners at five call sites and two did not exist at all. Now one module:
`saturation_stats`, `SE_sums` (F15's S, E and E/S), `crowding` (saturated
lines per unit ln lambda -- scale-free, so a forest at 4000 A and one at
12000 A with the same crowding land on the same point), `spacing_stats`
(wrapping the existing `nearest_neighbour_velocity_spacing`), and the new
`redistribution_range`, which is the one that matters: mean/median
|d ln lambda| per event plus the SAME-GROUP FRACTION. Ray-resolved S and E
already exist as `crossing_depths`; deliberately not reimplemented.

*`paper3/synthetic/forest.py`.* No new class needed -- `ForestAtom.__init__`
already builds `branch_lines`/`branch_cum`/`exit_cum` from `upper` and `A`
alone. Exit channels get f_osc = 0 and n_lower = 0, so tau = 0 keeps them out
of the opacity set while `beta_all` forces them to escape freely: the
`three_level` pattern, generalized to N lines.

*The dial is faithful, and there are two of them.* Measured against the
branch reference at 150 lines, tau = 5, span 0.25:

| dial dlnlam | measured &lt;\|dlnlam\|&gt; | same-group |
|---|---|---|
| n_exit = 1 | 0.0000 | 1.000 |
| 0.02 | 0.0167 | 0.166 |
| 0.15 | 0.1251 | 0.166 |
| 0.40 | 0.3339 | 0.165 |

n_exit = 1 is EXACTLY coherent scattering (`nu_out == nu_in` bit-for-bit, not
approximately) -- the zero of the range axis. The measured range is ~0.83x the
dial, which is not a leak but the A*beta weighting: the absorbing line has
beta = 0.199 at tau = 5, so 81% of exits take the free channel. And the
same-group fraction is flat in dlnlam while `f_return` moves it, so HOW FAR
exits reach and HOW OFTEN a photon returns to its own resonance are two
independent knobs. That is exactly the orthogonality the phase diagram needs.

Saturation (E/S) and crowding are invariant to the redistribution dial to
1e-12, pinned by test. 23 new tests.

## 9ac. E1: the answer is no, and the negative is better than the positive would have been (2026-08-30)

*The question.* P8, reframed: does the redistribution operator have only a few
macroscopic modes? If so that would EXPLAIN F25/F27's compressibility instead
of just recording it. Prediction going in was "La near rank 1-2, Ce several,
effective dimension <~ 3". All wrong.

*1. The dimension never saturates.* Photon-operator participation ratio grows
as N_g^0.64 (La), ^0.66 (Ce), ^0.75 (Nd), with PR/N_g falling 0.5 -> 0.2 and no
plateau anywhere from N_g = 4 to 128. Three ions with 949, 22,960 and 4,496
opacity lines give the same exponent to ~15%. There is no intrinsic mode count;
the operator looks the same at every resolution you examine it with.

*2. Low-rank approximation is simply bad.* NMF rank 8 of a 25-row operator
still misses row-L1 0.47 -- a 23% total-variation error on every row, against a
maximum possible 2.0 -- and the transport error tracks the reconstruction
error, so this is not transport being hypersensitive.

*3. Rank ANTI-correlates with compressibility.* Ce has the lowest
energy-operator dimension of the three (PR 1.60, sigma1 = 77%) and is the
hardest to compress. La has the highest (6.58) and compresses at four groups.
And the two operators disagree by a factor of seven on the SAME kernel (Ce: 1.60
by energy, 11.08 by photon count) because energy piles into one destination
while photons scatter everywhere. "The rank of the kernel" is not well posed.

*The comparison that settles it,* at matched parameter count on La II:

  coarsening to N_g = 4          16 params   1.62%
  rank-4 truncation at N_g = 32  228 params  76.11%
  rank-16 truncation at N_g = 32 912 params  11.30%

Sixteen numbers beat nine hundred and twelve by a factor of seven. Coarsening
AVERAGES neighbouring groups; truncation PROJECTS onto modes. Only the first
works, because the operator is LOCAL IN FREQUENCY and not LOW-RANK. A kernel
that varies smoothly with input frequency has high numerical rank on a fine
grid and coarse-grains perfectly -- smoothness is what F25/F27 were measuring
all along, and I had been calling it "few modes" without checking.

*Why the negative result is worth more.* It unifies the two halves of Paper
III into one statement. Redistribution is smooth at the group scale, so it
coarse-grains (F25, F27, F29). The opacity is a comb of resonances whose
ORDERING inside a bin decides a packet's fate, so it does not (F30, F31). Not
two unrelated facts about compressibility -- one fact about what is smooth at
the group scale.

*A method mistake, recorded.* I first truncated with SVD, which is ill-posed
for a stochastic matrix: it produces negatives, and clipping them to zero then
renormalizing destroys the distribution. It gave non-monotone transport errors
of several hundred per cent, which I nearly wrote up before noticing that a
709% error at rank 3 is a statement about my truncation and not about Ce. NMF
is well-posed and is also the right physical model (each input group as a
mixture of k archetypal exit distributions). The conclusion did not change --
which is luck, not vindication.

## 9ad. E2: the sequence is not the answer, and Nd breaks the density story (2026-08-30)

*Machinery.* `line_memory` is now a DEPTH, not a switch: a ring buffer of the
last m opacity-line INDICES per packet (indices not frequencies -- half the
memory, and no float-equality test at credit time), crediting each remembered
line still ahead. The ordering constraint matters: comoving frequency falls
monotonically along a leg, so only lines at or below the emission frequency can
still be reached, and crediting the ones already swept past would be wrong.
The credit draws no RNG, so every m shares one stream and m = 0 is
BIT-IDENTICAL to no memory -- pinned by test, exactly rather than
statistically, which is what makes the sweep a controlled comparison.

*Result 1: memory saturates immediately, and does nothing for Ce.*

| m | La binned | Nd binned | Ce binned |
|---|---|---|---|
| 1 | 6.50% | 8.41% | 116.31% |
| 4 | 6.39% | 7.98% | 116.86% |
| 16 | 6.40% | 7.95% | 116.79% |

La converges by m = 4, Nd by m = 8, each gaining a few tenths beyond m = 1.
Ce gains NOTHING across a sixteen-fold increase in remembered history.

So 9z's closing line -- "dense forests need the resonance sequence" -- is
WRONG and I am retracting it. It was a plausible reading of F31's two regimes
and the depth sweep kills it. Memory is a BETWEEN-STEP correction; the
dense-forest failure is a WITHIN-STEP one, about where in a crowded bin the
absorption happened. No amount of past history supplies that. I should have
been more careful calling an interpretation a conclusion.

*Result 2, which I did not expect at all: Nd's expansion closure WORKS.*
Running Nd through the full P11 decomposition for the first time:

| ion | opacity lines | sat. lines in 3800-3955 | sum tau in band | expansion + A*beta |
|---|---|---|---|---|
| Nd II | 4,496 | 1 | 11.4 | 1.79% |
| La II | 949 | 4 | 22.4 | 21.32% |
| Ce II | 22,960 | 24 | 89.6 | 112.86% |

Nd has 4.7x La's opacity lines and one twelfth its error. 1.79% with the exact
exit kernel, 2.15% as a full group closure -- a working closure on the ion with
the biggest line list of the three. So total line count does NOT order the
failure and F24's "density limit" is misnamed.

*What does order it, on three points:* band-local saturation. Saturated lines
inside the band that fails (1, 4, 24), or equivalently per transport bin
(0.001, 0.004, 0.025), or sum tau in the band (11.4, 22.4, 89.6). Forest-
averaged crowding also orders them (4, 15, 26 per unit ln lambda) but the
band-local version is the right scale -- the failure lives in one band while
the forest spans a decade and a half, and averaging over the whole span was
hiding the very thing that separates Nd from Ce.

Three points cannot fix an exponent; the successive log-log slopes are 1.8 and
0.9, which is not a law. But the phase diagram now has a TARGET ORDERING to
reproduce and a candidate axis that can be varied independently in synthetic
forests. That is a much better position than E3 would have been in yesterday.

## 9ae. E3: saturation controls it, redistribution does not, and the collapse is partial (2026-08-30)

*The calibration step mattered more than I expected.* Nobody had measured the
tau distribution INSIDE the failing band. Doing it: La 15 lines median tau 0.53
frac(tau>1) 0.267; Nd 163 lines median 0.005 frac 0.006; Ce 462 lines median
0.019 frac 0.052 -- and the ln-spread is 1.83 / 1.71 / 2.05, i.e. essentially
the same across ions differing 30x in line count. Real forests are mostly WEAK
lines with a saturated tail. My first synthetic grid used spread 0.4, near
monodisperse, which makes a forest either transparent or entirely black instead
of transmitting through the weak-line population -- every condition degenerate,
redistribution-only control at 155%. With spread 1.8 the control returns to
0.6-10.6% and the sweep means something. Worth remembering: calibrate the
synthetic against the real DISTRIBUTION, not just its summary statistics.

*Result 1, and the one I would keep if I could keep only one.* Across 96
conditions the Spearman correlations are +0.91 (sum tau in band), +0.86 (N_sat
in band) -- and +0.25 and -0.31 for the two REDISTRIBUTION axes. The
redistribution range and the same-group fraction barely matter. Coming from a
different direction, that is E2 again: memory depth did nothing, redistribution
range does nothing, so the grouped-closure failure is set by OPACITY structure
and not by redistribution structure at all. Two independent experiments now say
the same thing, which is more than either says alone.

*Result 2, the collapse.* dF = 0.162 N_sat^0.58, scatter x1.95, R2 = 0.64
(sum tau gives R2 = 0.77 with slightly worse scatter). Real atoms against that
fit: Nd 0.71, La 0.40, Ce 1.25 -- two of three inside the family's own scatter,
La the outlier at 2.5x, real slope 0.78 vs synthetic 0.58.

*A unit error of mine, caught.* I first compared the synthetic fit (a FRACTION)
against real errors quoted in PERCENT and got ratios of 13-125, and briefly
believed the real atoms missed the curve entirely. They do not. Check units
before believing a two-orders-of-magnitude disagreement -- especially one that
would have killed the PRL route on the spot.

*Two explanations tested, one eliminated.* The synthetic band is 0.19 of a
0.3-wide forest where the real band is 0.04 of ~3.4, so photons cannot enter
from as far away. Varying the span at fixed band crowding moves the error
83% -> 69% across a TENFOLD change in that ratio: geometry is not it. The exit
rule is: at matched N_sat = 13 the dialled forests give 78-94% and the LADDER
forests, whose branching emerges from real cascade structure, give 58.7-62.2%
against a real-atom interpolation of 62.1%. So the dial is qualitatively right
and quantitatively harsh, and the ladder is the better model -- which is exactly
what the "both" option was built to find out.

*Where this leaves the gate.* Partial. The controlling variable is identified
and the null axes are identified, but factor-2 scatter, one 2.5x outlier and
three real points are not "a general phase diagram and scaling law". The
cheapest strengthening is more real atoms: the GSI archive has 27 ions, all
already on disk, and three have been used.

## 9af. E3b: thirteen ions, a broken normalization, and a sign change (2026-08-30)

*The normalization had to be fixed before anything ran.* setup.py's recipe pins
tau_max = 5 INSIDE 3850-3950 A. That works for La/Ce/Nd because their strongest
lines sit in or near that window -- an accident. For Yb II the window holds only
weak lines, so reaching tau = 5 there wants n_ion = 1.7e12 cm^-3, which puts
tau_max at 1.7e8 and beta at 1.5e-8; a packet entering that resonance never
escapes and `run_mc` raised "re-absorption chain did not terminate". That crash
is a real physical statement about the recipe, not a bug. Switched the survey to
the ion's GLOBAL strongest line: scale-free, identical for every ion, beta >=
0.199 by construction.

Worth flagging beyond this experiment: that window recipe is load-bearing in
F24, F27, F30-F33. It does not invalidate them -- they are internally consistent
and all three ions are ones where the window happens to hold near-strongest
lines -- but no claim of UNIVERSALITY can rest on a normalization that silently
depends on which ion you picked. Ce II's band3800 binned error is +126.7% under
the window recipe and +12.2% under the global one, from a 25% density
difference. That should have been checked much earlier.

*Thirteen ions do not collapse.* La II S = 13.4 -> 6.55%; Pr II S = 13.8 ->
31.39%. Identical saturation, 5x the error. Ce II at S = 66.8 errs 12.20%, LESS
than Pr II at a fifth its saturation. The zero-saturation ions (Pr III, Tm III,
Yb II, Ce III, Er III) anchor the trivial limit at 0.01-1.25%, and that is where
most of the +0.75 rank correlation comes from.

*And on real ions the axes cannot be separated.* Saturation and redistribution
range correlate with the error equally (+0.75, +0.77) because in real atoms they
are correlated with EACH OTHER -- a bigger ion has more of both. This is exactly
the confounding the synthetic forests were built to break, and there, with the
axes decorrelated by construction, saturation wins 0.91 to 0.25. The real survey
cannot replace the synthetic experiment; it can only test transfer. Good to have
learned that the cheap-looking route (more real ions) does not substitute for
the expensive one.

*The thing that explains all of it: the error CHANGES SIGN.* Density scan on
Ce II, band 3800-3955:

  n_ion   2910  5821  8731  11641  17462  29103
  binned -33.4 -15.3 +21.4 +124.6  +94.7 +129.8 %

Too OPAQUE at low density, too TRANSPARENT at high, crossing zero between
S ~ 45 and 67 and then climbing +21 -> +125% across a density factor of 1.33.
La II stays negative over its whole range (+0.5% to -36.6%) and never crosses.

That is why nothing collapsed. A sign-changing error cannot be a power law, and
two ions at matched saturation can differ 5x simply by sitting on opposite sides
of the boundary. The two competing errors are ones already named -- exact-sum
binning over-counts interactions (too opaque), fluorescent refill and the
Poisson substitution add transparency -- and the locus where they cancel is a
genuine phase boundary. F34's power-law framing is superseded by its own
follow-up, one day later.

## 9ag. E4: the model has only one side of the boundary (2026-08-31)

*The plan.* Scan tau at fixed (n_lines, dlnlam) to BRACKET the Delta F = 0
crossing rather than grid blindly, then ask whether the crossing S moves with
redistribution range -- the question that decides 1-D vs 2-D.

*The answer: there is no crossing in the model.* All 36 conditions negative,
S from 2.7 to 1459, deepening monotonically to -99% and never turning.
dlnlam changes nothing structural (-4.6 -> -76.3% at 0.005 vs -5.4 -> -82.7%
at 0.15, same shape).

*Why, and the diagnostic that found it.* Reference transmission in the measured
band, at matched saturation:

  synthetic S = 54.1  ->  ref 0.229   dF -74.7%
  synthetic S = 87.5  ->  ref 0.251   dF -90.8%
  real Ce II S = 66.8 ->  ref 0.581   dF +12.2%
  real Tm II S =  5.8 ->  ref 1.049   dF  -6.0%

Tm II's band transmits MORE THAN THE CONTINUUM ENTERING IT. Dy III 0.970. Real
lanthanide forests FEED the measured band from outside; that is what keeps it
transparent as saturation rises, and eventually what makes a grouped closure
over-bright.

My synthetic forests cannot do that. Exit channels sit at a fixed +-dlnlam from
their OWN absorbing line and carry no opacity, so the model redistributes
locally and never delivers a net inflow to the band. The band only darkens. The
too-bright branch -- the entire reason the boundary exists -- is absent by
construction.

*So E4 is blocked on the model, not the measurement.* You cannot locate a
boundary in a model containing one side of it. I should have checked
`ref_core` against the real ions' band transmission when I calibrated the
forests in 9ae; I matched the tau DISTRIBUTION and never checked what the band
actually transmitted, which is the quantity the whole experiment is a ratio to.
Match the observable, not only the inputs.

*What the negative is worth.* In real atoms refill and saturation cannot be
separated -- bigger ions have more of both -- so the real-ion data implies
refill matters but cannot prove it. A model that lacks ONLY refill fails to
produce ONLY the bright branch. That is as close to a controlled demonstration
that fluorescent refill causes the sign change as this project has managed, and
it came from the model being wrong rather than right.

*The fix, for whoever runs it next:* exits distributed over the forest instead
of offset from their own line; enough channels per upper level that a photon
absorbed outside the band can reach it; and an acceptance test that the
reference band transmission approaches and sometimes exceeds unity, as Tm II
and Dy III do. Until then synthetic forests calibrate the too-opaque regime and
nothing else.

## 9ah. The fix: necessary, not sufficient, and the crossing runs backwards (2026-08-31)

*The change.* `synthetic_forest(delocalize=p)` puts a fraction p of exit
channels anywhere in the forest instead of at +-dlnlam from their own parent
line. That is the physical content of a real upper level decaying to lower
levels spread across the whole term structure, and it is what 9ag identified as
missing.

*It works, partially.* At S = 54.1, p = 0 -> 1 with six exit channels:
band transmission 0.226 -> 0.292, binned error -77.5% -> -25.3%. Factor of
three. But Ce II sits at 0.581, so the model is still about twice too opaque in
the band.

*And the crossing it produces runs BACKWARDS.* Scanning tau at p = 1,
n_exit = 6, the expansion leg goes

  S      8.1   20.3   54.1  135.1  337.8  810.8
  exp   +4.0   +1.5   -7.9  -20.8  -36.9  -47.3 %

Too bright at LOW saturation, too opaque at HIGH. Every real ion does the
reverse -- Ce II runs -33.4% at S = 22.3 to +124.6% at S = 89.6. So
delocalizing gives the model a boundary but not THE boundary.

*Why, most likely.* My exit channels terminate on sink levels with zero
population, so a photon leaving through one can never be re-absorbed and
cascaded again. Real forests refill the band MORE as saturation rises, because
more absorption elsewhere feeds more re-emission and the cascade continues;
mine refills LESS, because the exit is terminal. Giving exit lines their own
opacity is the next change and it is small.

*A bug this exposed, which matters more than the physics.* `crossing()` tested
only `a < 0 <= b`. It detected negative-to-positive crossings and silently
missed the opposite -- exactly the case the fixed model produces. It printed
"no crossing" for a row that plainly had one, and I nearly reported that. Now
detects either direction and says which. The earlier all-negative sweep is
unaffected, but the check had been one-sided since I wrote it, and a one-sided
test for a SIGN CHANGE is a bad piece of code to have written in an experiment
whose entire subject is the sign.

## 9ai. The one iteration: yes, and item 3 for free (2026-08-31)

*The binary question was: does adding recurrent opacity to the exit channels
recover the correct orientation of the phase boundary? Answer: yes.*

`exit_tau` gives exit lines their own optical depth on a small set of SHARED
POPULATED lower levels, as a real atom has. Absorption on an exit line then
returns the photon to that upper level and it branches again. Scanning tau at
delocalize = 1, n_exit = 6:

  exit_tau = 0    binned all negative, no crossing;
                  expansion crosses pos->neg at S = 23.7 (backwards)
  exit_tau = 0.5  binned -62.3 -> +2.5 -> +39.1%, crossing NEG->POS at S = 179
  exit_tau = 2.0  binned -53.7 -> +6.3 -> -2.7%, crossing NEG->POS at S = 688

So the boundary now exists in three logically independent settings: the Ce
density scan, the thirteen-ion survey, and a controlled forest where saturation
and redistribution range are separated -- which real atoms never allow. And it
required exactly the piece of physics 9ag predicted was missing, which is the
strongest form of confirmation a diagnosis gets.

*Stopping here, as agreed.* The crossing LOCATION is not universal (179 vs 688
as exit_tau goes 0.5 -> 2.0), the expansion leg still fails, and band
transmission (0.13-0.18) is still well below the real 0.58. Tuning further
would be fitting a toy to Ce II and would read as post hoc. Existence and
orientation are what the synthetic experiment was for; the real ions carry the
rest.

## 9aj. Item 3 came free, and it is better than expected

The three counterfactual legs already existed as modes -- sobolev_group (exact
opacity, grouped redistribution), expansion_branch (grouped opacity, exact
A*beta), expansion_group (both) -- so the audit data answered item 3 with no
new runs. Signed band 3800-3955:

  ion       S      A       B       C     C-(A+B)
  Nd II    8.6   +1.1%   -1.6%   -1.4%    -0.9%
  La II   13.4   -0.4%  +15.8%  +14.7%    -0.8%
  Pr II   13.8   +0.1%   +4.3%   -2.0%    -6.4%
  Ce II   66.8   +1.4%   +4.8%   +1.9%    -4.3%

|A| <= 1.4% everywhere: the redistribution approximation contributes almost
nothing, which is F30 reached a third way and this time MEASURED rather than
inferred by difference. B dominates. But at S >~ 14 the errors stop being
additive, and for Pr II the interaction term (-6.4%) is big enough to FLIP THE
SIGN: opacity alone leaves the band +4.3% too bright, both together -2.0% too
opaque.

That is the causal content of the boundary. The zero is not where the dominant
single-approximation error changes sign; it is where B plus the interaction
does. A closure whose two halves have each been validated separately can fail,
or appear to succeed, for reasons neither half shows alone -- the same lesson
as F36's memory correction taking Ce to +0.2% while moving La THROUGH zero, now
with a mechanism instead of a coincidence.

## 9ak. Item 5: the kilonova crosses at 1.2 days (2026-08-31)

*Setup.* Homologous ejecta, rho ~ t^-3 from 2e-17 g/cm3 at 1 d, X_lan = 0.1,
T ~ t^-1/2 from 5000 K, geometry expanding with t, 0.5-8 d. ASTROPHYSICAL
normalization (from_conditions), deliberately not the controlled standard --
different question. n_ion sweeps four orders of magnitude, so the trajectory
must sweep S.

*Ce II crosses at t = 1.17 d, S = 47.5.* +64.6% too bright at 0.5 d, zero at
1.2 d, -28.4% too opaque at 1.5 d. Ninety points across a factor of three in
time, straddling exactly the epoch kilonova spectra get taken.

The number that matters: S = 47.5 at the crossing, against a boundary located
independently at S ~ 50 by the density scan (44.8-67.2), the thirteen-ion
survey (13.8-66.8) and the synthetic forest. Four constructions, one value.

*Stated carefully:* in this trajectory the zero is the OPACITY error changing
sign (B crosses at 1.21 d, essentially with C), not two large errors
cancelling -- |A| <= 2.1% throughout. The cancellation mechanism shows up
separately at 2 d, where the binned closure reads -1.1% while its own opacity
piece reads -4.1%. Two different mechanisms with the same practical
consequence; the paper should not conflate them.

*La II, same history, does not cross -- and gives a better demonstration.* At
t = 0.75 d the EXPANSION closure reads +0.1% and the BINNED closure reads
-55.7%. Same epoch, same ejecta, same atom. The two differ only in whether a
bin carries sum(1-e^-tau) or sum tau, both defensible, one apparently exact and
one wrong by more than half. Anyone validating the first at this epoch learns
nothing about the second.

Ce shows a closure sweeping through zero as conditions evolve; La shows two
closures at ONE epoch 56 points apart with one sitting on zero. Neither
near-zero residual carries information about correctness. That is the thesis,
and La states it without needing the ejecta to evolve at all.

## 9al. §10: the error an observer sees is a colour, not a luminosity (2026-09-01)

*The gap.* Every result from 4.23 to 4.36 is a band ratio in a 155 A window.
Nobody observes that. A grep for filter / bandpass / magnitude / colour / light
curve over the whole repo returned nothing, and `trajectory.py` calls
`band_ratio` and throws away the emergent spectrum `run_mc` already hands it.
So §10 is mostly plumbing that should have existed from the start:
`sobolev/photometry.py` (absolute L_nu from the packet list, L_bol in the launch
window, AB magnitudes in seven top-hat bands, colours) and
`paper3/phase11_observables/observables.py` (F40's ejecta history, photometry
instead of band ratios).

*The normalization is checked twice* -- `planck_luminosity` against 4 pi r^2
sigma T^4 to 1e-4, and the MC spectrum divided by the analytic core continuum
against `band_ratio` to 0.2%. Getting an absolute scale wrong by 4 pi would have
been invisible in every plot and fatal in every number.

*Three changes from trajectory.py, each for a reason.* A fixed 1000-30000 A
launch window (F40 took it from each atom's own opacity extent: 1128 A to 36.8
um for Ce, so packets pile up in a far-IR tail no filter samples and no two
epochs compare). A cooling core, with the frozen 6000 K core kept as the control
that separates source reddening from opacity reddening. And `crossing_epoch`,
because F40's headline -- 1.17 d, S = 47.5 -- was interpolated by hand and
stored nowhere.

*Result 1: the redistribution approximation is invisible.* Worst |dm| over every
band and epoch: Ce 0.006, La 0.008, blend 0.008 mag, and 0.021 even at 0.2c.
F30/F38 said "the opacity
binds, not the redistribution" in band ratios; in magnitudes it says the kernel
compression of F25/F27 costs less than a photometric error bar. That is the
strongest form the compressible half of the hierarchy has taken.

*Result 2: the error is chromatic.* Ce II at 0.5 d is 0.14 mag too bright
bolometrically and 0.74 mag wrong in g-r -- five times larger -- because the
closure moves flux from r into g rather than creating or destroying it. The
four-ion blend does the same (0.737), so it is not a single-ion artefact.

*Result 3: and the band residual does not track it, in either direction.*
La II binned at 0.75 d: -59.5% in band, +0.007 mag bolometric. Ce II at 2.0 d:
-26.0% in band, every |dm| <= 0.006. Ce II at 0.5 d: +55% in band, 0.455 mag in
g. The band is a good diagnostic of the mechanism -- 4.32-4.36 were right to use
it -- and a bad proxy for an observable.

*Result 4: two defensible closures, opposite colours, one epoch.* La II at 1 d,
expansion gives d(g-r) = -0.057 and binned +0.096. Same ejecta, same atom,
differing only in whether a bin carries sum(1-e^-tau) or sum tau.

*Result 5: it is opacity, not the source.* Freezing the core at 6000 K instead
of cooling it changes the worst colour error by under 0.01 mag on both ions.

*Result 6: Paper I's geometry is a floor.* trajectory.py inherits a 1000-3000
km/s shell. tau does not care (tau ~ n t is velocity-free in homologous flow),
but the wavelength interval a packet sweeps is dv/c, and under worldline
transport the outer boundary recedes so packets stay in longer -- 2.07 vs 8.63
events per packet at 0.1c. `velocity.py` holds epoch, density and composition
fixed and moves only the shell velocity: |dm|max grows 34x (La) and 67x (Ce)
from 0.01c to 0.3c while dm_bol stays under 0.011 mag. Every magnitude in
Results 1-5 is a lower bound.

*Result 7, and the thing I nearly got wrong.* At v_out = 0.01c the worldline and
time-frozen treatments give reference band-3800 fluxes of 0.295 and 1.079 -- a
factor of 3.7, at a velocity where relativistic corrections are 1%. It is not a
bug: the frozen shell lets packets escape that the expanding one keeps. Rerunning
the whole history under worldline moves La II's 0.5 d band residual from +20.9%
to +89.3% -- and moves dm_bol from -0.032 to -0.055 and d(g-r) from -0.110 to
-0.185, the same signs, the same structure, the same 3.4x colour-to-bolometric
ratio. So: the band residual is fragile to the transport treatment and the
magnitudes are not. Had I not checked, I would have quoted velocity-scaled
magnitudes off a transport path I had not validated.

*The density arithmetic that should have been in 4.36.* rho(1 d) = 2e-17 is
tuned so Ce crosses inside the window; at 0.2c it corresponds to 5.9e-6 Msun,
which is not ejecta. But tau goes as n_ion, not rho, and in n_ion the tuned
value is within 5-20x of a LANTHANIDE-POOR component (M = 0.01 Msun, X_lan =
1e-3) and 5100x below a lanthanide-rich one. So the boundary lies where a blue
component lives, and red ejecta sit far past it for the whole observable window.
Whether a kilonova crosses is set by X_lan and M_ej -- the parameters people
infer from the spectra.

*Result 8: a physically normalized kilonova.* M_ej = 0.01 Msun, v = 0.05-0.2c,
X_lan = 1e-3, rho derived from the mass instead of tuned, worldline, 2-12 d.
The practical closure is 0.06 mag bolometric and 0.63 mag in g; the binned
variant reaches d(g-r) = -0.769 at 2 d and d(r-i) = +0.742 at 3 d, where its
bolometric error is +0.000. An order of magnitude between what a bolometric
check would see and what a colour would.

And the diagnostic band dies exactly where it matters: from 4 d on the reference
3800 A flux is below the level at which a ratio means anything, while dm is
still tenths of a magnitude. Anyone validating this closure on band residuals
runs out of signal before the interesting epochs.

*Noise discipline, and a mistake caught by it.* The `sobolev_group` leg doubles
as an empirical noise floor: it is <= 0.008 mag wherever the reference is well
sampled, so a run in which it reads 0.31 mag is reporting its own resolution.
The first blue-kilonova attempt did exactly that, and cranking the packet count
2.5x did not fix it -- because the cause was not sampling. The ejecta cool, and
by 12 d the g band carries 2.5e-6 of the reference bolometric luminosity; its
"0.75 mag error" was a handful of photons. The fix is a mask, not more packets:
`summary.py` drops any band carrying under 1% of the reference L_bol, computed
from the stored spectrum so nothing needs re-running. With it the blue run's
floor falls to 0.098 mag and the surviving entries are 6-8x above it.


## 9am. Paying the provenance debts, and a real filter curve (2026-09-02)

*Why first.* The action plan's Phase 1 starts with a source model and a grid,
but the audit (`findings_audit.md`) listed three findings whose numbers existed
only in this notebook: the Ce II density scan that first located the boundary,
the F38 A/B/C table, and the F39 exit_tau scan. Nothing new should be built on
top of numbers that cannot be regenerated, so these came first. Each took under
an hour of wall time and each turned up something.

*The density scan, three seeds.* `density_scan.py` reproduces §4.32's six
points within 7 points at five of them and within 9 at the third (+21.4 ->
+12.2%). That third point is where the curve climbs 115 points over a density
factor of 1.33, so a seed moves it a lot; the crossing itself is stable at
S = 55.7 (binned) and 64.6 (expansion), both inside the bracket F35 quotes.
Lesson filed: quote a crossing, not the nearest point to it.

*The F38 table.* Reads survey.json, not audit.json -- audit.json has no
sobolev_group leg, which is why the table was ever a manual join. The numbers
are exact to the digit (Pr II -6.4%, Ce II -4.3%).

*The exit_tau scan, and a parameter that does nothing.* The original F39 run
did not record dlnlam, so the driver scanned all three candidate values to find
which reproduced §4.34b. All three did, identically, to every printed digit.
At delocalize = 1 the exit line is placed anywhere in the forest and the
spacing parameter has no effect at all. Good to know; embarrassing that it took
a 30-minute scan to know it. Crossings 179.6 and 689.0 against the quoted 179
and 688.

*Real filters.* DECam g r i z and 2MASS J H Ks from SVO, integrated as
per-bin weights over the 200-bin spectral histogram rather than sampled at bin
centres (DECam g covers 19 bins; centre sampling would be a different
instrument). Gate 1 -- re-photometer the four committed F41 spectra, no packets
-- passes: worst colour errors move by -0.10 to +0.03 mag and stay 0.65-0.77
mag on Ce II and the blend where the top-hats gave 0.74-0.85. On the blue
kilonova the binned leg's worst colour moves from g-r to i-J at the same epoch
and about the same size. The redistribution leg stays at 0.008-0.009 mag. One
criterion had to be restated before it could be applied: "A_redist <= 0.02 mag"
is meaningless on the blue-kilonova file whose *top-hat* floor is already
0.098 (§9al), so the gate reads A_redist <= max(0.02, 1.5 x top-hat). Declared
in the driver's docstring, not silently.

## 10. Standing environment notes

- Everything SEDONA lives *outside* this repo: code `~/personal/pubsed`,
  deps `~/personal/sedona-deps`. Rebuild: `cd ~/personal/pubsed/src && bash
  install.sh wsl`. The `Makefile.wsl` that build needs is committed here at
  `docs/sedona/Makefile.wsl` -- copy it into `src/makefiles/` first. Full
  from-scratch setup for a new machine: `docs/sedona/SETUP.md`.
- Experiment drivers find SEDONA through `SEDONA_HOME` (default
  `~/personal/pubsed`) or `SEDONA_EXE`. Never hardcode a home directory.
- Figures are working outputs in `outputs/` (gitignored); the copies the
  report references live in `docs/figures/` (committed).
- SEDONA run directories under `experiments/**/run_*/` are gitignored;
  generators + param templates are committed, so every run is reproducible.
- MC noise at `core_n_emit = 2e6` is ~1–2% per band flux; raise it before
  chasing sub-percent effects.
