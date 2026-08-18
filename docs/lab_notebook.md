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
