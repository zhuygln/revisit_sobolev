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
