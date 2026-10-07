# Paper B — the effective redistribution operator

The PRL attempt: can detailed energy-conserving lanthanide fluorescence be
replaced by a compact effective operator R_ij, with the resolved Sobolev
opacity held fixed? The PI's statement is verbatim in
[`plan_review.md`](plan_review.md); the map is [`plan.md`](plan.md); the
first production gate is preregistered in [`prl_gate.md`](prl_gate.md) and
is not edited after the run.

## Gate G1 (the kill gate)

    .venv/bin/python paperB/gate1/run_gate1.py --ion 57LaII      # then 58CeII, 60NdII (~10-20 min each)
    .venv/bin/python paperB/gate1/analyse.py                      # readings B1-B3, gray first -> gate1_verdict.json
    .venv/bin/python paperB/gate1/figure.py                       # docs/figures/paperB/gate1_error_vs_complexity

`run_gate1.py` runs one ion alone on the P1 2 d photospheric state through
the energy-conserving downward macroatom reference (`R2`, evaluation seeds
1–3), R_ij at N_g = 2, 4, 8, 16, 32 (kernels built from a separate
reference run on build seeds 101–103, discrete within-group tables, energy
rows, transported on the evaluation seeds), the independent 128-group
matrix for the event-level metric, and ε = 0 … 1 in steps of 0.05, with
3×10⁵ energy packets per seed. `analyse.py` refuses a record that
is not the preregistered experiment. Tests: `tests/test_paperB_gate1.py`.

## The referee revision (2026-10-06): G2R, G3U, the minimax scalar, the costs

    .venv/bin/python paperB/gate2r/run_g2r.py --ion 58CeII --k 16   # G2R: NOT before the PI approves prl_gate.md's G2R section
    paperB/gate2r/run_all.sh                                        # Ce II then Nd II, N_g = 2 ... 32 (~9 h, sequential, retried)
    .venv/bin/python paperB/gate2r/analyse.py                        # rank of the physical ordering per (ion, N_g) -> g2r_verdict.json; --markdown; figure.py
    .venv/bin/python paperB/gate3u/analyse.py --affected             # the mechanically defined near-threshold G3 legs (33 in 18 records)
    .venv/bin/python paperB/gate3u/run_g3u.py --record gate3_P_P3d_60NdII.json   # G3U: one frozen record, 12 common seeds
    .venv/bin/python paperB/gate3u/analyse.py                        # paired bootstrap per case -> g3u_verdict.json; --markdown
    .venv/bin/python paperB/scalar/minimax.py                        # the joint minimax eps from the existing grids (no transport); --markdown
    .venv/bin/python paperB/cost/costs.py                            # offline / stored / online cost layers from the records; --markdown

**G2R ran 2026-10-06/07 (F71, §4.66):** the physical ordering reproduces
the light better than all 31 scrambled orderings in 8 of 10 (ion, block
count) cells and in every cell at or above the count where the ion passes
(p = 0.031 each); Ce II Green, Nd II Yellow (second at its two-block
passing count by 0.001 mag, inside the seed noise), **overall Yellow — the
PI's call**. On Nd II at 32 blocks the physical ordering fits the events
worse than most scrambled orderings and the light three times better than
the best: the inversion at identical rank and parameter count.

**G3U ran 2026-10-07 (F72, §4.67):** of 33 near-threshold decision legs
on twelve paired seeds (operators matched at 0.0), 5 decided pass, 13
decided fail, 13 within noise, 2 gray; the decisive readings are decided
(the Nd II trajectory interpolant fails, 0.118 [0.113, 0.123]), Ce II's
3×10⁵-packet readings near the threshold are noise-limited (7 point
estimates cross), and the frozen P1d Ce II reference at 10⁶ packets is not
reproduced seed for seed while its closures are (gray; the letter
unchanged; `gate3u/nonreproducible_p1d.md`). No frozen letter is re-read.

`gate2r/` is the frequency-adjacency ablation the PI preregistered on the
referee's reading (`operators_perm.py`: the 128 fine groups permuted, the
same block coarse-graining as G2's local family, the permutation undone;
32 orderings per ion and block count, read by the rank of the physical
ordering). `gate3u/` re-transports the near-threshold G3 decision legs
on twelve common seeds with the trained operators preserved and attaches
a paired-bootstrap interval to the joint decision statistic (it annotates
G3, it does not re-gate it). `scalar/` and `cost/` are post-hoc readings
of the frozen records. Tests: `tests/test_paperB_gate2r.py`,
`test_paperB_gate3u.py`, `test_paperB_minimax.py`, `test_paperB_cost.py`.

## After G1: R2M robustness, the exit-table audit, G2 (preregistered, not run)

    .venv/bin/python paperB/r2m/run_r2m.py --ion 57LaII          # then 58CeII, 60NdII (10 s / 2.5 min / 6 min at 3e5)
    .venv/bin/python paperB/r2m/analyse_r2m.py                    # survives? per ion -> r2m_verdict.json; --markdown for the report's tables
    .venv/bin/python paperB/audit/exit_tables.py --ion 57LaII     # representation audit (no transport scored); --figure; --markdown
    .venv/bin/python paperB/gate2/run_gate2.py --ion 58CeII       # G2: NOT before the PI approves prl_gate.md's G2 section
    .venv/bin/python paperB/gate2/analyse.py                      # H1/H2, gray first -> gate2_verdict.json; figure.py
    .venv/bin/python paperB/gate3/support.py                      # G3: freeze one support per ion (committed before the run)
    paperB/gate3/run_all.sh                                       # G3: Ce, Nd, La (anchor, then axes T D J P), then parts b and c (~10 h)
    .venv/bin/python paperB/gate3/analyse.py                      # A/B/C/D per state, C1-C5 -> gate3_verdict.json

`r2m/` is the post-pass robustness check of G1's own text (the
radiation-field-driven macroatom as the reference, the 8-group operator
rebuilt from its events; not a gate). `audit/` measures what the discrete
exit tables hold (distinct exit lines, exact rest frequencies, discovery
curve, energy concentration, what a truncation keeps, size split).
`gate2/` holds the operator families (`operators.py`: local on shared
128-group tables, rank-k NMF, exit-table truncation), the runner, the
readings H1/H2 and the figure. `gate3/` holds the state domain
(`states.py`: the T, D, J, P axes), the frozen supports
(`gate3_support.json`), the whole-operator interpolation and the support
diagnostic (`operators3.py`), the runner (transfer at R₁₆ from the θ₀
anchor, existence from the operators rebuilt at θ, the interior
interpolation from the endpoints' saved kernels), and the readings A/B/C/D
and C1–C5. Tests: `tests/test_paperB_r2m.py`, `test_paperB_audit.py`,
`test_paperB_gate2.py`, `test_paperB_gate3.py`, `test_paperB_freeze.py`.

## Where it stands

G1 was frozen at f9cfbff and run on 2026-09-22 (`docs/results_report.md`
§4.62, F67): N_g* = 2 (La II), 16 (Ce II), 2 (Nd II, after its
preregistered 10⁶-packet rerun cleared the precision gray); ε* is the
coherent limit and 1.1–1.4 mag off on the dense ions. **B1 Green, B2
Green, B3 Yellow → continue toward the PRL.** La II at its P1 partial
density is too thin to discriminate.

After G1 (§4.63, F68, 2026-09-22): **the compression survives the
radiation-field-driven macroatom** on all three ions (the 8-group
operator rebuilt from its events within 0.015–0.046 mag, the reference
itself moving by up to 2.1 mag; the downward-trained operator reproduces
the downward reference it came from). The exit tables are already the
exact aggregate over distinct exit lines; their size is the number of
exit lines the sample discovers, and 90 % of the exit energy sits in
5–26 % of them. G2 was preregistered, approved with the PI's H1 amendment, tagged
`paperB-g2-prereg` and run on 2026-09-22 (§4.64, F69): **H1 Green, H2
Yellow → write the mechanism claim.** K*_local = 2 / 16 / 2 against
K*_global = 2 / 32 / undefined (La II / Ce II / Nd II), so local
coarse-graining needs no more archetypal exit distributions than a global
non-negative factorisation — on Nd II two contiguous frequency blocks
against a rank-32 mixture that never reaches the threshold. At 32
archetypes the global operator fits the microscopic redistribution 2–4×
better and reproduces the light 2–20× worse. The exit spectrum compresses
to L* = 0.358 (Ce II) and 0.113 (Nd II) of its distinct lines, short of
the 10 % headline. Two preregistration amendments made during execution
(two implementation defects; the control firing and its declared remedy)
are recorded in `prl_gate.md`.

G3 was preregistered with the PI's amendment (transfer/interpolation at
N_g = 16, one frozen support per ion, whole-operator interpolation),
tagged `paperB-g3-prereg`, and run 2026-09-24-26 (§4.65, F70). **Lanthanide
fluorescence remains strongly compressible across physical state and
realistic composition, but the effective operator is state dependent:
one-dimensional interpolation predicts most state changes, while coupled
trajectories and species mixing expose nonlinearities that require
refitting or richer tabulation.** Transfer fails on every axis for Ce II
and on all but the source-spectrum axis for Nd II and La II;
whole-operator interpolation recovers every other
failure except one genuinely *different* outcome, not a variant of the
same one -- **Nd II's trajectory axis is outcome C**: direct transfer
fails, interpolation also fails at the interior state, and only a fresh
fit survives, so the operator manifold is compact but not separably
interpolable when several coordinates move together. Ce II is a distinct
qualification, not outcome D: its fresh 16-group operator fails at 7/14
states, always recovered at 32 -- compactness survives, a *fixed*
representation size does not. C4 (species mixing needs its own fit) reads
Yellow and is physically understandable once read against **C5 Green**:
the 13-ion blend passes at N_g = 4 while eps* misses by 1.69 mag, the
worst scalar failure in the program -- the effective-operator idea
survives a realistic blend even though the individual-ion state
dependence is nontrivial. No neural surrogate is warranted yet (one
failure of separable interpolation is not evidence against a better-chosen
low-dimensional table).

**The manuscript is a complete draft under referee revision**
(`docs/paperB/manuscript.tex`; the PRL reframing of 2026-10-05, merged as
PR #14; the referee revision of 2026-10-06): every section written on the
frozen G1–G3 record in the PI's progression, every number a macro from
`paperB/FROZEN.json`, the claim read-through (`make claims`) clean. The
referee's central criticism is accepted (the PI, `plan_review.md`
2026-10-06): the inversion is a 32-archetype result and does not by
itself separate adjacency from low rank; G2R tests that, G3U attaches
paired-seed intervals to the near-threshold G3 readings, and the wording
("by construction", "every run", "not separable", "best scalar", the
model equation) is corrected. **G4 is demoted** to a post-manuscript
robustness candidate and is not a gate. The Paper IV MNRAS submission is
independent.

Three implementation issues surfaced during the run, all disclosed in full
in §4.65 and `paperB/gate3/bug_snapshot/`: a confirmed, fixed, disclosed
bug in `thermal_sampler` (unreachable in any prior completed record); four
crashes checked and found to be a different, unresolved phenomenon; and a
third, energy-identity anomaly that the PI held PR #13 to audit before
merging -- closed 2026-09-27 (`energy_identity_anomaly.md`): all three
flagged legs reproduced clean on every seed (non-reproducible, not a logic
bug), the corrupted run's photometry was shown indistinguishable from a
clean rerun, and the real defect -- `read_blend()`'s C4 search not
excluding a leg its own `gray_checks()` had already named -- is fixed
(`leg_invalid`, gray-first, every exclusion disclosed). C4 reads Yellow
unchanged on the existing records; no family needed rerunning; C1-C3 and
C5 were not reopened.
