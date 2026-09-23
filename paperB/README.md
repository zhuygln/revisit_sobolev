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

## After G1: R2M robustness, the exit-table audit, G2 (preregistered, not run)

    .venv/bin/python paperB/r2m/run_r2m.py --ion 57LaII          # then 58CeII, 60NdII (10 s / 2.5 min / 6 min at 3e5)
    .venv/bin/python paperB/r2m/analyse_r2m.py                    # survives? per ion -> r2m_verdict.json; --markdown for the report's tables
    .venv/bin/python paperB/audit/exit_tables.py --ion 57LaII     # representation audit (no transport scored); --figure; --markdown
    .venv/bin/python paperB/gate2/run_gate2.py --ion 58CeII       # G2: NOT before the PI approves prl_gate.md's G2 section
    .venv/bin/python paperB/gate2/analyse.py                      # H1/H2, gray first -> gate2_verdict.json; figure.py

`r2m/` is the post-pass robustness check of G1's own text (the
radiation-field-driven macroatom as the reference, the 8-group operator
rebuilt from its events; not a gate). `audit/` measures what the discrete
exit tables hold (distinct exit lines, exact rest frequencies, discovery
curve, energy concentration, what a truncation keeps, size split).
`gate2/` holds the operator families (`operators.py`: local on shared
128-group tables, rank-k NMF, exit-table truncation), the runner, the
readings H1/H2 and the figure. Tests: `tests/test_paperB_r2m.py`,
`test_paperB_audit.py`, `test_paperB_gate2.py`.

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
are recorded in `prl_gate.md`. **G3 is untouched until the G2 verdict is
written up**, as the PI directed.
