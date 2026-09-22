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

## Where it stands

G1 was frozen at f9cfbff and run on 2026-09-22 (`docs/results_report.md`
§4.62, F67): **Yellow** — La II N_g* = 2 (too thin to discriminate), Ce II
N_g* = 16 with ε* 1.40 mag off, Nd II Gray on the precision rule (its 10⁶-
packet rerun pending). The PI decides with the numbers in hand; G2–G4 are
detailed in `prl_gate.md` only after that decision.
