# Paper IV — outline (draft skeleton, 2026-09-11; the PI decides whether this becomes the manuscript)

**Working title.** Fluorescence breaks the equivalence of coarse line-opacity treatments in kilonova radiative transfer

**Question.** How well do expansion and line-binned opacity closures survive when fluorescence is treated explicitly and energy is conserved?

**One-sentence result.** Under complete thermal redistribution both closures reproduce resolved Sobolev transport to a few tenths of a magnitude (line-binned the closer), as Fontes et al. (2020) found; with the single change to energy-conserving fluorescence the expansion closure keeps its 0.1–0.3 mag error while the line-binned closure develops colour errors of 0.3–1.7 mag in snapshots and 0.45 mag over a light curve, and on a strong-lined forest its transport does not terminate.

## Sections

1. **Introduction.** The two coarse-grainings of kilonova opacity (EP93 expansion opacity; area-preserving line binning), the Fontes et al. (2020) validation under ε = 1, ARTIS/Shingles line-by-line fluorescence without the controlled comparison, the SEDONA/SN Ia precedent; the gap this closes. The information-loss framing: a binned opacity does not retain the absorbing transition.
2. **Method.** The Sobolev Monte Carlo (indivisible energy packets, comoving energy conserved, the downward macroatom with β once, the exact energy scale); the three opacity treatments on one line list (resolved; expansion Σ(1 − e^{−τ}) per bin; line-binned Σ τ per bin, with the absorbing line restored within the bin); the two redistributions (complete thermal; downward macroatom); the radiation-field-driven macroatom (upward jumps under an imposed W B_ν(T)) as a robustness test, not an equilibrium solution; multi-shell zoned transport; the time-slab light curve (pause/resume, absolute heating energy, escape times). The inversion bug and its regression tests belong in the code section, briefly, with the split-invariance test as the diagnostic that found it.
3. **Benchmarks.** P1 (xkn secular, Gillanders Ye-0.21a pattern at X_lan 0.11, 1–5 d), P2 (AT2017gfo-like, Ye-0.29a rescaled), the Fontes Appendix C problem (pure Nd) as a snapshot and as a light curve; two line lists (GSI calibrated; Japan-Lithuania).
4. **Results.**
   4.1 F60/F61: the corrected closures on published states — Gate 2 Red everywhere; expansion within 0.2 mag; grid-invariant (Fig. F62 matrix; the multi-shell table).
   4.2 F62: ε = 1 vs fluorescence on P1/P2 — line-binned ≤ 0.1–0.3 under ε = 1, +0.3–0.8 under fluorescence; expansion unchanged (Fig. `fig_f62_matrix.png`).
   4.3 F64: the Fontes snapshot, two brackets, two line lists — line-binned colour error 0.0–0.2 → 1.4–1.7 mag; non-termination on the thinner forest (Fig. `fig_f64_snapshot.png`).
   4.4 F65: the light curve — bolometric ordering the same under both redistributions (resolved lowest, line-binned within 2–4 %, expansion +15–17 %); colour errors 0.13 → 0.45 (line-binned), 0.20 → 0.32 (expansion); the radiation-temperature variant (Fig. `fig_f65_lightcurve.png`; Table: peaks, E_rad, W, residuals).
   4.5 F63: the full-macroatom robustness test (Fig. `fig_f63_macroatom.png`).
5. **Discussion.** Why the closures fail differently (expansion: a nearly grey offset from Σ(1 − e^{−τ}) under-absorption; line-binned: over-absorption in saturated bins that complete redistribution forgives and fluorescence does not — the wrong upper level is activated); the information-loss interpretation of the non-termination; what the light curve adds (trapping time, expansion work: the resolved leg loses the most energy to work); relation to Fontes' 5–8 % and to Shingles' argument; what is not shown (no energy equation; T prescribed or radiation-fed; pure Nd; two line lists; no observer-time light curves beyond 12 d; g and r starved).
6. **Conclusions.**

## Figures and tables
- Fig. 1 `docs/figures/paper4/fig_f62_matrix.png` — the 2×2 closure/redistribution matrix on P1 (1–5 d), the robustness patterns and P2.
- Fig. 2 `fig_f64_snapshot.png` — the Fontes 4-d snapshot, both brackets, both line lists.
- Fig. 3 `fig_f65_lightcurve.png` — the six light curves, band and colour residuals.
- Fig. 4 `fig_f63_macroatom.png` — the full-macroatom check.
- Table 1 — the states and their saturation; Table 2 — the light-curve peaks, E_rad, W and residuals (§4.60); Table 3 — the grid-invariance sequence (§4.55–4.56).
- SI: the inversion bug (§4.55), the golden re-pin, the split-invariance and analytic-attenuation tests; the cost table; the dead-end and deposit variants; the Paper III correction as a note.

## Claims → findings (each with its gray condition in the report)
- "Expansion opacity within ~0.2 mag of resolved transport on published states" → F60, F61 (§4.55–4.56).
- "Line-binned ≤ 0.1–0.3 under ε = 1, 0.3–0.8 under fluorescence" → F62 (§4.57).
- "The Fontes ordering re-found; colour error 1.4–1.7 mag with fluorescence; non-termination on the thinner forest" → F64 (§4.59).
- "Over the light curve: bolometric ordering unchanged, colour errors 0.13 → 0.45 and 0.20 → 0.32" → F65 (§4.60; R1 Red on the amplitude of the expansion excess, R2 Gray on the cap).
- "Not an artefact of the downward table" → F63 (§4.58).

## Open before submission
- The line-binned fluorescence leg's cap (1.4 % of packets in the first slabs): a higher cap or a per-line escape treatment inside the bin to make R2 readable.
- A second light-curve model (P1's composition at its own heating) to show the light-curve result is not Nd-specific.
- The transient out-of-range level in one radiation-T slab (did not reproduce; a range check is in place).
