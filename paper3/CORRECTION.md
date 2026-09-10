# Paper III — scientific correction (2026-09-10)

> **Scientific correction — September 2026.** The grouped-opacity results
> reported here are **invalid** because of an implementation error in the
> inversion of cumulative binned opacity (`paper2/phase1/forest_mc.py::nu_of_G`,
> fixed in commit 38ebf30). The error caused packets in expansion-, line-binned-
> and dual-opacity modes to skip portions of non-smooth line forests. The
> resolved Sobolev calculations are unaffected. After correcting the transport
> and rerunning the analysis with energy-conserving transport on published
> ejecta benchmarks, the reported 1–3 mag grouped-opacity discrepancy is not
> reproduced. This manuscript is retained as a historical research record and
> is not a publication candidate. Corrected results: the Paper IV campaign
> (`docs/results_report.md` §4.55–4.59, F60–F64); the list of affected
> findings is in [`paper3/CORRECTION.md`](CORRECTION.md). The
> frozen tag `paper3-freeze` and the PDFs are left as they were.

## What the error was

`paper2/phase1/forest_mc.py::nu_of_G` (and its zoned copy) inverts the
cumulative expansion opacity G(ν) to find the frequency at which a packet's
drawn optical depth is reached. It located the bin as `nb − 1 − m` with
`m = searchsorted(G_edges[::-1], g, "right")`; the bin satisfying
`G_edges[b+1] ≤ g < G_edges[b]` is `nb − m`. One bin too low, the within-bin
fraction was formed with the neighbour bin's opacity: next to a thinner bin
the target overshot upward, landed above the packet's own comoving
frequency, was discarded as "behind", and the packet skipped the rest of the
forest. On smooth forests (every test toy, the synthetic forest of §4.31)
neighbouring bins have similar opacity and the error is a fraction of a bin;
on real lanthanide forests at the 4×10⁻⁵ bin width it leaked at every leg.
Found by multi-shell transport (the same physical state gave a different
closure answer with a different number of shell boundaries;
`docs/results_report.md` §4.55). Fixed in 38ebf30; the pre-fix golden
histories are kept in `tests/data/golden_run_mc_prefix_2026-09-10.json`; the
regression test is
`tests/test_zoned_run_mc.py::test_bin_legs_are_invariant_under_splitting_a_shell_on_a_spiky_forest`.

## Affected findings (README findings register)

Invalid — produced with the bug on real line forests; do not use
quantitatively or scientifically: F21, F22, F23, F24, F30, F31, F33, F35, F36, F38, F40, F41, F42, F43, F44, F45, F46, F47, F48, F49.
The manuscript headlines are F30, F43 and F44.

Affected in principle, effect unverified — smooth synthetic forests, where
the bug's effect is small but was not measured: F34, F37, F39.

Unaffected — Sobolev legs only (the redistribution results and the
references): F1–F20, F25–F29, F32; and F50 (Paper IV, Sobolev references).

Paper IV rows F51–F56 and F59 were produced with the same bug and are
superseded by F60–F64.

## What stands

The half of Paper III that compresses the *redistribution* physics (F25–F29,
F32: a small group-to-group kernel reproduces explicit A·β branching) used
identical exact line opacity in both legs and is unaffected. The corrected
picture of the *opacity* closures is in the Paper IV campaign: the expansion
closure is within ~0.2 mag of energy-conserving resolved transport on every
published state tried, the line-binned closure within 0.1–0.3 mag under
complete thermal redistribution and 0.3–1.7 mag off once fluorescence is
explicit (F60–F64).

## Provenance

`paper3-freeze` (2026-09-03) and `paper3-freeze-x4ln` (2026-09-09) are not
moved; `paper3/FROZEN.json` is not edited; `manuscript.pdf` and `si.pdf` are
the frozen artefacts and are not rebuilt. Every notice in this repository
is additive.
