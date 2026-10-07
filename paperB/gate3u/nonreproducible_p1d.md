# G3U found a second non-reproducible run: the P1d Ce II reference at 10⁶ packets (2026-10-07)

The G3U determinism check (prl_gate.md "G3U": the per-seed scatter over
seeds 1–3 of every rerun leg must equal the frozen record's
`mags_seed_std` to 10⁻⁶ mag) fired on one leg of one record:

| record | leg | frozen git | reproduced? |
|---|---|---|---|
| `gate3_P_P1d_58CeII_n1e6.json` | `R2` (the reference, seeds 1–3, 10⁶ packets) | 0b34f2f (2026-09-24) | **no** — r 24.0736 / i 23.4873 frozen against 24.0840 / 23.4805 rerun (seeds 1–3 pooled), seed scatter 0.0090 / 0.0007 against 0.0244 / 0.0111 |
| same | `R2build` (build seeds) | | yes — `n_interactions` 18,708,447 both, kernels identical at 0.0 |
| same | `Arec_ng16`, `Arec_ng32` | | yes — magnitudes and scatter identical to the printed digits, operators matched at 0.0 |
| `gate3_T_T5000_58CeII_n1e6.json` | `R2` | 0b34f2f (same day, same commit) | yes — scatter deviation 0.0 |
| every other rerun record (16) | every leg | | yes |

The frozen P1d `R2` closes its energy identity (3.7×10⁻¹⁷), was not gray,
and differs from a deterministic rerun of the same code, state, seeds and
packet count by up to 0.010 mag in a band. The build reference in the
same process reproduced exactly. This is the class the energy-identity
audit of 2026-09-27 (`paperB/gate3/bug_snapshot/energy_identity_anomaly.md`)
named — a run altered by the intermittent host fault without crashing —
now caught by a determinism check rather than by an exploding ledger.
The same chain's part (c) rerun died on its first attempt with the bit-60
`IndexError` (index 2⁶⁰ + 1,086,790 on an array of 375,912; retried,
completed clean).

**Consequence for the frozen readings: none.** The P1d Ce II state reads
D on the frozen record (fresh R16 0.119 / 0.151 fails, R32 passes) and
reads the same against the reproducible 12-seed reference (R16 0.137
[0.127, 0.149] decided fail; R32 0.075 [0.062, 0.089] decided pass). The
two P1d cases are reported GRAY under the preregistered rule, with their
12-seed numbers as diagnostics. No frozen record is edited.
