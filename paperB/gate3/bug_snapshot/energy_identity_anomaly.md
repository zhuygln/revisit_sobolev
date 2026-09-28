# The three energy-identity anomalies -- bounded audit (2026-09-27)

Three legs' internal energy identity blew up to 3.9e70 / 3.9e70 / 9.8e31
while every other accounting fraction and the emergent photometry stayed
ordinary (`docs/results_report.md` 4.65, finding (3)). The PI held PR #13
for a bounded, five-step audit before merging. This file records it.

## The three exact records

| # | record | leg | mode | seeds | build seeds | identity residual |
|---|---|---|---|---|---|---|
| 1 | `paperB/gate3/gate3_D_D0.1_60NdII.json` | `Arec_ng32` | sobolev_group (kernel from R2build) | 1,2,3 | 101,102,103 | 3.859736307910541e+70 |
| 2 | `paperB/gate3/gate3_T_T5000_60NdII.json` | `R2` | sobolev_dmacro | 1,2,3 | 101,102,103 | 3.8702947938624614e+70 |
| 3 | `paperB/gate3/gate3_partb_blend3.json` | `Amix_ng32` | sobolev_group (injected mixed kernel) | 1,2,3 | 101,102,103 | 9.773952155896828e+31 |

Every other field of every one of these three legs' accounting is
ordinary: `dep_lab_frac`/`work_frac`/`core_frac` sum close to 1 as
expected, `renorm_ratio` is a normal small number, and (checked directly
in `forest_mc.py`) `identity_residual = (E_esc + E_core + E_abs + E_capped
+ E_dep_lab + E_carried_out - E_inj) / E_inj` matches `esc_frac = E_esc /
E_inj` to full float precision in all three cases -- meaning `E_esc`
alone, the raw sum `np.sum(w[fate==1] * H * nu_final[fate==1])`, is the
corrupted quantity; every other term in the sum is negligible beside it.

## Step 2: reproduction (the deciding step)

Each of the three legs was rebuilt from scratch -- the same state, the
same atom, the same kernel (rebuilt identically for #1 and #3), the same
seeds -- and run per seed with `run_mc` called directly so each seed's
raw `accounting` dict could be read individually (the committed record
only stores the 3-seed average).

    seed 1: E_inj 1.558371e-06  E_esc ...  identity  1.4e-16   (leg 1)
    seed 2: E_inj 1.558893e-06  E_esc ...  identity  0.0
    seed 3: E_inj 1.559221e-06  E_esc ...  identity  2.7e-16

    seed 1: E_inj 1.558371e-06  E_esc ...  identity -1.4e-16   (leg 2)
    seed 2: E_inj 1.558893e-06  E_esc ...  identity  2.7e-16
    seed 3: E_inj 1.559221e-06  E_esc ...  identity  1.4e-16

    seed 1: E_inj 1.558371e-06  E_esc ...  identity -1.4e-16   (leg 3)
    seed 2: E_inj 1.558893e-06  E_esc ...  identity  1.4e-16
    seed 3: E_inj 1.559221e-06  E_esc ...  identity  2.7e-16

**All three reproduce clean, at machine precision, on every seed.** This
is the key discriminator the PI's audit asked for: a genuine transport or
accounting *logic* bug would reproduce identically given the same code and
the same seeds (as the `thermal_sampler` bug did, 4/4, in the companion
finding). None of these three did. The anomaly is a one-off,
non-reproducible corruption event, of the same general character as the
four crashes already characterised and left open elsewhere in this gate
(`docs/lab_notebook.md` 9bn) -- but this time landing in a floating-point
accumulator (a packet's `nu_final` or `w`) rather than an integer array
index, so it poisons a sum instead of crashing.

## Step 3: did the corrupted leg's photometry reach the verdict?

Yes -- `Amix_ng32`'s `mags` (used for C4's band/colour comparison) came
from the same corrupted run as its identity residual, and C4's reading
does use those numbers, not the identity. Reproducing the leg's full
photometry (not just its raw accounting) and comparing to the frozen
record:

| | g | r | i | z | J | H | K | max abs dm vs R2 |
|---|---|---|---|---|---|---|---|---|
| clean rerun | 25.6237 | 24.0743 | 23.3650 | 22.7870 | 22.2176 | 21.5217 | 20.0478 | 0.1114 |
| frozen (corrupted-identity) run | 25.6469 | 24.0811 | 23.3712 | 22.7898 | 22.2210 | 21.5223 | 20.0452 | 0.1108 |

**The corrupted run's photometry is statistically indistinguishable from
a clean rerun** (both ≈0.111 mag max error against the reference, the
difference well inside ordinary seed-to-seed Monte Carlo scatter). The
astronomical identity blow-up did not measurably affect the numbers C4
was actually built from. This rules out a rescue-by-omission concern: the
leg's own verdict-relevant numbers were fine; only its self-check
diagnostic was corrupted.

## Step 4: root cause and fix

Genuine transport nonconservation is ruled out by the clean, deterministic
reproduction (step 2). "Bad denominator, overflow, wrong field, stale
accumulator" was also ruled out as a *code* defect for the same reason: a
bug in `E_inj`/`identity_residual`'s formula (checked directly against
`forest_mc.py`, both terms are read from the correct fields with no
aliasing) would reproduce every time, not intermittently. **This is the
"transport (and its accounting) is sound; only the validity calculation
[i.e. the C4 analysis, not the transport code] is wrong" branch**: the
existing raw records are trustworthy, and the fix is entirely in
`paperB/gate3/analyse.py`'s reading of them, not in `forest_mc.py`.

**The real bug, precisely stated**: `read_blend()` (part b/c's K* search)
computed `gray = A.gray_checks(row, M)` -- which correctly names the
specific bad leg (`"3: Amix_ng32 identity residual 9.8e+31"`) -- but never
used that information to exclude the leg before searching for the
smallest passing N. G1's own per-leg gray check and G3's per-state
`read_state()` both already enforce this; only the blend readings did not.

## Fix

`leg_invalid(m)` (the same four thresholds `gray_checks` already applies:
identity, kernel validation, chain-capping, empty-row fallback) now gates
every entry before it can be returned by `k_star()`, and every excluded
leg is disclosed in `invalid_legs` and per-entry in `table[...]['invalid']`
-- never silently dropped. `Amix_ng32` is now excluded from part (b)'s
search regardless of what its face-value band/colour say.

## Step 3 (impact) confirmed on the real record

`Amix_ng32`'s face-value numbers (0.111 mag band, 0.100 mag colour) sit
almost exactly on the 0.10 mag threshold -- close enough that reading it
at face value COULD plausibly have flipped `k_mix`. It did not, because
Amix already fails at every other tested N (2, 4, 8, 16: 0.93, 0.96, 0.39,
0.14 mag) regardless of what N = 32 does. **`k_mix` is `None` and `k_direct`
is 4 both before and after the fix; C4 reads Yellow either way.** No
rerun of part (b), or of any other family, is warranted.

## What was not reopened

Per the PI's explicit instruction, C1, C2, C3 and C5 are untouched: the
same audit question (does the anomaly reach a verdict undetected) does
not apply to them, since G3's per-state `read_state()` already excludes
an entire state as GRAY the moment any of its legs' identity fails (a
more conservative rule than the one just added to `read_blend`), and no
identity anomaly was found in any part-(c) leg.
