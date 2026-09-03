# Paper III — the campaign directory

This directory holds the **code, run outputs and frozen record** of Paper III.
The **manuscript** lives in [`docs/paper3/`](../docs/paper3/) and is generated
from what is here: every number in it is a LaTeX macro produced by
`docs/paper3/latex_tables.py` from [`FROZEN.json`](FROZEN.json), and every
display item by `docs/paper3/display_items.py` from the JSONs this directory
records. The two are one paper, not two — same numbering, same tag
(`paper3-freeze`), split the way the repo splits every paper: `paperN/` is the
campaign, `docs/paperN/` is the write-up.

Manuscript title (Nature Astronomy Article format): *Coarse-grained line
opacity leaves a detectable chromatic signature in kilonovae that ejecta
parameters cannot mimic.*

## The question, and how it moved

The directory opened on the constructive question Paper II left: can a small
group-to-group redistribution matrix R_ij -- absorbed group in, re-emitted
group out, no atomic levels inspected after absorption -- reproduce explicit
A*beta branching? [`plan.md`](plan.md) is that program as received
(2026-08-29), kept verbatim as the historical input.

The answer came out yes (F25, F27): redistribution compresses. That turned the
campaign toward the approximation that does *not* compress -- the grouped
opacity (F30/F31) -- and then toward what its error does to an observer: the
kilonova trajectory (phase10), the observables (phase11), the (M_ej, v_ej,
X_lan) grid (phase12) and the observing scenarios (phase13). The manuscript is
the second half of that arc with the first half as its safe-closure control
leg, which is why its title is about opacity and `plan.md` is about
redistribution. Both halves are here.

## Design decisions taken in implementation

(Recorded here because the plan leaves them open.)

- The kernel is EVENT-level: one row draw per absorption event, where an
  "event" is the branch leg's full re-absorption chain collapsed to
  (nu_absorbed_cm -> nu_exit_rest), exactly what `run_mc` already resolves.
  Photon-number rows drive the sampling (each event re-emits one photon); the
  energy matrix R^E with q_dep = net comoving deposit per absorbed energy
  (which can be negative -- blueward fluorescence) is the
  conservation/validation object, and sum_j R^E_ij + q_dep_i = 1 is exact by
  construction per event. Core loss is transport's job, not the kernel's, so
  q_core lives in the reference metadata only.
- Within-group emission frequency: per-output-group sub-histogram (16
  sub-bins) marginalized over input groups, with a uniform fallback; the
  (i,j)-resolved sub-PDF is a second-order refinement left to Phase 3 analysis
  if morphology demands it.
- Rows never populated in the reference fall back to coherent scattering
  (keep nu_cm) and are counted; they carry ~0 of the absorbed energy.
- Paper II code is imported, not forked: `forest_mc.py` gains one mode
  ("sobolev_group") and an event collector, both inert for every existing leg
  (classical rng stream untouched).

## Layout

Directory numbers are execution order, not `plan.md`'s P-numbers; the mapping
is in the last column. There is no `phase2_`: the plan's P2 is the transport
leg itself, which lives in `paper2/phase1/forest_mc.py` as the
`sobolev_group` mode rather than in a directory of its own. The ions of P4
were folded into `phase1_groups/` (La II, Ce II, Nd II side by side) instead
of the `phase2_ions/` an earlier draft of this file announced.

```
redistribution/       the RedistributionKernel class (from_branching_mc,
                      sample_output_group, validate_energy, save/load)   plan §4
phase0_reference/     the frozen reference + Gate 0, La/Ce/Nd            P0
phase1_groups/        the N_g = 4..128 compression sweep and the ions    R1-R3, P1, P3, P4
phase3_temperature/   T_src and T_gas transfer, fixed vs recomputed      P5
phase4_epoch/         epoch dependence and the tau-collapse test         P6
phase5_mixture/       the opacity-weighted composition rule              P9
phase6_opacity/       kappa_grouped + R_ij: opacity vs redistribution    P11
phase7_rank/          how many independent redistribution modes (E1)     was P8
phase8_survey/        the collapse across thirteen real ions (E3b)       E3b
phase9_audit/         the normalization audit + counterfactual table     --
synthetic/            synthetic forests: crowding, saturation, spacing
                      and redistribution range dialled one at a time     E3/E4
phase10_kilonova/     does a real kilonova cross the cancellation
                      boundary                                           item 5
phase11_observables/  what the sign change does to an observer;
                      Gate 1 with real DECam + 2MASS passbands           §10
phase12_grid/         the (M_ej, v_ej, X_lan) heating-powered grid,
                      robustness, sensitivity, syserr, tscale, chain     P12
phase13_observability/ Gate 3: three pre-declared observing scenarios
figures/              working figures for docs/results_report.md (the
                      manuscript's own display items are in
                      docs/paper3/figures/)
freeze.py             regenerates every derived table and figure and
                      checks them against FROZEN.json
FROZEN.json           SHA-256 of every transport output, derived table
                      and figure, plus the commit and input tree hashes
plan.md               the program as received, kept unedited
```

Tests live in the repo-wide `tests/` so the suite runs them.

## Reproducing

```bash
python paper3/freeze.py --check --strict     # every derived quantity + figure
cd docs/paper3 && make                       # manuscript.pdf, si.pdf, structure check
```

The tag `paper3-freeze` marks the state the manuscript was built from. Because
the manuscript prose, `freeze.py`, `docs/paper3/display_items.py`,
`docs/paper3/check_structure.py` and that tag all name the literal path
`paper3/`, this directory's name is part of the frozen record: renaming it
means regenerating `FROZEN.json` and re-tagging.
