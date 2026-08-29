# Paper III — how much redistribution information is actually required?

Paper II established that no scalar thermalisation parameter reproduces
lanthanide fluorescence (F22), that the calibration is ion-dependent (F24),
and that the branching-aware Poisson closure is density-limited (F24).
This directory asks the constructive question: can a small group-to-group
redistribution matrix R_ij -- absorbed group in, re-emitted group out, no
atomic levels inspected after absorption -- reproduce explicit A*beta
branching? See `plan.md` for the full program and gates.

Design decisions taken in implementation (recorded here because the plan
leaves them open):
- The kernel is EVENT-level: one row draw per absorption event, where an
  "event" is the branch leg's full re-absorption chain collapsed to
  (nu_absorbed_cm -> nu_exit_rest), exactly what `run_mc` already resolves.
  Photon-number rows drive the sampling (each event re-emits one photon);
  the energy matrix R^E with q_dep = net comoving deposit per absorbed
  energy (which can be negative -- blueward fluorescence) is the
  conservation/validation object, and sum_j R^E_ij + q_dep_i = 1 is exact
  by construction per event. Core loss is transport's job, not the
  kernel's, so q_core lives in the reference metadata only.
- Within-group emission frequency: per-output-group sub-histogram (16
  sub-bins) marginalized over input groups, with a uniform fallback; the
  (i,j)-resolved sub-PDF is a second-order refinement left to Phase 3
  analysis if morphology demands it.
- Rows never populated in the reference fall back to coherent scattering
  (keep nu_cm) and are counted; they carry ~0 of the absorbed energy.
- Paper II code is imported, not forked: `forest_mc.py` gains one mode
  ("sobolev_group") and an event collector, both inert for every existing
  leg (classical rng stream untouched).

Layout: redistribution/ (kernel), phase0_reference/ (frozen reference +
Gate 0), phase1_groups/ (R1-R3 compression sweep), phase2_ions/ (R4 Ce II;
R5 Nd II blocked on data). Tests in the repo-wide tests/ so the suite runs
them.
