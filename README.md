# revisit_sobolev

How large is the error introduced by the **Sobolev approximation** in r-process
(kilonova) ejecta, and where in kilonova parameter space does that error matter?

## Phase 0 status: COMPLETE — GO

All four weekend deliverables from [docs/babystep_plan.md](docs/babystep_plan.md) are done:

- **0A** — `tau_exact` reproduces the analytic Sobolev optical depth to machine
  precision (E_Sob = 2×10⁻¹⁶) with a monotone convergence study.
- **0B / Figure 1** — a tanh population gradient breaks the approximation on a clean
  ε² power law (ε = v_D/v_scale), reaching ~5% error at ε = 1. Finding: a gradient
  symmetric about the resonance cancels at leading order, so the resonance must sit
  off-centre — gradient-only validity diagnostics are conservative.
- **0C / Figure 2** — real GSI La II data (3000 Å–3 μm): 11.3% of the broader subset
  has a neighbour within 10 km/s (3.3% within 3 km/s); even the high-confidence
  subset has a 1.2% tail. The thermal-width crowding regime exists.
- **0D** — SEDONA (public `dnkasen/pubsed`) built serially-invoked under WSL2 (deps
  via a micromamba env + Lua 5.1.5 source build + OpenMPI, custom
  `src/makefiles/Makefile.wsl`); the 1-D `spherical_lightbulb` example propagates
  2×10⁵ particles in ~1 s and writes a spectrum. Both line treatments needed later
  exist in this release: `opacity_line_expansion` (Sobolev) and
  `opacity_bound_bound` (resolved), `defaults/sedona_defaults.lua:115,118`.

## Week 1 status (babystep_plan.md §18)

- **Voigt profile** — implemented via the Faddeeva function; normalization,
  Gaussian-limit, and Lorentzian-wing tests.
- **Boltzmann populations** — partition function and level fractions straight from
  the GSI levels files (parsed by the same `load_gsi`); analytic two-level and
  temperature-limit tests, plus the per-statistical-weight monotonicity invariant
  on real La II levels.
- **SEDONA Type Ia pair** (`param_d20_lte_exp` vs `param_d20_lte_bb`, run in
  `~/personal/pubsed/examples/supernova/TypeIa/spectrum/run_{exp,bb}/`): the two
  configs differ **only** in the line treatment — the built-in controlled
  experiment this project needs. Results (Figure 3, `outputs/`):
  - cost: **14 s (Sobolev) vs 88 min (resolved)** for the same 4 iterations —
    a ~378× ratio, empirically confirming the §11 frequency-resolution cost and
    the necessity of narrow wavelength windows;
  - physics: the resolved treatment blankets the UV harder and redistributes
    ~3× more flux into the NIR (8000–12000 Å) — qualitatively the effect that
    matters for lanthanide-rich kilonova ejecta (toy model; qualitative only).

Next: Week 2 — a minimal 1-element/1-line SEDONA model compared against the
Python deterministic solver, then 2 → 20 lines.

---

This repo started as a **Phase 0** weekend feasibility check, not a paper: reproduce
the Sobolev limit in a toy problem, break it deliberately, and look at whether real
lanthanide data occupy the failure regime. Everything else is deferred.

## Layout

```
sobolev/     toy solver: constants, profiles, optical depths, atomic data
notebooks/   00 single line (Phase 0A+0B), 01 GSI line spacing (Phase 0C)
tests/       profile normalization + the Phase 0A Sobolev-limit gate
data/        raw atomic data, gitignored except its README
docs/        the two planning documents
outputs/     figures (gitignored)
```

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

All 7 tests pass. The La II data files are not committed — see `data/README.md` for
the Zenodo provenance; re-download and drop them in `data/` to run notebook 01.

Go/no-go after each phase. See [docs/babystep_plan.md](docs/babystep_plan.md) for the
full roadmap and [docs/research_requirements.md](docs/research_requirements.md) for the
data and code landscape.
