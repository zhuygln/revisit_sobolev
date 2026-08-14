# revisit_sobolev

How large is the error introduced by the **Sobolev approximation** in r-process
(kilonova) ejecta, and where in kilonova parameter space does that error matter?

This repo is at **Phase 0** — a weekend feasibility check, not a paper. The goal is
to reproduce the Sobolev limit in a toy problem, break it deliberately, and look at
whether real lanthanide data occupy the failure regime. Everything else is deferred.

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

`pytest` currently **fails** on `test_sobolev_limit.py`, by design — those tests are
the specification for Session 1, not regression cover. `tau_exact` is a stub with the
derivation recipe in its docstring; implementing it is the first task.

## Phase 0 deliverables

1. **Figure 1** — `E_Sob` vs `v_D / v_scale`, showing the Sobolev limit and its breakdown
2. **Figure 2** — `P(dv_nearest)` for La II, with reference widths at 1, 3, 10 km/s
3. **SEDONA smoke test** — run `spherical_lightbulb`; stop if the install eats the weekend

Go/no-go after each phase. See [docs/babystep_plan.md](docs/babystep_plan.md) for the
full roadmap and [docs/research_requirements.md](docs/research_requirements.md) for the
data and code landscape.
