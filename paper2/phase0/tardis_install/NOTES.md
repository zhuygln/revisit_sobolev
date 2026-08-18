# TARDIS install and smoke test (de-risking P2-0D)

**Verdict: installed and importable; the branching physics is present; an
end-to-end run is BLOCKED on atomic data, not on the code.**

Recorded 2026-08-18. TARDIS is Paper II's candidate *reference* leg because
the public SEDONA has no fluorescence at all (see `../sedona_source_audit/`).
It is not on the critical path for P2-0B, which was carried by the custom
branching Monte Carlo in `../three_level_atom/`.

## What is installed

```
prefix   ~/personal/tardis-env          (micromamba, Python 3.14)
version  TARDIS 2026.8.10.dev7+g3990eff76
```

Three install routes were tried, in order:

| route | outcome |
|---|---|
| `micromamba -c conda-forge tardis-sn` | **fails** — no such package on conda-forge |
| `pip install tardis-sn` | **fails** — the PyPI name is a ~2015 stub that still calls `ez_setup.py` and dies fetching a dead setuptools tarball |
| repo `conda-linux-64.lock` + `pip install git+https://github.com/tardis-sn/tardis.git` | **works** — this is the only supported path |

Reproduce with
`scratchpad/tardis_install2.sh` (lockfile from
`raw.githubusercontent.com/tardis-sn/tardis/master/conda-linux-64.lock`).

## The branching physics Paper II needs is there

```
tardis/io/configuration/schemas/plasma.yml   line_interaction_type:
                                               scatter | downbranch | macroatom
tardis/opacities/macro_atom/                 the branching implementation
tardis/io/atom_data/macro_atom_data.py       its atomic-data side
```

Note the package moved: it is `tardis.opacities.macro_atom`, **not**
`tardis.plasma.properties.macro_atom` as older documentation and notebooks
have it. An import check written against the old path reports a false
negative.

This is the qualitative difference from SEDONA and the reason TARDIS is worth
keeping in the plan: `downbranch` and `macroatom` are first-class,
configuration-level options, not modifications.

## What blocks an end-to-end run

Every mode — including plain `scatter` — fails at atomic-data load, before
any transport. The code was never reached, so **nothing here is evidence
about the modes themselves.**

1. **TARDIS's own downloader is broken upstream.** The URL it holds,
   `media.githubusercontent.com/media/tardis-sn/tardis-regression-data/main/atom_data/kurucz_cd23_chianti_H_He_latest.h5`,
   returns 404. `raw.githubusercontent.com` serves the Git-LFS *pointer*
   (134 bytes), and querying the LFS batch API with the pointer's oid returns

   ```
   {'code': 404, 'message': 'Object does not exist on the server'}
   ```

   The same is true for `chianti_He.h5` and `kurucz_atom_chianti_many.h5`, so
   the repository's LFS objects appear to be gone repo-wide, not just for one
   file. `tardis-refdata` (the older location named in the documentation) no
   longer exists either.

2. **The one file still retrievable is too old to load.** `tardis-atomdata`
   still has working LFS, and
   `kurucz_cd23_chianti_H_He.h5.zip` (27 MB) downloads via the LFS batch API.
   But it is `database_version = v0.9`: plain HDF5 datasets
   (`basic_atom_data`, `levels_data`, `lines_data`, `macro_atom_data`, ...),
   not the pandas-HDF tables the current `AtomData.from_hdf` reads. It has no
   pandas keys at all, so `store.select(...)` raises
   `TypeError: cannot create a storer ...`. TARDIS's own error message
   correctly guesses "old-format of the atomic database" and points at the
   repository whose LFS objects are missing.

3. **Carsus, which generates the modern format, is not installable by pip.**
   `pip install carsus` → "no matching distribution"; it is GitHub/conda only,
   and would additionally need the external Kurucz and CHIANTI sources.

## What this means for P2-0D

The blocker is a data-distribution problem upstream, not a defect in TARDIS
and not something about this machine. Two ways forward, and the second is
better:

1. Obtain a current carsus-generated atom file (install carsus from GitHub, or
   get a file from a TARDIS user). Unblocks the standard verysimple config,
   which is only a smoke test.

2. **Write the atom file directly.** Paper II never needed Kurucz H/He: it
   needs a three-level synthetic atom, then a La II mini-atom from the GSI
   database. Both must be authored either way, and the required format is a
   documented set of pandas tables (`atom_data`, `levels`, `lines`,
   `macro_atom_data`, `macro_atom_references`) listed in
   `tardis/io/atom_data/base.py:152`. This is the same task already solved
   once for SEDONA in `experiments/laII_forest/` — and there the lesson was
   that group attributes omitted from the atom file segfault the reader, so
   the format is worth building against a loader test rather than a run.

Path 2 also keeps TARDIS in the role the project can defend: an independent
cross-check on a controlled atom, not a production run whose atomic data we
did not author and cannot audit.

## Scope

Install and smoke test only, per plan. No three-level TARDIS atom was built
in this pass; P2-0B was delivered by the custom Monte Carlo instead.
