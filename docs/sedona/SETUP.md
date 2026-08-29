# Machine setup

Everything needed to take a fresh clone of this repo to "all 174 tests pass
and every experiment runs". Written 2026-08-29 while doing exactly that on a
second WSL2 machine; the paths below are the ones the code now defaults to.

Nothing here needs root.

## 1. Python

The package needs **Python >= 3.10** and **numpy >= 2**, so a distro python
3.8/3.9 will not do. The lab notebook's launch commands use an absolute
`.venv/bin/python`, so put the interpreter exactly there:

```bash
# if a >=3.10 python is already on PATH
python3 -m venv .venv

# otherwise (e.g. only python 3.9 present) any conda/micromamba will do:
conda create -y -p .venv python=3.12
```

Then:

```bash
.venv/bin/python -m pip install -e ".[dev]" h5py
.venv/bin/python -m pytest        # 174 passed
```

Tests that need atomic data `skip` cleanly when `data/` is empty, so a
data-less clone should report passes plus 5 skips and **no failures**.

## 2. Atomic data

`data/` is gitignored; re-download from the Zenodo record recorded in
`data/README.md` (875 MB for both archives):

```bash
cd data
curl -L -o GSI_lanthanides_calibrated_levels.zip \
  "https://zenodo.org/api/records/19335084/files/GSI_lanthanides_calibrated_levels.zip/content"
curl -L -o GSI_lanthanides_calibrated_transitions.zip \
  "https://zenodo.org/api/records/19335084/files/GSI_lanthanides_calibrated_transitions.zip/content"
```

Extract the per-ion files to `data/` **flat** (no subdirectory) — that is
where `load_gsi` looks. Both archives are Mac-made, so skip `__MACOSX/`
entries. All 27 lanthanide ions (La–Yb, II and III) are in there; the
experiments currently use:

| file | size |
|---|---|
| `57LaII_levels_calib.txt`, `57LaII_transitions_calib.txt` | 44 kB, 3.7 MB |
| `58CeII_levels_calib.txt`, `58CeII_transitions_calib.txt` | 258 kB, 84 MB |
| `58CeIII_*` | small, extracted for future use |
| `60NdII_levels_calib.txt`, `60NdII_transitions_calib.txt` | 910 kB, **687 MB** |

Provenance check after extracting: La II must read 472 levels / 17,743
transitions, Ce II 2,829 levels, Nd II 9,994 levels / 3,336,077 transitions.

Keeping the two zips costs 875 MB and saves re-downloading when another ion
is needed. Nd II parses in ~10 s at 1.9 GB peak RSS.

Some gitignored derived products must be regenerated once, because
`paper3/phase1_groups/compression.py` reads them:

```bash
cd paper3/phase0_reference
../../.venv/bin/python reference.py --ion laII   # ~6 s  -> reference_events_laII.npz
../../.venv/bin/python reference.py --ion ceII   # ~38 s -> reference_events_ceII.npz
```

Both print a Gate 0 check against the Paper II numbers; both should pass at
0.0 sigma.

## 3. SEDONA

Only `experiments/` needs it — everything in `paper2/` and `paper3/` is pure
Python. The build is the no-root WSL2 recipe from lab notebook §5, with
`Makefile.wsl` preserved in this directory (it lives in the pubsed clone,
which is outside this repo, so it was lost on the first machine move).

```bash
mkdir -p ~/personal && cd ~/personal

# 1. the code
git clone https://github.com/dnkasen/pubsed.git
cd pubsed && git config core.autocrlf input && git reset --hard HEAD
#   ^ CRLF trap: a global core.autocrlf=true checks the clone out with CRLF
#     line endings and every /bin/bash^M shebang fails.

# 2. deps, no root, via a static micromamba
mkdir -p ~/personal/sedona-deps && cd ~/personal/sedona-deps
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba create -y -p ~/personal/sedona-deps/env -c conda-forge \
    gsl hdf5 zlib openmpi openmpi-mpicxx

# 3. Lua 5.1.5 from source
curl -LO https://www.lua.org/ftp/lua-5.1.5.tar.gz
tar xzf lua-5.1.5.tar.gz && cd lua-5.1.5 && make generic

# 4. build
cp <this repo>/docs/sedona/Makefile.wsl ~/personal/pubsed/src/makefiles/
cd ~/personal/pubsed/src
SEDONA_HOME=~/personal/pubsed bash install.sh wsl     # -> src/sedona6.ex
```

Two things that are easy to lose:

- **MPI is not optional.** `src/main/sedona.h:10` hardcodes
  `#define MPI_PARALLEL 1` (the `#ifdef` guards around it are always true),
  so a serial `g++` build fails on `mpi.h`. Hence the mpicxx compile.
- **conda's `mpicxx` defaults to conda's own compiler**, which is not what we
  want here; `Makefile.wsl` forces the system g++ under it with `OMPI_CXX=g++`,
  and adds an `-Wl,-rpath` so the env's shared libs resolve at run time
  without `LD_LIBRARY_PATH`.

Running: **no shell setup is needed.** The drivers set `SEDONA_HOME`
themselves and exec `sedona6.ex` directly; MPI singleton-initialises to one
rank, and the `-Wl,-rpath` in `Makefile.wsl` resolves the conda libs without
`LD_LIBRARY_PATH`. Verified from a bare `PATH=/usr/bin:/bin`.

Put `~/personal/sedona-deps/env/bin` on `PATH` only if you want `mpirun` for
running SEDONA multi-rank by hand.

The experiment drivers locate the binary through the environment, defaulting
to `~/personal/pubsed`:

```bash
SEDONA_HOME=/somewhere/else/pubsed   # optional
SEDONA_EXE=/somewhere/else/sedona6.ex  # optional, overrides the above
```

(Before 2026-08-29 these were hardcoded to one machine's home directory,
which is why a clone elsewhere could not run any SEDONA experiment.)

Note the lab notebook §5 warning about passing a *trimmed* `env=` to
`subprocess`: `env={SEDONA_HOME, PATH}` alone breaks the OpenMPI runtime.
Pass the full environment plus your overrides.

## 4. LaTeX (only to rebuild the manuscript)

`docs/paper/manuscript.pdf` is committed, so this is needed only for
resubmission. conda-forge's `texlive-core` does **not** work — it ships
`mktexlsr` as a shell script where the TeX Live perl layer wants the
`mktexlsr.pl` module, so `fmtutil` cannot build `pdflatex.fmt`. Use TinyTeX,
which is the real TeX Live and designed for no-root:

```bash
curl -sSL https://yihui.org/tinytex/install-bin-unix.sh | sh
export PATH=~/.TinyTeX/bin/x86_64-linux:$PATH
tlmgr install amsmath amsfonts graphics booktabs natbib cm-super txfonts psnfss ec
cd docs/paper && make
```

`mnras.cls` and `mnras.bst` are committed, so no MNRAS package install is
needed.

## 5. Shell setup

None required for the experiments (see §3). The only thing worth adding to
`~/.bashrc`, and only if you rebuild the manuscript often:

```bash
export PATH=$HOME/.TinyTeX/bin/x86_64-linux:$PATH
```

## 6. Resource notes

Measured on the 2026-08-29 machine (24 GB RAM, 16 cores):

- full test suite 58 s
- Paper III La II reference 6 s; Ce II reference 38 s
- Paper III compression leg ~7 s per group count
- Paper II phase-1 forest driver at 2e5 packets: a few seconds per leg
- Nd II transition parse 10 s, 1.9 GB peak RSS

MC noise at `core_n_emit = 2e6` is ~1–2% per band flux; raise it before
chasing sub-percent effects.
