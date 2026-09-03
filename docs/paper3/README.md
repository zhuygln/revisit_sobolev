# Paper III — the manuscript

*Coarse-grained line opacity leaves a detectable chromatic signature in
kilonovae that ejecta parameters cannot mimic.* Nature Astronomy Article
format.

This is the write-up half of Paper III. The other half — the code, the run
outputs and the frozen record every number here is generated from — is
[`paper3/`](../../paper3/) at the repository root. Same paper, split the way
the repo splits every paper: `paperN/` is the campaign, `docs/paperN/` is the
manuscript.

| | |
|---|---|
| `manuscript.tex` / `.pdf` | the Article |
| `si.tex` / `si.pdf`, `si_tab_*.tex` | Supplementary Information |
| `numbers.tex`, `tab_*.tex` | generated — do not edit |
| `latex_tables.py` | numbers.tex and every table, from `paper3/FROZEN.json` |
| `display_items.py` | every figure, from the JSONs `paper3/freeze.py` records |
| `check_structure.py` | bans literal result numbers in the prose of both |
| `cover_letter.md`, `literature_notes.md`, `references.bib` | submission material |
| `figures/` | the manuscript display items (fig1-4, edfig1-4) — the working figures of the results report are in `paper3/figures/` |

```bash
make tables    # numbers.tex and every table fragment, from paper3/FROZEN.json
make figures   # every display item, from the frozen derived JSONs
make           # manuscript.pdf + si.pdf, then the structure check
make check     # structure check only
```

Nothing in this directory is transcribed by hand: rebuild from
`paper3/FROZEN.json`, whose state the tag `paper3-freeze` marks. The
affiliation and repository URL in `manuscript.tex` are still placeholders.
