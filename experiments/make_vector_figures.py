"""Re-emit every manuscript figure as vector PDF, without touching generators.

MNRAS prefers vector graphics for line art, and every figure here is line art.
The generators all write PNG at dpi=200; at the width they are placed that is
~310 dpi effective, below the journal's 400 dpi bar for raster.

Rather than edit thirteen scripts (and risk changing what they plot while
correcting how they are saved), this wrapper patches Figure.savefig so that
every PNG write is mirrored to a .pdf sibling. The generators are executed
unmodified, so the vector figure is by construction the same figure.

Usage:  python experiments/make_vector_figures.py [name ...]
"""

import runpy
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "docs/figures"

# (script, working directory) -- several read run outputs by relative path.
GENERATORS = [
    ("experiments/minimal_1line/compare.py", "experiments/minimal_1line"),
    ("experiments/line_ladder/compare.py", "experiments/line_ladder"),
    ("experiments/laII_forest/compare.py", "experiments/laII_forest"),
    ("experiments/laII_forest/fig7.py", "experiments/laII_forest"),
    ("experiments/laII_forest/fig8.py", "experiments/laII_forest"),
    ("experiments/multiion/compare.py", "experiments/multiion"),
    ("experiments/sobolev_proper/recompute_all.py", "experiments/sobolev_proper"),
    ("experiments/boundary/run.py", "experiments/boundary"),
    ("experiments/breadth/fig11.py", "experiments/breadth"),
    ("experiments/binwidth/fig13.py", "experiments/binwidth"),
]

_orig = Figure.savefig
written = []


def savefig(self, fname, *a, **kw):
    out = _orig(self, fname, *a, **kw)
    p = Path(str(fname))
    if p.suffix.lower() == ".png":
        kw.pop("dpi", None)
        pdf = p.with_suffix(".pdf")
        _orig(self, pdf, *a, **kw)
        written.append(pdf)
        # Keep docs/figures in step with outputs/.
        if pdf.parent.name != "figures":
            twin = FIGDIR / pdf.name
            if (FIGDIR / p.name).exists() or twin.exists():
                _orig(self, twin, *a, **kw)
                written.append(twin)
    return out


Figure.savefig = savefig

wanted = sys.argv[1:]
for script, cwd in GENERATORS:
    if wanted and not any(w in script for w in wanted):
        continue
    print(f"\n=== {script} ===", flush=True)
    here = Path.cwd()
    try:
        import os
        os.chdir(ROOT / cwd)
        sys.argv = [script]
        runpy.run_path(str(ROOT / script), run_name="__main__")
    except SystemExit:
        pass
    except Exception:
        print(f"  FAILED: {script}", flush=True)
        traceback.print_exc()
    finally:
        import os
        os.chdir(here)

print(f"\nwrote {len(written)} vector files")
for w in sorted(set(str(x) for x in written)):
    print("  ", w)
