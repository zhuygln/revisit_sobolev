"""Compatibility shim: the Paper IV figures live in display_items.py (PDF +
PNG into docs/paper4/figures/); this writes the PNGs to docs/figures/paper4/
for the results report."""
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import display_items  # noqa: E402

OUT = ROOT / "docs/figures/paper4"

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for p in display_items.main(HERE / "figures"):
        if p.suffix == ".png":
            shutil.copy(p, OUT / p.name)
    print("figures written to", OUT)
