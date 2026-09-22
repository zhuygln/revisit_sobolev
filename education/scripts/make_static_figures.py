#!/usr/bin/env python3
"""Static figures that belong to no notebook (chapter 0's chain diagram)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.visualization import chain_diagram, save_fig        # noqa: E402


def main(out_dir=rtedu.FIG_DIR):
    return save_fig(chain_diagram(), "ch00_chain", out_dir)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
