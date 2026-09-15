#!/usr/bin/env python3
"""Chapter 6 video: a few explicit cascades from the top level of the
five-level atom, one photon per downward step."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.atom import five_level_atom                          # noqa: E402
from rtedu.visualization import animate_cascades, save_animation   # noqa: E402


def main(out_dir=rtedu.VIDEO_DIR, n=4, fps=2):
    a = five_level_atom(); rng = np.random.default_rng(rtedu.SEEDS["ch06"])
    cascades = [a.cascade(rng, 4) for _ in range(n)]
    return save_animation(animate_cascades(a, cascades, fps=fps), "cascade", out_dir, fps=fps)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
