#!/usr/bin/env python3
"""Chapter 8 video: the macroatom as a Markov process on the level diagram,
internal jumps in blue and the de-activation in red."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.atom import five_level_atom                          # noqa: E402
from rtedu.macroatom import ToyMacroAtom                        # noqa: E402
from rtedu.visualization import animate_macroatom, save_animation   # noqa: E402


def walk_with_history(m, rng, level):
    hist = [level]
    while True:
        kind, x = m.step(rng, level)
        hist.append((kind, x))
        if kind == "deactivate":
            return hist
        level = x


def main(out_dir=rtedu.VIDEO_DIR, n=5, fps=2):
    a = five_level_atom(); m = ToyMacroAtom(a); rng = np.random.default_rng(rtedu.SEEDS["ch08"])
    walks = [walk_with_history(m, rng, 4) for _ in range(n)]
    return save_animation(animate_macroatom(a, walks, fps=fps), "macroatom", out_dir, fps=fps)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
