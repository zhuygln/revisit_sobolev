#!/usr/bin/env python3
"""Chapter 0/2 video: a few packets random-walking out of a disc of radius R
mean free paths (exponential step lengths, isotropic directions). GIF +
HTML fragment into education/videos/; deterministic (seeded)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.packets import walk_2d                              # noqa: E402
from rtedu.visualization import animate_packet_walk, save_animation   # noqa: E402


def main(out_dir=rtedu.VIDEO_DIR, radius=4.0, n_packets=5, fps=8):
    rng = np.random.default_rng(rtedu.SEEDS["ch00"])
    paths = [walk_2d(rng, radius, max_steps=60) for _ in range(n_packets)]
    anim = animate_packet_walk(paths, radius, fps=fps)
    return save_animation(anim, "mc_walk", out_dir, fps=fps)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
