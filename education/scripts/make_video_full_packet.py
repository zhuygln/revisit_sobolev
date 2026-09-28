#!/usr/bin/env python3
"""Chapter 13 video: one packet from creation to escape in the expanding
sphere, under the macroatom."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.atom import five_level_atom                          # noqa: E402
from rtedu.matrix import MacroatomRedistribution                # noqa: E402
from rtedu.transport3d import ToyEjecta, launch, trace_one      # noqa: E402
from rtedu.visualization import animate_full_packet, save_animation   # noqa: E402


def main(out_dir=rtedu.VIDEO_DIR, fps=2):
    a = five_level_atom(); ej = ToyEjecta(a); rng = np.random.default_rng(rtedu.SEEDS["ch13"] + 5)
    cache = {}
    def at(t):
        key = round(t / rtedu.DAY, 2)
        if key not in cache:
            cache[key] = MacroatomRedistribution(a, ej.tau(t))
        return cache[key]
    x, d, nu, t = launch(rng, ej, 40, 0, 3 * rtedu.DAY)
    best = None
    for i in range(40):                                   # the packet with the most interactions among forty
        tr = trace_one(np.random.default_rng(1000 + i), ej, x[i], d[i], nu[i], t[i], at)
        if best is None or len(tr) > len(best):
            best = tr
    return save_animation(animate_full_packet(best, 1e7 * a.lam_cm, ej.r_out, fps=fps), "full_packet", out_dir, fps=fps)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
