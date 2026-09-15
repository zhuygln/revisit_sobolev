#!/usr/bin/env python3
"""Chapter 11 video: R(T) at a sequence of temperatures."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.atom import five_level_atom                          # noqa: E402
from rtedu.matrix import MacroatomRedistribution, build_R, group_of_line   # noqa: E402
from rtedu.transport import run                                 # noqa: E402
from rtedu.visualization import animate_matrix_vs_T, save_animation   # noqa: E402


def main(out_dir=rtedu.VIDEO_DIR, fps=1.5):
    a = five_level_atom(); t = 2 * rtedu.DAY; r_out = 0.2 * rtedu.C * t; n_tot = 30.0
    g, ng = group_of_line(a.nu, 4)
    Ts = [2000.0, 2500.0, 3000.0, 4000.0, 5000.0, 7000.0, 10000.0]; Rs = []
    for i, T in enumerate(Ts):
        tau = a.line_list(T, n_tot, t); m = MacroatomRedistribution(a, tau)
        run(np.random.default_rng(rtedu.SEEDS["ch11"] + 100 + i), a.nu[0] * 1.001, 1500, a.nu, tau, r_out, t, m)
        Rs.append(build_R(m.events, g, ng))
    return save_animation(animate_matrix_vs_T(Ts, Rs, fps=fps), "matrix_T", out_dir, fps=fps)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
