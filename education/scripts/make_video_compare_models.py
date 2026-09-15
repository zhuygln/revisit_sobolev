#!/usr/bin/env python3
"""Chapter 9 video: one packet under epsilon, under R and under the
macroatom, the line it sits in after each interaction, side by side."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.atom import five_level_atom                          # noqa: E402
from rtedu.matrix import MacroatomRedistribution, build_R, group_of_line, MatrixRedistribution   # noqa: E402
from rtedu.redistribution import EpsilonRedistribution          # noqa: E402
from rtedu.transport import run, sweep                          # noqa: E402
from rtedu.visualization import animate_compare_models, save_animation   # noqa: E402


def main(out_dir=rtedu.VIDEO_DIR, fps=2):
    a = five_level_atom(); t = 2 * rtedu.DAY; r_out = 0.2 * rtedu.C * t; T, n_tot = 4000.0, 30.0
    tau = a.line_list(T, n_tot, t); emis = a.thermal_emissivity(T, n_tot); nm = 1e7 * a.lam_cm
    macro = MacroatomRedistribution(a, tau)
    run(np.random.default_rng(rtedu.SEEDS["ch09"] + 200), a.nu[0] * 1.001, 1500, a.nu, tau, r_out, t, macro)
    g, ng = group_of_line(a.nu, 4); R = build_R(macro.events, g, ng)
    models = [("epsilon = 0.5", EpsilonRedistribution(0.5, emis)), ("R, 4 groups", MatrixRedistribution(R, g, emis)), ("macroatom", MacroatomRedistribution(a, tau))]
    hist = []
    for k, (lab, m) in enumerate(models):
        rng = np.random.default_rng(rtedu.SEEDS["ch09"] + 300 + k)
        for _ in range(50):                                 # the first packet with at least three interactions
            nu, n_int, h = sweep(rng, a.nu[0] * 1.001, a.nu, tau, r_out, t, m)
            if len(h) >= 3:
                break
        hist.append([0] + [j for (_, j) in h])
    return save_animation(animate_compare_models(nm, hist, [l for l, _ in models], fps=fps), "compare_models", out_dir, fps=fps)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
