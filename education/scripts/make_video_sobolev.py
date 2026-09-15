#!/usr/bin/env python3
"""Chapter 4 video: one packet moving outward through homologous flow (top),
its comoving frequency sweeping down through the toy line list (bottom), a
ring where a resonance is crossed and a star where it interacted. The line
list is rtedu.sobolev.toy_line_list, shared with notebooks/04_sobolev.ipynb."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import rtedu                                                   # noqa: E402
from rtedu.sobolev import resonance_crossings                  # noqa: E402
from rtedu.visualization import animate_resonance, save_animation   # noqa: E402


def line_list():
    from rtedu.sobolev import toy_line_list
    L = toy_line_list()
    return L["nu_lab"], L["nu_lines"], L["tau"], L["r_out"], L["t"]


def main(out_dir=rtedu.VIDEO_DIR, fps=8):
    nu_lab, nu_lines, tau, r_out, t = line_list()
    rng = np.random.default_rng(rtedu.SEEDS["ch04"] + 7)
    crossings = resonance_crossings(nu_lab, nu_lines, tau, 0.0, r_out, t, rng=rng)
    anim = animate_resonance(nu_lab, nu_lines, tau, r_out, t, crossings, fps=fps, n_frames=60)
    return save_animation(anim, "sobolev_sweep", out_dir, fps=fps)


if __name__ == "__main__":
    for p in main():
        print("wrote", p)
