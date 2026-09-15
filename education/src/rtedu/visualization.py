"""Figures and animations for the book, written reproducibly.

save_fig: PNG + PDF with the metadata dates stripped (the docs/paper4 pattern).
save_animation: GIF via pillow (byte-stable, fixed fps) and an embeddable
HTML/JS fragment via matplotlib's to_jshtml (frames embedded; the fragment,
not a full document, so a chapter can include it directly). MP4 only when an
ffmpeg writer is available (it is not on the build machine).
"""
import hashlib
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from . import FIG_DIR, VIDEO_DIR

OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "red": "#D55E00",
      "purple": "#CC79A7", "sky": "#56B4E9", "yellow": "#F0E442", "black": "#000000"}
GROUP_COLOUR = {"B": OI["blue"], "V": OI["green"], "IR": OI["red"]}


def save_fig(fig, name, out_dir=FIG_DIR, dpi=150):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None}),
                      ("png", {"Software": None})):
        p = out_dir / f"{name}.{ext}"
        fig.savefig(p, metadata=meta, bbox_inches="tight", pad_inches=0.02, dpi=dpi)
        written.append(p)
    plt.close(fig)
    return written


def save_animation(anim, name, out_dir=VIDEO_DIR, fps=8, html_dpi=60):
    """GIF (pillow) + HTML fragment (to_jshtml); MP4 if ffmpeg exists."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    gif = out_dir / f"{name}.gif"
    anim.save(gif, writer=animation.PillowWriter(fps=fps))
    written.append(gif)
    html = out_dir / f"{name}.html"
    dpi0 = anim._fig.get_dpi(); anim._fig.set_dpi(html_dpi)        # smaller embedded frames for the page
    html.write_text(anim.to_jshtml(fps=fps, embed_frames=True, default_mode="loop"))
    anim._fig.set_dpi(dpi0)
    written.append(html)
    if animation.writers.is_available("ffmpeg"):
        mp4 = out_dir / f"{name}.mp4"
        anim.save(mp4, writer=animation.FFMpegWriter(fps=fps))
        written.append(mp4)
    plt.close(anim._fig)
    return written


def animation_frame_hashes(html_path):
    """The hashes of the frames embedded in a to_jshtml fragment: two
    fragments with the same frames are the same animation even if the
    generated element ids differ (PI amendment 6)."""
    text = Path(html_path).read_text()
    frames = re.findall(r'"data:image/png;base64,([^"]+)"', text)
    return [hashlib.sha256(f.encode()).hexdigest()[:16] for f in frames]


def chain_diagram():
    """Chapter 0: the kilonova chain with the RT box highlighted, and the
    RT loop underneath."""
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11, 4.6), gridspec_kw=dict(height_ratios=[1, 0.8]))
    chain = ["nuclear\nphysics", "abundances,\nheating", "ejecta\nstructure", "atomic\nphysics", "radiative\ntransfer", "SED", "light\ncurves"]
    x = np.linspace(0.06, 0.94, len(chain))
    for i, (xi, label) in enumerate(zip(x, chain)):
        rt = label.startswith("radiative")
        box = FancyBboxPatch((xi - 0.055, 0.3), 0.11, 0.4, boxstyle="round,pad=0.01",
                             fc=OI["orange"] if rt else "#f2f2f2", ec=OI["black"], lw=1.5 if rt else 0.8)
        ax.add_patch(box); ax.text(xi, 0.5, label, ha="center", va="center", fontsize=9, fontweight="bold" if rt else None)
        if i < len(chain) - 1:
            ax.add_patch(FancyArrowPatch((xi + 0.057, 0.5), (x[i + 1] - 0.057, 0.5), arrowstyle="-|>", mutation_scale=12, lw=1))
    ax.text(x[4], 0.86, "this book (Paper IV, Paper B, the emulator)", ha="center", fontsize=9, color=OI["red"])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    loop = ["propagate", "interact", "redistribute", "propagate again"]
    xl = np.linspace(0.15, 0.85, len(loop))
    for i, (xi, label) in enumerate(zip(xl, loop)):
        part = "E1: chapters 1-4" if label.startswith("propagate") or label == "interact" else "E2: chapters 5-8"
        box = FancyBboxPatch((xi - 0.09, 0.35), 0.18, 0.35, boxstyle="round,pad=0.01", fc="#e8f1fa" if "E1" in part else "#fbeee6", ec=OI["black"], lw=0.8)
        ax2.add_patch(box); ax2.text(xi, 0.525, label, ha="center", va="center", fontsize=9)
        ax2.text(xi, 0.18, part, ha="center", va="center", fontsize=7.5, color=OI["blue"] if "E1" in part else OI["red"])
        if i < len(loop) - 1:
            ax2.add_patch(FancyArrowPatch((xi + 0.092, 0.525), (xl[i + 1] - 0.092, 0.525), arrowstyle="-|>", mutation_scale=12, lw=1))
    ax2.add_patch(FancyArrowPatch((xl[-1], 0.35), (xl[0], 0.35), arrowstyle="-|>", mutation_scale=12, lw=1, connectionstyle="arc3,rad=0.25", ls="--"))
    ax2.text(0.5, 0.05, "the RT loop: E1 answers how a photon finds an atom; E2 what happens to its energy", ha="center", fontsize=8.5)
    ax2.set_xlim(0, 1); ax2.set_ylim(0, 1); ax2.axis("off")
    return fig


def animate_packet_walk(paths, radius, fps=8, title="packets in a disc of radius R mean free paths"):
    """Several packets' random walks drawn step by step until each escapes."""
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.add_patch(plt.Circle((0, 0), radius, fc="#fff3e0", ec=OI["black"], lw=1))
    ax.set_xlim(-1.3 * radius, 1.3 * radius); ax.set_ylim(-1.3 * radius, 1.3 * radius); ax.set_aspect("equal")
    ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    lines = [ax.plot([], [], "-", lw=1.2, color=c)[0] for c, _ in zip(list(OI.values()), paths)]
    dots = [ax.plot([], [], "o", ms=5, color=l.get_color())[0] for l in lines]
    n_frames = max(len(p) for p in paths)

    def draw(k):
        for p, l, d in zip(paths, lines, dots):
            j = min(k + 1, len(p))
            l.set_data(p[:j, 0], p[:j, 1]); d.set_data([p[j - 1, 0]], [p[j - 1, 1]])
        return lines + dots
    return animation.FuncAnimation(fig, draw, frames=n_frames, interval=1000 / fps, blit=True)


def animate_resonance(nu_lab, nu_lines, tau_lines, r_out, t, crossings, fps=8, n_frames=60, c=None):
    """Chapter 4: the packet's position on top, its comoving frequency
    against radius underneath with the lines as horizontal markers, and a
    flash where a resonance is crossed (a filled marker where it interacted)."""
    from . import C
    from .sobolev import comoving_frequency
    c = C if c is None else c
    r = np.linspace(0.0, r_out, n_frames)
    nu_c = comoving_frequency(nu_lab, r, t)
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7, 5.4), gridspec_kw=dict(height_ratios=[1, 2.2]))
    ax.set_xlim(0, r_out); ax.set_ylim(-1, 1); ax.set_yticks([]); ax.set_xlabel("r")
    ax.axvspan(0, r_out, color="#fff3e0"); ax.set_title("the packet moves outward through homologous flow v = r/t", fontsize=9)
    dot, = ax.plot([], [], "o", ms=9, color=OI["blue"])
    ax2.plot(r, nu_c / nu_lab, color=OI["blue"], lw=1.5, label=r"$\nu_{\rm com}(r)/\nu_{\rm lab}$")
    for k, (nl, tl) in enumerate(zip(nu_lines, tau_lines)):
        ax2.axhline(nl / nu_lab, color="grey", lw=0.6 + 0.8 * min(tl, 3) / 3, alpha=0.7)
    ax2.set_xlim(0, r_out); ax2.set_xlabel("r"); ax2.set_ylabel(r"$\nu/\nu_{\rm lab}$")
    ax2.set_ylim(min(nu_c.min() / nu_lab, min(nu_lines) / nu_lab) * 0.995, 1.005)
    ax2.legend(loc="upper right", fontsize=8); ax2.set_title("its comoving frequency sweeps down through the lines (grey: tau-weighted)", fontsize=9)
    trail, = ax2.plot([], [], "o", ms=7, color=OI["blue"])
    flashes = [ax2.plot([], [], "o", ms=14, mfc="none", mec=OI["red"], mew=2)[0] for _ in crossings]
    hits = [ax2.plot([], [], "*", ms=16, color=OI["red"])[0] for _ in crossings]

    def draw(k):
        dot.set_data([r[k]], [0.0]); trail.set_data([r[k]], [nu_c[k] / nu_lab])
        for cr, fl, ht in zip(crossings, flashes, hits):
            if r[k] >= cr["r_res"]:
                (ht if cr["interacted"] else fl).set_data([cr["r_res"]], [nu_lines[cr["line"]] / nu_lab])
            else:
                fl.set_data([], []); ht.set_data([], [])
        return [dot, trail] + flashes + hits
    return animation.FuncAnimation(fig, draw, frames=n_frames, interval=1000 / fps, blit=True)
