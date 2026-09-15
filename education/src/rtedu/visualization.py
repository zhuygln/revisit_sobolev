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
    text = anim.to_jshtml(fps=fps, embed_frames=True, default_mode="loop")
    anim._fig.set_dpi(dpi0)
    # matplotlib names the player's elements after a fresh uuid; replace it
    # by a name-derived id so the fragment is byte-stable across builds
    ids = set(re.findall(r"[0-9a-f]{32}", text))
    fixed = hashlib.sha256(name.encode()).hexdigest()[:32]
    for u in ids:
        text = text.replace(u, fixed)
    html.write_text(text)
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


def level_diagram(ax, atom, label_lines=True):
    """The energy-level diagram of a ToyAtom with every allowed transition
    as a faint arrow; returns the arrow artists keyed by line index."""
    E = atom.E
    for i, e in enumerate(E):
        ax.hlines(e, 0.1, 0.9, color=OI["black"], lw=1.5)
        ax.text(0.92, e, f"$E_{i}$", va="center", fontsize=9)
    arrows = {}
    xs = np.linspace(0.2, 0.8, atom.n_lines)
    for k in range(atom.n_lines):
        u, l = atom.upper[k], atom.lower[k]
        arrows[k] = FancyArrowPatch((xs[k], E[u]), (xs[k], E[l]), arrowstyle="-|>", mutation_scale=9, lw=0.8, color="grey", alpha=0.5)
        ax.add_patch(arrows[k])
        if label_lines:
            ax.text(xs[k], (E[u] + E[l]) / 2, f"{1e7 * atom.lam_cm[k]:.0f}", fontsize=6, ha="left", color="grey")
    ax.set_xlim(0, 1); ax.set_ylim(-0.05 * E[-1], 1.12 * E[-1]); ax.set_xticks([]); ax.set_ylabel("energy [cm$^{-1}$]")
    return arrows


def animate_cascades(atom, cascades, fps=2, title="explicit cascades from the top level"):
    """Chapter 6: one cascade after another on the level diagram; each step
    highlights the current level, draws the transition and names the
    emitted photon's wavelength."""
    fig, ax = plt.subplots(figsize=(5.2, 5))
    arrows = level_diagram(ax, atom); ax.set_title(title, fontsize=9)
    frames = []
    for c in cascades:
        level = int(atom.upper[c[0]])
        frames.append(("activate", level, None))
        for k in c:
            frames.append(("emit", int(atom.upper[k]), int(k)))
        frames.append(("done", 0, None))
    dot, = ax.plot([], [], "o", ms=12, color=OI["orange"])
    txt = ax.text(0.5, 1.06 * atom.E[-1], "", ha="center", fontsize=9, color=OI["red"])
    hi = [FancyArrowPatch((0, 0), (0, 0), arrowstyle="-|>", mutation_scale=14, lw=2.5, color=OI["red"]) for _ in range(1)]
    for h in hi: ax.add_patch(h); h.set_visible(False)

    def draw(i):
        kind, level, k = frames[i]
        dot.set_data([0.05], [atom.E[level]])
        if kind == "activate":
            txt.set_text(f"absorbed: activated at $E_{level}$"); hi[0].set_visible(False)
        elif kind == "emit":
            u, l = atom.upper[k], atom.lower[k]; x = np.linspace(0.2, 0.8, atom.n_lines)[k]
            hi[0].set_positions((x, atom.E[u]), (x, atom.E[l])); hi[0].set_visible(True)
            txt.set_text(f"emit {1e7 * atom.lam_cm[k]:.0f} nm  ($E_{u} \\to E_{l}$)")
        else:
            txt.set_text("back at the ground: the cascade is over"); hi[0].set_visible(False)
        return [dot, txt] + hi
    return animation.FuncAnimation(fig, draw, frames=len(frames), interval=1000 / fps, blit=True)


def animate_macroatom(atom, walks, fps=2):
    """Chapter 8: the macroatom as a Markov process on the level diagram:
    internal jumps in blue (the packet's energy stays inside), the
    de-activation in red (one packet leaves with all the energy)."""
    fig, ax = plt.subplots(figsize=(5.2, 5))
    level_diagram(ax, atom, label_lines=False); ax.set_title("the macroatom: jump (blue) or de-activate (red)", fontsize=9)
    frames = []
    for w in walks:                                    # w: list of ('jump', level) / ('deactivate', line), starting level first
        frames.append(("activate", w[0], None))
        cur = w[0]
        for kind, x in w[1:]:
            frames.append((kind, cur, x)); cur = x if kind == "jump" else cur
        frames.append(("done", 0, None))
    dot, = ax.plot([], [], "o", ms=12, color=OI["orange"])
    txt = ax.text(0.5, 1.06 * atom.E[-1], "", ha="center", fontsize=9)
    arr = FancyArrowPatch((0, 0), (0, 0), arrowstyle="-|>", mutation_scale=14, lw=2.5, color=OI["blue"]); ax.add_patch(arr); arr.set_visible(False)

    def draw(i):
        kind, level, x = frames[i]
        dot.set_data([0.05], [atom.E[level]]); arr.set_visible(False)
        if kind == "activate":
            txt.set_text(f"activated at $E_{level}$ with energy $E$"); txt.set_color(OI["black"])
        elif kind == "jump":
            arr.set_positions((0.5, atom.E[level]), (0.5, atom.E[x])); arr.set_color(OI["blue"]); arr.set_visible(True)
            txt.set_text(f"internal jump $E_{level} \\to E_{x}$: the energy stays in the atom"); txt.set_color(OI["blue"])
        elif kind == "deactivate":
            u, l = atom.upper[x], atom.lower[x]
            arr.set_positions((0.5, atom.E[u]), (0.5, atom.E[l])); arr.set_color(OI["red"]); arr.set_visible(True)
            txt.set_text(f"de-activate in {1e7 * atom.lam_cm[x]:.0f} nm: one packet leaves with all of $E$"); txt.set_color(OI["red"])
        else:
            txt.set_text(""); 
        return [dot, txt, arr]
    return animation.FuncAnimation(fig, draw, frames=len(frames), interval=1000 / fps, blit=True)


def animate_matrix_vs_T(Ts, Rs, fps=1.5, title="R(T): the redistribution matrix as the state changes"):
    """Chapter 11: R heat maps morphing with temperature."""
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    im = ax.imshow(np.asarray(Rs[0]), cmap="viridis", vmin=0, vmax=1); ax.set_xlabel("emitted group j"); ax.set_ylabel("absorbed group i")
    plt.colorbar(im, ax=ax, fraction=0.046); ttl = ax.set_title("", fontsize=9)
    n = np.asarray(Rs[0]).shape[0]
    texts = [[ax.text(j, i, "", ha="center", va="center", fontsize=8) for j in range(n)] for i in range(n)]

    def draw(k):
        R = np.asarray(Rs[k]); im.set_data(R); ttl.set_text(f"{title}\nT = {Ts[k]:.0f} K")
        for i in range(n):
            for j in range(n):
                texts[i][j].set_text(f"{R[i, j]:.2f}"); texts[i][j].set_color("w" if R[i, j] < 0.5 else "k")
        return [im, ttl] + [t for row in texts for t in row]
    return animation.FuncAnimation(fig, draw, frames=len(Ts), interval=1000 / fps, blit=False)


def animate_compare_models(nm, histories, labels, fps=2):
    """Chapter 9: three packets' frequency histories side by side, one per
    redistribution model: the line each packet sits in after each
    interaction, until it escapes."""
    fig, axes = plt.subplots(1, len(histories), figsize=(4 * len(histories), 3.6), sharey=True)
    n_frames = max(len(h) for h in histories)
    lines = []
    for ax, h, lab in zip(axes, histories, labels):
        ax.set_title(lab, fontsize=9); ax.set_xlabel("interaction"); ax.set_xlim(-0.5, n_frames - 0.5); ax.set_yscale("log")
        ax.set_ylim(min(nm) * 0.9, max(nm) * 1.1); ax.set_yticks(sorted(nm)); ax.set_yticklabels([f"{v:.0f}" for v in sorted(nm)], fontsize=6)
        for v in nm: ax.axhline(v, color="lightgrey", lw=0.5)
        lines.append(ax.plot([], [], "o-", color=OI["blue"])[0])
    axes[0].set_ylabel("line the packet is in [nm]")

    def draw(k):
        for l, h in zip(lines, histories):
            j = min(k + 1, len(h)); l.set_data(np.arange(j), [nm[x] for x in h[:j]])
        return lines
    return animation.FuncAnimation(fig, draw, frames=n_frames, interval=1000 / fps, blit=True)


def animate_full_packet(trace, nm, r_out_of_t, fps=2):
    """Chapter 13: one packet from creation to escape: its projected path
    in the expanding sphere (left) and the line it sits in against its own
    time (right)."""
    from . import DAY
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.2))
    xs = np.array([p[0] for p in trace]); ts = np.array([p[1] for p in trace]); ls = [p[2] for p in trace]
    R_end = r_out_of_t(ts[-1])
    ax.set_xlim(-1.1 * R_end, 1.1 * R_end); ax.set_ylim(-1.1 * R_end, 1.1 * R_end); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    circ = plt.Circle((0, 0), r_out_of_t(ts[0]), fc="#fff3e0", ec=OI["black"], lw=1); ax.add_patch(circ)
    path, = ax.plot([], [], "-", color=OI["blue"], lw=1.2); dot, = ax.plot([], [], "o", ms=7, color=OI["blue"])
    ax.set_title("the packet's path (x-y projection); the sphere grows with time", fontsize=8)
    ax2.set_yscale("log"); ax2.set_ylim(min(nm) * 0.9, max(nm) * 1.1); ax2.set_yticks(sorted(nm)); ax2.set_yticklabels([f"{v:.0f}" for v in sorted(nm)], fontsize=6)
    ax2.set_xlim(ts[0] / DAY, ts[-1] / DAY * 1.02); ax2.set_xlabel("the packet's own time [d]"); ax2.set_ylabel("line [nm]")
    steps, = ax2.plot([], [], "o-", color=OI["red"]); ax2.set_title("the line it was last emitted in, against its clock", fontsize=8)
    lam = [nm[l] if l >= 0 else np.nan for l in ls]
    cur = 0
    for i, l in enumerate(ls):
        if l >= 0: cur = nm[l]
        if np.isnan(lam[i]): lam[i] = cur if i > 0 else np.nan

    def draw(k):
        j = k + 1
        circ.set_radius(r_out_of_t(ts[k])); path.set_data(xs[:j, 0], xs[:j, 1]); dot.set_data([xs[k, 0]], [xs[k, 1]])
        tt = ts[:j] / DAY; ll = [v for v in lam[:j]]
        keep = [i for i in range(j) if not np.isnan(ll[i])]
        steps.set_data([tt[i] for i in keep], [ll[i] for i in keep])
        return [circ, path, dot, steps]
    return animation.FuncAnimation(fig, draw, frames=len(trace), interval=1000 / fps, blit=False)
