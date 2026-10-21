#!/usr/bin/env python3
"""Make the GIF animations of docs/images from finished cases.

    make_animations.py risingBubble <caseDir> [nx]
    make_animations.py zalesak <caseDir> [n]
    make_animations.py cavity <caseDir> [nx ny]

The cases must be written in ASCII on the structured (blockMesh) meshes of
validation/risingBubble, tutorials/ZalesakDisc/structured and
tests/regression/cases/thermocapillaryCavity. Needs numpy and matplotlib.
"""

import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402


def read_field(path, ncomp):
    s = open(path).read()
    s = s[s.index("internalField"):]
    m = re.match(r"internalField\s+uniform\s+(\S+?);", s)
    if m:
        raise ValueError("uniform field " + path)
    m = re.search(r"List<\w+>\s*(\d+)\s*\(", s)
    n = int(m.group(1))
    body = s[m.end():]
    if ncomp == 1:
        return np.array(body.split(None, n)[:n], dtype=float)
    vals = re.findall(r"\(([^()]*)\)", body[:body.index("\n)")])
    return np.array([[float(x) for x in v.split()] for v in vals[:n]])


def times(case):
    ts = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        f = os.path.join(case, d, "alpha.water")
        if os.path.exists(f) and "nonuniform" in open(f).read(2000):
            ts.append((t, d))
    return sorted(ts)


def save(anim, name, fps):
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    anim.save(out, writer=PillowWriter(fps=fps))
    print("wrote", out, "%.1f MB" % (os.path.getsize(out)/1e6))


def rising_bubble(case, nx=80):
    ny = 2*nx
    ts = times(case)
    fig, ax = plt.subplots(figsize=(2.6, 4.4), dpi=90)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.92)

    def frame(k):
        t, d = ts[k]
        a = read_field(os.path.join(case, d, "alpha.water"), 1).reshape(ny, nx)
        ax.clear()
        ax.imshow(a, origin="lower", extent=(0, 1, 0, 2), cmap="Blues",
                  vmin=-0.3, vmax=1.2, interpolation="bilinear")
        ax.contour(np.linspace(0.5/nx, 1 - 0.5/nx, nx),
                   np.linspace(0.5/nx, 2 - 0.5/nx, ny), a, [0.5],
                   colors="k", linewidths=1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Rising bubble  t = %.2f" % t, fontsize=9)

    save(FuncAnimation(fig, frame, frames=len(ts)), "risingBubble.gif", 12)


def zalesak(case, n=100):
    ts = times(case)
    a0 = None
    fig, ax = plt.subplots(figsize=(3.4, 3.6), dpi=90)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.92)
    x = np.linspace(0.5/n, 1 - 0.5/n, n)

    def frame(k):
        nonlocal a0
        t, d = ts[k]
        a = read_field(os.path.join(case, d, "alpha.water"), 1).reshape(n, n)
        if a0 is None:
            a0 = a
        ax.clear()
        ax.imshow(a, origin="lower", extent=(0, 1, 0, 1), cmap="Blues",
                  vmin=-0.3, vmax=1.2, interpolation="nearest")
        ax.contour(x, x, a0, [0.5], colors="0.4", linewidths=0.8,
                   linestyles="--")
        ax.contour(x, x, a, [0.5], colors="k", linewidths=1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Zalesak disc  %d/%d revolution" % (round(t/(2*np.pi)*32), 32),
                     fontsize=9)

    save(FuncAnimation(fig, frame, frames=len(ts)), "zalesakDisc.gif", 8)


def cavity(case, nx=50, ny=20):
    ts = times(case)[1::2]           # every second frame (file size)
    W, H = 5.0, 2.0
    x = np.linspace(W/nx/2, W - W/nx/2, nx)
    y = np.linspace(H/ny/2, H - H/ny/2, ny)
    fig, ax = plt.subplots(figsize=(6.4, 2.9), dpi=80)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.88)

    def frame(k):
        t, d = ts[k]
        a = read_field(os.path.join(case, d, "alpha.water"), 1).reshape(ny, nx)
        T = read_field(os.path.join(case, d, "T"), 1).reshape(ny, nx)
        U = read_field(os.path.join(case, d, "U"), 3).reshape(ny, nx, 3)
        ax.clear()
        ax.imshow(T, origin="lower", extent=(0, W, 0, H), cmap="coolwarm",
                  vmin=0, vmax=1, interpolation="bilinear", alpha=0.85)
        # 2D case in the x-y or in the x-z plane
        v = U[:, :, 1] if np.abs(U[:, :, 1]).max() > np.abs(U[:, :, 2]).max() else U[:, :, 2]
        ax.streamplot(x, y, U[:, :, 0], v, color="0.25",
                      density=1.2, linewidth=0.6, arrowsize=0.6)
        ax.contour(x, y, a, [0.5], colors="k", linewidths=1.5)
        ax.set_xlim(0, W)
        ax.set_ylim(0, H)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Thermocapillary cavity (hot left, cold right)  t = %.1f" % t,
                     fontsize=9)

    save(FuncAnimation(fig, frame, frames=len(ts)), "thermocapillaryCavity.gif", 8)


if __name__ == "__main__":
    kind, case = sys.argv[1], sys.argv[2]
    args = [int(v) for v in sys.argv[3:]]
    {"risingBubble": rising_bubble, "zalesak": zalesak, "cavity": cavity}[kind](case, *args)
