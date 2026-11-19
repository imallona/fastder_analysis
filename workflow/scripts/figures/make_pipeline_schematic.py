"""fastder pipeline schematic (Main Figure 1, panel A).

Five steps left to right: inputs, normalization and averaging, regions above
the coverage threshold, stitching across junctions, GTF output.
Output: fig_pipeline_schematic.pdf, .svg and .png in the working directory.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "text.color": "#222222",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

COVERAGE = "#222222"
REGION = "#e7a1c0"
JUNCTION = "#7fb8dc"
THRESHOLD = "#d6609a"
TEXT_SIZE = 7

# Three exons on a 0..1 axis; the same layout in every step.
EXONS = [(0.08, 0.22), (0.38, 0.50), (0.74, 0.92)]
STEP_WIDTH = 1.0
GAP = 0.45


def profile(x, scale=1.0, noise=0.0, rng=None):
    y = np.zeros_like(x)
    for start, end in EXONS:
        inside = (x >= start) & (x <= end)
        y[inside] = scale * (0.75 + 0.25 * np.sin((x[inside] - start) * 40) ** 2)
    if noise:
        y += noise * rng.random(x.size)
    return y


def coverage(ax, x0, y0, height, scale, noise, rng):
    x = np.linspace(0, 1, 400)
    y = profile(x, scale, noise, rng)
    ax.fill_between(x0 + x * STEP_WIDTH, y0, y0 + y * height, color=COVERAGE, lw=0)
    ax.plot([x0, x0 + STEP_WIDTH], [y0, y0], color=COVERAGE, lw=0.6)


def regions(ax, x0, y0, height=0.12, color=REGION):
    for start, end in EXONS:
        ax.add_patch(Rectangle((x0 + start * STEP_WIDTH, y0), (end - start) * STEP_WIDTH,
                               height, facecolor=color, edgecolor="none"))


def arcs(ax, x0, y0, rise=0.14):
    for (_, left_end), (right_start, _) in zip(EXONS[:-1], EXONS[1:]):
        a, b = x0 + left_end * STEP_WIDTH, x0 + right_start * STEP_WIDTH
        mid = (a + b) / 2
        ax.plot([a, mid, b], [y0, y0 + rise, y0], color=JUNCTION, lw=0.9)


def gtf_icon(ax, x0, y0):
    ax.add_patch(Rectangle((x0, y0), 0.32, 0.42, facecolor="white", edgecolor=COVERAGE, lw=0.7))
    for k in range(4):
        ax.plot([x0 + 0.06, x0 + 0.26], [y0 + 0.32 - 0.07 * k] * 2, color=COVERAGE, lw=0.6)


rng = np.random.default_rng(1)
fig, ax = plt.subplots(figsize=(6.3, 1.45))
starts = [k * (STEP_WIDTH + GAP) for k in range(5)]

coverage(ax, starts[0], 0.62, 0.22, 1.0, 0.05, rng)
coverage(ax, starts[0], 0.30, 0.22, 0.8, 0.05, rng)
arcs(ax, starts[0], 0.10)
regions(ax, starts[0], 0.07, 0.06, JUNCTION)

coverage(ax, starts[1], 0.30, 0.45, 0.9, 0.12, rng)

coverage(ax, starts[2], 0.36, 0.40, 0.9, 0.12, rng)
ax.plot([starts[2], starts[2] + STEP_WIDTH], [0.36 + 0.40 * 0.2] * 2,
        color=THRESHOLD, lw=0.8, ls=(0, (3, 2)))
regions(ax, starts[2], 0.10)

arcs(ax, starts[3], 0.40, rise=0.22)
regions(ax, starts[3], 0.34)

gtf_icon(ax, starts[4] + 0.3, 0.25)

labels = [
    "coverage per\nsample (bigWig),\njunctions (RR, MM)",
    "CPM, mean\nover samples",
    "regions above\nmin_coverage",
    "regions joined\nby a junction",
    "write GTF",
]
for k, (x0, text) in enumerate(zip(starts, labels)):
    ax.text(x0 + STEP_WIDTH / 2, -0.08, text, ha="center", va="top", fontsize=TEXT_SIZE)
    ax.text(x0, 1.0, str(k + 1), ha="left", va="top", fontsize=TEXT_SIZE, fontweight="bold")
    if k < 4:
        ax.add_patch(FancyArrowPatch((x0 + STEP_WIDTH + 0.05, 0.5), (starts[k + 1] - 0.05, 0.5),
                                     arrowstyle="-|>", mutation_scale=7, lw=0.8, color="#555555"))

ax.set_xlim(-0.05, starts[-1] + STEP_WIDTH + 0.05)
ax.set_ylim(-0.42, 1.02)
ax.axis("off")

for ext in ("pdf", "svg", "png"):
    fig.savefig(f"fig_pipeline_schematic.{ext}", bbox_inches="tight", dpi=200)
print("wrote fig_pipeline_schematic.pdf, .svg and .png")
