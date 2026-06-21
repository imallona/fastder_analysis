"""Sample-design schematics heading each worked-example block of Figure 2.

fig_tdp43_scheme: TDP-43 knockdown vs control.
fig_gtex_scheme:  GTEx four tissues, sub-group design.
Minimal and to the point: just the sample design, no fastder step (every panel
goes through fastder). Matches the other schematics (DejaVu Sans, blue/gold,
editable text).
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "text.color": "#222222",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

# WT blue, KD red, matching the coverage track and the TDP panels.
BLUE = "#4575b4"
RED = "#d73027"
GOLD = "#e0a52e"
GREY = "#555555"
LIGHT = "#e8edf2"


def dots(ax, x0, y, n, color):
    for i in range(n):
        ax.add_patch(plt.Circle((x0 + i * 0.05, y), 0.02, color=color, zorder=3))


# --- TDP-43: WT and KD groups, boxed like the GTEx scheme ---
fig, ax = plt.subplots(figsize=(3.8, 2.2))
ax.text(0.0, 0.95, "TDP-43 knockdown", fontsize=12, ha="left")
ax.text(0.0, 0.82, "recount3 SRP166282, motor neuron", fontsize=8.5, color=GREY, ha="left")
groups = [("TDP-43 WT", BLUE), ("TDP-43 KD", RED)]
for i, (g, c) in enumerate(groups):
    y = 0.60 - i * 0.20
    ax.add_patch(Rectangle((0.0, y), 0.40, 0.15, facecolor=LIGHT, edgecolor=c, lw=1.2))
    ax.text(0.20, y + 0.075, g, ha="center", va="center", fontsize=9)
ax.text(0.0, 0.05, "called separately", fontsize=7.5, color=GREY, ha="left", style="italic")
# Arrows from each group box to the evaluated cryptic-exon genes.
for y in (0.675, 0.475):
    ax.annotate("", xy=(0.56, 0.58), xytext=(0.40, y),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.1))
ax.text(0.58, 0.62, "cryptic-exon genes\nSTMN2, HDGFL2,\nELAVL3, CELF5, KCNQ2",
        fontsize=8, ha="left", va="center")
ax.text(0.58, 0.30, "1.0 CPM: STMN2 only\n0.02 CPM: full panel",
        fontsize=7.5, color=GREY, ha="left", va="center")
ax.set_xlim(0, 1.0)
ax.set_ylim(0.0, 1.02)
ax.axis("off")
fig.savefig("fig_tdp43_scheme.pdf", bbox_inches="tight")
fig.savefig("fig_tdp43_scheme.svg", bbox_inches="tight")
fig.savefig("fig_tdp43_scheme.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# --- GTEx: four tissues stacked vertically, sub-groups ---
# Same figsize and fonts as the TDP scheme so the two headers match in print.
fig, ax = plt.subplots(figsize=(3.8, 2.2))
ax.text(0.0, 0.95, "GTEx", fontsize=12, ha="left")
ax.text(0.0, 0.82, "recount3, genome-wide", fontsize=8.5, color=GREY, ha="left")
tissues = ["brain", "heart", "skeletal muscle", "whole blood"]
for i, t in enumerate(tissues):
    y = 0.66 - i * 0.145
    ax.add_patch(Rectangle((0.0, y), 0.46, 0.11, facecolor=LIGHT, edgecolor=GREY, lw=0.8))
    ax.text(0.23, y + 0.055, t, ha="center", va="center", fontsize=8.5)
# Arrow from the tissue stack to the fastder threshold, mirroring the TDP scheme.
ax.annotate("", xy=(0.62, 0.45), xytext=(0.46, 0.45),
            arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.1))
ax.text(0.64, 0.45, "fastder\n1.0 CPM", fontsize=8.5, ha="left", va="center")
ax.text(0.0, 0.04, "8 sub-groups of 5 samples -> 32 catalogs", fontsize=8, color=GREY, ha="left")
ax.set_xlim(0, 1.0)
ax.set_ylim(0.0, 1.02)
ax.axis("off")
fig.savefig("fig_gtex_scheme.pdf", bbox_inches="tight")
fig.savefig("fig_gtex_scheme.svg", bbox_inches="tight")
fig.savefig("fig_gtex_scheme.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("wrote fig_tdp43_scheme and fig_gtex_scheme (pdf, png)")
