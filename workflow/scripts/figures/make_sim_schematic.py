"""Simulation schematic (Main Figure 1, panel B).

Reference transcript plus one variant per ASimulatoR event class, and the two
expression scenarios on the exon-skipping example. Styled to match Martina's
pipeline schematic (Matplotlib, DejaVu Sans, rose expressed-region boxes, blue
for the alternative element, grey lines). Output: fig_sim_schematic.pdf/.png.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "text.color": "#222222",
    "axes.edgecolor": "#555555",
    "pdf.fonttype": 42,  # keep labels as editable text, not outlines
    "svg.fonttype": "none",
})

# Blue and gold, clear of the tool palette (teal/orange/periwinkle/pink) and the
# strand classes (grey/red/brown), so the schematic never reads as a data colour.
EXON = "#3b6fa0"        # reference exon
ALT = "#e0a52e"         # alternative exon
SKIP_EDGE = "#999999"   # skipped exon outline
LINE = "#555555"        # intron line and text rules
HEADER = "#555555"

exons = [(0.02, 0.10), (0.18, 0.26), (0.34, 0.42), (0.50, 0.58), (0.66, 0.74), (0.82, 0.95)]
exon_h = 0.5
LABEL_X = -0.03


def intron(ax, y, x0, x1):
    ax.plot([x0, x1], [y, y], color=LINE, lw=0.9, zorder=1)


def exon(ax, y, span, color):
    ax.add_patch(Rectangle((span[0], y - exon_h / 2), span[1] - span[0], exon_h,
                           facecolor=color, edgecolor="none", zorder=2))


# Skipped exon: empty box, dashed grey outline.
def skipped(ax, y, span):
    ax.add_patch(Rectangle((span[0], y - exon_h / 2), span[1] - span[0], exon_h,
                           facecolor="white", edgecolor=SKIP_EDGE, linestyle=(0, (3, 2)),
                           lw=0.9, zorder=2))


def transcript(ax, y, present, alt=None, skip=None, retained=None):
    alt, skip = alt or {}, skip or set()
    intron(ax, y, exons[0][0], exons[-1][1])
    if retained:
        exon(ax, y, retained, ALT)
    for i, span in enumerate(exons):
        if i in alt:
            exon(ax, y, alt[i], ALT)
        elif i in skip:
            skipped(ax, y, span)
        elif i in present:
            exon(ax, y, span, EXON)


# Wide and short: event classes on the left, expression scenarios on the right,
# so the panel sits as one method row across the full figure width.
fig, (ax, ax2) = plt.subplots(
    1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [5, 3], "wspace": 0.06})

SITE_SHIFT = 0.03  # how far an alternative splice site moves an exon edge
# label, exons drawn as reference, exons drawn as alternative, exons skipped,
# retained intron
rows = [
    ("reference", set(range(6)), None, set(), None),
    ("exon skipping (es)", {0, 1, 3, 4, 5}, None, {2}, None),
    ("multiple exon skipping (mes)", {0, 1, 4, 5}, None, {2, 3}, None),
    ("intron retention (ir)", set(range(6)), None, set(), (exons[2][1], exons[3][0])),
    ("alternative 3' splice site (a3)", {0, 1, 2, 4, 5},
     {3: (exons[3][0] + SITE_SHIFT, exons[3][1])}, set(), None),
    ("alternative 5' splice site (a5)", {0, 1, 3, 4, 5},
     {2: (exons[2][0], exons[2][1] - SITE_SHIFT)}, set(), None),
    ("mutually exclusive exons (mee)", {0, 1, 4, 5}, {3: exons[3]}, {2}, None),
    ("alternative first exon (afe)", {1, 2, 3, 4, 5}, {0: exons[0]}, set(), None),
    ("alternative last exon (ale)", {0, 1, 2, 3, 4}, {5: exons[5]}, set(), None),
]

TOP = len(rows) + 1.0  # shared top, so both headers line up
y = len(rows)
for label, present, alt, skip, retained in rows:
    transcript(ax, y, present, alt, skip, retained)
    ax.text(LABEL_X, y, label, ha="right", va="center", fontsize=10, color="#222222")
    y -= 1

ax.set_xlim(-0.42, 1.0)
ax.set_ylim(-0.2, TOP)
ax.axis("off")
ax.text(LABEL_X, len(rows) + 0.7, "Alternative-splicing event classes",
        fontsize=10.5, color=HEADER, ha="left")

# Legend strip below the transcripts.
lx = 0.0
ax.add_patch(Rectangle((lx, 0.12), 0.045, 0.16, facecolor=EXON))
ax.text(lx + 0.055, 0.2, "reference exon", fontsize=8.5, va="center", color="#222222")
ax.add_patch(Rectangle((lx + 0.22, 0.12), 0.045, 0.16, facecolor=ALT))
ax.text(lx + 0.275, 0.2, "alternative exon, splice site or retained intron", fontsize=8.5, va="center", color="#222222")
ax.add_patch(Rectangle((lx + 0.80, 0.12), 0.045, 0.16, facecolor="white",
                       edgecolor=SKIP_EDGE, linestyle=(0, (3, 2)), lw=0.9))
ax.text(lx + 0.855, 0.2, "skipped", fontsize=8.5, va="center", color="#222222")


def coverage(ax, y, heights):
    intron(ax, y, exons[0][0], exons[-1][1])
    for span, h in zip(exons, heights):
        if h > 0:
            ax.add_patch(Rectangle((span[0], y), span[1] - span[0], h * 0.55,
                                   facecolor=EXON, edgecolor="none"))


# Scenario labels sit above each track, since the right panel is narrow.
ax2.text(0.0, len(rows) + 0.7, "Expression scenarios (es example)",
         fontsize=10.5, color=HEADER, ha="left")
ax2.text(0.0, TOP - 2.45, "reference and variant", ha="left", va="bottom", fontsize=9.5, color="#222222")
coverage(ax2, TOP - 3.15, [1, 1, 0.5, 1, 1, 1])
ax2.text(0.0, TOP - 4.25, "variant", ha="left", va="bottom", fontsize=9.5, color="#222222")
coverage(ax2, TOP - 4.95, [1, 1, 0.0, 1, 1, 1])
ax2.set_xlim(-0.05, 1.0)
ax2.set_ylim(-0.2, TOP)
ax2.axis("off")

fig.savefig("fig_sim_schematic.pdf", bbox_inches="tight")
fig.savefig("fig_sim_schematic.svg", bbox_inches="tight")
fig.savefig("fig_sim_schematic.png", dpi=200, bbox_inches="tight")
print("wrote fig_sim_schematic.pdf and .png")
