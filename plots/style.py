"""One look for every figure: the grid, its labels, its colours, the rcParams.

Colour follows the entity. A platform or a phaser wears the same colour in every
panel of every figure, and the two families never share a hue:
  platforms  ONT purple (#3c1642), PacBio teal (#086375), a lighter tint for the
             older chemistry of each
  phasers    crimson, navy, amber, sky - checked pairwise under protan, deutan
             and tritan simulation (min OKLab dE 12), and each also has its own
             marker so identity never rests on colour alone
  callers    dark green, orchid, lime, separated by lightness as much as by hue
             (min dE 8); the GIAB truth is grey, it is the reference, not a caller
"""

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap

ROOT = "/autofs/bal36/dten/phasing_benchmark"

SAMPLE = "HG002"
REGION = "wgs"

REFERENCES = ["GRCh38", "T2T"]
PLATFORMS = ["ont4k", "ont5k", "hifi", "revio"]
DEPTHS = [10, 20, 30, 40, 50]
CALLERS = ["truth", "clair3", "deepvariant", "longcalld"]
PHASERS = ["whatshap", "longphase", "margin", "hapcut2", "hiphase"]
# a phaser run on some platforms only; every other one runs on all of them
PHASER_PLATFORMS = {"hiphase": ("hifi", "revio")}
# methylation-aware runs of meth_run.sh (truth input, 30x), next to the phasers:
# LongPhase with its modcall VCF, and MethPhaser on top of each phaser's VCF
METH_PHASERS = ["lpmod"] + ["mp" + p for p in PHASERS]

REFERENCE_LABEL = {"GRCh38": "GRCh38", "T2T": "T2T-CHM13"}
PLATFORM_LABEL = {
    "ont4k": "ONT 4 kHz",
    "ont5k": "ONT 5 kHz",
    "hifi": "PacBio Sequel II",
    "revio": "PacBio Revio",
}
CALLER_LABEL = {
    "truth": "GIAB truth",
    "clair3": "Clair3",
    "deepvariant": "DeepVariant",
    "longcalld": "longcallD",
    "truthindel": "GIAB truth + indels",
    "truthall": "GIAB truth + all HET variants",
}
PHASER_LABEL = {
    "whatshap": "WhatsHap",
    "longphase": "LongPhase",
    "margin": "Margin",
    "hapcut2": "HapCUT2",
    "hiphase": "HiPhase",
    "lpmod": "LongPhase + 5mC",
    "mpwhatshap": "WhatsHap + MethPhaser",
    "mplongphase": "LongPhase + MethPhaser",
    "mpmargin": "Margin + MethPhaser",
    "mphapcut2": "HapCUT2 + MethPhaser",
    "mphiphase": "HiPhase + MethPhaser",
}

PLATFORM_COLOR = {
    "ont4k": "#9B6BA8",
    "ont5k": "#3C1642",
    "hifi": "#5FB0BE",
    "revio": "#086375",
}
PLATFORM_MARKER = {"ont4k": "o", "ont5k": "s", "hifi": "^", "revio": "D"}

PHASER_COLOR = {
    "whatshap": "#D1495B",
    "longphase": "#1F4E79",
    "margin": "#E3A034",
    "hapcut2": "#6FA8DC",
    "hiphase": "#8C5A3C",
}
PHASER_MARKER = {"whatshap": "o", "longphase": "s", "margin": "^", "hapcut2": "D", "hiphase": "v"}


def applies(platform, phaser):
    """Is this phaser run on this platform at all (HiPhase: PacBio only)?"""
    base = phaser[2:] if phaser.startswith("mp") else phaser
    return platform in PHASER_PLATFORMS.get(base, (platform,))

CALLER_COLOR = {
    "truth": "#9A9A9A",
    "clair3": "#1E5631",
    "deepvariant": "#C2559D",
    "longcalld": "#A8B536",
}
CALLER_MARKER = {"truth": "o", "clair3": "v", "deepvariant": "P", "longcalld": "X"}

INK = "#222222"
MUTED = "#7A7A7A"
GRID = "#E6E6E6"
MISSING = "#D4D4D4"

# magnitude: one neutral hue, light to dark, so it reads as neither a platform
# nor a phaser
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "slate", ["#F4F6F8", "#B9C4CF", "#6D7F92", "#2C3A4B"]
)

# Every PDF of the paper is drawn at the size it is placed in the figure page
# (figures.pptx: 170 mm wide), so text is FONT_SIZE everywhere and never
# scaled; SMALL_FONT only for the numbers over bars, CELL_FONT in heatmap cells.
FONT_SIZE = 6
TITLE_FONT = 7  # the panel title over a figure
SMALL_FONT = 5
CELL_FONT = 4.5  # the numbers of the heatmap cells, 4 mm wide in Fig. 1b
MM = 1 / 25.4
PAGE_WIDTH = 170  # mm, the width of a figure page
FULL_WIDTH = PAGE_WIDTH  # mm, a figure across the page
PANEL_HEIGHT = 50  # mm, one row of panels
SMALL_WIDTH = 56  # mm, a third of the page (Fig. 1d, 1e)

# (width, height) in mm of each panel of the paper, title and legend included.
# Fig. 1: a the schematic; b two heatmaps one under the other beside c; d, e
# three in a row. Fig. 2: a beside b; c over d beside e over f; g across.
# Fig. 3: a, b, c, d one under the other, across the page.
PAPER_PANELS = {
    "heatmap": (99, 60),           # Fig. 1b, figS2b, one per reference
    "regions": (66, 86),           # figR2, beside the heatmap of the switch errors by region
    "small": (SMALL_WIDTH, 56),    # Fig. 1e fig6b/fig6d by platform
    "fig1d": (SMALL_WIDTH, 62),    # Fig. 1d fig6e, the phaser legend under it
    "fig3": (122, 80),             # Fig. 2a, both references
    "fig4": (48, 80),              # Fig. 2b
    "lines": (110, 44),            # Fig. 2c fig2a, no legend
    "lines_legend": (110, 50),     # Fig. 2d fig2b, the legend under it
    "ratio_by_technology": (60, 47),  # Fig. 2e, 2f
    "fig5d": (PAGE_WIDTH, 88),     # Fig. 2g
    "experiment": (PAGE_WIDTH, 62),         # Fig. 3a, 3c, no legend
    "experiment_legend": (PAGE_WIDTH, 68),  # Fig. 3b, 3d
    # the methylation figure: the MethPhaser schematic on top, 170 x 70 mm,
    # then B1, B3, B4 a little lower than the other experiment panels, so the
    # page stays within 230 mm: 70 + 50 + 50 + 56 and the gaps
    # Fig. 3: A1 over A3 (the legend under A3) on the left, the methylation
    # NG50 (figB3_side, horizontal bars) on their right; 130 + 4 + 36 = 170
    "experiment_narrow": (130, 46),         # figA1, a row of four panels
    "experiment_narrow_legend": (130, 50),  # figA3, the legend under it
    "experiment_grid": (130, 70),           # two platforms to a row (needs ~70 mm a metric)
    "experiment_grid_legend": (130, 76),
    "meth_side": (38, 96),                  # figB3_side, as high as A1 + A3
    # the assemblies against the phasers beside the chromosome figure of the
    # switch errors (assembly gaps)
    "asm_side": (52, 128),                  # figC1_side: switch | NG50, phasers on truth and on hifiasm SNPs, 4 platforms
    "asm_scatter": (56, 102),               # figC6 (Fig. 6c): beside the chromosome figure (6b, ~110 x 102 mm), as high as it
    # Fig. 2 (the page with the chromosome figure): b and c, side by side
    "samples_paper": (83, 34),              # fig9a/fig9b of ONT 5 kHz, the three samples
    # Fig. 1c: WhatsHap on T2T-CHM13, genome-wide to all difficult regions
    # Fig. 1d (the NG50 ratios by platform) beside 1e, between the heatmap
    # and the GIAB samples: as high as each other, their axes pinned to the
    # same top and bottom (ROW_1DE) so the two plots line up
    "ratio_by_platform": (45, 46),
    "r7_small": (52, 46),
    "meth_scatter": (46, 108),              # figB6: the platforms one over the other, as high as the region view + B5
    "meth_longphase": (121, 46),            # figB5, under the region view of the methylation figure
    # the read page, 225 mm: R5 | R3c, fig10d | R20 (~110 as placed), U2 | T3 (56), U1 (48), gaps
    "upset_crossref": (PAGE_WIDTH, 48),     # figU1: reads of a stratum on the other reference
    "upset_labels": (92, 56),               # figU2: the three parental labelling methods
    "runtime_mirror": (76, 56),             # figT3: run time up, memory down, one chromosome
    "experiment_meth": (PAGE_WIDTH, 50),         # figB1, figB3
    "experiment_meth_legend": (PAGE_WIDTH, 56),  # figB4, the legend under it
    "regions_profile": (84, 80),   # figR13a beside figR13b
    "bubbles": (84, 58),           # figR19a beside figR19b
    "enrichment": (103, 90),       # fig11a, under the bubbles
    # the read page: R5 | R3c over fig10d | R20 over R18 | the ideogram
    "donuts": (84, 38),            # figR5 of one platform, GRCh38 beside T2T
    "snps_per_read": (84, 38),     # figR3c, beside figR5
    "confidence": (70, 52),        # fig10d, under figR5
    "wrong_vs_snps": (70, 60),     # figR3d, the place of fig10d
    "transitions": (98, 52),       # figR20 of one platform, beside fig10d
    "r18": (100, 56),              # figR18, under fig10d, beside the ideogram
    "r18_fig1c": (62, 102),        # figR18 of read recall in bases, Fig. 1c: two columns, beside 1b and 1d/1e
    "rescued": (84, 62),           # figR21
}

# Fig. 1d and 1e: mm from the top and the bottom of the figure to the axes
ROW_1DE = (12, 17.2)

FONT_STACK = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]


def register_fonts(font_dir):
    """Arial is rarely installed on a Linux server. Any .ttf dropped into
    plots/fonts/ is picked up here; otherwise the closest available face is used
    and said so. Text stays editable either way (pdf.fonttype 42)."""
    if os.path.isdir(font_dir):
        for name in sorted(os.listdir(font_dir)):
            if name.lower().endswith((".ttf", ".otf")):
                try:
                    font_manager.fontManager.addfont(os.path.join(font_dir, name))
                except Exception:
                    pass
    installed = {f.name for f in font_manager.fontManager.ttflist}
    stack = [f for f in FONT_STACK if f in installed] or ["DejaVu Sans"]
    return stack


def apply_style(font_dir):
    stack = register_fonts(font_dir)
    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": stack,
        "mathtext.default": "regular",
        "font.size": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "legend.fontsize": FONT_SIZE,
        "legend.title_fontsize": FONT_SIZE,
        "figure.titlesize": TITLE_FONT,
        "figure.titleweight": "normal",
        "font.weight": "normal",
        "axes.titleweight": "normal",
        "axes.labelweight": "normal",
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlepad": 4,
        "axes.prop_cycle": cycler(color=[PHASER_COLOR[p] for p in PHASERS]),
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "lines.markeredgewidth": 0.5,
        "patch.linewidth": 0.6,
        "hatch.linewidth": 0.4,
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.columnspacing": 1.4,
        "legend.handletextpad": 0.5,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    return stack[0]
