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
    "truthall": "GIAB truth + all het variants",
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

FONT_SIZE = 7
MM = 1 / 25.4
FULL_WIDTH = 180  # mm, a two-column figure
PANEL_HEIGHT = 50  # mm, one row of panels
SMALL_WIDTH = 64  # mm, a narrow figure that sits beside others of its width (fig6 by platform, fig7g)

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
