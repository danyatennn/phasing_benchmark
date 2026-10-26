"""Grouped bars for the comparison figures: a panel per platform, a group of
bars per phaser (or caller, or method), a bar per thing compared inside the
group (input, read truth, sample). Shared by figures_extra.py and
figures_experiments.py so every comparison figure reads the same way.

Every bar carries its value; a bar compared with the first of its group also
carries the change (x6.7, +56%), because heights on a log axis are not ratios
one can read by eye. A rate that differs by orders of magnitude between groups
goes on a log axis, its bars starting at the bottom of the axis.
"""

import style as st
from figures import finite, log_ticks


def offsets(n, width=0.8):
    """x offsets of n bars side by side filling width, and the width of one"""
    w = width / max(n, 1)
    return [(i - (n - 1) / 2) * w for i in range(n)], w


def bar(ax, x, value, width, color, alpha=1.0, hatch=None):
    """One bar; None if there is nothing to draw."""
    if not finite(value):
        return None
    ax.bar(x, value, width=width * 0.9, color=color, alpha=alpha, hatch=hatch,
           edgecolor="white" if hatch else "none", linewidth=0, zorder=3)
    return value


def value_text(v):
    """2-3 significant digits, never a power of ten: 0.042, 0.23, 1.13, 7.6, 145"""
    a = abs(v)
    if a >= 100:
        return "%.0f" % v
    if a >= 10:
        return "%.1f" % v
    if a >= 0.1:
        return "%.2f" % v
    return "%.2g" % v


def change_text(new, old):
    """x6.7 for a large change, +56% / -12% for a small one."""
    if not (finite(new) and finite(old)) or old <= 0 or new <= 0:
        return ""
    r = new / old
    if 0.5 < r < 2:
        return "%+.0f%%" % (100 * (r - 1))
    return "×%.1f" % r if r >= 2 else "×%.2f" % r


def label_bar(ax, x, value, width, text=None):
    """The value over a bar (under it when negative), upright on a wide bar and
    turned on a narrow one; text replaces the plain value."""
    if not finite(value):
        return
    text = text if text is not None else value_text(value)
    upright = width >= 0.5
    below = value < 0
    ax.annotate(text, (x, value), xytext=(0, -1.5 if below else 1.5), textcoords="offset points",
                ha="center", va="top" if below else "bottom", rotation=0 if upright else 90,
                fontsize=st.FONT_SIZE - 2, color=st.INK, zorder=4, annotation_clip=False)


def group_ticks(ax, labels):
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.6, len(labels) - 0.4)


def scale_axes(axes, drawn, log, ylabel):
    """Shared y of a panel row, with room above the bars for their labels: on
    a log axis from half the smallest bar to a decade over the largest; on a
    linear one from 0 (or under the lowest negative bar)."""
    if log:
        values = [v for v in drawn if finite(v) and v > 0]
        lo, hi = min(values), max(values)
        axes[0].set_yscale("log")
        axes[0].set_ylim(lo / 2, hi * 10)
        log_ticks(axes[0].yaxis, lo / 2, hi * 10)
    else:
        values = [v for v in drawn if finite(v)]
        lo, hi = min(values + [0]), max(values + [0])
        axes[0].set_ylim(lo * 1.4 if lo < 0 else 0, hi * 1.35 if hi > 0 else 1)
    axes[0].set_ylabel(ylabel)
    for ax in axes:
        ax.grid(axis="y", color=st.GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
