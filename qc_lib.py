"""
qc_lib.py — shared publication-quality figure theme + helpers for HPV Nanopore QC.

Every figure produced via save_fig() is written as VECTOR (SVG + PDF) AND 300-dpi PNG,
with embedded fonts, a colourblind-safe palette, a consistent theme, and journal-column
sizing. Import this in qc_report.py and any custom plotting so all figures match.

Runtime env: hpv-py (matplotlib, seaborn; plotly+kaleido optional for interactive/extra).
"""
from __future__ import annotations
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl

# ---- amplicon-aware constants carried from the v5 design ----
AMPLICON_MIN_BP = 2000          # tiling amplicon design lower bound
AMPLICON_MAX_BP = 4000          # tiling amplicon design upper bound
HPV52_GAPFILL_BP = 1117         # HPV52 amp7 L1 gap-fill — MUST NOT be filtered out
FILTER_LOWER_BP = 900           # recommended length-filter lower bound (keeps HPV52 amp7)
FILTER_UPPER_BP = 4800          # recommended length-filter upper bound

HPV_TYPES = ["HPV16","HPV18","HPV31","HPV33","HPV35","HPV39",
             "HPV45","HPV52","HPV56","HPV58","HPV59","HPV66"]
CLADE = {"HPV16":"Alpha9","HPV31":"Alpha9","HPV33":"Alpha9","HPV35":"Alpha9","HPV52":"Alpha9",
         "HPV58":"Alpha9","HPV18":"Alpha7","HPV39":"Alpha7","HPV45":"Alpha7","HPV59":"Alpha7",
         "HPV56":"Alpha6","HPV66":"Alpha6"}

# Okabe–Ito colourblind-safe palette
OKABE_ITO = ["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00",
             "#56B4E9","#F0E442","#000000","#999999"]
POOL_COLORS = {"1": "#0072B2", "2": "#D55E00", 1: "#0072B2", 2: "#D55E00"}
CLADE_COLORS = {"Alpha9": "#0072B2", "Alpha7": "#D55E00", "Alpha6": "#009E73"}
# 12 distinct type colours (clade-family hues)
TYPE_COLORS = {
    "HPV16":"#08519c","HPV31":"#3182bd","HPV33":"#6baed6","HPV35":"#9ecae1","HPV52":"#2171b5","HPV58":"#c6dbef",
    "HPV18":"#a50f15","HPV39":"#de2d26","HPV45":"#fb6a4a","HPV59":"#fcae91",
    "HPV56":"#238b45","HPV66":"#74c476",
}
# colourblind-safe sequential map for log-depth heatmaps
DEPTH_CMAP = "viridis"

# journal sizing (inches)
COL_W = 3.46        # single column ~88 mm
FULL_W = 7.20       # full width ~183 mm


def set_theme():
    """Apply the shared publication theme. Fonts embedded for vector output."""
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,          # embed TrueType (editable, not outlined) in PDF
        "ps.fonttype": 42,
        "svg.fonttype": "none",      # keep text as text in SVG
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.titleweight": "bold",
        "axes.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "legend.frameon": False,
        "lines.linewidth": 1.2,
        "figure.constrained_layout.use": True,
    })


def save_fig(fig, outdir, basename, formats=("svg", "pdf", "png")):
    """Write fig to outdir/<basename>.{svg,pdf,png} (vector + 300-dpi PNG). Returns paths."""
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in formats:
        p = outdir / f"{basename}.{ext}"
        fig.savefig(p, dpi=300 if ext == "png" else None)
        paths.append(str(p))
    plt.close(fig)
    return paths


def annotate_length_window(ax, ymax=None):
    """Shade the expected amplicon length window and mark the HPV52 1.1 kb gap-fill."""
    ax.axvspan(AMPLICON_MIN_BP, AMPLICON_MAX_BP, color="#009E73", alpha=0.10,
               label="2–4 kb tiling amplicons", zorder=0)
    ax.axvline(HPV52_GAPFILL_BP, color="#D55E00", ls="--", lw=1.0,
               label=f"HPV52 gap-fill ({HPV52_GAPFILL_BP} bp)")
    ax.axvline(FILTER_LOWER_BP, color="#999999", ls=":", lw=0.9,
               label=f"filter lower bound ({FILTER_LOWER_BP} bp)")


def footer(fig, text):
    """Small provenance/footer note bottom-left (figure coords)."""
    fig.text(0.005, 0.002, text, fontsize=5.5, color="#666666", ha="left", va="bottom")
