#!/usr/bin/env python3
"""Render the HPV Nanopore pipeline runtime-analysis workflow schematic.

Diagram A — end-to-end data flow from instrument output to results/, faithful to
PIPELINE_STAGES.md, pipeline/run_pipeline.py + pipeline/hpvpipe.py, and the manual's
stage chapter. Pure matplotlib (already in the hpv-py stack); emits vector SVG + PDF
and a 300-dpi PNG. No network dependency.

Usage:  python docs/figures/make_workflow_figure.py
Output: docs/figures/pipeline_workflow.{svg,pdf,png}
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

# ---------------------------------------------------------------- palette
# Okabe-Ito colour-blind-safe hues; pale fill + saturated edge so each conda
# env is also separable in greyscale (distinct luminance) and labelled in-box.
ENV = {
    "hpv-py":    dict(fill="#D6E6F2", edge="#0072B2"),  # blue   - orchestration, typing, summaries
    "hpv-qc":    dict(fill="#CDEBDD", edge="#009E73"),  # green  - read QC, filter, coverage tools
    "hpv-artic": dict(fill="#F7DCC9", edge="#D55E00"),  # verm4  - ARTIC + Clair3 consensus
    "hpv-asm":   dict(fill="#EFD9E6", edge="#CC79A7"),  # purple - optional de novo path
}
INK = "#1A1A1A"
BAND = dict(fill="#ECECEC", edge="#5A5A5A")

plt.rcParams.update({
    "font.family": "DejaVu Sans",   # full glyph coverage -> 0 missing glyphs
    "svg.fonttype": "none",          # keep text as text in the SVG
    "pdf.fonttype": 42,
})

fig = plt.figure(figsize=(7.4, 10.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


def box(x, y, w, h, *, env=None, band=False, dashed=False, lw=1.6, z=2):
    style = BAND if band else ENV[env]
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                       linewidth=lw, edgecolor=style["edge"], facecolor=style["fill"],
                       linestyle=(0, (5, 3)) if dashed else "solid", zorder=z)
    ax.add_patch(p)
    return (x + w / 2, y + h / 2)


def vtext(cx, top, lines, *, gap=None):
    """Stack (text, size, weight, colour) tuples downward from `top`."""
    y = top
    for (s, sz, wt, col) in lines:
        ax.text(cx, y, s, ha="center", va="top", fontsize=sz, fontweight=wt,
                color=col, zorder=5)
        y -= (gap if gap else sz / 760.0)
    return y


def chip(cx, cy, label, env, size=7.5):
    st = ENV[env]
    ax.text(cx, cy, label, ha="center", va="center", fontsize=size, fontweight="bold",
            color=st["edge"], zorder=6,
            bbox=dict(boxstyle="round,pad=0.26", fc="white", ec=st["edge"], lw=1.1))


def down_arrow(cx, y0, y1, label=None):
    ax.add_patch(FancyArrowPatch((cx, y0), (cx, y1), arrowstyle="-|>",
                 mutation_scale=16, lw=1.8, color="#444444", zorder=1,
                 shrinkA=0, shrinkB=0))
    if label:
        ax.text(cx + 0.012, (y0 + y1) / 2, label, ha="left", va="center",
                fontsize=7.2, style="italic", color="#555555", zorder=5)


# ---------------------------------------------------------------- geometry
COL_X, COL_W = 0.045, 0.575          # main vertical flow column
COL_CX = COL_X + COL_W / 2
SB_X, SB_W = 0.655, 0.325            # right sidebar
H = 0.085                            # stage box height
GAP = 0.023                          # arrow gap between stages

# ---- instrument header band
itop, ih = 0.963, 0.083
box(COL_X, itop - ih, COL_W, ih, band=True, lw=1.7)
vtext(COL_CX, itop - 0.012, [
    ("INSTRUMENT  ·  GridION / PromethION P2 Solo", 10.5, "bold", INK),
    ("R10.4.1, 400 bps  ·  basecalling + demultiplexing (on-instrument)", 8.6, "normal", "#333333"),
    ("output:  fastq_pass/barcodeNN/*.fastq.gz  (+ sequencing_summary)", 8.3, "normal", "#333333"),
], gap=0.0205)

# ---- the seven runtime stages
stages = [
    ("1", "Ingest", "hpv-py",
     "Concatenate per-barcode FASTQ into one FASTQ per sample",
     "zcat / cat  ·  barcode -> sample_id via sample_sheet.csv"),
    ("2", "Read QC", "hpv-qc",
     "Per-read length & quality stats + QC figures",
     "NanoPlot / NanoComp  ·  qc_report.py   (pycoQC: optional)"),
    ("3", "Amplicon-aware length filter", "hpv-qc",
     "Length / quality window; keeps the ~1.1 kb HPV52 amp7",
     "chopper --minlength 900 --maxlength 4800 -q 9"),
    ("4", "Type assignment / binning", "hpv-py",
     "Competitive map vs 12 HPV refs; best-hit/read; split per type",
     "minimap2 -ax map-ont  ·  hpvpipe.py bin  (co-infection aware)"),
    ("5", "Per-type consensus", "hpv-artic",
     "ARTIC consensus/type -> consensus FASTA, VCF, trimmed BAM",
     "artic minion  ·  Clair3 2.x  ·  hpvpipe.py consensus"),
    ("6", "Coverage QC", "hpv-py",
     "Per-amplicon depth and dual-pool balance",
     "pysam / mosdepth-compatible  ·  qc_report.py"),
    ("7", "Summary + publication figures", "hpv-py",
     "Run summary, consensus QC and publication figure set",
     "run_summary · consensus_qc · junction_proof  ·  qc_report.py"),
]

cy = {}
# place stages, recording centre y of each
top = itop - ih - GAP
for (num, name, env, what, tool) in stages:
    ytop = top
    ybot = ytop - H
    box(COL_X, ybot, COL_W, H, env=env)
    # numbered disc
    ax.add_patch(plt.Circle((COL_X + 0.038, ytop - 0.021), 0.0175,
                 facecolor=ENV[env]["edge"], edgecolor="white", lw=1.3, zorder=6))
    ax.text(COL_X + 0.038, ytop - 0.021, num, ha="center", va="center",
            fontsize=9.6, fontweight="bold", color="white", zorder=7)
    # env chip in the top-right corner (kept clear of the body text below)
    chip(COL_X + COL_W - 0.050, ytop - 0.021, env, env)
    ax.text(COL_X + 0.066, ytop - 0.0125, name, ha="left", va="top",
            fontsize=9.8, fontweight="bold", color=INK, zorder=6)
    ax.text(COL_X + 0.022, ytop - 0.040, what, ha="left", va="top",
            fontsize=8.0, color="#222222", zorder=6)
    ax.text(COL_X + 0.022, ytop - 0.058, tool, ha="left", va="top",
            fontsize=7.4, color="#444444", fontstyle="italic", zorder=6)
    cy[num] = (ytop, ybot, (ytop + ybot) / 2)
    top = ybot - GAP

# arrows between stages
prev_bot = itop - ih
for num, _, _, _, _ in stages:
    ytop, ybot, _ = cy[num]
    down_arrow(COL_CX, prev_bot, ytop)
    prev_bot = ybot

# ---- results footer band
rtop = prev_bot - GAP
rh = 0.072
down_arrow(COL_CX, prev_bot, rtop)
box(COL_X, rtop - rh, COL_W, rh, band=True, lw=1.7)
vtext(COL_CX, rtop - 0.013, [
    ("results/  in the user's working directory", 10.2, "bold", INK),
    ("consensus FASTA · VCF · per-amplicon coverage", 8.3, "normal", "#333333"),
    ("QC + publication figures · run_summary.tsv", 8.3, "normal", "#333333"),
], gap=0.0195)

# ---------------------------------------------------------------- sidebar
# (a) conductor-model note, aligned with stages 1-3
note_top = cy["1"][0]
note_bot = cy["3"][1]
box(SB_X, note_bot, SB_W, note_top - note_bot, band=True, lw=1.4)
nx = SB_X + 0.018
ax.text(SB_X + SB_W / 2, note_top - 0.015, "Orchestration — the conductor model",
        ha="center", va="top", fontsize=9.0, fontweight="bold", color=INK, zorder=6)
ax.text(nx, note_top - 0.040,
        "The whole run is driven from a single\n"
        "orchestrator process in hpv-py\n"
        "(run_pipeline.py / hpvpipe.py).\n\n"
        "Each stage's tools run in their own\n"
        "conda env, invoked on demand via\n"
        "conda run -- e.g.\n"
        "HPV_ARTIC_RUN='conda run -n hpv-artic'.\n\n"
        "The envs are NOT all activated at\n"
        "once; the conductor calls into one\n"
        "env per stage.",
        ha="left", va="top", fontsize=7.9, color="#222222", zorder=6, linespacing=1.25)

# (b) origin-junction handling — INTEGRAL to step 5: it loops out for the rotated
#     amp3 run and the merged consensus returns to the main flow (solid arrows).
J_C = ENV["hpv-artic"]["edge"]
j_top = cy["5"][0] + 0.004
j_h = 0.150
j_bot = j_top - j_h
box(SB_X, j_bot, SB_W, j_h, env="hpv-artic", lw=1.6)
ax.text(SB_X + SB_W / 2, j_top - 0.014, "Origin-junction handling (amp3)",
        ha="center", va="top", fontsize=8.9, fontweight="bold", color=J_C, zorder=6)
ax.text(SB_X + 0.018, j_top - 0.036,
        "Part of step 5 — the origin-spanning\n"
        "amplicon, recovered with two ARTIC runs:\n"
        "• MAIN: 5 non-wrapping amplicons\n"
        "  on the canonical reference\n"
        "• JUNCTION: amp3 on a ROTATED ref\n"
        "  (no negative-length crash)\n"
        "rotate back + merge the seam INTO\n"
        "step 5's consensus => seam N = 0,\n"
        "continuous -> on to coverage QC",
        ha="left", va="top", fontsize=7.5, color="#222222", zorder=6, linespacing=1.25)
# IN: step 5 -> junction (solid)
ax.add_patch(FancyArrowPatch((COL_X + COL_W, cy["5"][0] - 0.024), (SB_X, j_top - 0.028),
             arrowstyle="-|>", mutation_scale=12, lw=1.7, color=J_C, zorder=3,
             connectionstyle="arc3,rad=-0.15"))
ax.text((COL_X + COL_W + SB_X) / 2, cy["5"][0] - 0.006, "amp3\njunction run",
        ha="center", va="bottom", fontsize=6.2, color=J_C, zorder=6, linespacing=1.0)
# RETURN: merged consensus -> back into the 5->6 flow (SOLID; integral, not a dead end)
merge_y = (cy["5"][1] + cy["6"][0]) / 2
ax.add_patch(FancyArrowPatch((SB_X, j_bot + 0.040), (COL_X + COL_W, merge_y),
             arrowstyle="-|>", mutation_scale=14, lw=1.9, color=J_C, zorder=3,
             connectionstyle="arc3,rad=0.28"))
ax.text((COL_X + COL_W + SB_X) / 2 + 0.004, merge_y - 0.010, "merged\nconsensus",
        ha="center", va="top", fontsize=6.2, fontstyle="italic", color=J_C,
        zorder=6, linespacing=1.0)

# (c) optional de novo branch off STEP 5 (DASHED, optional, NO return — distinct from junction)
D_C = ENV["hpv-asm"]["edge"]
dn_top = j_bot - 0.034
dn_h = 0.082
dn_bot = dn_top - dn_h
box(SB_X, dn_bot, SB_W, dn_h, env="hpv-asm", dashed=True, lw=1.5)
ax.text(SB_X + SB_W / 2, dn_top - 0.013, "Optional — de novo path (off step 5)",
        ha="center", va="top", fontsize=8.4, fontweight="bold", color=D_C, zorder=6)
ax.text(SB_X + 0.018, dn_top - 0.035,
        "Alternative for divergent / novel /\n"
        "co-infection cases: flye + medaka\n"
        "polish. Optional — not part of every\n"
        "run, and does not return to the core flow.",
        ha="left", va="top", fontsize=7.5, color="#222222", zorder=6, linespacing=1.22)
# dashed one-way arrow step 5 -> de novo (no return)
ax.add_patch(FancyArrowPatch((COL_X + COL_W, cy["5"][1] + 0.006), (SB_X, dn_top - 0.020),
             arrowstyle="-|>", mutation_scale=12, lw=1.5, color=D_C,
             linestyle=(0, (4, 2)), zorder=2, connectionstyle="arc3,rad=-0.30"))

# ---------------------------------------------------------------- legend
ly = dn_bot - 0.030
lh = 0.118
box(SB_X, ly - lh, SB_W, lh, band=True, lw=1.2)
ax.text(SB_X + SB_W / 2, ly - 0.013, "Legend — conda environment",
        ha="center", va="top", fontsize=9.0, fontweight="bold", color=INK, zorder=6)
order = [("hpv-py", "orchestration · typing · summaries"),
         ("hpv-qc", "read QC · length filter · coverage"),
         ("hpv-artic", "ARTIC minion + Clair3 2.x"),
         ("hpv-asm", "optional de novo (flye + medaka)")]
yy = ly - 0.036
for env, desc in order:
    st = ENV[env]
    ax.add_patch(FancyBboxPatch((SB_X + 0.020, yy - 0.010), 0.030, 0.016,
                 boxstyle="round,pad=0.002,rounding_size=0.004",
                 facecolor=st["fill"], edgecolor=st["edge"], lw=1.4, zorder=6))
    ax.text(SB_X + 0.058, yy - 0.002, env, ha="left", va="center",
            fontsize=7.8, fontweight="bold", color=st["edge"], zorder=6)
    ax.text(SB_X + 0.135, yy - 0.002, desc, ha="left", va="center",
            fontsize=6.9, color="#333333", zorder=6)
    yy -= 0.0205
ax.text(SB_X + 0.020, yy - 0.001,
        "dashed outline = optional stage / branch",
        ha="left", va="center", fontsize=6.9, fontstyle="italic", color="#555555", zorder=6)

# ---------------------------------------------------------------- title
ax.text(0.045, 0.997, "HPV Nanopore pipeline — runtime analysis workflow",
        ha="left", va="top", fontsize=12.4, fontweight="bold", color=INK)

out = Path(__file__).resolve().parent
for ext in ("svg", "pdf", "png"):
    fig.savefig(out / f"pipeline_workflow.{ext}", dpi=300,
                bbox_inches="tight", pad_inches=0.06,
                facecolor="white")
    print("wrote", out / f"pipeline_workflow.{ext}")
