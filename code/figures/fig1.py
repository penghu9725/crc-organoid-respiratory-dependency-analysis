# -*- coding: utf-8 -*-
"""FIG 1: study design — 1a workflow, 1b cohort 256, 1c screened 162, 1d universe, 1e consensus flow."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from style import (save, read_tsv, num, MM, PANELS, ASSEMBLED, MOD_COL, ROBUST_MODULES,
                   CRC, ESCA, GREY, BASE, NEUTRALS)

META = os.path.join(os.path.abspath(os.environ.get("ORGANOID_METADATA_ROOT", os.path.join(BASE, "data", "external_metadata"))), "ORGANOID_METADATA.tsv")
DARK = "#333333"

def _flow_box(ax, x, y, w, h, text, fc, tc="white", fs=6.6, ec="none", lw=0.5, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.04",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", wrap=True)

def _arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7,
                                 color="#555555", lw=0.9))

# ---- 1a workflow
def panel_1a(ax):
    ax.set_xlim(0, 20); ax.set_ylim(0, 4); ax.axis("off")
    steps = [("Source organoid biobank\n256 models; public release", "#4A7FB5"),
             ("CRISPR dependency data\n162 models \u00b7 2 libraries", "#6A9BC7"),
             ("1,412-gene metabolic\nuniverse (Tier 1+2)", CRC),
             ("CRC-vs-ESCA discovery\n(85 vs 59 organoids)", CRC),
             ("13 respiratory\nsubmodules", "#3D6E9E"),
             ("2D comparison \u00b7 patient RNA /\nprotein context", "#8AAECF")]
    x = 0.35
    for i, (txt, c) in enumerate(steps):
        w = 2.9 if i in (0, 1, 5) else 2.7
        _flow_box(ax, x, 1.1, w, 1.6, txt, c, fs=6.2)
        if i < len(steps) - 1:
            _arrow(ax, x + w + 0.08, 1.9, x + w + 0.42, 1.9)
        x += w + 0.5
    ax.text(0.35, 3.5, "No new patient or animal data generated; all analyses on public data",
            fontsize=6.2, style="italic", color="#666666")
    ax.set_title("a", loc="left", fontsize=10, fontweight="bold")

# ---- 1b / 1c cohort bars
def _cohort_ax(ax, meta, screened_only, title):
    from collections import Counter, OrderedDict
    order = ["Colorectal", "Oesophageal", "Ovarian", "Pancreatic", "Gastric"]
    if screened_only:
        rows = [r for r in meta if r.get("CRISPR_available", "").strip().lower() in ("yes", "true", "1")]
    else:
        rows = meta
    cnt = Counter(r["primary_tumour_type"].strip() for r in rows)
    n = [cnt.get(o, 0) for o in order]
    cols = [CRC, ESCA, NEUTRALS["Ovarian"], NEUTRALS["Pancreatic"], NEUTRALS["Gastric"]]
    y = np.arange(len(order))
    ax.barh(y, n, color=cols, height=0.62)
    for yi, ni in zip(y, n):
        ax.text(ni + 3, yi, str(ni), va="center", fontsize=7)
    ax.set_yticks(y); ax.set_yticklabels([o for o in order], fontsize=6.6)
    ax.set_xlim(0, 150)
    ax.set_xlabel("organoid models", fontsize=7)
    ax.set_title(title, fontsize=7.6, fontweight="bold")
    ax.tick_params(axis="x", labelsize=6.5)

def panel_1bc(ax_b, ax_c):
    meta = read_tsv(META)
    _cohort_ax(ax_b, meta, False, "b  256-model organoid biobank")
    _cohort_ax(ax_c, meta, True, "c  CRISPR-screened subset (n = 162)")

# ---- 1d universe
def panel_1d(ax):
    tiers = [("Tier 1\ncanonical membership", 2332, CRC),
             ("Tier 2\nregulatory/response", 716, "#6A9BC7"),
             ("Tier 3\nGO ontology (excluded)", 8614, "#D9D9D9")]
    y = np.arange(3)[::-1]
    for yi, (lab, n, c) in zip(y, tiers):
        ax.barh(yi, n, color=c, height=0.6)
        ax.text(n + 150, yi, f"{n:,} rows", va="center", fontsize=6.6)
    ax.set_yticks(y); ax.set_yticklabels([t[0] for t in tiers], fontsize=6.4)
    ax.set_xscale("log"); ax.set_xlim(100, 30000)
    ax.set_xlabel("membership rows (log)", fontsize=6.8)
    ax.text(100, 3.6, "1,412 unique Tier-1/Tier-2 genes dependency-tested; 27 metabolic categories; 3,734 unique genes",
            fontsize=6.2, color="#333333")
    ax.text(100, 3.15, "OXPHOS: 299 Tier-1 members (246 evaluable) \u2014 flagged aggregation-sensitive",
            fontsize=6.2, color="#8a5a00")
    ax.set_title("d", loc="left", fontsize=10, fontweight="bold")

# ---- 1e consensus flow
def panel_1e(ax):
    ax.set_xlim(0, 16); ax.set_ylim(0, 3.2); ax.axis("off")
    _flow_box(ax, 0.3, 0.9, 3.0, 1.4, "MinLib screen", "#B0B0B0", tc="white", fs=7)
    _flow_box(ax, 0.3 + 3.3, 0.9, 3.0, 1.4, "Yusa v1.1 screen", "#B0B0B0", tc="white", fs=7)
    ax.text(3.6, 0.55, "16 models measured with both libraries", fontsize=5.8, ha="center", color="#777777")
    _flow_box(ax, 7.1, 0.9, 4.2, 1.4, "gene-wise library-adjusted\nconsensus LFC", CRC, fs=6.6)
    _flow_box(ax, 11.7, 0.9, 3.9, 1.4, "organoid-level\nconsensus (n = 162)", "#3D6E9E", fs=6.6)
    _arrow(ax, 3.3, 1.6, 7.05, 1.6); _arrow(ax, 6.6, 1.6, 7.05, 1.6)
    _arrow(ax, 11.3, 1.6, 11.65, 1.6)
    ax.text(7.6, 2.7, "cross-library agreement: Spearman r = 0.363 (moderate) \u2014 adjustment prevents\n"
                      "16 dual-library models being double-counted", fontsize=5.9, color="#555555")
    ax.set_title("e", loc="left", fontsize=10, fontweight="bold")

def assemble():
    W, H = 180 * MM, 118 * MM
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.2, 0.8],
                          left=0.06, right=0.99, top=0.96, bottom=0.09,
                          hspace=0.55, wspace=0.55)
    ax_a = fig.add_subplot(gs[0, :]); panel_1a(ax_a)
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    panel_1bc(ax_b, ax_c)
    ax_d = fig.add_subplot(gs[1, 2]); panel_1d(ax_d)
    ax_e = fig.add_subplot(gs[2, :]); panel_1e(ax_e)
    save(fig, os.path.join(ASSEMBLED, "Figure1"))

def panels_only():
    specs = [("1a", panel_1a, (180 * MM, 30 * MM)),
             ("1b", lambda ax: _cohort_ax(ax, read_tsv(META), False, "b  256-model biobank"), (58 * MM, 40 * MM)),
             ("1c", lambda ax: _cohort_ax(ax, read_tsv(META), True, "c  screened subset (n = 162)"), (58 * MM, 40 * MM)),
             ("1d", panel_1d, (58 * MM, 40 * MM)),
             ("1e", panel_1e, (180 * MM, 26 * MM))]
    for name, fn, size in specs:
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([0.03, 0.03, 0.95, 0.92])
        fn(ax)
        save(fig, os.path.join(PANELS, f"Fig1_{name}"))

if __name__ == "__main__":
    panels_only()
    assemble()
    print("FIG1 DONE")
