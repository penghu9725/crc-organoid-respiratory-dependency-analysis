# -*- coding: utf-8 -*-
"""FIG 7: 7a paired CPTAC protein, 7b RNA-protein concordance, 7c integrated evidence matrix."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from style import (save, grayscale, read_tsv, num, MM, PANELS, ASSEMBLED, PREVIEW, MOD_COL, MOD_LABEL,
                   ROBUST_MODULES, CRC, ESCA, GREY, MOD_MTFAS, BASE, TUM, NORM)

T = os.path.abspath(os.environ.get("ORGANOID_DERIVED_ROOT", os.path.join(BASE, "data", "derived", "analysis_outputs")))

# frozen exact directions for the evidence matrix (FIGURE_SOURCE_MAP authority)
MATRIX = {
    # module: (organoid rB/FDR, TCGA, SC, CPTAC)  text rows
    "CI_assembly":   dict(org="\u2191  FDR 3.9\u00d710\u207b\u2074", tcga="\u2191  FDR 2.9\u00d710\u207b\u2076",
                          sc="\u2191*  (+0.482)", cptac="\u2191  FDR 3.2\u00d710\u207b\u2075", band="CROSS-PLATFORM CONVERGENCE", bandc=MOD_COL["CI_assembly"]),
    "CI_structural": dict(org="\u2191  FDR 1.1\u00d710\u207b\u00b3", tcga="\u2193  FDR 8.2\u00d710\u207b\u00b9\u00b9",
                          sc="\u2014  n.s.", cptac="\u2193  FDR 8.5\u00d710\u207b\u00b9\u00b3", band="DEPENDENCY\u2013ABUNDANCE DECOUPLING", bandc="#444444"),
    "CIV_structural": dict(org="\u2191  FDR 5.7\u00d710\u207b\u00b3", tcga="\u2193  FDR 4.3\u00d710\u207b\u2079",
                           sc="\u2193  (FDR 0.148)", cptac="\u2193  FDR 0.033*", band="DEPENDENCY\u2013ABUNDANCE DECOUPLING", bandc="#444444"),
    "CIV_assembly":  dict(org="\u2191  FDR 1.1\u00d710\u207b\u00b3", tcga="\u2193  FDR 8.2\u00d710\u207b\u00b9\u00b9",
                          sc="\u2014  n.s.", cptac="\u2014  FDR 0.090", band="", bandc=""),
    "FeS_biogenesis": dict(org="\u2191  FDR 5.7\u00d710\u207b\u00b3", tcga="\u2014  FDR 0.33",
                           sc="\u2191  (FDR 0.156)", cptac="\u2014  FDR 0.124", band="SUPPORTIVE (RNA-stage only)", bandc=MOD_COL["FeS_biogenesis"]),
}

# ============ 7a paired protein ============
def panel_7a(ax):
    d = read_tsv(os.path.join(T, "P2_4B_CPTAC_MODULE_SUMMARY.tsv"))
    mods = [m for m in ROBUST_MODULES if m != "mtFAS"]
    prot = [r for r in d if r["modality"] == "Protein"]
    y = np.arange(len(mods))[::-1]
    for yi, m in zip(y, mods):
        r = next((x for x in prot if x["module"] == m), None)
        if r is None:
            continue
        eff = num(r["paired_median_difference_tumor_minus_normal"])
        lo = num(r["CI95_low"]); hi = num(r["CI95_high"])
        fdr = num(r["BH_FDR"])
        ax.plot([lo, hi], [yi, yi], color="#888888", lw=1.8, zorder=2)
        ax.scatter([eff], [yi], s=34, facecolor=MOD_COL[m], edgecolor="none", zorder=3)
        ax.text(0.002, yi + 0.06, MOD_LABEL[m].replace("Complex ", "C"), fontsize=6.6, va="bottom")
        tag = "FDR {:.1e}".format(fdr) if fdr < 0.05 else "FDR {:.2f}".format(fdr)
        if m == "CIV_structural":
            tag = "FDR 0.033  (CI includes 0)"
        ax.text(0.052, yi, tag, fontsize=6.2, ha="left", va="center",
                color="#111111" if fdr < 0.05 or m == "CIV_structural" else "#777777")
    ax.axvline(0, color="#333333", lw=0.7)
    ax.set_yticks([])
    ax.set_xlim(-0.055, 0.075)
    ax.set_xticks([-0.04, -0.02, 0, 0.02, 0.04, 0.06])
    ax.set_xlabel("paired tumour\u2212normal protein score\n(positive = higher in tumours; n = 96 paired)")
    ax.set_title("a", loc="left", fontsize=10, fontweight="bold")

# ============ 7b concordance ============
def panel_7b(ax):
    d = read_tsv(os.path.join(T, "P2_4B_RNA_PROTEIN_CONCORDANCE.tsv"))
    mods = ROBUST_MODULES
    mod_rows = [r for r in d if r["level"] == "module"]
    y = np.arange(len(mods))[::-1]
    for yi, m in zip(y, mods):
        r = next((x for x in mod_rows if x["module"] == m), None)
        if r is None:
            continue
        rho = num(r["spearman_rho"]); fdr = num(r["BH_FDR"])
        sig = fdr < 0.05
        ax.scatter([rho], [yi], s=30, facecolor=MOD_COL[m] if sig else "white",
                   edgecolor=MOD_COL[m], linewidth=1.0, zorder=3)
        ax.text(-0.30, yi, MOD_LABEL[m].replace("Complex ", "C").replace("Fe-S ", "Fe-S "),
                fontsize=6.4, va="center", ha="right")
        ax.text(rho + 0.012, yi, "\u03c1 = {:.3f}".format(rho), fontsize=5.8, va="center",
                color="#666666")
    ax.axvline(0, color="#333333", lw=0.6)
    ax.set_yticks([])
    ax.set_xlim(-0.32, 0.35)
    ax.set_xlabel("module RNA\u2013protein Spearman \u03c1 (96 matched)")
    n_gene = sum(1 for r in d if r["level"] == "gene" and num(r["BH_FDR"]) < 0.05)
    ax.text(0.02, -0.9, f"gene level: {n_gene}/117 evaluable genes significant (all positive)",
            fontsize=6.4, color="#333333")
    ax.set_title("b", loc="left", fontsize=10, fontweight="bold")

# ============ 7c evidence matrix ============
def panel_7c(ax):
    ax.set_xlim(0, 12.6); ax.set_ylim(0, 10.2); ax.axis("off")
    mods = list(MATRIX.keys())
    cols = ["Organoid CRISPR", "TCGA RNA", "Single-cell", "CPTAC protein"]
    # column headers
    for j, c in enumerate(cols):
        ax.text(2.35 + j * 2.3, 9.7, c, ha="center", fontsize=6.8, fontweight="bold")
    ax.text(0.35, 9.7, "module", ha="center", fontsize=6.8, fontweight="bold")
    ax.plot([0, 12.4], [9.45, 9.45], color="#333333", lw=0.8)
    for i, m in enumerate(mods):
        yrow = 8.6 - i * 1.55
        # module label with colour chip
        col = MOD_COL[m]
        ax.add_patch(Rectangle((0.15, yrow - 0.42), 0.55, 0.84, fc=col, ec="none"))
        lab = MOD_LABEL[m].replace("Fe-S biogenesis", "Fe-S")
        ax.text(0.95, yrow, lab, fontsize=7.2, va="center", fontweight="bold")
        vals = MATRIX[m]
        for j, key in enumerate(["org", "tcga", "sc", "cptac"]):
            txt = vals[key]
            xc = 2.35 + j * 2.3
            star = "*" in txt
            ax.text(xc, yrow, txt, ha="center", va="center", fontsize=6.6,
                    color="#111111")
        band = vals["band"]
        if band:
            ax.text(12.05, yrow, band, ha="right", va="center", fontsize=6.4,
                    style="italic", color=vals["bandc"] if vals["bandc"] else "#555555")
    # exploratory note for sc column
    ax.text(2.35 + 2 * 2.3, 0.55, "* exploratory; n = 6 paired patients; not FDR-confirmed",
            ha="center", fontsize=5.9, style="italic", color="#666666")
    ax.text(0.35, 0.25, "CIV_structural protein: bootstrap CI includes zero (magnitude cautious).",
            fontsize=5.9, color="#888888")
    ax.text(12.1, 9.95, "", fontsize=1)
    ax.set_title("c", loc="left", fontsize=10, fontweight="bold")

def assemble():
    W, H = 180 * MM, 118 * MM
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], width_ratios=[1.0, 1.0],
                          left=0.06, right=0.995, top=0.97, bottom=0.10, hspace=0.7, wspace=0.4)
    ax_a = fig.add_subplot(gs[0, 0]); panel_7a(ax_a)
    ax_b = fig.add_subplot(gs[0, 1]); panel_7b(ax_b)
    ax_c = fig.add_subplot(gs[1, :]); panel_7c(ax_c)
    # footer
    fig.text(0.06, 0.015, "Patient RNA/protein data represent abundance programs, not dependency measurements.",
             fontsize=7, style="italic", color="#444444")
    save(fig, os.path.join(ASSEMBLED, "Figure7"))

def panels_only():
    specs = [("7a", panel_7a, (80 * MM, 62 * MM)),
             ("7b", panel_7b, (80 * MM, 62 * MM)),
             ("7c", panel_7c, (178 * MM, 66 * MM))]
    for name, fn, size in specs:
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([0.03, 0.03, 0.95, 0.94])
        fn(ax)
        save(fig, os.path.join(PANELS, f"Fig7_{name}"))

if __name__ == "__main__":
    panels_only()
    assemble()
    grayscale(os.path.join(ASSEMBLED, "Figure7.png"), os.path.join(PREVIEW, "Figure7_gray.png"))
    print("FIG7 DONE")
