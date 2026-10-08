# -*- coding: utf-8 -*-
"""FIG 4: 4a organoid vs 2D, 4b 10-gene shortlist, 4c common-essential filtering, 4d priority."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import (save, read_tsv, num, MM, PANELS, ASSEMBLED, MOD_COL, CRC, ESCA, GREY, BASE)

T = os.path.abspath(os.environ.get("ORGANOID_DERIVED_ROOT", os.path.join(BASE, "data", "derived", "analysis_outputs")))
CAT_COL = {"stronger_in_organoids": CRC, "similar_or_not_resolved": "#BBBBBB",
           "weaker_in_organoids": "#E2E2E2", "not_testable": "#AAAAAA"}
FLAG_COL = {"not_flagged": "#A6C9E2", "DepMap_inferred_common_essential": "#F0C36D",
            "Sanger_pan_cancer_core_fitness": "#E8A87C", "both_core_fitness_lists": "#C0504D"}

def panel_4a(ax):
    d = read_tsv(os.path.join(T, "P2_3A_DEPMAP_CANDIDATE_COMPARISON.tsv"))
    for r in d:
        x = num(r["organoid_effect_lineage_median"]); y = num(r["cellline_effect_matched_median"])
        cat = r["comparison_category"]
        if cat == "not_testable" or np.isnan(x) or np.isnan(y):
            continue
        ax.scatter([x], [y], s=11, c=CAT_COL[cat], alpha=0.9, linewidths=0)
    lim = [-1.35, 0.35]
    ax.plot(lim, lim, color="#999999", lw=0.7, ls="--", zorder=0)
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel("organoid dependency (median LFC, lineage-matched)")
    ax.set_ylabel("matched 2D cell-line dependency (median LFC)")
    ax.text(-1.33, 0.28, "stronger in organoids", fontsize=6.2, color=CRC, style="italic")
    ax.text(0.05, -1.28, "stronger in 2D lines", fontsize=6.2, color="#777777", style="italic")
    handles = [Patch(fc=CAT_COL["stronger_in_organoids"], label="organoid-enhanced (31)"),
               Patch(fc=CAT_COL["similar_or_not_resolved"], label="similar / unresolved (66)"),
               Patch(fc=CAT_COL["weaker_in_organoids"], edgecolor="#999999", label="weaker in organoids (25)")]
    ax.legend(handles=handles, loc="upper left", fontsize=6.0, frameon=False)
    ax.set_title("a", loc="left", fontsize=10, fontweight="bold")

def panel_4b(ax):
    s = read_tsv(os.path.join(T, "P2_3A_EXTERNAL_VALIDATION_SHORTLIST.tsv"))
    s = sorted(s, key=lambda r: (r["priority"] != "high", -num(r["organoid_vs_cell_line_difference"])))
    y = np.arange(len(s))[::-1]
    for yi, r in zip(y, s):
        oe = num(r["organoid_effect"]); ce = num(r["cell_line_effect"])
        hp = r["priority"] == "high"
        c = CRC if hp else "#8AAECF"
        ax.plot([oe, ce], [yi, yi], color="#CCCCCC", lw=1.0, zorder=1)
        ax.scatter([oe], [yi], s=30, facecolor=c, edgecolor="none", zorder=3)
        ax.scatter([ce], [yi], s=30, facecolor="white", edgecolor=c, linewidth=1.2, zorder=3)
        ax.text(-1.62, yi, r["gene"], fontsize=6.4, va="center",
                fontweight="bold" if hp else "normal")
    ax.set_yticks([])
    ax.set_xlim(-1.7, 0.45)
    ax.set_xticks([-1.5, -1.0, -0.5, 0, 0.5])
    ax.set_xlabel("median LFC (more-negative = stronger dependency)")
    handles = [Patch(fc=CRC, label="organoid (10/10 enhanced)"),
               Patch(fc="white", ec=CRC, label="matched 2D lines")]
    ax.legend(handles=handles, loc="lower right", fontsize=6.0)
    ax.set_title("b", loc="left", fontsize=10, fontweight="bold")

def panel_4c(ax):
    d = read_tsv(os.path.join(T, "P2_3A_DEPMAP_CANDIDATE_COMPARISON.tsv"))
    from collections import Counter
    order = ["not_flagged", "DepMap_inferred_common_essential", "Sanger_pan_cancer_core_fitness", "both_core_fitness_lists"]
    lab = {"not_flagged": "not flagged (74)", "DepMap_inferred_common_essential": "DepMap inferred (19)",
           "Sanger_pan_cancer_core_fitness": "Sanger core fitness (1)", "both_core_fitness_lists": "both lists (29)"}
    allc = Counter(r["common_essential_status"] for r in d)
    en = [r for r in d if r["comparison_category"] == "stronger_in_organoids"]
    enc = Counter(r["common_essential_status"] for r in en)
    x = np.arange(2)
    w = 0.5
    bottom = np.zeros(2)
    for k in order:
        v = [allc.get(k, 0), enc.get(k, 0)]
        ax.bar(x, v, w, bottom=bottom, color=FLAG_COL[k], label=lab[k])
        bottom += v
    ax.set_xticks(x)
    ax.set_xticklabels(["all 123 candidates", "organoid-enhanced (31)"], fontsize=6.6)
    ax.set_ylabel("candidates", fontsize=7)
    ax.set_ylim(0, 130)
    ax.legend(fontsize=5.4, loc="upper right", frameon=False)
    ax.text(0.0, 124, "shortlist: 10/10 not flagged common-essential", fontsize=6.2, color="#333333")
    ax.set_title("c", loc="left", fontsize=10, fontweight="bold")

def panel_4d(ax):
    d = read_tsv(os.path.join(T, "P2_3B_FINAL_CANDIDATE_RANKING.tsv"))
    d = sorted(d, key=lambda r: -num(r["evidence_score"]))
    y = np.arange(len(d))[::-1]
    for yi, r in zip(y, d):
        sc = num(r["evidence_score"])
        m = r["functional_module"]
        c = MOD_COL.get(m, GREY)
        ax.barh(yi, sc, height=0.62, color=c)
        ax.text(sc + 0.25, yi, f"{r['gene']}  ({r['final_priority']})", va="center", fontsize=6.2)
    ax.set_yticks([])
    ax.set_xlim(0, 14)
    ax.set_xlabel("integrated evidence score (frozen P2-3B ranking)")
    ax.set_title("d", loc="left", fontsize=10, fontweight="bold")

def assemble():
    W, H = 180 * MM, 140 * MM
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1.0, 0.6], width_ratios=[1.35, 1.0],
                          left=0.06, right=0.995, top=0.96, bottom=0.075, hspace=0.65, wspace=0.4)
    ax_a = fig.add_subplot(gs[0, :]); panel_4a(ax_a)
    ax_b = fig.add_subplot(gs[1, 0]); panel_4b(ax_b)
    ax_c = fig.add_subplot(gs[1, 1]); panel_4c(ax_c)
    ax_d = fig.add_subplot(gs[2, :]); panel_4d(ax_d)
    fig.text(0.06, 0.012, "organoid-enhanced comparison; not cell-line validation (matched 2D sets: 17 CRC / 5 ESCA lines)",
             fontsize=7, style="italic", color="#444444")
    save(fig, os.path.join(ASSEMBLED, "Figure4"))

def panels_only():
    specs = [("4a", panel_4a, (170 * MM, 55 * MM)),
             ("4b", panel_4b, (85 * MM, 45 * MM)),
             ("4c", panel_4c, (55 * MM, 45 * MM)),
             ("4d", panel_4d, (170 * MM, 30 * MM))]
    for name, fn, size in specs:
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([0.05, 0.06, 0.92, 0.9])
        fn(ax)
        save(fig, os.path.join(PANELS, f"Fig4_{name}"))

if __name__ == "__main__":
    panels_only()
    assemble()
    print("FIG4 DONE")
