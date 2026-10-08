# -*- coding: utf-8 -*-
"""FIG 6: 6a cell-type programs, 6b paired epithelial n=6 (EXPLORATORY)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import (save, read_tsv, num, MM, PANELS, ASSEMBLED, MOD_COL, MOD_LABEL,
                   ROBUST_MODULES, BASE, TUM, NORM)

T = os.path.abspath(os.environ.get("ORGANOID_DERIVED_ROOT", os.path.join(BASE, "data", "derived", "analysis_outputs")))

def panel_6a(ax):
    d = read_tsv(os.path.join(T, "P2_4A_SC_CELLTYPE_PROGRAM.tsv"))
    ct_ord = ["Tumour", "Normal"]
    types = sorted({r["Cell_type"] for r in d})
    # order cell types: tumour/normal epithelial first then others by median
    def keyf(ct):
        if ct == "Epithelial cells":
            return 0
        return 1
    cts = sorted(types, key=keyf)
    # prefer Class='Tumour' for 'Epithelial cells' display and other classes combined? use per Class+Cell_type rows
    cells = []  # (row_label, module, score)
    for r in d:
        cl = "Tumour" if r["Class"] == "Tumor" else r["Class"]; ct = r["Cell_type"]
        if ct == "Epithelial cells":
            lab = ("tumour\nepithelial" if cl == "Tumour" else "normal\nepithelial")
        else:
            lab = ct.replace(" cells", "").replace("_", " ")
            if cl == "Tumour" and ct not in ("Epithelial cells",):
                continue  # non-epithelial shown from normal class only (stromal/immune)
        cells.append((lab, r["module"], num(r["median_program_score"])))
    row_names = []
    for r in d:
        cl = "Tumour" if r["Class"] == "Tumor" else r["Class"]; ct = r["Cell_type"]
        if ct == "Epithelial cells":
            row_names.append(("tumour\nepithelial" if cl == "Tumour" else "normal\nepithelial"))
        elif cl == "Normal":
            row_names.append(ct.replace(" cells", "").replace("_", " "))
    uniq = []
    for x in row_names:
        if x not in uniq:
            uniq.append(x)
    # explicit row order: tumour/normal epithelial first, then normal-class immune/stromal types (stable)
    def row_key(x):
        if x == "tumour\nepithelial":
            return 0
        if x == "normal\nepithelial":
            return 1
        return 2
    uniq = sorted(uniq, key=lambda x: row_key(x))
    mods = ROBUST_MODULES
    M = np.full((len(uniq), len(mods)), np.nan)
    for lab, m, sc in cells:
        if lab in uniq and m in mods:
            M[uniq.index(lab), mods.index(m)] = sc
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto", rasterized=True)
    ax.set_yticks(range(len(uniq))); ax.set_yticklabels(uniq, fontsize=6.4)
    ax.set_xticks(range(len(mods)))
    ax.set_xticklabels([MOD_LABEL[m].replace("Fe-S biogenesis", "Fe-S").replace("Complex ", "C")
                        for m in mods], fontsize=5.8, rotation=30, ha="right")
    for i in range(len(uniq)):
        for j in range(len(mods)):
            v = M[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=5.0,
                        color="white" if abs(v) > 0.45 else "#333333")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cb.ax.tick_params(labelsize=5.5)
    cb.set_label("median patient-level program score", fontsize=5.8)
    ax.set_title("a", loc="left", fontsize=10, fontweight="bold")

def panel_6b(ax):
    d = read_tsv(os.path.join(T, "P2_4A_SC_PSEUDOBULK.tsv"))
    epi = [r for r in d if r["row_type"] == "patient_score"]
    tests = {r["module"]: r for r in d if r["row_type"] == "paired_test"}
    pats = sorted({r["Patient"] for r in epi if r["Class"] in ("Tumour", "Normal")})
    paired = set()
    for r in epi:
        if r["Class"] == "Tumour":
            paired.add(r["Patient"])
    mods = ROBUST_MODULES
    # collect only patients present in both classes
    both = []
    for p in pats:
        cls = {r["Class"] for r in epi if r["Patient"] == p}
        if cls >= {"Tumour", "Normal"}:
            both.append(p)
    both = sorted(both)
    n = len(both)
    for mi, m in enumerate(mods):
        axs = ax if len(mods) == 1 else None
    # draw 5 mini facets on one axis by offset; simpler: single axis grouped
    xs = []
    for mi, m in enumerate(mods):
        for pi, p in enumerate(both):
            vals = {r["Class"]: num(r["module_score"]) for r in epi
                    if r["Patient"] == p and r["module"] == m and r["Class"] in ("Tumour", "Normal")}
            if "Tumour" in vals and "Normal" in vals:
                x0 = mi * 4 + 0.4
                ax.plot([x0, x0 + 0.8], [vals["Normal"], vals["Tumour"]], color="#CCCCCC",
                        lw=0.7, zorder=1)
                ax.scatter([x0], [vals["Normal"]], s=14, facecolor=NORM, edgecolor="none", zorder=2)
                ax.scatter([x0 + 0.8], [vals["Tumour"]], s=14, facecolor=TUM, edgecolor="none", zorder=2)
        t = tests.get(m)
        if t:
            fdr = float(t["BH_FDR"]); rawp = float(t["raw_P"])
            ax.text(mi * 4 + 0.6, -4.6, "raw P {:.3f}\nFDR {:.3f}".format(rawp, fdr),
                    ha="center", fontsize=5.6, color="#333333")
    ax.set_xticks([mi * 4 + 0.6 for mi in range(len(mods))])
    ax.set_xticklabels([MOD_LABEL[m].replace("Fe-S biogenesis", "Fe-S").replace("Complex ", "C")
                        for m in mods], fontsize=6.0, rotation=25, ha="right")
    ax.set_ylabel("patient-level pseudobulk module score", fontsize=7)
    ax.axhline(0, color="#CCCCCC", lw=0.5)
    ax.text(0.4, 4.6, f"paired patients n = {n}  (tumour vs normal epithelial; \u2265100 cells each)",
            fontsize=6.4)
    handles = [Patch(fc=NORM, label="normal epithelial"), Patch(fc=TUM, label="tumour epithelial")]
    ax.legend(handles=handles, loc="upper right", fontsize=6.0)
    # big exploratory badge
    ax.text(0.99, 0.98, "EXPLORATORY", transform=ax.transAxes, ha="right", va="top",
            fontsize=9, fontweight="bold", color="#B8860B",
            bbox=dict(boxstyle="round,pad=0.35", fc="#FFF8DC", ec="#B8860B", lw=1.2))
    ax.set_title("b", loc="left", fontsize=10, fontweight="bold")
    ax.set_xlim(-1, len(mods) * 4 + 1)

def assemble():
    W, H = 89 * MM, 130 * MM
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.15], left=0.16, right=0.97,
                          top=0.97, bottom=0.07, hspace=0.5)
    ax_a = fig.add_subplot(gs[0]); panel_6a(ax_a)
    ax_b = fig.add_subplot(gs[1]); panel_6b(ax_b)
    fig.text(0.16, 0.012, "directional / exploratory; not single-cell validation (cells never treated as replicates)",
             fontsize=6.4, style="italic", color="#444444")
    save(fig, os.path.join(ASSEMBLED, "Figure6"))

def panels_only():
    specs = [("6a", panel_6a, (85 * MM, 55 * MM)),
             ("6b", panel_6b, (85 * MM, 62 * MM))]
    for name, fn, size in specs:
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([0.06, 0.06, 0.92, 0.9])
        fn(ax)
        save(fig, os.path.join(PANELS, f"Fig6_{name}"))

if __name__ == "__main__":
    panels_only()
    assemble()
    print("FIG6 DONE")
