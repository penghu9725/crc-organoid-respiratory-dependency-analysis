# -*- coding: utf-8 -*-
"""FIG 5: 5a TCGA module abundance T/N (COADREAD/COAD/READ), 5b FeS stage association."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import (save, read_tsv, num, MM, PANELS, ASSEMBLED, MOD_COL, MOD_LABEL,
                   ROBUST_MODULES, CRC, ESCA, GREY, BASE, TUM, NORM)

TC = os.path.join(os.path.abspath(os.environ.get("ORGANOID_P2_4A_TEMP", os.path.join(BASE, "data", "source", ".tmp_p2_4a"))), "tcga_scores.tsv")
MODS = ROBUST_MODULES
COMP_LABEL = {"COADREAD": "COADREAD combined", "COAD": "COAD", "READ": "READ"}

def _scores():
    import csv
    rows = list(csv.DictReader(open(TC, encoding="utf-8-sig"), delimiter="\t"))
    return rows

def _fmt_fdr(x):
    x = float(x)
    return "n.s." if x >= 0.05 else "{:.1e}".format(x)

def panel_5a(ax):
    rows = _scores()
    stats = {c: {} for c in ["COADREAD", "COAD", "READ"]}
for r in read_tsv(os.path.join(os.path.abspath(os.environ.get("ORGANOID_DERIVED_ROOT", os.path.join(BASE, "data", "derived", "analysis_outputs"))), "P2_4A_TCGA_MODULE_SUMMARY.tsv")):
        stats[r["comparison"]][r["module"]] = r
    comps = ["COADREAD", "COAD", "READ"]
    positions = {}
    xpos = 0
    for c in comps:
        positions[c] = {}
        for m in MODS:
            positions[c][m] = xpos
            xpos += 1
        xpos += 0.8
    xmax = xpos - 0.8
    data = {}
    for r in rows:
        m = r["module"]
        if m not in MODS:
            continue
        for comp in ["COADREAD", "COAD", "READ"]:
            ok = comp == "COADREAD" or r["project"] == comp
            if not ok:
                continue
            key = (comp, m, r["status"])
            data.setdefault(key, []).append(float(r["module_score"]))
    for comp in comps:
        for m in MODS:
            for si, status in enumerate(["Normal", "Tumor"]):
                key = (comp, m, status)
                if key not in data:
                    continue
                vals = np.array(data[key])
                x = positions[comp][m]
                w = 0.30
                xo = x - w / 2 if status == "Normal" else x + w / 2
                col = NORM if status == "Normal" else TUM
                vp = ax.violinplot([vals], [xo], widths=w, showextrema=False)
                for b in vp["bodies"]:
                    b.set_facecolor(col); b.set_alpha(0.55); b.set_edgecolor("none")
                bp = ax.boxplot([vals], positions=[xo], widths=w * 0.45, showfliers=False,
                                patch_artist=True, medianprops=dict(color="black", lw=0.6),
                                boxprops=dict(facecolor="white", edgecolor=col, lw=0.7),
                                whiskerprops=dict(color=col, lw=0.6))
    ax.axvline(4.9, color="#CCCCCC", lw=0.6)
    ax.axvline(10.7, color="#CCCCCC", lw=0.6)
    for j, c in enumerate(comps):
        grp = positions[c]["CI_structural"] + 2
        ax.text(grp, -4.2, COMP_LABEL[c], ha="center", fontsize=6.8, fontweight="bold")
        # FDR annotation under each module (from 15-test joint family, per comparison)
        for mi, m in enumerate(MODS):
            st = stats[c].get(m)
            if st is None:
                continue
            fdr = float(st["BH_FDR"])
            x = positions[c][m]
            ax.text(x, -2.35, "FDR " + _fmt_fdr(fdr), ha="center", fontsize=4.6,
                    color="#333333" if fdr < 0.05 else "#999999")
    ax.set_xticks([(positions["COADREAD"][m] + positions["READ"][m]) / 2 for m in MODS])
    ax.set_xticklabels([MOD_LABEL[m].replace("Fe-S biogenesis", "Fe-S").replace("Complex ", "C")
                        for m in MODS], fontsize=5.6, rotation=30, ha="right")
    ax.set_ylabel("within-sample percentile-rank module score\n(expression abundance)", fontsize=7)
    ax.set_xlim(-0.6, xmax + 1.6)
    handles = [Patch(fc=NORM, alpha=0.8, label="normal (n = 51)"),
               Patch(fc=TUM, alpha=0.8, label="tumour (n = 380)")]
    ax.legend(handles=handles, loc="upper right", fontsize=6.4)
    ax.text(xmax + 2.1, 0.0, "BH-FDR jointly across 15 module-by-\ncomparison tests (3 comparisons \u00d7 5 modules)",
            fontsize=5.6, color="#333333")
    ax.set_title("a", loc="left", fontsize=10, fontweight="bold")

def panel_5b(ax):
    rows = _scores()
    tum = [r for r in rows if r["status"] == "Tumor" and r["module"] == "FeS_biogenesis"
           and r.get("stage") and str(r["stage"]).strip() not in ("", "NA")]
    import re
    dat = []
    for r in tum:
        st = str(r["stage"]).strip()
        m = re.match(r"Stage ([IV]+)", st)
        if m:
            dat.append((m.group(1), float(r["module_score"])))
    order = ["I", "II", "III", "IV"]
    d = {o: [v for s, v in dat if s == o] for o in order}
    labels = [f"Stage {o}\n(n={len(d[o])})" for o in order]
    ax.boxplot([d[o] for o in order], showfliers=False, patch_artist=True,
               medianprops=dict(color="black", lw=0.7),
               boxprops=dict(facecolor=MOD_COL["FeS_biogenesis"], alpha=0.35, lw=0.7))
    for i, o in enumerate(order, 1):
        ax.scatter(np.full(len(d[o]), i) + np.random.RandomState(i).uniform(-0.08, 0.08, len(d[o])),
                   d[o], s=5, color=MOD_COL["FeS_biogenesis"], alpha=0.5, linewidths=0)
        ax.text(i, -0.35, f"n={len(d[o])}", ha="center", fontsize=5.8)
    ax.set_xticks(range(1, 5))
    ax.set_xticklabels(["I", "II", "III", "IV"], fontsize=7)
    ax.set_xlabel("pathologic stage", fontsize=7)
    ax.set_ylabel("Fe-S biogenesis module score", fontsize=7)
    ax.text(0.6, max(v for vals in d.values() for v in vals) * 1.02,
            "ordinal Spearman \u03c1 = \u22120.196, FDR 9.8\u00d710\u207b\u2074\nKruskal\u2013Wallis \u03b5\u00b2 = 0.041, FDR 2.8\u00d710\u207b\u00b3",
            fontsize=6.2, color="#333333")
    ax.set_title("b", loc="left", fontsize=10, fontweight="bold")

def assemble():
    W, H = 180 * MM, 110 * MM
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.8, 1.0], width_ratios=[2.2, 1.0],
                          left=0.07, right=0.99, top=0.96, bottom=0.08, hspace=0.55, wspace=0.35)
    ax_a = fig.add_subplot(gs[0, :]); panel_5a(ax_a)
    ax_b = fig.add_subplot(gs[1, 1]); panel_5b(ax_b)
    fig.text(0.07, 0.015, "expression programs are abundance measures; not CRISPR dependency",
             fontsize=7, style="italic", color="#444444")
    save(fig, os.path.join(ASSEMBLED, "Figure5"))

def panels_only():
    specs = [("5a", panel_5a, (175 * MM, 60 * MM)),
             ("5b", panel_5b, (85 * MM, 45 * MM))]
    for name, fn, size in specs:
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([0.05, 0.06, 0.93, 0.9])
        fn(ax)
        save(fig, os.path.join(PANELS, f"Fig5_{name}"))

if __name__ == "__main__":
    panels_only()
    assemble()
    print("FIG5 DONE")
