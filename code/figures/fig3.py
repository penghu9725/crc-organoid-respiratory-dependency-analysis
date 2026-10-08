# -*- coding: utf-8 -*-
"""FIG 3 (core): 3a schematic, 3b module effects, 3c CI/CIV heatmap, 3d drivers, 3e sensitivity."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
from style import (save, grayscale, read_tsv, num, MM, PANELS, ASSEMBLED, PREVIEW, MOD_COL, MOD_LABEL,
                   ROBUST_MODULES, CRC, ESCA, GREY, MOD_MTFAS, BASE)

T = os.path.abspath(os.environ.get("ORGANOID_DERIVED_ROOT", os.path.join(BASE, "data", "derived", "analysis_outputs")))
METADATA_ROOT = os.path.abspath(os.environ.get("ORGANOID_METADATA_ROOT", os.path.join(BASE, "data", "external_metadata")))
SOURCE_ROOT = os.path.abspath(os.environ.get("ORGANOID_SOURCE_ROOT", os.path.join(BASE, "data", "source")))
MAP = os.path.join(METADATA_ROOT, "OXPHOS_SUBMODULE_MAP.tsv")
CONS = os.path.join(SOURCE_ROOT, "02_crispr", "METABOLIC_ORGANOID_GENE_CONSENSUS.tsv")

# ============ 3a schematic ============
def panel_3a(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
    cx = [1.6, 4.2, 6.8, 9.4]  # complex positions (skip 5th slot handled below)
    # complexes I-IV large, V tucked; use two rows: row1 complexes I-IV? simpler: row of 5
    pos = {"CI": 1.35, "CII": 3.55, "CIII": 5.75, "CIV": 7.95}
    col = {"CI": MOD_COL["CI_structural"], "CII": GREY, "CIII": GREY, "CIV": MOD_COL["CIV_structural"]}
    # structural block (holoenzyme) and assembly block beneath
    y_str, y_asm = 3.4, 2.15
    for name, x in pos.items():
        c = col[name]
        ax.add_patch(FancyBboxPatch((x, y_str), 1.25, 1.5, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    fc=c, ec="none", alpha=0.85))
        ax.text(x + 0.625, y_str + 0.75, name, ha="center", va="center", fontsize=8,
                fontweight="bold", color="white")
        ax.add_patch(FancyBboxPatch((x + 0.12, y_asm), 1.0, 0.5, boxstyle="round,pad=0.01,rounding_size=0.05",
                                    fc="white", ec=c, lw=1.0))
    ax.text(1.95, 1.62, "assembly", ha="center", fontsize=6.5, color="#333333")
    ax.text(4.15, 1.62, "assembly", ha="center", fontsize=6.5, color="#333333")
    ax.text(6.35, 1.62, "assembly", ha="center", fontsize=6.5, color="#333333")
    ax.text(8.55, 1.62, "assembly", ha="center", fontsize=6.5, color="#333333")
    # supporting rows
    ax.text(0.05, 4.4, "respiratory chain", fontsize=7, style="italic", color="#555555")
    ax.text(0.05, 0.9, "supporting biosynthesis", fontsize=7, style="italic", color="#555555")
    sup = [("Fe-S biogenesis", MOD_COL["FeS_biogenesis"], "ROBUST"),
           ("mtFAS", MOD_MTFAS, "exploratory"),
           ("CoQ / mtDNA", GREY, "not supported")]
    x = 1.0
    for label, c, tag in sup:
        w = 2.2 if label != "CoQ / mtDNA" else 1.9
        fc = c if tag != "exploratory" else "white"
        ax.add_patch(FancyBboxPatch((x, 0.15), w, 0.6, boxstyle="round,pad=0.02,rounding_size=0.05",
                                    fc=fc, ec=c if tag != "exploratory" else MOD_MTFAS, lw=1.2,
                                    hatch="" if tag != "exploratory" else "///"))
        ax.text(x + w / 2, 0.45, label, ha="center", va="center", fontsize=6.8)
        x += w + 0.55
    ax.text(9.4, 4.75, "13 prespecified submodules", fontsize=7.5, ha="right", style="italic", color="#555555")
    ax.text(9.4, 4.35, "mtDNA-encoded genes excluded", fontsize=6.5, ha="right", color="#777777")

# ============ 3b module effects ============
def panel_3b(ax):
    d = read_tsv(os.path.join(T, "P2_3B_OXPHOS_SUBMODULE_RESULTS.tsv"))
    order = sorted(d, key=lambda r: num(r["effect_size_rank_biserial"]))
    y = np.arange(len(order))
    for r, yi in zip(order, y):
        mod = r["module"]
        eff = num(r["effect_size_rank_biserial"])
        lo = num(r["CI95_low"]); hi = num(r["CI95_high"])
        fdr = num(r["BH_FDR"]); rsfdr = num(r["rank_sensitivity_BH_FDR"])
        robust = mod in ROBUST_MODULES
        if mod == "mtFAS":
            fc, ec, hatch = "white", MOD_MTFAS, "///"
        elif robust:
            fc, ec, hatch = MOD_COL[mod], "none", None
        else:
            fc, ec, hatch = GREY, "none", None
        ax.plot([lo, hi], [yi, yi], color=fc if robust else "#888888", lw=1.6, zorder=2)
        ax.scatter([eff], [yi], s=26 if robust else 20, facecolor=fc, edgecolor=ec,
                   hatch=hatch, zorder=3, linewidth=0.8)
        label = MOD_LABEL.get(mod, mod)
        if mod == "mtFAS":
            label = "mtFAS (exploratory)"
        ax.text(-0.62, yi, label, ha="right", va="center", fontsize=6.6,
                color="#111111" if robust else "#555555")
        tag = "FDR {:.1e}".format(fdr) if robust else ("FDR {:.1e}*".format(fdr) if mod == "mtFAS" else "")
        ax.text(0.52, yi, tag, va="center", fontsize=6.2,
                color=MOD_COL[mod] if robust else ("#B8860B" if mod == "mtFAS" else "#777777"))
    ax.axvline(0, color="#333333", lw=0.7, ls="-")
    ax.set_yticks([])
    ax.set_xlim(-0.7, 0.7)
    ax.set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
    ax.set_xlabel("CRC-vs-ESCA rank-biserial effect\n(positive = stronger in CRC)")
    ax.set_title("b", loc="left", fontsize=10, fontweight="bold")
    # legend
    handles = [Line2D([], [], marker="o", ls="", ms=5, mfc=MOD_COL["CI_assembly"], mec="none", label="robust (FDR < 0.05)"),
               Line2D([], [], marker="o", ls="", ms=5, mfc="white", mec=MOD_MTFAS, label="exploratory (mtFAS)"),
               Line2D([], [], marker="o", ls="", ms=5, mfc=GREY, mec="none", label="not supported")]
    ax.legend(handles=handles, loc="lower right", fontsize=6.2, bbox_to_anchor=(1.0, -0.02))

# ============ 3c heatmap ============
def panel_3c(ax):
    import csv
    mods = ["CI_structural", "CI_assembly", "CIV_structural", "CIV_assembly"]
    gmap = {}
    with open(MAP, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["submodule_primary"] in mods:
                gmap.setdefault(r["submodule_primary"], []).append(r["gene_symbol"])
    keep = {}
    for m in mods:
        for g in gmap[m]:
            keep[g] = m
    # stream consensus: gene x sample adjusted LFC, CRC + ESCA only
    samples = {}
    vals = {}
    with open(CONS, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["cancer_type"] not in ("Colorectal", "Oesophageal"):
                continue
            g = r["gene"]
            if g not in keep:
                continue
            sid = r["sample_ID"]
            samples[sid] = r["cancer_type"]
            try:
                v = float(r["adjusted_consensus_LFC"])
            except ValueError:
                v = np.nan
            vals.setdefault(g, {})[sid] = v
    genes = [g for g, m in [(gg, mm) for mm in mods for gg in gmap[mm]] if g in vals]
    sid_order = sorted(samples, key=lambda s: (samples[s], s))
    X = np.full((len(genes), len(sid_order)), np.nan)
    for i, g in enumerate(genes):
        row = vals[g]
        for j, s in enumerate(sid_order):
            X[i, j] = row.get(s, np.nan)
    # gene-wise standardisation of -LFC (dependency direction) across organoids
    Y = -X
    mu = np.nanmean(Y, axis=1, keepdims=True); sd = np.nanstd(Y, axis=1, keepdims=True)
    sd[sd == 0] = 1
    Z = (Y - mu) / sd
    # order genes within module by CRC-mean minus ESCA-mean
    crc_cols = [j for j, s in enumerate(sid_order) if samples[s] == "Colorectal"]
    esc_cols = [j for j, s in enumerate(sid_order) if samples[s] == "Oesophageal"]
    new_order = []
    for m in mods:
        idx = [i for i, g in enumerate(genes) if keep[g] == m]
        dme = [np.nanmean(Z[i, crc_cols]) - np.nanmean(Z[i, esc_cols]) for i in idx]
        idx = [i for _, i in sorted(zip(dme, idx), reverse=True)]
        new_order += idx
    genes = [genes[i] for i in new_order]
    Z = Z[new_order, :]
    mod_of = [keep[g] for g in genes]
    im = ax.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                   interpolation="nearest", rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    # module strips
    prev = None
    ypos = []
    for i, m in enumerate(mod_of):
        if m != prev:
            ypos.append(i)
            prev = m
    ypos.append(len(genes))
    for k in range(len(ypos) - 1):
        y0, y1 = ypos[k], ypos[k + 1]
        ax.add_patch(Rectangle((-2.2, y0 - 0.5), 0.8, y1 - y0, fc=MOD_COL[mod_of[y0]], ec="none"))
        ax.text(-1.25, (y0 + y1) / 2 - 0.5, MOD_LABEL[mod_of[y0]].replace("Complex ", "C").replace("Fe-S ", "Fe-S\n"),
                ha="center", va="center", fontsize=5.2, rotation=90)
    # primary-driver labels on right
    drivers = {"NDUFA6": 1, "NDUFAF8": 1, "NDUFAF7": 1, "NDUFAF3": 1, "TIMMDC1": 1,
               "COX6B1": 1, "COX6A1": 1, "SCO1": 1}
    for i, g in enumerate(genes):
        if g in drivers:
            ax.text(Z.shape[1] + 0.3, i, g, fontsize=5.0, va="center", color="#111111")
    ax.set_xlim(-3.0, Z.shape[1] + 3.5)
    # CRC / ESCA spans
    n_crc = len(crc_cols); n_esc = len(esc_cols)
    ax.text((n_crc / 2) - 0.5, -1.3, f"CRC organoids (n = {n_crc})", ha="center", fontsize=6)
    ax.text(n_crc + n_esc / 2 - 0.5, -1.3, f"ESCA (n = {n_esc})", ha="center", fontsize=6)
    ax.axvline(n_crc - 0.5, color="#333333", lw=0.6)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.02, pad=0.30)
    cb.ax.tick_params(labelsize=5.5)
    cb.set_label("standardised dependency (-LFC z)", fontsize=5.8, labelpad=1)
    ax.set_title("c", loc="left", fontsize=10, fontweight="bold")

# ============ 3d driver structure ============
def panel_3d(ax):
    d = read_tsv(os.path.join(T, "P2_3A_OXPHOS_DRIVER_GENES.tsv"))
    prim = [r for r in d if r["driver_category"] == "primary_FDR_driver"]
    prim = sorted(prim, key=lambda r: -num(r["rank_biserial_CRC_vs_ESCA"]))
    n_support = sum(1 for r in d if r["driver_category"] == "supporting_aligned_driver")
    n_oppose = sum(1 for r in d if r["driver_category"] == "opposing_gene")
    n_weak = sum(1 for r in d if r["driver_category"] in ("weak_or_neutral", ""))
    ax.set_xlim(0, 10); ax.set_ylim(-4.5, 12); ax.axis("off")
    ax.text(0.1, 11.2, "primary FDR drivers (n = 9)", fontsize=7.5, fontweight="bold")
    for i, r in enumerate(prim):
        x = 0.55 + i * 1.02
        rb = num(r["rank_biserial_CRC_vs_ESCA"])
        ax.add_patch(Rectangle((x - 0.3, 0), 0.6, -rb * 8, fc=MOD_COL.get(r["submodule"], GREY),
                               ec="none"))
        ax.text(x, -rb * 8 - 0.5, r["gene"], ha="center", fontsize=6.6, fontweight="bold",
                rotation=0)
    ax.text(0.1, -1.2, "gene-level rank-biserial (CRC-stronger)", fontsize=6.2, color="#555555")
    counts = [("supporting aligned", n_support, "#B0B0B0"), ("opposing", n_oppose, "#D0D0D0"),
              ("weak / neutral", n_weak, "#E5E5E5")]
    x = 1.0
    for label, n, c in counts:
        ax.add_patch(Rectangle((x, -3.0), 2.3, 0.7, fc=c, ec="#888888", lw=0.5))
        ax.text(x + 1.15, -2.65, f"{label}: {n}", ha="center", va="center", fontsize=6.2)
        x += 2.75
    ax.text(0.1, -3.9, "module colour = functional membership; ISCA2 = Fe-S, TIMMDC1 = CI assembly/transport",
            fontsize=5.8, color="#666666")
    ax.set_title("d", loc="left", fontsize=10, fontweight="bold")

# ============ 3e sensitivity ============
def panel_3e(ax):
    d = read_tsv(os.path.join(T, "P2_3B_OXPHOS_SUBMODULE_RESULTS.tsv"))
    ax.axhline(0.05, color="#999999", lw=0.7, ls="--")
    ax.axvline(0.05, color="#999999", lw=0.7, ls="--")
    ax.text(0.06, 0.6, "FDR = 0.05", fontsize=5.6, color="#777777")
    for r in d:
        mod = r["module"]
        x = num(r["BH_FDR"]); y = num(r["rank_sensitivity_BH_FDR"])
        robust = mod in ROBUST_MODULES
        if mod == "mtFAS":
            ax.scatter([x], [y], s=30, facecolor="white", edgecolor=MOD_MTFAS, hatch="///",
                       linewidth=0.9, zorder=4)
            ax.annotate("mtFAS\n(exploratory)", (x, y), textcoords="offset points",
                        xytext=(6, -6), fontsize=5.6, color="#8a7a00")
        elif robust:
            ax.scatter([x], [y], s=32, facecolor=MOD_COL[mod], edgecolor="none", zorder=4)
        else:
            ax.scatter([x], [y], s=22, facecolor=GREY, edgecolor="none", alpha=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("primary module BH-FDR (13-module family)")
    ax.set_ylabel("fractional-rank sensitivity BH-FDR")
    ax.set_xlim(3e-5, 3); ax.set_ylim(3e-5, 3)
    ax.set_title("e", loc="left", fontsize=10, fontweight="bold")

# ============ assemble Fig3 ============
def assemble():
    W, H = 180 * MM, 175 * MM
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.15, 0.95], width_ratios=[1.0, 1.35],
                          left=0.05, right=0.99, top=0.97, bottom=0.06, hspace=0.55, wspace=0.35)
    ax_a = fig.add_subplot(gs[0, :]); panel_3a(ax_a)
    ax_b = fig.add_subplot(gs[1, 0]); panel_3b(ax_b)
    ax_c = fig.add_subplot(gs[1, 1]); panel_3c(ax_c)
    ax_d = fig.add_subplot(gs[2, 0]); panel_3d(ax_d)
    ax_e = fig.add_subplot(gs[2, 1]); panel_3e(ax_e)
    save(fig, os.path.join(ASSEMBLED, "Figure3"))

def panels_only():
    specs = [("3a", panel_3a, (180 * MM, 34 * MM)),
             ("3b", panel_3b, (62 * MM, 78 * MM)),
             ("3c", panel_3c, (110 * MM, 80 * MM)),
             ("3d", panel_3d, (62 * MM, 58 * MM)),
             ("3e", panel_3e, (62 * MM, 58 * MM))]
    for name, fn, size in specs:
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
        fn(ax)
        save(fig, os.path.join(PANELS, f"Fig3_{name}"))

if __name__ == "__main__":
    panels_only()
    assemble()
    grayscale(os.path.join(ASSEMBLED, "Figure3.png"), os.path.join(PREVIEW, "Figure3_gray.png"))
    print("FIG3 DONE")
