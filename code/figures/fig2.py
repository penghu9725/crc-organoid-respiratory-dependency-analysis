# -*- coding: utf-8 -*-
"""FIG 2: 2a landscape heatmap, 2b pathway forest, 2c volcano, 2d strength-MAD, 2e OXPHOS strip."""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from style import (save, read_tsv, num, MM, PANELS, ASSEMBLED, MOD_COL, MOD_LABEL,
                   CRC, ESCA, GREY, BASE, NEUTRALS)

T = os.path.abspath(os.environ.get("ORGANOID_DERIVED_ROOT", os.path.join(BASE, "data", "derived", "analysis_outputs")))
SOURCE_ROOT = os.path.abspath(os.environ.get("ORGANOID_SOURCE_ROOT", os.path.join(BASE, "data", "source")))
METADATA_ROOT = os.path.abspath(os.environ.get("ORGANOID_METADATA_ROOT", os.path.join(BASE, "data", "external_metadata")))
CONS = os.path.join(SOURCE_ROOT, "02_crispr", "METABOLIC_ORGANOID_GENE_CONSENSUS.tsv")
UNI = os.path.join(METADATA_ROOT, "METABOLIC_GENE_UNIVERSE_TIERED.tsv")
META = os.path.join(METADATA_ROOT, "ORGANOID_METADATA.tsv")
LI_COL = {"Colorectal": CRC, "Oesophageal": ESCA, "Ovarian": NEUTRALS["Ovarian"],
          "Pancreatic": NEUTRALS["Pancreatic"], "Gastric": NEUTRALS["Gastric"]}

# ============ 2a heatmap ============
def _load_matrix():
    # gene -> pathway (first occurrence in tiered universe), sample -> lineage
    gp = {}
    with open(UNI, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["tier"] in ("Tier 1", "Tier 2") and r["gene_symbol"] not in gp:
                gp[r["gene_symbol"]] = r["pathway_major"]
    sl = {}
    with open(META, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sl[r["sample_ID"].strip()] = r["primary_tumour_type"].strip()
    genes, sid, X = [], [], []
    vals = {}
    with open(CONS, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["gene"]
            if g not in gp:
                continue
            s = r["sample_ID"]
            try:
                v = float(r["adjusted_consensus_LFC"])
            except ValueError:
                v = np.nan
            vals.setdefault(g, {})[s] = v
    genes = list(vals.keys())
    sid = sorted(vals[genes[0]].keys(), key=lambda s: (sl.get(s, "Z"), s))
    X = np.full((len(genes), len(sid)), np.nan)
    for i, g in enumerate(genes):
        for j, s in enumerate(sid):
            X[i, j] = vals[g].get(s, np.nan)
    Y = -X
    mu = np.nanmean(Y, axis=1, keepdims=True); sd = np.nanstd(Y, axis=1, keepdims=True)
    sd[sd == 0] = 1
    return (Y - mu) / sd, genes, sid, gp, sl

def panel_2a(ax):
    Z, genes, sid, gp, sl = _load_matrix()
    pways = sorted({gp[g] for g in genes})
    pw_idx = {p: i for i, p in enumerate(pways)}
    order = sorted(range(len(genes)), key=lambda i: pw_idx[gp[genes[i]]])
    Z = Z[order, :]
    gord = [genes[i] for i in order]
    # pathway side strip (grayscale ramp)
    ramp = plt.cm.Greys(np.linspace(0.25, 0.9, len(pways)))
    strips = np.array([pw_idx[gp[g]] for g in gord]) / max(1, len(pways) - 1)
    ax.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
              interpolation="nearest", rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    # lineage colour strip on top
    lin = [sl.get(s, "Other") for s in sid]
    bounds = []
    prev = None
    for j, l in enumerate(lin):
        if l != prev:
            bounds.append(j); prev = l
    bounds.append(len(sid))
    for k in range(len(bounds) - 1):
        ax.add_patch(Rectangle((bounds[k] - 0.5, -2.6), bounds[k + 1] - bounds[k], 2.4,
                               fc=LI_COL.get(lin[bounds[k]], "#999999"), ec="none"))
        n = bounds[k + 1] - bounds[k]
        ax.text((bounds[k] + bounds[k + 1]) / 2 - 0.5, -1.4, f"{lin[bounds[k]]}\n(n={n})",
                ha="center", va="center", fontsize=4.6, color="white")
    cb = ax.figure.colorbar(ax.images[0], ax=ax, fraction=0.02, pad=0.12)
    cb.ax.tick_params(labelsize=5.5)
    cb.set_label("dependency z (-LFC)", fontsize=5.8)
    ax.text(0, len(genes) + 6, "1,412 Tier-1/Tier-2 metabolic genes \u00b7 rows ordered by pathway membership",
            fontsize=6.0, color="#444444")
    ax.set_title("a", loc="left", fontsize=10, fontweight="bold")
    ax.set_xlim(-0.5, len(sid) + 0.5); ax.set_ylim(len(genes) + 14, -4)

# ============ 2b pathway forest ============
def panel_2b(ax):
    d = read_tsv(os.path.join(T, "CRC_ESCA_PATHWAY_DIFFERENTIAL_DEPENDENCY.tsv"))
    order = sorted(d, key=lambda r: num(r["raw_P"]))
    y = np.arange(len(order))[::-1]
    for yi, r in zip(y, order):
        hl = num(r["hodges_lehmann_shift"]) if r.get("hodges_lehmann_shift") else num(r["median_difference_CRC_minus_ESCA"])
        lo = num(r["CI95_low"]); hi = num(r["CI95_high"])
        fdr = num(r["BH_FDR"])
        ox = r["pathway_major"] == "OXPHOS"
        c = CRC if ox else GREY
        ax.plot([lo, hi], [yi, yi], color=c, lw=1.5, alpha=0.9 if ox else 0.75)
        ax.scatter([hl], [yi], s=16, color="white", edgecolor=c, linewidth=1.0 if ox else 0.7, zorder=3)
        nm = r["pathway_major"].replace("_", " ").lower()
        ax.text(-0.24, yi, nm, ha="right", va="center", fontsize=5.6,
                color="#111111" if ox else "#666666")
        if ox:
            ax.text(0.015, yi + 0.06, "BH-FDR 6.2\u00d710\u207b\u2074  \u00b7  HL 0.098 (0.052\u20130.143)",
                    fontsize=5.8, color=CRC)
    ax.axvline(0, color="#333333", lw=0.7)
    ax.set_yticks([])
    ax.set_xlim(-0.26, 0.22)
    ax.set_xticks([-0.2, -0.1, 0, 0.1, 0.2])
    ax.set_xlabel("Hodges\u2013Lehmann shift (CRC \u2212 ESCA)\nwith 95% CI")
    ax.set_title("b", loc="left", fontsize=10, fontweight="bold")

# ============ 2c volcano ============
def panel_2c(ax):
    d = read_tsv(os.path.join(T, "CRC_ESCA_GENE_DIFFERENTIAL_DEPENDENCY.tsv"))
    labels9 = {"NDUFA6", "NDUFAF8", "NDUFAF7", "NDUFAF3", "TIMMDC1", "COX6B1", "COX6A1", "SCO1", "ISCA2"}
    lab2 = {"TYMS", "DHFR", "MSMO1"}
    xs, ys, cs = [], [], []
    for r in d:
        rb = num(r["rank_biserial_CRC_vs_ESCA"]); fdr = num(r["BH_FDR"])
        yv = -np.log10(max(fdr, 1e-300))
        if fdr < 0.05:
            c = CRC if rb < 0 else ESCA
        else:
            c = "#BBBBBB"
        ax.scatter([rb], [yv], s=4.5 if fdr >= 0.05 else 9, c=c, alpha=0.85, linewidths=0)
        if r["gene"] in labels9:
            ax.annotate(r["gene"], (rb, yv), fontsize=5.6, xytext=(3, 2),
                        textcoords="offset points", color="#111111", fontweight="bold")
        elif r["gene"] in lab2:
            ax.annotate(r["gene"], (rb, yv), fontsize=5.4, xytext=(3, 2),
                        textcoords="offset points", color="#555555")
    ax.axhline(-np.log10(0.05), color="#999999", lw=0.6, ls="--")
    ax.set_xlabel("rank-biserial (negative = CRC-stronger)")
    ax.set_ylabel("\u2212log\u2081\u2080 BH-FDR (1,412-gene family)")
    ax.set_title("c", loc="left", fontsize=10, fontweight="bold")
    handles = [Patch(fc=CRC, label="CRC-leaning, FDR < 0.05 (23)"),
               Patch(fc=ESCA, label="ESCA-leaning, FDR < 0.05 (10)"),
               Patch(fc="#BBBBBB", label="not significant")]
    ax.legend(handles=handles, loc="upper left", fontsize=5.8, frameon=False)

# ============ 2d strength-MAD ============
def panel_2d(ax):
    d = read_tsv(os.path.join(T, "SELECTIVE_METABOLIC_DEPENDENCIES.tsv"))
    cls = {"A_common_metabolic_essential": ("#7A7A7A", "A (44)"),
           "B_selective_metabolic_dependency": (CRC, "B (177)"),
           "C_weak_or_non_dependent": ("#DDDDDD", "C (1191)")}
    for r in d:
        c, lab = cls.get(r["dependency_class"], ("#DDDDDD", "C"))
        ax.scatter([-num(r["median_LFC"])], [num(r["MAD"])], s=4, c=c, alpha=0.75, linewidths=0)
    ax.set_xlabel("dependency strength (\u2212median LFC)")
    ax.set_ylabel("dispersion (MAD)")
    handles = [Patch(fc=c, label=lab) for c, lab in cls.values()]
    ax.legend(handles=handles, loc="upper left", fontsize=5.8, frameon=False)
    ax.set_title("d", loc="left", fontsize=10, fontweight="bold")

# ============ 2e OXPHOS gene strip ============
def panel_2e(ax):
    d = read_tsv(os.path.join(T, "CRC_ESCA_GENE_DIFFERENTIAL_DEPENDENCY.tsv"))
    drv = read_tsv(os.path.join(T, "P2_3A_OXPHOS_DRIVER_GENES.tsv"))
    cat = {r["gene"]: r["driver_category"] for r in drv}
    sub = {r["gene"]: r["submodule"] for r in drv}
    ox = [r for r in d if r["gene"] in cat or r["gene"] in sub]
    # order: driver categories then by effect
    def key(r):
        g = r["gene"]
        k = {"primary_FDR_driver": 0, "supporting_aligned_driver": 1,
             "opposing_gene": 2}.get(cat.get(g, "x"), 3)
        return (k, num(r["rank_biserial_CRC_vs_ESCA"]))
    ox = sorted(ox, key=key)
    y = np.arange(len(ox))
    for yi, r in zip(y, ox):
        g = r["gene"]; rb = num(r["rank_biserial_CRC_vs_ESCA"]); fdr = num(r["BH_FDR"])
        m = sub.get(g)
        c = MOD_COL.get(m, "#BBBBBB") if cat.get(g) != "opposing_gene" else "#999999"
        ax.plot([0, rb], [yi, yi], color="#DDDDDD", lw=0.5, zorder=1)
        ax.scatter([rb], [yi], s=7 if fdr < 0.05 else 4, c=c, linewidths=0, zorder=2)
        if cat.get(g) == "primary_FDR_driver":
            ax.text(rb, yi, " " + g, fontsize=5.2, va="center", color="#111111")
    ax.axvline(0, color="#333333", lw=0.6)
    ax.set_yticks([])
    ax.set_xlim(-0.62, 0.52)
    ax.set_xticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4])
    ax.set_xlabel("gene-level rank-biserial (CRC \u2212 ESCA)")
    ax.text(0.02, len(ox) + 1.5,
            "246 evaluable OXPHOS genes \u00b7 module colours per visual dictionary \u00b7 9 primary drivers labelled",
            fontsize=6.0, color="#444444")
    ax.set_title("e", loc="left", fontsize=10, fontweight="bold")

def assemble():
    W, H = 180 * MM, 155 * MM
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.55, 1.0, 0.5], width_ratios=[1.15, 1.0, 1.0],
                          left=0.06, right=0.995, top=0.97, bottom=0.07, hspace=0.6, wspace=0.5)
    ax_a = fig.add_subplot(gs[0, :]); panel_2a(ax_a)
    ax_b = fig.add_subplot(gs[1, 0]); panel_2b(ax_b)
    ax_c = fig.add_subplot(gs[1, 1]); panel_2c(ax_c)
    ax_d = fig.add_subplot(gs[1, 2]); panel_2d(ax_d)
    ax_e = fig.add_subplot(gs[2, :]); panel_2e(ax_e)
    save(fig, os.path.join(ASSEMBLED, "Figure2"))

def panels_only():
    specs = [("2a", panel_2a, (170 * MM, 60 * MM)),
             ("2b", panel_2b, (55 * MM, 60 * MM)),
             ("2c", panel_2c, (55 * MM, 60 * MM)),
             ("2d", panel_2d, (55 * MM, 60 * MM)),
             ("2e", panel_2e, (170 * MM, 30 * MM))]
    for name, fn, size in specs:
        fig = plt.figure(figsize=size)
        ax = fig.add_axes([0.04, 0.05, 0.93, 0.92])
        fn(ax)
        save(fig, os.path.join(PANELS, f"Fig2_{name}"))

if __name__ == "__main__":
    panels_only()
    assemble()
    print("FIG2 DONE")
