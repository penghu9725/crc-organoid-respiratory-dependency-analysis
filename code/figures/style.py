# -*- coding: utf-8 -*-
"""Frozen visual dictionary helpers for FINAL figure production."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---- frozen dictionary colours (Okabe-Ito based; see FIGURE_VISUAL_DICTIONARY.md)
CRC   = "#0072B2"
ESCA  = "#E69F00"
TUM   = "#D55E00"
NORM  = "#56B4E9"
MOD_CI_STR = "#0072B2"
MOD_CI_ASM = "#D55E00"
MOD_CIV_STR = "#009E73"
MOD_CIV_ASM = "#CC79A7"
MOD_FES = "#E69F00"
MOD_MTFAS = "#F0E442"
GREY = "#BBBBBB"
DARK = "#333333"
NEUTRALS = {"Ovarian": "#999999", "Pancreatic": "#8C8C8C", "Gastric": "#7A7A7A",
            "Oesophageal": ESCA, "Colorectal": CRC}
MOD_COL = {"CI_structural": MOD_CI_STR, "CI_assembly": MOD_CI_ASM,
           "CIV_structural": MOD_CIV_STR, "CIV_assembly": MOD_CIV_ASM,
           "FeS_biogenesis": MOD_FES, "mtFAS": MOD_MTFAS}
MOD_LABEL = {"CI_structural": "Complex I structural", "CI_assembly": "Complex I assembly",
             "CIV_structural": "Complex IV structural", "CIV_assembly": "Complex IV assembly",
             "FeS_biogenesis": "Fe-S biogenesis", "mtFAS": "mtFAS (exploratory, 2 genes)"}
ROBUST_MODULES = ["CI_structural", "CI_assembly", "CIV_structural", "CIV_assembly", "FeS_biogenesis"]

MM = 1 / 25.4

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8.5, "axes.labelsize": 8, "axes.linewidth": 0.6,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",  # text as text in SVG
    "pdf.fonttype": 42,      # TrueType embedded
})

BASE = os.path.abspath(os.environ.get("ORGANOID_PROJECT_ROOT", os.path.join(os.path.dirname(__file__), "../..")))
FIGURE_ROOT = os.path.abspath(os.environ.get("ORGANOID_FIGURE_ROOT", os.path.join(BASE, "figures", "generated")))
PANELS = os.path.join(FIGURE_ROOT, "panels")
ASSEMBLED = os.path.join(FIGURE_ROOT, "assembled")
PREVIEW = os.path.join(FIGURE_ROOT, "preview")
for d in (PANELS, ASSEMBLED, PREVIEW):
    os.makedirs(d, exist_ok=True)


def save(fig, path_stem, dpi=300):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(path_stem + "." + ext, dpi=dpi, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print("saved", os.path.basename(path_stem))


def grayscale(src_png, out_png):
    from PIL import Image
    im = Image.open(src_png).convert("L")
    im.save(out_png)
    print("grayscale", os.path.basename(out_png))


def read_tsv(path, skip=0):
    import csv
    with open(path, encoding="utf-8-sig", newline="") as f:
        rd = list(csv.DictReader(f, delimiter="\t", skipinitialspace=True))
    return rd


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")
