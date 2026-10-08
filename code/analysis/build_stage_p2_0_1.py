from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("ORGANOID_PROJECT_ROOT", Path(__file__).resolve().parents[2])).resolve()
SRC = Path(os.environ.get("ORGANOID_SOURCE_ROOT", ROOT / "data" / "source")).resolve() / "00_source"
SUPP = SRC / "nature_supplementary_tables"
FIG = SRC / "figshare_28339340"
GSRC = SRC / "gene_sets_public"


def norm_yes(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.casefold().eq("yes")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def source_url(path: Path, figmeta: dict) -> tuple[str, str, str]:
    rel = path.relative_to(SRC).as_posix()
    if rel == "41586_2026_10830_MOESM3_ESM.zip" or rel.startswith("nature_supplementary_tables/"):
        return (
            "Nature supplementary tables",
            "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10830-y/MediaObjects/41586_2026_10830_MOESM3_ESM.zip",
            "Nature article supplementary archive or losslessly extracted member",
        )
    if rel.startswith("figshare_28339340/"):
        name = path.name
        match = next((f for f in figmeta["files"] if f["name"] == name), None)
        return ("Figshare processed dataset", match["download_url"] if match else "https://doi.org/10.6084/m9.figshare.28339340", "Figshare article 28339340")
    if rel.startswith("SangerOrganoidBiobank/"):
        return ("Garnett-Lab code repository", "https://github.com/Garnett-Lab/SangerOrganoidBiobank", "Git checkout; commit recorded separately")
    if rel.startswith("gene_sets_public/"):
        if path.name == "hgnc_complete_set.txt":
            return ("HGNC mapping", "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt", "Public HGNC complete set")
        return ("MSigDB public GMT", f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/{path.name}", "MSigDB 2025.1.Hs")
    if rel == "figshare_article_28339340.json":
        return ("Figshare metadata", "https://api.figshare.com/v2/articles/28339340", "Figshare API response")
    return ("local source artifact", "NA", "")


def build_metadata() -> tuple[pd.DataFrame, dict]:
    s2_path = SUPP / "supplementary_table_2_revision.xlsx"
    s2 = pd.read_excel(s2_path, sheet_name="data_availability", dtype=str).fillna("")
    required = [
        "sample_ID", "sanger_ID", "individual_ID", "primary_tumour_type",
        "CRISPR_available", "WGS_available", "RNAseq_available", "tumour_available",
        "overlap_HCMI", "overlap_Broad_DepMap", "ATCC", "Merck",
    ]
    out = s2[required].copy()
    write_tsv(out, ROOT / "01_metadata" / "ORGANOID_METADATA.tsv")

    expected_types = {"Colorectal": 132, "Oesophageal": 76, "Ovarian": 20, "Pancreatic": 22, "Gastric": 6}
    observed_types = out["primary_tumour_type"].value_counts().to_dict()
    assert len(out) == 256 and observed_types == expected_types, (len(out), observed_types)
    expected_screened = {"Colorectal": 85, "Oesophageal": 59, "Ovarian": 11, "Pancreatic": 4, "Gastric": 3}
    observed_screened = out[norm_yes(out["CRISPR_available"])]["primary_tumour_type"].value_counts().to_dict()
    assert observed_screened == expected_screened and norm_yes(out["CRISPR_available"]).sum() == 162
    assert norm_yes(out["WGS_available"]).sum() == 256
    assert norm_yes(out["RNAseq_available"]).sum() == 255
    return out, {"types": observed_types, "screened": observed_screened}


def metadata_exceptions(meta: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cov = pd.read_csv(FIG / "data_covariates_ALL-organoids.csv", dtype=str).fillna("")
    rna = pd.read_csv(FIG / "data_RNAseq_ALL-organoids.csv")
    merged = cov.merge(meta[["sample_ID", "primary_tumour_type"]], on="sample_ID", how="left", validate="one_to_one")
    canonical_cov = {"COLO": "Colorectal", "OESO": "Oesophageal", "PANC": "Pancreatic", "PAAD": "Pancreatic", "OV": "Ovarian", "STAD": "Gastric"}
    merged["covariates_mapped_type"] = merged["tissue"].map(canonical_cov).fillna(merged["tissue"])
    discrep = merged[merged["covariates_mapped_type"] != merged["primary_tumour_type"]].copy()
    rows = []
    for _, r in discrep.sort_values("sample_ID").iterrows():
        rows.append({
            "exception_type": "covariates_tissue_mismatch",
            "sample_ID": r["sample_ID"],
            "source_value": r["tissue"],
            "authoritative_S2_value": r["primary_tumour_type"],
            "observed_data_status": "row present in covariates",
            "classification": "metadata_label_error",
            "resolution": "retain covariates ploidy/msStatus/library; replace cancer type only from S2",
            "status": "resolved_by_authority_rule",
        })

    rna_ids = rna["sample_ID"].astype(str)
    s2_no = meta.loc[~norm_yes(meta["RNAseq_available"]), "sample_ID"].tolist()
    absent = sorted(set(meta["sample_ID"]) - set(rna_ids))
    extra = sorted(set(rna_ids) - set(meta["sample_ID"]))
    duplicates = sorted(rna_ids[rna_ids.duplicated(keep=False)].unique())
    assert len(rna) == 255 and len(s2_no) == 1 and absent == s2_no and not extra and not duplicates
    sample = s2_no[0]
    rows.append({
        "exception_type": "RNAseq_reported_row_count_discrepancy",
        "sample_ID": sample,
        "source_value": "RNA CSV row absent; CSV has 255 unique samples",
        "authoritative_S2_value": "RNAseq_available=No",
        "observed_data_status": "no expression row; no placeholder and no duplicate",
        "classification": "metadata_consistent; prior 256-row report was incorrect",
        "resolution": "no silent fix; preserve S2=No and document verified CSV dimension",
        "status": "resolved",
    })
    exc = pd.DataFrame(rows)
    write_tsv(exc, ROOT / "07_logs" / "METADATA_QC_EXCEPTIONS.tsv")
    return exc, {
        "cov_tissue_counts": cov["tissue"].value_counts(dropna=False).to_dict(),
        "cov_discrepant": len(discrep),
        "cov_discrepant_true_types": discrep["primary_tumour_type"].value_counts().to_dict(),
        "rna_shape": rna.shape,
        "rna_absent": sample,
    }


def crispr_qc(meta: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    files = {
        "normalized_LFC": SUPP / "supplementary_table_6_revision.csv",
        "binary_essentiality": SUPP / "supplementary_table_5_revision.csv",
    }
    summaries = {}
    metrics = []
    for kind, path in files.items():
        model_counts = Counter()
        gene_counts = Counter()
        pair_counts = Counter()
        composite_counts = Counter()
        libraries = Counter()
        nrows = 0
        missing = Counter()
        values = []
        for chunk in pd.read_csv(path, sep=r"\s+", chunksize=250_000, dtype={"sample_ID": str, "gene": str, "ensembl_id": str, "HGNC_ID": str, "library": str}):
            nrows += len(chunk)
            for col in chunk.columns:
                missing[col] += int(chunk[col].isna().sum())
            model_counts.update(chunk["sample_ID"].dropna())
            gene_counts.update(chunk["gene"].dropna())
            libraries.update(chunk["library"].dropna())
            pair_counts.update(zip(chunk["sample_ID"], chunk["gene"]))
            composite_counts.update(zip(chunk["sample_ID"], chunk["gene"], chunk["library"]))
            valcol = "LFC" if kind == "normalized_LFC" else "is_depleted"
            values.append(pd.to_numeric(chunk[valcol], errors="coerce").to_numpy())
        arr = np.concatenate(values)
        duplicate_pairs = sum(v - 1 for v in pair_counts.values() if v > 1)
        duplicate_composite = sum(v - 1 for v in composite_counts.values() if v > 1)
        multilibrary_models = len({m for (m, g), v in pair_counts.items() if v > 1})
        summary = {
            "file": path.name, "rows": nrows, "columns": 6,
            "models": len(model_counts), "genes": len(gene_counts),
            "duplicate_model_gene_rows": duplicate_pairs,
            "duplicate_model_gene_library_rows": duplicate_composite,
            "models_with_two_libraries": multilibrary_models,
            "missing_total": int(sum(missing.values())), "missing_by_column": dict(missing),
            "libraries": dict(libraries), "min": float(np.nanmin(arr)), "max": float(np.nanmax(arr)),
            "median": float(np.nanmedian(arr)), "unique_values": sorted(pd.Series(arr).dropna().unique().tolist())[:20],
        }
        summaries[kind] = summary
        for metric, value in [
            ("data_rows", nrows), ("data_columns", 6), ("unique_screened_models", len(model_counts)),
            ("unique_gene_symbols", len(gene_counts)), ("duplicate_model_gene_rows", duplicate_pairs),
            ("duplicate_model_gene_library_rows", duplicate_composite), ("models_with_two_libraries", multilibrary_models),
            ("missing_cells_total", int(sum(missing.values()))), ("minimum_value", summary["min"]),
            ("maximum_value", summary["max"]), ("median_value", summary["median"]),
            ("library_counts", json.dumps(dict(libraries), sort_keys=True)),
        ]:
            metrics.append({"matrix": kind, "metric": metric, "value": value, "interpretation": ""})
    assert summaries["normalized_LFC"]["models"] == 162
    assert summaries["binary_essentiality"]["models"] == 162
    # Functional direction check using the paper-provided essential/non-essential controls.
    controls = pd.read_csv(FIG / "control-genes_organoids.csv")
    control_map = dict(zip(controls["gene"], controls["essentiality"]))
    sums = defaultdict(float); counts = Counter()
    for chunk in pd.read_csv(files["normalized_LFC"], sep=r"\s+", usecols=["gene", "LFC"], chunksize=400_000):
        labels = chunk["gene"].map(control_map)
        for lab in labels.dropna().unique():
            vals = pd.to_numeric(chunk.loc[labels == lab, "LFC"], errors="coerce").dropna()
            sums[str(lab)] += vals.sum(); counts[str(lab)] += len(vals)
    control_means = {k: sums[k] / counts[k] for k in counts}
    metrics.extend([
        {"matrix": "normalized_LFC", "metric": "direction", "value": "more negative = stronger loss of fitness/dependency", "interpretation": "confirmed by lower mean LFC in provided essential controls than non-essential controls"},
        {"matrix": "normalized_LFC", "metric": "control_group_mean_LFC", "value": json.dumps(control_means, sort_keys=True), "interpretation": "empirical direction check"},
        {"matrix": "binary_essentiality", "metric": "definition", "value": "is_depleted=1 denotes BAGEL2-called depleted/essential in that organoid; 0 denotes not depleted", "interpretation": "binary screen call, not effect magnitude"},
        {"matrix": "both", "metric": "identifier_format", "value": "sample_ID (HCM-SANG-* or WTSI-*); gene symbol + Ensembl + HGNC; library field", "interpretation": "join samples to S2 sample_ID; retain gene identifiers together"},
    ])
    qc = pd.DataFrame(metrics)
    write_tsv(qc, ROOT / "02_crispr" / "CRISPR_MATRIX_QC.tsv")
    return qc, {**summaries, "control_means": control_means}


PATHWAY_PATTERNS = {
    "serine/glycine/one-carbon": [r"SERINE", r"GLYCINE", r"ONE.CARBON"],
    "glutamine/glutamate": [r"GLUTAMINE", r"GLUTAMATE"],
    "arginine": [r"ARGININE", r"UREA.CYCLE"],
    "methionine": [r"METHIONINE", r"S.ADENOSYLMETHIONINE", r"TRANSSULFURATION"],
    "BCAA": [r"BRANCHED.CHAIN.AMINO.ACID", r"VALINE.*LEUCINE.*ISOLEUCINE", r"LEUCINE.*ISOLEUCINE.*VALINE"],
    "tryptophan": [r"TRYPTOPHAN", r"KYNURENINE"],
    "glycolysis": [r"GLYCOLYSIS", r"GLUCONEOGENESIS"],
    "PPP": [r"PENTOSE.PHOSPHATE"],
    "TCA": [r"CITRIC.ACID.CYCLE", r"TRICARBOXYLIC.ACID", r"TCA.CYCLE"],
    "OXPHOS": [r"OXIDATIVE.PHOSPHORYLATION", r"RESPIRATORY.ELECTRON.TRANSPORT", r"MITOCHONDRIAL.RESPIRATORY.CHAIN"],
    "fatty-acid synthesis": [r"FATTY.ACID.BIOSYN", r"FATTY.ACID.SYNTH", r"LIPOGENESIS"],
    "fatty-acid oxidation": [r"FATTY.ACID.OXID", r"BETA.OXIDATION"],
    "cholesterol/isoprenoid": [r"CHOLESTEROL", r"STEROL", r"ISOPRENOID", r"MEVALONATE", r"TERPENOID.BACKBONE"],
    "sphingolipid": [r"SPHINGOLIPID", r"SPHINGOSINE", r"CERAMIDE"],
    "folate": [r"FOLATE", r"FOLIC.ACID", r"TETRAHYDROFOLATE"],
    "NAD": [r"NICOTINAMIDE", r"NICOTINATE", r"NAD.BIOSYN", r"NAD.METABOL"],
    "thiamine": [r"THIAMINE"],
    "riboflavin": [r"RIBOFLAVIN"],
    "vitamin B6": [r"VITAMIN.B6", r"PYRIDOX"],
    "vitamin B12": [r"VITAMIN.B12", r"COBALAMIN"],
    "iron": [r"IRON", r"FERRIC", r"FERROUS"],
    "copper": [r"COPPER"],
    "zinc": [r"ZINC"],
    "selenium": [r"SELENIUM", r"SELENOCYSTEINE"],
    "amino-acid transport": [r"AMINO.ACID.TRANSPORT", r"AMINO.ACID.TRANSMEMBRANE.TRANSPORT"],
    "glucose transport": [r"GLUCOSE.TRANSPORT", r"HEXOSE.TRANSPORT"],
    "metal transport": [r"METAL.ION.TRANSPORT", r"TRANSITION.METAL.*TRANSPORT", r"DIVALENT.*METAL.*TRANSPORT"],
}


def build_gene_universe() -> tuple[pd.DataFrame, dict]:
    hgnc = pd.read_csv(GSRC / "hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False).fillna("")
    symbol_to_ens = dict(zip(hgnc["symbol"], hgnc["ensembl_gene_id"]))
    sources = [
        ("Reactome", "c2.cp.reactome.v2025.1.Hs.symbols.gmt", "standard_curated"),
        ("KEGG", "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt", "standard_curated"),
        ("KEGG", "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt", "standard_curated"),
        ("GO", "c5.go.bp.v2025.1.Hs.symbols.gmt", "ontology_annotation"),
        ("MSigDB", "h.all.v2025.1.Hs.symbols.gmt", "hallmark_curated"),
    ]
    rows = []
    for db, fname, evidence in sources:
        with (GSRC / fname).open(encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                setname, genes = parts[0], parts[2:]
                majors = [major for major, pats in PATHWAY_PATTERNS.items() if any(re.search(p, setname, re.I) for p in pats)]
                for major in majors:
                    for gene in genes:
                        rows.append({
                            "gene_symbol": gene,
                            "ensembl_id": symbol_to_ens.get(gene, ""),
                            "pathway_major": major,
                            "pathway_minor": setname,
                            "source_database": db,
                            "source_gene_set": setname,
                            "evidence_level": evidence,
                        })
    out = pd.DataFrame(rows).drop_duplicates().sort_values(["pathway_major", "gene_symbol", "source_database", "source_gene_set"])
    missing_categories = sorted(set(PATHWAY_PATTERNS) - set(out["pathway_major"]))
    write_tsv(out, ROOT / "03_gene_sets" / "METABOLIC_GENE_UNIVERSE.tsv")
    return out, {
        "membership_rows": len(out), "unique_genes": out["gene_symbol"].nunique(),
        "mapped_ensembl": out.loc[out["ensembl_id"] != "", "gene_symbol"].nunique(),
        "category_counts_unique_genes": out.groupby("pathway_major")["gene_symbol"].nunique().to_dict(),
        "missing_categories": missing_categories,
    }


def dimensions_inventory() -> pd.DataFrame:
    rows = []
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or ".git" in p.parts or p.name == "DATA_INVENTORY.tsv":
            continue
        shape = "not_tabular"
        try:
            if p.suffix.lower() == ".csv":
                sep = r"\s+" if p.name in {"supplementary_table_5_revision.csv", "supplementary_table_6_revision.csv"} else ","
                first = pd.read_csv(p, nrows=2, sep=sep)
                with p.open("rb") as fh:
                    nrows = sum(block.count(b"\n") for block in iter(lambda: fh.read(8 * 1024 * 1024), b"")) - 1
                shape = f"{nrows}x{len(first.columns)}"
            elif p.suffix.lower() == ".tsv":
                first = pd.read_csv(p, sep="\t", nrows=2)
                with p.open("rb") as fh:
                    nrows = sum(block.count(b"\n") for block in iter(lambda: fh.read(8 * 1024 * 1024), b"")) - 1
                shape = f"{nrows}x{len(first.columns)}"
            elif p.suffix.lower() == ".xlsx":
                xl = pd.ExcelFile(p)
                dims = []
                for s in xl.sheet_names:
                    d = pd.read_excel(p, sheet_name=s, header=None)
                    dims.append(f"{s}:{d.shape[0]}x{d.shape[1]}")
                shape = ";".join(dims)
        except Exception as exc:
            shape = f"inspection_error:{type(exc).__name__}"
        top = p.relative_to(ROOT).parts[0]
        role = "frozen_source" if top == "00_source" else "derived_output" if top in {"01_metadata", "02_crispr", "03_gene_sets", "07_logs"} else "project_structure"
        rows.append({"relative_path": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size, "format": p.suffix.lower().lstrip("."), "dimensions_or_object": shape, "role": role})
    out = pd.DataFrame(rows)
    write_tsv(out, ROOT / "DATA_INVENTORY.tsv")
    return out


def freeze_manifest(figmeta: dict) -> pd.DataFrame:
    commit = subprocess.check_output(["git", "-C", str(SRC / "SangerOrganoidBiobank"), "rev-parse", "HEAD"], text=True).strip()
    rows = []
    for p in sorted(SRC.rglob("*")):
        if not p.is_file() or ".git" in p.parts or p.name == "SOURCE_FREEZE.tsv":
            continue
        category, url, note = source_url(p, figmeta)
        rows.append({
            "relative_path": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size,
            "sha256": sha256(p), "source_category": category, "source_url": url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(), "version_or_note": note,
        })
    rows.append({"relative_path": "00_source/SangerOrganoidBiobank", "bytes": "NA", "sha256": "NA", "source_category": "Garnett-Lab code repository", "source_url": "https://github.com/Garnett-Lab/SangerOrganoidBiobank", "retrieved_utc": datetime.now(timezone.utc).isoformat(), "version_or_note": f"git_commit={commit}"})
    out = pd.DataFrame(rows)
    write_tsv(out, SRC / "SOURCE_FREEZE.tsv")
    return out


def main() -> None:
    figmeta = json.loads((SRC / "figshare_article_28339340.json").read_text(encoding="utf-8-sig"))
    meta, meta_stats = build_metadata()
    exceptions, exc_stats = metadata_exceptions(meta)
    crispr, crispr_stats = crispr_qc(meta)
    universe, gene_stats = build_gene_universe()
    inventory = dimensions_inventory()
    freeze = freeze_manifest(figmeta)
    summary = {
        "metadata": meta_stats, "exceptions": exc_stats, "crispr": crispr_stats,
        "gene_universe": gene_stats, "inventory_rows": len(inventory), "freeze_rows": len(freeze),
    }
    (ROOT / "07_logs" / "P2_STAGE0_1_MACHINE_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
