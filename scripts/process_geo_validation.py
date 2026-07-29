#!/usr/bin/env python
"""Process GSE13507/GSE31684 and apply fixed TCGA LASSO-Cox formula."""

from __future__ import annotations

import csv
import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_libs"))

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "external" / "raw"
PROCESSED = ROOT / "data" / "external" / "processed"
TABLES = ROOT / "results" / "tables"
MODEL = ROOT / "results" / "model" / "final_lasso_cox_model.tsv"


def read_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_meta: dict[str, list[str]] = {}
    rows = []
    header = None
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                header = [x.strip('"') for x in next(handle).rstrip("\n").split("\t")]
                continue
            if line == "!series_matrix_table_end":
                break
            if in_table:
                rows.append([x.strip('"') for x in line.split("\t")])
            elif line.startswith("!Sample_"):
                parts = [x.strip('"') for x in line.split("\t")]
                sample_meta[parts[0]] = parts[1:]
    if header is None:
        raise RuntimeError(f"No matrix table found in {path}")
    expr = pd.DataFrame(rows, columns=header)
    expr = expr.rename(columns={"ID_REF": "probe_id"}).set_index("probe_id")
    expr = expr.apply(pd.to_numeric, errors="coerce")

    meta = pd.DataFrame(
        {
            "geo_accession": sample_meta.get("!Sample_geo_accession", []),
            "title": sample_meta.get("!Sample_title", []),
            "source_name": sample_meta.get("!Sample_source_name_ch1", []),
            "platform": sample_meta.get("!Sample_platform_id", []),
        }
    )
    return expr, meta


def read_platform_mapping(platform: str) -> pd.DataFrame:
    path = RAW / "platforms" / f"{platform}.annot.gz"
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = None
        for line in handle:
            if line.startswith("ID\t"):
                reader = csv.DictReader(handle, fieldnames=line.rstrip("\n").split("\t"), delimiter="\t")
                break
        if reader is None:
            raise RuntimeError(f"No annotation table found in {path}")
        rows = []
        for rec in reader:
            probe = rec.get("ID")
            symbols = rec.get("Gene symbol") or rec.get("Gene Symbol") or rec.get("Symbol") or ""
            if not probe or not symbols:
                continue
            for symbol in re.split(r"///|;|,|\s//\s", symbols):
                symbol = symbol.strip()
                if symbol and symbol != "---":
                    rows.append({"probe_id": probe, "gene_symbol": symbol})
    mapping = pd.DataFrame(rows).drop_duplicates()
    return mapping


def collapse_to_gene(expr: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    common = mapping[mapping["probe_id"].isin(expr.index)].copy()
    merged = common.merge(expr, left_on="probe_id", right_index=True, how="inner")
    value_cols = [c for c in merged.columns if c not in {"probe_id", "gene_symbol"}]
    gene_expr = merged.groupby("gene_symbol")[value_cols].mean(numeric_only=True)
    return gene_expr


def gse13507_clinical(meta: pd.DataFrame) -> pd.DataFrame:
    clinical = pd.read_excel(RAW / "GSE13507" / "GSE13507_clinical_info.xls", sheet_name="Clinical.Info")
    primary = meta[meta["source_name"].eq("Primary bladder cancer")].copy()
    primary["sample_name"] = primary["title"].str.extract(r"(BT\d+)")
    out = primary.merge(clinical, left_on="sample_name", right_on="Sample name", how="inner")
    out["cohort"] = "GSE13507"
    out["sample_id"] = out["geo_accession"]
    out["patient_id"] = out["sample_name"]
    out["os_time_days"] = pd.to_numeric(out["survivalMonth"], errors="coerce") * 30.4375
    out["os_event"] = pd.to_numeric(out["overall survival"], errors="coerce").map(lambda x: 1 if x == 2 else (0 if x == 1 else np.nan))
    out["age"] = pd.to_numeric(out["AGE"], errors="coerce")
    out["gender"] = out["SEX"].map({"M": "male", "F": "female"})
    out["stage"] = out["TMN stage"]
    out["grade"] = out["Grade"]
    keep = ["cohort", "sample_id", "patient_id", "title", "source_name", "os_time_days", "os_event", "age", "gender", "stage", "grade"]
    return out[keep]


def gse31684_clinical(meta: pd.DataFrame) -> pd.DataFrame:
    clinical = pd.read_csv(RAW / "GSE31684" / "GSE31684_table_of_clinical_details.txt.gz", sep="\t")
    out = meta.merge(clinical, left_on="geo_accession", right_on="GEO", how="inner")
    out["cohort"] = "GSE31684"
    out["sample_id"] = out["geo_accession"]
    out["patient_id"] = out["ID"]
    out["os_time_days"] = pd.to_numeric(out["Survival Months"], errors="coerce") * 30.4375
    out["os_event"] = out["Last known status"].map(lambda x: 0 if x == "NED" else (1 if x in {"DOD", "DOC"} else np.nan))
    out["dss_event"] = out["Last known status"].map(lambda x: 1 if x == "DOD" else (0 if x in {"NED", "DOC"} else np.nan))
    out["age"] = pd.to_numeric(out["Age at RC"], errors="coerce")
    out["gender"] = out["Gender"]
    out["stage"] = out["RC_Stage"]
    out["grade"] = out["RC Grade"]
    keep = [
        "cohort",
        "sample_id",
        "patient_id",
        "title",
        "source_name",
        "os_time_days",
        "os_event",
        "dss_event",
        "age",
        "gender",
        "stage",
        "grade",
        "Last known status",
    ]
    return out[keep]


def apply_model(gene_expr: pd.DataFrame, clinical: pd.DataFrame, cohort: str) -> pd.DataFrame:
    model = pd.read_csv(MODEL, sep="\t")
    sample_cols = clinical["sample_id"].tolist()
    gene_expr = gene_expr.loc[:, [c for c in sample_cols if c in gene_expr.columns]]
    available = [g for g in model["gene"] if g in gene_expr.index]
    missing = [g for g in model["gene"] if g not in gene_expr.index]
    score = pd.Series(0.0, index=gene_expr.columns)
    used_rows = []
    for _, row in model.iterrows():
        gene = row["gene"]
        if gene not in gene_expr.index:
            continue
        z = (gene_expr.loc[gene, gene_expr.columns] - row["train_mean_log2_tpm"]) / row["train_sd_log2_tpm"]
        score += row["coefficient"] * z
        used_rows.append(row.to_dict())
    out = clinical[clinical["sample_id"].isin(score.index)].copy()
    out["tcga_fixed_risk_score"] = out["sample_id"].map(score.to_dict())
    cutoff = pd.read_csv(TABLES / "tcga_lasso_risk_scores.tsv", sep="\t")
    train_cutoff = cutoff.loc[cutoff["split"].eq("train"), "lasso_risk_score"].median()
    out["tcga_fixed_risk_group"] = np.where(out["tcga_fixed_risk_score"] >= train_cutoff, "High", "Low")
    out["n_model_genes_available"] = len(available)
    out["n_model_genes_missing"] = len(missing)
    out["missing_model_genes"] = ";".join(missing)
    pd.DataFrame(used_rows).to_csv(PROCESSED / f"{cohort}_tcga_model_genes_used.tsv", sep="\t", index=False)
    pd.DataFrame({"cohort": [cohort], "available_genes": [len(available)], "missing_genes": [len(missing)], "missing": [";".join(missing)]}).to_csv(
        TABLES / f"{cohort}_model_gene_availability.tsv", sep="\t", index=False
    )
    return out


def process_one(cohort: str, platform: str) -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    expr, meta = read_series_matrix(RAW / cohort / f"{cohort}_series_matrix.txt.gz")
    mapping = read_platform_mapping(platform)
    gene_expr = collapse_to_gene(expr, mapping)
    if cohort == "GSE13507":
        clinical = gse13507_clinical(meta)
    elif cohort == "GSE31684":
        clinical = gse31684_clinical(meta)
    else:
        raise ValueError(cohort)
    gene_expr = gene_expr.loc[:, [c for c in clinical["sample_id"] if c in gene_expr.columns]]
    gene_expr.to_csv(PROCESSED / f"{cohort}_gene_symbol_expression.tsv.gz", sep="\t", compression="gzip")
    clinical.to_csv(PROCESSED / f"{cohort}_clinical_survival.tsv", sep="\t", index=False)
    scored = apply_model(gene_expr, clinical, cohort)
    scored.to_csv(PROCESSED / f"{cohort}_tcga_fixed_risk_scores.tsv", sep="\t", index=False)
    print(
        f"{cohort}: genes={gene_expr.shape[0]} samples={gene_expr.shape[1]} "
        f"clinical={clinical.shape[0]} events={clinical['os_event'].sum()} "
        f"model_genes_available={scored['n_model_genes_available'].iloc[0]} missing={scored['n_model_genes_missing'].iloc[0]}"
    )


def main() -> int:
    process_one("GSE13507", "GPL6102")
    process_one("GSE31684", "GPL570")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
