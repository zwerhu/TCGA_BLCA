#!/usr/bin/env python
"""
Prepare downstream-ready TCGA-BLCA inputs from processed matrices.

Outputs are simple TSV/TSV.GZ files that R packages such as survival, DESeq2,
oncoPredict, pRRophetic, GSVA, and maftools can consume.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)

    clinical = pd.read_csv(PROCESSED / "blca_clinical_survival.tsv", sep="\t")
    sample_map = pd.read_csv(PROCESSED / "blca_rna_sample_map.tsv", sep="\t")
    annotation = pd.read_csv(PROCESSED / "blca_gene_annotation.tsv", sep="\t")
    tpm = pd.read_csv(PROCESSED / "blca_rna_tpm_matrix.tsv.gz", sep="\t", index_col=0)

    sample_map["patient_barcode"] = sample_map["sample_barcode"].str.slice(0, 12)
    sample_counts = (
        sample_map.groupby(["sample_type", "tissue_type"], dropna=False)
        .size()
        .reset_index(name="n_samples")
        .sort_values("n_samples", ascending=False)
    )
    sample_counts.to_csv(TABLES / "rna_sample_type_counts.tsv", sep="\t", index=False)

    os_summary = {
        "n_clinical_cases": int(clinical["patient_barcode"].nunique()),
        "n_with_os_time": int(clinical["os_time_days"].notna().sum()),
        "n_os_events": int(clinical["os_event"].fillna(0).astype(int).sum()),
        "median_os_time_days": float(clinical["os_time_days"].dropna().median()),
    }
    pd.DataFrame([os_summary]).to_csv(TABLES / "clinical_os_summary.tsv", sep="\t", index=False)

    tumor_cols = sample_map.loc[sample_map["sample_type"].eq("Primary Tumor"), "matrix_column"].tolist()
    normal_cols = sample_map.loc[sample_map["sample_type"].eq("Solid Tissue Normal"), "matrix_column"].tolist()
    tumor_cols = [c for c in tumor_cols if c in tpm.columns]
    normal_cols = [c for c in normal_cols if c in tpm.columns]

    expr_log2_tpm = np.log2(tpm + 1.0)
    gene_anno = annotation.set_index("gene_id").reindex(expr_log2_tpm.index)
    fallback_gene_id = pd.Series(expr_log2_tpm.index, index=expr_log2_tpm.index)
    gene_info = pd.DataFrame(
        {
            "gene_symbol": gene_anno["gene_name"].fillna(fallback_gene_id),
            "gene_type": gene_anno["gene_type"].fillna(""),
        },
        index=expr_log2_tpm.index,
    )
    log2_tpm = pd.concat([gene_info, expr_log2_tpm], axis=1)

    all_symbol = (
        log2_tpm.drop(columns=["gene_type"])
        .groupby("gene_symbol", sort=False)
        .mean(numeric_only=True)
    )
    all_symbol.to_csv(PROCESSED / "blca_log2_tpm_gene_symbol_all_samples.tsv.gz", sep="\t", compression="gzip")
    all_symbol[tumor_cols].to_csv(PROCESSED / "blca_log2_tpm_gene_symbol_tumor.tsv.gz", sep="\t", compression="gzip")
    all_symbol[tumor_cols].to_csv(PROCESSED / "blca_drug_prediction_expression_input.tsv.gz", sep="\t", compression="gzip")

    if tumor_cols and normal_cols:
        fc = pd.DataFrame(
            {
                "gene_id": log2_tpm.index,
                "gene_symbol": log2_tpm["gene_symbol"].values,
                "gene_type": log2_tpm["gene_type"].values,
                "tumor_mean_log2tpm": log2_tpm[tumor_cols].mean(axis=1).values,
                "normal_mean_log2tpm": log2_tpm[normal_cols].mean(axis=1).values,
            }
        )
        fc["exploratory_log2fc_tumor_vs_normal"] = fc["tumor_mean_log2tpm"] - fc["normal_mean_log2tpm"]
        fc = fc.sort_values("exploratory_log2fc_tumor_vs_normal", ascending=False)
        fc.to_csv(TABLES / "exploratory_log2fc_tumor_vs_normal.tsv", sep="\t", index=False)

    survival_sheet = sample_map[sample_map["matrix_column"].isin(tumor_cols)].merge(
        clinical,
        how="left",
        on="patient_barcode",
        suffixes=("", "_clinical"),
    )
    keep_cols = [
        "matrix_column",
        "sample_barcode",
        "patient_barcode",
        "sample_type",
        "os_time_days",
        "os_event",
        "age_at_diagnosis_years",
        "gender",
        "ajcc_pathologic_stage",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "ajcc_pathologic_m",
        "tumor_grade",
        "tobacco_smoking_status",
        "pack_years_smoked",
    ]
    survival_sheet[keep_cols].to_csv(PROCESSED / "blca_survival_expression_sample_sheet.tsv", sep="\t", index=False)

    print(f"tumor RNA samples: {len(tumor_cols)}")
    print(f"normal RNA samples: {len(normal_cols)}")
    print(f"gene-symbol expression matrix: {all_symbol.shape[0]} genes x {all_symbol.shape[1]} samples")
    print(f"survival sample sheet: {survival_sheet.shape[0]} tumor samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
