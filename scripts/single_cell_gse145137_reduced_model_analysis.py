#!/usr/bin/env python
"""Project the GEO-compatible reduced TCGA-BLCA model onto GSE145137 scRNA cells."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SC_DIR = ROOT / "data" / "external" / "single_cell" / "GSE145137"
OUT_DIR = ROOT / "results" / "single_cell"
TABLE_DIR = ROOT / "results" / "tables"
MODEL_FILE = ROOT / "results" / "model" / "reduced_geo_compatible_model.tsv"

PRIMARY_EXPR = SC_DIR / "GSM4307111_GEO_processed_BC159-T_3_log2TPM_matrix_final.txt.gz"
CELLTYPE_XLSX = SC_DIR / "GSE145137_BC159-T_3_All_SC_final_QC_Celltype_Information.xlsx"

MODEL_GENES = ["EMP1", "AHNAK", "TNFRSF14", "CLEC2D", "GSDMB"]


def read_expression_subset(path: Path, genes: list[str] | None = None) -> pd.DataFrame:
    """Read all genes or a requested gene subset from a gzipped GEO log2TPM matrix."""
    if genes is None:
        return pd.read_csv(path, sep="\t", compression="gzip", index_col=0)

    wanted = set(genes)
    rows: list[list[str]] = []
    header: list[str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in wanted:
                rows.append(parts)
    if not rows or header is None:
        return pd.DataFrame()
    df = pd.DataFrame(rows, columns=header).set_index("gene")
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def select_variable_expression(path: Path, cells: list[str], n_genes: int = 1500) -> pd.DataFrame:
    """Read the matrix and keep the most variable genes for lightweight PCA."""
    chunks = pd.read_csv(path, sep="\t", compression="gzip", index_col=0, chunksize=1000)
    variances: list[pd.Series] = []
    keep_chunks: list[pd.DataFrame] = []
    for chunk in chunks:
        chunk = chunk[cells].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        variances.append(chunk.var(axis=1))
        keep_chunks.append(chunk)
    var = pd.concat(variances).sort_values(ascending=False)
    top = set(var.head(n_genes).index)
    top_chunks = [chunk.loc[chunk.index.intersection(top)] for chunk in keep_chunks]
    expr = pd.concat(top_chunks, axis=0)
    return expr.loc[var.head(n_genes).index]


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    values = df.astype(float)
    means = values.mean(axis=1)
    sds = values.std(axis=1, ddof=0).replace(0, np.nan)
    return values.sub(means, axis=0).div(sds, axis=0).fillna(0.0)


def pca_from_expression(expr_gene_by_cell: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    x = expr_gene_by_cell.T.astype(float)
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0).replace(0, np.nan)
    x = x.fillna(0.0).to_numpy(dtype=np.float64).copy()
    x -= x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :2] * s[:2]
    var = s**2
    explained = var / var.sum()
    coords = pd.DataFrame(
        {
            "cell": expr_gene_by_cell.columns,
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
        }
    )
    return coords, {"PC1": float(explained[0]), "PC2": float(explained[1])}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    model = pd.read_csv(MODEL_FILE, sep="\t")
    model = model[model["gene"].isin(MODEL_GENES)].copy()
    coeff = dict(zip(model["gene"], model["coefficient"]))

    meta = pd.read_excel(CELLTYPE_XLSX, sheet_name="BC159-T#3")
    meta = meta.rename(columns={"Samples": "cell"})
    meta["cell"] = meta["cell"].astype(str)
    meta["cell_type"] = meta["cell_type"].astype(str)
    meta["cell_index"] = meta["cell_index"].astype(str)

    model_expr = read_expression_subset(PRIMARY_EXPR, MODEL_GENES)
    present = [g for g in MODEL_GENES if g in model_expr.index]
    missing = [g for g in MODEL_GENES if g not in model_expr.index]
    if missing:
        raise RuntimeError(f"Missing model genes in GSE145137 primary matrix: {missing}")

    common_cells = [cell for cell in meta["cell"] if cell in model_expr.columns]
    meta = meta[meta["cell"].isin(common_cells)].copy()
    model_expr = model_expr[common_cells]

    z = zscore_rows(model_expr.loc[present])
    score = pd.Series(0.0, index=common_cells)
    contribution_rows = []
    for gene in present:
        contribution = z.loc[gene] * coeff[gene]
        score += contribution
        for cell, value in contribution.items():
            contribution_rows.append({"cell": cell, "gene": gene, "contribution": value})

    meta["reduced_model_score"] = meta["cell"].map(score)
    meta["reduced_model_score_z"] = (
        (meta["reduced_model_score"] - meta["reduced_model_score"].mean())
        / meta["reduced_model_score"].std(ddof=0)
    )
    meta["risk_group_sc_median"] = np.where(
        meta["reduced_model_score"] >= meta["reduced_model_score"].median(), "High", "Low"
    )

    expr_long = (
        model_expr.loc[present]
        .T.reset_index(names="cell")
        .melt(id_vars="cell", var_name="gene", value_name="log2_tpm")
    )
    expr_long = expr_long.merge(meta[["cell", "cell_type", "cell_index", "reduced_model_score_z"]], on="cell")

    summary_rows = []
    for (cell_type, gene), sub in expr_long.groupby(["cell_type", "gene"], sort=False):
        vals = sub["log2_tpm"].astype(float)
        summary_rows.append(
            {
                "cell_type": cell_type,
                "gene": gene,
                "n_cells": int(vals.shape[0]),
                "mean_log2_tpm": float(vals.mean()),
                "median_log2_tpm": float(vals.median()),
                "pct_positive": float((vals > 0).mean() * 100.0),
            }
        )
    gene_summary = pd.DataFrame(summary_rows)
    gene_summary["mean_z_by_gene"] = gene_summary.groupby("gene")["mean_log2_tpm"].transform(
        lambda x: (x - x.mean()) / (x.std(ddof=0) if x.std(ddof=0) != 0 else 1.0)
    )

    risk_summary = (
        meta.groupby(["cell_type", "cell_index"], sort=False)
        .agg(
            n_cells=("cell", "size"),
            mean_score_z=("reduced_model_score_z", "mean"),
            median_score_z=("reduced_model_score_z", "median"),
            high_fraction=("risk_group_sc_median", lambda x: float((x == "High").mean())),
            mean_n_gene=("nGene", "mean"),
            mean_n_umi=("nUMI", "mean"),
        )
        .reset_index()
        .sort_values("median_score_z", ascending=False)
    )

    var_expr = select_variable_expression(PRIMARY_EXPR, common_cells, n_genes=1500)
    pca, pca_explained = pca_from_expression(var_expr)
    pca = pca.merge(meta, on="cell", how="left")

    meta_out = SC_DIR / "GSE145137_primary_cell_metadata_with_reduced_model_scores.tsv"
    expr_long_out = SC_DIR / "GSE145137_primary_reduced_model_gene_expression_long.tsv"
    gene_summary_out = TABLE_DIR / "single_cell_gse145137_model_gene_celltype_summary.tsv"
    risk_summary_out = TABLE_DIR / "single_cell_gse145137_risk_by_celltype_summary.tsv"
    pca_out = SC_DIR / "GSE145137_primary_pca_coordinates.tsv"
    contribution_out = SC_DIR / "GSE145137_primary_reduced_model_gene_contributions_long.tsv"
    manifest_out = OUT_DIR / "single_cell_gse145137_analysis_manifest.json"

    meta.to_csv(meta_out, sep="\t", index=False)
    expr_long.to_csv(expr_long_out, sep="\t", index=False)
    gene_summary.to_csv(gene_summary_out, sep="\t", index=False)
    risk_summary.to_csv(risk_summary_out, sep="\t", index=False)
    pca.to_csv(pca_out, sep="\t", index=False)
    pd.DataFrame(contribution_rows).to_csv(contribution_out, sep="\t", index=False)

    manifest = {
        "dataset": "GSE145137",
        "sample": "BC159-T#3 primary bladder cancer",
        "n_cells": int(meta.shape[0]),
        "n_genes_primary_matrix": 14233,
        "model_genes": present,
        "missing_model_genes": missing,
        "risk_score_transform": "cell-wise score = sum(TCGA coefficient * gene z-score across primary single cells)",
        "pca_top_variable_genes": int(var_expr.shape[0]),
        "pca_explained_variance": pca_explained,
        "outputs": {
            "cell_metadata_scores": str(meta_out),
            "model_gene_expression_long": str(expr_long_out),
            "gene_celltype_summary": str(gene_summary_out),
            "risk_celltype_summary": str(risk_summary_out),
            "pca_coordinates": str(pca_out),
            "gene_contributions": str(contribution_out),
        },
    }
    manifest_out.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(json.dumps(manifest, indent=2))
    print("\nTop cell types by median reduced-model score:")
    print(risk_summary.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
