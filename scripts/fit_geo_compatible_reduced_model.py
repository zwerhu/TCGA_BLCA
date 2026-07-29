#!/usr/bin/env python
"""
Fit a GEO-compatible reduced TCGA-BLCA Cox model.

The candidate genes are restricted to genes present in TCGA, GSE13507, and
GSE31684. The model is trained only in TCGA with a stratified train/test split.
GEO validation uses the fixed TCGA coefficients on cohort-wise z-scored
expression so microarray/RNA-seq platform offsets do not dominate the score.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fit_tcga_lasso_cox import (  # noqa: E402
    SEED,
    bh_fdr,
    concordance_index,
    cox_negloglik_grad,
    cox_score_screen,
    cv_lambdas,
    fit_l1_cox,
    stratified_split,
)

PROCESSED = ROOT / "data" / "processed"
EXT_PROCESSED = ROOT / "data" / "external" / "processed"
TABLES = ROOT / "results" / "tables"
MODEL_DIR = ROOT / "results" / "model"

TRAIN_FRACTION = 0.70
MAX_SCREEN_GENES = 10000
MAX_LASSO_CANDIDATES = 220
MIN_FINAL_GENES = 4
MAX_FINAL_GENES = 8


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)


def load_tcga() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(PROCESSED / "blca_log2_tpm_gene_symbol_tumor.tsv.gz", sep="\t", index_col=0)
    survival = pd.read_csv(PROCESSED / "blca_survival_expression_sample_sheet.tsv", sep="\t")
    survival = survival[survival["matrix_column"].isin(expr.columns)].copy()
    survival = survival[survival["os_time_days"].notna() & survival["os_event"].notna()].copy()
    survival = survival.drop_duplicates("matrix_column").sort_values("matrix_column")
    expr = expr.loc[:, survival["matrix_column"].tolist()]
    return expr, survival


def load_geo_gene_sets() -> tuple[set[str], set[str]]:
    g135 = pd.read_csv(EXT_PROCESSED / "GSE13507_gene_symbol_expression.tsv.gz", sep="\t", index_col=0, nrows=1)
    g316 = pd.read_csv(EXT_PROCESSED / "GSE31684_gene_symbol_expression.tsv.gz", sep="\t", index_col=0, nrows=1)
    # nrows=1 gives only first row; read index separately for correct gene universe.
    g135_full = pd.read_csv(EXT_PROCESSED / "GSE13507_gene_symbol_expression.tsv.gz", sep="\t", usecols=[0])
    g316_full = pd.read_csv(EXT_PROCESSED / "GSE31684_gene_symbol_expression.tsv.gz", sep="\t", usecols=[0])
    return set(g135_full.iloc[:, 0].astype(str)), set(g316_full.iloc[:, 0].astype(str))


def select_common_genes(tcga_expr: pd.DataFrame) -> list[str]:
    g135, g316 = load_geo_gene_sets()
    common = sorted(set(tcga_expr.index) & g135 & g316)
    annotation_path = PROCESSED / "blca_gene_annotation.tsv"
    if annotation_path.exists():
        anno = pd.read_csv(annotation_path, sep="\t")
        protein_symbols = set(anno.loc[anno["gene_type"].eq("protein_coding"), "gene_name"].astype(str))
        common = [g for g in common if g in protein_symbols]
    return common


def train_model(tcga_expr: pd.DataFrame, survival: pd.DataFrame, common_genes: list[str]) -> dict[str, object]:
    time = survival["os_time_days"].to_numpy(dtype=float)
    event = survival["os_event"].to_numpy(dtype=int)
    train_idx, test_idx = stratified_split(event, TRAIN_FRACTION, SEED + 17)
    split_label = np.array(["test"] * len(survival), dtype=object)
    split_label[train_idx] = "train"

    X_raw_all = tcga_expr.loc[common_genes].to_numpy(dtype=float).T
    train_raw = X_raw_all[train_idx]
    gene_mean = np.nanmean(train_raw, axis=0)
    gene_sd = np.nanstd(train_raw, axis=0, ddof=1)
    ok = np.isfinite(gene_mean) & np.isfinite(gene_sd) & (gene_mean > 0.2) & (gene_sd > 0.1)
    ok_idx = np.where(ok)[0]
    variable_order = np.argsort(-gene_sd[ok_idx])
    screen_idx = ok_idx[variable_order[: min(MAX_SCREEN_GENES, len(variable_order))]]
    screen_genes = np.array(common_genes)[screen_idx]
    X_screen_train_raw = X_raw_all[np.ix_(train_idx, screen_idx)]
    screen_mean = X_screen_train_raw.mean(axis=0)
    screen_sd = X_screen_train_raw.std(axis=0, ddof=1)
    X_screen_train = (X_screen_train_raw - screen_mean) / np.maximum(screen_sd, 1e-8)
    z, p, info = cox_score_screen(X_screen_train, time[train_idx], event[train_idx])
    screen_df = pd.DataFrame(
        {
            "gene": screen_genes,
            "z": z,
            "pvalue": p,
            "FDR": bh_fdr(p),
            "score_information": info,
            "train_mean_log2tpm": screen_mean,
            "train_sd_log2tpm": screen_sd,
        }
    ).sort_values(["pvalue", "FDR"])
    screen_df.to_csv(TABLES / "reduced_common_gene_train_cox_score_screen.tsv", sep="\t", index=False)

    candidate_genes = screen_df.head(MAX_LASSO_CANDIDATES)["gene"].tolist()
    candidate_idx = [common_genes.index(g) for g in candidate_genes]
    X_train_raw = X_raw_all[np.ix_(train_idx, candidate_idx)]
    train_mean = X_train_raw.mean(axis=0)
    train_sd = np.maximum(X_train_raw.std(axis=0, ddof=1), 1e-8)
    X_train = (X_train_raw - train_mean) / train_sd
    X_all = (X_raw_all[:, candidate_idx] - train_mean) / train_sd

    _, grad0 = cox_negloglik_grad(np.zeros(X_train.shape[1]), X_train, time[train_idx], event[train_idx])
    lambda_max = float(np.max(np.abs(grad0)))
    lambdas = np.exp(np.linspace(math.log(lambda_max), math.log(lambda_max * 0.035), 34))
    cv_df = cv_lambdas(X_train, time[train_idx], event[train_idx], lambdas, k=5)
    cv_df.to_csv(MODEL_DIR / "reduced_geo_compatible_lasso_cv.tsv", sep="\t", index=False)

    best = cv_df.loc[cv_df["mean_cv_negloglik"].idxmin()]
    threshold = best["mean_cv_negloglik"] + best["se_cv_negloglik"]
    lambda_order = cv_df[cv_df["mean_cv_negloglik"] <= threshold].sort_values("lambda", ascending=False)["lambda"].tolist()
    if not lambda_order:
        lambda_order = cv_df.sort_values("mean_cv_negloglik")["lambda"].tolist()

    chosen_lambda = float(best["lambda"])
    chosen_beta = None
    chosen_nz = None
    for lam in lambda_order + cv_df.sort_values("lambda", ascending=False)["lambda"].tolist():
        beta, _ = fit_l1_cox(X_train, time[train_idx], event[train_idx], float(lam), max_iter=700)
        nz = int(np.sum(np.abs(beta) > 1e-6))
        if MIN_FINAL_GENES <= nz <= MAX_FINAL_GENES:
            chosen_lambda = float(lam)
            chosen_beta = beta
            chosen_nz = nz
            break
    if chosen_beta is None:
        candidates = []
        for lam in cv_df["lambda"]:
            beta, _ = fit_l1_cox(X_train, time[train_idx], event[train_idx], float(lam), max_iter=700)
            nz = int(np.sum(np.abs(beta) > 1e-6))
            target_distance = min(abs(nz - MIN_FINAL_GENES), abs(nz - MAX_FINAL_GENES))
            candidates.append((target_distance, abs(float(lam) - float(best["lambda"])), float(lam), beta, nz))
        candidates.sort(key=lambda x: (x[0], x[1]))
        chosen_lambda = candidates[0][2]
        chosen_beta = candidates[0][3]
        chosen_nz = candidates[0][4]

    beta = chosen_beta
    nz_mask = np.abs(beta) > 1e-6
    final_genes = np.array(candidate_genes)[nz_mask]
    final_beta = beta[nz_mask]
    final_mean = train_mean[nz_mask]
    final_sd = train_sd[nz_mask]

    # If L1 gives too many because of numerical ties, keep strongest 8.
    if len(final_genes) > MAX_FINAL_GENES:
        strongest = np.argsort(-np.abs(final_beta))[:MAX_FINAL_GENES]
        final_genes = final_genes[strongest]
        final_beta = final_beta[strongest]
        final_mean = final_mean[strongest]
        final_sd = final_sd[strongest]
        nz_mask = np.isin(np.array(candidate_genes), final_genes)

    final_df = pd.DataFrame(
        {
            "gene": final_genes,
            "coefficient": final_beta,
            "tcga_train_mean_log2_tpm": final_mean,
            "tcga_train_sd_log2_tpm": final_sd,
            "geo_validation_transform": "cohort-wise z-score for this gene, then multiply by coefficient",
        }
    ).sort_values("coefficient", ascending=False)
    final_df.to_csv(MODEL_DIR / "reduced_geo_compatible_model.tsv", sep="\t", index=False)

    final_candidate_positions = [candidate_genes.index(g) for g in final_genes]
    risk_all = X_all[:, final_candidate_positions] @ final_beta
    cutoff = float(np.median(risk_all[train_idx]))
    risk_group = np.where(risk_all >= cutoff, "High", "Low")
    score_df = survival.copy().reset_index(drop=True)
    score_df["split"] = split_label
    score_df["reduced_risk_score"] = risk_all
    score_df["reduced_risk_group"] = risk_group
    score_df.to_csv(TABLES / "tcga_reduced_geo_compatible_risk_scores.tsv", sep="\t", index=False)

    metrics = []
    for label, idx in [("train", train_idx), ("test", test_idx), ("all", np.arange(len(survival)))]:
        metrics.append(
            {
                "cohort": f"TCGA-{label}",
                "n": len(idx),
                "events": int(event[idx].sum()),
                "c_index": concordance_index(time[idx], event[idx], risk_all[idx]),
                "risk_cutoff_train_median": cutoff,
                "lambda": chosen_lambda,
                "n_final_genes": len(final_genes),
            }
        )
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(TABLES / "tcga_reduced_geo_compatible_metrics.tsv", sep="\t", index=False)

    split_out = score_df[["matrix_column", "sample_barcode", "patient_barcode", "split", "os_time_days", "os_event"]]
    split_out.to_csv(MODEL_DIR / "reduced_geo_compatible_train_test_split.tsv", sep="\t", index=False)
    metadata = {
        "seed": SEED + 17,
        "train_fraction": TRAIN_FRACTION,
        "common_genes_total": len(common_genes),
        "common_genes_screened": len(screen_idx),
        "candidate_genes": len(candidate_genes),
        "lambda": chosen_lambda,
        "n_nonzero_before_optional_trim": chosen_nz,
        "risk_cutoff_tcga_train_median": cutoff,
        "final_genes": final_df.to_dict(orient="records"),
        "tcga_metrics": metrics,
        "note": "Trained only in TCGA on genes shared by TCGA/GSE13507/GSE31684. GEO validation uses cohort-wise z-scored expression.",
    }
    (MODEL_DIR / "reduced_geo_compatible_model.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Common genes: {len(common_genes)}")
    print(f"Screened genes: {len(screen_idx)}")
    print(f"Final genes ({len(final_genes)}): {', '.join(final_genes)}")
    print(metrics_df.to_string(index=False))
    return {
        "final_genes": final_genes.tolist(),
        "final_beta": final_beta.tolist(),
        "tcga_cutoff": cutoff,
        "tcga_metrics": metrics,
    }


def geo_score_one(cohort: str, final_model: pd.DataFrame) -> pd.DataFrame:
    expr = pd.read_csv(EXT_PROCESSED / f"{cohort}_gene_symbol_expression.tsv.gz", sep="\t", index_col=0)
    clinical = pd.read_csv(EXT_PROCESSED / f"{cohort}_clinical_survival.tsv", sep="\t")
    genes = final_model["gene"].tolist()
    missing = [g for g in genes if g not in expr.index]
    available = [g for g in genes if g in expr.index]
    expr = expr.loc[available, [c for c in clinical["sample_id"] if c in expr.columns]]
    means = expr.mean(axis=1)
    sds = expr.std(axis=1, ddof=1).replace(0, np.nan)
    z_expr = expr.sub(means, axis=0).div(sds, axis=0).fillna(0.0)
    coeff = final_model.set_index("gene").loc[available, "coefficient"]
    risk = z_expr.T @ coeff
    out = clinical[clinical["sample_id"].isin(risk.index)].copy()
    out["reduced_risk_score"] = out["sample_id"].map(risk.to_dict())
    cutoff = float(np.median(out["reduced_risk_score"]))
    out["reduced_risk_group"] = np.where(out["reduced_risk_score"] >= cutoff, "High", "Low")
    out["geo_median_cutoff"] = cutoff
    out["n_reduced_model_genes_available"] = len(available)
    out["n_reduced_model_genes_missing"] = len(missing)
    out["missing_reduced_model_genes"] = ";".join(missing)
    out.to_csv(EXT_PROCESSED / f"{cohort}_reduced_geo_compatible_risk_scores.tsv", sep="\t", index=False)
    return out


def score_geo() -> None:
    final_model = pd.read_csv(MODEL_DIR / "reduced_geo_compatible_model.tsv", sep="\t")
    rows = []
    for cohort in ["GSE13507", "GSE31684"]:
        out = geo_score_one(cohort, final_model)
        time = out["os_time_days"].to_numpy(dtype=float)
        event = out["os_event"].to_numpy(dtype=int)
        risk = out["reduced_risk_score"].to_numpy(dtype=float)
        rows.append(
            {
                "cohort": cohort,
                "n": len(out),
                "events": int(np.nansum(event)),
                "c_index": concordance_index(time, event, risk),
                "available_model_genes": out["n_reduced_model_genes_available"].iloc[0],
                "missing_model_genes": out["missing_reduced_model_genes"].iloc[0],
            }
        )
    pd.DataFrame(rows).to_csv(TABLES / "geo_reduced_compatible_python_metrics.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def main() -> int:
    ensure_dirs()
    tcga_expr, survival = load_tcga()
    common_genes = select_common_genes(tcga_expr)
    train_model(tcga_expr, survival, common_genes)
    score_geo()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
