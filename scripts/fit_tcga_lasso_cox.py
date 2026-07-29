#!/usr/bin/env python
"""
Train/test split plus L1-penalized Cox model for TCGA-BLCA.

This script avoids external ML dependencies so the model can be reproduced in a
clean R/Python installation. It performs:

1. Stratified train/test split by OS event.
2. Train-only univariate Cox score screening.
3. L1-penalized Cox optimization by proximal gradient.
4. Cross-validation over lambda values.
5. Final coefficient table and fixed risk score for all TCGA tumor samples.

For manuscript work, compare this output with glmnet when an R package
environment is available. The fixed formula produced here is suitable for
external validation because gene coefficients and training centering/scaling are
saved explicitly.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
MODEL_DIR = ROOT / "results" / "model"

SEED = 20260728
TRAIN_FRACTION = 0.70
MAX_SCREEN_GENES = 12000
MAX_LASSO_CANDIDATES = 250
MIN_FINAL_GENES = 4
MAX_FINAL_GENES = 15


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)


def stratified_split(events: np.ndarray, train_fraction: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    train_idx = []
    test_idx = []
    for event_value in sorted(np.unique(events)):
        idx = np.where(events == event_value)[0]
        rng.shuffle(idx)
        n_train = int(round(len(idx) * train_fraction))
        train_idx.extend(idx[:n_train])
        test_idx.extend(idx[n_train:])
    train_idx = np.array(sorted(train_idx), dtype=int)
    test_idx = np.array(sorted(test_idx), dtype=int)
    return train_idx, test_idx


def cox_score_screen(X: np.ndarray, time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cox score test at beta=0 for all columns of X."""
    order = np.argsort(-time)
    t = time[order]
    e = event[order].astype(bool)
    Xs = X[order, :]
    csum = np.cumsum(Xs, axis=0)
    csum2 = np.cumsum(Xs * Xs, axis=0)
    U = np.zeros(X.shape[1], dtype=float)
    I = np.zeros(X.shape[1], dtype=float)

    event_times = np.unique(t[e])
    for et in event_times:
        risk_end = np.searchsorted(-t, -et, side="right") - 1
        ev_mask = (t == et) & e
        d = int(ev_mask.sum())
        if d == 0 or risk_end < 0:
            continue
        risk_n = risk_end + 1
        risk_sum = csum[risk_end, :]
        risk_sum2 = csum2[risk_end, :]
        risk_mean = risk_sum / risk_n
        risk_var = np.maximum(risk_sum2 / risk_n - risk_mean * risk_mean, 1e-12)
        event_sum = Xs[ev_mask, :].sum(axis=0)
        U += event_sum - d * risk_mean
        I += d * risk_var

    z = U / np.sqrt(np.maximum(I, 1e-12))
    # two-sided normal p-value using erfc, vectorized without scipy.
    p = np.array([math.erfc(abs(float(v)) / math.sqrt(2.0)) for v in z])
    return z, p, I


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    n = len(vals)
    order = np.argsort(vals)
    ranked = vals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty_like(vals)
    tmp[order] = q
    out[ok] = tmp
    return out


def cox_negloglik_grad(beta: np.ndarray, X: np.ndarray, time: np.ndarray, event: np.ndarray) -> tuple[float, np.ndarray]:
    eta = X @ beta
    eta = np.clip(eta, -50, 50)
    exp_eta = np.exp(eta)
    unique_event_times = np.unique(time[event == 1])
    loglik = 0.0
    grad_eta = np.zeros(X.shape[0], dtype=float)

    for et in unique_event_times:
        event_mask = (time == et) & (event == 1)
        risk_mask = time >= et
        d = int(event_mask.sum())
        if d == 0:
            continue
        risk_sum = exp_eta[risk_mask].sum()
        loglik += eta[event_mask].sum() - d * math.log(max(risk_sum, 1e-300))
        grad_eta[event_mask] -= 1.0
        grad_eta[risk_mask] += d * exp_eta[risk_mask] / max(risk_sum, 1e-300)

    n_events = max(int(event.sum()), 1)
    f = -loglik / n_events
    grad = X.T @ grad_eta / n_events
    return f, grad


def soft_threshold(x: np.ndarray, threshold: float) -> np.ndarray:
    return np.sign(x) * np.maximum(np.abs(x) - threshold, 0.0)


def fit_l1_cox(
    X: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    lam: float,
    max_iter: int = 350,
    tol: float = 1e-6,
    beta_init: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
    beta = np.zeros(X.shape[1], dtype=float) if beta_init is None else beta_init.copy()
    L = 1.0
    previous_obj = np.inf
    for _ in range(max_iter):
        f, grad = cox_negloglik_grad(beta, X, time, event)
        while True:
            candidate = soft_threshold(beta - grad / L, lam / L)
            f_new, _ = cox_negloglik_grad(candidate, X, time, event)
            diff = candidate - beta
            rhs = f + grad.dot(diff) + 0.5 * L * diff.dot(diff)
            if f_new <= rhs + 1e-10 or L > 1e8:
                break
            L *= 2.0
        beta = candidate
        obj = f_new + lam * np.abs(beta).sum()
        L = max(L * 0.85, 1e-4)
        if abs(previous_obj - obj) < tol * (1.0 + abs(previous_obj)):
            break
        previous_obj = obj
    f_final, _ = cox_negloglik_grad(beta, X, time, event)
    return beta, f_final + lam * np.abs(beta).sum()


def make_folds(events: np.ndarray, k: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    folds = [[] for _ in range(k)]
    for event_value in sorted(np.unique(events)):
        idx = np.where(events == event_value)[0]
        rng.shuffle(idx)
        for i, sample_idx in enumerate(idx):
            folds[i % k].append(int(sample_idx))
    return [np.array(sorted(fold), dtype=int) for fold in folds]


def cv_lambdas(X: np.ndarray, time: np.ndarray, event: np.ndarray, lambdas: np.ndarray, k: int = 5) -> pd.DataFrame:
    folds = make_folds(event, k=k, seed=SEED + 101)
    rows = []
    for lam in lambdas:
        losses = []
        nonzeros = []
        for val_idx in folds:
            train_idx = np.setdiff1d(np.arange(X.shape[0]), val_idx)
            beta, _ = fit_l1_cox(X[train_idx], time[train_idx], event[train_idx], lam)
            val_loss, _ = cox_negloglik_grad(beta, X[val_idx], time[val_idx], event[val_idx])
            losses.append(val_loss)
            nonzeros.append(int(np.sum(np.abs(beta) > 1e-7)))
        rows.append(
            {
                "lambda": lam,
                "mean_cv_negloglik": float(np.mean(losses)),
                "se_cv_negloglik": float(np.std(losses, ddof=1) / math.sqrt(k)),
                "mean_nonzero": float(np.mean(nonzeros)),
            }
        )
        print(f"lambda={lam:.5g} cv_loss={rows[-1]['mean_cv_negloglik']:.5f} nz={rows[-1]['mean_nonzero']:.1f}", flush=True)
    return pd.DataFrame(rows)


def concordance_index(time: np.ndarray, event: np.ndarray, score: np.ndarray) -> float:
    permissible = 0
    concordant = 0.0
    n = len(time)
    for i in range(n):
        if event[i] != 1:
            continue
        for j in range(n):
            if time[j] <= time[i]:
                continue
            permissible += 1
            if score[i] > score[j]:
                concordant += 1
            elif score[i] == score[j]:
                concordant += 0.5
    return float(concordant / permissible) if permissible else float("nan")


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(PROCESSED / "blca_log2_tpm_gene_symbol_tumor.tsv.gz", sep="\t", index_col=0)
    survival = pd.read_csv(PROCESSED / "blca_survival_expression_sample_sheet.tsv", sep="\t")
    survival = survival[survival["matrix_column"].isin(expr.columns)].copy()
    survival = survival[survival["os_time_days"].notna() & survival["os_event"].notna()].copy()
    survival = survival.drop_duplicates("matrix_column")
    survival = survival.sort_values("matrix_column")
    expr = expr.loc[:, survival["matrix_column"].tolist()]
    return expr, survival


def main() -> int:
    ensure_dirs()
    expr_df, survival = load_data()
    time = survival["os_time_days"].to_numpy(dtype=float)
    event = survival["os_event"].to_numpy(dtype=int)
    train_idx, test_idx = stratified_split(event, TRAIN_FRACTION, SEED)

    split_rows = []
    split_label = np.array(["test"] * len(survival), dtype=object)
    split_label[train_idx] = "train"
    for i, row in survival.reset_index(drop=True).iterrows():
        split_rows.append(
            {
                "matrix_column": row["matrix_column"],
                "sample_barcode": row["sample_barcode"],
                "patient_barcode": row["patient_barcode"],
                "split": split_label[i],
                "os_time_days": row["os_time_days"],
                "os_event": row["os_event"],
            }
        )
    pd.DataFrame(split_rows).to_csv(MODEL_DIR / "tcga_train_test_split.tsv", sep="\t", index=False)

    X_all_raw = expr_df.to_numpy(dtype=float).T
    means_all = np.nanmean(X_all_raw[train_idx], axis=0)
    sds_all = np.nanstd(X_all_raw[train_idx], axis=0, ddof=1)
    ok = np.isfinite(means_all) & np.isfinite(sds_all) & (means_all > 0.2) & (sds_all > 0.1)
    variable_order = np.argsort(-sds_all[ok])
    ok_gene_indices = np.where(ok)[0][variable_order[:MAX_SCREEN_GENES]]
    genes_screen = expr_df.index.to_numpy()[ok_gene_indices]

    X_screen_train_raw = X_all_raw[np.ix_(train_idx, ok_gene_indices)]
    screen_mean = X_screen_train_raw.mean(axis=0)
    screen_sd = X_screen_train_raw.std(axis=0, ddof=1)
    X_screen_train = (X_screen_train_raw - screen_mean) / np.maximum(screen_sd, 1e-8)
    z, p, info = cox_score_screen(X_screen_train, time[train_idx], event[train_idx])
    fdr = bh_fdr(p)
    screen_df = pd.DataFrame(
        {
            "gene": genes_screen,
            "z": z,
            "pvalue": p,
            "FDR": fdr,
            "score_information": info,
            "train_mean_log2tpm": screen_mean,
            "train_sd_log2tpm": screen_sd,
        }
    ).sort_values(["pvalue", "FDR"])
    screen_df.to_csv(TABLES / "tcga_train_univariate_cox_score_screen.tsv", sep="\t", index=False)

    candidate_genes = screen_df.head(MAX_LASSO_CANDIDATES)["gene"].tolist()
    candidate_idx = [expr_df.index.get_loc(g) for g in candidate_genes]
    X_train_raw = X_all_raw[np.ix_(train_idx, candidate_idx)]
    X_test_raw = X_all_raw[np.ix_(test_idx, candidate_idx)]
    train_mean = X_train_raw.mean(axis=0)
    train_sd = X_train_raw.std(axis=0, ddof=1)
    train_sd = np.maximum(train_sd, 1e-8)
    X_train = (X_train_raw - train_mean) / train_sd
    X_test = (X_test_raw - train_mean) / train_sd
    X_all = (X_all_raw[:, candidate_idx] - train_mean) / train_sd

    f0, grad0 = cox_negloglik_grad(np.zeros(X_train.shape[1]), X_train, time[train_idx], event[train_idx])
    lambda_max = float(np.max(np.abs(grad0)))
    lambdas = np.exp(np.linspace(math.log(lambda_max), math.log(lambda_max * 0.025), 26))
    cv_df = cv_lambdas(X_train, time[train_idx], event[train_idx], lambdas, k=5)
    cv_df.to_csv(MODEL_DIR / "tcga_lasso_cox_cv.tsv", sep="\t", index=False)

    best = cv_df.loc[cv_df["mean_cv_negloglik"].idxmin()]
    threshold = best["mean_cv_negloglik"] + best["se_cv_negloglik"]
    eligible_lambdas = cv_df[cv_df["mean_cv_negloglik"] <= threshold].copy()
    eligible_lambdas = eligible_lambdas.sort_values("lambda", ascending=False)

    chosen_lambda = float(best["lambda"])
    chosen_beta = None
    for lam in eligible_lambdas["lambda"]:
        beta, _ = fit_l1_cox(X_train, time[train_idx], event[train_idx], float(lam), max_iter=600)
        nz = int(np.sum(np.abs(beta) > 1e-6))
        if MIN_FINAL_GENES <= nz <= MAX_FINAL_GENES:
            chosen_lambda = float(lam)
            chosen_beta = beta
            break
    if chosen_beta is None:
        # Fall back to the lambda with gene count closest to the desired range.
        beta_records = []
        for lam in cv_df["lambda"]:
            beta, _ = fit_l1_cox(X_train, time[train_idx], event[train_idx], float(lam), max_iter=600)
            nz = int(np.sum(np.abs(beta) > 1e-6))
            target = min(abs(nz - MIN_FINAL_GENES), abs(nz - MAX_FINAL_GENES))
            beta_records.append((target, abs(float(lam) - float(best["lambda"])), float(lam), beta, nz))
        beta_records.sort(key=lambda x: (x[0], x[1]))
        chosen_lambda = beta_records[0][2]
        chosen_beta = beta_records[0][3]

    beta = chosen_beta
    nz_mask = np.abs(beta) > 1e-6
    final_genes = np.array(candidate_genes)[nz_mask]
    final_beta = beta[nz_mask]
    final_mean = train_mean[nz_mask]
    final_sd = train_sd[nz_mask]

    model_df = pd.DataFrame(
        {
            "gene": final_genes,
            "coefficient": final_beta,
            "train_mean_log2_tpm": final_mean,
            "train_sd_log2_tpm": final_sd,
            "standardization": "(expression_log2_tpm - train_mean_log2_tpm) / train_sd_log2_tpm",
        }
    ).sort_values("coefficient", ascending=False)
    model_df.to_csv(MODEL_DIR / "final_lasso_cox_model.tsv", sep="\t", index=False)

    risk_all = X_all[:, nz_mask] @ final_beta
    risk_train = risk_all[train_idx]
    cutoff = float(np.median(risk_train))
    risk_group = np.where(risk_all >= cutoff, "High", "Low")
    score_df = survival.copy().reset_index(drop=True)
    score_df["split"] = split_label
    score_df["lasso_risk_score"] = risk_all
    score_df["lasso_risk_group"] = risk_group
    score_df.to_csv(TABLES / "tcga_lasso_risk_scores.tsv", sep="\t", index=False)

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
    pd.DataFrame(metrics).to_csv(TABLES / "tcga_lasso_model_metrics.tsv", sep="\t", index=False)

    formula = " + ".join(
        [
            f"({coef:.10g}) * (({gene} - {mean:.10g}) / {sd:.10g})"
            for gene, coef, mean, sd in zip(final_genes, final_beta, final_mean, final_sd)
        ]
    )
    metadata = {
        "seed": SEED,
        "train_fraction": TRAIN_FRACTION,
        "max_screen_genes": MAX_SCREEN_GENES,
        "max_lasso_candidates": MAX_LASSO_CANDIDATES,
        "lambda": chosen_lambda,
        "risk_cutoff_train_median": cutoff,
        "formula": f"risk_score = {formula}",
    }
    (MODEL_DIR / "final_lasso_cox_model.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(f"train/test: {len(train_idx)}/{len(test_idx)}")
    print(f"chosen lambda: {chosen_lambda:.6g}")
    print(f"final genes ({len(final_genes)}): {', '.join(final_genes)}")
    print(f"test c-index: {metrics[1]['c_index']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
