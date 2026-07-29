#!/usr/bin/env python3
"""Apply the fixed reduced TCGA-BLCA risk model to ICB cohorts.

Outputs harmonized sample-level tables for IMvigor210 and GSE176307, plus
response/survival/signature summary tables used by the R plotting script.
"""

from __future__ import annotations

import gzip
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
IMMUNE_DIR = ROOT / "data" / "external" / "immunotherapy"
PROCESSED_DIR = IMMUNE_DIR / "processed"
RESULT_TABLE_DIR = ROOT / "results" / "tables"
MODEL_PATH = ROOT / "results" / "model" / "reduced_geo_compatible_model.tsv"


MODULES = {
    "CD8_effector": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "IFNG", "CXCL9", "CXCL10"],
    "Checkpoint": ["CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1"],
    "TGFb_CAF": ["TGFB1", "TGFBI", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "FAP", "PDGFRA", "PDGFRB"],
    "EMT_stromal": ["EMP1", "AHNAK", "VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "COL1A1"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "CLDN5"],
    "Myeloid": ["LYZ", "AIF1", "CD68", "C1QA", "C1QB", "C1QC"],
    "B_cell": ["MS4A1", "CD79A"],
    "Basal_epithelial": ["KRT5", "KRT14"],
    "Luminal_epithelial": ["KRT20", "EPCAM"],
    "Proliferation": ["MKI67", "TOP2A"],
}


def ensure_dirs() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_TABLE_DIR.mkdir(parents=True, exist_ok=True)


def clean_key(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace({"NA": np.nan, "": np.nan, "nan": np.nan}), errors="coerce")


def load_model() -> pd.DataFrame:
    model = pd.read_csv(MODEL_PATH, sep="\t")
    if not {"gene", "coefficient"}.issubset(model.columns):
        raise ValueError(f"Unexpected model columns in {MODEL_PATH}: {list(model.columns)}")
    model = model[["gene", "coefficient"]].copy()
    model["coefficient"] = pd.to_numeric(model["coefficient"], errors="raise")
    return model


def aggregate_duplicate_genes(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = expr.index.astype(str)
    if expr.index.has_duplicates:
        expr = expr.groupby(expr.index).mean(numeric_only=True)
    return expr


def gene_zscore(expr: pd.DataFrame) -> pd.DataFrame:
    means = expr.mean(axis=1)
    sds = expr.std(axis=1, ddof=0).replace(0, np.nan)
    z = expr.sub(means, axis=0).div(sds, axis=0)
    return z.fillna(0.0)


def score_expression(expr: pd.DataFrame, model: pd.DataFrame) -> tuple[pd.Series, list[str], pd.DataFrame]:
    expr = aggregate_duplicate_genes(expr)
    present = [gene for gene in model["gene"] if gene in expr.index]
    if not present:
        raise ValueError("No model genes present in expression matrix")
    z = gene_zscore(expr.loc[present])
    coef = model.set_index("gene").loc[present, "coefficient"]
    score = z.mul(coef, axis=0).sum(axis=0)
    contrib = z.mul(coef, axis=0)
    return score, present, contrib


def module_scores(expr: pd.DataFrame, modules: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    expr = aggregate_duplicate_genes(expr)
    all_genes = sorted({gene for genes in modules.values() for gene in genes if gene in expr.index})
    z = gene_zscore(expr.loc[all_genes]) if all_genes else pd.DataFrame(index=[], columns=expr.columns)
    scores = {}
    present = {}
    for module, genes in modules.items():
        genes_present = [gene for gene in genes if gene in z.index]
        present[module] = genes_present
        if genes_present:
            scores[f"sig_{module}"] = z.loc[genes_present].mean(axis=0)
        else:
            scores[f"sig_{module}"] = pd.Series(np.nan, index=expr.columns)
    return pd.DataFrame(scores), present


def read_imvigor(model: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = IMMUNE_DIR / "IMvigor210" / "processed"
    clinical = pd.read_csv(base / "IMvigor210_clinical.tsv", sep="\t")
    counts = pd.read_csv(base / "IMvigor210_selected_gene_counts.tsv", sep="\t")
    counts = counts.rename(columns={counts.columns[0]: "gene"})
    counts = counts.set_index("gene").apply(pd.to_numeric, errors="coerce")

    clinical["sample_id"] = clinical["sample_id"].astype(str)
    clinical = clinical.set_index("sample_id", drop=False)
    counts = counts.loc[:, [col for col in counts.columns if col in clinical.index]]
    size_factor = numeric(clinical.loc[counts.columns, "sizeFactor"])
    size_factor = size_factor.fillna(size_factor.median())
    size_factor = size_factor.replace(0, size_factor.median())
    expr = np.log2(counts.div(size_factor, axis=1) + 1.0)

    risk, present, contrib = score_expression(expr, model)
    sigs, module_present = module_scores(expr, MODULES)

    out = pd.DataFrame(index=expr.columns)
    out["cohort"] = "IMvigor210"
    out["analysis_id"] = out.index
    out["sample_id"] = out.index
    out["patient_id"] = clinical.loc[out.index, "ANONPT_ID"].astype(str).values
    out["aliquot_id"] = out.index
    out["geo_accession"] = ""
    out["response_raw"] = clinical.loc[out.index, "Best Confirmed Overall Response"].astype(str).replace("nan", np.nan).values
    out["response_group"] = clinical.loc[out.index, "binaryResponse"].map({"CR/PR": "Responder", "SD/PD": "Non-responder"}).values
    out["response_binary"] = out["response_group"].map({"Responder": 1, "Non-responder": 0})
    out["nonresponse_binary"] = out["response_group"].map({"Responder": 0, "Non-responder": 1})
    out["os_time_months"] = numeric(clinical.loc[out.index, "os"]).values
    out["os_event"] = numeric(clinical.loc[out.index, "censOS"]).values
    out["pfs_time_months"] = np.nan
    out["pfs_event"] = np.nan
    out["age"] = numeric(clinical.loc[out.index, "Sample age"]).values
    out["sex"] = clinical.loc[out.index, "Sex"].astype(str).replace("nan", np.nan).values
    out["ecog"] = numeric(clinical.loc[out.index, "Baseline ECOG Score"]).values
    out["tmb"] = numeric(clinical.loc[out.index, "FMOne mutation burden per MB"]).values
    out["tmb_group"] = pd.Series(np.where(out["tmb"] >= 10, "TMB-high", "TMB-low"), index=out.index, dtype=object)
    out.loc[out["tmb"].isna(), "tmb_group"] = np.nan
    out["fgfr3_altered"] = np.nan
    out["treatment"] = "Atezolizumab"
    out["immune_phenotype"] = clinical.loc[out.index, "Immune phenotype"].astype(str).replace("nan", np.nan).values
    out["tcga_subtype"] = clinical.loc[out.index, "TCGA Subtype"].astype(str).replace("nan", np.nan).values
    out["risk_score"] = risk.loc[out.index].values
    out["model_genes_present"] = ",".join(present)
    out["n_model_genes_present"] = len(present)
    out["n_model_genes_total"] = len(model)
    out = pd.concat([out, sigs.loc[out.index]], axis=1)

    expr_out = expr.loc[sorted(set(model["gene"]).union(*MODULES.values()) & set(expr.index)), out.index].copy()
    expr_out.columns = [f"IMvigor210::{sample}" for sample in expr_out.columns]
    out.index = [f"IMvigor210::{sample}" for sample in out.index]
    return out.reset_index(drop=True), expr_out


def read_geo_series_samples(path: Path) -> pd.DataFrame:
    sample_rows: list[tuple[str, list[str]]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0]
            values = [part.strip().strip('"') for part in parts[1:]]
            sample_rows.append((key, values))

    title_rows = [values for key, values in sample_rows if key == "!Sample_title"]
    if not title_rows:
        raise ValueError(f"No !Sample_title row in {path}")
    n = len(title_rows[0])
    records = [dict() for _ in range(n)]
    direct_map = {
        "!Sample_title": "sample_title",
        "!Sample_geo_accession": "geo_accession",
        "!Sample_description": "aliquot_id",
        "!Sample_source_name_ch1": "source_name",
    }
    for key, values in sample_rows:
        if len(values) != n:
            continue
        if key == "!Sample_characteristics_ch1":
            for idx, value in enumerate(values):
                if ":" not in value:
                    continue
                k, v = value.split(":", 1)
                records[idx][clean_key(k)] = v.strip()
        elif key in direct_map:
            out_key = direct_map[key]
            for idx, value in enumerate(values):
                records[idx][out_key] = value

    df = pd.DataFrame(records)
    sample_ids = []
    for title in df["sample_title"].fillna(""):
        match = re.search(r"(BACI\d+)", title)
        sample_ids.append(match.group(1) if match else title.replace(" ", "_"))
    df["sample_id"] = sample_ids
    df["patient_id"] = df["sample_id"]

    counts = Counter(df["sample_id"])
    df["analysis_id"] = [
        sid if counts[sid] == 1 else f"{sid}_{aliquot}"
        for sid, aliquot in zip(df["sample_id"], df["aliquot_id"])
    ]
    return df


def read_selected_salmon_tpm(path: Path, genes: set[str]) -> pd.DataFrame:
    rows: dict[str, list[np.ndarray]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_cols = header[1:]
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            gene = parts[0]
            if gene not in genes:
                continue
            vals = pd.to_numeric(pd.Series(parts[1:]), errors="coerce").to_numpy(dtype=float)
            rows[gene].append(vals)

    if not rows:
        raise ValueError(f"No target genes found in {path}")
    data = {}
    for gene, matrices in rows.items():
        data[gene] = np.nanmean(np.vstack(matrices), axis=0)
    return pd.DataFrame(data, index=sample_cols).T


def read_gse176307(model: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = IMMUNE_DIR / "GSE176307"
    clinical = read_geo_series_samples(base / "GSE176307_series_matrix.txt.gz")
    target_genes = set(model["gene"])
    for genes in MODULES.values():
        target_genes.update(genes)
    tpm = read_selected_salmon_tpm(base / "GSE176307_salmon_tpm_gene.matrix.tsv.gz", target_genes)

    clinical = clinical[clinical["aliquot_id"].isin(tpm.columns)].copy()
    rs_to_analysis = dict(zip(clinical["aliquot_id"], clinical["analysis_id"]))
    keep_cols = [col for col in tpm.columns if col in rs_to_analysis]
    expr = np.log2(tpm.loc[:, keep_cols] + 1.0)
    expr = expr.rename(columns=rs_to_analysis)
    clinical = clinical.set_index("analysis_id", drop=False).loc[expr.columns]

    risk, present, contrib = score_expression(expr, model)
    sigs, module_present = module_scores(expr, MODULES)

    out = pd.DataFrame(index=expr.columns)
    out["cohort"] = "GSE176307"
    out["analysis_id"] = out.index
    out["sample_id"] = clinical.loc[out.index, "sample_id"].values
    out["patient_id"] = clinical.loc[out.index, "patient_id"].values
    out["aliquot_id"] = clinical.loc[out.index, "aliquot_id"].values
    out["geo_accession"] = clinical.loc[out.index, "geo_accession"].values
    out["response_raw"] = clinical.loc[out.index, "io_response"].astype(str).replace("nan", np.nan).values
    out["response_group"] = out["response_raw"].map({"CR": "Responder", "PR": "Responder", "SD": "Non-responder", "PD": "Non-responder"})
    out["response_binary"] = out["response_group"].map({"Responder": 1, "Non-responder": 0})
    out["nonresponse_binary"] = out["response_group"].map({"Responder": 0, "Non-responder": 1})
    out["os_time_months"] = numeric(clinical.loc[out.index, "overall_survival"]).values / 30.4375
    out["os_event"] = clinical.loc[out.index, "alive"].map({"No": 1, "Yes": 0}).values
    out["pfs_time_months"] = numeric(clinical.loc[out.index, "pfs"]).values / 30.4375
    out["pfs_event"] = clinical.loc[out.index, "progressed"].map({"Yes": 1, "No": 0}).values
    out["age"] = numeric(clinical.loc[out.index, "age"]).values
    out["sex"] = clinical.loc[out.index, "gender"].astype(str).replace("nan", np.nan).values
    out["ecog"] = numeric(clinical.loc[out.index, "ecog"]).values
    out["tmb"] = numeric(clinical.loc[out.index, "tmb"]).values
    out["tmb_group"] = clinical.loc[out.index, "tmb_interpretation"].astype(str).replace("nan", np.nan).values
    out["fgfr3_altered"] = clinical.loc[out.index, "fgfr_mutation_yn"].map({"Y": "FGFR3 altered", "N": "FGFR3 WT"}).values
    out["treatment"] = clinical.loc[out.index, "io_therapy"].astype(str).replace("nan", np.nan).values
    out["immune_phenotype"] = np.nan
    out["tcga_subtype"] = np.nan
    out["risk_score"] = risk.loc[out.index].values
    out["model_genes_present"] = ",".join(present)
    out["n_model_genes_present"] = len(present)
    out["n_model_genes_total"] = len(model)
    out = pd.concat([out, sigs.loc[out.index]], axis=1)

    expr_out = expr.loc[sorted(target_genes & set(expr.index)), out.index].copy()
    expr_out.columns = [f"GSE176307::{sample}" for sample in expr_out.columns]
    out.index = [f"GSE176307::{sample}" for sample in out.index]
    return out.reset_index(drop=True), expr_out


def auc_rank(y: np.ndarray, scores: np.ndarray) -> float:
    valid = (~pd.isna(y)) & (~pd.isna(scores))
    y = np.asarray(y[valid], dtype=float)
    scores = np.asarray(scores[valid], dtype=float)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = pd.Series(scores).rank(method="average").to_numpy()
    rank_sum_pos = ranks[y == 1].sum()
    return float((rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def wilcoxon_approx(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    x = pd.to_numeric(x, errors="coerce").dropna().to_numpy(dtype=float)
    y = pd.to_numeric(y, errors="coerce").dropna().to_numpy(dtype=float)
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan, np.nan
    combined = np.concatenate([x, y])
    ranks = pd.Series(combined).rank(method="average").to_numpy()
    u1 = ranks[:n1].sum() - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    tie_counts = pd.Series(combined).value_counts().to_numpy()
    tie_term = np.sum(tie_counts**3 - tie_counts)
    variance = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else np.nan
    if not np.isfinite(variance) or variance <= 0:
        return float(u1), np.nan
    z = (u1 - n1 * n2 / 2.0) / math.sqrt(variance)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return float(u1), float(p)


def harrell_c_index(time: pd.Series, event: pd.Series, score: pd.Series) -> tuple[float, int]:
    df = pd.DataFrame({"time": time, "event": event, "score": score}).dropna()
    if df.empty:
        return np.nan, 0
    t = df["time"].to_numpy(dtype=float)
    e = df["event"].to_numpy(dtype=float)
    s = df["score"].to_numpy(dtype=float)
    concordant = 0.0
    permissible = 0
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            if t[i] == t[j]:
                continue
            if t[i] < t[j] and e[i] == 1:
                permissible += 1
                concordant += 1.0 if s[i] > s[j] else 0.5 if s[i] == s[j] else 0.0
            elif t[j] < t[i] and e[j] == 1:
                permissible += 1
                concordant += 1.0 if s[j] > s[i] else 0.5 if s[i] == s[j] else 0.0
    if permissible == 0:
        return np.nan, 0
    return float(concordant / permissible), permissible


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, int]:
    df = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(df) < 3:
        return np.nan, len(df)
    rx = df["x"].rank(method="average")
    ry = df["y"].rank(method="average")
    return float(rx.corr(ry)), len(df)


def summarize(combined: pd.DataFrame, model: pd.DataFrame) -> None:
    response_rows = []
    survival_rows = []
    corr_rows = []
    count_rows = []
    availability_rows = []

    for cohort, sub in combined.groupby("cohort", sort=False):
        median_score = sub["risk_score"].median()
        combined.loc[sub.index, "risk_group"] = np.where(sub["risk_score"] >= median_score, "High", "Low")
        combined.loc[sub.index, "risk_group_median_rule"] = "cohort median: High >= median"

    for cohort, sub in combined.groupby("cohort", sort=False):
        evaluable = sub.dropna(subset=["nonresponse_binary", "risk_score"]).copy()
        nonresp = evaluable.loc[evaluable["nonresponse_binary"] == 1, "risk_score"]
        resp = evaluable.loc[evaluable["nonresponse_binary"] == 0, "risk_score"]
        u_stat, p_wilcox = wilcoxon_approx(nonresp, resp)
        auc_nonresponse = auc_rank(evaluable["nonresponse_binary"].to_numpy(), evaluable["risk_score"].to_numpy())
        rates = (
            evaluable.groupby("risk_group")["response_binary"]
            .agg(["count", "sum", "mean"])
            .rename(columns={"count": "n_response_evaluable", "sum": "n_responders", "mean": "responder_rate"})
        )
        response_rows.append(
            {
                "cohort": cohort,
                "response_evaluable_n": len(evaluable),
                "responders": int((evaluable["response_group"] == "Responder").sum()),
                "nonresponders": int((evaluable["response_group"] == "Non-responder").sum()),
                "risk_median_responder": float(resp.median()) if len(resp) else np.nan,
                "risk_median_nonresponder": float(nonresp.median()) if len(nonresp) else np.nan,
                "risk_delta_nonresponder_minus_responder": float(nonresp.median() - resp.median()) if len(resp) and len(nonresp) else np.nan,
                "auc_risk_predicts_nonresponse": auc_nonresponse,
                "wilcoxon_p_nonresponder_vs_responder": p_wilcox,
                "high_risk_responder_rate": rates.loc["High", "responder_rate"] if "High" in rates.index else np.nan,
                "low_risk_responder_rate": rates.loc["Low", "responder_rate"] if "Low" in rates.index else np.nan,
                "high_risk_response_evaluable_n": rates.loc["High", "n_response_evaluable"] if "High" in rates.index else np.nan,
                "low_risk_response_evaluable_n": rates.loc["Low", "n_response_evaluable"] if "Low" in rates.index else np.nan,
            }
        )

        for endpoint, time_col, event_col in [
            ("OS", "os_time_months", "os_event"),
            ("PFS", "pfs_time_months", "pfs_event"),
        ]:
            tmp = sub.dropna(subset=[time_col, event_col, "risk_score"])
            if len(tmp) == 0:
                continue
            c_index, comparable = harrell_c_index(tmp[time_col], tmp[event_col], tmp["risk_score"])
            survival_rows.append(
                {
                    "cohort": cohort,
                    "endpoint": endpoint,
                    "n": len(tmp),
                    "events": int(pd.to_numeric(tmp[event_col], errors="coerce").sum()),
                    "c_index_risk_higher_hazard": c_index,
                    "comparable_pairs": comparable,
                    "median_followup_or_time_months": float(pd.to_numeric(tmp[time_col], errors="coerce").median()),
                }
            )

        for col in [c for c in combined.columns if c.startswith("sig_")]:
            r, n = spearman(sub["risk_score"], sub[col])
            corr_rows.append(
                {
                    "cohort": cohort,
                    "signature": col.replace("sig_", ""),
                    "spearman_r_with_risk": r,
                    "n": n,
                }
            )

        for var in ["response_group", "immune_phenotype", "tcga_subtype", "fgfr3_altered", "tmb_group", "treatment"]:
            if var not in sub.columns:
                continue
            tmp = sub.dropna(subset=[var, "risk_group"])
            if tmp.empty:
                continue
            tab = tmp.groupby(["risk_group", var]).size().reset_index(name="n")
            tab.insert(0, "cohort", cohort)
            tab.insert(1, "variable", var)
            tab = tab.rename(columns={var: "category"})
            count_rows.append(tab)

        present = set(str(sub["model_genes_present"].dropna().iloc[0]).split(",")) if sub["model_genes_present"].notna().any() else set()
        for gene in model["gene"]:
            availability_rows.append(
                {
                    "cohort": cohort,
                    "gene": gene,
                    "coefficient": float(model.loc[model["gene"] == gene, "coefficient"].iloc[0]),
                    "present": gene in present,
                }
            )

    pd.DataFrame(response_rows).to_csv(RESULT_TABLE_DIR / "immunotherapy_response_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(survival_rows).to_csv(RESULT_TABLE_DIR / "immunotherapy_survival_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(corr_rows).to_csv(RESULT_TABLE_DIR / "immunotherapy_signature_correlations.tsv", sep="\t", index=False)
    if count_rows:
        pd.concat(count_rows, ignore_index=True).to_csv(RESULT_TABLE_DIR / "immunotherapy_clinical_association_counts.tsv", sep="\t", index=False)
    pd.DataFrame(availability_rows).to_csv(RESULT_TABLE_DIR / "immunotherapy_model_gene_availability.tsv", sep="\t", index=False)


def main() -> None:
    ensure_dirs()
    model = load_model()
    imvigor, imvigor_expr = read_imvigor(model)
    gse176307, gse_expr = read_gse176307(model)
    combined = pd.concat([imvigor, gse176307], ignore_index=True, sort=False)
    summarize(combined, model)

    score_path = PROCESSED_DIR / "immunotherapy_validation_scores.tsv"
    combined.to_csv(score_path, sep="\t", index=False)
    combined.to_csv(RESULT_TABLE_DIR / "immunotherapy_validation_scores.tsv", sep="\t", index=False)

    expr_combined = pd.concat([imvigor_expr, gse_expr], axis=1, sort=True)
    expr_combined.to_csv(PROCESSED_DIR / "immunotherapy_selected_log_expression.tsv.gz", sep="\t", compression="gzip")
    print(f"wrote {score_path}")
    print(combined.groupby("cohort").size().to_string())
    print(pd.read_csv(RESULT_TABLE_DIR / "immunotherapy_response_metrics.tsv", sep="\t").to_string(index=False))


if __name__ == "__main__":
    main()
