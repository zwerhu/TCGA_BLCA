#!/usr/bin/env python
"""Mutation/CNV/RPPA support analyses for TCGA fixed risk groups."""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"

NONSILENT_EXCLUDE = {"Silent", "Intron", "IGR", "3'UTR", "5'UTR", "RNA", "Targeted_Region", "Flank"}


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(vals)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty_like(vals)
    tmp[order] = q
    out[ok] = tmp
    return out


def p_from_z(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def compare_continuous(df: pd.DataFrame, value_col: str, group_col: str = "lasso_risk_group") -> dict[str, float]:
    high = pd.to_numeric(df.loc[df[group_col].eq("High"), value_col], errors="coerce").dropna()
    low = pd.to_numeric(df.loc[df[group_col].eq("Low"), value_col], errors="coerce").dropna()
    diff = high.mean() - low.mean()
    se = math.sqrt(high.var(ddof=1) / len(high) + low.var(ddof=1) / len(low)) if len(high) > 1 and len(low) > 1 else float("nan")
    z = diff / max(se, 1e-12) if math.isfinite(se) else float("nan")
    return {
        "high_mean": high.mean(),
        "low_mean": low.mean(),
        "diff_high_minus_low": diff,
        "z": z,
        "pvalue": p_from_z(z) if math.isfinite(z) else float("nan"),
        "n_high": len(high),
        "n_low": len(low),
    }


def proportion_z(a: int, n1: int, b: int, n2: int) -> float:
    p1 = a / n1 if n1 else 0
    p2 = b / n2 if n2 else 0
    p = (a + b) / (n1 + n2)
    se = math.sqrt(max(p * (1 - p) * (1 / n1 + 1 / n2), 1e-12))
    return (p1 - p2) / se


def load_risk() -> pd.DataFrame:
    risk = pd.read_csv(TABLES / "tcga_lasso_risk_scores.tsv", sep="\t")
    return risk[["patient_barcode", "sample_barcode", "matrix_column", "lasso_risk_score", "lasso_risk_group"]].drop_duplicates("patient_barcode")


def mutation_analysis(risk: pd.DataFrame) -> None:
    counts = pd.read_csv(TABLES / "mutation_counts_by_patient.tsv", sep="\t")
    mut = counts.merge(risk, on="patient_barcode", how="inner")
    burden_stats = compare_continuous(mut, "nonsilent_mutations")
    pd.DataFrame([{**{"feature": "nonsilent_mutations"}, **burden_stats}]).to_csv(
        TABLES / "tcga_fixed_risk_mutation_burden_stats.tsv", sep="\t", index=False
    )
    mut.to_csv(TABLES / "tcga_fixed_risk_mutation_burden_by_patient.tsv", sep="\t", index=False)

    maf = pd.read_csv(PROCESSED / "blca_somatic_mutations.maf.gz", sep="\t", low_memory=False)
    maf = maf[~maf["Variant_Classification"].isin(NONSILENT_EXCLUDE)].copy()
    maf["Patient_Barcode"] = maf["Patient_Barcode"].astype(str).str.slice(0, 12)
    risk_map = risk.set_index("patient_barcode")["lasso_risk_group"].to_dict()
    maf["risk_group"] = maf["Patient_Barcode"].map(risk_map)
    maf = maf[maf["risk_group"].notna()]
    top_genes = maf.drop_duplicates(["Patient_Barcode", "Hugo_Symbol"])["Hugo_Symbol"].value_counts().head(25).index.tolist()
    rows = []
    for gene in top_genes:
        pats = set(maf.loc[maf["Hugo_Symbol"].eq(gene), "Patient_Barcode"])
        high_pats = set(risk.loc[risk["lasso_risk_group"].eq("High"), "patient_barcode"])
        low_pats = set(risk.loc[risk["lasso_risk_group"].eq("Low"), "patient_barcode"])
        a = len(pats & high_pats)
        b = len(pats & low_pats)
        z = proportion_z(a, len(high_pats), b, len(low_pats))
        rows.append(
            {
                "gene": gene,
                "high_mutated": a,
                "high_total": len(high_pats),
                "low_mutated": b,
                "low_total": len(low_pats),
                "high_frequency": a / len(high_pats),
                "low_frequency": b / len(low_pats),
                "frequency_diff_high_minus_low": a / len(high_pats) - b / len(low_pats),
                "z": z,
                "pvalue": p_from_z(z),
            }
        )
    out = pd.DataFrame(rows)
    out["FDR"] = bh_fdr(out["pvalue"].to_numpy())
    out.to_csv(TABLES / "tcga_fixed_risk_top_mutated_gene_frequencies.tsv", sep="\t", index=False)


def cnv_analysis(risk: pd.DataFrame) -> None:
    cnv = pd.read_csv(PROCESSED / "blca_masked_cnv_segments.tsv.gz", sep="\t")
    cnv["segment_length"] = pd.to_numeric(cnv["End"], errors="coerce") - pd.to_numeric(cnv["Start"], errors="coerce") + 1
    cnv["abs_segment_mean"] = pd.to_numeric(cnv["Segment_Mean"], errors="coerce").abs()
    cnv["altered_length_abs_0_2"] = np.where(cnv["abs_segment_mean"] >= 0.2, cnv["segment_length"], 0)
    summary = (
        cnv.groupby("patient_barcode")
        .agg(
            cnv_segment_count=("Segment_Mean", "size"),
            mean_abs_segment_mean=("abs_segment_mean", "mean"),
            total_altered_length_abs_0_2=("altered_length_abs_0_2", "sum"),
        )
        .reset_index()
    )
    summary = summary.merge(risk, on="patient_barcode", how="inner")
    summary.to_csv(TABLES / "tcga_fixed_risk_cnv_burden_by_patient.tsv", sep="\t", index=False)
    rows = []
    for col in ["cnv_segment_count", "mean_abs_segment_mean", "total_altered_length_abs_0_2"]:
        rows.append({**{"feature": col}, **compare_continuous(summary, col)})
    pd.DataFrame(rows).to_csv(TABLES / "tcga_fixed_risk_cnv_burden_stats.tsv", sep="\t", index=False)


def rppa_analysis(risk: pd.DataFrame) -> None:
    rppa = pd.read_csv(PROCESSED / "blca_rppa_matrix.tsv.gz", sep="\t", index_col=0)
    smap = pd.read_csv(PROCESSED / "blca_rppa_sample_map.tsv", sep="\t")
    col_to_patient = smap.set_index("matrix_column")["patient_barcode"].to_dict()
    risk_map = risk.set_index("patient_barcode")["lasso_risk_group"].to_dict()
    cols = [c for c in rppa.columns if col_to_patient.get(c) in risk_map]
    rppa = rppa[cols].apply(pd.to_numeric, errors="coerce")
    groups = pd.Series({c: risk_map[col_to_patient[c]] for c in cols})
    rows = []
    for agid in rppa.index:
        vals = pd.DataFrame({"value": rppa.loc[agid, cols].to_numpy(dtype=float), "lasso_risk_group": groups.loc[cols].to_numpy()})
        rows.append({**{"protein_id": agid}, **compare_continuous(vals, "value")})
    out = pd.DataFrame(rows)
    out["FDR"] = bh_fdr(out["pvalue"].to_numpy())
    out = out.sort_values(["FDR", "pvalue"])
    out.to_csv(TABLES / "tcga_fixed_risk_rppa_differential.tsv", sep="\t", index=False)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    risk = load_risk()
    mutation_analysis(risk)
    cnv_analysis(risk)
    rppa_analysis(risk)
    print("mutation burden:")
    print(pd.read_csv(TABLES / "tcga_fixed_risk_mutation_burden_stats.tsv", sep="\t").to_string(index=False))
    print("CNV burden:")
    print(pd.read_csv(TABLES / "tcga_fixed_risk_cnv_burden_stats.tsv", sep="\t").to_string(index=False))
    print("Top RPPA:")
    print(pd.read_csv(TABLES / "tcga_fixed_risk_rppa_differential.tsv", sep="\t").head(10).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
