#!/usr/bin/env python
"""TCGA risk-group pathway enrichment and immune marker signatures."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
GENESETS = ROOT / "data" / "external" / "gene_sets"
TABLES = ROOT / "results" / "tables"

IMMUNE_SIGNATURES = {
    "CD8_T": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "CCL5"],
    "Cytotoxic": ["GZMB", "PRF1", "GNLY", "NKG7", "GZMA"],
    "NK_cell": ["NKG7", "KLRD1", "GNLY", "NCR1", "KLRF1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT"],
    "Macrophage": ["CD68", "CD163", "MSR1", "MRC1", "CSF1R"],
    "M1_like": ["IL1B", "CXCL9", "CXCL10", "CD86", "NOS2"],
    "M2_like": ["CD163", "MRC1", "MSR1", "CCL18", "IL10"],
    "Dendritic_cell": ["ITGAX", "HLA-DRA", "CD1C", "FCER1A", "CLEC10A"],
    "Neutrophil": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2"],
    "CAF": ["ACTA2", "COL1A1", "COL1A2", "FAP", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "CDH5", "ENG"],
    "Checkpoint": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
}


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    out = np.full_like(p, np.nan)
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


def normal_p(z: np.ndarray) -> np.ndarray:
    return np.array([math.erfc(abs(float(v)) / math.sqrt(2.0)) for v in z])


def expression_screen(expr: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    high_cols = groups[groups.eq("High")].index.tolist()
    low_cols = groups[groups.eq("Low")].index.tolist()
    high = expr[high_cols].to_numpy(dtype=float)
    low = expr[low_cols].to_numpy(dtype=float)
    mean_high = np.nanmean(high, axis=1)
    mean_low = np.nanmean(low, axis=1)
    var_high = np.nanvar(high, axis=1, ddof=1)
    var_low = np.nanvar(low, axis=1, ddof=1)
    se = np.sqrt(var_high / high.shape[1] + var_low / low.shape[1])
    z = (mean_high - mean_low) / np.maximum(se, 1e-8)
    p = normal_p(z)
    out = pd.DataFrame(
        {
            "gene": expr.index,
            "high_mean": mean_high,
            "low_mean": mean_low,
            "diff_high_minus_low": mean_high - mean_low,
            "z": z,
            "pvalue": p,
        }
    )
    out["FDR"] = bh_fdr(out["pvalue"].to_numpy())
    return out.sort_values(["FDR", "pvalue"])


def read_geneset_file(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term = parts[0]
        genes = {g.strip() for g in parts[2:] if g.strip()}
        if genes:
            out[term] = genes
    return out


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k_minus_1: int, N: int, K: int, n: int) -> float:
    max_k = min(K, n)
    vals = []
    for k in range(k_minus_1 + 1, max_k + 1):
        vals.append(log_choose(K, k) + log_choose(N - K, n - k) - log_choose(N, n))
    if not vals:
        return 1.0
    m = max(vals)
    return float(min(1.0, math.exp(m) * sum(math.exp(v - m) for v in vals)))


def enrich(gene_list: set[str], universe: set[str], genesets: dict[str, set[str]], direction: str, library: str) -> pd.DataFrame:
    rows = []
    gl = gene_list & universe
    N = len(universe)
    n = len(gl)
    for term, genes in genesets.items():
        gs = genes & universe
        if len(gs) < 5:
            continue
        overlap = gl & gs
        if len(overlap) < 2:
            continue
        p = hypergeom_sf(len(overlap) - 1, N, len(gs), n)
        rows.append(
            {
                "library": library,
                "direction": direction,
                "term": term,
                "overlap_n": len(overlap),
                "query_n": n,
                "term_n": len(gs),
                "pvalue": p,
                "overlap_genes": ";".join(sorted(overlap)),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["FDR"] = bh_fdr(df["pvalue"].to_numpy())
    return df.sort_values(["FDR", "pvalue"])


def immune_scores(expr: pd.DataFrame, groups: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_scores = pd.DataFrame(index=expr.columns)
    used_rows = []
    z_expr = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    for name, genes in IMMUNE_SIGNATURES.items():
        present = [g for g in genes if g in z_expr.index]
        if not present:
            continue
        sample_scores[name] = z_expr.loc[present].mean(axis=0)
        used_rows.append({"signature": name, "genes_used": ";".join(present), "n_genes_used": len(present)})
    sample_scores["risk_group"] = groups.reindex(sample_scores.index)
    rows = []
    for sig in [c for c in sample_scores.columns if c != "risk_group"]:
        high = sample_scores.loc[sample_scores["risk_group"].eq("High"), sig].dropna()
        low = sample_scores.loc[sample_scores["risk_group"].eq("Low"), sig].dropna()
        diff = high.mean() - low.mean()
        se = math.sqrt(high.var(ddof=1) / len(high) + low.var(ddof=1) / len(low))
        z = diff / max(se, 1e-8)
        rows.append(
            {
                "signature": sig,
                "high_mean": high.mean(),
                "low_mean": low.mean(),
                "diff_high_minus_low": diff,
                "z": z,
                "pvalue": math.erfc(abs(z) / math.sqrt(2.0)),
            }
        )
    stats = pd.DataFrame(rows)
    stats["FDR"] = bh_fdr(stats["pvalue"].to_numpy())
    stats = stats.merge(pd.DataFrame(used_rows), on="signature", how="left")
    return sample_scores, stats.sort_values(["FDR", "pvalue"])


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    expr = pd.read_csv(PROCESSED / "blca_log2_tpm_gene_symbol_tumor.tsv.gz", sep="\t", index_col=0)
    scores = pd.read_csv(TABLES / "tcga_lasso_risk_scores.tsv", sep="\t")
    groups = scores.set_index("matrix_column")["lasso_risk_group"]
    expr = expr.loc[:, [c for c in groups.index if c in expr.columns]]
    groups = groups.reindex(expr.columns)

    screen = expression_screen(expr, groups)
    screen.to_csv(TABLES / "tcga_fixed_risk_high_vs_low_expression_screen.tsv", sep="\t", index=False)
    up = set(screen.loc[(screen["FDR"] < 0.05) & (screen["diff_high_minus_low"] > 0.5), "gene"])
    down = set(screen.loc[(screen["FDR"] < 0.05) & (screen["diff_high_minus_low"] < -0.5), "gene"])
    universe = set(expr.index)

    enrichments = []
    for fname, lib in [("MSigDB_Hallmark_2020.txt", "Hallmark_2020"), ("KEGG_2021_Human.txt", "KEGG_2021_Human")]:
        genesets = read_geneset_file(GENESETS / fname)
        enrichments.append(enrich(up, universe, genesets, "High_risk_up", lib))
        enrichments.append(enrich(down, universe, genesets, "High_risk_down", lib))
    enrich_df = pd.concat([df for df in enrichments if not df.empty], ignore_index=True)
    enrich_df.to_csv(TABLES / "tcga_fixed_risk_pathway_enrichment.tsv", sep="\t", index=False)

    sample_scores, immune_stats = immune_scores(expr, groups)
    sample_scores.to_csv(TABLES / "tcga_fixed_risk_immune_signature_scores.tsv", sep="\t")
    immune_stats.to_csv(TABLES / "tcga_fixed_risk_immune_signature_stats.tsv", sep="\t", index=False)

    print(f"risk DE up genes: {len(up)}")
    print(f"risk DE down genes: {len(down)}")
    print("top pathway terms:")
    print(enrich_df.head(10)[["library", "direction", "term", "FDR"]].to_string(index=False))
    print("immune signatures:")
    print(immune_stats[["signature", "diff_high_minus_low", "FDR"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
