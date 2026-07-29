#!/usr/bin/env python
"""Validate the reduced TCGA-BLCA model in candidate GEO single-cell/spatial datasets.

The script intentionally avoids heavy single-cell dependencies. It extracts only
the five reduced-model genes plus marker genes needed for coarse cell-type
assignment, then computes dataset-wise z-scored reduced-model scores.
"""

from __future__ import annotations

import gzip
import io
import json
import math
import re
import shutil
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LOCAL_PY_LIBS = ROOT / "py_libs"
if LOCAL_PY_LIBS.exists():
    sys.path.append(str(LOCAL_PY_LIBS))


SC_ROOT = ROOT / "data" / "external" / "single_cell"
TABLE_DIR = ROOT / "results" / "tables"
RESULT_DIR = ROOT / "results" / "single_cell"
MODEL_FILE = ROOT / "results" / "model" / "reduced_geo_compatible_model.tsv"

MODEL_GENES = ["EMP1", "AHNAK", "TNFRSF14", "CLEC2D", "GSDMB"]

MARKERS: dict[str, list[str]] = {
    "Epithelial_luminal": ["EPCAM", "KRT8", "KRT18", "KRT19", "UPK1A", "UPK1B", "UPK2", "UPK3A", "KRT20"],
    "Epithelial_basal": ["KRT5", "KRT6A", "KRT14", "KRT17", "TP63", "DSG3", "S100A2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FAP", "COL6A1"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "CLDN5", "ENG", "ESAM", "RAMP2"],
    "T_NK": ["CD3D", "CD3E", "CD2", "TRAC", "NKG7", "GNLY", "PRF1", "KLRD1"],
    "Myeloid": ["LST1", "LYZ", "AIF1", "FCGR3A", "CD68", "MSR1", "C1QA", "C1QB", "C1QC"],
    "B_Plasma": ["MS4A1", "CD79A", "CD79B", "CD74", "MZB1", "JCHAIN", "IGHG1"],
    "Smooth_muscle_pericyte": ["ACTA2", "TAGLN", "MYH11", "DES", "RGS5", "MCAM", "CSPG4"],
    "Proliferating": ["MKI67", "TOP2A", "UBE2C", "STMN1", "HMGB2"],
}

EXTRA_GENES = ["CD274", "PDCD1LG2", "CDH12", "ALDH1A1"]
SELECTED_GENES = sorted(set(MODEL_GENES + EXTRA_GENES + [g for genes in MARKERS.values() for g in genes]))


def decode_array(values) -> list[str]:
    out = []
    for x in values:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8", "replace"))
        else:
            out.append(str(x))
    return out


def zscore(v: np.ndarray) -> np.ndarray:
    v = v.astype(float)
    mu = np.nanmean(v)
    sd = np.nanstd(v)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(v, dtype=float)
    return (v - mu) / sd


def normalize_counts(counts: dict[str, np.ndarray], libsize: np.ndarray) -> pd.DataFrame:
    denom = np.where(libsize > 0, libsize, np.nan)
    data = {}
    for gene in SELECTED_GENES:
        raw = counts.get(gene)
        if raw is None:
            raw = np.zeros_like(libsize, dtype=float)
        data[gene] = np.log1p(raw / denom * 10000.0)
    return pd.DataFrame(data).fillna(0.0)


def infer_sample_condition(dataset: str, sample: str) -> tuple[str, str]:
    s = sample.upper()
    patient = "unknown"
    condition = "Unknown"

    if dataset == "GSE222315":
        m = re.search(r"(P\d+)", s)
        patient = m.group(1).lower() if m else "unknown"
        if "BCA" in s:
            condition = "Tumor"
        elif "NAT" in s:
            condition = "NAT"
    elif dataset == "GSE293189":
        clean = s.replace("-", "_")
        m = re.search(r"(HY|PA|SG)_?([A-Z0-9]+)", clean)
        patient = m.group(0) if m else sample.split("_")[0]
        if "CIS" in clean:
            condition = "CIS"
        elif re.search(r"(^|_)N(\d|_|$)", clean) or re.search(r"\d+N($|_)", clean) or "_NX_" in clean:
            condition = "Normal"
        elif re.search(r"(^|_)T\d", clean) or "_MT" in clean or "_ST" in clean or "_HG_T" in clean or "_FG_T" in clean:
            condition = "Tumor"
        elif "_SN" in clean or "_MN" in clean:
            condition = "Tumor_variant_or_margin"
    elif dataset == "GSE171351":
        condition = "Spatial"
        patient = sample
    return patient, condition


def annotate_cells(dataset: str, meta: pd.DataFrame, log_expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = pd.read_csv(MODEL_FILE, sep="\t")
    coeff = dict(zip(model["gene"], model["coefficient"]))

    model_gene_present = [g for g in MODEL_GENES if g in log_expr.columns]
    score = np.zeros(log_expr.shape[0], dtype=float)
    for gene in MODEL_GENES:
        if gene in log_expr.columns:
            score += coeff[gene] * zscore(log_expr[gene].to_numpy(dtype=float))
    meta = meta.copy()
    meta["reduced_model_score"] = score
    meta["reduced_model_score_z"] = zscore(score)
    meta["risk_group_sc_median"] = np.where(score >= np.nanmedian(score), "High", "Low")

    marker_scores = {}
    for label, genes in MARKERS.items():
        present = [g for g in genes if g in log_expr.columns]
        if present:
            marker_scores[label] = log_expr[present].mean(axis=1).to_numpy(dtype=float)
        else:
            marker_scores[label] = np.zeros(log_expr.shape[0], dtype=float)
    marker_df = pd.DataFrame(marker_scores)
    best = marker_df.idxmax(axis=1)
    best_score = marker_df.max(axis=1)
    meta["inferred_cell_type"] = np.where(best_score >= 0.10, best, "Unknown")
    meta["inferred_cell_type_score"] = best_score
    meta = pd.concat([meta, marker_df.add_prefix("marker_score_")], axis=1)

    long_rows = []
    for gene in MODEL_GENES:
        vals = log_expr[gene].to_numpy(dtype=float) if gene in log_expr.columns else np.zeros(log_expr.shape[0])
        long_rows.append(
            pd.DataFrame(
                {
                    "dataset": dataset,
                    "cell": meta["cell"].to_numpy(),
                    "sample": meta["sample"].to_numpy(),
                    "condition": meta["condition"].to_numpy(),
                    "inferred_cell_type": meta["inferred_cell_type"].to_numpy(),
                    "gene": gene,
                    "log_norm_expr": vals,
                }
            )
        )
    gene_long = pd.concat(long_rows, ignore_index=True)
    return meta, gene_long


def expression_members_from_tar(tf: tarfile.TarFile, expr_name: str) -> tuple[list[str], dict[int, str]]:
    raw = tf.extractfile(expr_name)
    if raw is None:
        raise FileNotFoundError(expr_name)
    gz = gzip.GzipFile(fileobj=raw)
    txt = io.TextIOWrapper(gz, encoding="utf-8", errors="replace")
    header = txt.readline().rstrip("\n").split("\t")
    cells = header[2:]
    idx_to_gene: dict[int, str] = {}
    row_idx = 0
    for line in txt:
        row_idx += 1
        parts = line.split("\t", 2)
        if len(parts) < 2:
            continue
        gene = parts[1].strip()
        if gene in SELECTED_GENES:
            idx_to_gene[row_idx] = gene
    return cells, idx_to_gene


def parse_mtx_from_tar(
    tf: tarfile.TarFile,
    mtx_name: str,
    n_cells: int,
    idx_to_gene: dict[int, str],
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    counts = {gene: np.zeros(n_cells, dtype=float) for gene in SELECTED_GENES}
    libsize = np.zeros(n_cells, dtype=float)
    detected = np.zeros(n_cells, dtype=np.int32)
    raw = tf.extractfile(mtx_name)
    if raw is None:
        raise FileNotFoundError(mtx_name)
    gz = gzip.GzipFile(fileobj=raw)
    dims_seen = False
    for bline in gz:
        if bline.startswith(b"%"):
            continue
        if not dims_seen:
            dims_seen = True
            continue
        parts = bline.split()
        if len(parts) < 3:
            continue
        gi = int(parts[0])
        cj = int(parts[1]) - 1
        val = float(parts[2])
        libsize[cj] += val
        if val > 0:
            detected[cj] += 1
        gene = idx_to_gene.get(gi)
        if gene is not None:
            counts[gene][cj] += val
    return counts, libsize, detected


def sample_name_from_gse222_expr(expr_name: str) -> str:
    base = expr_name.rsplit("_expression.", 1)[0]
    base = re.sub(r"^GSM\d+_", "", base)
    return base


def process_gse222315() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    dataset = "GSE222315"
    tar_path = SC_ROOT / dataset / f"{dataset}_RAW.tar"
    metas = []
    exprs = []
    sample_stats = []
    with tarfile.open(tar_path, "r") as tf:
        names = tf.getnames()
        expr_names = [n for n in names if "_expression." in n and n.endswith(".gz")]
        for expr_name in expr_names:
            sample = sample_name_from_gse222_expr(expr_name)
            mtx_name = expr_name.replace("_expression.txt.gz", "_matrix.mtx.gz").replace("_expression.xls.gz", "_matrix.mtx.gz")
            print(f"GSE222315\t{sample}\treading genes/cells")
            cells, idx_to_gene = expression_members_from_tar(tf, expr_name)
            print(f"GSE222315\t{sample}\tcells={len(cells)}\tselected_indices={len(idx_to_gene)}")
            counts, libsize, detected = parse_mtx_from_tar(tf, mtx_name, len(cells), idx_to_gene)
            log_expr = normalize_counts(counts, libsize)
            patient, condition = infer_sample_condition(dataset, sample)
            keep = (libsize >= 500) & (detected >= 200)
            sample_meta = pd.DataFrame(
                {
                    "dataset": dataset,
                    "cell": [f"{sample}:{c}" for c in cells],
                    "barcode": cells,
                    "sample": sample,
                    "patient": patient,
                    "condition": condition,
                    "n_counts": libsize,
                    "n_genes": detected,
                    "qc_pass": keep,
                    "assay_unit": "single_cell",
                }
            )
            metas.append(sample_meta.loc[keep].reset_index(drop=True))
            exprs.append(log_expr.loc[keep].reset_index(drop=True))
            sample_stats.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "condition": condition,
                    "cells_raw": len(cells),
                    "cells_qc_pass": int(keep.sum()),
                    "selected_genes_detected": ";".join(sorted(set(idx_to_gene.values()))),
                }
            )
    meta = pd.concat(metas, ignore_index=True)
    log_expr = pd.concat(exprs, ignore_index=True)
    meta, gene_long = annotate_cells(dataset, meta, log_expr)
    return meta, gene_long, {"samples": sample_stats}


def sample_name_from_gse293_file(name: str) -> str:
    base = name.rsplit("_dge.txt.gz", 1)[0]
    return re.sub(r"^GSM\d+_", "", base)


def process_gse293189() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    dataset = "GSE293189"
    tar_path = SC_ROOT / dataset / f"{dataset}_RAW.tar"
    metas = []
    exprs = []
    sample_stats = []
    with tarfile.open(tar_path, "r") as tf:
        names = [n for n in tf.getnames() if n.endswith("_dge.txt.gz")]
        for dge_name in names:
            sample = sample_name_from_gse293_file(dge_name)
            raw = tf.extractfile(dge_name)
            if raw is None:
                continue
            gz = gzip.GzipFile(fileobj=raw)
            header = gz.readline().decode("utf-8", "replace").rstrip("\n").split("\t")
            cells = header[1:]
            n_cells = len(cells)
            counts = {gene: np.zeros(n_cells, dtype=float) for gene in SELECTED_GENES}
            libsize = np.zeros(n_cells, dtype=float)
            detected = np.zeros(n_cells, dtype=np.int32)
            selected_found = set()
            print(f"GSE293189\t{sample}\tcells={n_cells}")
            for bline in gz:
                if not bline:
                    continue
                gene_b, sep, rest = bline.partition(b"\t")
                if not sep:
                    continue
                gene = gene_b.decode("utf-8", "replace").strip()
                values = np.fromstring(rest.decode("utf-8", "replace").rstrip("\n"), sep="\t", dtype=float)
                if values.size != n_cells:
                    continue
                libsize += values
                detected += (values > 0)
                if gene in SELECTED_GENES:
                    counts[gene] += values
                    selected_found.add(gene)
            log_expr = normalize_counts(counts, libsize)
            patient, condition = infer_sample_condition(dataset, sample)
            keep = (libsize >= 100) & (detected >= 50)
            sample_meta = pd.DataFrame(
                {
                    "dataset": dataset,
                    "cell": [f"{sample}:{c}" for c in cells],
                    "barcode": cells,
                    "sample": sample,
                    "patient": patient,
                    "condition": condition,
                    "n_counts": libsize,
                    "n_genes": detected,
                    "qc_pass": keep,
                    "assay_unit": "single_cell",
                }
            )
            metas.append(sample_meta.loc[keep].reset_index(drop=True))
            exprs.append(log_expr.loc[keep].reset_index(drop=True))
            sample_stats.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "condition": condition,
                    "cells_raw": len(cells),
                    "cells_qc_pass": int(keep.sum()),
                    "selected_genes_detected": ";".join(sorted(selected_found)),
                }
            )
    meta = pd.concat(metas, ignore_index=True)
    log_expr = pd.concat(exprs, ignore_index=True)
    meta, gene_long = annotate_cells(dataset, meta, log_expr)
    return meta, gene_long, {"samples": sample_stats}


def process_gse171351() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    dataset = "GSE171351"
    import h5py  # noqa: PLC0415

    gz_path = SC_ROOT / dataset / f"{dataset}_combined_visium.h5ad.gz"
    h5ad_path = gz_path.with_suffix("")
    if not h5ad_path.exists():
        with gzip.open(gz_path, "rb") as fh, h5ad_path.open("wb") as out:
            shutil.copyfileobj(fh, out)

    with h5py.File(h5ad_path, "r") as f:
        genes = decode_array(f["var"]["_index"][()])
        gene_to_idx = {g: i for i, g in enumerate(genes)}
        selected_idx = {gene_to_idx[g]: g for g in SELECTED_GENES if g in gene_to_idx}
        obs_index = decode_array(f["obs"]["_index"][()])

        obs = {}
        for key in ["Patient", "sampleID", "array_row", "array_col", "in_tissue", "total_counts", "n_genes_by_counts"]:
            if key not in f["obs"]:
                continue
            arr = f["obs"][key][()]
            if key in ("Patient", "sampleID") and "__categories" in f["obs"] and key in f["obs"]["__categories"]:
                cats = decode_array(f["obs"]["__categories"][key][()])
                obs[key] = [cats[int(x)] for x in arr]
            else:
                obs[key] = arr

        n_spots = len(obs_index)
        counts = {gene: np.zeros(n_spots, dtype=float) for gene in SELECTED_GENES}
        libsize = np.asarray(obs.get("total_counts", np.zeros(n_spots)), dtype=float)
        detected = np.asarray(obs.get("n_genes_by_counts", np.zeros(n_spots)), dtype=float)
        indptr = f["X"]["indptr"][()]
        indices = f["X"]["indices"][()]
        data = f["X"]["data"][()]
        for row in range(n_spots):
            start = indptr[row]
            end = indptr[row + 1]
            cols = indices[start:end]
            vals = data[start:end]
            if libsize[row] == 0:
                libsize[row] = vals.sum()
            for col, val in zip(cols, vals):
                gene = selected_idx.get(int(col))
                if gene is not None:
                    counts[gene][row] += float(val)

        log_expr = normalize_counts(counts, libsize)
        samples = obs.get("sampleID", ["unknown"] * n_spots)
        patients = obs.get("Patient", samples)
        meta = pd.DataFrame(
            {
                "dataset": dataset,
                "cell": obs_index,
                "barcode": obs_index,
                "sample": samples,
                "patient": patients,
                "condition": "Spatial",
                "n_counts": libsize,
                "n_genes": detected,
                "qc_pass": True,
                "assay_unit": "spatial_spot",
                "array_row": obs.get("array_row", np.repeat(np.nan, n_spots)),
                "array_col": obs.get("array_col", np.repeat(np.nan, n_spots)),
                "in_tissue": obs.get("in_tissue", np.repeat(np.nan, n_spots)),
            }
        )
        meta, gene_long = annotate_cells(dataset, meta, log_expr)
        stats = {
            "spots": int(n_spots),
            "selected_genes_detected": ";".join(sorted(set(selected_idx.values()))),
            "note": "Visium spatial transcriptomics spots, not single cells.",
        }
    return meta, gene_long, stats


def summarize(meta_all: pd.DataFrame, gene_long: pd.DataFrame, status_rows: list[dict[str, object]]) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    processed_dir = SC_ROOT / "processed_candidate_sets"
    processed_dir.mkdir(parents=True, exist_ok=True)

    meta_out = processed_dir / "candidate_geo_single_cell_spatial_reduced_model_scores.tsv.gz"
    gene_long_out = processed_dir / "candidate_geo_model_gene_expression_long.tsv.gz"
    meta_all.to_csv(meta_out, sep="\t", index=False, compression="gzip")
    gene_long.to_csv(gene_long_out, sep="\t", index=False, compression="gzip")

    celltype_summary = (
        meta_all.groupby(["dataset", "assay_unit", "inferred_cell_type"], sort=False)
        .agg(
            n_units=("cell", "size"),
            median_score_z=("reduced_model_score_z", "median"),
            mean_score_z=("reduced_model_score_z", "mean"),
            high_fraction=("risk_group_sc_median", lambda x: float((x == "High").mean())),
            mean_n_counts=("n_counts", "mean"),
            mean_n_genes=("n_genes", "mean"),
        )
        .reset_index()
        .sort_values(["dataset", "median_score_z"], ascending=[True, False])
    )
    celltype_summary.to_csv(TABLE_DIR / "candidate_geo_reduced_model_by_inferred_celltype.tsv", sep="\t", index=False)

    sample_summary = (
        meta_all.groupby(["dataset", "assay_unit", "sample", "patient", "condition"], sort=False)
        .agg(
            n_units=("cell", "size"),
            median_score_z=("reduced_model_score_z", "median"),
            mean_score_z=("reduced_model_score_z", "mean"),
            high_fraction=("risk_group_sc_median", lambda x: float((x == "High").mean())),
            epithelial_fraction=("inferred_cell_type", lambda x: float(np.mean([str(v).startswith("Epithelial") for v in x]))),
            stromal_fraction=("inferred_cell_type", lambda x: float(np.mean([str(v) in ["Fibroblast", "Endothelial", "Smooth_muscle_pericyte"] for v in x]))),
            immune_fraction=("inferred_cell_type", lambda x: float(np.mean([str(v) in ["T_NK", "Myeloid", "B_Plasma"] for v in x]))),
        )
        .reset_index()
        .sort_values(["dataset", "median_score_z"], ascending=[True, False])
    )
    sample_summary.to_csv(TABLE_DIR / "candidate_geo_reduced_model_by_sample.tsv", sep="\t", index=False)

    gene_summary = (
        gene_long.groupby(["dataset", "inferred_cell_type", "gene"], sort=False)
        .agg(
            n_units=("cell", "size"),
            mean_log_norm_expr=("log_norm_expr", "mean"),
            median_log_norm_expr=("log_norm_expr", "median"),
            pct_positive=("log_norm_expr", lambda x: float((x > 0).mean() * 100.0)),
        )
        .reset_index()
    )
    gene_summary.to_csv(TABLE_DIR / "candidate_geo_model_gene_by_inferred_celltype.tsv", sep="\t", index=False)

    dataset_summary = (
        meta_all.groupby(["dataset", "assay_unit"], sort=False)
        .agg(
            n_units=("cell", "size"),
            n_samples=("sample", "nunique"),
            median_score_z=("reduced_model_score_z", "median"),
            top_inferred_cell_type=("inferred_cell_type", lambda x: x.value_counts().index[0] if len(x) else ""),
        )
        .reset_index()
    )
    status = pd.DataFrame(status_rows)
    dataset_summary = dataset_summary.merge(status, on="dataset", how="outer")
    dataset_summary.to_csv(TABLE_DIR / "candidate_geo_dataset_status_and_summary.tsv", sep="\t", index=False)

    manifest = {
        "cell_or_spot_scores": str(meta_out),
        "model_gene_expression_long": str(gene_long_out),
        "celltype_summary": str(TABLE_DIR / "candidate_geo_reduced_model_by_inferred_celltype.tsv"),
        "sample_summary": str(TABLE_DIR / "candidate_geo_reduced_model_by_sample.tsv"),
        "gene_summary": str(TABLE_DIR / "candidate_geo_model_gene_by_inferred_celltype.tsv"),
        "dataset_status": str(TABLE_DIR / "candidate_geo_dataset_status_and_summary.tsv"),
    }
    (RESULT_DIR / "candidate_geo_single_cell_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    metas = []
    gene_longs = []
    status_rows: list[dict[str, object]] = []

    for dataset, func in [
        ("GSE222315", process_gse222315),
        ("GSE293189", process_gse293189),
        ("GSE171351", process_gse171351),
    ]:
        try:
            meta, gene_long, stats = func()
            metas.append(meta)
            gene_longs.append(gene_long)
            status_rows.append(
                {
                    "dataset": dataset,
                    "analysis_status": "analyzed",
                    "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset}",
                    "notes": json.dumps(stats, ensure_ascii=False),
                }
            )
            meta.to_csv(SC_ROOT / dataset / f"{dataset}_reduced_model_scores.tsv.gz", sep="\t", index=False, compression="gzip")
        except Exception as exc:  # noqa: BLE001
            status_rows.append(
                {
                    "dataset": dataset,
                    "analysis_status": "failed",
                    "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset}",
                    "notes": repr(exc),
                }
            )
            print(f"FAILED\t{dataset}\t{exc}")

    status_rows.append(
        {
            "dataset": "GSE269877",
            "analysis_status": "deferred_large_rds",
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269877",
            "notes": "Supplementary files are GSE269877_dta_cancer.submission.rds.gz (13G) and GSE269877_dta_normal.submission.rds.gz (2.4G); not downloaded for this lightweight test.",
        }
    )

    if metas:
        summarize(pd.concat(metas, ignore_index=True), pd.concat(gene_longs, ignore_index=True), status_rows)
    else:
        pd.DataFrame(status_rows).to_csv(TABLE_DIR / "candidate_geo_dataset_status_and_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
