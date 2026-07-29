#!/usr/bin/env python
"""
TCGA-BLCA open-data pipeline using the GDC API.

The script intentionally uses only Python standard-library networking for
downloads. pandas/numpy are used for matrix building because they make the
large expression tables much less error-prone.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

GDC_API = "https://api.gdc.cancer.gov"
PROJECT = "TCGA-BLCA"

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
META_DIR = DATA_DIR / "metadata"
PROCESSED_DIR = DATA_DIR / "processed"
RESULT_TABLE_DIR = ROOT / "results" / "tables"
LOG_DIR = ROOT / "logs"

DATASETS: dict[str, dict[str, Any]] = {
    "rna_star_counts": {
        "raw_subdir": "rna_star_counts",
        "filters": [
            ("data_type", ["Gene Expression Quantification"]),
            ("analysis.workflow_type", ["STAR - Counts"]),
            ("data_format", ["TSV"]),
        ],
    },
    "somatic_maf": {
        "raw_subdir": "somatic_maf",
        "filters": [
            ("data_type", ["Masked Somatic Mutation"]),
            ("data_format", ["MAF"]),
        ],
    },
    "masked_cnv_segments": {
        "raw_subdir": "masked_cnv_segments",
        "filters": [
            ("data_type", ["Masked Copy Number Segment"]),
        ],
    },
    "rppa_protein": {
        "raw_subdir": "rppa_protein",
        "filters": [
            ("data_type", ["Protein Expression Quantification"]),
            ("data_format", ["TSV"]),
        ],
    },
    "methylation_beta": {
        "raw_subdir": "methylation_beta",
        "filters": [
            ("data_type", ["Methylation Beta Value"]),
        ],
    },
}

CORE_DATASETS = ["rna_star_counts", "somatic_maf", "masked_cnv_segments", "rppa_protein"]


def ensure_dirs() -> None:
    for path in [RAW_DIR, META_DIR, PROCESSED_DIR, RESULT_TABLE_DIR, LOG_DIR]:
        path.mkdir(parents=True, exist_ok=True)
    for spec in DATASETS.values():
        (RAW_DIR / spec["raw_subdir"]).mkdir(parents=True, exist_ok=True)


def log(message: str) -> None:
    print(message, flush=True)


def gdc_request(endpoint: str, params: dict[str, Any] | None = None, timeout: int = 120) -> Any:
    query = urllib.parse.urlencode(params or {})
    url = f"{GDC_API}/{endpoint.lstrip('/')}"
    if query:
        url += "?" + query
    last_error: Exception | None = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            wait = min(60, 2**attempt)
            log(f"[WARN] GDC request failed on attempt {attempt}: {exc}; retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"GDC request failed after retries: {url}") from last_error


def filter_content(project: str, dataset: str | None = None) -> dict[str, Any]:
    content: list[dict[str, Any]] = [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": [project]}},
        {"op": "=", "content": {"field": "access", "value": ["open"]}},
    ]
    if dataset:
        for field, value in DATASETS[dataset]["filters"]:
            content.append({"op": "=", "content": {"field": field, "value": value}})
    return {"op": "and", "content": content}


def fetch_files(project: str, dataset: str) -> list[dict[str, Any]]:
    fields = ",".join(
        [
            "file_id",
            "file_name",
            "file_size",
            "md5sum",
            "data_category",
            "data_type",
            "data_format",
            "access",
            "analysis.workflow_type",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.sample_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
            "cases.samples.sample_type_id",
            "cases.samples.tissue_type",
        ]
    )
    params = {
        "filters": json.dumps(filter_content(project, dataset)),
        "fields": fields,
        "expand": "cases.samples",
        "format": "JSON",
        "size": "5000",
        "sort": "file_name:asc",
    }
    data = gdc_request("files", params)
    hits = data["data"]["hits"]
    if len(hits) >= 5000:
        raise RuntimeError("Increase pagination support: more than 5000 files returned.")
    return hits


def first_case(hit: dict[str, Any]) -> dict[str, Any]:
    cases = hit.get("cases") or []
    return cases[0] if cases else {}


def first_sample(hit: dict[str, Any]) -> dict[str, Any]:
    case = first_case(hit)
    samples = case.get("samples") or []
    return samples[0] if samples else {}


def flatten_file_hit(hit: dict[str, Any], dataset: str) -> dict[str, Any]:
    case = first_case(hit)
    sample = first_sample(hit)
    analysis = hit.get("analysis") or {}
    return {
        "dataset": dataset,
        "file_id": hit.get("file_id") or hit.get("id"),
        "file_name": hit.get("file_name"),
        "file_size": hit.get("file_size"),
        "md5sum": hit.get("md5sum"),
        "data_category": hit.get("data_category"),
        "data_type": hit.get("data_type"),
        "data_format": hit.get("data_format"),
        "workflow_type": analysis.get("workflow_type"),
        "case_id": case.get("case_id"),
        "patient_barcode": case.get("submitter_id"),
        "sample_id": sample.get("sample_id"),
        "sample_barcode": sample.get("submitter_id"),
        "sample_type": sample.get("sample_type"),
        "sample_type_id": sample.get("sample_type_id"),
        "tissue_type": sample.get("tissue_type"),
    }


def write_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def query_files(project: str, datasets: list[str]) -> None:
    ensure_dirs()
    summary_rows = []
    for dataset in datasets:
        log(f"[QUERY] {project} {dataset}")
        hits = fetch_files(project, dataset)
        rows = [flatten_file_hit(hit, dataset) for hit in hits]
        json_path = META_DIR / f"{dataset}_files.json"
        tsv_path = META_DIR / f"{dataset}_metadata.tsv"
        manifest_path = META_DIR / f"{dataset}_gdc_manifest.tsv"
        json_path.write_text(json.dumps(hits, indent=2), encoding="utf-8")
        write_tsv(rows, tsv_path)
        manifest_rows = [{"id": row["file_id"], "filename": row["file_name"]} for row in rows]
        write_tsv(manifest_rows, manifest_path)
        total_size = sum(int(row["file_size"] or 0) for row in rows)
        summary_rows.append(
            {
                "dataset": dataset,
                "files": len(rows),
                "total_size_gb": f"{total_size / 1024**3:.3f}",
                "metadata_tsv": str(tsv_path),
                "raw_subdir": str(RAW_DIR / DATASETS[dataset]["raw_subdir"]),
            }
        )
    write_tsv(summary_rows, META_DIR / "dataset_summary.tsv")
    log(f"[DONE] Wrote metadata to {META_DIR}")


def download_one(file_id: str, file_name: str, out_path: Path, expected_size: int | None = None) -> tuple[str, bool, str]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and expected_size and out_path.stat().st_size == expected_size:
        return file_id, False, "exists"
    if out_path.exists() and expected_size is None:
        return file_id, False, "exists"

    tmp = out_path.with_suffix(out_path.suffix + ".part")
    url = f"{GDC_API}/data/{file_id}"
    last_error: Exception | None = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(url, timeout=240) as response, tmp.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
            tmp.replace(out_path)
            return file_id, True, "downloaded"
        except Exception as exc:  # noqa: BLE001 - network retries need broad catch
            last_error = exc
            if tmp.exists():
                tmp.unlink(missing_ok=True)
            wait = min(90, 2**attempt)
            time.sleep(wait)
    return file_id, False, f"failed: {last_error}"


def download_dataset(dataset: str, threads: int = 4) -> None:
    metadata_path = META_DIR / f"{dataset}_metadata.tsv"
    rows = read_tsv(metadata_path)
    if not rows:
        raise FileNotFoundError(f"Missing metadata: {metadata_path}. Run query first.")
    raw_subdir = RAW_DIR / DATASETS[dataset]["raw_subdir"]
    jobs = []
    for row in rows:
        file_id = row["file_id"]
        file_name = row["file_name"]
        expected_size = int(row["file_size"]) if row.get("file_size") else None
        out_path = raw_subdir / file_id / file_name
        jobs.append((file_id, file_name, out_path, expected_size))

    log(f"[DOWNLOAD] {dataset}: {len(jobs)} files -> {raw_subdir}")
    status_rows = []
    downloaded = skipped = failed = 0
    with ThreadPoolExecutor(max_workers=max(1, threads)) as executor:
        futures = [executor.submit(download_one, *job) for job in jobs]
        for index, future in enumerate(as_completed(futures), start=1):
            file_id, did_download, status = future.result()
            if status.startswith("failed"):
                failed += 1
            elif did_download:
                downloaded += 1
            else:
                skipped += 1
            status_rows.append({"dataset": dataset, "file_id": file_id, "status": status})
            if index % 25 == 0 or index == len(jobs):
                log(
                    f"  {index}/{len(jobs)} complete "
                    f"(downloaded={downloaded}, skipped={skipped}, failed={failed})"
                )
    write_tsv(status_rows, LOG_DIR / f"{dataset}_download_status.tsv")
    if failed:
        raise RuntimeError(f"{dataset}: {failed} downloads failed. See logs.")


def download_files(datasets: list[str], threads: int = 4) -> None:
    ensure_dirs()
    for dataset in datasets:
        download_dataset(dataset, threads=threads)


def as_int(value: Any) -> int | None:
    if value in [None, "", "not reported", "Not Reported", "unknown", "Unknown"]:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def clean_value(value: Any) -> Any:
    if value in [None, "", "not reported", "Not Reported", "unknown", "Unknown", "--"]:
        return None
    return value


def primary_diagnosis(case: dict[str, Any]) -> dict[str, Any]:
    diagnoses = case.get("diagnoses") or []
    primaries = [d for d in diagnoses if d.get("diagnosis_is_primary_disease") is True]
    if primaries:
        return sorted(primaries, key=lambda x: x.get("days_to_diagnosis") or 0)[0]
    if diagnoses:
        return sorted(diagnoses, key=lambda x: x.get("days_to_diagnosis") or 0)[0]
    return {}


def max_follow_up_days(case: dict[str, Any]) -> int | None:
    values = []
    for follow_up in case.get("follow_ups") or []:
        day = as_int(follow_up.get("days_to_follow_up"))
        if day is not None:
            values.append(day)
    diagnosis_day = as_int(primary_diagnosis(case).get("days_to_last_follow_up"))
    if diagnosis_day is not None:
        values.append(diagnosis_day)
    return max(values) if values else None


def query_clinical(project: str) -> None:
    ensure_dirs()
    filters = {"op": "=", "content": {"field": "project.project_id", "value": [project]}}
    params = {
        "filters": json.dumps(filters),
        "expand": "diagnoses,demographic,exposures,follow_ups",
        "format": "JSON",
        "size": "2000",
    }
    data = gdc_request("cases", params)
    cases = data["data"]["hits"]
    (META_DIR / "clinical_cases.json").write_text(json.dumps(cases, indent=2), encoding="utf-8")

    rows = []
    for case in cases:
        demographic = case.get("demographic") or {}
        diagnosis = primary_diagnosis(case)
        days_to_death = as_int(demographic.get("days_to_death"))
        days_to_last_follow_up = max_follow_up_days(case)
        vital_status = clean_value(demographic.get("vital_status"))
        os_event = 1 if str(vital_status).lower() == "dead" else 0
        os_time = days_to_death if os_event and days_to_death is not None else days_to_last_follow_up
        exposure = next((e for e in case.get("exposures") or [] if e.get("exposure_type") == "Tobacco"), {})
        rows.append(
            {
                "patient_barcode": case.get("submitter_id"),
                "case_id": case.get("case_id"),
                "vital_status": vital_status,
                "os_event": os_event,
                "os_time_days": os_time,
                "gender": clean_value(demographic.get("gender") or demographic.get("sex_at_birth")),
                "race": clean_value(demographic.get("race")),
                "ethnicity": clean_value(demographic.get("ethnicity")),
                "age_at_index": clean_value(demographic.get("age_at_index")),
                "age_at_diagnosis_years": (
                    round(as_int(diagnosis.get("age_at_diagnosis")) / 365.25, 2)
                    if as_int(diagnosis.get("age_at_diagnosis")) is not None
                    else None
                ),
                "year_of_diagnosis": clean_value(diagnosis.get("year_of_diagnosis")),
                "primary_diagnosis": clean_value(diagnosis.get("primary_diagnosis")),
                "tumor_grade": clean_value(diagnosis.get("tumor_grade")),
                "ajcc_pathologic_stage": clean_value(diagnosis.get("ajcc_pathologic_stage")),
                "ajcc_pathologic_t": clean_value(diagnosis.get("ajcc_pathologic_t")),
                "ajcc_pathologic_n": clean_value(diagnosis.get("ajcc_pathologic_n")),
                "ajcc_pathologic_m": clean_value(diagnosis.get("ajcc_pathologic_m")),
                "tobacco_smoking_status": clean_value(exposure.get("tobacco_smoking_status")),
                "pack_years_smoked": clean_value(exposure.get("pack_years_smoked")),
            }
        )
    write_tsv(rows, PROCESSED_DIR / "blca_clinical_survival.tsv")
    log(f"[DONE] clinical cases={len(rows)} -> {PROCESSED_DIR / 'blca_clinical_survival.tsv'}")


def require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:
        raise RuntimeError("pandas is required for matrix building.") from exc
    return pd


def metadata_by_file_id(dataset: str) -> dict[str, dict[str, str]]:
    return {row["file_id"]: row for row in read_tsv(META_DIR / f"{dataset}_metadata.tsv")}


def unique_column_names(rows: list[dict[str, str]]) -> dict[str, str]:
    base_names = []
    for row in rows:
        sample = row.get("sample_barcode") or row.get("patient_barcode") or row["file_id"]
        base_names.append(sample)
    counts = Counter(base_names)
    seen = Counter()
    mapping = {}
    for row, base in zip(rows, base_names):
        if counts[base] == 1:
            name = base
        else:
            seen[base] += 1
            name = f"{base}__{seen[base]}_{row['file_id'][:8]}"
        mapping[row["file_id"]] = name
    return mapping


def build_expression_matrices() -> None:
    pd = require_pandas()
    rows = read_tsv(META_DIR / "rna_star_counts_metadata.tsv")
    if not rows:
        log("[SKIP] No RNA metadata found.")
        return
    file_to_col = unique_column_names(rows)
    count_series = []
    tpm_series = []
    fpkm_uq_series = []
    annotations = None
    sample_map = []
    for idx, row in enumerate(rows, start=1):
        file_id = row["file_id"]
        path = RAW_DIR / "rna_star_counts" / file_id / row["file_name"]
        if not path.exists():
            log(f"[WARN] Missing RNA file, skip: {path}")
            continue
        df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
        if "gene_id" not in df.columns:
            log(f"[WARN] Unexpected RNA format, skip: {path}")
            continue
        df = df[df["gene_id"].astype(str).str.startswith("ENSG")].copy()
        df["gene_id_clean"] = df["gene_id"].astype(str).str.replace(r"\..*$", "", regex=True)
        if annotations is None:
            annotations = df[["gene_id", "gene_id_clean", "gene_name", "gene_type"]].drop_duplicates("gene_id")
        col = file_to_col[file_id]
        count_series.append(df.set_index("gene_id")["unstranded"].rename(col))
        if "tpm_unstranded" in df.columns:
            tpm_series.append(df.set_index("gene_id")["tpm_unstranded"].rename(col))
        if "fpkm_uq_unstranded" in df.columns:
            fpkm_uq_series.append(df.set_index("gene_id")["fpkm_uq_unstranded"].rename(col))
        sample_map.append({**row, "matrix_column": col})
        if idx % 50 == 0 or idx == len(rows):
            log(f"  parsed RNA {idx}/{len(rows)}")

    if annotations is not None:
        annotations.to_csv(PROCESSED_DIR / "blca_gene_annotation.tsv", sep="\t", index=False)
    if count_series:
        pd.concat(count_series, axis=1).to_csv(PROCESSED_DIR / "blca_rna_counts_matrix.tsv.gz", sep="\t", compression="gzip")
    if tpm_series:
        pd.concat(tpm_series, axis=1).to_csv(PROCESSED_DIR / "blca_rna_tpm_matrix.tsv.gz", sep="\t", compression="gzip")
    if fpkm_uq_series:
        pd.concat(fpkm_uq_series, axis=1).to_csv(
            PROCESSED_DIR / "blca_rna_fpkm_uq_matrix.tsv.gz", sep="\t", compression="gzip"
        )
    write_tsv(sample_map, PROCESSED_DIR / "blca_rna_sample_map.tsv")
    log("[DONE] RNA matrices built.")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def gzip_writer(path: Path):
    return gzip.open(path, "wt", encoding="utf-8", newline="")


def build_maf() -> None:
    rows = read_tsv(META_DIR / "somatic_maf_metadata.tsv")
    if not rows:
        log("[SKIP] No MAF metadata found.")
        return
    out_maf = PROCESSED_DIR / "blca_somatic_mutations.maf.gz"
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    header: list[str] | None = None
    nonsilent_exclude = {
        "Silent",
        "Intron",
        "IGR",
        "3'UTR",
        "5'UTR",
        "RNA",
        "Targeted_Region",
        "Flank",
    }
    with gzip_writer(out_maf) as out:
        writer = None
        for idx, row in enumerate(rows, start=1):
            path = RAW_DIR / "somatic_maf" / row["file_id"] / row["file_name"]
            if not path.exists():
                continue
            with open_text(path) as handle:
                reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
                if reader.fieldnames is None:
                    continue
                if header is None:
                    header = list(reader.fieldnames) + ["Patient_Barcode"]
                    writer = csv.DictWriter(out, fieldnames=header, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                assert writer is not None
                for rec in reader:
                    sample = rec.get("Tumor_Sample_Barcode") or row.get("sample_barcode") or ""
                    patient = sample[:12] if sample else row.get("patient_barcode")
                    rec["Patient_Barcode"] = patient
                    writer.writerow({key: rec.get(key, "") for key in header})
                    classification = rec.get("Variant_Classification") or "Unknown"
                    counts[patient]["all_mutations"] += 1
                    if classification not in nonsilent_exclude:
                        counts[patient]["nonsilent_mutations"] += 1
                    if rec.get("Hugo_Symbol"):
                        counts[patient][f"gene::{rec['Hugo_Symbol']}"] += 1
            if idx % 100 == 0 or idx == len(rows):
                log(f"  merged MAF {idx}/{len(rows)}")

    count_rows = []
    for patient, counter in sorted(counts.items()):
        count_rows.append(
            {
                "patient_barcode": patient,
                "all_mutations": counter["all_mutations"],
                "nonsilent_mutations": counter["nonsilent_mutations"],
            }
        )
    write_tsv(count_rows, RESULT_TABLE_DIR / "mutation_counts_by_patient.tsv")
    log(f"[DONE] MAF merged -> {out_maf}")


def build_cnv_segments() -> None:
    rows = read_tsv(META_DIR / "masked_cnv_segments_metadata.tsv")
    if not rows:
        log("[SKIP] No CNV metadata found.")
        return
    out_path = PROCESSED_DIR / "blca_masked_cnv_segments.tsv.gz"
    header_written = False
    total_records = 0
    with gzip_writer(out_path) as out:
        writer = None
        for idx, row in enumerate(rows, start=1):
            path = RAW_DIR / "masked_cnv_segments" / row["file_id"] / row["file_name"]
            if not path.exists():
                continue
            with open_text(path) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames is None:
                    continue
                extra = ["patient_barcode", "sample_barcode", "source_file_id", "source_file_name"]
                if not header_written:
                    writer = csv.DictWriter(
                        out, fieldnames=list(reader.fieldnames) + extra, delimiter="\t", lineterminator="\n"
                    )
                    writer.writeheader()
                    header_written = True
                assert writer is not None
                for rec in reader:
                    rec.update(
                        {
                            "patient_barcode": row.get("patient_barcode"),
                            "sample_barcode": row.get("sample_barcode"),
                            "source_file_id": row.get("file_id"),
                            "source_file_name": row.get("file_name"),
                        }
                    )
                    writer.writerow(rec)
                    total_records += 1
            if idx % 200 == 0 or idx == len(rows):
                log(f"  merged CNV {idx}/{len(rows)}")
    log(f"[DONE] CNV segments records={total_records} -> {out_path}")


def build_rppa_matrix() -> None:
    pd = require_pandas()
    rows = read_tsv(META_DIR / "rppa_protein_metadata.tsv")
    if not rows:
        log("[SKIP] No RPPA metadata found.")
        return
    file_to_col = unique_column_names(rows)
    series = []
    sample_map = []
    for row in rows:
        path = RAW_DIR / "rppa_protein" / row["file_id"] / row["file_name"]
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
        cols = list(df.columns)
        id_col = next((c for c in cols if c.lower() in {"composite.element.ref", "protein", "protein_name", "target"}), cols[0])
        value_col = next((c for c in cols if "expression" in c.lower() or "value" in c.lower()), cols[-1])
        col = file_to_col[row["file_id"]]
        series.append(df.set_index(id_col)[value_col].rename(col))
        sample_map.append({**row, "matrix_column": col, "protein_id_column": id_col, "protein_value_column": value_col})
    if series:
        pd.concat(series, axis=1).to_csv(PROCESSED_DIR / "blca_rppa_matrix.tsv.gz", sep="\t", compression="gzip")
        write_tsv(sample_map, PROCESSED_DIR / "blca_rppa_sample_map.tsv")
        log("[DONE] RPPA matrix built.")


def build_all() -> None:
    ensure_dirs()
    query_clinical(PROJECT)
    build_expression_matrices()
    build_maf()
    build_cnv_segments()
    build_rppa_matrix()


def parse_datasets(value: str) -> list[str]:
    if value == "core":
        return CORE_DATASETS
    if value == "all":
        return list(DATASETS)
    datasets = [v.strip() for v in value.split(",") if v.strip()]
    invalid = [d for d in datasets if d not in DATASETS]
    if invalid:
        raise SystemExit(f"Unknown dataset(s): {', '.join(invalid)}. Choices: {', '.join(DATASETS)}")
    return datasets


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="TCGA-BLCA open-data download and processing pipeline.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_query = sub.add_parser("query", help="Query GDC metadata and write manifests.")
    p_query.add_argument("--project", default=PROJECT)
    p_query.add_argument("--datasets", default="core", help="core, all, or comma-separated dataset names")

    p_download = sub.add_parser("download", help="Download files listed in metadata manifests.")
    p_download.add_argument("--datasets", default="core", help="core, all, or comma-separated dataset names")
    p_download.add_argument("--threads", type=int, default=4)

    sub.add_parser("clinical", help="Download/flatten GDC clinical data.")
    sub.add_parser("build", help="Build processed matrices/tables from downloaded raw files.")
    sub.add_parser("all", help="Run query, download core datasets, and build processed outputs.")

    args = parser.parse_args(argv)
    ensure_dirs()
    if args.command == "query":
        query_files(args.project, parse_datasets(args.datasets))
        query_clinical(args.project)
    elif args.command == "download":
        download_files(parse_datasets(args.datasets), threads=args.threads)
    elif args.command == "clinical":
        query_clinical(PROJECT)
    elif args.command == "build":
        build_all()
    elif args.command == "all":
        query_files(PROJECT, CORE_DATASETS)
        query_clinical(PROJECT)
        download_files(CORE_DATASETS, threads=4)
        build_all()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
