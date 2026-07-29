#!/usr/bin/env python
"""Download and inspect GEO series matrix files for external validation."""

from __future__ import annotations

import gzip
import json
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "external" / "raw"
META = ROOT / "data" / "external" / "metadata"

COHORTS = ["GSE13507", "GSE31684"]


def geo_series_family(accession: str) -> str:
    digits = accession[3:]
    return f"GSE{digits[:2]}nnn"


def download(url: str, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 0:
        print(f"[SKIP] {out}")
        return
    tmp = out.with_suffix(out.suffix + ".part")
    print(f"[GET] {url}")
    with urllib.request.urlopen(url, timeout=300) as response, tmp.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    tmp.replace(out)


def parse_series_metadata(matrix_gz: Path) -> dict[str, object]:
    metadata: dict[str, list[str]] = {}
    sample_metadata: dict[str, list[str]] = {}
    table_header = None
    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                table_header = next(handle).rstrip("\n").split("\t")
                break
            if not line.startswith("!"):
                continue
            parts = [p.strip('"') for p in line.split("\t")]
            key = parts[0]
            values = parts[1:]
            if key.startswith("!Sample_"):
                sample_metadata[key] = values
            else:
                metadata[key] = values
    samples = sample_metadata.get("!Sample_geo_accession", [])
    platforms = sorted(set(metadata.get("!Series_platform_id", []) + sample_metadata.get("!Sample_platform_id", [])))
    characteristics = sample_metadata.get("!Sample_characteristics_ch1", [])
    char_preview = characteristics[: min(12, len(characteristics))]
    return {
        "matrix_file": str(matrix_gz),
        "n_samples": len(samples),
        "samples_preview": samples[:5],
        "platforms": platforms,
        "sample_metadata_keys": sorted(sample_metadata),
        "series_metadata_keys": sorted(metadata),
        "characteristics_preview": char_preview,
        "table_header_preview": table_header[:8] if table_header else None,
    }


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    META.mkdir(parents=True, exist_ok=True)
    summaries = {}
    for accession in COHORTS:
        family = geo_series_family(accession)
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{family}/{accession}/matrix/{accession}_series_matrix.txt.gz"
        out = RAW / accession / f"{accession}_series_matrix.txt.gz"
        download(url, out)
        summaries[accession] = parse_series_metadata(out)
        print(json.dumps({accession: summaries[accession]}, indent=2)[:5000])
    (META / "geo_series_download_summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
