#!/usr/bin/env python
"""Download and catalogue selected BLCA single-cell/spatial GEO datasets."""

from __future__ import annotations

import json
import ssl
import time
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SC_ROOT = ROOT / "data" / "external" / "single_cell"

DATASETS = {
    "GSE222315": {
        "type": "scRNA",
        "series": "GSE222nnn",
        "files": {
            "RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222315/suppl/GSE222315_RAW.tar",
            "filelist.txt": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222315/suppl/filelist.txt",
        },
    },
    "GSE171351": {
        "type": "spatial_visium",
        "series": "GSE171nnn",
        "files": {
            "combined_visium.h5ad.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171351/suppl/GSE171351_combined_visium.h5ad.gz",
            "RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171351/suppl/GSE171351_RAW.tar",
            "filelist.txt": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171351/suppl/filelist.txt",
        },
    },
    "GSE269877": {
        "type": "scRNA_large_rds",
        "series": "GSE269nnn",
        "files": {
            "directory_listing.html": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269877/suppl/",
        },
        "deferred_large_files": {
            "GSE269877_dta_cancer.submission.rds.gz": "13G",
            "GSE269877_dta_normal.submission.rds.gz": "2.4G",
        },
    },
    "GSE293189": {
        "type": "scRNA_dge",
        "series": "GSE293nnn",
        "files": {
            "RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293189/suppl/GSE293189_RAW.tar",
            "filelist.txt": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293189/suppl/filelist.txt",
        },
    },
}


def ssl_context() -> ssl.SSLContext:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def url_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "TCGA-BLCA-analysis/1.0"})
    try:
        with urllib.request.urlopen(req, context=ssl_context(), timeout=60) as r:
            length = r.headers.get("Content-Length")
            return int(length) if length is not None else None
    except Exception:
        return None


def download(url: str, path: Path, tries: int = 4) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    expected = url_size(url)
    if path.exists() and path.stat().st_size > 0:
        if expected is None or path.stat().st_size == expected:
            print(f"exists\t{path}")
            return {"file": str(path), "downloaded": False, "size_bytes": path.stat().st_size, "expected_bytes": expected}

    tmp = path.with_suffix(path.suffix + ".part")
    headers = {"User-Agent": "TCGA-BLCA-analysis/1.0"}
    for attempt in range(1, tries + 1):
        try:
            print(f"download\t{path.name}\tattempt={attempt}")
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, context=ssl_context(), timeout=120) as r, tmp.open("wb") as fh:
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    fh.write(chunk)
            tmp.replace(path)
            print(f"done\t{path.name}\t{path.stat().st_size}")
            return {"file": str(path), "downloaded": True, "size_bytes": path.stat().st_size, "expected_bytes": expected}
        except Exception as exc:  # noqa: BLE001
            print(f"retry\t{path.name}\t{exc}")
            if tmp.exists():
                tmp.unlink()
            if attempt == tries:
                raise
            time.sleep(2 * attempt)
    raise RuntimeError(f"Failed to download {url}")


def main() -> None:
    SC_ROOT.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {}
    for acc, info in DATASETS.items():
        print(f"\n### {acc}")
        ddir = SC_ROOT / acc
        ddir.mkdir(parents=True, exist_ok=True)
        outputs: dict[str, object] = {}
        for label, url in info["files"].items():
            if label == "directory_listing.html":
                fname = "supplementary_directory_listing.html"
            elif label == "RAW.tar":
                fname = f"{acc}_RAW.tar"
            else:
                fname = f"{acc}_{label}" if not label.startswith(acc) else label
            outputs[label] = download(url, ddir / fname)
        summary[acc] = {
            "type": info["type"],
            "geo": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
            "outputs": outputs,
            "deferred_large_files": info.get("deferred_large_files", {}),
        }
    out = SC_ROOT / "candidate_geo_download_summary.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nsummary\t{out}")


if __name__ == "__main__":
    main()
