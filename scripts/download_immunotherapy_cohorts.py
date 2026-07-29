#!/usr/bin/env python
"""Download IMvigor210 and GSE176307 immunotherapy validation resources."""

from __future__ import annotations

import json
import ssl
import time
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "external" / "immunotherapy"

FILES = {
    "IMvigor210": {
        "IMvigor210CoreBiologies_1.0.1.tar.gz": (
            "http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/"
            "IMvigor210CoreBiologies_1.0.1.tar.gz"
        ),
    },
    "GSE176307": {
        "GSE176307_series_matrix.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/matrix/"
            "GSE176307_series_matrix.txt.gz"
        ),
        "GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/suppl/"
            "GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz"
        ),
        "GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/suppl/"
            "GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz"
        ),
        "GSE176307_salmon_tpm_gene.matrix.tsv.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/suppl/"
            "GSE176307_salmon_tpm_gene.matrix.tsv.gz"
        ),
        "GSE176307_baci_rsem_RS_BACI_headers_tab.txt.gz": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176307/suppl/"
            "GSE176307_baci_rsem_RS_BACI_headers_tab.txt.gz"
        ),
    },
}


def ssl_context() -> ssl.SSLContext:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def download(url: str, path: Path, tries: int = 4) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        print(f"exists\t{path.name}\t{path.stat().st_size}")
        return {"file": str(path), "downloaded": False, "size_bytes": path.stat().st_size, "url": url}

    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(1, tries + 1):
        try:
            print(f"download\t{path.name}\tattempt={attempt}")
            req = urllib.request.Request(url, headers={"User-Agent": "TCGA-BLCA-analysis/1.0"})
            with urllib.request.urlopen(req, context=ssl_context(), timeout=120) as resp, tmp.open("wb") as fh:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    fh.write(chunk)
            tmp.replace(path)
            print(f"done\t{path.name}\t{path.stat().st_size}")
            return {"file": str(path), "downloaded": True, "size_bytes": path.stat().st_size, "url": url}
        except Exception as exc:  # noqa: BLE001
            print(f"retry\t{path.name}\t{exc}")
            if tmp.exists():
                tmp.unlink()
            if attempt == tries:
                raise
            time.sleep(2 * attempt)
    raise RuntimeError(f"Failed to download {url}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {}
    for cohort, files in FILES.items():
        cdir = OUT / cohort
        cdir.mkdir(parents=True, exist_ok=True)
        manifest[cohort] = {}
        for filename, url in files.items():
            manifest[cohort][filename] = download(url, cdir / filename)
    manifest_path = OUT / "immunotherapy_download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"manifest\t{manifest_path}")


if __name__ == "__main__":
    main()
