#!/usr/bin/env python
"""Download GEO single-cell files for the BLCA GSE145137 pilot analysis."""

from __future__ import annotations

import ssl
import sys
import time
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "external" / "single_cell" / "GSE145137"

FILES = {
    "celltype_xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145137/suppl/"
        "GSE145137_BC159-T_3_All_SC_final_QC_Celltype_Information.xlsx"
    ),
    "human_log2tpm": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145137/suppl/"
        "GSE145137_GEO_processed_BC159-T_3_PDX_human_log2TPM_matrix_final.txt.gz"
    ),
    "mouse_log2tpm": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145137/suppl/"
        "GSE145137_GEO_processed_BC159-T_3_PDX_mouse_log2TPM_matrix_final.txt.gz"
    ),
    "raw_tar": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE145137&format=file",
}


def download(url: str, path: Path, tries: int = 4) -> None:
    if path.exists() and path.stat().st_size > 0:
        print(f"exists\t{path.name}\t{path.stat().st_size}")
        return

    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(1, tries + 1):
        try:
            print(f"download\t{path.name}\tattempt={attempt}")
            req = urllib.request.Request(url, headers={"User-Agent": "TCGA-BLCA-analysis/1.0"})
            with urllib.request.urlopen(req, context=ctx, timeout=90) as r, tmp.open("wb") as fh:
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    fh.write(chunk)
            tmp.replace(path)
            print(f"done\t{path.name}\t{path.stat().st_size}")
            return
        except Exception as exc:  # noqa: BLE001
            print(f"retry\t{path.name}\t{exc}", file=sys.stderr)
            if tmp.exists():
                tmp.unlink()
            if attempt == tries:
                raise
            time.sleep(2 * attempt)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for key, url in FILES.items():
        suffix = "GSE145137_RAW.tar" if key == "raw_tar" else url.rsplit("/", 1)[-1]
        download(url, OUT / suffix)

    manifest = OUT / "download_manifest.tsv"
    with manifest.open("w", encoding="utf-8") as fh:
        fh.write("key\tfile\turl\tsize_bytes\n")
        for key, url in FILES.items():
            path = OUT / ("GSE145137_RAW.tar" if key == "raw_tar" else url.rsplit("/", 1)[-1])
            fh.write(f"{key}\t{path.name}\t{url}\t{path.stat().st_size if path.exists() else ''}\n")
    print(f"manifest\t{manifest}")


if __name__ == "__main__":
    main()
