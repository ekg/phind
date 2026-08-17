#!/usr/bin/env python3
"""
sync_repo.py — copy panel summaries + provenance into the repository.

The bulky FASTA and the full manifest stay in external storage
(/mnt/nvme3n1/.../external_validation/panel/); the repo gets:

  REPORT.md, PROVENANCE.md            (verbatim copies of the build outputs)
  summaries/manifest_panel.tsv.gz     (full manifest, gzipped ~500 kB)
  summaries/crosswalk.tsv.gz
  summaries/exclusions.tsv            (small, plain)
  summaries/source_summary.tsv
  summaries/host_distribution.tsv
  summaries/CHECKSUMS.sha256          (sha256 of every panel artifact, incl. FASTA)
  summaries/panel_files.txt           (names, sizes of external panel files)

Usage: python3 sync_repo.py [--data-dir DIR] [--repo-dir auto]
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_DIR = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation"
REPO_EXTERNAL_DIR = os.path.join(os.path.dirname(HERE), "summaries")

PLAIN = ["exclusions.tsv", "source_summary.tsv", "host_distribution.tsv"]
DOCS = ["REPORT.md", "PROVENANCE.md"]


def gzip_copy(src: str, dst: str) -> None:
    import gzip as _gz
    data = open(src, "rb").read()
    with open(dst, "wb") as fh:
        with _gz.GzipFile(fileobj=fh, mode="wb", compresslevel=9, mtime=0) as fo:
            fo.write(data)


def sha256_file(path: str) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    ap.add_argument("--repo-dir", default=REPO_EXTERNAL_DIR)
    args = ap.parse_args(argv)
    panel = os.path.join(args.data_dir, "panel")
    os.makedirs(args.repo_dir, exist_ok=True)

    for name in DOCS:
        shutil.copyfile(os.path.join(panel, name),
                        os.path.join(args.repo_dir, name))
    for name in PLAIN:
        shutil.copyfile(os.path.join(panel, name),
                        os.path.join(args.repo_dir, name))
    for name in ["manifest_panel.tsv", "crosswalk.tsv"]:
        gzip_copy(os.path.join(panel, name),
                  os.path.join(args.repo_dir, name + ".gz"))

    # checksums of every external panel artifact
    lines = []
    for fn in sorted(os.listdir(panel)):
        p = os.path.join(panel, fn)
        if os.path.isfile(p):
            lines.append(f"{sha256_file(p)}  panel/{fn}")
    open(os.path.join(args.repo_dir, "CHECKSUMS.sha256"), "w").write(
        "\n".join(lines) + "\n")
    open(os.path.join(args.repo_dir, "panel_files.txt"), "w").write(
        "\n".join(f"{os.path.getsize(os.path.join(panel, fn))}\t{fn}"
                  for fn in sorted(os.listdir(panel))
                  if os.path.isfile(os.path.join(panel, fn))) + "\n")
    print(f"synced {len(lines)} panel artifacts -> {args.repo_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
