#!/usr/bin/env python3
"""Emit NVMe artifact manifest (RUN_PLAN §10): every bulky artifact under the
NVMe run root with size + sha256 (downloads optional, size-dominated), so the
repo stays small while everything remains traceable.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="NVMe run root")
    ap.add_argument("--out", required=True, help="repo-side TSV")
    ap.add_argument("--hash-downloads", action="store_true",
                    help="hash the downloads/ tree too (slow; ledger already "
                         "has sha256 per file)")
    args = ap.parse_args()
    rows = []
    for dirpath, _dirnames, filenames in os.walk(args.root):
        if "env" in dirpath.split(os.sep):  # conda env: versioned via VERSIONS.txt
            continue
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, args.root)
            if rel.startswith("downloads/") and not args.hash_downloads:
                sha = "see_download_ledger"
            elif rel.startswith(("screen/scratch", "confirm/scratch")):
                continue  # transient scratch
            elif rel == "downloads/download_ledger.jsonl":
                sha = sha256_file(p)
            else:
                sha = sha256_file(p)
            rows.append((rel, os.path.getsize(p), sha))
    rows.sort()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["relative_path", "size_bytes", "sha256_or_ref"])
        w.writerows(rows)
    print(f"{len(rows)} artifacts -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
