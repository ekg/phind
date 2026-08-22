#!/usr/bin/env python3
"""Emit the full-provenance confirmed-hits TSV (RUN_PLAN §10.3).

Every row: accession -> S3 URL -> zst sha256 -> bait -> kmer_cov -> contig id
+ contig length + target span + strand + CIGAR + identity + classification.
Reads confirm jsonl + availability.tsv. Deterministic order.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--confirm-jsonl", nargs="+", required=True)
    ap.add_argument("--availability", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sha_by_acc, url_by_acc = {}, {}
    for r in csv.DictReader(open(args.availability), delimiter="\t"):
        if r.get("sha256") and r.get("status") == "200":
            sha_by_acc.setdefault(r["accession"], r["sha256"])
        url_by_acc.setdefault(r["accession"],
                              f'https://s3.amazonaws.com/logan-pub/c/{r["accession"]}/'
                              f'{r["accession"]}.contigs.fa.zst')

    rows = []
    for path in args.confirm_jsonl:
        for line in open(path):
            rec = json.loads(line)
            acc = rec["accession"]
            for b, info in rec["baits"].items():
                aln = info.get("minimap2") or {}
                best = (aln.get("best") or {})
                rows.append({
                    "accession": acc, "tier": rec.get("tier", ""),
                    "scientific_name": rec.get("scientific_name", ""),
                    "bioproject": rec.get("bioproject", ""),
                    "s3_url": rec.get("s3_url") or url_by_acc.get(acc, ""),
                    "zst_sha256": rec.get("zst_sha256") or sha_by_acc.get(acc, ""),
                    "bait": b, "kmer_cov": info["kmer_cov"],
                    "qcov": info["qcov"], "identity": info["identity"],
                    "contig": info["best_contig"],
                    "contig_len": info["contig_len"],
                    "target_span": info["target_span"],
                    "strand": (best or {}).get("strand", ""),
                    "cigar": (best or {}).get("cigar", ""),
                    "classification": rec["classification"],
                    "ntm_strict_rule": rec.get("ntm_strict_rule_near_complete", False),
                })
    rows.sort(key=lambda r: (r["accession"], r["bait"]))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    cols = list(rows[0].keys()) if rows else ["accession"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
