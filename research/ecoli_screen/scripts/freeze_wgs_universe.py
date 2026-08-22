#!/usr/bin/env python3
"""Freeze tier manifests for the E. coli local screen (RUN_PLAN §4, §5, §6).

Filter (RUN_PLAN §4, frozen):
  - library_strategy == "WGS"
  - scientific_name startswith "Escherichia coli"
Output: metadata/ecoli_wgs_frozen.tsv (run, scientific_name, library_strategy,
platform, bioproject, uid) + freeze receipt in the universe manifest json.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe-manifest", required=True)
    ap.add_argument("--runinfo", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = {}
    n_rows_seen = 0
    for r in csv.DictReader(open(args.runinfo), delimiter="\t"):
        n_rows_seen += 1
        sci = (r.get("scientific_name") or "")
        if (r.get("library_strategy") or "").strip().upper() != "WGS":
            continue
        if not sci.startswith("Escherichia coli"):
            continue
        for run in (r.get("runs") or "").split(","):
            run = run.strip()
            if not run:
                continue
            e = rows.setdefault(run, {"run": run, "scientific_name": sci,
                                      "library_strategy": r["library_strategy"],
                                      "platform": r.get("platform", ""),
                                      "bioproject": r.get("bioproject", ""),
                                      "uids": set()})
            e["uids"].add(r["uid"])

    kept = sorted(rows.values(), key=lambda e: e["run"])
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "scientific_name", "library_strategy", "platform",
                    "bioproject", "uids"])
        for e in kept:
            w.writerow([e["run"], e["scientific_name"], e["library_strategy"],
                        e["platform"], e["bioproject"],
                        ",".join(sorted(e["uids"]))])
    frozen = open(args.out, "rb").read()
    receipt = {
        "freeze_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_runinfo_rows": n_rows_seen, "n_runs_kept": len(kept),
        "frozen_tsv_sha256": hashlib.sha256(frozen).hexdigest(),
        "filters": ("library_strategy == WGS AND scientific_name startswith "
                    "'Escherichia coli' (RUN_PLAN §4, frozen before any "
                    "tier design or download)"),
    }
    m = json.load(open(args.universe_manifest))
    m.setdefault("freeze", []).append(receipt)
    with open(args.universe_manifest, "w") as fh:
        json.dump(m, fh, indent=1, sort_keys=True)
    print(json.dumps(receipt))
    return 0


if __name__ == "__main__":
    sys.exit(main())
