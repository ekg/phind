#!/usr/bin/env python3
"""Freeze a tier manifest: uid list + runinfo -> final run set with filters.

Tier A (per RUN_PLAN §5 + AMENDMENTS.md A1):
  - NTM filter on scientific_name for all queries
  - A2bare uids additionally restricted to M. smegmatis (preregistered A2
    organism clause, evaluated locally per Amendment A1)
Tier B: organism-only queries, no extra filter.

Output: frozen TSV (run, scientific_name, queries, library_strategy,
platform, bioproject) + freeze receipt block appended to manifest json.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys

NTM_EXCLUDE = re.compile(r"tuberculos|bovis|africanum|canettii|leprae|"
                         r"caprae|pinnipedii|mungi|orycis", re.I)
# modern NTM genera after the Mycobacterium sensu lato split
NTM_GENERA = ("Mycobacterium", "Mycolicibacterium", "Mycobacteroides",
              "Mycolicibacillus")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--runinfo", required=True)
    ap.add_argument("--tier", choices=["A", "B"], required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    m = json.load(open(args.manifest))
    uid_queries = {u: set(qs) for u, qs in m["uids"].items()}

    rows = {}
    for r in csv.DictReader(open(args.runinfo), delimiter="\t"):
        qs_for_uid = set(uid_queries.get(r["uid"], []))
        for run in (r.get("runs") or "").split(","):
            if not run:
                continue
            e = rows.setdefault(run, {"run": run,
                                      "scientific_name": r["scientific_name"],
                                      "library_strategy": r["library_strategy"],
                                      "platform": r["platform"],
                                      "bioproject": r["bioproject"], "queries": set(),
                                      "uids": set()})
            e["queries"] |= qs_for_uid
            e["uids"].add(r["uid"])

    kept, dropped = [], []
    for e in rows.values():
        sci = e["scientific_name"] or ""
        keep = True
        if not sci.startswith(NTM_GENERA):
            keep = False
        elif NTM_EXCLUDE.search(sci):
            keep = False
        elif args.tier == "A" and "A2bare" in e["queries"] and "smegmatis" not in sci.lower():
            # preregistered A2 organism clause, evaluated locally (Amendment A1)
            keep = False
        (kept if keep else dropped).append(e)

    kept.sort(key=lambda e: e["run"])
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "scientific_name", "queries", "library_strategy",
                    "platform", "bioproject", "uids"])
        for e in kept:
            w.writerow([e["run"], e["scientific_name"],
                        ",".join(sorted(e["queries"])), e["library_strategy"],
                        e["platform"], e["bioproject"],
                        ",".join(sorted(e["uids"]))])
    frozen = open(args.out, "rb").read()
    receipt = {
        "freeze_utc": __import__("time").strftime("%Y-%m-%dT%H:%M:%SZ",
                                                  __import__("time").gmtime()),
        "n_runs_runinfo": len(rows), "n_runs_kept": len(kept),
        "n_runs_dropped_filter": len(dropped),
        "frozen_tsv_sha256": hashlib.sha256(frozen).hexdigest(),
        "filters": ("NTM: scientific_name genus in Mycobacterium|Mycolicibacterium|"
                    "Mycobacteroides|Mycolicibacillus and not match tuberculos|bovis|"
                    "africanum|canettii|leprae|caprae|pinnipedii|mungi|orycis"
                    + ("; A2bare restricted to smegmatis (Amendment A1)"
                       if args.tier == "A" else "")),
    }
    m.setdefault("freeze", []).append(receipt)
    with open(args.manifest, "w") as fh:
        json.dump(m, fh, indent=1, sort_keys=True)
    print(json.dumps(receipt))
    return 0


if __name__ == "__main__":
    sys.exit(main())
