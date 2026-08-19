#!/usr/bin/env python3
"""uid -> run accession (+ metadata) via EUtils efetch runinfo.

Rate limit <= 1 req/s, retries with backoff. Output TSV:
uid, run, scientific_name, library_strategy, library_source, platform, bioproject
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def http_get(url, retries=3, backoff=(5, 20, 60)):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "local-ntm-logan/1.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < retries:
                time.sleep(backoff[min(attempt, len(backoff) - 1)])
    raise RuntimeError(f"GET failed after {retries} retries: {url}: {last}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uids", required=True, help="file with one numeric uid per line")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--log", required=True)
    ap.add_argument("--chunk", type=int, default=400)
    args = ap.parse_args()

    uids = [u.strip() for u in open(args.uids) if u.strip()]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    done_uids = set()
    if os.path.exists(args.out):
        with open(args.out) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            done_uids = {r["uid"] for r in rd}
    todo = [u for u in uids if u not in done_uids]
    print(f"{len(uids)} uids, {len(done_uids)} already done, {len(todo)} to fetch")

    cols = ["uid", "run", "scientific_name", "library_strategy", "library_source",
            "platform", "bioproject"]
    write_header = not os.path.exists(args.out)
    with open(args.out, "a", newline="") as out, open(args.log, "a") as log:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=cols, lineterminator="\n")
        if write_header:
            w.writeheader()
        for i in range(0, len(todo), args.chunk):
            batch = todo[i:i + args.chunk]
            q = urllib.parse.urlencode({"db": "sra", "id": ",".join(batch),
                                        "rettype": "runinfo", "retmode": "csv"})
            body = http_get(f"{EUTILS}/efetch.fcgi?{q}")
            txt = body.decode("utf-8", "replace")
            rows = list(csv.DictReader(io.StringIO(txt)))
            log.write(json.dumps({"event": "efetch_runinfo", "n_uids": len(batch),
                                  "n_rows": len(rows),
                                  "sha256": __import__("hashlib").sha256(body).hexdigest(),
                                  "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}) + "\n")
            seen_uids = set()
            for r in rows:
                run = (r.get("Run") or "").strip()
                if not run:
                    continue
                w.writerow({
                    "uid": r.get("eFetchResult_num" ) or "",
                    "run": run,
                    "scientific_name": (r.get("ScientificName") or "").strip(),
                    "library_strategy": (r.get("LibraryStrategy") or "").strip(),
                    "library_source": (r.get("LibrarySource") or "").strip(),
                    "platform": (r.get("Platform") or "").strip(),
                    "bioproject": (r.get("BioProject") or "").strip(),
                })
                seen_uids.add(r.get("eFetchResult_num"))
            out.flush()
            time.sleep(1.1)
            if i % (args.chunk * 10) == 0:
                print(f"  {i + len(batch)}/{len(todo)} uids -> {len(rows)} rows", flush=True)
    print(f"done -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
