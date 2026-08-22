#!/usr/bin/env python3
"""Fetch the E. coli WGS SRA universe manifest (RUN_PLAN §3, §4).

Frozen query (calibrated 2026-08-22, receipts in logs/eutils_calibration.log):
    "Escherichia coli"[Organism] AND wgs[strategy]

esearch with usehistory + retstart pagination (retmax 100000; deep paging
verified at retstart=500000 pre-preregistration). Saves the uid list,
computes sha256 of the frozen list, writes metadata/ecoli_universe_manifest.json.
Rate limit <= 1 req/s, 3 retries, 5/20/60 s backoff.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Amendment-free frozen query (RUN_PLAN §3)
UNIVERSE_TERM = '"Escherichia coli"[Organism] AND wgs[strategy]'
PAGE = 100000


def http_get(url, retries=3, backoff=(5, 20, 60)):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "local-ecoli-logan/1.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < retries:
                time.sleep(backoff[min(attempt, len(backoff) - 1)])
    raise RuntimeError(f"GET failed after {retries} retries: {url}: {last}")


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="metadata dir (NVMe)")
    ap.add_argument("--log", required=True, help="append-only eutils log (repo)")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    uid_path = os.path.join(args.out_dir, "ecoli_universe_uids.txt")
    done = []
    if os.path.exists(uid_path):
        done = [u.strip() for u in open(uid_path) if u.strip()]
        print(f"resuming: {len(done)} uids already saved")

    with open(args.log, "a") as log:
        if not done:
            q = urllib.parse.urlencode({"db": "sra", "term": UNIVERSE_TERM,
                                        "retmax": 0, "retmode": "json",
                                        "usehistory": "y"})
            body = http_get(f"{EUTILS}/esearch.fcgi?{q}")
            head = json.loads(body)["esearchresult"]
            count = int(head["count"])
            webenv, qkey = head["webenv"], head["querykey"]
            log.write(json.dumps({"event": "esearch_head", "term": UNIVERSE_TERM,
                                  "count": count, "webenv": webenv[:32] + "…",
                                  "query_key": qkey, "sha256": sha256_bytes(body),
                                  "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}) + "\n")
            print(f"universe count: {count}")
        else:
            # resume: re-open history for pagination
            q = urllib.parse.urlencode({"db": "sra", "term": UNIVERSE_TERM,
                                        "retmax": 0, "retmode": "json",
                                        "usehistory": "y"})
            body = http_get(f"{EUTILS}/esearch.fcgi?{q}")
            head = json.loads(body)["esearchresult"]
            count = int(head["count"])
            webenv, qkey = head["webenv"], head["querykey"]

        uids = done
        retstart = len(done)
        while retstart < count:
            q = urllib.parse.urlencode({"db": "sra", "query_key": qkey,
                                        "WebEnv": webenv, "retstart": retstart,
                                        "retmax": PAGE, "retmode": "json"})
            body = http_get(f"{EUTILS}/esearch.fcgi?{q}")
            page = json.loads(body)["esearchresult"].get("idlist", [])
            log.write(json.dumps({"event": "esearch_page", "retstart": retstart,
                                  "returned": len(page), "sha256": sha256_bytes(body),
                                  "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}) + "\n")
            if not page:
                print(f"no more results at retstart={retstart}; stopping")
                break
            uids.extend(page)
            with open(uid_path, "a") as fh:
                fh.write("\n".join(page) + "\n")
            retstart += len(page)
            print(f"  {retstart}/{count}", flush=True)
            time.sleep(1.1)

    uids = sorted(set(uids), key=int)
    uid_bytes = ("\n".join(uids) + "\n").encode()
    manifest = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "db": "sra",
        "universe_query": UNIVERSE_TERM,
        "reported_count": count,
        "frozen_uid_count": len(uids),
        "uid_list_sha256": sha256_bytes(uid_bytes),
        "note": "RUN_PLAN §3 frozen universe query; uid list frozen before any "
                "tier design or download; pagination receipts in eutils log",
    }
    mpath = os.path.join(args.out_dir, "ecoli_universe_manifest.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    # rewrite the uid file deduplicated + sorted (freeze order)
    with open(uid_path, "w") as fh:
        fh.write("\n".join(uids) + "\n")
    print(json.dumps(manifest))
    return 0


if __name__ == "__main__":
    sys.exit(main())
