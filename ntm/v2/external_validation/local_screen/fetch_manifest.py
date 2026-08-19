#!/usr/bin/env python3
"""Fetch preregistered EUtils SRA manifests for the local NTM screen.

RUN_PLAN §5 (Tier A) and §6 (Tier B). Executes the exact preregistered
esearch terms, saves uid lists, computes sha256 of the frozen uid list,
and (optionally, --runinfo) pulls organism + layout metadata for filtering
(attribution only — the uid set itself is frozen at query time).

Rate limit: <= 1 req/s (no API key), 3 retries with backoff per request.
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

TIER_A = [
    ("A1", '"Mycobacterium smegmatis"[Organism] AND mitomycin[All Fields]'),
    ("A2", '"Mycobacterium smegmatis"[Organism] AND prophage[All Fields]'),
    ("A3", '"Mycobacterium smegmatis"[Organism] AND temperate[All Fields]'),
    ("A4", '"Mycobacterium"[Organism] AND mitomycin[All Fields]'),
]
# Amendment A1 (AMENDMENTS.md): compound prophage x organism queries currently
# return 0 on EUtils; A2 is executed as the bare term + LOCAL smegmatis filter.
TIER_A_BARE = [("A2bare", "prophage[All Fields]")]
# Amendment A3 (AMENDMENTS.md): stable induction terms completing Tier A
TIER_A_AMEND3 = [
    ("A5", '"Mycobacterium smegmatis"[Organism] AND ciprofloxacin[All Fields]'),
    ("A6", '"Mycobacterium smegmatis"[Organism] AND induction[All Fields]'),
    ("A7", '"Mycobacterium smegmatis"[Organism] AND induced[All Fields]'),
    ("A8", '"Mycobacterium smegmatis"[Organism] AND phage[All Fields]'),
]
TIER_B = [
    ("B1", '"Mycobacterium abscessus"[Organism]'),
    ("B2", '"Mycobacterium avium"[Organism]'),
    ("B3", '"Mycobacterium smegmatis"[Organism]'),
]


def http_get(url, retries=3, backoff=(5, 20, 60)):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "local-ntm-logan/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read(), dict(r.headers)
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < retries:
                time.sleep(backoff[min(attempt, len(backoff) - 1)])
    raise RuntimeError(f"GET failed after {retries} retries: {url}: {last}")


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def esearch(term, log):
    q = urllib.parse.urlencode({"db": "sra", "term": term, "retmax": 100000,
                                "retmode": "json", "sort": "accession"})
    body, hdrs = http_get(f"{EUTILS}/esearch.fcgi?{q}")
    data = json.loads(body)["esearchresult"]
    uids = data.get("idlist", [])
    log.write(json.dumps({"event": "esearch", "term": term, "count": data.get("count"),
                          "returned": len(uids), "sha256": sha256_bytes(body),
                          "querytranslation": data.get("querytranslation")}) + "\n")
    time.sleep(1.1)  # rate limit
    return uids, data.get("count")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", choices=["A", "B"], required=True)
    ap.add_argument("--out", required=True, help="manifest json path")
    ap.add_argument("--log", required=True, help="append-only eutils log")
    args = ap.parse_args()

    queries = TIER_A if args.tier == "A" else TIER_B
    if args.tier == "A":
        queries = queries + TIER_A_BARE + TIER_A_AMEND3  # Amendments A1 + A3
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    manifest = {
        "tier": args.tier,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "db": "sra",
        "queries": [],
        "note": "preregistered in RUN_PLAN.md before any download; flat queries only "
                "(nested-paren terms return 0 on db=sra — see logs/eutils_calibration.log)",
    }
    union = {}
    with open(args.log, "a") as log:
        for qid, term in queries:
            uids, count = esearch(term, log)
            uid_bytes = ("\n".join(uids) + "\n").encode()
            manifest["queries"].append({
                "id": qid, "term": term, "reported_count": count,
                "returned_uids": len(uids), "uid_list_sha256": sha256_bytes(uid_bytes),
                "status": "amended_see_AMENDMENTS_A1" if qid == "A2" else "ok",
            })
            for u in uids:
                union.setdefault(u, []).append(qid)
        union_bytes = json.dumps(union, sort_keys=True).encode()
        manifest["union_size"] = len(union)
        manifest["union_sha256"] = sha256_bytes(union_bytes)
        manifest["uids"] = {u: qs for u, qs in sorted(union.items(), key=lambda kv: int(kv[0]))}
    with open(args.out, "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    print(json.dumps({k: manifest[k] for k in ("tier", "union_size", "union_sha256")}))
    # also drop a plain accession file for the screen driver
    acc_path = args.out.replace(".json", "_uids.txt")
    with open(acc_path, "w") as fh:
        fh.write("\n".join(sorted(union, key=int)) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
