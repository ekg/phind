#!/usr/bin/env python3
"""uid -> run accessions + metadata via EUtils esummary (db=sra, JSON).

esummary is used instead of efetch runinfo because:
  - exact uid -> run mapping (efetch runinfo CSV is headerless, no uid col)
  - organism/strategy live in expxml; runs blob lists Run accessions
Rate limit <= 1 req/s, retries with backoff. Output TSV (one row per uid):
uid, runs, scientific_name, library_strategy, library_source, platform,
bioproject, n_runs
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RUN_ACC = re.compile(r'Run acc="([A-Z]+[0-9.]+)"')
SCI_NAME = re.compile(r'ScientificName="([^"]+)"')
LIB_STRAT = re.compile(r"<LIBRARY_STRATEGY>([^<]+)<")
LIB_SRC = re.compile(r"<LIBRARY_SOURCE>([^<]+)<")
PLATFORM = re.compile(r"<Platform[^>]*>([A-Z]+)</Platform>")
STUDY = re.compile(r'<Study acc="([^"]+)"')


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


def parse_uid_entry(e):
    expxml = e.get("expxml", "") or ""
    runs = RUN_ACC.findall(e.get("runs", "") or "")
    g = lambda rx: (rx.search(expxml).group(1) if rx.search(expxml) else "")  # noqa: E731
    return runs, g(SCI_NAME), g(LIB_STRAT), g(LIB_SRC), g(PLATFORM), g(STUDY)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uids", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--chunk", type=int, default=250)
    args = ap.parse_args()

    uids = [u.strip() for u in open(args.uids) if u.strip()]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    done_uids = set()
    if os.path.exists(args.out):
        with open(args.out) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            done_uids = {r["uid"] for r in rd if r.get("uid")}
    todo = [u for u in uids if u not in done_uids]
    print(f"{len(uids)} uids, {len(done_uids)} done, {len(todo)} to fetch")

    cols = ["uid", "runs", "scientific_name", "library_strategy",
            "library_source", "platform", "bioproject", "n_runs"]
    write_header = not os.path.exists(args.out)
    with open(args.out, "a", newline="") as out, open(args.log, "a") as log:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=cols, lineterminator="\n")
        if write_header:
            w.writeheader()
        for i in range(0, len(todo), args.chunk):
            batch = todo[i:i + args.chunk]
            q = urllib.parse.urlencode({"db": "sra", "id": ",".join(batch),
                                        "retmode": "json"})
            body = http_get(f"{EUTILS}/esummary.fcgi?{q}")
            data = json.loads(body).get("result", {})
            n_rows = 0
            for uid in batch:
                e = data.get(uid)
                if not e:
                    continue
                runs, sci, strat, src, plat, study = parse_uid_entry(e)
                if not runs:
                    continue
                w.writerow({"uid": uid, "runs": ",".join(runs),
                            "scientific_name": sci, "library_strategy": strat,
                            "library_source": src, "platform": plat,
                            "bioproject": study, "n_runs": len(runs)})
                n_rows += 1
            out.flush()
            log.write(json.dumps({"event": "esummary", "n_uids": len(batch),
                                  "n_rows": n_rows,
                                  "sha256": hashlib.sha256(body).hexdigest(),
                                  "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                      time.gmtime())}) + "\n")
            time.sleep(1.1)
            if i and i % (args.chunk * 8) == 0:
                print(f"  {i + len(batch)}/{len(todo)}", flush=True)
    print(f"done -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
