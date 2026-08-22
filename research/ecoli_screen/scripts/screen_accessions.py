#!/usr/bin/env python3
"""Batch S3 logan-pub download + 31-mer bait screen driver (RUN_PLAN §3, §7).

Per accession (RUN accession, e.g. DRR000016):
  1. HEAD https://s3.amazonaws.com/logan-pub/c/<acc>/<acc>.contigs.fa.zst
     -> availability.tsv (200 + size | 404 not_in_logan)
  2. download (<=4 concurrent, 3 retries, 5/20/60 s backoff) + sha256 ledger
  3. back_to_sequences --in-kmers panel.fa --in-sequences <acc>.contigs.fa.zst
     --out-kmers <scratch>   (canonical 31-mers, both strands)
  4. per-bait kmer_coverage = |distinct bait 31-mers found| / |distinct bait
     31-mers|  -> screen/<tier>_results.jsonl

Budgets (RUN_PLAN §7): 400 GB downloaded, 25k accessions, 48 h wall,
concurrency <= 4. Resumable: ledger + results drive skipping.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

S3_BASE = "https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst"

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def canon(k):
    r = k.translate(COMP)[::-1]
    return k if k <= r else r


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_baits(panel_fa):
    """Return {bait_id: {canonical 31-mer set}} and reverse index
    {canonical kmer: set(bait_ids)}."""
    baits, name = {}, None
    for line in open(panel_fa):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            baits[name] = []
        elif line and name:
            baits[name].append(line.upper())
    baits = {k: "".join(v) for k, v in baits.items()}
    sets = {}
    for bid, seq in baits.items():
        ks = set()
        for i in range(len(seq) - 30):
            k = seq[i:i + 31]
            if set(k) <= {"A", "C", "G", "T"}:
                ks.add(canon(k))
        sets[bid] = ks
    index = {}
    for bid, ks in sets.items():
        for k in ks:
            index.setdefault(k, []).append(bid)
    return sets, index


def http_head(url, timeout=60):
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, {k: r.headers.get(k) for k in
                              ("Content-Length", "ETag", "Last-Modified")}
    except urllib.error.HTTPError as e:
        return e.code, {}
    except Exception as e:  # noqa: BLE001
        return -1, {"err": str(e)}


def download(url, dest, log):
    t0 = time.time()
    last = None
    for attempt in range(4):
        try:
            tmp = dest + ".part"
            with urllib.request.urlopen(urllib.request.Request(url), timeout=900) as r, \
                    open(tmp, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            os.replace(tmp, dest)
            size = os.path.getsize(dest)
            log({"event": "download", "url": url, "dest": dest, "size_bytes": size,
                 "elapsed_s": round(time.time() - t0, 2), "attempts": attempt + 1})
            return size, True
        except Exception as e:  # noqa: BLE001
            last = e
            try:
                os.remove(dest + ".part")
            except OSError:
                pass
            if attempt < 3:
                time.sleep((5, 20, 60)[attempt])
    log({"event": "download_error", "url": url, "err": str(last)[:300]})
    return 0, False


def screen_one(acc, zst, b2s, scratch_dir, sets, index, log):
    """Run b2s and return {bait: cov} (all baits) + n_found."""
    kmers_out = os.path.join(scratch_dir, f"{acc}.kmers.txt")
    t0 = time.time()
    cmd = [b2s, "--in-kmers", PANEL_PATH, "--in-sequences", zst,
           "--out-kmers", kmers_out, "--output-kmer-positions"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log({"event": "b2s_error", "accession": acc,
             "err": (r.stderr or r.stdout)[-300:]})
        return None
    found = {}
    with open(kmers_out) as fh:
        for line in fh:
            # found kmers carry "(read_id,pos,strand)" annotations; zero-count
            # kmers are written bare (verified on DRR000016 probe + D29/L5 refs)
            if "(" in line:
                k = canon(line.split()[0])
                bids = index.get(k)
                if bids:
                    for b in bids:
                        found[b] = found.get(b, 0) + 1
    try:
        os.remove(kmers_out)
    except OSError:
        pass
    cov = {b: round(found.get(b, 0) / len(ks), 4) for b, ks in sets.items()}
    n_found_total = len(found)
    log({"event": "screen", "accession": acc, "elapsed_s": round(time.time() - t0, 2),
         "n_baits_hit": sum(1 for v in cov.values() if v > 0)})
    return cov, n_found_total


PANEL_PATH = None  # set in main (b2s needs a file path)


def main():
    global PANEL_PATH
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-tsv", required=True,
                    help="TSV uid/run/scientific_name from fetch_runinfo.py")
    ap.add_argument("--baits", required=True)
    ap.add_argument("--run-root", required=True)
    ap.add_argument("--tier", required=True)
    ap.add_argument("--b2s", required=True)
    ap.add_argument("--concurrency", type=int, default=4)
    ap.add_argument("--max-bytes", type=int, default=400 * 1024**3)
    ap.add_argument("--max-accessions", type=int, default=25000)
    ap.add_argument("--deadline", required=True, help="ISO ts, e.g. 2026-08-21T14:08:52Z")
    args = ap.parse_args()

    root = args.run_root
    dl_dir = os.path.join(root, "downloads")
    scratch = os.path.join(root, "screen", "scratch")
    results_path = os.path.join(root, "screen", f"{args.tier}_results.jsonl")
    ledger_path = os.path.join(dl_dir, "download_ledger.jsonl")
    avail_path = os.path.join(root, "downloads", "availability.tsv")
    PANEL_PATH = os.path.abspath(args.baits)
    for d in (dl_dir, scratch, os.path.dirname(results_path)):
        os.makedirs(d, exist_ok=True)

    deadline = time.mktime(time.strptime(args.deadline, "%Y-%m-%dT%H:%M:%SZ")) - time.timezone
    log_lock = threading.Lock()

    def log(ev):
        ev["ts"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with log_lock:
            with open(ledger_path, "a") as fh:
                fh.write(json.dumps(ev, sort_keys=True) + "\n")

    sets, index = load_baits(PANEL_PATH)
    n_baits = len(sets)
    print(f"panel: {n_baits} baits, {len(index)} distinct canonical 31-mers")

    # run accessions (dedup, keep order); skip already-screened
    import csv as _csv
    runs = {}
    for r in _csv.DictReader(open(args.runs_tsv), delimiter="\t"):
        if r.get("run"):
            runs.setdefault(r["run"], r.get("scientific_name", ""))
    done = set()
    if os.path.exists(results_path):
        for line in open(results_path):
            try:
                done.add(json.loads(line)["accession"])
            except Exception:  # noqa: BLE001
                pass
    todo = [a for a in runs if a not in done]
    print(f"{len(runs)} runs ({len(done)} already screened), {len(todo)} to do")

    budget = {"bytes": 0, "accessions": 0}
    budget_lock = threading.Lock()
    stopped = {"reason": None}

    avail_cols = ["accession", "status", "size_bytes", "etag", "last_modified",
                  "sha256", "scientific_name", "tier"]
    avail_exists = os.path.exists(avail_path)
    avail_fh = open(avail_path, "a", newline="")
    import csv as _csv2
    avail_w = _csv2.DictWriter(avail_fh, delimiter="\t", fieldnames=avail_cols,
                               lineterminator="\n")
    if not avail_exists:
        avail_w.writeheader()

    def process(acc):
        if stopped["reason"]:
            return None
        url = S3_BASE.format(acc=acc)
        zst = os.path.join(dl_dir, f"{acc}.contigs.fa.zst")
        have = os.path.exists(zst)
        if not have:
            status, meta = http_head(url)
            with log_lock:
                avail_w.writerow({"accession": acc, "status": status,
                                  "size_bytes": meta.get("Content-Length", ""),
                                  "etag": meta.get("ETag", ""),
                                  "last_modified": meta.get("Last-Modified", ""),
                                  "sha256": "", "scientific_name": runs[acc],
                                  "tier": args.tier})
                avail_fh.flush()
            if status != 200:
                return None
            size = int(meta.get("Content-Length") or 0)
            with budget_lock:
                if budget["bytes"] + size > args.max_bytes:
                    stopped["reason"] = f"byte budget {args.max_bytes}"
                    return None
                budget["bytes"] += size
                budget["accessions"] += 1
                if budget["accessions"] > args.max_accessions:
                    stopped["reason"] = f"accession budget {args.max_accessions}"
                    return None
            got, ok = download(url, zst, log)
            if not ok:
                with budget_lock:
                    budget["bytes"] -= size
                return None
            sha = sha256_file(zst)
            with log_lock:
                avail_w.writerow({"accession": acc, "status": 200,
                                  "size_bytes": got, "etag": meta.get("ETag", ""),
                                  "last_modified": meta.get("Last-Modified", ""),
                                  "sha256": sha, "scientific_name": runs[acc],
                                  "tier": args.tier})
                avail_fh.flush()
            log({"event": "sha256", "accession": acc, "sha256": sha,
                 "size_bytes": got})
        res = screen_one(acc, zst, args.b2s, scratch, sets, index, log)
        if res is None:
            return None
        cov, n_found = res
        row = {"accession": acc, "tier": args.tier,
               "scientific_name": runs[acc],
               "n_found_panel_kmers": n_found,
               "per_bait_kmer_cov": cov,
               "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        with log_lock:
            with open(results_path, "a") as fh:
                fh.write(json.dumps(row, sort_keys=True) + "\n")
        hits = [b for b, v in cov.items() if v >= 0.7]
        return acc, hits

    t0 = time.time()
    n_ok = n_avail = n_hits = 0
    with ThreadPoolExecutor(max_workers=args.concurrency) as ex:
        futs = {ex.submit(process, a): a for a in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                out = fut.result()
            except Exception as e:  # noqa: BLE001
                log({"event": "worker_error", "accession": futs[fut],
                     "err": str(e)[:300]})
                out = None
            if out is not None:
                n_ok += 1
                if out[1]:
                    n_hits += 1
                    print(f"  HIT {out[0]}: {len(out[1])} baits >= 0.7 "
                          f"{sorted(out[1])[:4]}", flush=True)
            if stopped["reason"]:
                print(f"STOP: {stopped['reason']}")
                break
            if i % 50 == 0:
                el = time.time() - t0
                print(f"  {i}/{len(todo)} processed ({n_ok} screened, "
                      f"{n_hits} hit-accessions) {el:.0f}s "
                      f"budget={budget['bytes']/1e9:.1f}GB", flush=True)
            if time.time() > deadline:
                stopped["reason"] = "wall clock deadline"
                print("STOP: deadline")
                break
    avail_fh.close()
    print(json.dumps({"processed": n_ok, "hit_accessions": n_hits,
                      "budget": budget, "stopped": stopped["reason"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
