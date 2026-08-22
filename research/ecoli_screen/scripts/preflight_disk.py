#!/usr/bin/env python3
"""Tier B disk preflight: stratified S3 HEAD sample + budget projection.

Samples N accessions uniformly (deterministic stride) from a frozen tier
TSV, records size/404 ratio, and extrapolates the projected download volume
against RUN_PLAN §7 budgets. Writes results/<tier>_preflight.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
import statistics
import sys
import time
import urllib.error
import urllib.request


def head(url, timeout=45):
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, int(r.headers.get("Content-Length") or 0)
    except urllib.error.HTTPError as e:
        return e.code, 0
    except Exception as e:  # noqa: BLE001
        return -1, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen", required=True)
    ap.add_argument("--tier", required=True)
    ap.add_argument("--sample", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260819)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-bytes", type=float, default=400 * 1024**3)
    args = ap.parse_args()

    runs = [r["run"] for r in csv.DictReader(open(args.frozen), delimiter="\t")]
    rng = random.Random(args.seed)
    sample = sorted(rng.sample(runs, min(args.sample, len(runs))))

    sizes, miss404, err = [], 0, 0
    for i, acc in enumerate(sample):
        s, n = head(f"https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst")
        if s == 200:
            sizes.append(n)
        elif s == 404:
            miss404 += 1
        else:
            err += 1
        if i % 50 == 49:
            print(f"  {i+1}/{len(sample)} ok={len(sizes)} 404={miss404} err={err}", flush=True)
        time.sleep(0.05)

    frac_ok = len(sizes) / len(sample) if sample else 0
    mean = statistics.mean(sizes) if sizes else 0
    projected_total = mean * frac_ok * len(runs)
    sizes_sorted = sorted(sizes)
    res = {
        "tier": args.tier, "n_runs_frozen": len(runs),
        "n_sampled": len(sample), "seed": args.seed,
        "n_ok": len(sizes), "n_404": miss404, "n_err": err,
        "frac_in_logan": round(frac_ok, 4),
        "size_mb": {"min": round(sizes_sorted[0] / 1e6, 3),
                    "median": round(statistics.median(sizes) / 1e6, 3) if sizes else None,
                    "mean": round(mean / 1e6, 3), "p90": round(sizes_sorted[int(.9 * len(sizes))] / 1e6, 3) if sizes else None,
                    "max": round(sizes_sorted[-1] / 1e6, 3)},
        "projected_total_gb": round(projected_total / 1e9, 1),
        "max_bytes_gb": round(args.max_bytes / 1e9, 1),
        "within_budget": projected_total < args.max_bytes,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, f"{args.tier}_preflight.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res))
    return 0 if res["within_budget"] else 1


if __name__ == "__main__":
    sys.exit(main())
