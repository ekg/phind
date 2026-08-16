#!/usr/bin/env python3
"""ENA probe at scale: do collaborator run-accessions (ERR/SRR/DRR) have
retrievable submitted assemblies in ENA?

For each sampled run accession:
  1. run -> sample pivot : result=read_run, fields accession,sample_accession
  2. sample -> analysis  : result=analysis, fields analysis_accession,submitted_ftp

A "hit" = at least one analysis record with a non-empty submitted_ftp
(meaning a retrievable FASTA exists on ENA's FTP mirrors).

Deterministic sampling (fixed seed) so the probe is reproducible.

Usage:
  python ntm/v2/scripts/probe_ena_runs.py \
      --manifest /mnt/.../NTM_QC_passed_prophage_master_manifest.tsv \
      --n 100 --seed 42 \
      --out ntm/v2/run_assemblies/ena_probe_results.tsv
"""
import argparse
import csv
import json
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "https://www.ebi.ac.uk/ena/portal/api/search"
UA = "Mozilla/5.0 (phind NTM v2 run-assembly ENA probe; erik@hypervolu.me)"
TIMEOUT = 30


def api_search(result, query, fields):
    params = {
        "result": result,
        "query": query,
        "fields": fields,
        "format": "tsv",
        "limit": "1000",
    }
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    rows = []
    for line in raw.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        rows.append(parts)
    return rows


def probe_run(acc, delay):
    """Return dict with pivot + analysis results for one run accession."""
    time.sleep(delay)
    rec = {
        "run": acc,
        "sample_accession": "",
        "pivot_ok": False,
        "n_analysis": 0,
        "n_hits": 0,
        "analysis_accessions": "",
        "submitted_ftp": "",
        "error": "",
    }
    try:
        pivot = api_search("read_run", f"accession={acc}",
                           "accession,sample_accession")
        rec["pivot_ok"] = bool(pivot and len(pivot) > 1)
        if not pivot or len(pivot) < 2:
            rec["error"] = "no read_run record"
            return rec
        # header + rows; take first data row
        sample = pivot[1][1] if len(pivot[1]) > 1 else ""
        rec["sample_accession"] = sample
        if not sample:
            rec["error"] = "no sample_accession"
            return rec
        anl = api_search("analysis", f"sample_accession={sample}",
                         "analysis_accession,submitted_ftp")
        data_rows = anl[1:] if anl and anl[0][0] == "analysis_accession" else anl
        rec["n_analysis"] = len(data_rows)
        hits = [r for r in data_rows if len(r) > 1 and r[1].strip()]
        rec["n_hits"] = len(hits)
        rec["analysis_accessions"] = ";".join(r[0] for r in data_rows)
        rec["submitted_ftp"] = ";".join(r[1] for r in data_rows)
        return rec
    except urllib.error.HTTPError as e:
        rec["error"] = f"HTTP {e.code}"
        return rec
    except Exception as e:  # noqa: BLE001 - probe must not die on one run
        rec["error"] = f"{type(e).__name__}: {e}"
        return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="NTM_QC_passed_prophage_master_manifest.tsv path")
    ap.add_argument("--n", type=int, default=100, help="sample size")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", required=True, help="TSV output path")
    ap.add_argument("--summary", help="JSON summary output path")
    ap.add_argument("--save-sample", help="save sampled accession list to this path")
    args = ap.parse_args()

    # unique run accessions from the manifest
    by_acc = {}
    with open(args.manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = r["accession"]
            if acc[:3] not in ("ERR", "SRR", "DRR"):
                continue
            if acc not in by_acc:
                by_acc[acc] = {
                    "has_prophage": r["has_prophage"] == "True",
                    "assembly_contigs": r["assembly_contigs"],
                    "assembly_length_bp": r["assembly_length_bp"],
                }
            elif r["has_prophage"] == "True":
                by_acc[acc]["has_prophage"] = True

    all_accs = sorted(by_acc)
    rng = random.Random(args.seed)
    sampled = rng.sample(all_accs, min(args.n, len(all_accs)))

    if args.save_sample:
        with open(args.save_sample, "w") as fh:
            for a in sampled:
                fh.write(a + "\n")

    delay = 0.35  # base delay per request; ENA guidance ~15 req/s max
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(probe_run, acc, delay * (i % args.workers)):
                acc for i, acc in enumerate(sampled)}
        done = 0
        for fut in as_completed(futs):
            acc = futs[fut]
            rec = fut.result()
            rec["has_prophage"] = by_acc[acc]["has_prophage"]
            rec["assembly_contigs"] = by_acc[acc]["assembly_contigs"]
            rec["assembly_length_bp"] = by_acc[acc]["assembly_length_bp"]
            results.append(rec)
            done += 1
            sys.stderr.write(f"\rprobed {done}/{len(sampled)} {acc} "
                             f"hits={rec['n_hits']}  ")
            sys.stderr.flush()
    sys.stderr.write("\n")

    results.sort(key=lambda r: r["run"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(results[0].keys()),
                           lineterminator="\n")
        w.writeheader()
        w.writerows(results)

    n_hit = sum(1 for r in results if r["n_hits"] > 0)
    n_analysis = sum(1 for r in results if r["n_analysis"] > 0)
    summary = {
        "n_probed": len(results),
        "n_pivot_ok": sum(1 for r in results if r["pivot_ok"]),
        "n_with_any_analysis": n_analysis,
        "n_with_retrievable_ftp": n_hit,
        "hit_rate_pct": round(100.0 * n_hit / len(results), 2),
        "hit_runs": [r["run"] for r in results if r["n_hits"] > 0],
        "errors": [{"run": r["run"], "error": r["error"]}
                   for r in results if r["error"]],
    }
    if args.summary:
        with open(args.summary, "w") as fh:
            json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"details -> {args.out}")


if __name__ == "__main__":
    main()
