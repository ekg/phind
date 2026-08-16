#!/usr/bin/env python3
"""
NTM v2 — validation + report generation for the PanSN download.

Reads the downloader's progress.json plus the on-disk canonical_objects and
emits:

  {output_dir}/gcf_only.tsv        GCF-only accessions (one per line)
  {output_dir}/failed_accessions.txt   accession<TAB>reason (if any)
  {output_dir}/download_report.md  counts, total bp, failures

Validation performed (mirrors the task's ## Validation section):
  * resolved == downloaded + linked + failed == cohort total
  * samtools faidx re-runs cleanly on every downloaded genome (.fai valid)
  * gzip -t passes on >= 200 random files (or all if fewer)
  * spot-check 3 genomes: .fai contig names vs manifest prophage_contig
    (NZ_ stripped on both sides for GCF-only accessions)

Usage:
  python3 validate_download.py [--data-dir .../ntm/v2/inputs] \
      [--output-dir .../ntm/v2/genomes] [--v1-dir .../ntm/v1]
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import os
import random
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs")
    ap.add_argument("--output-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/genomes")
    ap.add_argument("--v1-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1")
    ap.add_argument("--report-path", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/download_report.md")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--random-check", type=int, default=200)
    ap.add_argument("--spotcheck", type=int, default=3)
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.output_dir)
    co_dir = out_dir / "canonical_objects"
    progress = json.loads((out_dir / "progress.json").read_text())

    completed = progress.get("completed", {})
    failed = progress.get("failed", {})
    total = progress.get("total", 0)

    downloaded = {a: v for a, v in completed.items() if v.get("status") == "COMPLETE"}
    linked = {a: v for a, v in completed.items() if v.get("status") == "LINKED"}
    already = {a: v for a, v in completed.items() if v.get("status") == "ALREADY_COMPLETE"}
    # ALREADY_COMPLETE files are on disk; count them as downloaded-or-linked by
    # whether their source was v1 (link) — we infer by .pansn.fa.gz is symlink.
    already_linked = {}
    already_downloaded = {}
    for a, v in already.items():
        p = co_dir / a / f"{a}.pansn.fa.gz"
        if p.is_symlink():
            already_linked[a] = v
        else:
            already_downloaded[a] = v

    n_downloaded = len(downloaded) + len(already_downloaded)
    n_linked = len(linked) + len(already_linked)
    n_failed = len(failed)
    resolved = n_downloaded + n_linked + n_failed

    print(f"cohort total={total} downloaded={n_downloaded} linked={n_linked} "
          f"failed={n_failed} resolved={resolved}")

    # --- fail if unresolved -------------------------------------------------
    checks: list[tuple[str, bool]] = []
    checks.append(("resolved == cohort total", resolved == total))

    # --- samtools faidx re-index on every downloaded genome ------------------
    fai_failures: list[str] = []
    def check_faidx(acc: str) -> tuple[str, str]:
        p = co_dir / acc / f"{acc}.pansn.fa.gz"
        try:
            r = subprocess.run(["samtools", "faidx", str(p)],
                               capture_output=True, timeout=300)
            if r.returncode != 0:
                return acc, r.stderr.decode(errors="replace")[:200]
            fai = Path(str(p) + ".fai")
            if not fai.exists() or fai.stat().st_size == 0:
                return acc, "fai missing/empty"
            return acc, ""
        except Exception as e:  # noqa: BLE001
            return acc, str(e)

    all_downloaded = list(downloaded) + list(already_downloaded)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for acc, err in ex.map(check_faidx, all_downloaded):
            if err:
                fai_failures.append(f"{acc}\t{err}")
    checks.append(("samtools faidx ok on all downloaded", not fai_failures))
    print(f"  faidx failures: {len(fai_failures)}")

    # --- gzip -t on random sample --------------------------------------------
    all_files = [co_dir / a / f"{a}.pansn.fa.gz" for a in list(completed)]
    sample = random.sample(all_files, min(args.random_check, len(all_files)))
    gzip_failures: list[str] = []
    for p in sample:
        try:
            with gzip.open(p, "rb") as fh:
                fh.read(1024)
            r = subprocess.run(["gzip", "-t", str(p)], capture_output=True)
            if r.returncode != 0:
                gzip_failures.append(str(p))
        except Exception as e:  # noqa: BLE001
            gzip_failures.append(f"{p}\t{e}")
    checks.append(("gzip -t passes on sampled files", not gzip_failures))
    print(f"  gzip -t failures on {len(sample)} sampled: {len(gzip_failures)}")

    # --- spot-check 3 genomes vs manifest prophage_contig ---------------------
    manifest_rows: dict[str, list[dict[str, str]]] = {}
    with open(data_dir / "NTM_QC_passed_prophage_master_manifest.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("has_prophage") == "True":
                manifest_rows.setdefault(row["accession"], []).append(row)

    # prefer accessions that have prophage rows and are in our completed set
    candidates = [a for a in completed if a in manifest_rows]
    spot = candidates[: args.spotcheck]
    spot_results: list[str] = []
    spot_ok = True
    for acc in spot:
        fai_path = co_dir / acc / f"{acc}.pansn.fa.gz.fai"
        fai_contigs = set()
        if fai_path.exists():
            for line in open(fai_path):
                fai_contigs.add(line.split("\t")[0].split("#")[-1])
        for row in manifest_rows[acc][:3]:
            contig = row["prophage_contig"]
            norm = contig[3:] if contig.startswith("NZ_") else contig  # NZ_ strip
            hit = norm in fai_contigs
            if not hit:
                spot_ok = False
            spot_results.append(
                f"- {acc}: manifest prophage_contig {contig} -> "
                f"normalized {norm} {'MATCH' if hit else 'NO MATCH'} "
                f"(fai has {len(fai_contigs)} contigs)"
            )
    checks.append(("spot-check contig names match manifest", spot_ok))
    print(f"  spot-checked: {spot}")

    # --- total bp --------------------------------------------------------------
    total_bp = 0
    per_acc_bp: dict[str, int] = {}
    for a, v in completed.items():
        tb = v.get("total_bases")
        if tb is None:
            fai = co_dir / a / f"{a}.pansn.fa.gz.fai"
            if fai.exists():
                tb = sum(int(l.split("\t")[1]) for l in open(fai))
        if tb:
            total_bp += tb
            per_acc_bp[a] = tb

    # --- write gcf_only.tsv -----------------------------------------------------
    gcf_only = [a for a, v in completed.items() if v.get("source") == "GCF"]
    gcf_only += [a for a, v in failed.items() if v.get("source") == "GCF"]
    gcf_only_path = out_dir / "gcf_only.tsv"
    with open(gcf_only_path, "w") as fh:
        fh.write("accession\n")
        for a in sorted(set(gcf_only)):
            fh.write(f"{a}\n")
    print(f"gcf_only.tsv: {len(set(gcf_only))} accessions")

    # --- write failed_accessions.txt ----------------------------------------------
    failed_path = out_dir / "failed_accessions.txt"
    if failed:
        with open(failed_path, "w") as fh:
            for a in sorted(failed):
                err = failed[a].get("error", "unknown")
                fh.write(f"{a}\t{err}\n")
    else:
        failed_path.write_text("")

    # --- write report -------------------------------------------------------------
    lines: list[str] = []
    lines.append("# NTM v2 — PanSN genome download report")
    lines.append("")
    lines.append(f"Generated: {utcnow()}")
    lines.append("")
    lines.append("## Counts")
    lines.append("")
    lines.append(f"| metric | count |")
    lines.append(f"|---|---|")
    lines.append(f"| cohort total (unique assemblies) | {total} |")
    lines.append(f"| downloaded | {n_downloaded} |")
    lines.append(f"| linked from v1 | {n_linked} |")
    lines.append(f"| failed | {n_failed} |")
    lines.append(f"| resolved (downloaded+linked+failed) | {resolved} |")
    lines.append(f"| total bp | {total_bp:,} |")
    lines.append(f"| GCF-only accessions (gcf_only.tsv) | {len(set(gcf_only))} |")
    lines.append("")
    lines.append("## Validation")
    lines.append("")
    for label, ok in checks:
        lines.append(f"- [{'x' if ok else ' '}] {label}")
    lines.append("")
    lines.append("## Failures")
    lines.append("")
    if failed:
        lines.append("| accession | reason |")
        lines.append("|---|---|")
        for a in sorted(failed):
            lines.append(f"| {a} | {failed[a].get('error', 'unknown')} |")
    else:
        lines.append("None.")
    lines.append("")
    if fai_failures:
        lines.append("## faidx failures")
        lines.append("")
        for f in fai_failures:
            lines.append(f"- {f}")
        lines.append("")
    if gzip_failures:
        lines.append("## gzip -t failures")
        lines.append("")
        for f in gzip_failures:
            lines.append(f"- {f}")
        lines.append("")
    lines.append("## Spot-check (contig names vs manifest prophage_contig)")
    lines.append("")
    lines.extend(spot_results)
    lines.append("")
    lines.append("## GCF-only accessions (recorded in gcf_only.tsv)")
    lines.append("")
    lines.append(f"{len(set(gcf_only))} accessions — see `gcf_only.tsv`.")

    report_path = Path(args.report_path)
    report_path.write_text("\n".join(lines))
    print(f"report written: {report_path}")

    # also write a copy in the output dir
    (out_dir / "download_report.md").write_text("\n".join(lines))

    all_ok = all(ok for _, ok in checks)
    print("ALL CHECKS PASS" if all_ok else "CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
