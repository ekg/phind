#!/usr/bin/env python3
"""
NTM v3 — S0b: build canonical PanSN bgzip genome objects for delivered SRA runs.

Scope (compute-determined, do not widen): exactly the 17,920 genomes whose
master-table `source == ASSEMBLY` (SRA runs ERR/SRR/DRR). Non-run genomes
already have canonical objects; this driver must never build or overwrite
anything else, and never deletes a pre-existing object.

Input  : {raw_dir}/{run}.fasta                      (bare-name FASTA headers)
Output : {out_dir}/{run}/{run}.pansn.fa.gz{,.fai,.gzi}

Contig renaming and bgzip/faidx follow acquire_v3_genomes.py exactly:
`contig_name_from_header` + `rename_fasta_to_pansn` + `write_bgzip_faidx`.

Resume-first: a run is reused (skipped) when {run}.pansn.fa.gz.fai exists and
is non-empty. Use --no-resume to rebuild regardless.

Usage:
  python3 ntm/v3/scripts/build_v3_pansn_objects.py [--jobs 16] [--limit N] \
      [--run DRR015955 --run DRR015956] [--dry-run] [--no-resume]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "ntm/v3/scripts"))
from acquire_v3_genomes import (  # noqa: E402
    contig_name_from_header,
    rename_fasta_to_pansn,
    write_bgzip_faidx,
)

NVME = Path("/mnt/nvme3n1/erikg/phind-genome-work")
DEFAULT_RAW_DIR = NVME / "ntm/v3/genomes/delivery_raw/sra_assembled"
DEFAULT_OUT_DIR = NVME / "ntm/v3/genomes/canonical_objects"
DEFAULT_MANIFEST = REPO / "NTM_Data/NTM_master_table_complete_20260922.tsv"


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def load_runs(manifest: Path) -> list[str]:
    """Canonical run list: genome_id of every `source == ASSEMBLY` row."""
    runs: list[str] = []
    seen: set[str] = set()
    with open(manifest, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("source") != "ASSEMBLY":
                continue
            run = row["genome_id"]
            if run and run not in seen:
                seen.add(run)
                runs.append(run)
    return sorted(runs)


def build_one(run: str, raw_dir: Path, out_dir: Path) -> dict:
    """Build a single run's object. Returns a status record (never raises)."""
    rec: dict = {"run": run, "status": "failed", "contigs": 0, "total_bp": 0,
                 "object_bytes": 0, "error": ""}
    src = raw_dir / f"{run}.fasta"
    obj_dir = out_dir / run
    fai = obj_dir / f"{run}.pansn.fa.gz.fai"
    try:
        if not src.exists():
            rec["error"] = f"missing input {src}"
            return rec
        fasta = src.read_bytes()
        pansn = rename_fasta_to_pansn(fasta, run)
        obj_dir.mkdir(parents=True, exist_ok=True)
        log: dict = {}
        stats = write_bgzip_faidx(pansn, obj_dir, run, log)
        if not fai.exists() or fai.stat().st_size == 0:
            rec["error"] = "faidx produced no .fai"
            return rec
        rec.update({"status": "built", "contigs": stats["contigs"],
                    "total_bp": stats["total_bp"],
                    "object_bytes": log.get("canonical_bgzf_bytes", 0)})
    except Exception as e:  # noqa: BLE001
        rec["error"] = str(e)[:500]
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--run", action="append", default=[],
                    help="restrict to this run id (repeatable)")
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--summary", type=Path, default=None,
                    help="summary path (default {out_dir}/../pansn_build_summary)")
    ap.add_argument("--no-resume", action="store_true",
                    help="rebuild even when a valid object already exists")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    summary_base = args.summary or (args.out_dir.parent / "pansn_build_summary")

    runs = load_runs(args.manifest)
    if args.run:
        wanted = set(args.run)
        runs = [r for r in runs if r in wanted]
        missing = sorted(wanted - set(runs))
        if missing:
            print(f"warning: --run ids not in ASSEMBLY manifest: {missing}",
                  file=sys.stderr)
    if args.limit:
        runs = runs[:args.limit]
    print(f"runs to consider: {len(runs)}", flush=True)

    todo: list[str] = []
    reused = 0
    for run in runs:
        fai = args.out_dir / run / f"{run}.pansn.fa.gz.fai"
        if not args.no_resume and fai.exists() and fai.stat().st_size > 0:
            reused += 1
            continue
        todo.append(run)

    if args.dry_run:
        n_missing = sum(1 for r in todo if not (args.raw_dir / f"{r}.fasta").exists())
        print(f"[dry-run] would build {len(todo)} (missing input: {n_missing}); "
              f"reused {reused}")
        return 0

    t0 = time.time()
    records: list[dict] = []
    built = failed = 0
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(build_one, r, args.raw_dir, args.out_dir): r
                for r in todo}
        done = 0
        for fut in as_completed(futs):
            rec = fut.result()
            records.append(rec)
            done += 1
            if rec["status"] == "built":
                built += 1
            else:
                failed += 1
                print(f"FAIL {rec['run']}: {rec['error']}", file=sys.stderr,
                      flush=True)
            if done % 200 == 0:
                print(f"  progress {done}/{len(todo)} built={built} "
                      f"failed={failed}", flush=True)

    records.sort(key=lambda r: r["run"])
    summary = {
        "generated_utc": utcnow(),
        "manifest": str(args.manifest),
        "raw_dir": str(args.raw_dir),
        "out_dir": str(args.out_dir),
        "runs_considered": len(runs),
        "built": built,
        "reused": reused,
        "failed": failed,
        "elapsed_s": round(time.time() - t0, 1),
        "records": records,
    }
    json_path = Path(str(summary_base) + ".json")
    tsv_path = Path(str(summary_base) + ".tsv")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=1))
    with open(tsv_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["run", "status", "contigs", "total_bp", "object_bytes", "error"])
        for r in records:
            w.writerow([r["run"], r["status"], r["contigs"], r["total_bp"],
                        r["object_bytes"], r["error"]])
    print(f"summary: {json_path} ; {tsv_path}")

    print(f"BUILD_DONE built={built} reused={reused} failed={failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
