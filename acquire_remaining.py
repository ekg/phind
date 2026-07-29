#!/usr/bin/env python3
"""
Acquire remaining ~25k E. coli genomes for the full 26k prophage set.

Downloads from NCBI Datasets v2 API in parallel, converts to PanSN bgzip format.

Usage:
  python acquire_remaining.py --accessions need_to_download.txt \\
      --workers 10 --output-dir /mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/acquire-remaining-25k

Output structure:
  {output_dir}/
    canonical_objects/{accession}/{accession}.pansn.fa.gz
    canonical_objects/{accession}/{accession}.pansn.fa.gz.fai
    canonical_objects/{accession}/{accession}.pansn.fa.gz.gzi
    progress.json
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class RateLimiter:
    """Simple token-bucket rate limiter for API calls."""

    def __init__(self, calls_per_second: float = 2.0):
        self.min_interval = 1.0 / calls_per_second
        self._lock = threading.Lock()
        self._last_call = 0.0

    def acquire(self) -> None:
        """Block until the next allowed call time."""
        with self._lock:
            now = time.time()
            elapsed = now - self._last_call
            if elapsed < self.min_interval:
                sleep_time = self.min_interval - elapsed
                time.sleep(sleep_time)
            self._last_call = time.time()

DOWNLOAD_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
USER_AGENT = "phind-acquire-remaining-25k/1.0"
INCLUDE_TYPES = ("GENOME_FASTA",)
ACCESSION_RE = re.compile(r"^GC[AF]_[0-9]{9}\.[1-9][0-9]*$")
CONTIG_HEADER_RE = re.compile(r"^>(\S+)")
RETRY_DELAY = 2.0
MAX_RETRIES = 3


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_url(accession: str) -> str:
    """Build the NCBI Datasets v2 download URL for a single accession."""
    query = urllib.parse.urlencode(
        {"include_annotation_type": list(INCLUDE_TYPES), "filename": f"{accession}.zip"},
        doseq=True,
    )
    return DOWNLOAD_BASE.format(accession=accession) + "?" + query


def download_zip(accession: str, dst: Path, rate_limiter: RateLimiter | None = None) -> dict[str, Any]:
    """Download a single accession ZIP from NCBI Datasets v2 API.

    Returns download metadata including HTTP response headers and local SHA256.
    """
    url = download_url(accession)
    headers = {"User-Agent": USER_AGENT, "Accept": "application/zip"}
    attempt = 0
    last_error = None

    while attempt < MAX_RETRIES:
        attempt += 1
        try:
            # Rate limit before making the request
            if rate_limiter is not None:
                rate_limiter.acquire()

            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=300) as response:
                data = response.read()
                info = {
                    "status": response.status,
                    "url": response.geturl(),
                    "content_length": response.headers.get("Content-Length", "."),
                    "etag": response.headers.get("ETag", "."),
                    "last_modified": response.headers.get("Last-Modified", "."),
                }

            # Handle rate limiting (HTTP 429)
            if info["status"] == 429:
                last_error = f"HTTP 429 Rate Limited"
                retry_after = response.headers.get("Retry-After", "5")
                try:
                    sleep_time = float(retry_after)
                except ValueError:
                    sleep_time = 5.0
                time.sleep(sleep_time)
                continue

            # Write to temporary path, then atomically rename
            tmp = dst.with_suffix(".zip.tmp")
            tmp.write_bytes(data)
            os.rename(tmp, dst)

            info["local_sha256"] = sha_file(dst)
            info["local_bytes"] = dst.stat().st_size
            info["accession"] = accession
            info["attempts"] = attempt
            info["downloaded_at_utc"] = utcnow()
            return info

        except urllib.error.HTTPError as e:
            if e.code == 429:
                last_error = f"HTTP 429 Rate Limited"
                retry_after = e.headers.get("Retry-After", "5")
                try:
                    sleep_time = float(retry_after)
                except ValueError:
                    sleep_time = 5.0
                time.sleep(sleep_time)
                continue
            last_error = f"HTTP {e.code}: {e.reason}"
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY * attempt)
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last_error = str(e)
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY * attempt)

    raise RuntimeError(f"Failed to download {accession} after {MAX_RETRIES} attempts: {last_error}")


def extract_fasta_from_zip(zip_path: Path) -> tuple[bytes, str, str]:
    """Extract the genomic FASTA from a Datasets ZIP.

    Returns (uncompressed_fasta_bytes, fasta_member_path, accession_assembly_name).
    The FASTA member is typically like:
      ncbi_dataset/data/GCF_XXXXX.Y/GCF_XXXXX.Y_ASMXXX_genomic.fna
    """
    with zipfile.ZipFile(zip_path, "r") as zf:
        # Find the genomic FASTA member
        fasta_member = None
        for name in zf.namelist():
            if name.endswith("_genomic.fna"):
                fasta_member = name
                break
        if fasta_member is None:
            raise ValueError(f"No genomic FASTA found in {zip_path}")

        fasta_bytes = zf.read(fasta_member)

        # Extract assembly name from the member path
        # e.g. ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna
        assembly_name = Path(fasta_member).stem
        # Remove trailing _genomic
        if assembly_name.endswith("_genomic"):
            assembly_name = assembly_name[:-8]

        return fasta_bytes, fasta_member, assembly_name


def rename_fasta_to_pansn(fasta_bytes: bytes, accession: str) -> bytes:
    """Rename FASTA contig headers to PanSN format.

    Original: >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
    PanSN:    >{accession}#1#{contig_accession_version}

    Also handles alternative contig accessions like NZ_* or unnamed contigs.
    """
    lines = fasta_bytes.split(b"\n")
    output = []
    for line in lines:
        if line.startswith(b">"):
            # Extract the first token as the contig accession
            header = line[1:].decode("utf-8", errors="replace")
            first_token = header.split()[0] if header.split() else "unknown"
            # Sanitize: remove any characters that might cause issues
            safe_contig = first_token.replace(" ", "_")
            pansn_name = f"{accession}#1#{safe_contig}"
            output.append(f">{pansn_name}".encode("utf-8"))
        else:
            output.append(line)
    return b"\n".join(output)


def convert_to_bgzip(fasta_bytes: bytes, pansn_fasta_path: Path, bgzip_bin: str = "bgzip", threads: int = 2) -> dict[str, Any]:
    """Convert plain FASTA bytes to bgzip-compressed PanSN FASTA.

    Writes to pansn_fasta_path, creates .fai and .gzi indexes.
    """
    # Write renamed FASTA to temp file
    # pansn_fasta_path is like: .../canonical_objects/GCF_X/GCF_X.pansn.fa.gz
    tmp_fa = pansn_fasta_path.with_name(pansn_fasta_path.name.replace(".pansn.fa.gz", ".fa.tmp"))
    tmp_fa.write_bytes(fasta_bytes)

    bgz_tmp = pansn_fasta_path.with_name(pansn_fasta_path.name + ".tmp")

    # Run bgzip
    result = subprocess.run(
        [bgzip_bin, "-c", "-@", str(threads), str(tmp_fa)],
        capture_output=True,
        check=True,
    )
    bgz_tmp.write_bytes(result.stdout)

    # Atomically rename
    os.rename(bgz_tmp, pansn_fasta_path)

    # Clean up temp FASTA
    tmp_fa.unlink(missing_ok=True)

    info = {
        "canonical_bgzf_bytes": pansn_fasta_path.stat().st_size,
        "canonical_bgzf_sha256": sha_file(pansn_fasta_path),
    }

    # Create .fai index using samtools faidx on the bgzip file
    try:
        result = subprocess.run(
            ["samtools", "faidx", str(pansn_fasta_path)],
            capture_output=True, check=True, timeout=120,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, TimeoutError):
        pass  # Indexing is optional

    fai_path = pansn_fasta_path.with_name(pansn_fasta_path.name + ".fai")
    if fai_path.exists():
        info["fai_sha256"] = sha_file(fai_path)
        info["fai_bytes"] = fai_path.stat().st_size
    else:
        info["fai_sha256"] = "."
        info["fai_bytes"] = 0

    # Create .gzi index (bgzip index)
    gzi_path = pansn_fasta_path.with_name(pansn_fasta_path.name + ".gzi")
    try:
        result = subprocess.run(
            [bgzip_bin, "-r", str(pansn_fasta_path)],
            capture_output=True, check=True, timeout=120,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, TimeoutError):
        pass

    # After bgzip -r, the .gzi is created alongside the .gz file
    actual_gzi = pansn_fasta_path.with_name(pansn_fasta_path.name + ".gzi")
    if actual_gzi.exists():
        info["gzi_sha256"] = sha_file(actual_gzi)
        info["gzi_bytes"] = actual_gzi.stat().st_size
    else:
        info["gzi_sha256"] = "."
        info["gzi_bytes"] = 0

    return info


def process_accession(
    accession: str,
    output_dir: Path,
    scratch_dir: Path,
    bgzip_bin: str,
    bgzip_threads: int,
    force: bool = False,
    rate_limiter: RateLimiter | None = None,
) -> dict[str, Any]:
    """Download and convert a single accession to PanSN bgzip format.

    Returns a dict with processing metadata.
    """
    result: dict[str, Any] = {
        "accession": accession,
        "status": "PENDING",
        "started_at_utc": utcnow(),
    }

    canonical_dir = output_dir / "canonical_objects" / accession
    canonical_dir.mkdir(parents=True, exist_ok=True)

    pansn_fa_gz = canonical_dir / f"{accession}.pansn.fa.gz"

    # Skip if already complete
    if pansn_fa_gz.exists() and not force:
        result["status"] = "ALREADY_COMPLETE"
        result["canonical_bgzf_bytes"] = pansn_fa_gz.stat().st_size
        result["canonical_bgzf_sha256"] = sha_file(pansn_fa_gz)
        result["completed_at_utc"] = utcnow()
        return result

    try:
        # Download ZIP to scratch
        zip_path = scratch_dir / f"{accession}.zip"
        download_info = download_zip(accession, zip_path, rate_limiter)

        # Extract FASTA from ZIP
        fasta_bytes, fasta_member, assembly_name = extract_fasta_from_zip(zip_path)

        # Contig count (count > headers)
        contig_count = fasta_bytes.count(b"\n>") + 1 if fasta_bytes.startswith(b">") else 0
        total_bases = sum(len(line) for line in fasta_bytes.split(b"\n") if line and not line.startswith(b">"))

        # Rename contigs to PanSN format
        pansn_fasta = rename_fasta_to_pansn(fasta_bytes, accession)

        # Write renamed FASTA to temp and bgzip compress
        conversion_info = convert_to_bgzip(pansn_fasta, pansn_fa_gz, bgzip_bin, bgzip_threads)

        # Clean up ZIP
        zip_path.unlink(missing_ok=True)

        result.update({
            "status": "COMPLETE",
            "source_package_bytes": download_info.get("local_bytes", 0),
            "source_package_sha256": download_info.get("local_sha256", "."),
            "source_fasta_member": fasta_member,
            "source_fasta_bytes": len(fasta_bytes),
            "contig_count": contig_count,
            "total_bases": total_bases,
            "download_attempts": download_info.get("attempts", 1),
            "download_etag": download_info.get("etag", "."),
            "download_status": download_info.get("status", 0),
        })
        result.update(conversion_info)
        result["completed_at_utc"] = utcnow()

    except Exception as e:
        result["status"] = "FAILED"
        result["error"] = str(e)
        result["completed_at_utc"] = utcnow()

    return result


def load_progress(progress_path: Path) -> dict[str, Any]:
    """Load progress state from file."""
    if progress_path.exists():
        try:
            return json.loads(progress_path.read_text())
        except (json.JSONDecodeError, OSError):
            pass
    return {"completed": {}, "failed": {}, "started_at": utcnow(), "total": 0, "done": 0}


def save_progress(progress_path: Path, state: dict[str, Any]) -> None:
    """Save progress state atomically."""
    tmp = progress_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True))
    os.rename(tmp, progress_path)


def main():
    parser = argparse.ArgumentParser(
        description="Acquire remaining E. coli genomes for full 26k prophage set"
    )
    parser.add_argument(
        "--accessions", required=True,
        help="File with one accession per line (e.g., need_to_download.txt)"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for canonical objects"
    )
    parser.add_argument(
        "--scratch-dir",
        default=None,
        help="Scratch directory for temp files (default: {output_dir}/.scratch)"
    )
    parser.add_argument(
        "--workers", type=int, default=10,
        help="Number of parallel download workers (default: 10)"
    )
    parser.add_argument(
        "--bgzip", default="bgzip",
        help="Path to bgzip binary (default: bgzip)"
    )
    parser.add_argument(
        "--bgzip-threads", type=int, default=2,
        help="Threads per bgzip call (default: 2)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download even if already complete"
    )
    parser.add_argument(
        "--batch-size", type=int, default=100,
        help="Progress save interval in accessions (default: 100)"
    )
    parser.add_argument(
        "--max-failures", type=int, default=100,
        help="Maximum allowed failures before aborting (default: 100)"
    )
    parser.add_argument(
        "--rate-limit", type=float, default=2.0,
        help="NCBI API calls per second (default: 2.0, max ~3.0)"
    )
    args = parser.parse_args()

    # Read accessions
    accessions = []
    with open(args.accessions) as f:
        for line in f:
            line = line.strip()
            if line and ACCESSION_RE.match(line):
                accessions.append(line)

    print(f"Loaded {len(accessions)} accessions to process", flush=True)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scratch_dir = Path(args.scratch_dir) if args.scratch_dir else output_dir / ".scratch"
    scratch_dir.mkdir(parents=True, exist_ok=True)

    progress_path = output_dir / "progress.json"
    state = load_progress(progress_path)

    if "total" not in state or state["total"] == 0:
        state["total"] = len(accessions)
        state["started_at"] = utcnow()
        state["completed"] = {}
        state["failed"] = {}
        state["done"] = 0
        save_progress(progress_path, state)

    # Filter out already completed
    to_process = []
    for acc in accessions:
        if acc in state["completed"] and not args.force:
            continue
        if acc in state["failed"] and not args.force:
            continue
        to_process.append(acc)

    print(f"Already completed: {len(accessions) - len(to_process)}", flush=True)
    print(f"Remaining to process: {len(to_process)}", flush=True)

    if not to_process:
        print("All accessions already processed!", flush=True)
        return

    # Process in parallel with rate limiting
    failure_count = len(state["failed"])
    # NCBI rate limit is ~3 req/s; use user-specified rate to be safe
    rate_limiter = RateLimiter(calls_per_second=args.rate_limit)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_acc = {
            executor.submit(
                process_accession,
                acc,
                output_dir,
                scratch_dir,
                args.bgzip,
                args.bgzip_threads,
                args.force,
                rate_limiter,
            ): acc
            for acc in to_process
        }

        batch_count = 0
        for future in concurrent.futures.as_completed(future_to_acc):
            acc = future_to_acc[future]
            try:
                result = future.result()
            except Exception as e:
                result = {
                    "accession": acc,
                    "status": "FAILED",
                    "error": str(e),
                    "completed_at_utc": utcnow(),
                }

            if result["status"] == "COMPLETE" or result["status"] == "ALREADY_COMPLETE":
                state["completed"][acc] = result
                state["done"] += 1
            else:
                state["failed"][acc] = result
                failure_count += 1

            # Print progress
            done = state["done"]
            failed = len(state["failed"])
            total = state["total"]
            elapsed = time.time() - 0  # Will be approximate
            pct = 100.0 * (done + failed) / total if total > 0 else 0
            print(
                f"[{done + failed}/{total} {pct:.1f}%] {acc}: {result['status']} "
                f"(failures: {failed})",
                flush=True,
            )

            # Save progress periodically
            batch_count += 1
            if batch_count % args.batch_size == 0:
                state["updated_at"] = utcnow()
                save_progress(progress_path, state)

            # Check failure threshold
            if failure_count >= args.max_failures:
                print(
                    f"ABORTING: reached {failure_count} failures (max {args.max_failures})",
                    flush=True,
                )
                executor.shutdown(wait=False, cancel_futures=True)
                break

    # Final save
    state["updated_at"] = utcnow()
    state["finished_at"] = utcnow()
    save_progress(progress_path, state)

    # Summary
    done = len(state["completed"])
    failed = len(state["failed"])
    total = state["total"]
    print(f"\n{'='*60}", flush=True)
    print(f"SUMMARY: {done + failed}/{total} processed", flush=True)
    print(f"  Completed: {done}", flush=True)
    print(f"  Failed:    {failed}", flush=True)

    if failed > 0:
        failed_file = output_dir / "failed_accessions.txt"
        with open(failed_file, "w") as f:
            for acc in state["failed"]:
                err = state["failed"][acc].get("error", "unknown")
                f.write(f"{acc}\t{err}\n")
        print(f"  Failed list: {failed_file}", flush=True)
        sys.exit(1 if failed > 0 else 0)


if __name__ == "__main__":
    main()