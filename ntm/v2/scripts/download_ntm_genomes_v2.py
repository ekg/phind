#!/usr/bin/env python3
"""
NTM v2 — download missing NCBI assemblies as PanSN bgzip canonical objects.

Cohort: QC-passed NCBI assemblies from the v2 accession list, deduplicated by
numeric id, GCA-preferred (GenBank) with GCF fallback for GCF-only numerics.
Genomes already present in the v1 canonical_objects tree are linked (symlink)
instead of re-downloaded.

Layout (mirrors v1):
  {output_dir}/canonical_objects/{acc}/{acc}.pansn.fa.gz
  {output_dir}/canonical_objects/{acc}/{acc}.pansn.fa.gz.fai
  {output_dir}/canonical_objects/{acc}/{acc}.pansn.fa.gz.gzi

PanSN header format: >{acc}#1#{contig_name}, where contig_name is the first
token of the source FASTA header (bare WGS accession for GenBank fna, e.g.
CP000325.1; NZ_/NC_-prefixed for RefSeq fna). This matches the v1 download
path convention (see scripts/acquire_remaining.py: rename_fasta_to_pansn).

Download sources (in order):
  1. NCBI FTP direct: https://ftp.ncbi.nlm.nih.gov/genomes/all/{GCA|GCF}/{d1}/{d2}/{d3}/{acc}_{name}/{acc}_{name}_genomic.fna.gz
     (directory name taken from the manifest genome_id, which mirrors the
     ncbi_dataset fna filename). Canonical headers, fast.
  2. NCBI Datasets v2 API ZIP (single-accession download), as used in v1.
  3. ENA browser API fasta (gzip), headers parsed to bare contig accessions.

Usage:
  python3 download_ntm_genomes_v2.py \
      [--data-dir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs] \
      [--output-dir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/genomes] \
      [--v1-dir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v1] \
      [--workers 24] [--link-mode symlink|hardlink]
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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

ACCESSION_RE = re.compile(r"^GC[AF]_([0-9]{9})\.([0-9]+)$")
USER_AGENT = "phind-ntm-v2-download/1.0"
MAX_ATTEMPTS = 2  # retry once per task spec
RETRY_DELAY = 2.0
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
ENA_FASTA_BASE = "https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}?download=true&gzip=true"


class RateLimiter:
    """Token-bucket rate limiter for API calls."""

    def __init__(self, calls_per_second: float = 2.0):
        self.min_interval = 1.0 / calls_per_second
        self._lock = threading.Lock()
        self._last_call = 0.0

    def acquire(self) -> None:
        with self._lock:
            now = time.time()
            elapsed = now - self._last_call
            if elapsed < self.min_interval:
                time.sleep(self.min_interval - elapsed)
            self._last_call = time.time()


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def http_get_bytes(url: str, timeout: int = 300, rate_limiter: RateLimiter | None = None,
                 retries: int = 3) -> tuple[bytes, dict[str, Any]]:
    """GET a URL and return (body, info). Retries transient 429/5xx errors."""
    last_info: dict[str, Any] = {}
    for attempt in range(retries + 1):
        if rate_limiter is not None:
            rate_limiter.acquire()
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read()
                info = {
                    "status": resp.status,
                    "url": resp.geturl(),
                    "content_length": resp.headers.get("Content-Length", "."),
                    "etag": resp.headers.get("ETag", "."),
                }
            if info["status"] in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(RETRY_DELAY * (attempt + 1) * 5)
                continue
            return body, info
        except urllib.error.HTTPError as e:
            last_info = {"status": e.code, "error": str(e.reason)}
            if e.code in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(RETRY_DELAY * (attempt + 1) * 5)
                continue
            raise
    return b"", last_info


# ---------------------------------------------------------------------------
# cohort construction
# ---------------------------------------------------------------------------

def build_cohort(data_dir: Path) -> tuple[list[dict[str, Any]], dict[str, str]]:
    """Build the chosen-accession cohort.

    Returns (cohort, genome_id_by_acc) where cohort is a list of dicts with
    keys: numeric, accession, source (GCA|GCF), genome_id, and
    genome_id_by_acc maps accession -> genome_id (from the manifest).
    """
    acc_file = data_dir / "NTM_QC_passed_accession_list.tsv"
    man_file = data_dir / "NTM_QC_passed_prophage_master_manifest.tsv"

    # numeric id (digits + version) -> list of accessions from the accession list.
    # Grouping by the FULL numeric id (e.g. "000194015.2") keeps distinct GCA/GCF
    # versions as separate assemblies, matching the task cohort of 13,122.
    by_numeric: dict[str, list[str]] = {}
    with open(acc_file, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            m = ACCESSION_RE.match(acc)
            if m:
                by_numeric.setdefault(m.group(1) + "." + m.group(2), []).append(acc)

    # accession -> genome_id from manifest (first occurrence)
    genome_id_by_acc: dict[str, str] = {}
    with open(man_file, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            gid = (row.get("genome_id") or "").strip()
            if acc and gid and acc not in genome_id_by_acc:
                genome_id_by_acc[acc] = gid

    cohort: list[dict[str, Any]] = []
    for numeric, accs in sorted(by_numeric.items()):
        gca = sorted(a for a in accs if a.startswith("GCA_"))
        gcf = sorted(a for a in accs if a.startswith("GCF_"))
        if gca:
            chosen, source = gca[0], "GCA"
        else:
            chosen, source = gcf[0], "GCF"
        gid = genome_id_by_acc.get(chosen, "")
        m = ACCESSION_RE.match(chosen)
        cohort.append({
            "numeric": numeric,
            "accession": chosen,
            "source": source,
            "genome_id": gid,
            "digits": m.group(1) if m else "",
        })
    return cohort, genome_id_by_acc


# ---------------------------------------------------------------------------
# URL construction
# ---------------------------------------------------------------------------

def ftp_url(entry: dict[str, Any]) -> str | None:
    """NCBI FTP URL for the assembly fna.gz built from the manifest genome_id."""
    digits = entry.get("digits") or entry["accession"].split("_", 1)[1].split(".", 1)[0]
    gid = entry["genome_id"]
    acc = entry["accession"]
    if not digits or not gid:
        return None
    prefix = acc.split("_", 1)[0]  # GCA or GCF
    if not gid.endswith("_genomic"):
        return None
    dirname = gid[: -len("_genomic")]
    return f"{FTP_BASE}/{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{dirname}/{dirname}_genomic.fna.gz"


def datasets_zip_url(acc: str) -> str:
    query = urllib.parse.urlencode(
        {"include_annotation_type": ["GENOME_FASTA"], "filename": f"{acc}.zip"},
        doseq=True,
    )
    return DATASETS_BASE.format(acc=acc) + "?" + query


def ena_fasta_url(acc: str) -> str:
    return ENA_FASTA_BASE.format(acc=acc)


# ---------------------------------------------------------------------------
# source fetch -> raw FASTA bytes
# ---------------------------------------------------------------------------

def fetch_ftp(entry: dict[str, Any], rate_limiter: RateLimiter | None = None) -> tuple[bytes, str]:
    """Download the genomic fna.gz from NCBI FTP; returns (fasta_bytes, source)."""
    url = ftp_url(entry)
    if url is None:
        raise RuntimeError("no FTP URL (missing genome_id)")
    data, info = http_get_bytes(url, rate_limiter=rate_limiter)
    if info["status"] != 200:
        raise RuntimeError(f"HTTP {info['status']} from NCBI FTP")
    return gunzip_bytes(data), f"NCBI FTP {url}"


def gunzip_bytes(data: bytes) -> bytes:
    """Decompress gzip bytes; raise on invalid gzip."""
    if data[:2] != b"\x1f\x8b":
        raise RuntimeError("not a gzip stream")
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as fh:
        return fh.read()


def fetch_datasets(acc: str, rate_limiter: RateLimiter | None = None) -> tuple[bytes, str]:
    """Download single-accession ZIP from NCBI Datasets v2; returns (fasta_bytes, source)."""
    url = datasets_zip_url(acc)
    data, info = http_get_bytes(url, rate_limiter=rate_limiter)
    if info["status"] != 200:
        raise RuntimeError(f"HTTP {info['status']} from NCBI datasets")
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        fasta_member = None
        for name in zf.namelist():
            if name.endswith("_genomic.fna"):
                fasta_member = name
                break
        if fasta_member is None:
            raise RuntimeError(f"no _genomic.fna in datasets ZIP ({len(zf.namelist())} members)")
        return zf.read(fasta_member), f"NCBI datasets ZIP {url}"


def fetch_ena(acc: str, rate_limiter: RateLimiter | None = None) -> tuple[bytes, str]:
    """Download FASTA from ENA browser API (gzip); returns (fasta_bytes, source)."""
    url = ena_fasta_url(acc)
    data, info = http_get_bytes(url, rate_limiter=rate_limiter)
    if info["status"] != 200:
        raise RuntimeError(f"HTTP {info['status']} from ENA fasta API")
    return gunzip_bytes(data), f"ENA fasta API {url}"


def fetch_fasta(entry: dict[str, Any], rate_limiter: RateLimiter | None) -> tuple[bytes, str]:
    """Fetch raw FASTA bytes via the source chain (FTP -> datasets -> ENA)."""
    acc = entry["accession"]
    last_err: Exception | None = None
    sources = [("ftp", lambda: fetch_ftp(entry, rate_limiter)),
               ("datasets", lambda: fetch_datasets(acc, rate_limiter)),
               ("ena", lambda: fetch_ena(acc, rate_limiter))]
    for name, fn in sources:
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            last_err = e
    raise RuntimeError(f"all sources failed: {last_err}")


# ---------------------------------------------------------------------------
# PanSN conversion
# ---------------------------------------------------------------------------

def contig_name_from_header(header: str) -> str:
    """Extract the contig name from a FASTA header line (first token).

    NCBI fna: '>CP000325.1 Mycobacterium ...' -> 'CP000325.1'
    ENA API:  '>ENA|CP000325|CP000325.1 Mycobacterium ...' -> 'CP000325.1'
    """
    token = header.strip().split()[0] if header.strip().split() else "unknown"
    if "|" in token:
        token = token.split("|")[-1]
    return token.replace(" ", "_")


def rename_fasta_to_pansn(fasta_bytes: bytes, accession: str) -> bytes:
    """Rewrite FASTA headers to PanSN format: >{acc}#1#{contig_name}."""
    out = []
    for line in fasta_bytes.split(b"\n"):
        if line.startswith(b">"):
            header = line[1:].decode("utf-8", errors="replace")
            name = contig_name_from_header(header)
            out.append(f">{accession}#1#{name}".encode("utf-8"))
        else:
            out.append(line)
    return b"\n".join(out)


# ---------------------------------------------------------------------------
# bgzip + index
# ---------------------------------------------------------------------------

def convert_to_bgzip(fasta_bytes: bytes, pansn_path: Path, bgzip_bin: str = "bgzip", threads: int = 2) -> dict[str, Any]:
    """Compress renamed FASTA with bgzip; create .fai and .gzi indexes."""
    tmp_fa = pansn_path.with_name(pansn_path.name.replace(".pansn.fa.gz", ".fa.tmp"))
    tmp_fa.write_bytes(fasta_bytes)
    bgz_tmp = pansn_path.with_name(pansn_path.name + ".tmp")
    try:
        result = subprocess.run(
            [bgzip_bin, "-c", "-@", str(threads), str(tmp_fa)],
            capture_output=True, check=True,
        )
        bgz_tmp.write_bytes(result.stdout)
        os.rename(bgz_tmp, pansn_path)
    finally:
        tmp_fa.unlink(missing_ok=True)

    info: dict[str, Any] = {
        "canonical_bgzf_bytes": pansn_path.stat().st_size,
        "canonical_bgzf_sha256": sha_file(pansn_path),
    }

    # .fai via samtools faidx
    try:
        subprocess.run(["samtools", "faidx", str(pansn_path)], capture_output=True, check=True, timeout=180)
    except (subprocess.CalledProcessError, FileNotFoundError, TimeoutError) as e:
        info["fai_error"] = str(e)
    fai_path = pansn_path.with_name(pansn_path.name + ".fai")
    info["fai_exists"] = fai_path.exists()
    if fai_path.exists():
        info["fai_sha256"] = sha_file(fai_path)
        info["contig_count"] = sum(1 for _ in open(fai_path))
        info["total_bases"] = sum(int(line.split("\t")[1]) for line in open(fai_path))

    # .gzi via bgzip -r (reindex)
    try:
        subprocess.run([bgzip_bin, "-r", str(pansn_path)], capture_output=True, check=True, timeout=180)
    except (subprocess.CalledProcessError, FileNotFoundError, TimeoutError) as e:
        info["gzi_error"] = str(e)
    gzi_path = pansn_path.with_name(pansn_path.name + ".gzi")
    info["gzi_exists"] = gzi_path.exists()
    if gzi_path.exists():
        info["gzi_sha256"] = sha_file(gzi_path)

    return info


# ---------------------------------------------------------------------------
# per-accession processing
# ---------------------------------------------------------------------------

def process_accession(
    entry: dict[str, Any],
    output_dir: Path,
    scratch_dir: Path,
    link_mode: str,
    bgzip_bin: str,
    bgzip_threads: int,
    rate_limiter: RateLimiter | None = None,
    force: bool = False,
) -> dict[str, Any]:
    acc = entry["accession"]
    result: dict[str, Any] = {
        "accession": acc,
        "numeric": entry["numeric"],
        "source": entry["source"],
        "status": "PENDING",
        "started_at_utc": utcnow(),
    }

    canonical_dir = output_dir / "canonical_objects" / acc
    canonical_dir.mkdir(parents=True, exist_ok=True)
    pansn_path = canonical_dir / f"{acc}.pansn.fa.gz"

    if pansn_path.exists() and not force:
        result["status"] = "ALREADY_COMPLETE"
        result["completed_at_utc"] = utcnow()
        return result

    # 1) link from v1 if present
    v1_path = entry.get("v1_path")
    if v1_path is not None:
        for suffix in ("", ".fai", ".gzi"):
            src = Path(str(v1_path) + suffix)
            dst = Path(str(pansn_path) + suffix)
            if src.exists() and not dst.exists():
                if link_mode == "hardlink":
                    os.link(src, dst)
                else:
                    os.symlink(src, dst)
        result.update({
            "status": "LINKED",
            "v1_path": str(v1_path),
            "canonical_bgzf_bytes": pansn_path.stat().st_size if pansn_path.exists() else 0,
        })
        # total bp from .fai
        fai = Path(str(pansn_path) + ".fai")
        if fai.exists():
            result["total_bases"] = sum(int(l.split("\t")[1]) for l in open(fai))
            result["contig_count"] = sum(1 for _ in open(fai))
        result["completed_at_utc"] = utcnow()
        return result

    # 2) download + convert
    attempt = 0
    last_error = None
    while attempt < MAX_ATTEMPTS:
        attempt += 1
        try:
            fasta_bytes, source = fetch_fasta(entry, rate_limiter)
            pansn_fasta = rename_fasta_to_pansn(fasta_bytes, acc)
            conv = convert_to_bgzip(pansn_fasta, pansn_path, bgzip_bin, bgzip_threads)
            result.update({
                "status": "COMPLETE",
                "download_source": source,
                "download_attempts": attempt,
                "source_fasta_bytes": len(fasta_bytes),
            })
            result.update(conv)
            result["completed_at_utc"] = utcnow()
            return result
        except Exception as e:  # noqa: BLE001
            last_error = f"{type(e).__name__}: {e}"
            # clean partial output
            for suffix in ("", ".fai", ".gzi"):
                Path(str(pansn_path) + suffix).unlink(missing_ok=True)
            if attempt < MAX_ATTEMPTS:
                time.sleep(RETRY_DELAY * attempt)

    result["status"] = "FAILED"
    result["error"] = last_error
    result["completed_at_utc"] = utcnow()
    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def load_progress(path: Path) -> dict[str, Any]:
    if path.exists():
        try:
            return json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            pass
    return {"completed": {}, "failed": {}, "started_at": utcnow()}


def save_progress(path: Path, state: dict[str, Any]) -> None:
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True))
    os.rename(tmp, path)


def main() -> int:
    ap = argparse.ArgumentParser(description="NTM v2 PanSN downloader")
    ap.add_argument("--data-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs")
    ap.add_argument("--output-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/genomes")
    ap.add_argument("--v1-dir", default="/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--bgzip", default="bgzip")
    ap.add_argument("--bgzip-threads", type=int, default=4)
    ap.add_argument("--link-mode", choices=["symlink", "hardlink"], default="symlink")
    ap.add_argument("--rate-limit", type=float, default=0.0,
                    help="API calls/sec for non-FTP sources (0 = unlimited)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--retry-failed", action="store_true",
                    help="re-process accessions currently recorded as failed")
    ap.add_argument("--limit", type=int, default=0, help="process at most N accessions (testing)")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    v1_dir = Path(args.v1_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir = output_dir / ".scratch"
    scratch_dir.mkdir(parents=True, exist_ok=True)

    cohort, _ = build_cohort(data_dir)
    print(f"Cohort: {len(cohort)} unique assemblies "
          f"({sum(1 for c in cohort if c['source']=='GCA')} GCA, "
          f"{sum(1 for c in cohort if c['source']=='GCF')} GCF-only)", flush=True)

    # v1 canonical objects available for linking
    v1_co = v1_dir / "genomes" / "canonical_objects"
    v1_available: dict[str, Path] = {}
    if v1_co.is_dir():
        for d in os.listdir(v1_co):
            p = v1_co / d / f"{d}.pansn.fa.gz"
            if p.exists():
                v1_available[d] = p
    print(f"v1 canonical objects with pansn.fa.gz: {len(v1_available)}", flush=True)

    for entry in cohort:
        p = v1_available.get(entry["accession"])
        entry["v1_path"] = p if p is not None else None

    progress_path = output_dir / "progress.json"
    state = load_progress(progress_path)
    if not state.get("total"):
        state["total"] = len(cohort)
        state["started_at"] = utcnow()

    to_process = []
    for entry in cohort:
        acc = entry["accession"]
        skip = False
        if acc in state["completed"] and not args.force:
            skip = True
        if acc in state["failed"] and not (args.force or args.retry_failed):
            skip = True
        if skip:
            continue
        to_process.append(entry)
    if args.limit:
        to_process = to_process[: args.limit]
    print(f"Already resolved: {len(cohort) - len(to_process)}; to process: {len(to_process)}", flush=True)

    rate_limiter = RateLimiter(calls_per_second=args.rate_limit) if args.rate_limit > 0 else None

    done = len([a for a in state["completed"]])
    failed = len(state["failed"])
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_accession, e, output_dir, scratch_dir,
                          args.link_mode, args.bgzip, args.bgzip_threads,
                          rate_limiter, args.force): e["accession"] for e in to_process}
        batch = 0
        for fut in concurrent.futures.as_completed(futs):
            acc = futs[fut]
            try:
                res = fut.result()
            except Exception as e:  # noqa: BLE001
                res = {"accession": acc, "status": "FAILED", "error": f"{type(e).__name__}: {e}"}
            if res["status"] in ("COMPLETE", "LINKED", "ALREADY_COMPLETE"):
                state["completed"][acc] = res
            else:
                state["failed"][acc] = res
            batch += 1
            if batch % 50 == 0:
                state["updated_at"] = utcnow()
                save_progress(progress_path, state)
                n_ok = len(state["completed"]); n_fail = len(state["failed"])
                print(f"[{n_ok + n_fail}/{state['total']}] ok={n_ok} failed={n_fail} "
                      f"last={acc}:{res['status']}", flush=True)

    state["updated_at"] = utcnow()
    state["finished_at"] = utcnow()
    save_progress(progress_path, state)

    n_ok = len(state["completed"]); n_fail = len(state["failed"])
    print(f"\nSUMMARY: ok={n_ok} failed={n_fail} total={state['total']}", flush=True)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
