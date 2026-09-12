#!/usr/bin/env python3
"""
NTM v3 — acquire the cohort: link local v1/v2 objects, download the delta.

Delta (from build_v3_acquisition_manifest.py):
  * 16,696 link_v1v2 rows  -> symlink existing v1/v2 canonical objects into
                             ntm/v3/genomes/canonical_objects/{acc}/ (v2-style
                             file-level symlinks; frozen v1/v2 bytes untouched)
  * 70   download_bvbrc    -> BV-BRC API genome_sequence (contig accessions ==
                             coordinates scaffolds minus "accn|")
  * 25   download_ncbi     -> NCBI FTP (ASM dirname) / Datasets v2 ZIP / ENA
                             fallback — v2 download chain, exact export version
  * 18,055 blocked_run      -> terminal state, no public coordinate-compatible
                             source (collaborator delivery pending)

Layout (mirrors v1/v2):
  {out}/canonical_objects/{acc}/{acc}.pansn.fa.gz{,.fai,.gzi}
  {out}/acquisition_log.jsonl   per-genome ok/fail record
  {out}/progress.json           counters (resumable; existing objects skipped)

PanSN headers: >{acc}#1#{contig}, bgzip + samtools faidx (workflow/ntm
convention).

Usage:
  python3 ntm/v3/scripts/acquire_v3_genomes.py \
      [--work-dir /mnt/nvme3n1/erikg/phind-genome-work] [--force] [--workers 8]
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
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
MANIFEST_TSV = REPO / "ntm/v3/inputs/v3_acquisition_manifest.tsv"
GC_RE = re.compile(r"^GC([AF])_([0-9]{9})\.([0-9]+)$")
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
ENA_FASTA_BASE = "https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}?download=true&gzip=true"
BVBRC_API = "https://www.bv-brc.org/api"
USER_AGENT = "phind-ntm-v3-acquire/1.0"


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


class RateLimiter:
    def __init__(self, calls_per_second: float):
        self.min_interval = 1.0 / calls_per_second
        self._lock = threading.Lock()
        self._last = 0.0

    def acquire(self):
        with self._lock:
            now = time.time()
            elapsed = now - self._last
            if elapsed < self.min_interval:
                time.sleep(self.min_interval - elapsed)
            self._last = time.time()


NCBI_LIMITER = RateLimiter(2.0)
BVBRC_LIMITER = RateLimiter(3.0)
LOG_LOCK = threading.Lock()


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def http_get(url: str, limiter: RateLimiter | None, retries: int = 3,
             timeout: int = 300) -> bytes:
    last = None
    for attempt in range(retries + 1):
        if limiter:
            limiter.acquire()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            body = e.read()[:200]
            last = f"HTTP {e.code} {body!r}"
            if e.code in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(5.0 * (attempt + 1))
                continue
            raise RuntimeError(f"{url}: {last}") from e
        except Exception as e:  # noqa: BLE001
            last = str(e)
            if attempt < retries:
                time.sleep(5.0 * (attempt + 1))
                continue
            raise RuntimeError(f"{url}: {last}") from e
    raise RuntimeError(f"{url}: {last}")


def http_json(url: str, limiter: RateLimiter | None, retries: int = 4,
              timeout: int = 120):
    for attempt in range(retries):
        try:
            if limiter:
                limiter.acquire()
            req = urllib.request.Request(
                url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read())
        except Exception as e:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            time.sleep(5.0 * (attempt + 1))


def gunzip(data: bytes) -> bytes:
    if data[:2] != b"\x1f\x8b":
        raise RuntimeError("not a gzip stream")
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as fh:
        return fh.read()


# ------------------------------------------------------------------ FASTA
def contig_name_from_header(header: str) -> str:
    token = header.strip().split()[0] if header.strip().split() else "unknown"
    if "|" in token:
        token = token.split("|")[-1]
    return token.replace(" ", "_")


def rename_fasta_to_pansn(fasta_bytes: bytes, accession: str) -> bytes:
    out = []
    for line in fasta_bytes.split(b"\n"):
        if line.startswith(b">"):
            name = contig_name_from_header(line[1:].decode("utf-8", errors="replace"))
            out.append(f">{accession}#1#{name}".encode())
        else:
            out.append(line)
    return b"\n".join(out)


def wrap_seq(seq: str, width: int = 60) -> list[str]:
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def write_bgzip_faidx(fasta_bytes: bytes, canonical_dir: Path, acc: str,
                      log: dict) -> dict:
    """Write PanSN fasta, bgzip it, faidx + gzi. Returns object stats."""
    gz_path = canonical_dir / f"{acc}.pansn.fa.gz"
    tmp_fa = canonical_dir / f"{acc}.fa.tmp"
    tmp_gz = canonical_dir / f"{acc}.pansn.fa.gz.tmp"
    tmp_fa.write_bytes(fasta_bytes)
    try:
        res = subprocess.run(["bgzip", "-c", "-@2", str(tmp_fa)],
                             capture_output=True, check=True)
        tmp_gz.write_bytes(res.stdout)
        os.replace(tmp_gz, gz_path)
    finally:
        tmp_fa.unlink(missing_ok=True)
        tmp_gz.unlink(missing_ok=True)
    subprocess.run(["samtools", "faidx", str(gz_path)],
                   capture_output=True, check=True, timeout=300)
    subprocess.run(["bgzip", "-r", str(gz_path)],
                   capture_output=True, check=True, timeout=300)
    fai = gz_path.with_suffix(".gz.fai") if False else canonical_dir / f"{acc}.pansn.fa.gz.fai"
    contigs = 0
    total_bp = 0
    with open(fai) as fh:
        for line in fh:
            contigs += 1
            total_bp += int(line.split("\t")[1])
    log.update({"canonical_bgzf_bytes": gz_path.stat().st_size,
                "contigs": contigs, "total_bp": total_bp})
    return {"contigs": contigs, "total_bp": total_bp}


# ------------------------------------------------------------- BV-BRC fetch
def fetch_bvbrc(genome_id: str) -> tuple[bytes, dict]:
    """Fetch all sequences of a BV-BRC genome; return (pansn_fasta_bytes, info)."""
    url = (f"{BVBRC_API}/genome_sequence/?eq(genome_id,"
           f"{urllib.parse.quote(genome_id)})"
           "&select(sequence_id,accession,length,sequence,sequence_md5,sequence_type)")
    seqs = http_json(url, BVBRC_LIMITER)
    if not seqs:
        raise RuntimeError(f"BV-BRC returned 0 sequences for {genome_id}")
    out = []
    md5_ok = 0
    md5_bad = []
    for s in seqs:
        acc = s.get("accession") or s.get("sequence_id")
        seq = (s.get("sequence") or "").strip()
        if not seq:
            raise RuntimeError(f"BV-BRC sequence {acc} of {genome_id} empty")
        if s.get("sequence_md5"):
            import hashlib as _h
            if _h.md5(seq.encode()).hexdigest() == s["sequence_md5"]:
                md5_ok += 1
            else:
                md5_bad.append(acc)
        if len(seq) != s.get("length"):
            raise RuntimeError(
                f"length mismatch for {acc}: got {len(seq)} expected {s.get('length')}")
        out.append(f">{genome_id}#1#{acc}")
        out.extend(wrap_seq(seq.upper()))
    info = {"n_seqs": len(seqs), "md5_ok": md5_ok, "md5_bad": md5_bad,
            "total_bp": sum(s.get("length", 0) for s in seqs),
            "seq_accessions": ",".join(
                (s.get("accession") or s.get("sequence_id")) for s in seqs)}
    return ("\n".join(out) + "\n").encode(), info


# ------------------------------------------------------------- NCBI fetch
def ftp_url_for(acc: str, genome_key: str) -> str | None:
    m = GC_RE.match(acc)
    if not m:
        return None
    digits = m.group(2)
    if not genome_key.endswith("_genomic"):
        return None
    dirname = genome_key[: -len("_genomic")]
    return (f"{FTP_BASE}/{m.group(1)}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
            f"{dirname}/{dirname}_genomic.fna.gz")


def fetch_ncbi(acc: str, genome_key: str) -> tuple[bytes, str]:
    """Source chain: NCBI FTP -> Datasets ZIP -> ENA. Returns (fasta, source)."""
    errors = []
    url = ftp_url_for(acc, genome_key)
    if url:
        try:
            return gunzip(http_get(url, NCBI_LIMITER)), f"NCBI FTP {url}"
        except Exception as e:  # noqa: BLE001
            errors.append(f"ftp: {e}")
    try:
        q = urllib.parse.urlencode(
            {"include_annotation_type": ["GENOME_FASTA"], "filename": f"{acc}.zip"},
            doseq=True)
        data = http_get(DATASETS_BASE.format(acc=acc) + "?" + q, NCBI_LIMITER)
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            member = next((n for n in zf.namelist() if n.endswith("_genomic.fna")), None)
            if member is None:
                raise RuntimeError(f"no _genomic.fna in datasets ZIP")
            return zf.read(member), f"NCBI datasets ZIP {acc}"
    except Exception as e:  # noqa: BLE001
        errors.append(f"datasets: {e}")
    try:
        return gunzip(http_get(ENA_FASTA_BASE.format(acc=acc), NCBI_LIMITER)), \
            f"ENA fasta API {acc}"
    except Exception as e:  # noqa: BLE001
        errors.append(f"ena: {e}")
    raise RuntimeError("; ".join(errors))


# ------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", default="/mnt/nvme3n1/erikg/phind-genome-work")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    out_dir = Path(args.work_dir) / "ntm/v3/genomes"
    co_dir = out_dir / "canonical_objects"
    co_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "acquisition_log.jsonl"

    rows = list(csv.DictReader(open(MANIFEST_TSV, newline=""), delimiter="\t"))
    print(f"manifest rows: {len(rows)}")

    # prune orphan objects from previous runs (only pure-symlink dirs)
    keep = {r["canonical_acc"] for r in rows if r["canonical_acc"]}
    n_pruned = 0
    for d in sorted(os.listdir(co_dir)):
        if d in keep:
            continue
        p = co_dir / d
        if p.is_dir() and all(os.path.islink(str(f)) or f.is_symlink()
                               for f in p.iterdir()):
            for f in p.iterdir():
                f.unlink()
            p.rmdir()
            n_pruned += 1
    if n_pruned:
        print(f"pruned {n_pruned} orphan symlink objects from earlier runs")

    log_fh = open(log_path, "a")

    def log(rec: dict):
        rec["ts"] = utcnow()
        with LOG_LOCK:
            log_fh.write(json.dumps(rec) + "\n")
            log_fh.flush()

    # ------------------------------------------------------- phase 1: link
    links: dict[str, dict] = {}  # canonical_acc -> row
    for r in rows:
        if r["plan"] == "link_v1v2":
            links.setdefault(r["canonical_acc"], r)
    print(f"link targets: {len(links)} objects")
    n_linked = n_linked_existing = 0
    link_errors = []
    for acc, r in sorted(links.items()):
        rec = {"genome_key": r["genome_key"], "canonical_acc": acc,
               "plan": "link_v1v2"}
        dst_dir = co_dir / acc
        dst_dir.mkdir(parents=True, exist_ok=True)
        src_gz = Path(r["link_gz"])
        if not src_gz.exists():
            rec.update({"status": "fail", "error": f"source missing {src_gz}"})
            log(rec)
            link_errors.append(acc)
            continue
        dst_gz = dst_dir / f"{acc}.pansn.fa.gz"
        if dst_gz.exists() and not args.force:
            n_linked_existing += 1
            rec.update({"status": "ok", "note": "link already present"})
            log(rec)
            continue
        for suffix in ("", ".fai", ".gzi"):
            src = Path(str(src_gz) + suffix)
            dst = dst_dir / f"{acc}.pansn.fa.gz{suffix}"
            if dst.is_symlink() or dst.exists():
                if args.force:
                    dst.unlink()
                else:
                    continue
            if src.exists():
                os.symlink(src, dst)
            elif suffix:  # index missing at source: regenerate below via faidx
                pass
        # ensure indexes exist (regenerated locally if the source lacks them)
        fai = dst_dir / f"{acc}.pansn.fa.gz.fai"
        gzi = dst_dir / f"{acc}.pansn.fa.gz.gzi"
        try:
            if not fai.exists():
                subprocess.run(["samtools", "faidx", "-o", str(fai), str(dst_gz)],
                               capture_output=True, check=True, timeout=300)
            if not gzi.exists():
                subprocess.run(["bgzip", "-r", str(dst_gz)],
                               capture_output=True, check=True, timeout=300)
            rec.update({"status": "ok", "link_target": str(src_gz)})
            n_linked += 1
        except subprocess.CalledProcessError as e:
            rec.update({"status": "fail", "error": f"faidx: {e.stderr[:200]}"})
            link_errors.append(acc)
        log(rec)
    print(f"linked {n_linked} objects ({n_linked_existing} already present); "
          f"errors {len(link_errors)}")

    # ------------------------------------------- phase 2/3: downloads
    dl_rows: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        if r["plan"] in ("download_bvbrc", "download_ncbi"):
            dl_rows[r["canonical_acc"]].append(r)
    print(f"download targets: {len(dl_rows)} objects")

    def do_download(pair):
        acc, rlist = pair
        r0 = rlist[0]
        rec = {"genome_key": r0["genome_key"], "canonical_acc": acc,
               "genome_keys": [x["genome_key"] for x in rlist],
               "plan": r0["plan"]}
        dst_dir = co_dir / acc
        dst_gz = dst_dir / f"{acc}.pansn.fa.gz"
        if dst_gz.exists() and (dst_dir / f"{acc}.pansn.fa.gz.fai").exists() \
                and not args.force:
            rec.update({"status": "ok", "note": "already downloaded"})
            log(rec)
            return rec
        try:
            dst_dir.mkdir(parents=True, exist_ok=True)
            if r0["plan"] == "download_bvbrc":
                fasta, info = fetch_bvbrc(acc)
                stats = write_bgzip_faidx(fasta, dst_dir, acc, rec)
                rec.update({"status": "ok", "source": "BV-BRC API",
                            "bvbrc_seqs": info["n_seqs"],
                            "bvbrc_md5_ok": info["md5_ok"],
                            "bvbrc_md5_bad": ",".join(info["md5_bad"]),
                            **stats})
            else:
                fasta, source = fetch_ncbi(acc, r0["genome_key"])
                pansn = rename_fasta_to_pansn(fasta, acc)
                stats = write_bgzip_faidx(pansn, dst_dir, acc, rec)
                rec.update({"status": "ok", "source": source, **stats})
        except Exception as e:  # noqa: BLE001
            rec.update({"status": "fail", "error": str(e)[:500]})
        log(rec)
        return rec

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for rec in ex.map(do_download, sorted(dl_rows.items())):
            results.append(rec)
    dl_ok = sum(1 for r in results if r.get("status") == "ok")
    dl_fail = [r["canonical_acc"] for r in results if r.get("status") == "fail"]
    print(f"downloads ok {dl_ok}, failed {len(dl_fail)}: {dl_fail[:10]}")

    # ------------------------------------------------ phase 4: blocked runs
    blocked = [r for r in rows if r["plan"] == "blocked_run"]
    print(f"blocked run assemblies (terminal state): {len(blocked)}")

    progress = {
        "updated_utc": utcnow(),
        "manifest_rows": len(rows),
        "link_targets": len(links),
        "linked_ok": n_linked,
        "linked_existing": n_linked_existing,
        "link_errors": link_errors,
        "download_targets": len(dl_rows),
        "downloaded_ok": dl_ok,
        "download_fail": dl_fail,
        "blocked_runs": len(blocked),
    }
    (out_dir / "progress.json").write_text(json.dumps(progress, indent=1))
    print(json.dumps({k: v for k, v in progress.items() if k != "link_errors"},
                     indent=1))
    log_fh.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
