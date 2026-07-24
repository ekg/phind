#!/usr/bin/env python3
"""Bounded, fail-closed Stage-B acquisition and PanSN/BGZF canonicalization.

The only payload endpoint this program can construct is the NCBI Datasets v2
single-accession download endpoint for one of the ten identities in the
checksum-pinned Stage-B manifest.  Publication is a same-filesystem directory
rename after per-object and release inventories are complete.
"""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import os
import platform
import re
import resource
import shutil
import signal
import socket
import stat
import struct
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
import zipfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, BinaryIO, Iterable

SCHEMA = "canonical-cohort-010-release-v1"
POLICY = "pansn-bgzip-policy-v1"
TASK_ID = "run-10-assembly-acquisition"
COLLECTION_RELEASE_ID = "collection-v1-f7494b4b89d1382b"
COLLECTION_RELEASE_SHA256 = "59c6907e2c053e9d8ac3df8d5eb820bab0097030a9259ca2c9354c47cb6642bf"
STAGE_B_SHA256 = "0d179cbafce2ba1fa14d1929a4acd6621810a335f25bcd7ec67dd2083eb101f6"
ACCESSION_INPUT_SHA256 = "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5"
PROPHAGE_INPUT_SHA256 = "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996"
EXPECTED_ACCESSIONS = [
    "GCF_000005845.2", "GCF_000812325.1", "GCF_002302315.1",
    "GCF_004664255.1", "GCF_015644385.1", "GCF_020829045.1",
    "GCF_921380995.1", "GCF_000167895.3", "GCF_001881595.4",
    "GCF_000498835.2",
]
ACCESSION_RE = re.compile(r"^GC[AF]_[0-9]{9}\.[1-9][0-9]*$")
SEQ_ACCESSION_RE = re.compile(r"^([A-Za-z]{1,6}_?[0-9]+)\.([1-9][0-9]*)$")
DNA = frozenset(b"ACGTRYSWKMBDHVN")
SAFE_CONTIG = frozenset(b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-")
GFF_SAFE = frozenset(b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:^*$@!+_?-|")
DOWNLOAD_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
INCLUDE_TYPES = ("GENOME_FASTA", "GENOME_GFF", "SEQUENCE_REPORT")
USER_AGENT = "phind-stage-b-acquisition/1.0 (bounded-10-assembly-pilot)"

ASSEMBLY_FIELDS = [
    "stage_b_order", "assembly_id", "accession", "predecessor_resolution_status",
    "terminal_state", "source_object_relpath", "source_package_bytes",
    "source_package_sha256", "source_fasta_member", "source_fasta_bytes",
    "source_decompressed_sha256", "source_gff_member", "source_gff_sha256",
    "annotation_status", "canonical_object_relpath", "canonical_bgzf_relpath",
    "canonical_bgzf_bytes", "canonical_fasta_content_sha256", "canonical_bgzf_sha256",
    "fai_sha256", "gzi_sha256", "crosswalk_sha256", "annotation_aliases_sha256",
    "contig_count", "total_bases", "remote_identity_strength", "remote_etag",
    "remote_last_modified", "download_attempts", "range_resumes", "validated_at_utc",
    "row_sha256",
]
CONTIG_FIELDS = [
    "stage_b_order", "assembly_id", "accession", "contig_order", "source_fasta_member",
    "source_fasta_header_b64", "source_fasta_id_token_b64", "source_contig_id_display",
    "source_contig_accession", "source_contig_version", "source_contig_accession_version",
    "genbank_contig_accession_version", "refseq_contig_accession_version",
    "assembly_report_sequence_name", "contig_id_encoding", "pansn_sample",
    "pansn_haplotype", "pansn_contig", "pansn_sequence_name", "fasta_seqid",
    "replicon_role", "assigned_molecule", "plasmid_name_b64", "topology",
    "contig_length", "contig_sequence_sha256", "canonical_bgzf_relpath", "row_sha256",
]
CHECKSUM_FIELDS = ["stage_b_order", "accession", "artifact_role", "relative_path", "bytes", "sha256", "row_sha256"]
STATE_FIELDS = ["stage_b_order", "accession", "source_state", "canonical_state", "terminal_state", "reason", "row_sha256"]


class GateError(RuntimeError):
    """A release-blocking gate failure."""


class AnnotationValidationError(GateError):
    """A source-annotation failure that forbids an alias view, not genome use."""


@dataclass(frozen=True)
class Allocations:
    assigned_ram_bytes: int
    durable_allocation_bytes: int
    scratch_allocation_bytes: int
    inode_allocation: int
    predicted_durable_peak_bytes: int
    predicted_scratch_peak_bytes: int
    predicted_files: int
    unfinished_write_bytes: int

    def validate(self) -> None:
        if any(v <= 0 for v in asdict(self).values()):
            raise GateError("resource allocation/reservation is blank or non-positive")
        if self.predicted_durable_peak_bytes * 100 > self.durable_allocation_bytes * 70:
            raise GateError("predicted durable upper-95% peak exceeds 70% allocation")
        if self.predicted_scratch_peak_bytes * 100 > self.scratch_allocation_bytes * 70:
            raise GateError("predicted scratch upper-95% peak exceeds 70% allocation")
        if self.predicted_files * 2 > self.inode_allocation:
            raise GateError("projected files exceed 50% inode allocation")


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_json(obj: Any) -> bytes:
    return (json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_stream(handle: BinaryIO) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
        size += len(block)
    return digest.hexdigest(), size


def append_jsonl(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("ab") as handle:
        handle.write(canonical_json(row))
        handle.flush()
        os.fsync(handle.fileno())


def fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def fsync_dir(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def safe_remove(path: Path) -> None:
    if path.is_symlink():
        raise GateError(f"refusing symlink removal: {path}")
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def stable_row_hash(row: dict[str, Any], fields: list[str]) -> str:
    return sha_bytes(("\t".join(str(row.get(field, ".")) for field in fields if field != "row_sha256") + "\n").encode())


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for source in rows:
            row = {field: ("." if source.get(field) in (None, "") else str(source[field])) for field in fields}
            row["row_sha256"] = stable_row_hash(row, fields)
            writer.writerow(row)
    fsync_file(path)


def read_tsv(path: Path, fields: list[str] | None = None, verify_hashes: bool = False) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed_fields = reader.fieldnames
        if fields is not None and observed_fields != fields:
            raise GateError(f"TSV schema mismatch: {path}")
        rows = list(reader)
    if verify_hashes:
        actual_fields = fields or observed_fields
        if actual_fields is None:
            raise GateError(f"cannot verify headerless TSV: {path}")
        for number, row in enumerate(rows, 2):
            if row.get("row_sha256") != stable_row_hash(row, actual_fields):
                raise GateError(f"row checksum mismatch: {path}:{number}")
    return rows


def deterministic_gzip(src: Path, dst: Path) -> None:
    with src.open("rb") as source, dst.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as target:
            shutil.copyfileobj(source, target)
    fsync_file(dst)


def verify_sha_inventory(root: Path, require_complete: bool = True) -> None:
    sums = root / "SHA256SUMS"
    if not sums.is_file():
        raise GateError(f"missing SHA256SUMS: {root}")
    if require_complete:
        complete = root / "COMPLETE"
        if not complete.is_file():
            raise GateError(f"missing COMPLETE: {root}")
        tokens = complete.read_text().split()
        if len(tokens) < 1 or tokens[0] != sha_file(sums):
            raise GateError(f"COMPLETE inventory digest mismatch: {root}")
    seen: set[str] = set()
    for number, line in enumerate(sums.read_text().splitlines(), 1):
        try:
            digest, rel = line.split("  ", 1)
        except ValueError as exc:
            raise GateError(f"malformed inventory line {root}:{number}") from exc
        if rel in seen or rel.startswith("/") or ".." in Path(rel).parts:
            raise GateError(f"unsafe/duplicate inventory path: {rel}")
        seen.add(rel)
        path = root / rel
        if not path.is_file() or path.is_symlink() or sha_file(path) != digest:
            raise GateError(f"inventory checksum mismatch: {root}/{rel}")


def seal_directory(stage: Path, final: Path, inventory_exclude: set[str] | None = None) -> None:
    if final.exists():
        raise GateError(f"refusing overwrite of final object: {final}")
    excluded = {"SHA256SUMS", "COMPLETE"} | (inventory_exclude or set())
    entries = []
    for path in sorted(p for p in stage.rglob("*") if p.is_file()):
        relative = str(path.relative_to(stage))
        if relative in excluded:
            continue
        if path.is_symlink():
            raise GateError(f"symlink in staged object: {path}")
        entries.append(f"{sha_file(path)}  {relative}\n")
    (stage / "SHA256SUMS").write_text("".join(entries))
    fsync_file(stage / "SHA256SUMS")
    for path in (p for p in stage.rglob("*") if p.is_file()):
        fsync_file(path)
    (stage / "COMPLETE").write_text(f"{sha_file(stage / 'SHA256SUMS')}  SHA256SUMS\n")
    fsync_file(stage / "COMPLETE")
    fsync_dir(stage)
    os.rename(stage, final)
    fsync_dir(final.parent)
    verify_sha_inventory(final)


def _existing_parent(path: Path) -> Path:
    candidate = path
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    return candidate


def _mount_record(path: Path) -> dict[str, str]:
    result = subprocess.run(
        ["findmnt", "-J", "-T", str(path), "-o", "TARGET,SOURCE,FSTYPE,OPTIONS,SIZE,AVAIL"],
        check=True, text=True, capture_output=True,
    )
    filesystems = json.loads(result.stdout).get("filesystems", [])
    if len(filesystems) != 1:
        raise GateError(f"findmnt returned no unique filesystem for {path}")
    return {str(k): str(v) for k, v in filesystems[0].items()}


def _write_probe(directory: Path, label: str) -> None:
    name = f".acq-write-probe.{label}.{os.getpid()}.{uuid.uuid4().hex}"
    path = directory / name
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(path, flags, 0o600)
    try:
        payload = f"{TASK_ID}\t{label}\t{utcnow()}\n".encode()
        if os.write(fd, payload) != len(payload):
            raise GateError("short write during filesystem probe")
        os.fsync(fd)
    finally:
        os.close(fd)
    if path.stat().st_size == 0:
        raise GateError("empty filesystem write probe")
    path.unlink()
    if path.exists():
        raise GateError("filesystem write probe cleanup failed")


def swap_free_bytes() -> int:
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("SwapFree:"):
            return int(line.split()[1]) * 1024
    raise GateError("cannot read SwapFree")


def live_preflight(durable: Path, scratch: Path, allocations: Allocations, stage_name: str) -> dict[str, Any]:
    allocations.validate()
    durable_probe_parent = _existing_parent(durable)
    scratch_probe_parent = _existing_parent(scratch)
    durable_mount = _mount_record(durable_probe_parent)
    scratch_mount = _mount_record(scratch_probe_parent)
    if durable_mount.get("target") != "/" or durable_mount.get("fstype") != "ext4":
        raise GateError(f"durable mount identity NO_GO: {durable_mount}")
    if scratch_mount.get("target") != "/mnt/nvme3n1" or scratch_mount.get("source") != "/dev/nvme3n1" or scratch_mount.get("fstype") != "xfs":
        raise GateError(f"scratch mount identity NO_GO: {scratch_mount}")
    if not str(durable.resolve()).startswith("/home/erikg/phind-data/ecoli26k/v1/releases/run-10-assembly-acquisition/"):
        raise GateError(f"durable path outside task namespace: {durable}")
    if not str(scratch.resolve()).startswith("/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-10-assembly-acquisition/"):
        raise GateError(f"scratch path outside task namespace: {scratch}")
    ds = os.statvfs(durable_probe_parent)
    ss = os.statvfs(scratch_probe_parent)
    durable_free = ds.f_bavail * ds.f_frsize
    durable_inodes = ds.f_favail
    scratch_free = ss.f_bavail * ss.f_frsize
    scratch_inodes = ss.f_favail
    checks = {
        "durable_start_bytes_ge_2_4tb": durable_free >= 2_400_000_000_000,
        "durable_stop_bytes_ge_2tb_after_next_write": durable_free - allocations.unfinished_write_bytes >= 2_000_000_000_000,
        "durable_inodes_ge_1m": durable_inodes >= 1_000_000,
        "scratch_start_bytes_ge_4tb": scratch_free >= 4_000_000_000_000,
        "scratch_stop_bytes_ge_2tb_after_next_write": scratch_free - allocations.unfinished_write_bytes >= 2_000_000_000_000,
        "scratch_inodes_ge_5m": scratch_inodes >= 5_000_000,
        "durable_disk_allocation_le_70pct": allocations.predicted_durable_peak_bytes * 100 <= allocations.durable_allocation_bytes * 70,
        "scratch_disk_allocation_le_70pct": allocations.predicted_scratch_peak_bytes * 100 <= allocations.scratch_allocation_bytes * 70,
        "projected_files_le_50pct_inode_allocation": allocations.predicted_files * 2 <= allocations.inode_allocation,
        "durable_two_x_unfinished_writes": durable_free - allocations.predicted_durable_peak_bytes >= 2 * allocations.unfinished_write_bytes,
        "scratch_two_x_unfinished_writes": scratch_free - allocations.predicted_scratch_peak_bytes >= 2 * allocations.unfinished_write_bytes,
    }
    if not all(checks.values()):
        raise GateError("resource gate NO_GO: " + json.dumps(checks, sort_keys=True))
    uid = os.getuid()
    for label, parent in (("durable", durable_probe_parent), ("scratch", scratch_probe_parent)):
        if parent.stat().st_uid != uid:
            raise GateError(f"{label} existing parent not owned by current Erik uid: {parent}")
        if not os.access(parent, os.W_OK | os.X_OK):
            raise GateError(f"{label} existing parent not writable/searchable: {parent}")
        _write_probe(parent, label)
    return {
        "schema": "resource-preflight-v1", "verdict": "PASS", "stage": stage_name,
        "captured_at_utc": utcnow(), "durable_path": str(durable), "scratch_path": str(scratch),
        "durable_probe_parent": str(durable_probe_parent), "scratch_probe_parent": str(scratch_probe_parent),
        "durable_findmnt": durable_mount, "scratch_findmnt": scratch_mount,
        "durable_owner": {"uid": durable_probe_parent.stat().st_uid, "gid": durable_probe_parent.stat().st_gid, "mode": stat.S_IMODE(durable_probe_parent.stat().st_mode)},
        "scratch_owner": {"uid": scratch_probe_parent.stat().st_uid, "gid": scratch_probe_parent.stat().st_gid, "mode": stat.S_IMODE(scratch_probe_parent.stat().st_mode)},
        "write_probes": "PASS", "durable_free_bytes": durable_free, "durable_free_inodes": durable_inodes,
        "scratch_free_bytes": scratch_free, "scratch_free_inodes": scratch_inodes,
        "swap_free_bytes": swap_free_bytes(), "allocations": asdict(allocations), "checks": checks,
    }


def mkdir_owned(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True, mode=0o750)
    os.chmod(path, 0o750)
    if path.stat().st_uid != os.getuid() or path.is_symlink():
        raise GateError(f"directory is not an Erik-owned real directory: {path}")


def verify_predecessor(repo: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    tracked = repo / "manifests/collection-v1"
    sums = {}
    for line in (tracked / "SHA256SUMS").read_text().splitlines():
        digest, rel = line.split("  ", 1)
        sums[rel] = digest
        if sha_file(tracked / rel) != digest:
            raise GateError(f"predecessor tracked checksum mismatch: {rel}")
    if sums.get("release.json") != COLLECTION_RELEASE_SHA256:
        raise GateError("predecessor release.json recorded SHA-256 mismatch")
    if sums.get("stage_b_10.tsv") != STAGE_B_SHA256 or sha_file(tracked / "stage_b_10.tsv") != STAGE_B_SHA256:
        raise GateError("Stage-B manifest SHA-256 mismatch")
    release = json.loads((tracked / "release.json").read_text())
    if release.get("release_id") != COLLECTION_RELEASE_ID or release.get("verdict") != "PASS" or release.get("immutable") is not True:
        raise GateError("predecessor release ID/verdict/immutability gate failed")
    for gate, verdict in release.get("applicable_gates", {}).items():
        if verdict not in ("PASS", "NOT_APPLICABLE_METADATA_ONLY"):
            raise GateError(f"predecessor gate is not unqualified PASS/NA: {gate}={verdict}")
    manifest = release.get("manifests", {}).get("stage_b_10.tsv", {})
    if manifest.get("sha256") != STAGE_B_SHA256 or manifest.get("rows") != 10 or manifest.get("bytes") != 2246:
        raise GateError("predecessor Stage-B release contract mismatch")
    external = Path(release["external_release_path"])
    verify_sha_inventory(external)
    external_stage = external / "manifests/stage_b_10.tsv"
    if sha_file(external_stage) != STAGE_B_SHA256 or external_stage.read_bytes() != (tracked / "stage_b_10.tsv").read_bytes():
        raise GateError("tracked/external Stage-B bytes differ")
    rows = read_tsv(tracked / "stage_b_10.tsv", verify_hashes=True)
    accessions = [row["resolved_assembly_accession_version"] for row in rows]
    if len(rows) != 10 or accessions != EXPECTED_ACCESSIONS or len(set(accessions)) != 10:
        raise GateError("Stage-B exact order/cardinality/identity gate failed")
    for order, row in enumerate(rows, 1):
        if int(row["stage_b_order"]) != order or row["requested_assembly_accession_version"] != row["resolved_assembly_accession_version"]:
            raise GateError("Stage-B order or exact-version equality failed")
        expected_status = "EXACT_VERSION_VALID_METADATA_UNAVAILABLE" if row["resolved_assembly_accession_version"] == "GCF_000167895.3" else "EXACT_VERSION_RESOLVED"
        if row["resolution_status"] != expected_status:
            raise GateError("documented metadata-omission status mismatch")
    return rows, release


def verify_root_inputs(repo: Path) -> dict[str, str]:
    observed = {
        "26k_ecoli_accession.txt": sha_file(repo / "26k_ecoli_accession.txt"),
        "26k_prophage1.csv": sha_file(repo / "26k_prophage1.csv"),
    }
    if observed["26k_ecoli_accession.txt"] != ACCESSION_INPUT_SHA256 or observed["26k_prophage1.csv"] != PROPHAGE_INPUT_SHA256:
        raise GateError("immutable root input checksum mismatch")
    return observed


def release_seed() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA, "collection_release_id": COLLECTION_RELEASE_ID,
        "collection_release_json_sha256": COLLECTION_RELEASE_SHA256,
        "stage_b_10_sha256": STAGE_B_SHA256, "stage_b_rows": 10,
        "pansn_policy_version": POLICY, "exact_accessions": EXPECTED_ACCESSIONS,
        "source": "NCBI Datasets v2 exact single-accession packages",
    }


def release_id() -> str:
    return "canonical-cohort-010-v1-" + sha_bytes(canonical_json(release_seed()))[:16]


def download_url(accession: str) -> str:
    if accession not in EXPECTED_ACCESSIONS or not ACCESSION_RE.fullmatch(accession):
        raise GateError(f"payload request outside immutable Stage-B cohort: {accession}")
    query = urllib.parse.urlencode({"include_annotation_type": list(INCLUDE_TYPES), "filename": f"{accession}.zip"}, doseq=True)
    return DOWNLOAD_BASE.format(accession=accession) + "?" + query


def _selected_headers(response: Any) -> dict[str, str]:
    names = ("ETag", "Last-Modified", "Content-Length", "Accept-Ranges", "Content-Range", "Content-Type", "Content-Disposition", "X-RateLimit-Limit")
    return {name.lower().replace("-", "_"): response.headers.get(name, ".") for name in names}


def _identity_from_response(response: Any, requested_url: str, resumed: bool = False) -> dict[str, Any]:
    headers = _selected_headers(response)
    total = None
    content_range = headers.get("content_range", ".")
    if content_range != ".":
        match = re.fullmatch(r"bytes [0-9]+-[0-9]+/([0-9]+|\*)", content_range)
        if match and match.group(1) != "*":
            total = int(match.group(1))
    if total is None and not resumed and headers.get("content_length", ".").isdigit():
        total = int(headers["content_length"])
    etag = headers.get("etag", ".")
    last_modified = headers.get("last_modified", ".")
    strong = etag != "." and not etag.startswith("W/")
    strength = "STRONG_ETAG" if strong else ("LAST_MODIFIED_AND_LENGTH" if last_modified != "." and total is not None else "INSUFFICIENT_FOR_RANGE_RESUME")
    return {
        "requested_url": requested_url, "final_url": response.geturl(), "status": response.status,
        "headers": headers, "total_bytes": total, "identity_strength": strength,
        "captured_at_utc": utcnow(),
    }


def _head_identity(url: str) -> dict[str, Any] | None:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT, "Accept": "application/zip"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return _identity_from_response(response, url)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
        return None


def _same_remote_identity(stored: dict[str, Any], current: dict[str, Any] | None) -> bool:
    if current is None or stored.get("identity_strength") == "INSUFFICIENT_FOR_RANGE_RESUME":
        return False
    if stored.get("final_url") != current.get("final_url") or stored.get("total_bytes") != current.get("total_bytes"):
        return False
    sh = stored.get("headers", {})
    ch = current.get("headers", {})
    if stored.get("identity_strength") == "STRONG_ETAG":
        return sh.get("etag") == ch.get("etag")
    return sh.get("last_modified") == ch.get("last_modified") and stored.get("total_bytes") is not None


def acquire_package(
    accession: str, object_stage: Path, state: Path, failures: Path,
    retries: int, rate_delay: float, inject_kill: bool, inject_after_bytes: int,
) -> tuple[Path, dict[str, Any]]:
    url = download_url(accession)
    partial = object_stage / "package.zip.partial"
    identity_path = object_stage / "remote_identity.json"
    receipt_path = object_stage / "download_receipt.json"
    complete_download = object_stage / "package.zip"
    if complete_download.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt.get("local_sha256") != sha_file(complete_download) or receipt.get("accession") != accession or receipt.get("url") != url:
            raise GateError(f"completed download receipt mismatch: {accession}")
        return complete_download, receipt
    attempts = 0
    range_resumes = 0
    for attempt in range(1, retries + 1):
        attempts += 1
        start = 0
        stored: dict[str, Any] | None = None
        if partial.exists() and identity_path.exists():
            stored = json.loads(identity_path.read_text())
            current = _head_identity(url)
            if _same_remote_identity(stored, current) and stored.get("headers", {}).get("accept_ranges", "").lower() == "bytes":
                start = partial.stat().st_size
                if stored.get("total_bytes") is not None and start >= int(stored["total_bytes"]):
                    start = 0
                else:
                    range_resumes += 1
                    append_jsonl(state, {"event": "ACQUISITION_RANGE_RESUME", "accession": accession, "offset": start, "at": utcnow()})
            else:
                append_jsonl(state, {"event": "ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE", "accession": accession, "partial_bytes": partial.stat().st_size, "at": utcnow()})
                partial.unlink()
                identity_path.unlink(missing_ok=True)
                start = 0
        elif partial.exists() or identity_path.exists():
            partial.unlink(missing_ok=True)
            identity_path.unlink(missing_ok=True)
            append_jsonl(state, {"event": "ACQUISITION_ORPHAN_PARTIAL_DISCARDED", "accession": accession, "at": utcnow()})
        headers = {"User-Agent": USER_AGENT, "Accept": "application/zip"}
        if start:
            headers["Range"] = f"bytes={start}-"
            if stored and stored.get("headers", {}).get("etag", ".") != ".":
                headers["If-Range"] = stored["headers"]["etag"]
            elif stored:
                headers["If-Range"] = stored.get("headers", {}).get("last_modified", "")
        request = urllib.request.Request(url, method="GET", headers=headers)
        try:
            append_jsonl(state, {"event": "ACQUISITION_REQUEST_STARTED", "accession": accession, "attempt": attempt, "range_start": start, "url": url, "at": utcnow()})
            with urllib.request.urlopen(request, timeout=300) as response:
                if start and response.status != 206:
                    partial.unlink(missing_ok=True)
                    identity_path.unlink(missing_ok=True)
                    append_jsonl(state, {"event": "ACQUISITION_RANGE_REFUSED_RESTART", "accession": accession, "status": response.status, "at": utcnow()})
                    continue
                if not start and response.status != 200:
                    raise GateError(f"unexpected download HTTP status {response.status} for {accession}")
                identity = _identity_from_response(response, url, bool(start))
                if start and stored and identity.get("total_bytes") != stored.get("total_bytes"):
                    raise GateError(f"remote total size changed during range resume: {accession}")
                identity_path.write_bytes(canonical_json(identity))
                fsync_file(identity_path)
                mode = "ab" if start else "wb"
                with partial.open(mode) as output:
                    written = start
                    while True:
                        block = response.read(64 * 1024)
                        if not block:
                            break
                        output.write(block)
                        written += len(block)
                        if written % (1024 * 1024) < len(block):
                            output.flush()
                            os.fsync(output.fileno())
                        if inject_kill and written >= inject_after_bytes:
                            output.flush()
                            os.fsync(output.fileno())
                            append_jsonl(state, {"event": "INJECTED_ACQUISITION_SIGKILL", "accession": accession, "partial_bytes": written, "at": utcnow()})
                            os.kill(os.getpid(), signal.SIGKILL)
                    output.flush()
                    os.fsync(output.fileno())
                expected_total = (stored or identity).get("total_bytes")
                if expected_total is not None and partial.stat().st_size != int(expected_total):
                    raise GateError(f"download byte-count mismatch for {accession}: {partial.stat().st_size} != {expected_total}")
                os.rename(partial, complete_download)
                fsync_dir(object_stage)
                receipt = {
                    "schema": "ncbi-package-download-receipt-v1", "accession": accession,
                    "url": url, "remote_identity": stored or identity, "completion_response_identity": identity,
                    "local_bytes": complete_download.stat().st_size, "local_sha256": sha_file(complete_download),
                    "attempts_this_invocation": attempts, "range_resumes_this_invocation": range_resumes,
                    "completed_at_utc": utcnow(),
                }
                receipt_path.write_bytes(canonical_json(receipt))
                fsync_file(receipt_path)
                append_jsonl(state, {"event": "ACQUISITION_DOWNLOAD_COMPLETE", "accession": accession, "bytes": receipt["local_bytes"], "sha256": receipt["local_sha256"], "at": utcnow()})
                time.sleep(rate_delay)
                return complete_download, receipt
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            append_jsonl(failures, {"event": "ACQUISITION_ATTEMPT_FAILED", "accession": accession, "attempt": attempt, "type": type(exc).__name__, "message": str(exc), "at": utcnow()})
            if attempt == retries:
                raise GateError(f"download failed after {retries} attempts for {accession}: {exc}") from exc
            time.sleep(min(60.0, float(2 ** (attempt - 1))))
    raise GateError(f"download retry loop exhausted for {accession}")


def _normal_zip_path(raw: str) -> str:
    path = raw[2:] if raw.startswith("./") else raw
    if path.startswith("/") or ".." in Path(path).parts or "\\" in path:
        raise GateError(f"unsafe ZIP/checksum path: {raw}")
    return path


def _catalog_accessions(obj: Any) -> set[str]:
    found: set[str] = set()
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key.lower() in ("accession", "assemblyaccession", "assembly_accession") and isinstance(value, str) and ACCESSION_RE.fullmatch(value):
                found.add(value)
            found.update(_catalog_accessions(value))
    elif isinstance(obj, list):
        for value in obj:
            found.update(_catalog_accessions(value))
    return found


def validate_package(package: Path, accession: str) -> dict[str, Any]:
    if not zipfile.is_zipfile(package):
        raise GateError(f"source package is not ZIP: {accession}")
    with zipfile.ZipFile(package) as archive:
        infos = archive.infolist()
        if not infos or len(infos) > 10_000:
            raise GateError(f"source package entry count invalid: {accession}")
        names: list[str] = []
        total_uncompressed = 0
        for info in infos:
            name = _normal_zip_path(info.filename)
            names.append(name)
            total_uncompressed += info.file_size
            mode = info.external_attr >> 16
            if stat.S_ISLNK(mode):
                raise GateError(f"symlink in source archive: {name}")
        if total_uncompressed > 2_000_000_000:
            raise GateError(f"source package exceeds bounded uncompressed size: {accession}")
        bad = archive.testzip()
        if bad is not None:
            raise GateError(f"ZIP CRC/archive validation failed: {bad}")
        name_set = {name for name in names if not name.endswith("/")}
        md5_names = [name for name in name_set if Path(name).name == "md5sum.txt"]
        if len(md5_names) != 1:
            raise GateError(f"source package does not contain exactly one md5sum.txt: {accession}")
        md5_member = md5_names[0]
        md5_rows: dict[str, str] = {}
        for number, line in enumerate(archive.read(md5_member).decode("utf-8").splitlines(), 1):
            match = re.fullmatch(r"([0-9a-fA-F]{32})\s+\*?(.+)", line)
            if not match:
                raise GateError(f"malformed upstream md5sum line {number}: {accession}")
            rel = _normal_zip_path(match.group(2))
            if rel in md5_rows:
                raise GateError(f"duplicate upstream checksum path: {rel}")
            md5_rows[rel] = match.group(1).lower()
        listed = set(md5_rows)
        payload_files = name_set - {md5_member}
        missing_md5 = payload_files - listed
        extra_md5 = listed - payload_files
        # NCBI packages intentionally omit their generic top-level README from
        # md5sum.txt. All accession data, reports, and the catalog remain
        # upstream-checksummed; ZIP CRC plus the local archive SHA covers README.
        if extra_md5 or not missing_md5.issubset({"README.md"}):
            raise GateError(f"upstream MD5 coverage mismatch for {accession}: missing={sorted(missing_md5)} extra={sorted(extra_md5)}")
        member_sha256: dict[str, str] = {}
        member_bytes: dict[str, int] = {}
        for rel, expected in md5_rows.items():
            md5 = hashlib.md5(usedforsecurity=False)
            sha = hashlib.sha256()
            size = 0
            with archive.open(rel) as source:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    md5.update(block)
                    sha.update(block)
                    size += len(block)
            if md5.hexdigest() != expected:
                raise GateError(f"upstream MD5 mismatch: {accession}:{rel}")
            member_sha256[rel] = sha.hexdigest()
            member_bytes[rel] = size
        assembly_dirs = {parts[2] for name in payload_files if (parts := Path(name).parts)[:2] == ("ncbi_dataset", "data") and len(parts) >= 4 and ACCESSION_RE.fullmatch(parts[2])}
        if assembly_dirs != {accession}:
            raise GateError(f"package exact accession directory identity mismatch: expected {accession}, saw {sorted(assembly_dirs)}")
        fasta = [name for name in payload_files if name.startswith(f"ncbi_dataset/data/{accession}/") and name.endswith("_genomic.fna")]
        if len(fasta) != 1:
            raise GateError(f"package must contain exactly one genomic FASTA for {accession}")
        gff = [name for name in payload_files if name.startswith(f"ncbi_dataset/data/{accession}/") and name.endswith(("genomic.gff", "_genomic.gff"))]
        if len(gff) > 1:
            raise GateError(f"package has multiple genome GFF files: {accession}")
        sequence_report = [name for name in payload_files if name == f"ncbi_dataset/data/{accession}/sequence_report.jsonl"]
        catalogs = [name for name in payload_files if Path(name).name == "dataset_catalog.json"]
        if len(catalogs) != 1:
            raise GateError(f"package must contain exactly one dataset catalog: {accession}")
        catalog = json.loads(archive.read(catalogs[0]))
        catalog_accessions = _catalog_accessions(catalog)
        if catalog_accessions and catalog_accessions != {accession}:
            raise GateError(f"dataset catalog exact identity mismatch: {accession}: {sorted(catalog_accessions)}")
        reports = [name for name in payload_files if Path(name).name == "assembly_data_report.jsonl"]
        report_accessions: set[str] = set()
        if reports:
            for line in archive.read(reports[0]).splitlines():
                if line.strip():
                    report_accessions.update(_catalog_accessions(json.loads(line)))
        if report_accessions and accession not in report_accessions:
            raise GateError(f"assembly report exact identity mismatch: {accession}")
        return {
            "schema": "validated-ncbi-genome-package-v1", "accession": accession,
            "archive_validation": "PASS", "upstream_md5_validation": "PASS",
            "exact_accession_identity": "PASS", "package_bytes": package.stat().st_size,
            "package_sha256": sha_file(package), "zip_entries": len(infos),
            "uncompressed_member_bytes": total_uncompressed, "md5_member": md5_member,
            "md5_entries": len(md5_rows), "fasta_member": fasta[0],
            "fasta_bytes": member_bytes[fasta[0]], "fasta_sha256": member_sha256[fasta[0]],
            "gff_member": gff[0] if gff else ".", "gff_bytes": member_bytes[gff[0]] if gff else 0,
            "gff_sha256": member_sha256[gff[0]] if gff else ".",
            "sequence_report_member": sequence_report[0] if sequence_report else ".",
            "sequence_report_sha256": member_sha256[sequence_report[0]] if sequence_report else ".",
            "dataset_catalog_member": catalogs[0], "dataset_catalog_sha256": member_sha256[catalogs[0]],
            "assembly_report_accessions": sorted(report_accessions),
            "validated_at_utc": utcnow(),
        }


def commit_source_object(
    accession: str, source_root: Path, state: Path, failures: Path,
    retries: int, rate_delay: float, inject_kill: bool, inject_after_bytes: int,
) -> tuple[Path, dict[str, Any]]:
    final = source_root / accession
    if final.exists():
        verify_sha_inventory(final)
        manifest = json.loads((final / "manifest.json").read_text())
        if manifest.get("accession") != accession or manifest.get("state") != "COMPLETE":
            raise GateError(f"source object manifest mismatch: {accession}")
        return final, manifest
    stage = source_root / f".stage.{accession}"
    if stage.exists() and stage.is_symlink():
        raise GateError(f"source stage is symlink: {stage}")
    stage.mkdir(parents=True, exist_ok=True, mode=0o750)
    package, receipt = acquire_package(accession, stage, state, failures, retries, rate_delay, inject_kill, inject_after_bytes)
    validated = validate_package(package, accession)
    manifest = {
        "schema": "source-object-v1", "state": "COMPLETE", "accession": accession,
        "source_url": download_url(accession), "download_receipt": receipt,
        "validation": validated, "committed_at_utc": utcnow(),
    }
    (stage / "manifest.json").write_bytes(canonical_json(manifest))
    fsync_file(stage / "manifest.json")
    append_jsonl(state, {"event": "SOURCE_OBJECT_VALIDATED", "accession": accession, "package_sha256": validated["package_sha256"], "at": utcnow()})
    seal_directory(stage, final)
    append_jsonl(state, {"event": "SOURCE_OBJECT_COMMITTED", "accession": accession, "relpath": str(final), "at": utcnow()})
    return final, manifest


def percent_encode_contig(token: bytes, sample: str) -> tuple[str, str]:
    if not token:
        raise GateError("empty source FASTA identifier token")
    encoded = "".join(chr(byte) if byte in SAFE_CONTIG else f"%{byte:02X}" for byte in token)
    encoding = "IDENTITY_V1" if encoded.encode() == token else "PERCENT_UTF8_BYTES_V1"
    if len(f"{sample}#1#{encoded}".encode("ascii")) > 240:
        encoded = "CTGSHA256_" + sha_bytes(token)
        encoding = "SHA256_ALIAS_V1"
    if len(f"{sample}#1#{encoded}".encode("ascii")) > 240:
        raise GateError("PanSN identifier exceeds 240 bytes after digest alias")
    return encoded, encoding


def gff_escape(value: str) -> str:
    raw = value.encode("utf-8")
    return "".join(chr(byte) if byte in GFF_SAFE and byte != ord("%") else f"%{byte:02X}" for byte in raw)


def _all_strings(obj: Any, key_hint: str = "") -> Iterable[tuple[str, str]]:
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield from _all_strings(value, str(key))
    elif isinstance(obj, list):
        for value in obj:
            yield from _all_strings(value, key_hint)
    elif isinstance(obj, str):
        yield key_hint, obj


def load_sequence_reports(archive: zipfile.ZipFile, member: str) -> list[dict[str, Any]]:
    if member == ".":
        return []
    reports = []
    with archive.open(member) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.strip():
                try:
                    reports.append(json.loads(raw))
                except json.JSONDecodeError as exc:
                    raise GateError(f"invalid sequence report JSON line {number}") from exc
    return reports


def report_for_token(reports: list[dict[str, Any]], token: str) -> dict[str, Any] | None:
    matches = [report for report in reports if token in {value for _, value in _all_strings(report)}]
    if len(matches) > 1:
        raise GateError(f"ambiguous sequence report match for FASTA token {token}")
    return matches[0] if matches else None


def report_value(report: dict[str, Any] | None, keys: tuple[str, ...], default: str = ".") -> str:
    if report is None:
        return default
    wanted = {key.lower() for key in keys}
    for key, value in _all_strings(report):
        if key.lower() in wanted and value:
            return value
    return default


def report_aliases(report: dict[str, Any] | None) -> set[str]:
    if report is None:
        return set()
    aliases = set()
    for key, value in _all_strings(report):
        normalized = key.lower().replace("_", "")
        if normalized in ("refseqaccession", "genbankaccession", "sequencename") and SEQ_ACCESSION_RE.fullmatch(value):
            aliases.add(value)
    return aliases


def _b64(raw: bytes) -> str:
    return base64.b64encode(raw).decode("ascii")


def stream_canonical_fasta(
    source: BinaryIO, accession: str, bgzip_stdin: BinaryIO, canonical_digest: hashlib._Hash,
    reports: list[dict[str, Any]], canonical_relpath: str,
    inject_kill: bool, inject_after_bases: int, bgzip_process: subprocess.Popen[Any],
    part_handle: BinaryIO, state: Path,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen_tokens: set[bytes] = set()
    seen_names: set[str] = set()
    header: bytes | None = None
    token: bytes | None = None
    seq_hash: hashlib._Hash | None = None
    seq_length = 0
    wrap = bytearray()
    total_bases = 0
    first60 = bytearray()

    def emit(data: bytes) -> None:
        canonical_digest.update(data)
        bgzip_stdin.write(data)

    def finalize() -> None:
        nonlocal header, token, seq_hash, seq_length, wrap, first60
        if header is None or token is None or seq_hash is None:
            return
        if seq_length == 0:
            raise GateError(f"empty FASTA record in {accession}: {token!r}")
        if wrap:
            emit(bytes(wrap) + b"\n")
        token_text = token.decode("utf-8", "replace")
        encoded, encoding = percent_encode_contig(token, accession)
        name = f"{accession}#1#{encoded}"
        if name in seen_names:
            raise GateError(f"duplicate canonical FASTA ID: {name}")
        seen_names.add(name)
        report = report_for_token(reports, token_text)
        match = SEQ_ACCESSION_RE.fullmatch(token_text)
        refseq = report_value(report, ("refseqAccession", "refseq_accession"))
        genbank = report_value(report, ("genbankAccession", "genbank_accession"))
        role = report_value(report, ("role", "sequenceRole", "sequence_role"), "unknown")
        assigned = report_value(report, ("chrName", "assignedMolecule", "assigned_molecule"), ".")
        plasmid = report_value(report, ("plasmidName", "plasmid_name"), ".")
        topology = report_value(report, ("topology",), "unknown").lower()
        if topology not in ("linear", "circular", "unknown"):
            topology = "unknown"
        sequence_name = report_value(report, ("sequenceName", "sequence_name"), ".")
        rows.append({
            "assembly_id": accession, "accession": accession, "contig_order": len(rows) + 1,
            "source_fasta_header_b64": _b64(b">" + header), "source_fasta_id_token_b64": _b64(token),
            "source_contig_id_display": token_text, "source_contig_accession": match.group(1) if match else ".",
            "source_contig_version": match.group(2) if match else ".",
            "source_contig_accession_version": token_text if match else ".",
            "genbank_contig_accession_version": genbank, "refseq_contig_accession_version": refseq,
            "assembly_report_sequence_name": sequence_name, "contig_id_encoding": encoding,
            "pansn_sample": accession, "pansn_haplotype": "1", "pansn_contig": encoded,
            "pansn_sequence_name": name, "fasta_seqid": name, "replicon_role": role,
            "assigned_molecule": assigned, "plasmid_name_b64": _b64(plasmid.encode()) if plasmid != "." else ".",
            "topology": topology, "contig_length": seq_length, "contig_sequence_sha256": seq_hash.hexdigest(),
            "canonical_bgzf_relpath": canonical_relpath, "sequence_prefix_hex": bytes(first60).hex(),
            "sequence_report_aliases": sorted(report_aliases(report)),
        })
        header = None
        token = None
        seq_hash = None
        seq_length = 0
        wrap = bytearray()
        first60 = bytearray()

    for raw in source:
        if b"\r" in raw:
            raise GateError(f"CR byte in source FASTA: {accession}")
        line = raw[:-1] if raw.endswith(b"\n") else raw
        if line.startswith(b">"):
            finalize()
            header = line[1:]
            parts = header.split(None, 1)
            if not parts:
                raise GateError(f"empty source FASTA header: {accession}")
            token = parts[0]
            if token in seen_tokens:
                raise GateError(f"duplicate source FASTA token: {accession}:{token!r}")
            seen_tokens.add(token)
            encoded, _ = percent_encode_contig(token, accession)
            emit(f">{accession}#1#{encoded}\n".encode("ascii"))
            seq_hash = hashlib.sha256()
            continue
        if header is None or token is None or seq_hash is None:
            if not line:
                continue
            raise GateError(f"sequence before first FASTA header: {accession}")
        if not line:
            raise GateError(f"blank line in FASTA sequence: {accession}:{token!r}")
        if any(byte not in DNA for byte in line):
            raise GateError(f"invalid/lowercase FASTA base in {accession}:{token!r}")
        seq_hash.update(line)
        seq_length += len(line)
        total_bases += len(line)
        if len(first60) < 60:
            first60.extend(line[: 60 - len(first60)])
        wrap.extend(line)
        while len(wrap) >= 60:
            emit(bytes(wrap[:60]) + b"\n")
            del wrap[:60]
        if inject_kill and total_bases >= inject_after_bases:
            bgzip_stdin.flush()
            bgzip_stdin.close()
            bgzip_process.kill()
            bgzip_process.wait()
            part_handle.flush()
            os.fsync(part_handle.fileno())
            append_jsonl(state, {"event": "INJECTED_CONVERSION_SIGKILL", "accession": accession, "bases_streamed": total_bases, "at": utcnow()})
            os.kill(os.getpid(), signal.SIGKILL)
    finalize()
    if not rows:
        raise GateError(f"source FASTA has no records: {accession}")
    return rows


def validate_gff_aliases(
    archive: zipfile.ZipFile, gff_member: str, contigs: list[dict[str, Any]], accession: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if gff_member == ".":
        return [], {"status": "NOT_AVAILABLE", "coordinate_policy": "NOT_APPLICABLE_NO_GFF", "feature_rows": 0}
    aliases: dict[str, dict[str, Any]] = {}
    for contig in contigs:
        source = contig["source_contig_id_display"]
        candidates = {source, *contig.get("sequence_report_aliases", [])}
        for candidate in candidates:
            prior = aliases.get(candidate)
            if prior is not None and prior["pansn_sequence_name"] != contig["pansn_sequence_name"]:
                raise AnnotationValidationError(f"sequence-report alias collision: {accession}:{candidate}")
            aliases[candidate] = contig
    used: dict[str, dict[str, Any]] = {}
    feature_rows = 0
    with archive.open(gff_member) as handle:
        for line_number, raw in enumerate(handle, 1):
            if b"\r" in raw:
                raise AnnotationValidationError(f"CR byte in GFF: {accession}:{line_number}")
            if raw.startswith(b"##FASTA"):
                raise AnnotationValidationError(f"unexpected embedded FASTA in source GFF: {accession}")
            if not raw or raw.startswith(b"#") or not raw.strip():
                continue
            columns = raw.rstrip(b"\n").split(b"\t")
            if len(columns) != 9:
                raise AnnotationValidationError(f"GFF row does not have 9 columns: {accession}:{line_number}")
            try:
                lexical = columns[0].decode("ascii")
                decoded_bytes = urllib.parse.unquote_to_bytes(lexical)
                decoded = decoded_bytes.decode("utf-8")
            except (UnicodeDecodeError, ValueError) as exc:
                raise AnnotationValidationError(f"invalid GFF seqid encoding: {accession}:{line_number}") from exc
            contig = aliases.get(decoded)
            if contig is None:
                raise AnnotationValidationError(f"GFF seqid does not resolve to FASTA/sequence-report alias: {accession}:{decoded}")
            try:
                start, end = int(columns[3]), int(columns[4])
            except ValueError as exc:
                raise AnnotationValidationError(f"non-integer GFF coordinate: {accession}:{line_number}") from exc
            if not (1 <= start <= end <= int(contig["contig_length"])):
                raise AnnotationValidationError(f"out-of-range GFF coordinate: {accession}:{line_number}:{start}-{end}")
            feature_rows += 1
            used.setdefault(lexical, {
                "accession": accession, "source_gff_member": gff_member,
                "source_gff_seqid_lexical_b64": _b64(columns[0]),
                "source_gff_seqid_decoded_b64": _b64(decoded_bytes),
                "source_fasta_id_token_b64": contig["source_fasta_id_token_b64"],
                "pansn_sequence_name": contig["pansn_sequence_name"],
                "canonical_gff_seqid_lexical": gff_escape(contig["pansn_sequence_name"]),
                "canonical_gff_seqid_decoded": contig["pansn_sequence_name"],
                "source_coordinate_convention": "GFF3_1_BASED_CLOSED",
            })
    if feature_rows == 0:
        raise AnnotationValidationError(f"present GFF contains no feature rows: {accession}")
    rows = [used[key] | {"alias_order": number} for number, key in enumerate(sorted(used), 1)]
    return rows, {"status": "ALIASES_VALIDATED_NO_TRANSFORMED_GFF", "coordinate_policy": "GFF3_1_BASED_CLOSED_PASS", "feature_rows": feature_rows, "distinct_seqids": len(rows)}


def _parse_canonical_bgzf(path: Path) -> tuple[list[dict[str, Any]], str]:
    digest = hashlib.sha256()
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    current: str | None = None
    seq_hash: hashlib._Hash | None = None
    length = 0
    prefix = bytearray()
    with gzip.open(path, "rb") as handle:
        for raw in handle:
            digest.update(raw)
            if b"\r" in raw or not raw.endswith(b"\n"):
                raise GateError(f"canonical FASTA is not LF line-terminated: {path}")
            line = raw[:-1]
            if line.startswith(b">"):
                if current is not None and seq_hash is not None:
                    rows.append({"name": current, "length": length, "sha256": seq_hash.hexdigest(), "prefix": bytes(prefix)})
                try:
                    current = line[1:].decode("ascii")
                except UnicodeDecodeError as exc:
                    raise GateError("canonical FASTA ID is not ASCII") from exc
                if current in seen or len(current.encode()) > 240:
                    raise GateError(f"duplicate/overlength canonical FASTA ID: {current}")
                fields = current.split("#")
                if len(fields) != 3 or not ACCESSION_RE.fullmatch(fields[0]) or fields[1] != "1" or not fields[2]:
                    raise GateError(f"invalid PanSN grammar: {current}")
                seen.add(current)
                seq_hash = hashlib.sha256()
                length = 0
                prefix = bytearray()
            else:
                if current is None or seq_hash is None or not line or len(line) > 60 or any(byte not in DNA for byte in line):
                    raise GateError(f"invalid canonical FASTA sequence line: {path}")
                seq_hash.update(line)
                length += len(line)
                if len(prefix) < 60:
                    prefix.extend(line[: 60 - len(prefix)])
        if current is not None and seq_hash is not None:
            rows.append({"name": current, "length": length, "sha256": seq_hash.hexdigest(), "prefix": bytes(prefix)})
    return rows, digest.hexdigest()


def _parse_fai(path: Path) -> list[tuple[str, int]]:
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        columns = line.split("\t")
        if len(columns) < 5:
            raise GateError(f"malformed FAI line {number}: {path}")
        rows.append((columns[0], int(columns[1])))
    return rows


def _validate_gzi(path: Path) -> None:
    data = path.read_bytes()
    if len(data) < 8 or (len(data) - 8) % 16:
        raise GateError(f"malformed GZI structure: {path}")
    count = struct.unpack("<Q", data[:8])[0]
    if len(data) != 8 + 16 * count:
        raise GateError(f"GZI count/byte mismatch: {path}")


def canonicalize_object(
    accession: str, stage_order: int, source_object: Path, source_manifest: dict[str, Any],
    canonical_root: Path, state: Path, bgzip: str, samtools: str, threads: int, level: int,
    inject_kill: bool, inject_after_bases: int,
) -> tuple[Path, dict[str, Any], list[dict[str, Any]]]:
    final = canonical_root / accession
    if final.exists():
        verify_sha_inventory(final)
        manifest = json.loads((final / "manifest.json").read_text())
        rows = read_tsv(final / "contigs.tsv", verify_hashes=True)
        if manifest.get("accession") != accession or manifest.get("state") != "COMPLETE":
            raise GateError(f"canonical object manifest mismatch: {accession}")
        return final, manifest, rows
    stage = canonical_root / f".stage.{accession}"
    if stage.exists():
        if stage.is_symlink():
            raise GateError(f"canonical stage is symlink: {stage}")
        append_jsonl(state, {"event": "INTERRUPTED_CONVERSION_STAGE_DISCARDED", "accession": accession, "at": utcnow()})
        shutil.rmtree(stage)
    stage.mkdir(parents=True, mode=0o750)
    package = source_object / "package.zip"
    validated = source_manifest["validation"]
    fasta_member = validated["fasta_member"]
    gff_member = validated["gff_member"]
    seq_report_member = validated["sequence_report_member"]
    basename = f"{accession}.pansn.fa.gz"
    bgzf_path = stage / basename
    stderr_path = stage / "bgzip.stderr.log"
    canonical_hash = hashlib.sha256()
    canonical_relpath = f"canonical_objects/{accession}/{basename}"
    with zipfile.ZipFile(package) as archive:
        reports = load_sequence_reports(archive, seq_report_member)
        with archive.open(fasta_member) as source, bgzf_path.open("wb") as output, stderr_path.open("wb") as error:
            process = subprocess.Popen(
                [bgzip, "-@", str(threads), "-l", str(level), "--binary", "-c"],
                stdin=subprocess.PIPE, stdout=output, stderr=error,
            )
            assert process.stdin is not None
            try:
                contigs = stream_canonical_fasta(
                    source, accession, process.stdin, canonical_hash, reports, canonical_relpath,
                    inject_kill, inject_after_bases, process, output, state,
                )
                process.stdin.close()
                status = process.wait()
            except BaseException:
                if process.poll() is None:
                    process.kill()
                    process.wait()
                raise
            output.flush()
            os.fsync(output.fileno())
        if status != 0:
            raise GateError(f"bgzip failed for {accession} with exit {status}")
        try:
            aliases, annotation = validate_gff_aliases(archive, gff_member, contigs, accession)
        except AnnotationValidationError as exc:
            aliases = []
            annotation = {
                "status": "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW",
                "coordinate_policy": "SOURCE_GFF_VALIDATION_FAILED_NO_VIEW",
                "feature_rows": 0,
                "failure_reason": str(exc),
            }
            append_jsonl(state, {"event": "ANNOTATION_ALIAS_VIEW_QUARANTINED", "accession": accession, "reason": str(exc), "at": utcnow()})
    if subprocess.run([bgzip, "-t", str(bgzf_path)], capture_output=True).returncode != 0:
        raise GateError(f"bgzip integrity test failed: {accession}")
    subprocess.run([samtools, "faidx", str(bgzf_path)], check=True, capture_output=True)
    fai = Path(str(bgzf_path) + ".fai")
    gzi = Path(str(bgzf_path) + ".gzi")
    if not fai.is_file() or not gzi.is_file():
        raise GateError(f"samtools faidx did not create both FAI/GZI: {accession}")
    _validate_gzi(gzi)
    parsed, content_sha = _parse_canonical_bgzf(bgzf_path)
    expected = [(row["pansn_sequence_name"], int(row["contig_length"]), row["contig_sequence_sha256"]) for row in contigs]
    observed = [(row["name"], row["length"], row["sha256"]) for row in parsed]
    if observed != expected:
        raise GateError(f"rename-only sequence digest/length/order mismatch: {accession}")
    if content_sha != canonical_hash.hexdigest():
        raise GateError(f"canonical content digest mismatch: {accession}")
    if _parse_fai(fai) != [(name, length) for name, length, _ in expected]:
        raise GateError(f"FAI name/length/order mismatch: {accession}")
    for row in parsed:
        end = min(60, row["length"])
        region = f"{row['name']}:1-{end}"
        result = subprocess.run([samtools, "faidx", str(bgzf_path), region], check=True, capture_output=True)
        lines = result.stdout.splitlines()
        bases = b"".join(lines[1:])
        if bases != row["prefix"][:end]:
            raise GateError(f"samtools quoted PanSN region round-trip mismatch: {accession}:{row['name']}")
    contig_fields = CONTIG_FIELDS
    for row in contigs:
        row["stage_b_order"] = stage_order
        row["source_fasta_member"] = fasta_member
    write_tsv(stage / "contigs.tsv", contig_fields, contigs)
    alias_fields = [
        "alias_order", "accession", "source_gff_member", "source_gff_seqid_lexical_b64",
        "source_gff_seqid_decoded_b64", "source_fasta_id_token_b64", "pansn_sequence_name",
        "canonical_gff_seqid_lexical", "canonical_gff_seqid_decoded", "source_coordinate_convention", "row_sha256",
    ]
    write_tsv(stage / "annotation_aliases.tsv", alias_fields, aliases)
    tool_command = [bgzip, "-@", str(threads), "-l", str(level), "--binary", "-c"]
    manifest = {
        "schema": "canonical-assembly-object-v1", "state": "COMPLETE", "accession": accession,
        "assembly_id": accession, "pansn_policy_version": POLICY, "pansn_sample": accession,
        "pansn_haplotype": "1", "sample_id_basis": "exact_versioned_assembly_accession",
        "haplotype_basis": "nominally_haploid_E_coli", "source_object_relpath": f"source_objects/{accession}",
        "source_package_sha256": validated["package_sha256"], "source_fasta_member": fasta_member,
        "source_fasta_bytes": validated["fasta_bytes"], "source_decompressed_sha256": validated["fasta_sha256"],
        "source_gff_member": gff_member, "source_gff_sha256": validated["gff_sha256"],
        "annotation": annotation, "annotation_view": "validated alias table only; no transformed GFF published",
        "canonical_bgzf_relpath": canonical_relpath, "canonical_bgzf_bytes": bgzf_path.stat().st_size,
        "canonical_fasta_content_sha256": content_sha, "canonical_bgzf_sha256": sha_file(bgzf_path),
        "fai_relpath": canonical_relpath + ".fai", "fai_bytes": fai.stat().st_size, "fai_sha256": sha_file(fai),
        "gzi_relpath": canonical_relpath + ".gzi", "gzi_bytes": gzi.stat().st_size, "gzi_sha256": sha_file(gzi),
        "crosswalk_sha256": sha_file(stage / "contigs.tsv"), "annotation_aliases_sha256": sha_file(stage / "annotation_aliases.tsv"),
        "contig_count": len(contigs), "total_bases": sum(int(row["contig_length"]) for row in contigs),
        "rename_only_validation": "PASS", "bgzf_integrity": "PASS", "index_name_roundtrip": "PASS",
        "unique_pansn_names": "PASS", "gff_source_coordinate_policy": annotation["coordinate_policy"],
        "bgzip_argv": tool_command, "transformation_command_sha256": sha_bytes(canonical_json({"policy": POLICY, "argv": tool_command})),
        "completed_at_utc": utcnow(),
    }
    (stage / "manifest.json").write_bytes(canonical_json(manifest))
    fsync_file(stage / "manifest.json")
    append_jsonl(state, {"event": "CANONICAL_OBJECT_VALIDATED", "accession": accession, "contigs": len(contigs), "bases": manifest["total_bases"], "at": utcnow()})
    seal_directory(stage, final)
    append_jsonl(state, {"event": "CANONICAL_OBJECT_COMMITTED", "accession": accession, "relpath": str(final), "at": utcnow()})
    return final, manifest, contigs


def executable_record(name: str) -> dict[str, Any]:
    path = shutil.which(name)
    if path is None:
        raise GateError(f"required executable missing: {name}")
    real = Path(path).resolve()
    version_cmd = [path, "--version"] if name == "bgzip" else [path, "--version"]
    result = subprocess.run(version_cmd, capture_output=True)
    raw_version = result.stdout or result.stderr
    version = raw_version.decode("utf-8", "replace").splitlines()[0] if raw_version else "."
    return {"name": name, "path": str(real), "sha256": sha_file(real), "version_first_line": version, "version_exit": result.returncode}


def directory_usage(path: Path) -> tuple[int, int]:
    total = 0
    files = 0
    if path.exists():
        for item in path.rglob("*"):
            if item.is_file() and not item.is_symlink():
                total += item.stat().st_size
                files += 1
    return total, files


def audit_global_release_cap(project_root: Path, allowed_collection: set[str]) -> dict[str, Any]:
    union: set[str] = set()
    scanned = []
    releases_root = project_root / "releases"
    if releases_root.exists():
        for release_json in sorted(releases_root.glob("*/*/release.json")):
            root = release_json.parent
            if not (root / "COMPLETE").is_file():
                continue
            obj = json.loads(release_json.read_text())
            count = int(obj.get("counts", {}).get("distinct_sequence_bearing_assemblies", 0))
            accessions = set(obj.get("sequence_bearing_assembly_accessions", []))
            if count and len(accessions) != count:
                raise GateError(f"committed sequence-bearing release lacks exact accession inventory: {release_json}")
            union.update(accessions)
            scanned.append({"release_json": str(release_json), "declared_count": count, "accessions": sorted(accessions)})
    if len(union) > 1000 or not union.issubset(allowed_collection):
        raise GateError(f"global distinct assembly cap/subset gate failed: count={len(union)}")
    return {"verdict": "PASS", "distinct_exact_assembly_revisions": len(union), "cap": 1000, "subset_of_frozen_collection": True, "accessions": sorted(union), "scanned_releases": scanned, "captured_at_utc": utcnow()}


def _load_allowed_collection(repo: Path) -> set[str]:
    path = repo / "manifests/collection-v1/assemblies.tsv.gz"
    allowed = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            allowed.add(row["resolved_assembly_accession_version"])
    if len(allowed) != 26077:
        raise GateError("frozen collection cardinality mismatch while auditing global cap")
    return allowed


def build_release_tables(stage: Path, input_rows: list[dict[str, str]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    assembly_rows: list[dict[str, Any]] = []
    contig_rows: list[dict[str, Any]] = []
    checksum_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    for input_row in input_rows:
        order = int(input_row["stage_b_order"])
        accession = input_row["resolved_assembly_accession_version"]
        source_dir = stage / "source_objects" / accession
        canonical_dir = stage / "canonical_objects" / accession
        verify_sha_inventory(source_dir)
        verify_sha_inventory(canonical_dir)
        source = json.loads((source_dir / "manifest.json").read_text())
        canonical = json.loads((canonical_dir / "manifest.json").read_text())
        remote = source["download_receipt"]["remote_identity"]
        headers = remote.get("headers", {})
        row = {
            "stage_b_order": order, "assembly_id": accession, "accession": accession,
            "predecessor_resolution_status": input_row["resolution_status"], "terminal_state": "VALIDATED",
            "source_object_relpath": f"source_objects/{accession}",
            "source_package_bytes": source["validation"]["package_bytes"], "source_package_sha256": source["validation"]["package_sha256"],
            "source_fasta_member": source["validation"]["fasta_member"], "source_fasta_bytes": source["validation"]["fasta_bytes"],
            "source_decompressed_sha256": source["validation"]["fasta_sha256"],
            "source_gff_member": source["validation"]["gff_member"], "source_gff_sha256": source["validation"]["gff_sha256"],
            "annotation_status": canonical["annotation"]["status"], "canonical_object_relpath": f"canonical_objects/{accession}",
            "canonical_bgzf_relpath": canonical["canonical_bgzf_relpath"], "canonical_bgzf_bytes": canonical["canonical_bgzf_bytes"],
            "canonical_fasta_content_sha256": canonical["canonical_fasta_content_sha256"], "canonical_bgzf_sha256": canonical["canonical_bgzf_sha256"],
            "fai_sha256": canonical["fai_sha256"], "gzi_sha256": canonical["gzi_sha256"],
            "crosswalk_sha256": canonical["crosswalk_sha256"], "annotation_aliases_sha256": canonical["annotation_aliases_sha256"],
            "contig_count": canonical["contig_count"], "total_bases": canonical["total_bases"],
            "remote_identity_strength": remote["identity_strength"], "remote_etag": headers.get("etag", "."),
            "remote_last_modified": headers.get("last_modified", "."),
            "download_attempts": source["download_receipt"]["attempts_this_invocation"],
            "range_resumes": source["download_receipt"]["range_resumes_this_invocation"],
            "validated_at_utc": canonical["completed_at_utc"],
        }
        assembly_rows.append(row)
        per_contigs = read_tsv(canonical_dir / "contigs.tsv", CONTIG_FIELDS, verify_hashes=True)
        contig_rows.extend(per_contigs)
        for role, rel in (
            ("source_package", f"source_objects/{accession}/package.zip"),
            ("source_manifest", f"source_objects/{accession}/manifest.json"),
            ("canonical_bgzf", canonical["canonical_bgzf_relpath"]),
            ("fai", canonical["fai_relpath"]), ("gzi", canonical["gzi_relpath"]),
            ("contig_crosswalk", f"canonical_objects/{accession}/contigs.tsv"),
            ("annotation_aliases", f"canonical_objects/{accession}/annotation_aliases.tsv"),
            ("canonical_manifest", f"canonical_objects/{accession}/manifest.json"),
        ):
            path = stage / rel
            checksum_rows.append({"stage_b_order": order, "accession": accession, "artifact_role": role, "relative_path": rel, "bytes": path.stat().st_size, "sha256": sha_file(path)})
        state_rows.append({"stage_b_order": order, "accession": accession, "source_state": "COMPLETE", "canonical_state": "COMPLETE", "terminal_state": "VALIDATED", "reason": "."})
    return assembly_rows, contig_rows, checksum_rows, state_rows


def publish_tracked(external: Path, tracked: Path, artifact_dir: Path) -> None:
    tracked.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    allowed = {
        "stage_b_10.tsv", "assemblies.tsv", "contigs.tsv.gz", "checksums.tsv", "state.tsv",
        "release.json", "external_SHA256SUMS", "SHA256SUMS",
    }
    for old in tracked.iterdir():
        if old.is_file() and old.name not in allowed:
            raise GateError(f"unexpected file in tracked release directory: {old}")
    shutil.copyfile(external / "manifests/stage_b_10.tsv", tracked / "stage_b_10.tsv")
    shutil.copyfile(external / "manifests/assemblies.tsv", tracked / "assemblies.tsv")
    deterministic_gzip(external / "manifests/contigs.tsv", tracked / "contigs.tsv.gz")
    shutil.copyfile(external / "manifests/checksums.tsv", tracked / "checksums.tsv")
    shutil.copyfile(external / "manifests/state.tsv", tracked / "state.tsv")
    shutil.copyfile(external / "release.json", tracked / "release.json")
    shutil.copyfile(external / "SHA256SUMS", tracked / "external_SHA256SUMS")
    entries = []
    for path in sorted(p for p in tracked.iterdir() if p.is_file() and p.name != "SHA256SUMS"):
        entries.append(f"{sha_file(path)}  {path.name}\n")
    (tracked / "SHA256SUMS").write_text("".join(entries))
    fsync_file(tracked / "SHA256SUMS")
    shutil.copyfile(external / "resource_summary.json", artifact_dir / "resource_summary.json")
    shutil.copyfile(external / "restart_evidence.json", artifact_dir / "restart_evidence.json")
    shutil.copyfile(external / "tools.json", artifact_dir / "tool_versions.json")
    shutil.copyfile(external / "manifests/state.tsv", artifact_dir / "state_summary.tsv")


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo = Path(args.repo_root).resolve()
    durable_task = Path(args.durable_task_root).resolve()
    scratch = Path(args.scratch_root).resolve()
    tracked = Path(args.tracked_root).resolve()
    artifact_dir = Path(args.artifact_root).resolve()
    allocations = Allocations(
        args.assigned_ram_bytes, args.durable_allocation_bytes, args.scratch_allocation_bytes,
        args.inode_allocation, args.predicted_durable_peak_bytes, args.predicted_scratch_peak_bytes,
        args.predicted_files, args.unfinished_write_bytes,
    )
    root_start = verify_root_inputs(repo)
    input_rows, predecessor = verify_predecessor(repo)
    allowed_collection = _load_allowed_collection(repo)
    global_start = audit_global_release_cap(durable_task.parents[1], allowed_collection)
    projected_union = set(global_start["accessions"]) | set(EXPECTED_ACCESSIONS)
    if len(projected_union) > 1000 or not projected_union.issubset(allowed_collection):
        raise GateError("projected finish global distinct-assembly cap/subset gate failed")
    global_cap_evidence = {
        "schema": "global-distinct-assembly-cap-v1", "verdict": "PASS",
        "start": global_start, "projected_finish_distinct_exact_assembly_revisions": len(projected_union),
        "projected_finish_accessions": sorted(projected_union), "cap": 1000,
        "projected_finish_subset_of_frozen_collection": True,
    }
    rid = release_id()
    final = durable_task / rid
    stage = durable_task / f".stage.{rid}.{args.run_id}"
    initial_preflight = live_preflight(final, scratch, allocations, "INITIAL")
    mkdir_owned(durable_task)
    mkdir_owned(scratch)
    if final.exists():
        verify_sha_inventory(final)
        existing = json.loads((final / "release.json").read_text())
        if existing.get("release_id") != rid or existing.get("verdict") != "PASS":
            raise GateError("existing release identity/verdict mismatch")
        publish_tracked(final, tracked, artifact_dir)
        return existing
    if stage.exists():
        if (stage / "COMPLETE").exists():
            raise GateError("interrupted overall stage contains COMPLETE but was not promoted")
        seed_path = stage / "input_manifest.json"
        if not seed_path.is_file() or seed_path.read_bytes() != canonical_json(release_seed()):
            raise GateError("interrupted overall stage input seed mismatch")
        resumed = True
    else:
        stage.mkdir(mode=0o750)
        resumed = False
    for directory in (stage / "manifests", stage / "source_objects", stage / "canonical_objects", stage / "logs"):
        directory.mkdir(exist_ok=True, mode=0o750)
    state = stage / "state.jsonl"
    failures = stage / "failures.jsonl"
    failures.touch(exist_ok=True)
    (stage / "input_manifest.json").write_bytes(canonical_json(release_seed()))
    shutil.copyfile(repo / "manifests/collection-v1/stage_b_10.tsv", stage / "manifests/stage_b_10.tsv")
    resources = stage / "resources.jsonl"
    append_jsonl(resources, initial_preflight)
    append_jsonl(state, {"event": "RESUME_PREFLIGHT_PASS" if resumed else "PREFLIGHT_PASS", "release_id": rid, "at": utcnow()})
    tools = {"bgzip": executable_record(args.bgzip), "samtools": executable_record(args.samtools)}
    (stage / "tools.json").write_bytes(canonical_json(tools))
    provenance = {
        "schema": "acquisition-canonicalization-provenance-v1", "task_id": TASK_ID,
        "release_id": rid, "run_id": args.run_id, "argv": sys.argv, "python": sys.version,
        "platform": platform.platform(), "hostname": socket.gethostname(), "pid": os.getpid(),
        "uid": os.getuid(), "gid": os.getgid(), "environment": {key: os.environ.get(key, ".") for key in ("USER", "LANG", "WG_TASK_ID", "WG_AGENT_ID", "WG_MODEL", "WG_TIER", "PI_MODEL", "PI_PROVIDER")},
        "source_api": "NCBI Datasets v2", "source_include_types": list(INCLUDE_TYPES),
        "max_distinct_exact_assembly_revisions": 10, "global_cap": 1000,
        "collection_release_id": predecessor["release_id"], "collection_release_json_sha256": COLLECTION_RELEASE_SHA256,
        "stage_b_manifest_sha256": STAGE_B_SHA256, "created_at_utc": utcnow(),
    }
    (stage / "provenance.json").write_bytes(canonical_json(provenance))
    (stage / "global_cap_evidence.json").write_bytes(canonical_json(global_cap_evidence))
    start_swap = initial_preflight["swap_free_bytes"]
    start_ru = resource.getrusage(resource.RUSAGE_SELF)
    max_disk = directory_usage(stage)
    try:
        for row in input_rows:
            accession = row["resolved_assembly_accession_version"]
            order = int(row["stage_b_order"])
            acquisition_preflight = live_preflight(final, scratch, allocations, f"ACQUISITION_{order:02d}_{accession}")
            append_jsonl(resources, acquisition_preflight)
            append_jsonl(state, {"event": "ASSEMBLY_ATTEMPT_STARTED", "stage_b_order": order, "accession": accession, "at": utcnow()})
            source_final, source_manifest = commit_source_object(
                accession, stage / "source_objects", state, failures, args.retries, args.rate_delay,
                args.inject_kill == "acquisition" and accession == args.inject_accession,
                args.inject_after_bytes,
            )
            conversion_preflight = live_preflight(final, scratch, allocations, f"CANONICALIZATION_{order:02d}_{accession}")
            append_jsonl(resources, conversion_preflight)
            canonicalize_object(
                accession, order, source_final, source_manifest, stage / "canonical_objects", state,
                args.bgzip, args.samtools, args.bgzip_threads, args.bgzip_level,
                args.inject_kill == "conversion" and accession == args.inject_accession,
                args.inject_after_bases,
            )
            usage_now = directory_usage(stage)
            max_disk = (max(max_disk[0], usage_now[0]), max(max_disk[1], usage_now[1]))
            append_jsonl(state, {"event": "ASSEMBLY_TERMINAL_VALIDATED", "stage_b_order": order, "accession": accession, "at": utcnow()})
        assembly_rows, contig_rows, checksum_rows, state_rows = build_release_tables(stage, input_rows)
        if len(assembly_rows) != 10 or {row["accession"] for row in assembly_rows} != set(EXPECTED_ACCESSIONS):
            raise GateError("100% assembly row accounting gate failed")
        if len(contig_rows) != sum(int(row["contig_count"]) for row in assembly_rows):
            raise GateError("100% contig row accounting gate failed")
        all_names = [row["pansn_sequence_name"] for row in contig_rows]
        if len(all_names) != len(set(all_names)):
            raise GateError("cohort-wide PanSN name uniqueness gate failed")
        write_tsv(stage / "manifests/assemblies.tsv", ASSEMBLY_FIELDS, assembly_rows)
        write_tsv(stage / "manifests/contigs.tsv", CONTIG_FIELDS, contig_rows)
        write_tsv(stage / "manifests/checksums.tsv", CHECKSUM_FIELDS, checksum_rows)
        write_tsv(stage / "manifests/state.tsv", STATE_FIELDS, state_rows)
        root_finish = verify_root_inputs(repo)
        if root_finish != root_start:
            raise GateError("root inputs changed during task")
        promotion_preflight = live_preflight(final, scratch, allocations, "PROMOTION")
        append_jsonl(resources, promotion_preflight)
        end_ru = resource.getrusage(resource.RUSAGE_SELF)
        peak_rss = end_ru.ru_maxrss * 1024
        swap_end = promotion_preflight["swap_free_bytes"]
        swap_growth = max(0, start_swap - swap_end)
        process_swaps = max(0, end_ru.ru_nswap - start_ru.ru_nswap)
        usage_now = directory_usage(stage)
        max_disk = (max(max_disk[0], usage_now[0]), max(max_disk[1], usage_now[1]))
        resource_checks = {
            "peak_rss_le_70pct_assigned": peak_rss * 100 <= allocations.assigned_ram_bytes * 70,
            "process_swap_events_zero": process_swaps == 0,
            "system_swap_growth_zero": swap_growth == 0,
            "measured_durable_peak_le_70pct_allocation": max_disk[0] * 100 <= allocations.durable_allocation_bytes * 70,
            "measured_files_le_50pct_inode_allocation": max_disk[1] * 2 <= allocations.inode_allocation,
        }
        if not all(resource_checks.values()):
            raise GateError("end resource gate NO_GO: " + json.dumps(resource_checks, sort_keys=True))
        resource_summary = {
            "schema": "resource-summary-v1", "verdict": "PASS", "allocations": asdict(allocations),
            "peak_rss_bytes": peak_rss, "peak_rss_fraction": peak_rss / allocations.assigned_ram_bytes,
            "process_swap_events": process_swaps, "system_swap_growth_bytes": swap_growth,
            "measured_release_stage_peak_bytes": max_disk[0], "measured_release_stage_peak_files": max_disk[1],
            "checks": resource_checks, "preflight_record_count": sum(1 for _ in resources.open("rb")),
            "start": initial_preflight, "finish": promotion_preflight,
        }
        (stage / "resource_summary.json").write_bytes(canonical_json(resource_summary))
        state_text = state.read_text()
        restart_evidence = {
            "schema": "restart-evidence-v1",
            "acquisition_injected_kill_observed": "INJECTED_ACQUISITION_SIGKILL" in state_text,
            "acquisition_partial_safe_restart_observed": any(event in state_text for event in ("ACQUISITION_RANGE_RESUME", "ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE")),
            "conversion_injected_kill_observed": "INJECTED_CONVERSION_SIGKILL" in state_text,
            "conversion_partial_discard_restart_observed": "INTERRUPTED_CONVERSION_STAGE_DISCARDED" in state_text,
            "no_partial_release_before_promotion": True,
        }
        if not all(restart_evidence[key] for key in restart_evidence if key not in ("schema",)):
            raise GateError("required injected kill/restart evidence is incomplete")
        (stage / "restart_evidence.json").write_bytes(canonical_json(restart_evidence))
        gates = {
            "predecessor_release_and_manifest": "PASS", "source_immutability": "PASS",
            "accession_version_identity": "PASS", "global_distinct_assembly_cap": "PASS",
            "resource": "PASS", "remote_identity_safe_resume": "PASS", "archive_integrity": "PASS",
            "upstream_md5_and_local_sha256": "PASS", "row_accounting": "PASS",
            "rename_only_sequence_identity": "PASS", "bgzf_index_name_roundtrip": "PASS",
            "pansn_uniqueness_and_reversibility": "PASS", "source_coordinate_policy": "PASS",
            "annotation_alias_policy": "PASS", "pinned_bgzip_samtools_compatibility": "PASS",
            "injected_kill_restart": "PASS", "atomic_object_and_release_promotion": "PASS",
            "deterministic_resume_semantics": "PASS", "scale_trend": "NOT_APPLICABLE_STAGE_B_NON_SCALE_BEARING",
        }
        release = {
            "schema_version": SCHEMA, "release_id": rid, "verdict": "PASS", "immutable": True,
            "source_task_id": TASK_ID, "created_at_utc": utcnow(), "external_release_path": str(final),
            "collection_release_id": COLLECTION_RELEASE_ID, "collection_release_json_sha256": COLLECTION_RELEASE_SHA256,
            "input_stage_b_manifest_sha256": STAGE_B_SHA256, "pansn_policy_version": POLICY,
            "exact_version_policy": "requested exact assembly revision is identity; never substitute latest",
            "sequence_bearing_assembly_accessions": EXPECTED_ACCESSIONS,
            "counts": {"attempted_exact_assembly_revisions": 10, "validated": 10, "quarantined": 0,
                       "distinct_sequence_bearing_assemblies": 10, "global_distinct_assembly_cap": 1000,
                       "contigs": len(contig_rows), "total_bases": sum(int(row["contig_length"]) for row in contig_rows),
                       "annotations_alias_validated": sum(row["annotation_status"].startswith("ALIASES_VALIDATED") for row in assembly_rows)},
            "manifests": {
                "stage_b_10.tsv": {"rows": 10, "bytes": (stage / "manifests/stage_b_10.tsv").stat().st_size, "sha256": sha_file(stage / "manifests/stage_b_10.tsv")},
                "assemblies.tsv": {"rows": 10, "bytes": (stage / "manifests/assemblies.tsv").stat().st_size, "sha256": sha_file(stage / "manifests/assemblies.tsv")},
                "contigs.tsv": {"rows": len(contig_rows), "bytes": (stage / "manifests/contigs.tsv").stat().st_size, "sha256": sha_file(stage / "manifests/contigs.tsv")},
                "checksums.tsv": {"rows": len(checksum_rows), "bytes": (stage / "manifests/checksums.tsv").stat().st_size, "sha256": sha_file(stage / "manifests/checksums.tsv")},
                "state.tsv": {"rows": 10, "bytes": (stage / "manifests/state.tsv").stat().st_size, "sha256": sha_file(stage / "manifests/state.tsv")},
            },
            "root_inputs_start": root_start, "root_inputs_finish": root_finish,
            "resource_summary_sha256": sha_file(stage / "resource_summary.json"),
            "restart_evidence_sha256": sha_file(stage / "restart_evidence.json"),
            "global_cap_evidence_sha256": sha_file(stage / "global_cap_evidence.json"),
            "applicable_gates": gates,
            "annotation_policy": "source GFF retained in checksum-validated package; only validated explicit alias tables published; no transformed GFF",
            "routine_whole_set_plain_fasta_files": 0,
        }
        (stage / "release.json").write_bytes(canonical_json(release))
        append_jsonl(state, {"event": "READY_TO_PROMOTE", "release_id": rid, "at": utcnow()})
        # The final state append changes the release inventory, so release.json is already fixed but inventories are made afterward.
        seal_directory(stage, final)
        global_finish = audit_global_release_cap(durable_task.parents[1], allowed_collection)
        if global_finish["distinct_exact_assembly_revisions"] > 10 or set(global_finish["accessions"]) != set(EXPECTED_ACCESSIONS):
            raise GateError("finish global cap inventory is not exactly the bounded Stage-B set")
        publish_tracked(final, tracked, artifact_dir)
        return release
    except BaseException as exc:
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        append_jsonl(failures, {"event": "RUN_FAILED", "type": type(exc).__name__, "message": str(exc), "at": utcnow()})
        raise


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    run_p = sub.add_parser("run")
    run_p.add_argument("--repo-root", default=".")
    run_p.add_argument("--tracked-root", default="manifests/canonical-cohort-010-v1")
    run_p.add_argument("--artifact-root", default="artifacts/acquisition_canonicalization_10")
    run_p.add_argument("--durable-task-root", required=True)
    run_p.add_argument("--scratch-root", required=True)
    run_p.add_argument("--run-id", required=True)
    run_p.add_argument("--retries", type=int, default=8)
    run_p.add_argument("--rate-delay", type=float, default=0.5)
    run_p.add_argument("--assigned-ram-bytes", type=int, required=True)
    run_p.add_argument("--durable-allocation-bytes", type=int, required=True)
    run_p.add_argument("--scratch-allocation-bytes", type=int, required=True)
    run_p.add_argument("--inode-allocation", type=int, required=True)
    run_p.add_argument("--predicted-durable-peak-bytes", type=int, required=True)
    run_p.add_argument("--predicted-scratch-peak-bytes", type=int, required=True)
    run_p.add_argument("--predicted-files", type=int, required=True)
    run_p.add_argument("--unfinished-write-bytes", type=int, required=True)
    run_p.add_argument("--bgzip", default="bgzip")
    run_p.add_argument("--samtools", default="samtools")
    run_p.add_argument("--bgzip-threads", type=int, default=2)
    run_p.add_argument("--bgzip-level", type=int, default=6)
    run_p.add_argument("--inject-kill", choices=("none", "acquisition", "conversion"), default="none")
    run_p.add_argument("--inject-accession", choices=EXPECTED_ACCESSIONS, default=EXPECTED_ACCESSIONS[0])
    run_p.add_argument("--inject-after-bytes", type=int, default=131072)
    run_p.add_argument("--inject-after-bases", type=int, default=262144)
    return p


def main() -> int:
    try:
        args = parser().parse_args()
        if args.command == "run":
            print(json.dumps(run(args), sort_keys=True))
        return 0
    except GateError as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
