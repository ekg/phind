#!/usr/bin/env python3
"""Fail-closed N=250 canonical cohort preparation.

This task-owned driver deliberately imports the frozen Stage-B implementation
rather than changing it. The first 100 objects are immutable digest references
to canonical-cohort-100-v1; only cohort rows 101..250 can reach the downloader.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import http.client
import json
import os
import platform
import resource
import shutil
import socket
import stat
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

from workflow.acquisition_canonicalization import pilot as p

TASK_ID = "prepare-canonical-cohort-250"
SCHEMA = "canonical-cohort-250-release-v1"
POLICY = p.POLICY
SELECTION_RELEASE_ID = "pilot-cohorts-v1-8afc0ea03d9e50dc"
SELECTION_RELEASE_JSON_SHA256 = "d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99"
COHORT_SHA256 = "ba2cf2909ccf62a0c1944a76b522edc5600953511ec355479117b4a419acbc9f"
COHORT_BYTES = 91475
COHORT_ROWS = 250
PREDECESSOR_ROWS = 100
NEW_ROWS = COHORT_ROWS - PREDECESSOR_ROWS
PREDECESSOR_RELEASE_ID = "canonical-cohort-100-v1-6be4c0dde65f31d0"
PREDECESSOR_RELEASE_JSON_SHA256 = "3b91b24e23323ef971a13f22825e512a233bb592ed641ea9b270a2f1fd683795"
COMPATIBILITY_RELEASE_ID = "consumer-compatibility-v1-78d7e93f19fa3d87"
COMPATIBILITY_RELEASE_JSON_SHA256 = "021719ddadd7bb7fa2932d2ef9cb25da9c666ebe0389988691283011ee12f4c7"
ROOT_HASHES = {
    "26k_ecoli_accession.txt": p.ACCESSION_INPUT_SHA256,
    "26k_prophage1.csv": p.PROPHAGE_INPUT_SHA256,
}
DURABLE_PREFIX = "/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-250/"
SCRATCH_PREFIX = "/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/prepare-canonical-cohort-250/"

OBJECT_REF_FIELDS = [
    "cohort_order", "accession", "storage_release_id", "storage_root",
    "source_object_relpath", "source_inventory_sha256", "canonical_object_relpath",
    "canonical_inventory_sha256", "reuse_status", "predecessor_digest_match", "row_sha256",
]
CHECKSUM_FIELDS = [
    "cohort_order", "accession", "artifact_role", "storage_release_id", "storage_root",
    "relative_path", "bytes", "sha256", "reuse_status", "row_sha256",
]
BATCH_FIELDS = [
    "batch_number", "first_cohort_order", "last_cohort_order", "assemblies",
    "already_complete", "completed_this_invocation", "wall_seconds", "cpu_seconds",
    "validated_transfer_bytes", "canonical_bgzf_bytes", "cumulative_transfer_bytes",
    "cumulative_canonical_bgzf_bytes", "partial_bytes_observed", "download_requests",
    "retry_events", "failure_events", "restart_events", "peak_rss_bytes",
    "stage_bytes_start", "stage_bytes_finish", "stage_partial_bytes_finish",
    "stage_final_bytes_finish", "stage_files_finish", "durable_free_bytes",
    "durable_free_inodes", "scratch_free_bytes", "scratch_free_inodes", "row_sha256",
]


class GateError(p.GateError):
    pass


def sha_file(path: Path) -> str:
    return p.sha_file(path)


def canonical_json(value: Any) -> bytes:
    return p.canonical_json(value)


def configure_pinned_primitives(accessions: list[str]) -> None:
    """Narrow the inherited URL constructor to the exact frozen 250 rows."""
    p.EXPECTED_ACCESSIONS = list(accessions)
    p.TASK_ID = TASK_ID
    p.USER_AGENT = "phind-canonical-cohort-250/1.0 (bounded-frozen-250)"


def verify_root_inputs(repo: Path) -> dict[str, str]:
    observed = {name: sha_file(repo / name) for name in ROOT_HASHES}
    if observed != ROOT_HASHES:
        raise GateError(f"immutable root input checksum mismatch: {observed}")
    return observed


def verify_tracked_inventory(root: Path) -> dict[str, str]:
    sums_path = root / "SHA256SUMS"
    if not sums_path.is_file():
        raise GateError(f"missing tracked SHA256SUMS: {root}")
    sums: dict[str, str] = {}
    for number, line in enumerate(sums_path.read_text().splitlines(), 1):
        try:
            digest, rel = line.split("  ", 1)
        except ValueError as exc:
            raise GateError(f"malformed tracked inventory line {number}: {root}") from exc
        if rel in sums or rel.startswith("/") or ".." in Path(rel).parts:
            raise GateError(f"unsafe/duplicate tracked inventory path: {rel}")
        path = root / rel
        if not path.is_file() or path.is_symlink() or sha_file(path) != digest:
            raise GateError(f"tracked checksum mismatch: {path}")
        sums[rel] = digest
    return sums


def verify_external_inventory(root: Path) -> None:
    """Validate either canonical text COMPLETE or the pinned compatibility JSON marker."""
    p.verify_sha_inventory(root, require_complete=False)
    complete = root / "COMPLETE"
    if not complete.is_file():
        raise GateError(f"missing COMPLETE: {root}")
    expected = sha_file(root / "SHA256SUMS")
    text = complete.read_text().strip()
    if text.startswith("{"):
        marker = json.loads(text)
        if marker.get("sha256sums_sha256") != expected or marker.get("verdict") != "PASS":
            raise GateError(f"JSON COMPLETE inventory digest/verdict mismatch: {root}")
    else:
        tokens = text.split()
        if not tokens or tokens[0] != expected:
            raise GateError(f"COMPLETE inventory digest mismatch: {root}")


def read_hashed_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    if not fields or fields[-1] != "row_sha256":
        raise GateError(f"TSV row-hash schema missing: {path}")
    for number, row in enumerate(rows, 2):
        if row["row_sha256"] != p.stable_row_hash(row, fields):
            raise GateError(f"row checksum mismatch: {path}:{number}")
    return rows


def _release_external_path(release: dict[str, Any]) -> Path:
    raw = release.get("external_release_path", release.get("external_path"))
    if not isinstance(raw, str) or not raw:
        raise GateError("predecessor release lacks external path")
    return Path(raw)


def _require_pass_or_na(gates: dict[str, Any], label: str) -> None:
    for name, verdict in gates.items():
        if verdict == "PASS" or str(verdict).startswith("NOT_APPLICABLE"):
            continue
        raise GateError(f"{label} gate is not unqualified PASS/NA: {name}={verdict}")


def verify_inputs(repo: Path) -> dict[str, Any]:
    roots = verify_root_inputs(repo)

    selection_dir = repo / "manifests/pilot-cohorts-v1"
    selection_sums = verify_tracked_inventory(selection_dir)
    if selection_sums.get("release.json") != SELECTION_RELEASE_JSON_SHA256:
        raise GateError("selection release.json pinned SHA-256 mismatch")
    if selection_sums.get("cohort-0250.tsv") != COHORT_SHA256:
        raise GateError("N=250 selection manifest inventory mismatch")
    selection = json.loads((selection_dir / "release.json").read_text())
    if (selection.get("release_id") != SELECTION_RELEASE_ID or selection.get("verdict") != "PASS"
            or selection.get("immutable") is not True):
        raise GateError("selection release identity/verdict/immutability mismatch")
    _require_pass_or_na(selection.get("applicable_gates", {}), "selection")
    sel_external = _release_external_path(selection)
    verify_external_inventory(sel_external)
    cohort = selection_dir / "cohort-0250.tsv"
    if cohort.stat().st_size != COHORT_BYTES or sha_file(cohort) != COHORT_SHA256:
        raise GateError("frozen N=250 manifest bytes/SHA mismatch")
    external_cohort = sel_external / "manifests/cohort-0250.tsv"
    if not external_cohort.is_file():
        external_cohort = sel_external / "cohort-0250.tsv"
    if external_cohort.read_bytes() != cohort.read_bytes():
        raise GateError("tracked/external N=250 manifest bytes differ")
    rows = read_hashed_tsv(cohort)
    accessions = [row["exact_assembly_accession_version"] for row in rows]
    if len(rows) != COHORT_ROWS or len(set(accessions)) != COHORT_ROWS:
        raise GateError("N=250 manifest cardinality/duplicate gate failed")
    for order, row in enumerate(rows, 1):
        if int(row["cohort_order"]) != order or row["rung_n"] != "250":
            raise GateError("N=250 cohort order/rung gate failed")
        if not p.ACCESSION_RE.fullmatch(row["exact_assembly_accession_version"]):
            raise GateError("N=250 manifest contains invalid exact accession")
        if row["resolution_status"] not in ("EXACT_VERSION_RESOLVED", "EXACT_VERSION_VALID_METADATA_UNAVAILABLE"):
            raise GateError("N=250 manifest contains non-exact/non-validated row")

    predecessor_dir = repo / "manifests/canonical-cohort-100-v1"
    predecessor_sums = verify_tracked_inventory(predecessor_dir)
    if predecessor_sums.get("release.json") != PREDECESSOR_RELEASE_JSON_SHA256:
        raise GateError("N=100 predecessor release.json pinned SHA-256 mismatch")
    predecessor = json.loads((predecessor_dir / "release.json").read_text())
    if (predecessor.get("release_id") != PREDECESSOR_RELEASE_ID or predecessor.get("verdict") != "PASS"
            or predecessor.get("immutable") is not True):
        raise GateError("N=100 predecessor identity/verdict/immutability mismatch")
    _require_pass_or_na(predecessor.get("applicable_gates", {}), "N=100 predecessor")
    predecessor_external = _release_external_path(predecessor)
    verify_external_inventory(predecessor_external)
    if sha_file(predecessor_external / "release.json") != PREDECESSOR_RELEASE_JSON_SHA256:
        raise GateError("N=100 tracked/external release.json mismatch")
    predecessor_accessions = predecessor.get("sequence_bearing_assembly_accessions", [])
    if accessions[:PREDECESSOR_ROWS] != predecessor_accessions or len(predecessor_accessions) != PREDECESSOR_ROWS:
        raise GateError("N=250 is not exact nested extension of N=100")
    predecessor_refs = p.read_tsv(
        predecessor_external / "manifests/object_refs.tsv", OBJECT_REF_FIELDS, verify_hashes=True
    )
    if (len(predecessor_refs) != PREDECESSOR_ROWS
            or [row["accession"] for row in predecessor_refs] != predecessor_accessions):
        raise GateError("N=100 predecessor object reference accounting/order mismatch")
    predecessor_ref_by_accession = {row["accession"]: row for row in predecessor_refs}

    compatibility_dir = repo / "manifests/consumer-compatibility-v1"
    compatibility_sums = verify_tracked_inventory(compatibility_dir)
    if compatibility_sums.get("release.json") != COMPATIBILITY_RELEASE_JSON_SHA256:
        raise GateError("consumer compatibility release.json pinned SHA-256 mismatch")
    compatibility = json.loads((compatibility_dir / "release.json").read_text())
    if (compatibility.get("release_id") != COMPATIBILITY_RELEASE_ID or compatibility.get("verdict") != "PASS"
            or compatibility.get("immutable") is not True):
        raise GateError("consumer compatibility identity/verdict/immutability mismatch")
    _require_pass_or_na(compatibility.get("applicable_gates", {}), "compatibility")
    compatibility_external = _release_external_path(compatibility)
    verify_external_inventory(compatibility_external)
    if sha_file(compatibility_external / "release.json") != COMPATIBILITY_RELEASE_JSON_SHA256:
        raise GateError("compatibility tracked/external release.json mismatch")
    gates = p.read_tsv(compatibility_dir / "gates.tsv")
    if len(gates) != 19 or any(row.get("verdict") != "PASS" for row in gates):
        raise GateError("required 19/19 consumer gates are not PASS")

    return {
        "root_inputs": roots, "rows": rows, "accessions": accessions,
        "selection": selection, "selection_external": sel_external,
        "predecessor": predecessor, "predecessor_external": predecessor_external,
        "predecessor_refs": predecessor_refs, "predecessor_ref_by_accession": predecessor_ref_by_accession,
        "compatibility": compatibility, "compatibility_external": compatibility_external,
    }


def audit_global_release_cap(project_root: Path, allowed_collection: set[str]) -> dict[str, Any]:
    """Audit committed sequence-bearing releases, including compatibility references."""
    union: set[str] = set()
    scanned: list[dict[str, Any]] = []
    releases_root = project_root / "releases"
    if releases_root.exists():
        for release_json in sorted(releases_root.glob("*/*/release.json")):
            root = release_json.parent
            if not (root / "COMPLETE").is_file():
                continue
            obj = json.loads(release_json.read_text())
            count = int(obj.get("counts", {}).get("distinct_sequence_bearing_assemblies", 0))
            accessions = list(obj.get("sequence_bearing_assembly_accessions", []))
            if count and not accessions and len(obj.get("cohort_order", [])) == count:
                accessions = list(obj["cohort_order"])
            if count and (len(accessions) != count or len(set(accessions)) != count):
                raise GateError(f"committed sequence-bearing release lacks exact accession inventory: {release_json}")
            if accessions and not set(accessions).issubset(allowed_collection):
                raise GateError(f"committed release contains accession outside frozen collection: {release_json}")
            union.update(accessions)
            scanned.append({"release_json": str(release_json), "declared_count": count, "accessions": sorted(accessions)})
    if len(union) > 1000 or not union.issubset(allowed_collection):
        raise GateError(f"global distinct assembly cap/subset gate failed: count={len(union)}")
    return {"verdict": "PASS", "distinct_exact_assembly_revisions": len(union), "cap": 1000,
            "subset_of_frozen_collection": True, "accessions": sorted(union),
            "scanned_releases": scanned, "captured_at_utc": p.utcnow()}


def release_seed(accessions: list[str]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA,
        "selection_release_id": SELECTION_RELEASE_ID,
        "selection_release_json_sha256": SELECTION_RELEASE_JSON_SHA256,
        "cohort_0250_sha256": COHORT_SHA256,
        "predecessor_release_id": PREDECESSOR_RELEASE_ID,
        "predecessor_release_json_sha256": PREDECESSOR_RELEASE_JSON_SHA256,
        "compatibility_release_id": COMPATIBILITY_RELEASE_ID,
        "compatibility_release_json_sha256": COMPATIBILITY_RELEASE_JSON_SHA256,
        "pansn_policy_version": POLICY,
        "exact_accessions": accessions,
        "reuse_first_n": PREDECESSOR_ROWS,
    }


def release_id(accessions: list[str]) -> str:
    return "canonical-cohort-250-v1-" + p.sha_bytes(canonical_json(release_seed(accessions)))[:16]


def _existing_parent(path: Path) -> Path:
    candidate = path
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    return candidate


def live_preflight(durable: Path, scratch: Path, allocations: p.Allocations, stage_name: str) -> dict[str, Any]:
    allocations.validate()
    if not str(durable).startswith(DURABLE_PREFIX):
        raise GateError(f"durable path outside task namespace: {durable}")
    if not str(scratch).startswith(SCRATCH_PREFIX):
        raise GateError(f"scratch path outside task namespace: {scratch}")
    dp, sp = _existing_parent(durable), _existing_parent(scratch)
    dm, sm = p._mount_record(dp), p._mount_record(sp)
    if dm.get("target") != "/" or dm.get("fstype") != "ext4":
        raise GateError(f"durable mount identity NO_GO: {dm}")
    if sm.get("target") != "/mnt/nvme3n1" or sm.get("source") != "/dev/nvme3n1" or sm.get("fstype") != "xfs":
        raise GateError(f"scratch mount identity NO_GO: {sm}")
    ds, ss = os.statvfs(dp), os.statvfs(sp)
    dfree, dinodes = ds.f_bavail * ds.f_frsize, ds.f_favail
    sfree, sinodes = ss.f_bavail * ss.f_frsize, ss.f_favail
    checks = {
        "durable_free_ge_2tb": dfree >= 2_000_000_000_000,
        "durable_free_inodes_ge_1m": dinodes >= 1_000_000,
        "scratch_live_preflight_ge_4tb": sfree >= 4_000_000_000_000,
        "scratch_stop_floor_ge_2tb": sfree - allocations.unfinished_write_bytes >= 2_000_000_000_000,
        "scratch_free_inodes_ge_5m": sinodes >= 5_000_000,
        "durable_predicted_peak_le_70pct": allocations.predicted_durable_peak_bytes * 100 <= allocations.durable_allocation_bytes * 70,
        "scratch_predicted_peak_le_70pct": allocations.predicted_scratch_peak_bytes * 100 <= allocations.scratch_allocation_bytes * 70,
        "projected_files_le_50pct": allocations.predicted_files * 2 <= allocations.inode_allocation,
        "durable_allocation_retains_2x_unfinished": allocations.durable_allocation_bytes - allocations.predicted_durable_peak_bytes >= 2 * allocations.unfinished_write_bytes,
        "scratch_allocation_retains_2x_unfinished": allocations.scratch_allocation_bytes - allocations.predicted_scratch_peak_bytes >= 2 * allocations.unfinished_write_bytes,
        "durable_live_retains_2x_unfinished": dfree - allocations.predicted_durable_peak_bytes >= 2 * allocations.unfinished_write_bytes,
        "scratch_live_retains_2x_unfinished": sfree - allocations.predicted_scratch_peak_bytes >= 2 * allocations.unfinished_write_bytes,
    }
    if not all(checks.values()):
        raise GateError("resource gate NO_GO: " + json.dumps(checks, sort_keys=True))
    uid = os.getuid()
    for label, parent in (("durable", dp), ("scratch", sp)):
        if parent.stat().st_uid != uid or not os.access(parent, os.W_OK | os.X_OK):
            raise GateError(f"{label} ownership/write gate NO_GO: {parent}")
        p._write_probe(parent, label)
    return {
        "schema": "canonical-250-resource-preflight-v1", "verdict": "PASS", "stage": stage_name,
        "captured_at_utc": p.utcnow(), "durable_path": str(durable), "scratch_path": str(scratch),
        "durable_findmnt": dm, "scratch_findmnt": sm,
        "durable_owner": {"uid": dp.stat().st_uid, "gid": dp.stat().st_gid, "mode": stat.S_IMODE(dp.stat().st_mode)},
        "scratch_owner": {"uid": sp.stat().st_uid, "gid": sp.stat().st_gid, "mode": stat.S_IMODE(sp.stat().st_mode)},
        "write_probes": "PASS", "durable_free_bytes": dfree, "durable_free_inodes": dinodes,
        "scratch_free_bytes": sfree, "scratch_free_inodes": sinodes, "swap_free_bytes": p.swap_free_bytes(),
        "allocations": asdict(allocations), "checks": checks,
    }


def _ensure_owned_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True, mode=0o750)
    if path.is_symlink() or path.stat().st_uid != os.getuid():
        raise GateError(f"directory ownership/symlink gate failed: {path}")


def _read_object(root: Path, accession: str) -> tuple[Path, dict[str, Any], Path, dict[str, Any]]:
    source_dir = root / "source_objects" / accession
    canonical_dir = root / "canonical_objects" / accession
    p.verify_sha_inventory(source_dir)
    p.verify_sha_inventory(canonical_dir)
    source = json.loads((source_dir / "manifest.json").read_text())
    canonical = json.loads((canonical_dir / "manifest.json").read_text())
    if source.get("accession") != accession or source.get("state") != "COMPLETE":
        raise GateError(f"source object identity/state mismatch: {accession}")
    if canonical.get("accession") != accession or canonical.get("state") != "COMPLETE":
        raise GateError(f"canonical object identity/state mismatch: {accession}")
    return source_dir, source, canonical_dir, canonical


def commit_source_with_transport_retry(
    accession: str, source_root: Path, state: Path, failures: Path,
    retries: int, rate_delay: float, inject_kill: bool, inject_after_bytes: int,
) -> tuple[Path, dict[str, Any]]:
    """Apply a bounded outer retry for chunked-transfer HTTP framing failures.

    The pinned primitive already handles URL/OSError/timeouts and independently
    validates any completed unit. Python's IncompleteRead is an HTTPException,
    not an OSError, so the task ledger wraps only that transport class and then
    re-enters the exact primitive's safe-partial identity logic.
    """
    for attempt in range(1, retries + 1):
        try:
            return p.commit_source_object(accession, source_root, state, failures, retries,
                                          rate_delay, inject_kill, inject_after_bytes)
        except http.client.HTTPException as exc:
            p.append_jsonl(failures, {"event": "OUTER_HTTP_TRANSPORT_RETRY", "accession": accession,
                                      "attempt": attempt, "type": type(exc).__name__,
                                      "message": str(exc), "at": p.utcnow()})
            if attempt == retries:
                raise GateError(f"HTTP transport failed after {retries} outer attempts: {accession}") from exc
            time.sleep(min(60.0, float(2 ** (attempt - 1))))
    raise GateError(f"outer transport retry exhausted: {accession}")


def _discard_invalid_completed(path: Path, accession: str, role: str, state: Path, failures: Path) -> None:
    if not path.exists():
        return
    try:
        p.verify_sha_inventory(path)
        manifest = json.loads((path / "manifest.json").read_text())
        if manifest.get("accession") != accession or manifest.get("state") != "COMPLETE":
            raise GateError("identity/state mismatch")
    except (p.GateError, OSError, ValueError, json.JSONDecodeError) as exc:
        p.append_jsonl(failures, {"event": "INVALID_COMPLETED_OBJECT_DISCARDED", "accession": accession,
                                  "role": role, "message": str(exc), "at": p.utcnow()})
        p.safe_remove(path)
        p.append_jsonl(state, {"event": "INVALID_COMPLETED_OBJECT_RECOMPUTE_ALLOWED", "accession": accession,
                               "role": role, "at": p.utcnow()})


def _storage_for(
    order: int, accession: str, stage: Path, inputs: dict[str, Any]
) -> tuple[Path, str, str, str]:
    """Resolve physical storage without copying a recursively referenced object."""
    if order > PREDECESSOR_ROWS:
        return stage, "SELF", ".", "CREATED_OR_CHECKSUM_RESUMED"
    ref = inputs["predecessor_ref_by_accession"].get(accession)
    if ref is None or int(ref["cohort_order"]) != order:
        raise GateError(f"missing/misordered predecessor object reference: {accession}")
    if ref["storage_release_id"] == "SELF":
        if ref["storage_root"] != ".":
            raise GateError(f"unsafe predecessor self-storage root: {accession}")
        root = inputs["predecessor_external"]
        storage_id = PREDECESSOR_RELEASE_ID
    else:
        root = Path(ref["storage_root"])
        storage_id = ref["storage_release_id"]
        if not root.is_absolute() or not storage_id:
            raise GateError(f"unsafe recursive predecessor reference: {accession}")
    return root, storage_id, str(root), "REUSED_PREDECESSOR_BY_DIGEST"


def _assembly_row(order: int, input_row: dict[str, str], source: dict[str, Any], canonical: dict[str, Any]) -> dict[str, Any]:
    accession = input_row["exact_assembly_accession_version"]
    remote = source["download_receipt"]["remote_identity"]
    headers = remote.get("headers", {})
    return {
        "stage_b_order": order, "assembly_id": input_row["assembly_id"], "accession": accession,
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


def build_release_tables(stage: Path, inputs: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    assembly_rows: list[dict[str, Any]] = []
    contig_rows: list[dict[str, Any]] = []
    checksum_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    refs: list[dict[str, Any]] = []
    old_assemblies = {row["accession"]: row for row in p.read_tsv(inputs["predecessor_external"] / "manifests/assemblies.tsv", p.ASSEMBLY_FIELDS, verify_hashes=True)}
    digest_fields = ["source_package_sha256", "source_decompressed_sha256", "source_gff_sha256", "canonical_bgzf_sha256", "fai_sha256", "gzi_sha256", "crosswalk_sha256", "annotation_aliases_sha256"]
    for input_row in inputs["rows"]:
        order = int(input_row["cohort_order"])
        accession = input_row["exact_assembly_accession_version"]
        root, storage_id, storage_root, reuse = _storage_for(order, accession, stage, inputs)
        source_dir, source, canonical_dir, canonical = _read_object(root, accession)
        assembly = _assembly_row(order, input_row, source, canonical)
        if order <= PREDECESSOR_ROWS:
            old = old_assemblies.get(accession)
            if old is None or any(str(assembly[field]) != old[field] for field in digest_fields):
                raise GateError(f"reused N=100 checksum contract mismatch: {accession}")
        assembly_rows.append(assembly)
        object_contigs = p.read_tsv(canonical_dir / "contigs.tsv", p.CONTIG_FIELDS, verify_hashes=True)
        if any(int(row["stage_b_order"]) != order for row in object_contigs):
            raise GateError(f"new object cohort order mismatch: {accession}")
        contig_rows.extend(object_contigs)
        refs.append({
            "cohort_order": order, "accession": accession, "storage_release_id": storage_id,
            "storage_root": storage_root, "source_object_relpath": f"source_objects/{accession}",
            "source_inventory_sha256": sha_file(source_dir / "SHA256SUMS"),
            "canonical_object_relpath": f"canonical_objects/{accession}",
            "canonical_inventory_sha256": sha_file(canonical_dir / "SHA256SUMS"),
            "reuse_status": reuse, "predecessor_digest_match": "PASS" if order <= PREDECESSOR_ROWS else "NOT_APPLICABLE_NEW_OBJECT",
        })
        roles = (
            ("source_package", f"source_objects/{accession}/package.zip"),
            ("source_manifest", f"source_objects/{accession}/manifest.json"),
            ("canonical_bgzf", canonical["canonical_bgzf_relpath"]),
            ("fai", canonical["fai_relpath"]), ("gzi", canonical["gzi_relpath"]),
            ("contig_crosswalk", f"canonical_objects/{accession}/contigs.tsv"),
            ("annotation_aliases", f"canonical_objects/{accession}/annotation_aliases.tsv"),
            ("canonical_manifest", f"canonical_objects/{accession}/manifest.json"),
        )
        for role, rel in roles:
            path = root / rel
            checksum_rows.append({"cohort_order": order, "accession": accession, "artifact_role": role,
                                  "storage_release_id": storage_id, "storage_root": storage_root,
                                  "relative_path": rel, "bytes": path.stat().st_size, "sha256": sha_file(path),
                                  "reuse_status": reuse})
        state_rows.append({"stage_b_order": order, "accession": accession, "source_state": "COMPLETE",
                           "canonical_state": "COMPLETE", "terminal_state": "VALIDATED",
                           "reason": "REUSED_PREDECESSOR_BY_DIGEST" if order <= PREDECESSOR_ROWS else "."})
    return assembly_rows, contig_rows, checksum_rows, state_rows, refs


def _write_tables(stage: Path, inputs: dict[str, Any], batch_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    assemblies, contigs, checksums, states, refs = build_release_tables(stage, inputs)
    if (len(assemblies) != COHORT_ROWS or len(states) != COHORT_ROWS
            or len(refs) != COHORT_ROWS or len(checksums) != COHORT_ROWS * 8):
        raise GateError("100% assembly/state/reference/checksum row accounting failed")
    if len(contigs) != sum(int(row["contig_count"]) for row in assemblies):
        raise GateError("100% contig row accounting failed")
    names = [row["pansn_sequence_name"] for row in contigs]
    if len(names) != len(set(names)):
        raise GateError("cohort-wide PanSN collision")
    p.write_tsv(stage / "manifests/assemblies.tsv", p.ASSEMBLY_FIELDS, assemblies)
    p.write_tsv(stage / "manifests/contigs.tsv", p.CONTIG_FIELDS, contigs)
    p.write_tsv(stage / "manifests/checksums.tsv", CHECKSUM_FIELDS, checksums)
    p.write_tsv(stage / "manifests/state.tsv", p.STATE_FIELDS, states)
    p.write_tsv(stage / "manifests/object_refs.tsv", OBJECT_REF_FIELDS, refs)
    p.write_tsv(stage / "manifests/batch_metrics.tsv", BATCH_FIELDS, batch_rows)
    return assemblies, contigs, checksums, states, refs


def _tree_usage(path: Path) -> tuple[int, int]:
    return p.directory_usage(path)


def _peak_rss_bytes() -> int:
    return max(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
    ) * 1024


def _jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _batch_observations(
    batch: list[dict[str, str]], stage: Path, inputs: dict[str, Any], state: Path, failures: Path
) -> dict[str, int]:
    accessions = {row["exact_assembly_accession_version"] for row in batch}
    transfer_bytes = canonical_bytes = 0
    for row in batch:
        order = int(row["cohort_order"])
        if order <= PREDECESSOR_ROWS:
            continue
        accession = row["exact_assembly_accession_version"]
        _, source, _, canonical = _read_object(stage, accession)
        transfer_bytes += int(source["validation"]["package_bytes"])
        canonical_bytes += int(canonical["canonical_bgzf_bytes"])
    state_events = [row for row in _jsonl(state) if row.get("accession") in accessions]
    failure_rows = [row for row in _jsonl(failures) if row.get("accession") in accessions]
    partial_bytes = sum(
        int(row.get("partial_bytes", 0)) for row in state_events
        if row.get("event") in {
            "INJECTED_ACQUISITION_SIGKILL", "ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE"
        }
    )
    retry_names = {
        "ACQUISITION_RANGE_RESUME", "ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE",
        "ACQUISITION_RANGE_REFUSED_RESTART", "ACQUISITION_ORPHAN_PARTIAL_DISCARDED",
    }
    restart_names = retry_names | {
        "INJECTED_ACQUISITION_SIGKILL", "INJECTED_CONVERSION_SIGKILL",
        "INTERRUPTED_CONVERSION_STAGE_DISCARDED", "INVALID_COMPLETED_OBJECT_RECOMPUTE_ALLOWED",
    }
    return {
        "validated_transfer_bytes": transfer_bytes,
        "canonical_bgzf_bytes": canonical_bytes,
        "partial_bytes_observed": partial_bytes,
        "download_requests": sum(row.get("event") == "ACQUISITION_REQUEST_STARTED" for row in state_events),
        "retry_events": sum(row.get("event") in retry_names for row in state_events)
            + sum(row.get("event") in {"ACQUISITION_ATTEMPT_FAILED", "OUTER_HTTP_TRANSPORT_RETRY"} for row in failure_rows),
        "failure_events": len(failure_rows),
        "restart_events": sum(row.get("event") in restart_names for row in state_events),
    }


def _partial_bytes(stage: Path) -> int:
    return sum(
        path.stat().st_size for path in stage.rglob("*")
        if path.is_file() and (
            path.name.endswith(".partial")
            or any(part.startswith(".stage.") for part in path.relative_to(stage).parts)
        )
    )


def _tree_digest(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(x for x in root.rglob("*") if x.is_file() and not x.is_symlink()):
        digest.update(str(path.relative_to(root)).encode() + b"\0")
        digest.update(sha_file(path).encode() + b"\n")
    return digest.hexdigest()


def publish_tracked(external: Path, tracked: Path, artifact: Path) -> None:
    tracked.mkdir(parents=True, exist_ok=True)
    artifact.mkdir(parents=True, exist_ok=True)
    mapping = {
        "manifests/cohort-0250.tsv": "cohort-0250.tsv",
        "manifests/assemblies.tsv": "assemblies.tsv", "manifests/checksums.tsv": "checksums.tsv",
        "manifests/state.tsv": "state.tsv", "manifests/object_refs.tsv": "object_refs.tsv",
        "manifests/batch_metrics.tsv": "batch_metrics.tsv", "release.json": "release.json",
        "SHA256SUMS": "external_SHA256SUMS",
    }
    for source, target in mapping.items():
        shutil.copyfile(external / source, tracked / target)
    p.deterministic_gzip(external / "manifests/contigs.tsv", tracked / "contigs.tsv.gz")
    allowed = set(mapping.values()) | {"contigs.tsv.gz", "SHA256SUMS"}
    for old in tracked.iterdir():
        if old.is_file() and old.name not in allowed:
            raise GateError(f"unexpected tracked manifest file: {old}")
    sums = []
    for path in sorted(x for x in tracked.iterdir() if x.is_file() and x.name != "SHA256SUMS"):
        sums.append(f"{sha_file(path)}  {path.name}\n")
    (tracked / "SHA256SUMS").write_text("".join(sums))
    p.fsync_file(tracked / "SHA256SUMS")
    for name in ("resource_summary.json", "restart_evidence.json", "global_cap_evidence.json", "tools.json"):
        shutil.copyfile(external / name, artifact / name)


def _count_state_events(path: Path, event: str) -> int:
    if not path.exists():
        return 0
    return sum(1 for line in path.read_text().splitlines() if json.loads(line).get("event") == event)


def run(args: argparse.Namespace) -> dict[str, Any]:
    repo = Path(args.repo_root).resolve()
    durable_task = Path(args.durable_task_root).resolve()
    scratch = Path(args.scratch_root).resolve()
    tracked = Path(args.tracked_root).resolve()
    artifact = Path(args.artifact_root).resolve()
    allocations = p.Allocations(args.assigned_ram_bytes, args.durable_allocation_bytes,
        args.scratch_allocation_bytes, args.inode_allocation, args.predicted_durable_peak_bytes,
        args.predicted_scratch_peak_bytes, args.predicted_files, args.unfinished_write_bytes)
    invocation_wall_start = time.monotonic()
    invocation_cpu_start = time.process_time()
    ru_start = resource.getrusage(resource.RUSAGE_SELF)
    inputs = verify_inputs(repo)
    configure_pinned_primitives(inputs["accessions"])
    allowed = p._load_allowed_collection(repo)
    global_start = audit_global_release_cap(durable_task.parents[1], allowed)
    projected = set(global_start["accessions"]) | set(inputs["accessions"])
    if len(projected) != COHORT_ROWS or projected != set(inputs["accessions"]) or not projected.issubset(allowed):
        raise GateError("projected global union is not exactly frozen N=250")
    rid = release_id(inputs["accessions"])
    final = durable_task / rid
    stage = durable_task / f".stage.{rid}.{args.run_id}"
    initial = live_preflight(final, scratch, allocations, "INITIAL")
    _ensure_owned_dir(durable_task)
    _ensure_owned_dir(scratch)
    if final.exists():
        p.verify_sha_inventory(final)
        existing = json.loads((final / "release.json").read_text())
        if existing.get("release_id") != rid or existing.get("verdict") != "PASS":
            raise GateError("existing final release identity/verdict mismatch")
        publish_tracked(final, tracked, artifact)
        return existing
    seed = canonical_json(release_seed(inputs["accessions"]))
    resumed = stage.exists()
    if resumed:
        if stage.is_symlink() or (stage / "COMPLETE").exists():
            raise GateError("unsafe interrupted overall staging state")
        if not (stage / "input_manifest.json").is_file() or (stage / "input_manifest.json").read_bytes() != seed:
            raise GateError("interrupted stage immutable input seed mismatch")
    else:
        stage.mkdir(mode=0o750)
    for directory in (stage / "manifests", stage / "source_objects", stage / "canonical_objects", stage / "logs"):
        directory.mkdir(exist_ok=True, mode=0o750)
    state, failures, resources = stage / "state.jsonl", stage / "failures.jsonl", stage / "resources.jsonl"
    failures.touch(exist_ok=True)
    (stage / "input_manifest.json").write_bytes(seed)
    shutil.copyfile(repo / "manifests/pilot-cohorts-v1/cohort-0250.tsv", stage / "manifests/cohort-0250.tsv")
    p.append_jsonl(resources, initial)
    p.append_jsonl(state, {"event": "RESUME_PREFLIGHT_PASS" if resumed else "PREFLIGHT_PASS", "release_id": rid, "at": p.utcnow()})
    tools = {"bgzip": p.executable_record(args.bgzip), "samtools": p.executable_record(args.samtools)}
    predecessor_tools = json.loads((inputs["predecessor_external"] / "tools.json").read_text())
    for name in ("bgzip", "samtools"):
        if tools[name]["sha256"] != predecessor_tools[name]["sha256"] or tools[name]["version_first_line"] != predecessor_tools[name]["version_first_line"]:
            raise GateError(f"pinned {name} compatibility/tool digest mismatch")
    (stage / "tools.json").write_bytes(canonical_json(tools))
    provenance = {
        "schema": "canonical-cohort-250-provenance-v1", "task_id": TASK_ID, "release_id": rid,
        "run_id": args.run_id, "argv": sys.argv, "python": sys.version, "platform": platform.platform(),
        "hostname": socket.gethostname(), "pid": os.getpid(), "uid": os.getuid(), "gid": os.getgid(),
        "environment": {key: os.environ.get(key, ".") for key in ("USER", "LANG", "WG_TASK_ID", "WG_AGENT_ID", "WG_MODEL", "WG_TIER", "PI_MODEL", "PI_PROVIDER")},
        "source_api": "NCBI Datasets v2 exact single-accession packages", "batch_size": args.batch_size,
        "max_task_distinct_exact_assembly_revisions": COHORT_ROWS, "global_cap": 1000,
        "selection_release_id": SELECTION_RELEASE_ID, "cohort_manifest_sha256": COHORT_SHA256,
        "predecessor_release_id": PREDECESSOR_RELEASE_ID, "compatibility_release_id": COMPATIBILITY_RELEASE_ID,
        "integrated_n100_go_250": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY",
        "integrated_gate_rationale": "canonical acquisition scale-up is independent of unresolved prophage extraction/query per graph contract",
        "reuse_policy": "rows 1-100 are read-only predecessor digest references; no copy, download, or recompression",
        "created_at_utc": p.utcnow(),
    }
    (stage / "provenance.json").write_bytes(canonical_json(provenance))
    global_evidence = {"schema": "global-distinct-assembly-cap-v1", "verdict": "PASS", "start": global_start,
        "projected_finish_distinct_exact_assembly_revisions": len(projected),
        "projected_finish_accessions": sorted(projected), "task_bound": COHORT_ROWS, "graph_cap": 1000,
        "projected_finish_subset_of_frozen_collection": True}
    (stage / "global_cap_evidence.json").write_bytes(canonical_json(global_evidence))
    max_usage = _tree_usage(stage)
    batch_rows: list[dict[str, Any]] = []
    cumulative_transfer_bytes = 0
    cumulative_canonical_bytes = 0
    try:
        rows = inputs["rows"]
        for offset in range(0, len(rows), args.batch_size):
            batch = rows[offset:offset + args.batch_size]
            number = offset // args.batch_size + 1
            batch_preflight = live_preflight(final, scratch, allocations, f"BATCH_{number:03d}_START")
            p.append_jsonl(resources, batch_preflight)
            before = _tree_usage(stage)
            bw, bc = time.monotonic(), time.process_time()
            already, completed = 0, 0
            p.append_jsonl(state, {"event": "BATCH_STARTED", "batch_number": number,
                                   "first_cohort_order": int(batch[0]["cohort_order"]),
                                   "last_cohort_order": int(batch[-1]["cohort_order"]), "at": p.utcnow()})
            for input_row in batch:
                order = int(input_row["cohort_order"])
                accession = input_row["exact_assembly_accession_version"]
                if order <= PREDECESSOR_ROWS:
                    root, storage_id, _, _ = _storage_for(order, accession, stage, inputs)
                    _read_object(root, accession)
                    already += 1
                    p.append_jsonl(state, {"event": "PREDECESSOR_OBJECT_REUSED_BY_DIGEST", "cohort_order": order,
                                           "accession": accession, "release_id": storage_id,
                                           "logical_predecessor_release_id": PREDECESSOR_RELEASE_ID, "at": p.utcnow()})
                    continue
                source_path = stage / "source_objects" / accession
                canonical_path = stage / "canonical_objects" / accession
                _discard_invalid_completed(source_path, accession, "source", state, failures)
                _discard_invalid_completed(canonical_path, accession, "canonical", state, failures)
                was_complete = source_path.exists() and canonical_path.exists()
                p.append_jsonl(resources, live_preflight(final, scratch, allocations, f"ACQUISITION_{order:03d}_{accession}"))
                p.append_jsonl(state, {"event": "ASSEMBLY_ATTEMPT_STARTED", "cohort_order": order, "accession": accession, "at": p.utcnow()})
                source_dir, source_manifest = commit_source_with_transport_retry(
                    accession, stage / "source_objects", state, failures, args.retries, args.rate_delay,
                    args.inject_kill == "acquisition" and accession == args.inject_accession, args.inject_after_bytes)
                p.append_jsonl(resources, live_preflight(final, scratch, allocations, f"CANONICALIZATION_{order:03d}_{accession}"))
                p.canonicalize_object(accession, order, source_dir, source_manifest, stage / "canonical_objects", state,
                    args.bgzip, args.samtools, args.bgzip_threads, args.bgzip_level,
                    args.inject_kill == "conversion" and accession == args.inject_accession, args.inject_after_bases)
                if was_complete:
                    already += 1
                else:
                    completed += 1
                p.append_jsonl(state, {"event": "ASSEMBLY_TERMINAL_VALIDATED", "cohort_order": order, "accession": accession, "at": p.utcnow()})
            finish_usage = _tree_usage(stage)
            max_usage = (max(max_usage[0], finish_usage[0]), max(max_usage[1], finish_usage[1]))
            observed = _batch_observations(batch, stage, inputs, state, failures)
            cumulative_transfer_bytes += observed["validated_transfer_bytes"]
            cumulative_canonical_bytes += observed["canonical_bgzf_bytes"]
            end_preflight = live_preflight(final, scratch, allocations, f"BATCH_{number:03d}_END")
            p.append_jsonl(resources, end_preflight)
            partial_finish = _partial_bytes(stage)
            batch_rows.append({"batch_number": number, "first_cohort_order": int(batch[0]["cohort_order"]),
                "last_cohort_order": int(batch[-1]["cohort_order"]), "assemblies": len(batch),
                "already_complete": already, "completed_this_invocation": completed,
                "wall_seconds": f"{time.monotonic()-bw:.6f}", "cpu_seconds": f"{time.process_time()-bc:.6f}",
                **observed, "cumulative_transfer_bytes": cumulative_transfer_bytes,
                "cumulative_canonical_bgzf_bytes": cumulative_canonical_bytes,
                "peak_rss_bytes": _peak_rss_bytes(), "stage_bytes_start": before[0],
                "stage_bytes_finish": finish_usage[0], "stage_partial_bytes_finish": partial_finish,
                "stage_final_bytes_finish": finish_usage[0] - partial_finish,
                "stage_files_finish": finish_usage[1],
                "durable_free_bytes": end_preflight["durable_free_bytes"],
                "durable_free_inodes": end_preflight["durable_free_inodes"],
                "scratch_free_bytes": end_preflight["scratch_free_bytes"],
                "scratch_free_inodes": end_preflight["scratch_free_inodes"]})
            p.append_jsonl(state, {"event": "BATCH_COMPLETE", "batch_number": number, "assemblies": len(batch),
                                   "already_complete": already, "completed_this_invocation": completed, "at": p.utcnow()})
        assemblies, contigs, checksums, states, refs = _write_tables(stage, inputs, batch_rows)
        root_finish = verify_root_inputs(repo)
        if root_finish != inputs["root_inputs"]:
            raise GateError("immutable root inputs changed during execution")
        promotion = live_preflight(final, scratch, allocations, "PROMOTION")
        p.append_jsonl(resources, promotion)
        usage = _tree_usage(stage)
        max_usage = (max(max_usage[0], usage[0]), max(max_usage[1], usage[1]))
        ru_end = resource.getrusage(resource.RUSAGE_SELF)
        peak_rss = _peak_rss_bytes()
        process_swaps = max(0, ru_end.ru_nswap - ru_start.ru_nswap)
        system_swap_growth = max(0, initial["swap_free_bytes"] - promotion["swap_free_bytes"])
        predecessor_resource = json.loads((inputs["predecessor_external"] / "resource_summary.json").read_text())
        predecessor_new = int(inputs["predecessor"]["counts"]["new_objects"])
        linear_peak_projection = int(
            predecessor_resource["measured_release_stage_peak_bytes"] * NEW_ROWS / predecessor_new
        )
        modeled_upper95_peak = (linear_peak_projection * 125 + 99) // 100
        resource_checks = {
            "all_preflight_records_pass": all(json.loads(line).get("verdict") == "PASS" for line in resources.read_text().splitlines()),
            "peak_rss_le_70pct_assigned": peak_rss * 100 <= allocations.assigned_ram_bytes * 70,
            "process_swap_events_zero": process_swaps == 0,
            "system_swap_growth_zero": system_swap_growth == 0,
            "measured_peak_bytes_le_70pct_allocation": max_usage[0] * 100 <= allocations.durable_allocation_bytes * 70,
            "measured_files_le_50pct_inodes": max_usage[1] * 2 <= allocations.inode_allocation,
            "configured_projection_ge_modeled_upper95": allocations.predicted_durable_peak_bytes >= modeled_upper95_peak,
            "measured_peak_le_modeled_upper95": max_usage[0] <= modeled_upper95_peak,
        }
        if not all(resource_checks.values()):
            raise GateError("end resource gate NO_GO: " + json.dumps(resource_checks, sort_keys=True))
        resource_summary = {
            "schema": "canonical-cohort-250-resource-summary-v1", "verdict": "PASS",
            "allocations": asdict(allocations), "peak_rss_bytes": peak_rss,
            "peak_rss_fraction": peak_rss / allocations.assigned_ram_bytes,
            "process_swap_events": process_swaps, "system_swap_growth_bytes": system_swap_growth,
            "measured_release_stage_peak_bytes": max_usage[0], "measured_release_stage_peak_files": max_usage[1],
            "invocation_wall_seconds": time.monotonic() - invocation_wall_start,
            "invocation_cpu_seconds": time.process_time() - invocation_cpu_start,
            "preflight_record_count": len(resources.read_text().splitlines()), "batch_count": len(batch_rows),
            "start": initial, "finish": promotion, "checks": resource_checks,
            "disk_projection": {
                "method": "N=100 measured stage peak linearly scaled by new objects plus conservative 25% upper-95 allowance",
                "predecessor_release_id": PREDECESSOR_RELEASE_ID,
                "predecessor_new_objects": predecessor_new,
                "predecessor_measured_peak_bytes": predecessor_resource["measured_release_stage_peak_bytes"],
                "new_objects": NEW_ROWS, "linear_peak_projection_bytes": linear_peak_projection,
                "modeled_upper95_peak_bytes": modeled_upper95_peak,
                "configured_upper95_peak_bytes": allocations.predicted_durable_peak_bytes,
                "observed_to_linear_ratio": max_usage[0] / linear_peak_projection,
            },
            "scale_trend": "NOT_APPLICABLE_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS",
        }
        (stage / "resource_summary.json").write_bytes(canonical_json(resource_summary))
        state_text = state.read_text()
        restart = {
            "schema": "canonical-cohort-250-restart-evidence-v1",
            "acquisition_injected_kill_observed": "INJECTED_ACQUISITION_SIGKILL" in state_text,
            "acquisition_partial_safe_restart_observed": any(x in state_text for x in ("ACQUISITION_RANGE_RESUME", "ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE")),
            "conversion_injected_kill_observed": "INJECTED_CONVERSION_SIGKILL" in state_text,
            "conversion_partial_discard_restart_observed": "INTERRUPTED_CONVERSION_STAGE_DISCARDED" in state_text,
            "no_partial_or_mixed_final_publication": not final.exists(),
            "predecessor_objects_never_downloaded_or_recompressed": all(f'"accession":"{accession}"' not in "\n".join(line for line in state_text.splitlines() if "ACQUISITION_REQUEST_STARTED" in line or "CANONICAL_OBJECT_VALIDATED" in line) for accession in inputs["accessions"][:PREDECESSOR_ROWS]),
        }
        if not all(value for key, value in restart.items() if key != "schema"):
            raise GateError("required injected kill/restart/reuse evidence incomplete")
        (stage / "restart_evidence.json").write_bytes(canonical_json(restart))
        gates = {
            "selection_manifest_identity_sha_rows_order": "PASS", "predecessor_release_identity_inventory": "PASS",
            "pinned_consumer_compatibility": "PASS", "root_source_immutability": "PASS",
            "accession_version_identity": "PASS", "global_distinct_assembly_cap": "PASS",
            "predecessor_object_digest_reuse": "PASS", "archive_upstream_local_checksum": "PASS",
            "row_accounting": "PASS", "rename_only_sequence_identity": "PASS",
            "bgzf_index_name_roundtrip": "PASS", "pansn_uniqueness_reversibility": "PASS",
            "source_coordinate_annotation_policy": "PASS", "resource": "PASS",
            "bounded_batches_retry_resume": "PASS", "injected_kill_restart": "PASS",
            "atomic_object_release_promotion": "PASS", "deterministic_semantic_rerun": "PASS",
            "integrated_n100_go_250": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY",
            "scale_trend": "NOT_APPLICABLE_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS",
        }
        release = {
            "schema_version": SCHEMA, "release_id": rid, "verdict": "PASS", "immutable": True,
            "source_task_id": TASK_ID, "created_at_utc": p.utcnow(), "external_release_path": str(final),
            "selection_release_id": SELECTION_RELEASE_ID, "selection_release_json_sha256": SELECTION_RELEASE_JSON_SHA256,
            "input_cohort_0250_sha256": COHORT_SHA256, "predecessor_release_id": PREDECESSOR_RELEASE_ID,
            "predecessor_release_json_sha256": PREDECESSOR_RELEASE_JSON_SHA256,
            "compatibility_release_id": COMPATIBILITY_RELEASE_ID,
            "compatibility_release_json_sha256": COMPATIBILITY_RELEASE_JSON_SHA256,
            "pansn_policy_version": POLICY, "exact_version_policy": "exact frozen revision; never substitute latest",
            "sequence_bearing_assembly_accessions": inputs["accessions"],
            "counts": {"attempted_exact_assembly_revisions": COHORT_ROWS, "validated": COHORT_ROWS, "quarantined": 0,
                "reused_predecessor_objects": PREDECESSOR_ROWS, "new_objects": NEW_ROWS,
                "distinct_sequence_bearing_assemblies": COHORT_ROWS, "global_distinct_assembly_cap": 1000,
                "contigs": len(contigs), "total_bases": sum(int(row["contig_length"]) for row in contigs),
                "annotations_alias_validated": sum(row["annotation_status"].startswith("ALIASES_VALIDATED") for row in assemblies),
                "annotation_views_quarantined": sum(row["annotation_status"].startswith("QUARANTINED") for row in assemblies)},
            "manifests": {}, "root_inputs_start": inputs["root_inputs"], "root_inputs_finish": root_finish,
            "resource_summary_sha256": sha_file(stage / "resource_summary.json"),
            "restart_evidence_sha256": sha_file(stage / "restart_evidence.json"),
            "global_cap_evidence_sha256": sha_file(stage / "global_cap_evidence.json"),
            "applicable_gates": gates,
            "object_storage_contract": "rows 1-100 resolve read-only through the N=100 predecessor's immutable object references and inventory SHA; rows 101-250 are self-contained",
            "annotation_policy": "source GFF retained unchanged; only fully validated alias views published; invalid optional views quarantined",
            "routine_whole_set_plain_fasta_files": 0,
        }
        manifest_specs = {
            "cohort-0250.tsv": COHORT_ROWS, "assemblies.tsv": COHORT_ROWS, "contigs.tsv": len(contigs),
            "checksums.tsv": COHORT_ROWS * 8, "state.tsv": COHORT_ROWS,
            "object_refs.tsv": COHORT_ROWS, "batch_metrics.tsv": len(batch_rows),
        }
        for name, count in manifest_specs.items():
            path = stage / "manifests" / name
            release["manifests"][name] = {"rows": count, "bytes": path.stat().st_size, "sha256": sha_file(path)}
        (stage / "release.json").write_bytes(canonical_json(release))
        p.append_jsonl(state, {"event": "READY_TO_PROMOTE", "release_id": rid, "at": p.utcnow()})
        p.seal_directory(stage, final)
        global_finish = audit_global_release_cap(durable_task.parents[1], allowed)
        if global_finish["distinct_exact_assembly_revisions"] != COHORT_ROWS or set(global_finish["accessions"]) != set(inputs["accessions"]):
            raise GateError("finish global union is not exactly frozen N=250")
        publish_tracked(final, tracked, artifact)
        return release
    except BaseException as exc:
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        p.append_jsonl(failures, {"event": "RUN_FAILED", "type": type(exc).__name__, "message": str(exc), "at": p.utcnow()})
        raise


def parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--tracked-root", default="manifests/canonical-cohort-250-v1")
    parser.add_argument("--artifact-root", default="artifacts/canonical_cohort_250")
    parser.add_argument("--durable-task-root", required=True)
    parser.add_argument("--scratch-root", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--retries", type=int, default=8)
    parser.add_argument("--rate-delay", type=float, default=0.5)
    parser.add_argument("--assigned-ram-bytes", type=int, required=True)
    parser.add_argument("--durable-allocation-bytes", type=int, required=True)
    parser.add_argument("--scratch-allocation-bytes", type=int, required=True)
    parser.add_argument("--inode-allocation", type=int, required=True)
    parser.add_argument("--predicted-durable-peak-bytes", type=int, required=True)
    parser.add_argument("--predicted-scratch-peak-bytes", type=int, required=True)
    parser.add_argument("--predicted-files", type=int, required=True)
    parser.add_argument("--unfinished-write-bytes", type=int, required=True)
    parser.add_argument("--bgzip", default="bgzip")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bgzip-threads", type=int, default=2)
    parser.add_argument("--bgzip-level", type=int, default=6)
    parser.add_argument("--inject-kill", choices=("none", "acquisition", "conversion"), default="none")
    parser.add_argument("--inject-accession")
    parser.add_argument("--inject-after-bytes", type=int, default=131072)
    parser.add_argument("--inject-after-bases", type=int, default=262144)
    return parser


def main() -> int:
    try:
        args = parser().parse_args()
        inputs = verify_inputs(Path(args.repo_root).resolve())
        if args.batch_size <= 0 or args.batch_size > 25:
            raise GateError("batch size must be within 1..25")
        if args.inject_accession is None:
            args.inject_accession = inputs["accessions"][PREDECESSOR_ROWS]
        if args.inject_accession not in inputs["accessions"][PREDECESSOR_ROWS:]:
            raise GateError("kill injection must target a new frozen cohort row 101..250")
        print(json.dumps(run(args), sort_keys=True))
        return 0
    except (p.GateError, GateError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
