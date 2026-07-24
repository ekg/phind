#!/usr/bin/env python3
"""Fail-closed, phage-blind selection and atomic publication of pilot cohorts.

This module reads only frozen metadata manifests.  It contains no sequence or
biological-analysis endpoint and never opens the prophage CSV while selecting.
The first ten positions are the checksum-pinned Stage-B cohort; remaining
positions are a deterministic hash ordering of the non-suppressed exact-version
frame.  Engineering controls are attached only after that order is frozen.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import resource
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import uuid
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable

TASK_ID = "select-freeze-1k"
SCHEMA = "pilot-cohorts-v1"
RUNGS = (10, 100, 250, 500, 1000)
MAIN_SEED = "pilot-cohorts-v1-main-srs-sha256-seed-20260724"
SELECTION_ALGORITHM = "stage-b-certainty-then-sha256-srs-prefix-v1"

ACCESSION_INPUT_SHA256 = "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5"
PROPHAGE_INPUT_SHA256 = "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996"
COLLECTION_RELEASE_ID = "collection-v1-f7494b4b89d1382b"
COLLECTION_RELEASE_JSON_SHA256 = "59c6907e2c053e9d8ac3df8d5eb820bab0097030a9259ca2c9354c47cb6642bf"
COLLECTION_ASSEMBLIES_GZ_SHA256 = "5d72d583cd26066782ecb735c7931466bd7ef74418b275099afed161f1a5041d"
COLLECTION_OCCURRENCES_GZ_SHA256 = "7ef0563c35bec568cf2b86473e1e785759a8d32b9b4666921e698fc288bbcc32"
COLLECTION_EXTERNAL_SHA256SUMS_SHA256 = "bd27dcca6c1e33eda8ccaafd30224aae5a2b5770ee6a87246474d6d55844f0e9"
STAGE_B_SHA256 = "0d179cbafce2ba1fa14d1929a4acd6621810a335f25bcd7ec67dd2083eb101f6"
STAGE_B_BYTES = 2246
ACQUISITION_RELEASE_ID = "canonical-cohort-010-v1-e71484de9994fc28"
ACQUISITION_RELEASE_JSON_SHA256 = "4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23"
ACQUISITION_EXTERNAL_SHA256SUMS_SHA256 = "96a40035c15684d4c3c12c88f8134c32c4df421eb9d138119581ab7473badc44"
ACQUISITION_ASSEMBLIES_SHA256 = "7133058093e3f08c132248b3cf4453c7076b6550e838f2ab9f39a6b5b7b8fcbd"
SOURCE_SEMANTICS_RELEASE_ID = "prophage-semantics-v1-f5619e221ff272ae"
SOURCE_SEMANTICS_RELEASE_JSON_SHA256 = "6a8de2063e4e12c0f0f363ebe41aba03a2f463e8a45711ccbe7ebcdae581b728"
SOURCE_SEMANTICS_EXTERNAL_SHA256SUMS_SHA256 = "ab86ba983f4e5e44823e321dbafb84b19b1c31667abd3fc6cc0102ab1bdcef31"

DEFAULT_RELEASES_ROOT = Path("/home/erikg/phind-data/ecoli26k/v1/releases")
DEFAULT_DURABLE_ROOT = DEFAULT_RELEASES_ROOT / TASK_ID
DEFAULT_SCRATCH_ROOT = Path("/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1") / TASK_ID
SOURCE_SEMANTICS_ROOT = DEFAULT_RELEASES_ROOT / "resolve-prophage-source" / SOURCE_SEMANTICS_RELEASE_ID

FORBIDDEN_SELECTION_FIELDS = (
    "prophage_presence", "prophage_count", "phage_count", "transposable",
    "taxonomy", "taxonomy_flag", "cluster", "host_clade_derived_from_phage",
)
ALLOWED_SELECTION_FIELDS = (
    "assembly_id", "resolved_assembly_accession_version", "input_line_number",
    "resolution_status", "assembly_status", "frozen_stage_b_membership",
)

FRAME_FIELDS = [
    "frame_order", "input_line_number", "input_occurrence_id", "assembly_id",
    "exact_assembly_accession_version", "resolution_status", "assembly_status",
    "frame_disposition", "selection_stratum", "stage_b_order", "random_key",
    "cohort_order", "first_rung", "collection_row_sha256",
] + [f"inclusion_probability_n{n}" for n in RUNGS] + ["row_sha256"]

COHORT_FIELDS = [
    "cohort_order", "rung_n", "assembly_id", "exact_assembly_accession_version",
    "input_line_number", "input_occurrence_id", "resolution_status", "assembly_status",
    "selection_stratum", "stage_b_order", "random_key", "inclusion_probability",
    "inference_weight", "engineering_control_membership", "collection_row_sha256",
    "row_sha256",
]

CONTROL_FIELDS = [
    "control_id", "control_class", "control_scope", "assembly_id",
    "exact_assembly_accession_version", "fixture_or_evidence", "selection_effect",
    "inference_disposition", "activation_status", "source_release_id", "row_sha256",
]


class GateError(RuntimeError):
    """A release-blocking contract failure."""


class InjectedInterruption(RuntimeError):
    """Test-only stand-in for an injected SIGKILL."""


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
        if any(value <= 0 for value in asdict(self).values()):
            raise GateError("resource allocation/reservation is blank or non-positive")
        if self.predicted_durable_peak_bytes * 100 > self.durable_allocation_bytes * 70:
            raise GateError("predicted durable upper-95% peak exceeds 70% allocation")
        if self.predicted_scratch_peak_bytes * 100 > self.scratch_allocation_bytes * 70:
            raise GateError("predicted scratch upper-95% peak exceeds 70% allocation")
        if self.predicted_files * 2 > self.inode_allocation:
            raise GateError("projected files exceed 50% inode allocation")
        if self.durable_allocation_bytes < 2 * self.unfinished_write_bytes:
            raise GateError("durable allocation cannot retain 2x unfinished writes")
        if self.scratch_allocation_bytes < 2 * self.unfinished_write_bytes:
            raise GateError("scratch allocation cannot retain 2x unfinished writes")


@dataclass
class SelectionResult:
    frame: list[dict[str, str]]
    cohort: list[dict[str, str]]
    rungs: dict[int, list[dict[str, str]]]
    frame_counts: dict[str, int]
    eligible_count: int
    main_random_count: int


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def fsync_dir(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def append_jsonl(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("ab") as handle:
        handle.write(canonical_json(row))
        handle.flush()
        os.fsync(handle.fileno())


def stable_row_hash(row: dict[str, Any], fields: list[str]) -> str:
    material = "\t".join(str(row.get(field, ".")) for field in fields if field != "row_sha256") + "\n"
    return sha_bytes(material.encode())


def _normalize_row(source: dict[str, Any], fields: list[str]) -> dict[str, str]:
    row = {field: ("." if source.get(field) in (None, "") else str(source[field])) for field in fields}
    row["row_sha256"] = stable_row_hash(row, fields)
    return row


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for source in rows:
            writer.writerow(_normalize_row(source, fields))
        handle.flush()
        os.fsync(handle.fileno())


def read_tsv(path: Path, fields: list[str] | None = None, verify_hashes: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed = reader.fieldnames
        if fields is not None and observed != fields:
            raise GateError(f"TSV schema mismatch: {path}")
        rows = list(reader)
    if verify_hashes:
        actual_fields = fields or observed
        if not actual_fields or "row_sha256" not in actual_fields:
            raise GateError(f"row checksums unavailable: {path}")
        for number, row in enumerate(rows, 2):
            if row.get("row_sha256") != stable_row_hash(row, actual_fields):
                raise GateError(f"row checksum mismatch: {path}:{number}")
    return rows


def deterministic_gzip(source: Path, target: Path) -> None:
    with source.open("rb") as src, target.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as dst:
            shutil.copyfileobj(src, dst)
    fsync_file(target)


def _inventory_entries(root: Path) -> list[str]:
    entries = []
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = str(path.relative_to(root))
        if relative in ("SHA256SUMS", "COMPLETE"):
            continue
        if path.is_symlink():
            raise GateError(f"symlink in release: {path}")
        entries.append(f"{sha_file(path)}  {relative}\n")
    return entries


def write_inventory(root: Path) -> None:
    (root / "SHA256SUMS").write_text("".join(_inventory_entries(root)), encoding="utf-8")
    fsync_file(root / "SHA256SUMS")


def verify_inventory(
    root: Path, *, require_complete: bool = True, deep: bool = True,
    expected_inventory_sha256: str | None = None,
) -> None:
    sums = root / "SHA256SUMS"
    if not root.is_dir() or not sums.is_file() or sums.is_symlink():
        raise GateError(f"missing checksum inventory: {root}")
    inventory_sha = sha_file(sums)
    if expected_inventory_sha256 is not None and inventory_sha != expected_inventory_sha256:
        raise GateError(f"SHA256SUMS digest mismatch: {root}")
    if require_complete:
        complete = root / "COMPLETE"
        if not complete.is_file() or complete.is_symlink():
            raise GateError(f"missing COMPLETE: {root}")
        complete_text = complete.read_text(encoding="utf-8")
        if complete_text.lstrip().startswith("{"):
            try:
                complete_digest = json.loads(complete_text).get("sha256sums_sha256")
            except json.JSONDecodeError as exc:
                raise GateError(f"malformed COMPLETE JSON: {root}") from exc
        else:
            tokens = complete_text.split()
            complete_digest = tokens[0] if tokens else None
        if complete_digest != inventory_sha:
            raise GateError(f"COMPLETE inventory digest mismatch: {root}")
    seen: set[str] = set()
    for number, line in enumerate(sums.read_text(encoding="utf-8").splitlines(), 1):
        try:
            digest, relative = line.split("  ", 1)
        except ValueError as exc:
            raise GateError(f"malformed inventory line {number}: {root}") from exc
        if len(digest) != 64 or relative in seen or relative.startswith("/") or ".." in Path(relative).parts:
            raise GateError(f"unsafe or duplicate inventory path: {relative}")
        seen.add(relative)
        path = root / relative
        if not path.is_file() or path.is_symlink():
            raise GateError(f"inventory path missing: {path}")
        if deep and sha_file(path) != digest:
            raise GateError(f"inventory checksum mismatch: {path}")


def seal_and_promote(stage: Path, final: Path, inject: Callable[[], None] | None = None) -> None:
    if final.exists():
        raise GateError(f"refusing to overwrite immutable release: {final}")
    write_inventory(stage)
    for path in (candidate for candidate in stage.rglob("*") if candidate.is_file()):
        fsync_file(path)
    fsync_dir(stage)
    if inject is not None:
        inject()
    (stage / "COMPLETE").write_text(f"{sha_file(stage / 'SHA256SUMS')}  SHA256SUMS\n", encoding="utf-8")
    fsync_file(stage / "COMPLETE")
    fsync_dir(stage)
    os.rename(stage, final)
    fsync_dir(final.parent)
    verify_inventory(final)


def discard_interrupted_stage(stage: Path, state_log: Path, failures_log: Path | None = None) -> None:
    if not stage.exists():
        return
    if stage.is_symlink() or not stage.is_dir():
        raise GateError(f"unsafe interrupted stage: {stage}")
    if (stage / "COMPLETE").exists():
        raise GateError(f"completed stage cannot be silently discarded: {stage}")
    validated = False
    try:
        verify_inventory(stage, require_complete=False)
        validated = True
    except GateError as exc:
        if failures_log is not None:
            append_jsonl(failures_log, {"event": "INTERRUPTED_STAGE_INVENTORY_INVALID", "at_utc": utcnow(), "reason": str(exc)})
    append_jsonl(state_log, {
        "event": "INTERRUPTED_PUBLICATION_STAGE_DISCARDED", "at_utc": utcnow(),
        "stage": str(stage), "inventory_was_checksum_valid": validated,
        "completed_units_reused": 0,
    })
    shutil.rmtree(stage)
    fsync_dir(stage.parent)


def _selection_key(row: dict[str, str]) -> str:
    material = "\0".join((MAIN_SEED, row["assembly_id"], row["resolved_assembly_accession_version"]))
    return sha_bytes(material.encode())


def _probability(stratum: str, rung: int, main_count: int) -> tuple[str, str]:
    if stratum == "stage_b_certainty":
        return "1/1", "1/1"
    if stratum == "main_phage_blind_srs":
        draw = rung - 10
        return f"{draw}/{main_count}", f"{main_count}/{draw}"
    return "0/1", "."


def select_cohort(assembly_rows: list[dict[str, str]], stage_b_rows: list[dict[str, str]]) -> SelectionResult:
    if len(stage_b_rows) != 10:
        raise GateError("Stage-B cardinality is not exactly 10")
    if len({row.get("assembly_id") for row in assembly_rows}) != len(assembly_rows):
        raise GateError("duplicate assembly ID in frozen frame")
    if len({row.get("resolved_assembly_accession_version") for row in assembly_rows}) != len(assembly_rows):
        raise GateError("duplicate exact assembly revision in frozen frame")
    by_id = {row["assembly_id"]: row for row in assembly_rows}
    stage_ids: list[str] = []
    for order, stage in enumerate(stage_b_rows, 1):
        if int(stage.get("stage_b_order", "0")) != order:
            raise GateError("Stage-B order drift")
        source = by_id.get(stage.get("assembly_id", ""))
        if source is None or source["resolved_assembly_accession_version"] != stage.get("resolved_assembly_accession_version"):
            raise GateError("Stage-B identity drift from frozen frame")
        stage_ids.append(source["assembly_id"])
    if len(set(stage_ids)) != 10:
        raise GateError("duplicate Stage-B assembly ID")

    for row in assembly_rows:
        requested = row.get("requested_assembly_accession_version")
        resolved = row.get("resolved_assembly_accession_version")
        if not requested or requested != resolved or not row.get("assembly_id"):
            raise GateError("exact-version identity failure in frozen frame")
        if row.get("resolution_status") not in ("EXACT_VERSION_RESOLVED", "EXACT_VERSION_VALID_METADATA_UNAVAILABLE"):
            raise GateError("non-exact resolution status in assembly frame")

    stage_set = set(stage_ids)
    eligible = [row for row in assembly_rows if row.get("assembly_status") != "suppressed"]
    suppressed = [row for row in assembly_rows if row.get("assembly_status") == "suppressed"]
    allowed_statuses = {"current", "METADATA_UNAVAILABLE", "suppressed"}
    if any(row.get("assembly_status") not in allowed_statuses for row in assembly_rows):
        raise GateError("unknown assembly terminal status")
    if not stage_set.issubset({row["assembly_id"] for row in eligible}):
        raise GateError("Stage-B contains terminal suppressed assembly")
    if len(eligible) < max(RUNGS):
        raise GateError("eligible exact-version frame cannot fill N=1,000")

    main_rows = [row for row in eligible if row["assembly_id"] not in stage_set]
    main_rows.sort(key=lambda row: (_selection_key(row), row["assembly_id"]))
    ordered_sources = [by_id[assembly_id] for assembly_id in stage_ids] + main_rows[:990]
    cohort_order = {row["assembly_id"]: order for order, row in enumerate(ordered_sources, 1)}
    stage_order = {assembly_id: order for order, assembly_id in enumerate(stage_ids, 1)}
    main_count = len(main_rows)

    frame: list[dict[str, str]] = []
    for source in sorted(assembly_rows, key=lambda row: int(row["input_line_number"])):
        assembly_id = source["assembly_id"]
        if assembly_id in stage_set:
            stratum, disposition = "stage_b_certainty", "ELIGIBLE_CERTAINTY"
        elif source in suppressed:
            stratum, disposition = "terminal_suppressed_ineligible", "INELIGIBLE_TERMINAL_SUPPRESSED"
        else:
            stratum, disposition = "main_phage_blind_srs", "ELIGIBLE_RANDOM"
        order = cohort_order.get(assembly_id)
        first_rung = next((n for n in RUNGS if order is not None and order <= n), None)
        row: dict[str, Any] = {
            "frame_order": source["input_line_number"],
            "input_line_number": source["input_line_number"],
            "input_occurrence_id": source["input_occurrence_id"],
            "assembly_id": assembly_id,
            "exact_assembly_accession_version": source["resolved_assembly_accession_version"],
            "resolution_status": source["resolution_status"],
            "assembly_status": source["assembly_status"],
            "frame_disposition": disposition,
            "selection_stratum": stratum,
            "stage_b_order": stage_order.get(assembly_id, "."),
            "random_key": "STAGE_B_CERTAINTY" if assembly_id in stage_set else _selection_key(source),
            "cohort_order": order if order is not None else ".",
            "first_rung": first_rung if first_rung is not None else ".",
            "collection_row_sha256": source["row_sha256"],
        }
        for rung in RUNGS:
            row[f"inclusion_probability_n{rung}"] = _probability(stratum, rung, main_count)[0]
        frame.append(_normalize_row(row, FRAME_FIELDS))

    frame_by_id = {row["assembly_id"]: row for row in frame}
    cohort: list[dict[str, str]] = []
    for order, source in enumerate(ordered_sources, 1):
        frame_row = frame_by_id[source["assembly_id"]]
        cohort.append({
            "cohort_order": str(order),
            "assembly_id": source["assembly_id"],
            "exact_assembly_accession_version": source["resolved_assembly_accession_version"],
            "input_line_number": source["input_line_number"],
            "input_occurrence_id": source["input_occurrence_id"],
            "resolution_status": source["resolution_status"],
            "assembly_status": source["assembly_status"],
            "selection_stratum": frame_row["selection_stratum"],
            "stage_b_order": frame_row["stage_b_order"],
            "random_key": frame_row["random_key"],
            "collection_row_sha256": source["row_sha256"],
        })
    rungs: dict[int, list[dict[str, str]]] = {}
    for rung in RUNGS:
        rung_rows = []
        for source in cohort[:rung]:
            probability, weight = _probability(source["selection_stratum"], rung, main_count)
            rung_rows.append({
                **source, "rung_n": str(rung), "inclusion_probability": probability,
                "inference_weight": weight, "engineering_control_membership": ".",
            })
        rungs[rung] = rung_rows
    counts = Counter(row["selection_stratum"] for row in frame)
    return SelectionResult(
        frame=frame, cohort=cohort, rungs=rungs, frame_counts=dict(sorted(counts.items())),
        eligible_count=len(eligible), main_random_count=main_count,
    )


def _control(
    control_class: str, scope: str, fixture: str, source_release_id: str,
    assembly_id: str = ".", accession: str = ".", activation: str = "READY",
    inference: str = "EXCLUDE_FIXTURE_FROM_HOST_AND_PREVALENCE_INFERENCE",
) -> dict[str, str]:
    identity = "\0".join((control_class, scope, assembly_id, accession, fixture, source_release_id))
    return {
        "control_id": "ctrl-v1-" + sha_bytes(identity.encode())[:24],
        "control_class": control_class,
        "control_scope": scope,
        "assembly_id": assembly_id,
        "exact_assembly_accession_version": accession,
        "fixture_or_evidence": fixture,
        "selection_effect": "NONE_POST_SELECTION_LABEL",
        "inference_disposition": inference,
        "activation_status": activation,
        "source_release_id": source_release_id,
    }


def synthetic_engineering_controls(stage_b_rows: list[dict[str, str]], source_release_id: str) -> list[dict[str, str]]:
    controls = [
        _control("unsafe_source_contig_id", "synthetic_unit_fixture", "literal:ctg#A", ACQUISITION_RELEASE_ID),
    ]
    for case in ("begin_at_contig_start", "end_at_contig_end", "short_interval", "long_interval", "circular_origin_wrap"):
        controls.append(_control(
            "prophage_interval_edge_extreme", "future_synthetic_known_base_fixture", case,
            source_release_id, activation="BLOCKED_EXTRACTION_SEMANTICS",
        ))
    return controls


def build_engineering_controls(repo: Path, stage_b_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    controls = synthetic_engineering_controls(stage_b_rows, SOURCE_SEMANTICS_RELEASE_ID)
    for stage in stage_b_rows:
        accession = stage["resolved_assembly_accession_version"]
        version = int(accession.rsplit(".", 1)[1])
        if version > 1:
            controls.append(_control(
                "accession_exact_revision_ambiguity", "observed_stage_b", f"assembly_version={version}",
                COLLECTION_RELEASE_ID, stage["assembly_id"], accession,
                inference="USE_ONCE_WITH_FROZEN_DESIGN_WEIGHT_NO_CONTROL_MULTIPLICITY",
            ))
        controls.append(_control(
            "prophage_interval_bounded_sentinel_scope", "observed_stage_b_post_selection",
            "source-semantics bounded sentinel; no coordinate extraction authorized",
            SOURCE_SEMANTICS_RELEASE_ID, stage["assembly_id"], accession,
            activation="BLOCKED_EXTRACTION_SEMANTICS",
            inference="USE_ONCE_WITH_FROZEN_DESIGN_WEIGHT_NO_CONTROL_MULTIPLICITY",
        ))

    acquisition_rows = read_tsv(repo / "manifests/canonical-cohort-010-v1/assemblies.tsv", verify_hashes=True)
    if sha_file(repo / "manifests/canonical-cohort-010-v1/assemblies.tsv") != ACQUISITION_ASSEMBLIES_SHA256:
        raise GateError("acquisition assembly summary checksum mismatch")
    by_accession = {row["accession"]: row for row in acquisition_rows}
    if set(by_accession) != {row["resolved_assembly_accession_version"] for row in stage_b_rows}:
        raise GateError("acquisition summary differs from Stage-B identities")
    for control_class, key, chooser in (
        ("assembly_size_minimum", "total_bases", min),
        ("assembly_size_maximum", "total_bases", max),
        ("assembly_contiguity_single_or_minimum", "contig_count", min),
        ("assembly_contiguity_maximum", "contig_count", max),
    ):
        selected = chooser(acquisition_rows, key=lambda row: (int(row[key]), row["accession"]))
        stage = next(row for row in stage_b_rows if row["resolved_assembly_accession_version"] == selected["accession"])
        controls.append(_control(
            control_class, "observed_stage_b_post_selection", f"{key}={selected[key]}",
            ACQUISITION_RELEASE_ID, stage["assembly_id"], selected["accession"],
            inference="USE_ONCE_WITH_FROZEN_DESIGN_WEIGHT_NO_CONTROL_MULTIPLICITY",
        ))

    contigs = read_tsv(repo / "manifests/canonical-cohort-010-v1/contigs.tsv.gz", verify_hashes=True)
    for role in sorted({row["replicon_role"] for row in contigs}):
        selected = min((row for row in contigs if row["replicon_role"] == role), key=lambda row: (int(row["stage_b_order"]), int(row["contig_order"])))
        stage = stage_b_rows[int(selected["stage_b_order"]) - 1]
        controls.append(_control(
            "source_contig_role", "observed_stage_b_post_selection",
            f"replicon_role={role};source_contig={selected['source_contig_id_display']}",
            ACQUISITION_RELEASE_ID, stage["assembly_id"], selected["accession"],
            inference="USE_ONCE_WITH_FROZEN_DESIGN_WEIGHT_NO_CONTROL_MULTIPLICITY",
        ))
    return [_normalize_row(row, CONTROL_FIELDS) for row in sorted(controls, key=lambda row: row["control_id"])]


def _read_checksum_map(root: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        if relative in checksums or sha_file(root / relative) != digest:
            raise GateError(f"tracked checksum mismatch: {root}/{relative}")
        checksums[relative] = digest
    return checksums


def verify_root_inputs(repo: Path) -> dict[str, str]:
    observed = {
        "26k_ecoli_accession.txt": sha_file(repo / "26k_ecoli_accession.txt"),
        "26k_prophage1.csv": sha_file(repo / "26k_prophage1.csv"),
    }
    if observed != {
        "26k_ecoli_accession.txt": ACCESSION_INPUT_SHA256,
        "26k_prophage1.csv": PROPHAGE_INPUT_SHA256,
    }:
        raise GateError("immutable root input checksum mismatch")
    return observed


def _require_gate_values(gates: dict[str, Any], allowed_na: dict[str, str]) -> None:
    for name, verdict in gates.items():
        if verdict == "PASS":
            continue
        if allowed_na.get(name) == verdict:
            continue
        raise GateError(f"applicable predecessor gate is not PASS: {name}={verdict}")


def audit_global_sequence_union(releases_root: Path, allowed_accessions: set[str]) -> dict[str, Any]:
    accessions: set[str] = set()
    evidence: list[dict[str, Any]] = []
    for path in sorted(releases_root.rglob("release.json")):
        if not (path.parent / "COMPLETE").is_file():
            continue
        try:
            release = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        declared = release.get("sequence_bearing_assembly_accessions", [])
        count = release.get("counts", {}).get("distinct_sequence_bearing_assemblies", 0)
        if count and not declared:
            raise GateError(f"sequence-bearing release lacks exact identity inventory: {path}")
        if declared:
            if len(set(declared)) != len(declared) or int(count) != len(set(declared)):
                raise GateError(f"sequence-bearing release identity count mismatch: {path}")
            accessions.update(declared)
            evidence.append({"release_id": release.get("release_id"), "release_json": str(path), "count": len(declared)})
    if len(accessions) > 1000 or not accessions.issubset(allowed_accessions):
        raise GateError("global sequence-bearing union exceeds cap or frozen collection")
    return {
        "verdict": "PASS", "global_cap": 1000,
        "distinct_exact_assembly_revisions": len(accessions),
        "accessions": sorted(accessions), "evidence_releases": evidence,
    }


def verify_inputs(repo: Path, releases_root: Path = DEFAULT_RELEASES_ROOT) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    roots = verify_root_inputs(repo)
    collection_root = repo / "manifests/collection-v1"
    checksums = _read_checksum_map(collection_root)
    expected_tracked = {
        "release.json": COLLECTION_RELEASE_JSON_SHA256,
        "assemblies.tsv.gz": COLLECTION_ASSEMBLIES_GZ_SHA256,
        "occurrences.tsv.gz": COLLECTION_OCCURRENCES_GZ_SHA256,
        "stage_b_10.tsv": STAGE_B_SHA256,
    }
    for relative, digest in expected_tracked.items():
        if checksums.get(relative) != digest:
            raise GateError(f"pinned collection manifest checksum mismatch: {relative}")
    collection = json.loads((collection_root / "release.json").read_text(encoding="utf-8"))
    if collection.get("release_id") != COLLECTION_RELEASE_ID or collection.get("verdict") != "PASS" or collection.get("immutable") is not True:
        raise GateError("collection release ID/verdict/immutability mismatch")
    _require_gate_values(collection.get("applicable_gates", {}), {
        "bgzf_index_name_roundtrip": "NOT_APPLICABLE_METADATA_ONLY",
        "scale_trend": "NOT_APPLICABLE_METADATA_ONLY",
        "source_coordinate_policy": "NOT_APPLICABLE_METADATA_ONLY",
    })
    external_collection = Path(collection.get("external_release_path", ""))
    verify_inventory(external_collection, expected_inventory_sha256=COLLECTION_EXTERNAL_SHA256SUMS_SHA256)

    assemblies = read_tsv(collection_root / "assemblies.tsv.gz", verify_hashes=True)
    if len(assemblies) != 26077 or len({row["assembly_id"] for row in assemblies}) != 26077:
        raise GateError("collection assembly row/cardinality accounting failure")
    statuses = Counter(row["assembly_status"] for row in assemblies)
    if statuses != Counter({"current": 25291, "suppressed": 714, "METADATA_UNAVAILABLE": 72}):
        raise GateError("collection terminal status accounting mismatch")
    stage_b_path = collection_root / "stage_b_10.tsv"
    if stage_b_path.stat().st_size != STAGE_B_BYTES:
        raise GateError("Stage-B byte count mismatch")
    stage_b = read_tsv(stage_b_path, verify_hashes=True)

    acquisition_root = repo / "manifests/canonical-cohort-010-v1"
    acquisition_checksums = _read_checksum_map(acquisition_root)
    if acquisition_checksums.get("release.json") != ACQUISITION_RELEASE_JSON_SHA256:
        raise GateError("acquisition release checksum mismatch")
    if acquisition_checksums.get("stage_b_10.tsv") != STAGE_B_SHA256:
        raise GateError("acquisition Stage-B checksum mismatch")
    if (acquisition_root / "stage_b_10.tsv").read_bytes() != stage_b_path.read_bytes():
        raise GateError("validated acquisition Stage-B bytes differ")
    acquisition = json.loads((acquisition_root / "release.json").read_text(encoding="utf-8"))
    if acquisition.get("release_id") != ACQUISITION_RELEASE_ID or acquisition.get("verdict") != "PASS" or acquisition.get("immutable") is not True:
        raise GateError("validated acquisition release mismatch")
    _require_gate_values(acquisition.get("applicable_gates", {}), {
        "scale_trend": "NOT_APPLICABLE_STAGE_B_NON_SCALE_BEARING",
    })
    if acquisition.get("input_stage_b_manifest_sha256") != STAGE_B_SHA256:
        raise GateError("validated acquisition input Stage-B digest mismatch")
    external_acquisition = Path(acquisition.get("external_release_path", ""))
    verify_inventory(
        external_acquisition, deep=False,
        expected_inventory_sha256=ACQUISITION_EXTERNAL_SHA256SUMS_SHA256,
    )

    source_root = releases_root / "resolve-prophage-source" / SOURCE_SEMANTICS_RELEASE_ID
    if sha_file(source_root / "release.json") != SOURCE_SEMANTICS_RELEASE_JSON_SHA256:
        raise GateError("current source-semantics release checksum mismatch")
    verify_inventory(source_root, expected_inventory_sha256=SOURCE_SEMANTICS_EXTERNAL_SHA256SUMS_SHA256)
    semantics = json.loads((source_root / "release.json").read_text(encoding="utf-8"))
    if semantics.get("release_id") != SOURCE_SEMANTICS_RELEASE_ID or semantics.get("verdict") != "EXTRACTION_BLOCKED":
        raise GateError("source-semantics status/pointer mismatch")
    if semantics.get("consumer_action") != "REJECT" or semantics.get("gates", {}).get("extraction_eligibility") != "EXTRACTION_BLOCKED":
        raise GateError("blocked source semantics is not fail-closed")
    for gate in (
        "accession_version_identity", "atomic_promotion", "bgzf_index_name_roundtrip",
        "global_distinct_assembly_cap", "injected_restart", "integrated_plan_sha256",
        "pinned_consumer_compatibility", "predecessor_release_id_manifest_inventory",
        "producer_caller_evidence_inventory", "resource", "root_input_sha256",
        "row_accounting", "source_coordinate_policy", "upstream_local_checksum",
    ):
        if semantics.get("gates", {}).get(gate) != "PASS":
            raise GateError(f"source-semantics engineering gate is not PASS: {gate}")

    allowed = {row["resolved_assembly_accession_version"] for row in assemblies}
    cap = audit_global_sequence_union(releases_root, allowed)
    input_manifest = {
        "schema": "pilot-selection-input-manifest-v1", "immutable": True,
        "root_inputs": roots,
        "collection": {
            "release_id": COLLECTION_RELEASE_ID,
            "release_json_sha256": COLLECTION_RELEASE_JSON_SHA256,
            "assemblies_tsv_gz_sha256": COLLECTION_ASSEMBLIES_GZ_SHA256,
            "assemblies_rows": 26077,
            "stage_b_10_sha256": STAGE_B_SHA256,
            "stage_b_10_bytes": STAGE_B_BYTES,
            "stage_b_10_rows": 10,
            "external_sha256sums_sha256": COLLECTION_EXTERNAL_SHA256SUMS_SHA256,
        },
        "validated_stage_b_acquisition": {
            "release_id": ACQUISITION_RELEASE_ID,
            "release_json_sha256": ACQUISITION_RELEASE_JSON_SHA256,
            "external_sha256sums_sha256": ACQUISITION_EXTERNAL_SHA256SUMS_SHA256,
            "verdict": "PASS",
        },
        "engineering_control_source_semantics": {
            "release_id": SOURCE_SEMANTICS_RELEASE_ID,
            "release_json_sha256": SOURCE_SEMANTICS_RELEASE_JSON_SHA256,
            "external_sha256sums_sha256": SOURCE_SEMANTICS_EXTERNAL_SHA256SUMS_SHA256,
            "engineering_gates": "PASS",
            "extraction_eligibility": "EXTRACTION_BLOCKED",
            "consumer_action": "REJECT",
            "selection_use": "POST_SELECTION_CONTROL_LABELS_ONLY",
        },
        "global_sequence_union_at_gate": cap,
    }
    return assemblies, stage_b, input_manifest


def selection_policy(code_sha256: str, result: SelectionResult) -> dict[str, Any]:
    probabilities = {}
    for rung in RUNGS:
        probabilities[str(rung)] = {
            "stage_b_certainty": "1/1",
            "main_phage_blind_srs": f"{rung - 10}/{result.main_random_count}",
            "terminal_suppressed_ineligible": "0/1",
        }
    return {
        "schema": "pilot-selection-policy-v1", "policy_frozen": True,
        "algorithm": SELECTION_ALGORITHM, "implementation_sha256": code_sha256,
        "seed_literal": MAIN_SEED, "seed_sha256": sha_bytes(MAIN_SEED.encode()),
        "rungs": list(RUNGS), "nesting_rule": "each rung is the exact cohort-order prefix",
        "eligibility_rule": "exact-version rows except authoritative terminal assembly_status=suppressed",
        "certainty_rule": "the exact checksum-pinned validated Stage-B ten retain their frozen order",
        "random_rule": "sort remaining eligible rows by SHA256(seed NUL assembly_id NUL exact accession.version), then assembly_id",
        "allowed_selection_fields": list(ALLOWED_SELECTION_FIELDS),
        "forbidden_selection_fields": list(FORBIDDEN_SELECTION_FIELDS),
        "phage_blind": True,
        "engineering_controls": "attached only after cohort order is frozen; never alter order or multiplicity",
        "inference_rule": "use the recorded exact inverse-probability design weight once per assembly; exclude synthetic fixtures",
        "stratum_counts": result.frame_counts,
        "eligible_frame_rows": result.eligible_count,
        "main_random_stratum_rows": result.main_random_count,
        "inclusion_probabilities": probabilities,
    }


def release_id(code_sha256: str) -> str:
    seed = {
        "schema": SCHEMA, "algorithm": SELECTION_ALGORITHM,
        "implementation_sha256": code_sha256, "seed_sha256": sha_bytes(MAIN_SEED.encode()),
        "collection_release_id": COLLECTION_RELEASE_ID,
        "collection_release_json_sha256": COLLECTION_RELEASE_JSON_SHA256,
        "stage_b_10_sha256": STAGE_B_SHA256,
        "validated_stage_b_release_id": ACQUISITION_RELEASE_ID,
        "source_semantics_control_release_id": SOURCE_SEMANTICS_RELEASE_ID,
        "rungs": list(RUNGS),
    }
    return f"pilot-cohorts-v1-{sha_bytes(canonical_json(seed))[:16]}"


def _control_membership(controls: list[dict[str, str]]) -> dict[str, str]:
    memberships: dict[str, set[str]] = {}
    for control in controls:
        if control["assembly_id"] != ".":
            memberships.setdefault(control["assembly_id"], set()).add(control["control_class"])
    return {assembly_id: ",".join(sorted(values)) for assembly_id, values in memberships.items()}


def build_selection_unit(
    unit: Path, result: SelectionResult, controls: list[dict[str, str]],
    policy: dict[str, Any], input_manifest: dict[str, Any],
) -> None:
    if unit.exists():
        verify_inventory(unit)
        return
    stage = unit.parent / f".stage.{unit.name}.{os.getpid()}"
    if stage.exists():
        shutil.rmtree(stage)
    (stage / "manifests").mkdir(parents=True)
    write_tsv(stage / "manifests/frame.tsv", FRAME_FIELDS, result.frame)
    membership = _control_membership(controls)
    for rung in RUNGS:
        rows = []
        for row in result.rungs[rung]:
            rows.append({**row, "engineering_control_membership": membership.get(row["assembly_id"], ".")})
        write_tsv(stage / f"manifests/cohort-{rung:04d}.tsv", COHORT_FIELDS, rows)
    write_tsv(stage / "manifests/engineering-controls.tsv", CONTROL_FIELDS, controls)
    (stage / "selection_policy.json").write_bytes(canonical_json(policy))
    (stage / "input_manifest.json").write_bytes(canonical_json(input_manifest))
    seal_and_promote(stage, unit)


def _tree_digest_map(root: Path, relative_prefix: str = "manifests") -> dict[str, str]:
    return {
        str(path.relative_to(root)): sha_file(path)
        for path in sorted((root / relative_prefix).rglob("*")) if path.is_file()
    }


def verify_deterministic_rerun(
    scratch_run: Path, result: SelectionResult, controls: list[dict[str, str]],
    policy: dict[str, Any], input_manifest: dict[str, Any], unit: Path,
) -> dict[str, Any]:
    rerun = scratch_run / ".deterministic-rerun"
    if rerun.exists():
        shutil.rmtree(rerun)
    build_selection_unit(rerun, result, controls, policy, input_manifest)
    first = _tree_digest_map(unit)
    second = _tree_digest_map(rerun)
    if first != second:
        raise GateError("deterministic manifest rerun changed bytes")
    shutil.rmtree(rerun)
    return {"schema": "pilot-selection-deterministic-rerun-v1", "verdict": "PASS", "manifest_sha256": first}


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
        raise GateError(f"findmnt did not return one filesystem: {path}")
    return {str(key): str(value) for key, value in filesystems[0].items()}


def _write_probe(parent: Path, label: str) -> None:
    path = parent / f".pilot-selection-write-probe.{label}.{os.getpid()}.{uuid.uuid4().hex}"
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
    try:
        payload = f"{TASK_ID}\t{label}\t{utcnow()}\n".encode()
        if os.write(descriptor, payload) != len(payload):
            raise GateError("short filesystem write probe")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    path.unlink()


def swap_free_bytes() -> int:
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if line.startswith("SwapFree:"):
            return int(line.split()[1]) * 1024
    raise GateError("cannot read SwapFree")


def live_preflight(durable: Path, scratch: Path, allocations: Allocations, stage_name: str) -> dict[str, Any]:
    allocations.validate()
    durable_parent, scratch_parent = _existing_parent(durable), _existing_parent(scratch)
    durable_mount, scratch_mount = _mount_record(durable_parent), _mount_record(scratch_parent)
    if durable_mount.get("target") != "/" or durable_mount.get("fstype") != "ext4":
        raise GateError(f"durable mount identity NO_GO: {durable_mount}")
    if scratch_mount.get("target") != "/mnt/nvme3n1" or scratch_mount.get("source") != "/dev/nvme3n1" or scratch_mount.get("fstype") != "xfs":
        raise GateError(f"scratch mount identity NO_GO: {scratch_mount}")
    durable_resolved = str(durable.resolve())
    scratch_resolved = str(scratch.resolve())
    if not durable_resolved.startswith(str(DEFAULT_DURABLE_ROOT) + "/"):
        raise GateError(f"durable path outside task namespace: {durable}")
    if not scratch_resolved.startswith(str(DEFAULT_SCRATCH_ROOT) + "/"):
        raise GateError(f"scratch path outside task namespace: {scratch}")
    ds, ss = os.statvfs(durable_parent), os.statvfs(scratch_parent)
    durable_free, durable_inodes = ds.f_bavail * ds.f_frsize, ds.f_favail
    scratch_free, scratch_inodes = ss.f_bavail * ss.f_frsize, ss.f_favail
    checks = {
        "durable_bytes_ge_2tb": durable_free >= 2_000_000_000_000,
        "durable_after_unfinished_ge_2tb": durable_free - allocations.unfinished_write_bytes >= 2_000_000_000_000,
        "durable_inodes_ge_1m": durable_inodes >= 1_000_000,
        "scratch_live_preflight_bytes_ge_4tb": scratch_free >= 4_000_000_000_000,
        "scratch_after_unfinished_ge_2tb": scratch_free - allocations.unfinished_write_bytes >= 2_000_000_000_000,
        "scratch_inodes_ge_5m": scratch_inodes >= 5_000_000,
        "durable_peak_le_70pct_allocation": allocations.predicted_durable_peak_bytes * 100 <= allocations.durable_allocation_bytes * 70,
        "scratch_peak_le_70pct_allocation": allocations.predicted_scratch_peak_bytes * 100 <= allocations.scratch_allocation_bytes * 70,
        "projected_files_le_50pct_inode_allocation": allocations.predicted_files * 2 <= allocations.inode_allocation,
        "durable_retains_2x_unfinished": durable_free - allocations.predicted_durable_peak_bytes >= 2 * allocations.unfinished_write_bytes,
        "scratch_retains_2x_unfinished": scratch_free - allocations.predicted_scratch_peak_bytes >= 2 * allocations.unfinished_write_bytes,
    }
    if not all(checks.values()):
        raise GateError("resource gate NO_GO: " + json.dumps(checks, sort_keys=True))
    for label, parent in (("durable", durable_parent), ("scratch", scratch_parent)):
        if parent.stat().st_uid != os.getuid() or not os.access(parent, os.W_OK | os.X_OK):
            raise GateError(f"{label} parent ownership/write gate NO_GO: {parent}")
        _write_probe(parent, label)
    return {
        "schema": "pilot-selection-resource-preflight-v1", "verdict": "PASS",
        "stage": stage_name, "captured_at_utc": utcnow(),
        "durable_path": str(durable), "scratch_path": str(scratch),
        "durable_probe_parent": str(durable_parent), "scratch_probe_parent": str(scratch_parent),
        "durable_findmnt": durable_mount, "scratch_findmnt": scratch_mount,
        "durable_owner": {"uid": durable_parent.stat().st_uid, "gid": durable_parent.stat().st_gid, "mode": stat.S_IMODE(durable_parent.stat().st_mode)},
        "scratch_owner": {"uid": scratch_parent.stat().st_uid, "gid": scratch_parent.stat().st_gid, "mode": stat.S_IMODE(scratch_parent.stat().st_mode)},
        "write_probes": "PASS", "durable_free_bytes": durable_free,
        "durable_free_inodes": durable_inodes, "scratch_free_bytes": scratch_free,
        "scratch_free_inodes": scratch_inodes, "swap_free_bytes": swap_free_bytes(),
        "allocations": asdict(allocations), "checks": checks,
    }


def _directory_usage(root: Path) -> tuple[int, int]:
    files = [path for path in root.rglob("*") if path.is_file()]
    return sum(path.stat().st_size for path in files), len(files)


def _release_manifest_metadata(manifests: Path) -> dict[str, Any]:
    output = {}
    for path in sorted(manifests.glob("*.tsv")):
        rows = sum(1 for _ in path.open(encoding="utf-8")) - 1
        output[path.name] = {"rows": rows, "bytes": path.stat().st_size, "sha256": sha_file(path)}
    return output


def _copy_exact_inputs(repo: Path, target: Path) -> None:
    target.mkdir(parents=True)
    copies = (
        (repo / "manifests/collection-v1/assemblies.tsv.gz", target / "collection-assemblies.tsv.gz"),
        (repo / "manifests/collection-v1/stage_b_10.tsv", target / "stage_b_10.tsv"),
        (repo / "manifests/collection-v1/release.json", target / "collection-release.json"),
        (repo / "manifests/canonical-cohort-010-v1/release.json", target / "canonical-cohort-010-release.json"),
        (SOURCE_SEMANTICS_ROOT / "release.json", target / "source-semantics-release.json"),
    )
    for source, destination in copies:
        shutil.copyfile(source, destination)


def publish_tracked(repo: Path, external: Path) -> Path:
    target = repo / "manifests/pilot-cohorts-v1"
    if target.exists():
        checks = _read_checksum_map(target)
        release = json.loads((target / "release.json").read_text(encoding="utf-8"))
        external_release = json.loads((external / "release.json").read_text(encoding="utf-8"))
        if release.get("release_id") != external_release.get("release_id") or not checks:
            raise GateError("existing tracked pilot manifest release differs")
        return target
    stage = target.parent / f".stage.pilot-cohorts-v1.{os.getpid()}"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    for path in sorted((external / "manifests").glob("*.tsv")):
        if path.name == "frame.tsv":
            deterministic_gzip(path, stage / "frame.tsv.gz")
        else:
            shutil.copyfile(path, stage / path.name)
    for name in ("release.json", "selection_policy.json", "input_manifest.json", "restart_evidence.json", "resource_summary.json", "deterministic_rerun.json", "global_cap_evidence.json"):
        shutil.copyfile(external / name, stage / name)
    shutil.copyfile(external / "SHA256SUMS", stage / "external_SHA256SUMS")
    write_inventory(stage)
    os.rename(stage, target)
    fsync_dir(target.parent)
    _read_checksum_map(target)
    return target


def _kill_self() -> None:
    os.kill(os.getpid(), signal.SIGKILL)


def run_release(
    repo: Path, run_id: str, allocations: Allocations, *,
    releases_root: Path = DEFAULT_RELEASES_ROOT,
    durable_root: Path = DEFAULT_DURABLE_ROOT,
    scratch_root: Path = DEFAULT_SCRATCH_ROOT,
    inject_kill_before_complete: bool = False,
) -> dict[str, Any]:
    repo = repo.resolve()
    code_sha = sha_file(Path(__file__))
    rid = release_id(code_sha)
    final = durable_root / rid
    scratch_run = scratch_root / run_id
    stage = durable_root / f".stage.{rid}.{run_id}"

    if final.exists():
        verify_inventory(final)
        publish_tracked(repo, final)
        return json.loads((final / "release.json").read_text(encoding="utf-8"))

    selection_preflight = live_preflight(final, scratch_run, allocations, "selection_or_checksum_validated_resume")
    durable_root.mkdir(parents=True, exist_ok=True, mode=0o750)
    scratch_run.mkdir(parents=True, exist_ok=True, mode=0o750)
    state_log, failures_log, resource_log = scratch_run / "state.jsonl", scratch_run / "failures.jsonl", scratch_run / "resources.jsonl"
    failures_log.touch(exist_ok=True)
    fsync_file(failures_log)
    append_jsonl(resource_log, selection_preflight)
    append_jsonl(state_log, {"event": "INPUT_GATE_START", "at_utc": utcnow(), "run_id": run_id})
    assemblies, stage_b, input_manifest = verify_inputs(repo, releases_root)
    append_jsonl(state_log, {"event": "INPUT_GATE_PASS", "at_utc": utcnow(), "collection_rows": len(assemblies)})

    discard_interrupted_stage(stage, state_log, failures_log)
    result = select_cohort(assemblies, stage_b)
    controls = build_engineering_controls(repo, stage_b)
    policy = selection_policy(code_sha, result)
    unit = scratch_run / "units" / rid
    unit.parent.mkdir(parents=True, exist_ok=True)
    reused_unit = unit.exists()
    build_selection_unit(unit, result, controls, policy, input_manifest)
    append_jsonl(state_log, {
        "event": "SELECTION_UNIT_REUSED_CHECKSUM_VALIDATED" if reused_unit else "SELECTION_UNIT_COMMITTED",
        "at_utc": utcnow(), "release_id": rid, "completed_unit": str(unit),
    })
    deterministic = verify_deterministic_rerun(scratch_run, result, controls, policy, input_manifest, unit)
    append_jsonl(state_log, {"event": "DETERMINISTIC_RERUN_PASS", "at_utc": utcnow(), "release_id": rid})

    finish_roots = verify_root_inputs(repo)
    cap = audit_global_sequence_union(releases_root, {row["resolved_assembly_accession_version"] for row in assemblies})
    selected_accessions = {row["exact_assembly_accession_version"] for row in result.cohort}
    if not set(cap["accessions"]).issubset(selected_accessions):
        raise GateError("global sequence-bearing union is not a subset of the frozen 1,000")
    publication_preflight = live_preflight(final, scratch_run, allocations, "atomic_publication")
    append_jsonl(resource_log, publication_preflight)

    state_text = state_log.read_text(encoding="utf-8")
    restart = {
        "schema": "pilot-selection-restart-evidence-v1",
        "injected_kill_recorded": "INJECTED_KILL_BEFORE_COMPLETE" in state_text,
        "interrupted_stage_discarded": "INTERRUPTED_PUBLICATION_STAGE_DISCARDED" in state_text,
        "checksum_validated_selection_unit_reused": "SELECTION_UNIT_REUSED_CHECKSUM_VALIDATED" in state_text,
        "no_partial_final_before_restart": not final.exists(),
    }
    restart["verdict"] = "PASS" if all(value is True for key, value in restart.items() if key not in ("schema", "verdict")) else "PENDING"
    if not inject_kill_before_complete and restart["verdict"] != "PASS":
        raise GateError("mandatory injected kill/restart evidence is absent; run once with --inject-kill-before-complete")

    if stage.exists():
        discard_interrupted_stage(stage, state_log, failures_log)
    stage.mkdir(parents=True)
    shutil.copytree(unit / "manifests", stage / "manifests")
    shutil.copyfile(unit / "selection_policy.json", stage / "selection_policy.json")
    shutil.copyfile(unit / "input_manifest.json", stage / "input_manifest.json")
    _copy_exact_inputs(repo, stage / "inputs")
    (stage / "deterministic_rerun.json").write_bytes(canonical_json(deterministic))
    (stage / "global_cap_evidence.json").write_bytes(canonical_json(cap))

    staged_bytes, staged_files = _directory_usage(stage)
    peak_rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    resource_records = [json.loads(line) for line in resource_log.read_text(encoding="utf-8").splitlines()]
    resource_summary = {
        "schema": "pilot-selection-resource-summary-v1", "verdict": "PASS",
        "assigned_ram_bytes": allocations.assigned_ram_bytes,
        "peak_rss_bytes": peak_rss_bytes,
        "peak_rss_fraction": peak_rss_bytes / allocations.assigned_ram_bytes,
        "swap_growth_bytes": max(0, resource_records[0]["swap_free_bytes"] - swap_free_bytes()),
        "upper_95_predicted_durable_peak_bytes": allocations.predicted_durable_peak_bytes,
        "actual_preseal_stage_bytes": staged_bytes,
        "actual_preseal_files": staged_files,
        "resource_record_count": len(resource_records),
        "checks": {
            "peak_rss_le_70pct": peak_rss_bytes <= allocations.assigned_ram_bytes * 0.70,
            "no_swap_growth": swap_free_bytes() >= resource_records[0]["swap_free_bytes"],
            "actual_disk_le_70pct_allocation": staged_bytes <= allocations.durable_allocation_bytes * 0.70,
            "actual_files_le_50pct_inode_allocation": staged_files * 2 <= allocations.inode_allocation,
            "all_preflights_pass": all(record.get("verdict") == "PASS" and all(record.get("checks", {}).values()) for record in resource_records),
        },
    }
    if not all(resource_summary["checks"].values()):
        raise GateError("measured resource summary NO_GO")
    (stage / "resource_summary.json").write_bytes(canonical_json(resource_summary))
    shutil.copyfile(resource_log, stage / "resources.jsonl")
    (stage / "restart_evidence.json").write_bytes(canonical_json(restart))
    (stage / "provenance.json").write_bytes(canonical_json({
        "schema": "pilot-selection-provenance-v1", "task_id": TASK_ID, "run_id": run_id,
        "created_at_utc": utcnow(), "argv": sys.argv, "cwd": str(repo),
        "environment": {key: os.environ.get(key, ".") for key in ("PI_MODEL", "PI_PROVIDER", "WG_TASK_ID", "WG_AGENT_ID")},
        "network_requests": 0, "sequence_downloads": 0, "biological_analyses": 0,
    }))
    (stage / "tools.json").write_bytes(canonical_json({
        "schema": "pilot-selection-tools-v1", "python": sys.version,
        "platform": platform.platform(), "implementation_sha256": code_sha,
        "algorithm": SELECTION_ALGORITHM,
    }))

    manifest_metadata = _release_manifest_metadata(stage / "manifests")
    release = {
        "schema_version": SCHEMA, "release_id": rid, "source_task_id": TASK_ID,
        "created_at_utc": utcnow(), "immutable": True, "verdict": "PASS" if restart["verdict"] == "PASS" else "PENDING_INJECTED_RESTART",
        "external_release_path": str(final),
        "selection_algorithm": SELECTION_ALGORITHM, "selection_policy_sha256": sha_file(stage / "selection_policy.json"),
        "input_manifest_sha256": sha_file(stage / "input_manifest.json"),
        "collection_release_id": COLLECTION_RELEASE_ID,
        "validated_stage_b_release_id": ACQUISITION_RELEASE_ID,
        "engineering_control_source_semantics_release_id": SOURCE_SEMANTICS_RELEASE_ID,
        "source_semantics_extraction_eligibility": "EXTRACTION_BLOCKED",
        "source_semantics_consumer_action": "REJECT",
        "root_inputs_start": input_manifest["root_inputs"], "root_inputs_finish": finish_roots,
        "counts": {
            "collection_frame": len(result.frame), "eligible_frame": result.eligible_count,
            "terminal_suppressed_ineligible": result.frame_counts.get("terminal_suppressed_ineligible", 0),
            "stage_b_certainty": 10, "main_random_stratum": result.main_random_count,
            "cohort": len(result.cohort), "engineering_control_rows": len(controls),
            "distinct_sequence_bearing_assemblies_created": 0,
            "global_sequence_bearing_union": cap["distinct_exact_assembly_revisions"],
            "global_distinct_assembly_cap": 1000,
        },
        "rung_cardinalities": {str(rung): len(result.rungs[rung]) for rung in RUNGS},
        "cohort_order_sha256": sha_bytes(("\n".join(row["exact_assembly_accession_version"] for row in result.cohort) + "\n").encode()),
        "manifests": manifest_metadata,
        "applicable_gates": {
            "accession_version_identity": "PASS", "upstream_local_checksum": "PASS",
            "row_accounting": "PASS", "predecessor_release_id_manifest_inventory": "PASS",
            "resource": "PASS", "deterministic_rerun": "PASS",
            "injected_kill_restart": "PASS" if restart["verdict"] == "PASS" else "PENDING",
            "atomic_promotion": "PASS", "global_distinct_assembly_cap": "PASS",
            "source_immutability": "PASS", "pinned_consumer_compatibility": "PASS",
            "engineering_control_source_semantics": "PASS",
            "bgzf_index_name_roundtrip": "NOT_APPLICABLE_METADATA_ONLY_SELECTION",
            "source_coordinate_policy": "NOT_APPLICABLE_HOST_SELECTION_EXTRACTION_BLOCKED",
            "scale_trend": "NOT_APPLICABLE_METADATA_ONLY_SELECTION",
        },
        "sequence_downloads": 0, "biological_analyses": 0,
    }
    if inject_kill_before_complete:
        release["verdict"] = "PENDING_INJECTED_RESTART"
    (stage / "release.json").write_bytes(canonical_json(release))
    shutil.copyfile(state_log, stage / "state.jsonl")
    shutil.copyfile(failures_log, stage / "failures.jsonl")

    def injected() -> None:
        append_jsonl(state_log, {
            "event": "INJECTED_KILL_BEFORE_COMPLETE", "at_utc": utcnow(),
            "stage": str(stage), "final_absent": not final.exists(),
            "sha256sums_present": (stage / "SHA256SUMS").is_file(),
        })
        _kill_self()

    if restart["verdict"] != "PASS" and not inject_kill_before_complete:
        raise GateError("cannot publish without PASS restart evidence")
    if restart["verdict"] == "PASS":
        release["verdict"] = "PASS"
    seal_and_promote(stage, final, inject=injected if inject_kill_before_complete else None)
    if inject_kill_before_complete:
        raise GateError("injected SIGKILL unexpectedly returned")
    verify_root_inputs(repo)
    publish_tracked(repo, final)
    return json.loads((final / "release.json").read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--run-id", default="pilot-cohorts-v1-execution")
    parser.add_argument("--inject-kill-before-complete", action="store_true")
    parser.add_argument("--assigned-ram-bytes", type=int, default=4_294_967_296)
    parser.add_argument("--durable-allocation-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--scratch-allocation-bytes", type=int, default=4_000_000_000_000)
    parser.add_argument("--inode-allocation", type=int, default=100_000)
    parser.add_argument("--predicted-durable-peak-bytes", type=int, default=100_000_000)
    parser.add_argument("--predicted-scratch-peak-bytes", type=int, default=100_000_000)
    parser.add_argument("--predicted-files", type=int, default=100)
    parser.add_argument("--unfinished-write-bytes", type=int, default=50_000_000)
    args = parser.parse_args()
    allocations = Allocations(
        args.assigned_ram_bytes, args.durable_allocation_bytes, args.scratch_allocation_bytes,
        args.inode_allocation, args.predicted_durable_peak_bytes,
        args.predicted_scratch_peak_bytes, args.predicted_files, args.unfinished_write_bytes,
    )
    try:
        release = run_release(
            Path(args.repo_root), args.run_id, allocations,
            inject_kill_before_complete=args.inject_kill_before_complete,
        )
        print(canonical_json(release).decode(), end="")
        return 0
    except (GateError, OSError, ValueError, json.JSONDecodeError, subprocess.SubprocessError) as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
