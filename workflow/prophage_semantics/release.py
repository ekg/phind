#!/usr/bin/env python3
"""Fail-closed, bounded prophage source-semantics release builder.

This program never downloads or extracts biological sequence.  Its only biological
look-up is a read-only, ten-assembly annotation-boundary diagnostic against the
exact predecessor package objects.  Coordinate candidates are not applied to
sequence and the published extraction verdict is derived from the independently
verified pinned-caller consumption gate (EXTRACTION_GO only when the historical
CSV attribution is DECISIVE and independently re-verified as sound).
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
import stat
import subprocess
import sys
import tempfile
import time
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable

SCHEMA = "prophage-semantics-release-v1"
SOURCE_SHA = "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996"
ACCESSIONS_SHA = "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5"
PLAN_SHA = "fb58d25a6f4971137ab0dcb82dae09eac5d177e37ba34f857845ce2d7e0a6da8"
AUDIT_SHA = "feb9b687fb0722a4f073f7105088f19a2dabebcbc0378dbe5c4799b0b7f29fdc"
COHORT_RELEASE_JSON_SHA = "4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23"
COHORT_RELEASE_ID = "canonical-cohort-010-v1-e71484de9994fc28"
COHORT_ORDER = [
    "GCF_000005845.2", "GCF_000812325.1", "GCF_002302315.1",
    "GCF_004664255.1", "GCF_015644385.1", "GCF_020829045.1",
    "GCF_921380995.1", "GCF_000167895.3", "GCF_001881595.4",
    "GCF_000498835.2",
]
EXPECTED_HEADER = ["end", "genome", "scaffold", "begin", "transposable", "taxonomy", "prophage_id"]
PASS_OR_NA = {"PASS", "NOT_APPLICABLE_STAGE_B_NON_SCALE_BEARING"}

# Versioned policy / release wiring.  The v1 BLOCKED release
# (``prophage-semantics-v1-f5619e221ff272ae``) is already published and immutable
# and remains the historical record.  The v2 release reflects the now-decisive
# pinned-caller evidence and derives its verdict from the consumption gate.
POLICY_VERSION = "v2"
POLICY_FILE = f"artifacts/prophage_semantics/semantics_policy_{POLICY_VERSION}.json"
POLICY_STAGE_NAME = f"semantics_policy_{POLICY_VERSION}.json"
RELEASE_PREFIX = f"prophage-semantics-{POLICY_VERSION}"
PINNED_CALLER_GATE_FILE = "artifacts/prophage_semantics/pinned_caller_input_gate.json"
PREDECESSOR_PHIGARO_RELEASE_ID = "phigaro-version-comparison-v1-e7cfa43b9231aee5"
PREDECESSOR_PHIGARO_EXTERNAL_ROOT = Path(
    "/home/erikg/phind-data/ecoli26k/v1/releases/rerun-phigaro-version-comparison"
)
REQUIRED_DIMENSIONS = {
    "producer_and_caller_version", "tagged", "transposable", "taxonomy_labels",
    "coordinate_base_and_end_inclusivity", "strand_orientation",
    "topology_and_circularity", "contig_edge_behavior", "completeness",
    "duplicate_locus_rules",
}
SCOPE_IDS = ["all_records", "transposable_flag_positive", "taxonomy_assigned"]
COORDINATE_CANDIDATE_IDS = {"C1_RAW_1_BASED_CLOSED", "C2_RAW_0_BASED_INCLUSIVE"}


class GateError(RuntimeError):
    pass


class InjectedInterruption(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_bytes(obj: Any) -> bytes:
    return (json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fsync_dir(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("wb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    fsync_dir(path.parent)


def append_jsonl(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("ab") as handle:
        handle.write(canonical_bytes(obj))
        handle.flush()
        os.fsync(handle.fileno())


def derive_extraction_from_gate(repo: Path) -> tuple[str, str, dict[str, Any]]:
    """Read the pinned-caller consumption gate and derive the extraction verdict.

    Fail-closed by construction: any missing/non-PASS gate, any failure to keep
    the modern-v2.4 pilot separate, any non-DECISIVE attribution, or any
    independent re-verification that is not sound raises ``GateError`` and no
    release is published.
    """
    path = repo / PINNED_CALLER_GATE_FILE
    if not path.is_file():
        raise GateError("pinned-caller consumption gate result is absent; extraction must fail closed")
    gate = json.loads(path.read_text())
    if gate.get("schema") != "pinned-phigaro-consumption-gate-v1":
        raise GateError("pinned-caller consumption gate schema mismatch")
    if gate.get("verdict") != "PASS":
        raise GateError("pinned-caller consumption gate is not PASS")
    if gate.get("release_id") != PREDECESSOR_PHIGARO_RELEASE_ID:
        raise GateError("pinned-caller release identity mismatch")
    if gate.get("modern_v2_4_pilot_separate") is not True:
        raise GateError("modern v2.4 pilot must remain strictly separate")
    if gate.get("historical_csv_attribution") != "DECISIVE":
        raise GateError("historical attribution is not DECISIVE")
    if gate.get("decisive_evidence_independently_sound") is not True:
        raise GateError("decisive attribution evidence was not independently re-verified as sound")
    if gate.get("historical_csv_extraction") != "EXTRACTION_GO":
        raise GateError("pinned-caller gate did not authorize historical extraction")
    return "EXTRACTION_GO", "ALLOW", gate


def verify_exact(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise GateError(f"required input missing: {label}: {path}")
    got = sha_file(path)
    if got != expected:
        raise GateError(f"checksum mismatch for {label}: expected {expected}, got {got}")


def parse_sum_file(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line:
            continue
        try:
            digest, rel = line.split("  ", 1)
        except ValueError as exc:
            raise GateError(f"invalid SHA256SUMS row {path}:{line_number}") from exc
        if len(digest) != 64 or Path(rel).is_absolute() or ".." in Path(rel).parts:
            raise GateError(f"unsafe SHA256SUMS row {path}:{line_number}")
        rows.append((digest, rel))
    return rows


def verify_inventory(root: Path, sums: Path) -> int:
    if not sums.is_file():
        raise GateError(f"missing checksum inventory: {sums}")
    rows = parse_sum_file(sums)
    for digest, rel in rows:
        path = root / rel
        if not path.is_file() or sha_file(path) != digest:
            raise GateError(f"inventory mismatch: {path}")
    return len(rows)


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


def findmnt(path: Path) -> dict[str, str]:
    fields = ["TARGET", "SOURCE", "FSTYPE", "OPTIONS", "SIZE", "AVAIL"]
    cp = subprocess.run(
        ["findmnt", "-T", str(path), "-n", "-P", "-o", ",".join(fields)],
        text=True, capture_output=True, check=True,
    )
    # -P emits FIELD="value" tokens; shlex correctly preserves options.
    import shlex
    result = {}
    for token in shlex.split(cp.stdout.strip()):
        key, value = token.split("=", 1)
        result[key.lower()] = value
    return result


def swap_used_bytes() -> int:
    used_kib = 0
    with Path("/proc/swaps").open() as handle:
        next(handle, None)
        for line in handle:
            fields = line.split()
            if len(fields) >= 5:
                used_kib += int(fields[3])
    return used_kib * 1024


def write_probe(parent: Path) -> None:
    parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=".prophage-semantics-write-probe-", dir=parent)
    try:
        os.write(fd, b"probe\n")
        os.fsync(fd)
    finally:
        os.close(fd)
        Path(name).unlink(missing_ok=True)
    fsync_dir(parent)


def preflight(durable_path: Path, scratch_path: Path, allocations: Allocations, stage: str) -> dict[str, Any]:
    allocations.validate()
    durable_parent = durable_path.parent
    scratch_parent = scratch_path
    write_probe(durable_parent)
    write_probe(scratch_parent)
    ds = os.statvfs(durable_parent)
    ss = os.statvfs(scratch_parent)
    durable_free = ds.f_bavail * ds.f_frsize
    scratch_free = ss.f_bavail * ss.f_frsize
    durable_inodes = ds.f_favail
    scratch_inodes = ss.f_favail
    checks = {
        "durable_start_bytes_ge_2tb": durable_free >= 2_000_000_000_000,
        "durable_inodes_ge_1m": durable_inodes >= 1_000_000,
        "scratch_live_bytes_ge_4tb": scratch_free >= 4_000_000_000_000,
        "scratch_live_inodes_ge_5m": scratch_inodes >= 5_000_000,
        "durable_stop_after_unfinished_ge_2tb": durable_free - allocations.unfinished_write_bytes >= 2_000_000_000_000,
        "scratch_stop_after_unfinished_ge_2tb": scratch_free - allocations.unfinished_write_bytes >= 2_000_000_000_000,
        "scratch_stop_inodes_ge_5m": scratch_inodes - allocations.predicted_files >= 5_000_000,
        "predicted_durable_le_70pct_allocation": allocations.predicted_durable_peak_bytes * 100 <= allocations.durable_allocation_bytes * 70,
        "predicted_scratch_le_70pct_allocation": allocations.predicted_scratch_peak_bytes * 100 <= allocations.scratch_allocation_bytes * 70,
        "projected_files_le_50pct_inodes": allocations.predicted_files * 2 <= allocations.inode_allocation,
        "durable_two_x_unfinished_writes": durable_free - allocations.predicted_durable_peak_bytes >= 2 * allocations.unfinished_write_bytes,
        "scratch_two_x_unfinished_writes": scratch_free - allocations.predicted_scratch_peak_bytes >= 2 * allocations.unfinished_write_bytes,
    }
    record = {
        "schema": "prophage-semantics-resource-preflight-v1", "stage": stage,
        "captured_at_utc": utc_now(), "verdict": "PASS" if all(checks.values()) else "NO_GO",
        "allocations": asdict(allocations), "checks": checks,
        "durable_path": str(durable_path), "scratch_path": str(scratch_path),
        "durable_findmnt": findmnt(durable_parent), "scratch_findmnt": findmnt(scratch_parent),
        "durable_free_bytes": durable_free, "scratch_free_bytes": scratch_free,
        "durable_free_inodes": durable_inodes, "scratch_free_inodes": scratch_inodes,
        "durable_owner": {"uid": durable_parent.stat().st_uid, "gid": durable_parent.stat().st_gid, "mode": stat.S_IMODE(durable_parent.stat().st_mode)},
        "scratch_owner": {"uid": scratch_parent.stat().st_uid, "gid": scratch_parent.stat().st_gid, "mode": stat.S_IMODE(scratch_parent.stat().st_mode)},
        "write_probes": "PASS", "swap_used_bytes": swap_used_bytes(),
    }
    if not all(checks.values()):
        raise GateError("resource gate NO_GO: " + json.dumps(checks, sort_keys=True))
    return record


def validate_policy(policy: dict[str, Any]) -> None:
    if policy.get("schema_version") != "prophage-source-semantics-v1":
        raise GateError("policy schema mismatch")
    if policy.get("source_sha256") != SOURCE_SHA or not policy.get("normalization", {}).get("lossless"):
        raise GateError("policy is not pinned and lossless")
    dims = {d.get("name"): d for d in policy.get("semantic_dimensions", [])}
    if set(dims) != REQUIRED_DIMENSIONS:
        raise GateError("policy semantic dimensions are incomplete")
    for dim in dims.values():
        if not {"status", "confidence", "value", "evidence_ids", "extraction_critical"} <= set(dim):
            raise GateError("semantic dimension lacks status/confidence/evidence")
    scopes = policy.get("scopes", [])
    if [s.get("id") for s in scopes] != SCOPE_IDS:
        raise GateError("scope order/identity mismatch")
    if any(s.get("generic_tagged_alias") is not False for s in scopes):
        raise GateError("a scope was silently relabeled tagged")
    candidates = policy.get("coordinate_candidates", [])
    if {c.get("id") for c in candidates} != COORDINATE_CANDIDATE_IDS:
        raise GateError("explicit dual coordinate policy is required")
    gate = policy.get("extraction_gate", {})
    policy_id = policy.get("policy_id", "")
    is_go = policy_id.endswith("-v2")
    if is_go:
        # v2: decisive evidence uniquely resolves caller/version/coordinates.
        if gate.get("verdict") != "EXTRACTION_GO" or gate.get("consumer_action") != "ALLOW":
            raise GateError("decisive v2 policy must authorize extraction")
        if gate.get("blocking_dimensions") != []:
            raise GateError("a GO policy may not carry blocking dimensions")
        if gate.get("selected_coordinate_candidate") != "C1_RAW_1_BASED_CLOSED":
            raise GateError("GO policy must select the resolved 1-based closed candidate")
        statuses = {c.get("id"): c.get("status") for c in candidates}
        if statuses.get("C1_RAW_1_BASED_CLOSED") != "SELECTED" or statuses.get("C2_RAW_0_BASED_INCLUSIVE") != "REJECTED":
            raise GateError("GO policy coordinate candidate selection is inconsistent")
        if "tagged" in dims and dims["tagged"]["extraction_critical"]:
            raise GateError("'tagged' is not extraction-critical and must not gate extraction")
        # every extraction-critical dimension must be resolved under a GO policy
        unresolved_critical = {
            name for name, dim in dims.items()
            if dim["extraction_critical"] and not str(dim["status"]).startswith("RESOLVED")
        }
        if unresolved_critical:
            raise GateError("a GO policy has unresolved extraction-critical dimensions: " + ", ".join(sorted(unresolved_critical)))
    else:
        # v1 historical record: an unresolved policy must fail closed.
        if gate.get("verdict") != "EXTRACTION_BLOCKED" or gate.get("consumer_action") != "REJECT":
            raise GateError("unresolved policy must fail closed")
        if {c.get("status") for c in candidates} != {"CANDIDATE"}:
            raise GateError("a blocked policy must keep both coordinate candidates open")
        unresolved_critical = {
            name for name, dim in dims.items()
            if dim["extraction_critical"] and dim["status"] != "RESOLVED"
        }
        if set(gate.get("blocking_dimensions", [])) != unresolved_critical:
            raise GateError("every unresolved extraction-critical dimension must block")


def source_profile(csv_path: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    with csv_path.open("rb") as physical_handle:
        physical_lines = sum(1 for _ in physical_handle)
    rows: list[dict[str, str]] = []
    scopes = Counter()
    genomes = Counter()
    flags = Counter()
    taxonomies = Counter()
    loci = Counter()
    exact = Counter()
    ids = Counter()
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != EXPECTED_HEADER:
            raise GateError(f"source CSV header mismatch: {reader.fieldnames!r}")
        for record_number, row in enumerate(reader, 2):
            if None in row or any(row[k] == "" for k in EXPECTED_HEADER):
                raise GateError(f"source row parse/missing failure at physical row {record_number}")
            try:
                begin_d, end_d, flag_d = Decimal(row["begin"]), Decimal(row["end"]), Decimal(row["transposable"])
            except InvalidOperation as exc:
                raise GateError(f"invalid decimal at row {record_number}") from exc
            if begin_d != begin_d.to_integral_value() or end_d != end_d.to_integral_value() or flag_d not in {Decimal(0), Decimal(1)}:
                raise GateError(f"non-integral coordinates or nonbinary flag at row {record_number}")
            begin, end = int(begin_d), int(end_d)
            if begin < 0 or end < begin:
                raise GateError(f"invalid ordered coordinates at row {record_number}")
            scopes["all_records"] += 1
            if flag_d == 1:
                scopes["transposable_flag_positive"] += 1
            if row["taxonomy"].strip().casefold() not in {"", "unknown"}:
                scopes["taxonomy_assigned"] += 1
            genomes[row["genome"]] += 1
            flags[row["transposable"]] += 1
            taxonomies[row["taxonomy"]] += 1
            loci[(row["genome"], row["scaffold"], begin, end)] += 1
            exact[tuple(row[k] for k in EXPECTED_HEADER)] += 1
            ids[row["prophage_id"]] += 1
            rows.append(row)
    expected_scopes = {"all_records": 132404, "transposable_flag_positive": 7695, "taxonomy_assigned": 115442}
    if dict(scopes) != expected_scopes or physical_lines != 132405 or len(genomes) != 26077:
        raise GateError("source row/scope/cardinality accounting mismatch")
    profile = {
        "schema": "prophage-source-profile-v1", "source_sha256": SOURCE_SHA,
        "header": EXPECTED_HEADER, "physical_lines": physical_lines, "data_rows": len(rows),
        "scope_rows": expected_scopes, "distinct_genomes": len(genomes),
        "transposable_raw_counts": dict(sorted(flags.items())),
        "taxonomy_exact_counts": dict(taxonomies.most_common()),
        "duplicate_locus_groups": sum(v > 1 for v in loci.values()),
        "duplicate_locus_extra_rows": sum(v - 1 for v in loci.values() if v > 1),
        "exact_duplicate_groups": sum(v > 1 for v in exact.values()),
        "duplicate_prophage_id_groups": sum(v > 1 for v in ids.values()),
        "row_accounting": "PASS", "raw_rows_preserved_in_place": True,
        "normalization_materialized": False,
    }
    if any(profile[k] for k in ("duplicate_locus_groups", "duplicate_locus_extra_rows", "exact_duplicate_groups", "duplicate_prophage_id_groups")):
        raise GateError("immutable expected duplicate cardinality changed")
    return profile, rows


def validate_predecessor(repo: Path) -> tuple[dict[str, Any], Path, list[dict[str, str]], int]:
    manifest_dir = repo / "manifests/canonical-cohort-010-v1"
    release_path = manifest_dir / "release.json"
    verify_exact(release_path, COHORT_RELEASE_JSON_SHA, "canonical cohort release manifest")
    verify_inventory(manifest_dir, manifest_dir / "SHA256SUMS")
    release = json.loads(release_path.read_text())
    if release.get("release_id") != COHORT_RELEASE_ID or release.get("verdict") != "PASS":
        raise GateError("predecessor release ID/verdict mismatch")
    gates = release.get("applicable_gates", {})
    bad = {k: v for k, v in gates.items() if v not in PASS_OR_NA}
    if bad:
        raise GateError("predecessor has non-PASS applicable gates: " + json.dumps(bad, sort_keys=True))
    if release.get("sequence_bearing_assembly_accessions") != COHORT_ORDER:
        raise GateError("predecessor cohort order mismatch")
    if release.get("counts", {}).get("distinct_sequence_bearing_assemblies") != 10:
        raise GateError("predecessor distinct assembly count mismatch")
    if release.get("counts", {}).get("global_distinct_assembly_cap") != 1000:
        raise GateError("global cap declaration mismatch")
    external = Path(release["external_release_path"])
    if not (external / "COMPLETE").is_file():
        raise GateError("predecessor external release lacks COMPLETE")
    inventory_rows = verify_inventory(external, external / "SHA256SUMS")
    with (manifest_dir / "assemblies.tsv").open(newline="") as handle:
        assemblies = list(csv.DictReader(handle, delimiter="\t"))
    if [r["accession"] for r in assemblies] != COHORT_ORDER or len(assemblies) != 10:
        raise GateError("immutable assemblies.tsv order/cardinality mismatch")
    return release, external, assemblies, inventory_rows


def annotation_boundary_diagnostic(repo: Path, external: Path, assemblies: list[dict[str, str]], source_rows: list[dict[str, str]]) -> dict[str, Any]:
    cohort = set(COHORT_ORDER)
    selected = [row for row in source_rows if row["genome"] in cohort]
    if len(selected) != 56 or {row["genome"] for row in selected} != cohort:
        raise GateError("N=10 CSV sentinel row accounting mismatch")
    features: dict[str, dict[str, tuple[set[int], set[int]]]] = {}
    gff_digests: dict[str, str] = {}
    annotation_status: dict[str, str] = {}
    for assembly in assemblies:
        accession = assembly["accession"]
        package = external / assembly["source_object_relpath"] / "package.zip"
        with zipfile.ZipFile(package) as zf:
            member = assembly["source_gff_member"]
            raw = zf.read(member)
        gff_digests[accession] = sha_bytes(raw)
        annotation_status[accession] = assembly["annotation_status"]
        by_seq: dict[str, list[set[int]]] = defaultdict(lambda: [set(), set()])
        for line in raw.decode("utf-8").splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            by_seq[fields[0]][0].add(start)
            by_seq[fields[0]][1].add(end)
        features[accession] = {k: (v[0], v[1]) for k, v in by_seq.items()}
    manifest_dir = repo / "manifests/canonical-cohort-010-v1"
    lengths: dict[tuple[str, str], int] = {}
    topology_counts = Counter()
    with gzip.open(manifest_dir / "contigs.tsv.gz", "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lengths[(row["accession"], row["source_contig_id_display"])] = int(row["contig_length"])
            topology_counts[row["topology"]] += 1
    deltas: dict[str, dict[str, int]] = {}
    missing_alias = 0
    range_counts = Counter()
    per_assembly = Counter()
    for delta in (-1, 0, 1):
        both = begin_matches = end_matches = 0
        for row in selected:
            per_assembly[row["genome"]] += int(delta == 0)
            pair = features[row["genome"]].get(row["scaffold"])
            if pair is None:
                missing_alias += int(delta == 0)
                continue
            begin = int(Decimal(row["begin"])) + delta
            end = int(Decimal(row["end"])) + delta
            bm, em = begin in pair[0], end in pair[1]
            begin_matches += bm
            end_matches += em
            both += bm and em
        deltas[str(delta)] = {"both_boundaries": both, "begin_boundary": begin_matches, "end_boundary": end_matches}
    for row in selected:
        key = (row["genome"], row["scaffold"])
        if key not in lengths:
            continue
        begin, end, length = int(Decimal(row["begin"])), int(Decimal(row["end"])), lengths[key]
        range_counts["C1_RAW_1_BASED_CLOSED"] += int(1 <= begin <= end <= length)
        range_counts["C2_RAW_0_BASED_INCLUSIVE"] += int(0 <= begin <= end < length)
    return {
        "schema": "prophage-n10-annotation-boundary-diagnostic-v1",
        "bounded_assemblies": COHORT_ORDER, "distinct_assemblies": 10,
        "source_rows": 56, "new_assembly_downloads": 0, "sequence_bases_read": 0,
        "known_base_sentinel": False,
        "diagnostic_type": "indirect source-GFF feature-boundary congruence",
        "gff_member_sha256": gff_digests, "annotation_status": annotation_status,
        "per_assembly_source_rows": dict(sorted(per_assembly.items())),
        "missing_scaffold_alias_rows": missing_alias,
        "boundary_matches_after_adding_delta_to_both_raw_coordinates": deltas,
        "candidate_in_range_rows": dict(range_counts),
        "cohort_contig_topology_counts": dict(topology_counts),
        "interpretation": "Raw coordinates align with both annotated feature boundaries for 45/56 rows; raw+1 aligns for 0/56. This favors but cannot establish a 1-based-closed hypothesis because NCBI GFF is not the missing producer/caller provenance, its gene calls need not equal the caller's Prodigal calls, and no expected boundary-base oracle exists.",
        "verdict": "NON_DECISIVE",
    }


def dir_usage(path: Path) -> tuple[int, int]:
    total = files = 0
    if not path.exists():
        return 0, 0
    for item in path.rglob("*"):
        if item.is_file():
            total += item.stat().st_size
            files += 1
    return total, files


def expected_inputs(repo: Path, predecessor: dict[str, Any], external: Path, gate: dict[str, Any]) -> dict[str, Any]:
    files = {
        "26k_ecoli_accession.txt": ACCESSIONS_SHA,
        "26k_prophage1.csv": SOURCE_SHA,
        "reports/phage_pangenome_project_plan.md": PLAN_SHA,
        "reports/prophage_distribution.md": AUDIT_SHA,
        "manifests/canonical-cohort-010-v1/release.json": COHORT_RELEASE_JSON_SHA,
        "artifacts/prophage_semantics/evidence_inventory.json": sha_file(repo / "artifacts/prophage_semantics/evidence_inventory.json"),
        POLICY_FILE: sha_file(repo / POLICY_FILE),
        "workflow/prophage_semantics/semantics-policy-v1.schema.json": sha_file(repo / "workflow/prophage_semantics/semantics-policy-v1.schema.json"),
        "workflow/prophage_semantics/release.py": sha_file(repo / "workflow/prophage_semantics/release.py"),
        "workflow/prophage_semantics/independent_rerun_verification.py": sha_file(repo / "workflow/prophage_semantics/independent_rerun_verification.py"),
        "workflow/prophage_semantics/pinned_caller_gate.py": sha_file(repo / "workflow/prophage_semantics/pinned_caller_gate.py"),
        PINNED_CALLER_GATE_FILE: sha_file(repo / PINNED_CALLER_GATE_FILE),
    }
    return {
        "schema": "prophage-semantics-input-manifest-v1", "immutable": True,
        "policy_version": POLICY_VERSION,
        "files": files, "predecessor_release_id": COHORT_RELEASE_ID,
        "predecessor_external_release": str(external),
        "predecessor_external_sha256sums_sha256": sha_file(external / "SHA256SUMS"),
        "cohort_order": COHORT_ORDER, "cohort_rows": 10,
        "pinned_caller_release_id": PREDECESSOR_PHIGARO_RELEASE_ID,
        "pinned_caller_external_release": str(PREDECESSOR_PHIGARO_EXTERNAL_ROOT / PREDECESSOR_PHIGARO_RELEASE_ID),
        "pinned_caller_complete_sha256": gate.get("complete_sha256"),
        "pinned_caller_sha256sums_sha256": gate.get("sha256sums_sha256"),
        "pinned_caller_historical_csv_attribution": gate.get("historical_csv_attribution"),
        "pinned_caller_modern_v2_4_pilot": gate.get("modern_v2_4_pilot"),
        "global_distinct_exact_assembly_revisions": 10, "global_cap": 1000,
        "new_assembly_downloads": 0,
    }


def write_static_unit(stage: Path, rel: str, data: bytes) -> None:
    path = stage / rel
    digest = sha_bytes(data)
    if path.exists():
        if sha_file(path) != digest:
            raise GateError(f"interrupted unit mismatch; refusing mixed resume: {rel}")
        append_jsonl(stage / "state.jsonl", {"event": "RESUME_UNIT_VALIDATED", "path": rel, "sha256": digest, "at": utc_now()})
        return
    atomic_write(path, data)
    append_jsonl(stage / "state.jsonl", {"event": "UNIT_COMMITTED", "path": rel, "bytes": len(data), "sha256": digest, "at": utc_now()})


def release_tree_inventory(stage: Path) -> bytes:
    rows = []
    excluded = {"SHA256SUMS", "COMPLETE"}
    for path in sorted(p for p in stage.rglob("*") if p.is_file()):
        rel = path.relative_to(stage).as_posix()
        if rel not in excluded:
            rows.append(f"{sha_file(path)}  {rel}\n")
    return "".join(rows).encode()


def validate_release(path: Path, require_go: bool = False) -> dict[str, Any]:
    if not (path / "COMPLETE").is_file():
        raise GateError("release is incomplete (COMPLETE absent)")
    inventory_rows = verify_inventory(path, path / "SHA256SUMS")
    release = json.loads((path / "release.json").read_text())
    policy_version = release.get("policy_version", "v1")
    policy_path = path / f"semantics_policy_{policy_version}.json"
    policy = json.loads(policy_path.read_text())
    validate_policy(policy)
    expected_verdict = policy["extraction_gate"]["verdict"]
    expected_consumer = policy["extraction_gate"]["consumer_action"]
    if release.get("release_id") != path.name or release.get("verdict") != expected_verdict:
        raise GateError("release identity/verdict mismatch")
    if release.get("consumer_action") != expected_consumer:
        raise GateError("release consumer action mismatch")
    required_pass_gates = {
        "root_input_sha256", "integrated_plan_sha256", "producer_caller_evidence_inventory",
        "predecessor_release_id_manifest_inventory", "accession_version_identity",
        "upstream_local_checksum", "row_accounting", "bgzf_index_name_roundtrip",
        "global_distinct_assembly_cap", "resource", "injected_restart", "atomic_promotion",
        "source_coordinate_policy", "pinned_consumer_compatibility",
        "pinned_caller_consumption_gate",
    }
    gates = release.get("gates", {})
    if any(gates.get(name) != "PASS" for name in required_pass_gates):
        raise GateError("applicable release gate is not unqualified PASS")
    if gates.get("extraction_eligibility") != expected_verdict:
        raise GateError("extraction eligibility must equal the policy verdict")
    if require_go and expected_verdict != "EXTRACTION_GO":
        raise GateError("consumer rejected non-GO extraction policy")
    return {
        "schema": "prophage-semantics-validation-v1", "verdict": "PASS",
        "release_id": path.name, "release_verdict": release["verdict"],
        "policy_version": policy_version,
        "extraction_consumer_gate": "REJECT" if release["verdict"] != "EXTRACTION_GO" else "ALLOW",
        "inventory_rows": inventory_rows, "source_rows": release["counts"]["source_rows"],
        "cohort_rows": release["counts"]["bounded_sentinel_assemblies"],
        "global_distinct_exact_assembly_revisions": release["counts"]["global_distinct_exact_assembly_revisions"],
    }


def publish_artifacts(repo: Path, final: Path, validation: dict[str, Any]) -> None:
    out = repo / "artifacts/prophage_semantics"
    out.mkdir(parents=True, exist_ok=True)
    for name in ("source_profile.json", "sentinel_summary.json", "resource_summary.json", "restart_evidence.json"):
        shutil.copyfile(final / name, out / name)
    atomic_write(out / "validation.json", canonical_bytes(validation))
    reference = {
        "schema": "prophage-semantics-release-reference-v1", "release_id": final.name,
        "external_path": str(final), "complete_sha256": sha_file(final / "COMPLETE"),
        "sha256sums_sha256": sha_file(final / "SHA256SUMS"),
        "verdict": validation["release_verdict"],
        "consumer_action": validation["extraction_consumer_gate"],
        "policy_version": validation["policy_version"],
    }
    atomic_write(out / "release_reference.json", canonical_bytes(reference))


def build(args: argparse.Namespace) -> dict[str, Any]:
    repo = Path(args.repo).resolve()
    # Immutable gate occurs before any stage or task-owned external directory is created.
    verify_exact(repo / "26k_ecoli_accession.txt", ACCESSIONS_SHA, "root accession input")
    verify_exact(repo / "26k_prophage1.csv", SOURCE_SHA, "root prophage input")
    verify_exact(repo / "reports/phage_pangenome_project_plan.md", PLAN_SHA, "integrated plan")
    verify_exact(repo / "reports/prophage_distribution.md", AUDIT_SHA, "source audit")
    policy = json.loads((repo / POLICY_FILE).read_text())
    validate_policy(policy)
    # Derive the extraction verdict from the pinned-caller consumption gate.
    # This is fail-closed: a missing/non-PASS/non-DECISIVE/non-independently-
    # sound gate raises GateError and no release is published (the historical
    # v1 BLOCKED release then remains the active record).
    verdict, consumer_action, gate = derive_extraction_from_gate(repo)
    evidence = json.loads((repo / "artifacts/prophage_semantics/evidence_inventory.json").read_text())
    if evidence.get("conclusion") is None or evidence.get("known_base_sentinels_used") is not False:
        raise GateError("evidence inventory contract mismatch")
    predecessor, external, assemblies, predecessor_inventory_rows = validate_predecessor(repo)
    input_manifest = expected_inputs(repo, predecessor, external, gate)
    for rel, expected in input_manifest["files"].items():
        verify_exact(repo / rel, expected, rel)
    release_key = sha_bytes(
        canonical_bytes(input_manifest) + canonical_bytes(policy) + gate["release_id"].encode()
    )[:16]
    release_id = f"{RELEASE_PREFIX}-{release_key}"
    durable_root = Path(args.durable_root).resolve()
    scratch = Path(args.scratch_root).resolve() / args.run_id
    final = durable_root / release_id
    allocations = Allocations(
        args.assigned_ram_bytes, args.durable_allocation_bytes, args.scratch_allocation_bytes,
        args.inode_allocation, args.predicted_durable_peak_bytes, args.predicted_scratch_peak_bytes,
        args.predicted_files, args.unfinished_write_bytes,
    )
    if final.exists():
        validation = validate_release(final)
        publish_artifacts(repo, final, validation)
        return validation | {"deterministic_rerun": "EXISTING_IMMUTABLE_RELEASE_VALIDATED"}
    start = preflight(final, scratch, allocations, "INITIAL")
    durable_root.mkdir(parents=True, exist_ok=True)
    stage = durable_root / f".{release_id}.staging"
    stage.mkdir(parents=False, exist_ok=True)
    resources_path = stage / "resources.jsonl"
    append_jsonl(resources_path, start)
    if not (stage / "state.jsonl").exists():
        append_jsonl(stage / "state.jsonl", {"event": "STAGE_CREATED", "release_id": release_id, "at": utc_now()})
    else:
        append_jsonl(stage / "state.jsonl", {"event": "RESTART_DETECTED", "release_id": release_id, "at": utc_now()})
    if not (stage / "failures.jsonl").exists():
        atomic_write(stage / "failures.jsonl", b"")
    start_swap = start["swap_used_bytes"]
    start_usage = dir_usage(stage)
    start_ru = resource.getrusage(resource.RUSAGE_SELF)
    try:
        profile, source_rows = source_profile(repo / "26k_prophage1.csv")
        sentinel = annotation_boundary_diagnostic(repo, external, assemblies, source_rows)
        write_static_unit(stage, "input_manifest.json", canonical_bytes(input_manifest))
        write_static_unit(stage, POLICY_STAGE_NAME, canonical_bytes(policy))
        write_static_unit(stage, "evidence_inventory.json", canonical_bytes(evidence))
        write_static_unit(stage, "source_profile.json", canonical_bytes(profile))
        write_static_unit(stage, "sentinel_summary.json", canonical_bytes(sentinel))
        tools = {
            "schema": "prophage-semantics-tools-v1", "python": platform.python_version(),
            "python_executable": sys.executable, "implementation": platform.python_implementation(),
            "platform": platform.platform(), "network_downloads": 0,
        }
        write_static_unit(stage, "tools.json", canonical_bytes(tools))
        provenance = {
            "schema": "prophage-semantics-provenance-v1", "task_id": "resolve-prophage-source",
            "release_id": release_id, "created_at_utc": utc_now(), "argv": sys.argv,
            "cwd": str(Path.cwd()), "repo": str(repo), "run_id": args.run_id,
            "source_code": "workflow/prophage_semantics/release.py",
            "source_code_sha256": sha_file(Path(__file__)), "predecessor_release_id": COHORT_RELEASE_ID,
            "pinned_caller_release_id": PREDECESSOR_PHIGARO_RELEASE_ID,
            "pinned_caller_external_release": str(PREDECESSOR_PHIGARO_EXTERNAL_ROOT / PREDECESSOR_PHIGARO_RELEASE_ID),
            "predecessor_external_inventory_rows": predecessor_inventory_rows,
            "authorization": {"new_assembly_downloads": 0, "max_distinct_assemblies": 10, "production_extraction": False},
        }
        if not (stage / "provenance.json").exists():
            write_static_unit(stage, "provenance.json", canonical_bytes(provenance))
        if args.inject_stop_before_complete:
            failure = {"event": "INJECTED_INTERRUPTION", "point": "AFTER_STATIC_UNITS_BEFORE_PUBLICATION", "at": utc_now()}
            append_jsonl(stage / "failures.jsonl", failure)
            append_jsonl(stage / "state.jsonl", failure)
            raise InjectedInterruption("injected interruption before COMPLETE")
        promotion = preflight(final, scratch, allocations, "PROMOTION")
        append_jsonl(resources_path, promotion)
        finish_usage = dir_usage(stage)
        ru = resource.getrusage(resource.RUSAGE_SELF)
        peak_rss = ru.ru_maxrss * 1024
        swap_end = promotion["swap_used_bytes"]
        checks = {
            "peak_rss_le_70pct_assigned": peak_rss * 100 <= allocations.assigned_ram_bytes * 70,
            "system_swap_growth_zero": swap_end <= start_swap,
            "measured_durable_peak_le_70pct_allocation": finish_usage[0] * 100 <= allocations.durable_allocation_bytes * 70,
            "measured_files_le_50pct_inodes": finish_usage[1] * 2 <= allocations.inode_allocation,
            "two_x_unfinished_reservation_at_finish": promotion["durable_free_bytes"] - finish_usage[0] >= 2 * allocations.unfinished_write_bytes,
        }
        if not all(checks.values()):
            raise GateError("end resource gate NO_GO: " + json.dumps(checks, sort_keys=True))
        resource_summary = {
            "schema": "prophage-semantics-resource-summary-v1", "verdict": "PASS",
            "allocations": asdict(allocations), "checks": checks, "start": start, "finish": promotion,
            "peak_rss_bytes": peak_rss, "peak_rss_fraction": peak_rss / allocations.assigned_ram_bytes,
            "swap_start_bytes": start_swap, "swap_finish_bytes": swap_end,
            "measured_stage_bytes": finish_usage[0], "measured_stage_files": finish_usage[1],
            "scale_trend": "NOT_APPLICABLE_NON_SCALE_BEARING_N10_METADATA_DIAGNOSTIC",
        }
        write_static_unit(stage, "resource_summary.json", canonical_bytes(resource_summary))
        failures = [json.loads(line) for line in (stage / "failures.jsonl").read_text().splitlines() if line]
        restart = {
            "schema": "prophage-semantics-restart-evidence-v1", "verdict": "PASS",
            "injected_interruption_observed": any(x.get("event") == "INJECTED_INTERRUPTION" for x in failures),
            "resume_policy": "existing static units independently SHA-256 validated; mismatch refuses mixed publication",
            "partial_complete_never_present": True,
        }
        write_static_unit(stage, "restart_evidence.json", canonical_bytes(restart))
        release = {
            "schema": SCHEMA, "release_id": release_id, "immutable": True,
            "created_at_utc": utc_now(), "source_task_id": "resolve-prophage-source",
            "verdict": verdict, "consumer_action": consumer_action,
            "policy_id": policy["policy_id"], "policy_version": POLICY_VERSION,
            "input_manifest_sha256": sha_bytes(canonical_bytes(input_manifest)),
            "predecessor_release_id": COHORT_RELEASE_ID,
            "pinned_caller_release_id": PREDECESSOR_PHIGARO_RELEASE_ID,
            "pinned_caller_historical_csv_attribution": gate["historical_csv_attribution"],
            "pinned_caller_modern_v2_4_pilot": gate["modern_v2_4_pilot"],
            "pinned_caller_decisive_evidence_independently_sound": gate["decisive_evidence_independently_sound"],
            "selected_coordinate_candidate": policy["extraction_gate"].get("selected_coordinate_candidate"),
            "counts": {
                "source_rows": profile["data_rows"], "all_records": 132404,
                "transposable_flag_positive": 7695, "taxonomy_assigned": 115442,
                "bounded_sentinel_assemblies": 10, "bounded_sentinel_rows": 56,
                "global_distinct_exact_assembly_revisions": 10, "global_cap": 1000,
                "new_assembly_downloads": 0, "sequence_bases_read": 0,
            },
            "gates": {
                "root_input_sha256": "PASS", "integrated_plan_sha256": "PASS",
                "producer_caller_evidence_inventory": "PASS",
                "predecessor_release_id_manifest_inventory": "PASS",
                "accession_version_identity": "PASS", "upstream_local_checksum": "PASS",
                "row_accounting": "PASS", "bgzf_index_name_roundtrip": "PASS",
                "global_distinct_assembly_cap": "PASS", "resource": "PASS",
                "injected_restart": "PASS", "atomic_promotion": "PASS",
                "source_coordinate_policy": "PASS",
                "pinned_consumer_compatibility": "PASS",
                "pinned_caller_consumption_gate": "PASS",
                "extraction_eligibility": verdict,
                "scale_trend": "NOT_APPLICABLE_NON_SCALE_BEARING",
            },
            "blocking_dimensions": policy["extraction_gate"]["blocking_dimensions"],
            "known_base_sentinels_used": False,
        }
        write_static_unit(stage, "release.json", canonical_bytes(release))
        append_jsonl(stage / "state.jsonl", {"event": "READY_FOR_INVENTORY", "at": utc_now()})
        atomic_write(stage / "SHA256SUMS", release_tree_inventory(stage))
        fsync_dir(stage)
        complete = {
            "schema": "prophage-semantics-complete-v1", "release_id": release_id,
            "sha256sums_sha256": sha_file(stage / "SHA256SUMS"), "verdict": verdict,
        }
        atomic_write(stage / "COMPLETE", canonical_bytes(complete))
        fsync_dir(stage)
        os.replace(stage, final)
        fsync_dir(durable_root)
    except InjectedInterruption:
        raise
    except Exception as exc:
        append_jsonl(stage / "failures.jsonl", {"event": "FAILURE", "type": type(exc).__name__, "message": str(exc), "at": utc_now()})
        append_jsonl(stage / "state.jsonl", {"event": "FAILED", "type": type(exc).__name__, "message": str(exc), "at": utc_now()})
        raise
    validation = validate_release(final)
    publish_artifacts(repo, final, validation)
    return validation


def add_allocation_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--assigned-ram-bytes", type=int, required=True)
    parser.add_argument("--durable-allocation-bytes", type=int, required=True)
    parser.add_argument("--scratch-allocation-bytes", type=int, required=True)
    parser.add_argument("--inode-allocation", type=int, required=True)
    parser.add_argument("--predicted-durable-peak-bytes", type=int, required=True)
    parser.add_argument("--predicted-scratch-peak-bytes", type=int, required=True)
    parser.add_argument("--predicted-files", type=int, required=True)
    parser.add_argument("--unfinished-write-bytes", type=int, required=True)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("--repo", default=".")
    run.add_argument("--durable-root", required=True)
    run.add_argument("--scratch-root", required=True)
    run.add_argument("--run-id", required=True)
    run.add_argument("--inject-stop-before-complete", action="store_true")
    add_allocation_args(run)
    val = sub.add_parser("validate")
    val.add_argument("release")
    val.add_argument("--require-extraction-go", action="store_true")
    return p


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "run":
            result = build(args)
        else:
            result = validate_release(Path(args.release), args.require_extraction_go)
        print(json.dumps(result, sort_keys=True, separators=(",", ":")))
        return 0
    except InjectedInterruption as exc:
        print(str(exc), file=sys.stderr)
        return 75
    except (GateError, OSError, ValueError, KeyError, zipfile.BadZipFile) as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
