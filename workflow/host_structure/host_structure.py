#!/usr/bin/env python3
"""Shared fail-closed primitives for the frozen N=1,000 host-only run.

Only standard-library code lives here so the independent validator does not
need the biological tool environment.  The workflow hashes but never parses
the immutable prophage root file; no phage-derived field is an allowed input.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
import resource
import shutil
import stat
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

PASS = "PASS"
NA_HOST = "NOT_APPLICABLE_HOST_ONLY_ANALYSIS_EXTRACTION_BLOCKED"
NA_RUNG = "NOT_APPLICABLE_ONE_RUNG_HOST_ONLY_TASK"
TASK_ID = "run-host-structure-1000"
SCHEMA = "host-structure-1000-release-v1"
EXPECTED_N = 1000
EXPECTED_PAIRS = 499_500
ROOT_HASHES = {
    "26k_ecoli_accession.txt": "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5",
    "26k_prophage1.csv": "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996",
}
SELECTION_RELEASE_ID = "pilot-cohorts-v1-8afc0ea03d9e50dc"
SELECTION_RELEASE_JSON_SHA256 = "d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99"
COHORT_SHA256 = "265a1e7784a4d5db3ea3577892feba8173290518b6c621f7e5091dbad66bfe77"
CANONICAL_RELEASE_ID = "canonical-cohort-1000-v1-4bc3e029e6e0be44"
CANONICAL_RELEASE_JSON_SHA256 = "14a39b424f2a23de6fa52c173b00e03b167e897baf3a9dbcd9876e31e999740c"
CANONICAL_EXTERNAL_PATH = Path("/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-1000/canonical-cohort-1000-v1-4bc3e029e6e0be44")
CANONICAL_TRACKED_INVENTORY_SHA256 = "7a819c283c49a18be0ad5e3ec1138b4e70e35524da73d0c27725933806cad01e"
GLOBAL_CAP_EVIDENCE_SHA256 = "8de39ac8cd6df912736085b5299e345242724a88891ca117d786f2d83eb9c769"
COMPATIBILITY_RELEASE_ID = "consumer-compatibility-v1-78d7e93f19fa3d87"
COMPATIBILITY_RELEASE_JSON_SHA256 = "021719ddadd7bb7fa2932d2ef9cb25da9c666ebe0389988691283011ee12f4c7"
COMPATIBILITY_TRACKED_INVENTORY_SHA256 = "fa503af287f9361d2d350177ce8e54e3859fb255b6d33be10d155d4707960823"
HOST_LOCK_SHA256 = "0b8db59d01eed5762db2bcb52e581e66b61750988c6db98cc36cc4721e53ccd4"
HOST_PACKAGE_INVENTORY_SHA256 = "b473c14076ebf01c35e3a496a4fde403ce7bf90a874bb20ec4c17afb01c8c34a"
MINIMAP2_SHA256 = "b6c81294dc0b68b2f54f8e2f6f3ad6be71a40bccc47a6e244ecebc73bae9501d"
PINNED_TOOL_SHA256 = {
    "mash": "83c85c063118c8c12659baa5d990aa8821f65cde40d1c82c8eaedd144dfec205",
    "rapidnj": "73bbf9615f3592d540084241634b4f09f1205e75107bc529efc4359453dd0208",
    "skani": "b1d20cb7170fe40a964526eadeb6fcdf61eefa748b88ed701ec3ff8dfbc07f5f",
    "minimap2": MINIMAP2_SHA256,
}
DURABLE_PREFIX = Path("/home/erikg/phind-data/ecoli26k/v1/releases/run-host-structure-1000")
SCRATCH_PREFIX = Path("/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-host-structure-1000")

# Biological commands may receive only these host-only sources.  The prophage
# root appears solely in ROOT_HASHES and is read as opaque bytes by root audit.
HOST_ONLY_INPUT_ALLOWLIST = (
    "manifests/pilot-cohorts-v1/",
    "manifests/canonical-cohort-1000-v1/",
    "manifests/consumer-compatibility-v1/",
    "workflow/compatibility/environment-linux-64.explicit.lock",
    "workflow/compatibility/environment-package-sha256.tsv",
)
FORBIDDEN_BIOLOGICAL_TOKENS = (
    "prophage", "phage", "viral", "cluster_id", "integrase", "att_site",
    "taxonomy", "prophage_count", "coordinates",
)
ACCESSION_RE = re.compile(r"GC[AF]_\d{9}\.\d+")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


class GateError(RuntimeError):
    """A hard, non-waivable gate failure."""


def sha256_file(path: Path, chunk: int = 4 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while data := handle.read(chunk):
            h.update(data)
    return h.hexdigest()


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode() + b"\n"


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json(value))


def append_jsonl(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("ab") as handle:
        handle.write(canonical_json(value))
        handle.flush()
        os.fsync(handle.fileno())


def stable_row_hash(row: dict[str, str], fields: Sequence[str]) -> str:
    return hashlib.sha256(("\t".join(str(row.get(f, ".")) for f in fields if f != "row_sha256") + "\n").encode()).hexdigest()


def write_hashed_tsv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    if not fields or fields[-1] != "row_sha256":
        raise GateError("hashed TSV schema must end in row_sha256")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        for raw in rows:
            row = {f: str(raw.get(f, "")) for f in fields}
            row["row_sha256"] = stable_row_hash(row, fields)
            writer.writerow(row)


def read_hashed_tsv(path: Path, expected_fields: Sequence[str] | None = None) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    if not fields or fields[-1] != "row_sha256":
        raise GateError(f"row checksum column missing: {path}")
    if expected_fields is not None and list(fields) != list(expected_fields):
        raise GateError(f"unexpected TSV schema: {path}")
    for number, row in enumerate(rows, 2):
        if row.get("row_sha256") != stable_row_hash(row, fields):
            raise GateError(f"row checksum mismatch: {path}:{number}")
    return rows


def verify_tracked_inventory(root: Path, expected_digest: str | None = None) -> dict[str, str]:
    sums_path = root / "SHA256SUMS"
    if expected_digest and sha256_file(sums_path) != expected_digest:
        raise GateError(f"tracked inventory digest mismatch: {root}")
    result: dict[str, str] = {}
    for number, line in enumerate(sums_path.read_text().splitlines(), 1):
        try:
            digest, rel = line.split("  ", 1)
        except ValueError as exc:
            raise GateError(f"malformed inventory line {number}: {sums_path}") from exc
        p = Path(rel)
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or p.is_absolute() or ".." in p.parts or rel in result:
            raise GateError(f"unsafe inventory line {number}: {sums_path}")
        target = root / p
        if not target.is_file() or target.is_symlink() or sha256_file(target) != digest:
            raise GateError(f"tracked checksum mismatch: {target}")
        result[rel] = digest
    return result


def verify_external_inventory(root: Path) -> dict[str, str]:
    sums_path = root / "SHA256SUMS"
    complete = root / "COMPLETE"
    if not sums_path.is_file() or not complete.is_file() or complete.is_symlink():
        raise GateError(f"external release lacks immutable completion contract: {root}")
    result: dict[str, str] = {}
    for number, line in enumerate(sums_path.read_text().splitlines(), 1):
        try:
            digest, rel = line.split("  ", 1)
        except ValueError as exc:
            raise GateError(f"malformed external inventory line {number}: {root}") from exc
        p = Path(rel)
        if p.is_absolute() or ".." in p.parts or rel in result:
            raise GateError(f"unsafe external inventory path: {rel}")
        target = root / p
        if not target.is_file() or target.is_symlink() or sha256_file(target) != digest:
            raise GateError(f"external checksum mismatch: {target}")
        result[rel] = digest
    marker = complete.read_text().strip()
    sums_digest = sha256_file(sums_path)
    if marker.startswith("{"):
        value = json.loads(marker)
        if value.get("sha256sums_sha256") != sums_digest or value.get("verdict") != PASS:
            raise GateError(f"COMPLETE marker mismatch: {root}")
    elif not marker.split() or marker.split()[0] != sums_digest:
        raise GateError(f"COMPLETE marker mismatch: {root}")
    return result


def require_pass_or_explicit_na(gates: dict[str, Any], label: str) -> None:
    for gate, verdict in gates.items():
        text = str(verdict)
        if text == PASS or text.startswith("NOT_APPLICABLE"):
            continue
        raise GateError(f"{label} applicable gate not unqualified PASS/NA: {gate}={verdict}")


def verify_root_hashes(repo: Path) -> dict[str, str]:
    # Opaque hashing only.  The prophage CSV is never opened as text or parsed.
    observed = {name: sha256_file(repo / name) for name in ROOT_HASHES}
    if observed != ROOT_HASHES:
        raise GateError(f"immutable root input checksum mismatch: {observed}")
    return observed


def _release_path(value: dict[str, Any]) -> Path:
    raw = value.get("external_release_path", value.get("external_path"))
    if not isinstance(raw, str) or not raw:
        raise GateError("release does not declare an external path")
    return Path(raw)


def verify_inputs(repo: Path) -> dict[str, Any]:
    roots = verify_root_hashes(repo)

    selection_root = repo / "manifests/pilot-cohorts-v1"
    selection_inventory = verify_tracked_inventory(selection_root)
    if selection_inventory.get("release.json") != SELECTION_RELEASE_JSON_SHA256:
        raise GateError("selection release pin mismatch")
    if selection_inventory.get("cohort-1000.tsv") != COHORT_SHA256:
        raise GateError("selection cohort pin mismatch")
    selection = json.loads((selection_root / "release.json").read_text())
    if selection.get("release_id") != SELECTION_RELEASE_ID or selection.get("verdict") != PASS or selection.get("immutable") is not True:
        raise GateError("selection release identity/verdict/immutability mismatch")
    require_pass_or_explicit_na(selection.get("applicable_gates", {}), "selection")
    verify_external_inventory(_release_path(selection))

    canonical_root = repo / "manifests/canonical-cohort-1000-v1"
    canonical_inventory = verify_tracked_inventory(canonical_root, CANONICAL_TRACKED_INVENTORY_SHA256)
    if canonical_inventory.get("release.json") != CANONICAL_RELEASE_JSON_SHA256:
        raise GateError("canonical release pin mismatch")
    if canonical_inventory.get("cohort-1000.tsv") != COHORT_SHA256:
        raise GateError("canonical cohort pin mismatch")
    canonical = json.loads((canonical_root / "release.json").read_text())
    if canonical.get("release_id") != CANONICAL_RELEASE_ID or canonical.get("immutable") is not True:
        raise GateError("canonical release identity/immutability mismatch")
    require_pass_or_explicit_na(canonical.get("applicable_gates", {}), "canonical")
    counts = canonical.get("counts", {})
    if counts.get("validated") != EXPECTED_N or counts.get("distinct_sequence_bearing_assemblies") != EXPECTED_N:
        raise GateError("canonical release is not exact N=1,000")
    canonical_external = _release_path(canonical)
    if canonical_external != CANONICAL_EXTERNAL_PATH:
        raise GateError("canonical external release path pin mismatch")
    verify_external_inventory(canonical_external)
    if sha256_file(canonical_external / "release.json") != CANONICAL_RELEASE_JSON_SHA256:
        raise GateError("tracked/external canonical release mismatch")

    compatibility_root = repo / "manifests/consumer-compatibility-v1"
    compatibility_inventory = verify_tracked_inventory(compatibility_root, COMPATIBILITY_TRACKED_INVENTORY_SHA256)
    if compatibility_inventory.get("release.json") != COMPATIBILITY_RELEASE_JSON_SHA256:
        raise GateError("consumer compatibility release pin mismatch")
    compatibility = json.loads((compatibility_root / "release.json").read_text())
    if (compatibility.get("release_id") != COMPATIBILITY_RELEASE_ID
            or compatibility.get("verdict") != PASS or compatibility.get("immutable") is not True):
        raise GateError("consumer compatibility identity/verdict/immutability mismatch")
    require_pass_or_explicit_na(compatibility.get("applicable_gates", {}), "consumer compatibility")
    compatibility_external = _release_path(compatibility)
    verify_external_inventory(compatibility_external)
    if sha256_file(compatibility_external / "release.json") != COMPATIBILITY_RELEASE_JSON_SHA256:
        raise GateError("tracked/external compatibility release mismatch")
    with (compatibility_root / "gates.tsv").open(newline="") as gate_handle:
        gates = list(csv.DictReader(gate_handle, delimiter="\t"))
    needed = {"bgzip", "samtools-faidx", "mash", "rapidnj", "skani"}
    observed_needed = {r["gate_id"] for r in gates if r.get("verdict") == PASS}
    if not needed <= observed_needed or any(r.get("verdict") != PASS for r in gates):
        raise GateError("required host consumer compatibility gates are not all PASS")

    if sha256_file(repo / "workflow/compatibility/environment-linux-64.explicit.lock") != HOST_LOCK_SHA256:
        raise GateError("host environment lock digest mismatch")
    if sha256_file(repo / "workflow/compatibility/environment-package-sha256.tsv") != HOST_PACKAGE_INVENTORY_SHA256:
        raise GateError("host package inventory digest mismatch")

    cohort = read_hashed_tsv(canonical_root / "cohort-1000.tsv")
    assemblies = read_hashed_tsv(canonical_root / "assemblies.tsv")
    refs = read_hashed_tsv(canonical_root / "object_refs.tsv")
    if not (len(cohort) == len(assemblies) == len(refs) == EXPECTED_N):
        raise GateError("canonical row accounting is not 100%")
    accessions = [r["exact_assembly_accession_version"] for r in cohort]
    if len(set(accessions)) != EXPECTED_N or accessions != [r["accession"] for r in assemblies] or accessions != [r["accession"] for r in refs]:
        raise GateError("canonical accession order/identity accounting mismatch")
    cap_path = repo / "artifacts/canonical_cohort_1000/global_cap_evidence.json"
    if sha256_file(cap_path) != GLOBAL_CAP_EVIDENCE_SHA256:
        raise GateError("global distinct-assembly cap evidence pin mismatch")
    cap = json.loads(cap_path.read_text())
    projected = cap.get("projected_finish_accessions", [])
    if (cap.get("verdict") != PASS or cap.get("graph_cap") != EXPECTED_N
            or cap.get("projected_finish_distinct_exact_assembly_revisions") != EXPECTED_N
            or set(projected) != set(accessions) or len(projected) != EXPECTED_N):
        raise GateError("global distinct-assembly cap is not exact frozen N=1,000 PASS")
    for order, (c, a, ref) in enumerate(zip(cohort, assemblies, refs), 1):
        if int(c["cohort_order"]) != order or int(a["stage_b_order"]) != order or int(ref["cohort_order"]) != order:
            raise GateError("canonical cohort order is not immutable 1..1,000")
        if (not ACCESSION_RE.fullmatch(accessions[order - 1]) or a["terminal_state"] != "VALIDATED"
                or ref["predecessor_digest_match"] not in (PASS, "NOT_APPLICABLE_NEW_OBJECT")):
            raise GateError("canonical exact-version/validation/reference gate failed")

    return {
        "root_hashes": roots,
        "selection": selection,
        "canonical": canonical,
        "canonical_external": canonical_external,
        "compatibility": compatibility,
        "compatibility_external": compatibility_external,
        "cohort": cohort,
        "assemblies": assemblies,
        "refs": refs,
        "accessions": accessions,
    }


def canonical_object_path(ref: dict[str, str], accession: str) -> Path:
    root = Path(ref["storage_root"])
    # SELF is intentionally represented as '.' in the immutable object-ref
    # table and resolves only against the checksum-pinned canonical N=1,000
    # external release, never against the caller's working directory.
    if ref.get("storage_release_id") == "SELF":
        if root != Path("."):
            raise GateError("SELF object reference has noncanonical storage root")
        root = CANONICAL_EXTERNAL_PATH
    elif not root.is_absolute():
        raise GateError("predecessor object reference is not absolute")
    obj = root / ref["canonical_object_relpath"]
    return obj / f"{accession}.pansn.fa.gz"


def assert_host_only_manifest(value: Any, path: str = "$") -> None:
    """Reject phage-trait keys/paths while permitting the opaque root hash label."""
    if isinstance(value, dict):
        for key, child in value.items():
            lower = str(key).lower()
            if any(token in lower for token in FORBIDDEN_BIOLOGICAL_TOKENS):
                opaque_root = key == "26k_prophage1.csv" and child == ROOT_HASHES["26k_prophage1.csv"]
                explicit_block = key in {"prophage_source_coordinate_policy", "prophage_extraction_semantics"} and child == NA_HOST
                if not (opaque_root or explicit_block):
                    raise GateError(f"phage-derived key absent gate failed at {path}.{key}")
            assert_host_only_manifest(child, f"{path}.{key}")
    elif isinstance(value, list):
        for i, child in enumerate(value):
            assert_host_only_manifest(child, f"{path}[{i}]")
    elif isinstance(value, str):
        lower = value.lower()
        looks_like_artifact = ("/" in value or "\\" in value or Path(value).suffix.lower() in {".tsv", ".csv", ".fa", ".fasta", ".gz", ".json", ".msh"})
        if (value != NA_HOST and "26k_prophage1.csv" not in lower and looks_like_artifact
                and any(token in Path(lower).name for token in ("prophage", "phage", "viral"))):
            raise GateError(f"phage artifact absent gate failed at {path}")


@dataclass(frozen=True)
class Allocations:
    assigned_ram_bytes: int
    durable_allocation_bytes: int
    scratch_allocation_bytes: int
    inode_allocation: int
    predicted_durable_upper95_bytes: int
    predicted_scratch_upper95_bytes: int
    predicted_files: int
    unfinished_write_bytes: int

    def validate(self) -> None:
        values = asdict(self)
        if any(not isinstance(v, int) or v <= 0 for v in values.values()):
            raise GateError(f"blank/invalid allocation is NO_GO: {values}")
        if self.predicted_durable_upper95_bytes > 0.70 * self.durable_allocation_bytes:
            raise GateError("durable upper-95 exceeds 70% allocation")
        if self.predicted_scratch_upper95_bytes > 0.70 * self.scratch_allocation_bytes:
            raise GateError("scratch upper-95 exceeds 70% allocation")
        if self.predicted_files > 0.50 * self.inode_allocation:
            raise GateError("projected files exceed 50% inode allocation")
        if self.durable_allocation_bytes - self.predicted_durable_upper95_bytes < 2 * self.unfinished_write_bytes:
            raise GateError("durable allocation lacks 2x unfinished-write reserve")
        if self.scratch_allocation_bytes - self.predicted_scratch_upper95_bytes < 2 * self.unfinished_write_bytes:
            raise GateError("scratch allocation lacks 2x unfinished-write reserve")


def _mount_record(path: Path) -> dict[str, str]:
    proc = subprocess.run(
        ["findmnt", "-J", "-T", str(path), "-o", "TARGET,SOURCE,FSTYPE,OPTIONS,SIZE,AVAIL"],
        check=True, text=True, capture_output=True,
    )
    fs = json.loads(proc.stdout)["filesystems"][0]
    return {k: str(v) for k, v in fs.items()}


def _write_probe(parent: Path) -> None:
    parent.mkdir(parents=True, exist_ok=True)
    if parent.stat().st_uid != os.getuid():
        raise GateError(f"task path is not owned by current uid: {parent}")
    fd, name = tempfile.mkstemp(prefix=".write-probe-", dir=parent)
    try:
        os.write(fd, b"host-structure-write-probe\n")
        os.fsync(fd)
    finally:
        os.close(fd)
        Path(name).unlink(missing_ok=True)


def cgroup_memory_snapshot() -> dict[str, int]:
    rel = "/"
    for line in Path("/proc/self/cgroup").read_text().splitlines():
        if line.startswith("0::"):
            rel = line.split("::", 1)[1]
    root = Path("/sys/fs/cgroup") / rel.lstrip("/")
    result: dict[str, int] = {}
    for name in ("memory.current", "memory.peak", "memory.swap.current"):
        try:
            result[name.replace(".", "_")] = int((root / name).read_text().strip())
        except (FileNotFoundError, ValueError):
            result[name.replace(".", "_")] = 0
    try:
        for line in (root / "memory.events").read_text().splitlines():
            key, value = line.split()
            result[f"memory_events_{key}"] = int(value)
    except FileNotFoundError:
        pass
    return result


def system_swap_used() -> int:
    total = 0
    lines = Path("/proc/swaps").read_text().splitlines()[1:]
    for line in lines:
        fields = line.split()
        if len(fields) >= 4:
            total += int(fields[3]) * 1024
    return total


def live_preflight(durable_parent: Path, scratch_root: Path, allocations: Allocations, stage: str) -> dict[str, Any]:
    allocations.validate()
    durable_parent.mkdir(parents=True, exist_ok=True)
    scratch_root.mkdir(parents=True, exist_ok=True)
    _write_probe(durable_parent)
    _write_probe(scratch_root)
    ds, ss = os.statvfs(durable_parent), os.statvfs(scratch_root)
    durable_free = ds.f_bavail * ds.f_frsize
    scratch_free = ss.f_bavail * ss.f_frsize
    checks = {
        "durable_free_ge_2tb": durable_free >= 2_000_000_000_000,
        "durable_free_inodes_ge_1m": ds.f_favail >= 1_000_000,
        "scratch_live_preflight_ge_4tb": scratch_free >= 4_000_000_000_000,
        "scratch_free_inodes_ge_5m": ss.f_favail >= 5_000_000,
        "scratch_stop_floor_ge_2tb": scratch_free >= 2_000_000_000_000,
        "durable_live_retains_2x_unfinished": durable_free - allocations.predicted_durable_upper95_bytes >= 2 * allocations.unfinished_write_bytes,
        "scratch_live_retains_2x_unfinished": scratch_free - allocations.predicted_scratch_upper95_bytes >= 2 * allocations.unfinished_write_bytes,
        "durable_predicted_le_70pct": allocations.predicted_durable_upper95_bytes <= 0.70 * allocations.durable_allocation_bytes,
        "scratch_predicted_le_70pct": allocations.predicted_scratch_upper95_bytes <= 0.70 * allocations.scratch_allocation_bytes,
        "projected_files_le_50pct": allocations.predicted_files <= 0.50 * allocations.inode_allocation,
    }
    record = {
        "schema": "host-structure-resource-preflight-v1", "stage": stage,
        "captured_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "allocations": asdict(allocations), "durable_path": str(durable_parent),
        "scratch_path": str(scratch_root), "durable_findmnt": _mount_record(durable_parent),
        "scratch_findmnt": _mount_record(scratch_root), "durable_free_bytes": durable_free,
        "scratch_free_bytes": scratch_free, "durable_free_inodes": ds.f_favail,
        "scratch_free_inodes": ss.f_favail, "durable_owner": durable_parent.stat().st_uid,
        "scratch_owner": scratch_root.stat().st_uid, "write_probes": PASS,
        "checks": checks, "cgroup": cgroup_memory_snapshot(), "system_swap_used_bytes": system_swap_used(),
        "verdict": PASS if all(checks.values()) else "NO_GO",
    }
    if record["verdict"] != PASS:
        raise GateError(f"resource gate NO_GO before {stage}: {record}")
    return record


def tree_usage(root: Path) -> tuple[int, int]:
    total = files = 0
    if not root.exists():
        return 0, 0
    for base, dirs, names in os.walk(root):
        files += len(names)
        for name in names:
            try:
                total += (Path(base) / name).stat().st_size
            except FileNotFoundError:
                pass
    return total, files


def seal_directory(stage: Path, final: Path) -> None:
    if final.exists():
        raise GateError(f"refusing to overwrite immutable release: {final}")
    if not stage.is_dir() or (stage / "COMPLETE").exists():
        raise GateError("invalid staging directory for promotion")
    lines: list[str] = []
    for path in sorted(p for p in stage.rglob("*") if p.is_file() and p.name not in {"SHA256SUMS", "COMPLETE"}):
        rel = path.relative_to(stage).as_posix()
        lines.append(f"{sha256_file(path)}  {rel}\n")
    sums = stage / "SHA256SUMS"
    sums.write_text("".join(lines))
    with sums.open("rb") as handle:
        os.fsync(handle.fileno())
    marker = {"schema": "host-structure-complete-v1", "sha256sums_sha256": sha256_file(sums), "verdict": PASS}
    complete = stage / "COMPLETE"
    complete.write_bytes(canonical_json(marker))
    with complete.open("rb") as handle:
        os.fsync(handle.fileno())
    dirfd = os.open(stage, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(dirfd)
    finally:
        os.close(dirfd)
    os.replace(stage, final)
    parentfd = os.open(final.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(parentfd)
    finally:
        os.close(parentfd)


def parse_mash_triangle(path: Path, expected_labels: Sequence[str] | None = None) -> tuple[list[str], list[list[float]]]:
    lines = path.read_text().splitlines()
    if not lines:
        raise GateError(f"empty Mash triangle: {path}")
    try:
        n = int(lines[0].strip())
    except ValueError as exc:
        raise GateError(f"invalid Mash triangle count: {path}") from exc
    if len(lines) != n + 1:
        raise GateError(f"Mash triangle row count mismatch: {path}")
    labels: list[str] = []
    matrix = [[0.0] * n for _ in range(n)]
    pairs = 0
    for i, line in enumerate(lines[1:]):
        fields = line.split()
        if len(fields) != i + 1:
            raise GateError(f"Mash triangle row {i} has {len(fields)-1}, expected {i} distances")
        label = Path(fields[0]).name
        if label.endswith(".fa"):
            label = label[:-3]
        labels.append(label)
        for j, token in enumerate(fields[1:]):
            value = float(token)
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise GateError(f"nonfinite/out-of-range Mash distance at {i},{j}")
            matrix[i][j] = matrix[j][i] = value
            pairs += 1
    if n != EXPECTED_N or pairs != EXPECTED_PAIRS or len(set(labels)) != n:
        raise GateError(f"Mash cardinality/pair gate failed: n={n}, pairs={pairs}")
    if expected_labels is not None and labels != list(expected_labels):
        raise GateError("Mash labels/order differ from frozen cohort")
    return labels, matrix


def write_full_phylip(path: Path, labels: Sequence[str], matrix: Sequence[Sequence[float]]) -> None:
    with path.open("w") as handle:
        handle.write(f"{len(labels)}\n")
        for label, row in zip(labels, matrix):
            if any(c.isspace() or c in "(),:;'[]" for c in label):
                raise GateError(f"unsafe PHYLIP label: {label}")
            handle.write(label + " " + " ".join(f"{x:.10g}" for x in row) + "\n")


def validate_directed_mash(path: Path, labels: Sequence[str], matrix: Sequence[Sequence[float]], tolerance: float = 5e-7) -> dict[str, Any]:
    index = {label: i for i, label in enumerate(labels)}
    seen: dict[tuple[int, int], float] = {}
    diagonal = 0
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise GateError(f"malformed mash dist record {number}")
            a, b = Path(fields[0]).name, Path(fields[1]).name
            a = a[:-3] if a.endswith(".fa") else a
            b = b[:-3] if b.endswith(".fa") else b
            if a not in index or b not in index:
                raise GateError("mash dist emitted unknown tip")
            i, j = index[a], index[b]
            value = float(fields[2])
            if (i, j) in seen or abs(value - matrix[i][j]) > tolerance:
                raise GateError(f"mash dist duplicate/matrix mismatch at {a},{b}")
            seen[i, j] = value
            if i == j:
                diagonal += 1
                if value != 0.0:
                    raise GateError("mash dist diagonal is not zero")
    expected = len(labels) ** 2
    if len(seen) != expected or diagonal != len(labels):
        raise GateError(f"directed Mash pair count mismatch: {len(seen)} != {expected}")
    for (i, j), value in seen.items():
        if abs(value - seen[j, i]) > tolerance:
            raise GateError("mash dist symmetry gate failed")
    return {"directed_records": len(seen), "diagonal_records": diagonal,
            "unordered_off_diagonal_pairs": len(labels) * (len(labels) - 1) // 2,
            "symmetry": PASS, "triangle_exact_match": PASS}


@dataclass
class TreeNode:
    name: str | None = None
    length: float | None = None
    children: list["TreeNode"] | None = None

    def is_leaf(self) -> bool:
        return not self.children


def parse_newick(text: str) -> TreeNode:
    """Parse the conservative RapidNJ/Newick subset used by this workflow."""
    text = text.strip()
    if not text.endswith(";"):
        raise GateError("Newick lacks terminal semicolon")
    i = 0

    def token() -> str:
        nonlocal i
        start = i
        while i < len(text) and text[i] not in ",():;":
            i += 1
        return text[start:i].strip().strip("'\"")

    def subtree() -> TreeNode:
        nonlocal i
        if text[i] == "(":
            i += 1
            children = [subtree()]
            while text[i] == ",":
                i += 1
                children.append(subtree())
            if text[i] != ")":
                raise GateError("malformed Newick internal node")
            i += 1
            name = token() or None
            node = TreeNode(name=name, children=children)
        else:
            name = token()
            if not name:
                raise GateError("empty Newick leaf")
            node = TreeNode(name=name, children=[])
        if i < len(text) and text[i] == ":":
            i += 1
            raw = token()
            try:
                node.length = float(raw)
            except ValueError as exc:
                raise GateError("invalid Newick branch length") from exc
        return node

    root = subtree()
    if i != len(text) - 1:
        raise GateError(f"unparsed Newick content at {i}")
    return root


def leaf_names(root: TreeNode) -> list[str]:
    if root.is_leaf():
        return [root.name or ""]
    result: list[str] = []
    for child in root.children or []:
        result.extend(leaf_names(child))
    return result


def canonical_split(side: Iterable[str], all_tips: frozenset[str]) -> frozenset[str] | None:
    a = frozenset(side)
    b = all_tips - a
    if len(a) < 2 or len(b) < 2:
        return None
    if len(a) < len(b):
        return a
    if len(b) < len(a):
        return b
    return min(a, b, key=lambda s: tuple(sorted(s)))


def tree_splits(root: TreeNode) -> set[frozenset[str]]:
    tips = frozenset(leaf_names(root))
    result: set[frozenset[str]] = set()

    def walk(node: TreeNode) -> frozenset[str]:
        if node.is_leaf():
            return frozenset([node.name or ""])
        here = frozenset().union(*(walk(c) for c in node.children or []))
        split = canonical_split(here, tips)
        if split:
            result.add(split)
        return here

    walk(root)
    return result


def newick_string(root: TreeNode) -> str:
    def emit(n: TreeNode) -> str:
        if n.is_leaf():
            base = n.name or ""
        else:
            base = "(" + ",".join(emit(c) for c in n.children or []) + ")" + (n.name or "")
        if n.length is not None:
            base += f":{max(0.0, n.length):.10g}"
        return base
    return emit(root) + ";\n"


def collapse_unsupported(root: TreeNode, support: dict[frozenset[str], float], threshold: float) -> TreeNode:
    all_tips = frozenset(leaf_names(root))

    def visit(node: TreeNode) -> tuple[TreeNode, frozenset[str]]:
        if node.is_leaf():
            return TreeNode(node.name, node.length, []), frozenset([node.name or ""])
        rebuilt: list[tuple[TreeNode, frozenset[str]]] = [visit(c) for c in node.children or []]
        children: list[TreeNode] = []
        here = frozenset().union(*(tips for _, tips in rebuilt))
        for child, child_tips in rebuilt:
            split = canonical_split(child_tips, all_tips)
            if split is not None and support.get(split, 0.0) < threshold and not child.is_leaf():
                for grand in child.children or []:
                    children.append(grand)
            else:
                children.append(child)
        return TreeNode(node.name, node.length, children), here

    return visit(root)[0]


def nearest_neighbors(labels: Sequence[str], matrix: Sequence[Sequence[float]]) -> dict[str, tuple[float, tuple[str, ...]]]:
    result: dict[str, tuple[float, tuple[str, ...]]] = {}
    for i, label in enumerate(labels):
        best = min(matrix[i][j] for j in range(len(labels)) if i != j)
        ties = tuple(labels[j] for j in range(len(labels)) if i != j and abs(matrix[i][j] - best) <= 1e-12)
        result[label] = best, ties
    return result


def rank_correlation_sample(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], limit: int = 50_000) -> float:
    pairs: list[tuple[float, float]] = []
    n = len(a)
    stride = max(1, (n * (n - 1) // 2) // limit)
    k = 0
    for i in range(n):
        for j in range(i):
            if k % stride == 0:
                pairs.append((a[i][j], b[i][j]))
            k += 1
    if len(pairs) < 2:
        return 1.0
    # Deterministic average ranks; ties are common for near-identical sketches.
    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda x: values[x])
        out = [0.0] * len(values)
        p = 0
        while p < len(order):
            q = p + 1
            while q < len(order) and values[order[q]] == values[order[p]]:
                q += 1
            rank = (p + q - 1) / 2.0
            for idx in order[p:q]:
                out[idx] = rank
            p = q
        return out
    ra, rb = ranks([x for x, _ in pairs]), ranks([y for _, y in pairs])
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db = math.sqrt(sum((y - mb) ** 2 for y in rb))
    return num / (da * db) if da and db else 1.0


def union_find_groups(labels: Sequence[str], matrix: Sequence[Sequence[float]], threshold: float) -> list[list[int]]:
    parent = list(range(len(labels)))
    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[max(a, b)] = min(a, b)
    for i in range(len(labels)):
        for j in range(i):
            if matrix[i][j] <= threshold:
                union(i, j)
    groups: dict[int, list[int]] = defaultdict(list)
    for i in range(len(labels)):
        groups[find(i)].append(i)
    return sorted(groups.values(), key=lambda g: g[0])


def farthest_first_partition(labels: Sequence[str], matrix: Sequence[Sequence[float]], k: int) -> tuple[list[int], list[int]]:
    if not 1 <= k <= len(labels):
        raise GateError("invalid medoid partition k")
    medoids = [0]
    while len(medoids) < k:
        candidate = max((i for i in range(len(labels)) if i not in medoids),
                        key=lambda i: (min(matrix[i][m] for m in medoids), -i))
        medoids.append(candidate)
    for _ in range(5):
        assignment = [min(range(k), key=lambda c: (matrix[i][medoids[c]], c)) for i in range(len(labels))]
        changed = False
        for c in range(k):
            members = [i for i, a in enumerate(assignment) if a == c]
            if not members:
                continue
            new = min(members, key=lambda x: (sum(matrix[x][y] for y in members), x))
            if new != medoids[c]:
                medoids[c] = new
                changed = True
        if not changed:
            break
    assignment = [min(range(k), key=lambda c: (matrix[i][medoids[c]], c)) for i in range(len(labels))]
    order = sorted(range(k), key=lambda c: medoids[c])
    remap = {old: new for new, old in enumerate(order)}
    return [medoids[c] for c in order], [remap[a] for a in assignment]


def select_high_fidelity_cases(
    labels: Sequence[str], matrix: Sequence[Sequence[float]], medoids: Sequence[int],
    assignment: Sequence[int], assemblies: Sequence[dict[str, str]], per_lineage: int = 12,
) -> dict[int, list[int]]:
    selected: dict[int, list[int]] = {}
    for lineage, medoid in enumerate(medoids):
        members = [i for i, a in enumerate(assignment) if a == lineage]
        chosen: list[int] = [medoid]
        # Diverse/boundary, then QC extremes, then maximin genetic coverage.
        if len(members) > 1:
            chosen.append(max(members, key=lambda i: (matrix[i][medoid], -i)))
            chosen.append(min(members, key=lambda i: (matrix[i][medoid], i)))
            chosen.append(max(members, key=lambda i: (int(assemblies[i]["contig_count"]), -i)))
            chosen.append(min(members, key=lambda i: (int(assemblies[i]["contig_count"]), i)))
        chosen = list(dict.fromkeys(chosen))
        while len(chosen) < min(per_lineage, len(members)):
            nxt = max((i for i in members if i not in chosen),
                      key=lambda i: (min(matrix[i][j] for j in chosen), -i))
            chosen.append(nxt)
        selected[lineage] = chosen
    return selected


def read_fasta(path: Path) -> tuple[list[str], dict[str, bytes]]:
    opener = gzip.open if path.suffix == ".gz" else open
    order: list[str] = []
    seqs: dict[str, bytearray] = {}
    current: str | None = None
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in seqs:
                    raise GateError(f"invalid/duplicate FASTA name in {path}")
                order.append(current)
                seqs[current] = bytearray()
            else:
                if current is None:
                    raise GateError(f"FASTA sequence precedes header: {path}")
                seqs[current].extend(line.strip().upper().encode())
    return order, {k: bytes(v) for k, v in seqs.items()}


def reverse_complement(seq: bytes) -> bytes:
    return seq.translate(bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan"))[::-1]


def build_reference_calls(
    paf: Path, reference_fasta: Path, query_fastas: dict[str, Path], samples: Sequence[str],
    min_mapq: int = 20, min_block: int = 1000, min_identity: float = 0.90,
) -> tuple[list[bytearray], dict[str, Any], list[tuple[str, int, int]]]:
    ref_order, ref_seqs = read_fasta(reference_fasta)
    offsets: dict[str, int] = {}
    coordinate: list[tuple[str, int, int]] = []
    total = 0
    for name in ref_order:
        offsets[name] = total
        coordinate.append((name, total, len(ref_seqs[name])))
        total += len(ref_seqs[name])
    sample_index = {s: i for i, s in enumerate(samples)}
    calls = [bytearray(b"N") * total for _ in samples]
    ambiguous = [bytearray(total) for _ in samples]
    query_cache: dict[str, tuple[list[str], dict[str, bytes]]] = {}
    accepted = rejected = 0
    mapped_bases = Counter()
    for line in paf.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            raise GateError("malformed minimap2 PAF")
        qname, qlen, qstart, qend, strand = fields[0], int(fields[1]), int(fields[2]), int(fields[3]), fields[4]
        tname, tstart = fields[5], int(fields[7])
        matches, block, mapq = int(fields[9]), int(fields[10]), int(fields[11])
        sample = qname.split("#", 1)[0]
        cg = next((x[5:] for x in fields[12:] if x.startswith("cg:Z:")), None)
        if (sample not in sample_index or tname not in offsets or cg is None or mapq < min_mapq
                or block < min_block or matches / max(1, block) < min_identity):
            rejected += 1
            continue
        if sample not in query_cache:
            query_cache[sample] = read_fasta(query_fastas[sample])
        qseqs = query_cache[sample][1]
        if qname not in qseqs or len(qseqs[qname]) != qlen:
            raise GateError(f"PAF query/name length mismatch: {qname}")
        qseq = qseqs[qname]
        if strand == "+":
            oriented = qseq
            qpos = qstart
        elif strand == "-":
            oriented = reverse_complement(qseq)
            qpos = qlen - qend
        else:
            raise GateError("invalid PAF strand")
        tpos = offsets[tname] + tstart
        idx = sample_index[sample]
        consumed_q = consumed_t = 0
        cigar_parts = CIGAR_RE.findall(cg)
        if "".join(n + op for n, op in cigar_parts) != cg:
            raise GateError("unsupported/malformed CIGAR")
        for raw_n, op in cigar_parts:
            n = int(raw_n)
            if op in "M=X":
                segment = oriented[qpos:qpos+n]
                if len(segment) != n:
                    raise GateError("CIGAR exceeds query")
                for p, base in enumerate(segment):
                    target = tpos + p
                    if base not in b"ACGT":
                        continue
                    prior = calls[idx][target]
                    if prior == ord("N"):
                        calls[idx][target] = base
                    elif prior != base:
                        ambiguous[idx][target] = 1
                        calls[idx][target] = ord("N")
                qpos += n; tpos += n; consumed_q += n; consumed_t += n
            elif op in "IS":
                qpos += n; consumed_q += n
            elif op in "DN":
                tpos += n; consumed_t += n
            elif op in "HP":
                continue
            else:
                raise GateError(f"unsupported CIGAR operation: {op}")
        mapped_bases[sample] += consumed_t
        accepted += 1
    # The reference itself is exact and callable at every A/C/G/T position.
    reference_sample = ref_order[0].split("#", 1)[0] if ref_order else ""
    if reference_sample in sample_index:
        idx = sample_index[reference_sample]
        for name, offset, length in coordinate:
            seq = ref_seqs[name]
            calls[idx][offset:offset+length] = bytes(base if base in b"ACGT" else ord("N") for base in seq)
    for idx in range(len(samples)):
        for p, flag in enumerate(ambiguous[idx]):
            if flag:
                calls[idx][p] = ord("N")
    stats = {"reference_bases": total, "paf_records_accepted": accepted,
             "paf_records_rejected": rejected, "mapped_target_bases_by_sample": dict(mapped_bases)}
    return calls, stats, coordinate


def core_alignment(calls: Sequence[bytearray], min_sample_fraction: float = 0.95) -> tuple[list[bytes], list[int], dict[str, Any]]:
    if not calls or len({len(x) for x in calls}) != 1:
        raise GateError("invalid reference-call arrays")
    n, length = len(calls), len(calls[0])
    required = math.ceil(n * min_sample_fraction)
    positions = [p for p in range(length) if sum(c[p] in b"ACGT" for c in calls) >= required]
    alignment = [bytes(c[p] for p in positions) for c in calls]
    sample_callable = [sum(base in b"ACGT" for base in c) / max(1, length) for c in calls]
    missing = [sum(base not in b"ACGT" for base in seq) / max(1, len(seq)) for seq in alignment]
    variable = 0
    informative = 0
    for col in zip(*alignment):
        counts = Counter(x for x in col if x in b"ACGT")
        if len(counts) >= 2:
            variable += 1
        if sum(1 for value in counts.values() if value >= 2) >= 2:
            informative += 1
    stats = {"reference_bases": length, "core_callable_sites": len(positions),
             "core_fraction": len(positions) / max(1, length), "required_samples": required,
             "sample_reference_callable_fraction": sample_callable,
             "sample_core_missing_fraction": missing, "mean_core_missing_fraction": sum(missing) / n,
             "variable_sites": variable, "parsimony_informative_sites": informative}
    return alignment, positions, stats


def recombination_candidate_mask(alignment: Sequence[bytes], positions: Sequence[int], z: float = 6.0, window: int = 1000) -> tuple[list[bool], dict[str, Any]]:
    if not alignment:
        raise GateError("empty alignment")
    variable_indices: list[int] = []
    for i, col in enumerate(zip(*alignment)):
        if len({x for x in col if x in b"ACGT"}) >= 2:
            variable_indices.append(i)
    bins: Counter[int] = Counter(positions[i] // window for i in variable_indices)
    values = sorted(bins.values()) or [0]
    median = values[len(values)//2]
    deviations = sorted(abs(x - median) for x in values)
    mad = deviations[len(deviations)//2]
    threshold = max(10, math.floor(median + z * max(1, mad)))
    hot = {b for b, count in bins.items() if count > threshold}
    mask = [(pos // window) in hot for pos in positions]
    stats = {"method": "host-core SNP-density candidate mask; not a Gubbins claim",
             "window_bases": window, "median_variable_sites_per_nonempty_window": median,
             "mad": mad, "z_multiplier": z, "hot_window_threshold": threshold,
             "hot_windows": len(hot), "masked_sites": sum(mask),
             "masked_fraction": sum(mask) / max(1, len(mask)),
             "variable_sites_before": len(variable_indices)}
    return mask, stats


def p_distance_matrix(alignment: Sequence[bytes], keep: Sequence[bool] | None = None) -> list[list[float]]:
    n = len(alignment)
    if n == 0 or len({len(x) for x in alignment}) != 1:
        raise GateError("invalid alignment for p-distance")
    length = len(alignment[0])
    if keep is None:
        keep = [True] * length
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            mismatch = compared = 0
            for p in range(length):
                if not keep[p]:
                    continue
                a, b = alignment[i][p], alignment[j][p]
                if a in b"ACGT" and b in b"ACGT":
                    compared += 1
                    mismatch += a != b
            matrix[i][j] = matrix[j][i] = mismatch / compared if compared else 1.0
    return matrix


def neighbor_joining(labels: Sequence[str], matrix: Sequence[Sequence[float]]) -> TreeNode:
    """Deterministic Saitou-Nei NJ for small high-fidelity lineage trees."""
    active: dict[int, TreeNode] = {i: TreeNode(name=label, children=[]) for i, label in enumerate(labels)}
    dist: dict[tuple[int, int], float] = {}
    for i in active:
        for j in active:
            if i < j:
                dist[i, j] = float(matrix[i][j])
    next_id = len(active)
    while len(active) > 2:
        ids = sorted(active)
        n = len(ids)
        sums = {i: sum(dist[min(i, j), max(i, j)] for j in ids if j != i) for i in ids}
        i, j = min(((i, j) for x, i in enumerate(ids) for j in ids[x+1:]),
                   key=lambda p: ((n - 2) * dist[min(p), max(p)] - sums[p[0]] - sums[p[1]], p))
        dij = dist[min(i, j), max(i, j)]
        li = 0.5 * dij + (sums[i] - sums[j]) / (2 * (n - 2))
        lj = dij - li
        active[i].length = max(0.0, li)
        active[j].length = max(0.0, lj)
        node = TreeNode(children=[active[i], active[j]])
        others = [x for x in ids if x not in (i, j)]
        for x in others:
            dix = dist[min(i, x), max(i, x)]
            djx = dist[min(j, x), max(j, x)]
            dist[min(next_id, x), max(next_id, x)] = max(0.0, 0.5 * (dix + djx - dij))
        active.pop(i); active.pop(j)
        active[next_id] = node
        for key in list(dist):
            if i in key or j in key:
                del dist[key]
        next_id += 1
    i, j = sorted(active)
    d = dist[min(i, j), max(i, j)] if dist else 0.0
    active[i].length = active[j].length = max(0.0, d / 2)
    return TreeNode(children=[active[i], active[j]])


def bootstrap_splits(labels: Sequence[str], alignment: Sequence[bytes], keep: Sequence[bool], replicates: int, seed: int = 104729) -> dict[frozenset[str], float]:
    import random
    usable = [i for i, value in enumerate(keep) if value]
    if not usable:
        raise GateError("no unmasked sites for bootstrap")
    rng = random.Random(seed)
    counts: Counter[frozenset[str]] = Counter()
    for _ in range(replicates):
        sampled = [usable[rng.randrange(len(usable))] for _ in usable]
        sampled_alignment = [bytes(seq[p] for p in sampled) for seq in alignment]
        tree = neighbor_joining(labels, p_distance_matrix(sampled_alignment))
        counts.update(tree_splits(tree))
    return {split: count / replicates for split, count in counts.items()}


def split_concordance(a: TreeNode, b: TreeNode) -> float:
    sa, sb = tree_splits(a), tree_splits(b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


def write_inventory(root: Path) -> None:
    lines = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p.name != "SHA256SUMS"):
        lines.append(f"{sha256_file(path)}  {path.relative_to(root).as_posix()}\n")
    (root / "SHA256SUMS").write_text("".join(lines))
