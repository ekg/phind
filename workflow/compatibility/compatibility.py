#!/usr/bin/env python3
"""Shared fail-closed contracts for the ten-assembly consumer certification.

The module deliberately uses only the Python standard library.  Biological
payloads are read from the checksum-validated predecessor release and temporary
sequence-bearing views are permitted only below the caller-supplied scratch
root.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Iterator, Mapping, Sequence
from urllib.parse import quote_from_bytes, unquote_to_bytes

ROOT_HASHES = {
    "26k_ecoli_accession.txt": "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5",
    "26k_prophage1.csv": "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996",
}
PREDECESSOR_RELEASE_JSON_SHA256 = "4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23"
PREDECESSOR_RELEASE_ID = "canonical-cohort-010-v1-e71484de9994fc28"
EXPECTED_ASSEMBLIES = [
    "GCF_000005845.2", "GCF_000812325.1", "GCF_002302315.1",
    "GCF_004664255.1", "GCF_015644385.1", "GCF_020829045.1",
    "GCF_921380995.1", "GCF_000167895.3", "GCF_001881595.4",
    "GCF_000498835.2",
]
PASS = "PASS"


class GateError(RuntimeError):
    """A fail-closed validation error."""


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while data := fh.read(chunk):
            h.update(data)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json(value))


def append_jsonl(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("ab") as fh:
        fh.write(canonical_json(value))
        fh.flush()
        os.fsync(fh.fileno())


def secure_relpath(raw: str) -> Path:
    p = PurePosixPath(raw)
    if p.is_absolute() or not p.parts or any(x in {"", ".", ".."} for x in p.parts):
        raise GateError(f"unsafe relative path: {raw!r}")
    return Path(*p.parts)


def read_sha256sums(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line_no, line in enumerate(path.read_text().splitlines(), 1):
        if not line:
            continue
        m = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not m:
            raise GateError(f"invalid SHA256SUMS line {line_no} in {path}")
        secure_relpath(m.group(2))
        rows.append((m.group(1), m.group(2)))
    if len({rel for _, rel in rows}) != len(rows):
        raise GateError(f"duplicate SHA256SUMS path in {path}")
    return rows


def verify_inventory(root: Path, inventory: Path, *, exact: bool = False,
                     excluded: Iterable[str] = ("SHA256SUMS", "COMPLETE")) -> int:
    rows = read_sha256sums(inventory)
    for expected, rel in rows:
        p = root / secure_relpath(rel)
        if not p.is_file() or p.is_symlink():
            raise GateError(f"missing/non-regular inventory object: {p}")
        observed = sha256_file(p)
        if observed != expected:
            raise GateError(f"checksum mismatch: {rel}: {observed} != {expected}")
    if exact:
        excluded_set = set(excluded)
        actual = {
            p.relative_to(root).as_posix() for p in root.rglob("*")
            if p.is_file() and not p.is_symlink() and p.relative_to(root).as_posix() not in excluded_set
        }
        listed = {rel for _, rel in rows}
        if actual != listed:
            raise GateError(f"inventory coverage mismatch missing={sorted(actual-listed)} extra={sorted(listed-actual)}")
    return len(rows)


def write_inventory(root: Path, *, excluded: Iterable[str] = ("SHA256SUMS", "COMPLETE")) -> int:
    excluded_set = set(excluded)
    paths = sorted(
        p for p in root.rglob("*")
        if p.is_file() and not p.is_symlink() and p.relative_to(root).as_posix() not in excluded_set
    )
    text = "".join(f"{sha256_file(p)}  {p.relative_to(root).as_posix()}\n" for p in paths)
    (root / "SHA256SUMS").write_text(text)
    return len(paths)


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def verify_root_hashes(repo: Path) -> dict[str, str]:
    observed = {name: sha256_file(repo / name) for name in ROOT_HASHES}
    if observed != ROOT_HASHES:
        raise GateError(f"immutable root digest mismatch: {observed}")
    return observed


@dataclass(frozen=True)
class Assembly:
    order: int
    accession: str
    bgzf: Path
    fai: Path
    gzi: Path
    bgzf_sha256: str
    content_sha256: str
    contig_count: int
    total_bases: int


def _verify_tracked_manifest_dir(repo: Path) -> None:
    root = repo / "manifests/canonical-cohort-010-v1"
    verify_inventory(root, root / "SHA256SUMS", exact=True, excluded=("SHA256SUMS",))


def _bgzf_content_sha(path: Path) -> str:
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        while data := fh.read(1024 * 1024):
            h.update(data)
    return h.hexdigest()


def _fai_names_lengths(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            raise GateError(f"bad FAI row in {path}")
        rows.append((fields[0], int(fields[1])))
    return rows


def verify_predecessor(repo: Path) -> tuple[dict[str, Any], list[Assembly], dict[str, Any]]:
    """Verify every tracked/external predecessor digest and semantic invariant."""
    verify_root_hashes(repo)
    _verify_tracked_manifest_dir(repo)
    release_path = repo / "manifests/canonical-cohort-010-v1/release.json"
    observed_release_sha = sha256_file(release_path)
    if observed_release_sha != PREDECESSOR_RELEASE_JSON_SHA256:
        raise GateError(f"predecessor release.json mismatch: {observed_release_sha}")
    release = json.loads(release_path.read_text())
    if release.get("release_id") != PREDECESSOR_RELEASE_ID or release.get("verdict") != PASS:
        raise GateError("predecessor release identity/verdict is not exact PASS")
    gates = release.get("applicable_gates", {})
    for name, verdict in gates.items():
        if verdict != PASS and not (name == "scale_trend" and verdict == "NOT_APPLICABLE_STAGE_B_NON_SCALE_BEARING"):
            raise GateError(f"predecessor gate {name} is {verdict}")
    counts = release.get("counts", {})
    if counts.get("validated") != 10 or counts.get("distinct_sequence_bearing_assemblies") != 10:
        raise GateError("predecessor is not the exact validated ten-assembly cohort")
    if counts.get("contigs") != 1223 or counts.get("total_bases") != 51731662:
        raise GateError("predecessor contig/base cardinality mismatch")
    if release.get("sequence_bearing_assembly_accessions") != EXPECTED_ASSEMBLIES:
        raise GateError("predecessor cohort order changed")
    if counts.get("global_distinct_assembly_cap") != 1000:
        raise GateError("global distinct-assembly cap changed")

    external = Path(release["external_release_path"])
    if external.name != PREDECESSOR_RELEASE_ID or not (external / "COMPLETE").is_file():
        raise GateError("predecessor external release absent/incomplete")
    expected_external_inventory = (repo / "manifests/canonical-cohort-010-v1/external_SHA256SUMS").read_bytes()
    if (external / "SHA256SUMS").read_bytes() != expected_external_inventory:
        raise GateError("external SHA256SUMS bytes differ from tracked immutable copy")
    external_files = verify_inventory(external, external / "SHA256SUMS")
    external_release = json.loads((external / "release.json").read_text())
    if external_release.get("release_id") != PREDECESSOR_RELEASE_ID or external_release.get("verdict") != PASS:
        raise GateError("external predecessor release identity/verdict mismatch")

    checksums = load_tsv(repo / "manifests/canonical-cohort-010-v1/checksums.tsv")
    if len(checksums) != 80:
        raise GateError(f"expected 80 checksum rows, got {len(checksums)}")
    for row in checksums:
        p = external / secure_relpath(row["relative_path"])
        if not p.is_file() or p.stat().st_size != int(row["bytes"]) or sha256_file(p) != row["sha256"]:
            raise GateError(f"checksum manifest object mismatch: {row['relative_path']}")

    rows = load_tsv(repo / "manifests/canonical-cohort-010-v1/assemblies.tsv")
    if len(rows) != 10 or [r["accession"] for r in rows] != EXPECTED_ASSEMBLIES:
        raise GateError("assembly manifest order/cardinality mismatch")
    assemblies: list[Assembly] = []
    all_names: set[str] = set()
    contigs = bases = 0
    for expected_order, row in enumerate(rows, 1):
        if int(row["stage_b_order"]) != expected_order or row["terminal_state"] != "VALIDATED":
            raise GateError(f"row order/state mismatch for {row['accession']}")
        bgzf = external / secure_relpath(row["canonical_bgzf_relpath"])
        fai = Path(str(bgzf) + ".fai")
        gzi = Path(str(bgzf) + ".gzi")
        checks = [(bgzf, row["canonical_bgzf_sha256"]), (fai, row["fai_sha256"]), (gzi, row["gzi_sha256"])]
        for p, digest in checks:
            if not p.is_file() or sha256_file(p) != digest:
                raise GateError(f"canonical object/index mismatch: {p}")
        if _bgzf_content_sha(bgzf) != row["canonical_fasta_content_sha256"]:
            raise GateError(f"decompressed canonical checksum mismatch: {bgzf}")
        fai_rows = _fai_names_lengths(fai)
        if len(fai_rows) != int(row["contig_count"]):
            raise GateError(f"FAI contig count mismatch: {fai}")
        if sum(n for _, n in fai_rows) != int(row["total_bases"]):
            raise GateError(f"FAI base count mismatch: {fai}")
        for name, _ in fai_rows:
            parts = name.split("#")
            if len(parts) != 3 or parts[0] != row["accession"] or parts[1] != "1" or not parts[2]:
                raise GateError(f"non-canonical PanSN path: {name}")
            if name in all_names:
                raise GateError(f"duplicate PanSN path: {name}")
            all_names.add(name)
        contigs += len(fai_rows)
        bases += sum(n for _, n in fai_rows)
        assemblies.append(Assembly(
            expected_order, row["accession"], bgzf, fai, gzi,
            row["canonical_bgzf_sha256"], row["canonical_fasta_content_sha256"],
            int(row["contig_count"]), int(row["total_bases"]),
        ))
    if contigs != 1223 or bases != 51731662 or len(all_names) != 1223:
        raise GateError("100% contig/base/name accounting failed")
    evidence = {
        "verdict": PASS,
        "release_json_sha256": observed_release_sha,
        "release_id": release["release_id"],
        "external_release": str(external),
        "external_inventory_sha256": sha256_file(external / "SHA256SUMS"),
        "external_inventory_rows": external_files,
        "assemblies": len(assemblies), "contigs": contigs, "bases": bases,
        "cohort_order": EXPECTED_ASSEMBLIES, "global_distinct_assembly_cap": 1000,
    }
    return release, assemblies, evidence


@dataclass(frozen=True)
class ResourceRequest:
    assigned_ram_bytes: int
    durable_allocation_bytes: int
    scratch_allocation_bytes: int
    inode_allocation: int
    predicted_durable_peak_bytes: int
    predicted_scratch_peak_bytes: int
    predicted_files: int
    unfinished_write_bytes: int

    def validate(self) -> None:
        for key, value in asdict(self).items():
            if not isinstance(value, int) or value <= 0:
                raise GateError(f"blank/nonpositive allocation: {key}")
        if self.predicted_durable_peak_bytes > self.durable_allocation_bytes * 0.70:
            raise GateError("predicted durable upper-95 peak exceeds 70% allocation")
        if self.predicted_scratch_peak_bytes > self.scratch_allocation_bytes * 0.70:
            raise GateError("predicted scratch upper-95 peak exceeds 70% allocation")
        if self.predicted_files > self.inode_allocation * 0.50:
            raise GateError("projected files exceed 50% inode allocation")
        if self.durable_allocation_bytes - self.predicted_durable_peak_bytes < 2 * self.unfinished_write_bytes:
            raise GateError("durable allocation lacks 2x unfinished-write reserve")
        if self.scratch_allocation_bytes - self.predicted_scratch_peak_bytes < 2 * self.unfinished_write_bytes:
            raise GateError("scratch allocation lacks 2x unfinished-write reserve")


def _mount_record(path: Path) -> dict[str, Any]:
    proc = subprocess.run(
        ["findmnt", "-T", str(path), "-n", "-o", "TARGET,SOURCE,FSTYPE,OPTIONS"],
        text=True, capture_output=True,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        raise GateError(f"findmnt failed for {path}: {proc.stderr.strip()}")
    parts = proc.stdout.strip().split(None, 3)
    if len(parts) != 4:
        raise GateError(f"incomplete findmnt record for {path}")
    return dict(zip(("target", "source", "fstype", "options"), parts))


def _fs_record(path: Path) -> dict[str, int]:
    st = os.statvfs(path)
    return {
        "free_bytes": st.f_bavail * st.f_frsize,
        "total_bytes": st.f_blocks * st.f_frsize,
        "free_inodes": st.f_favail,
        "total_inodes": st.f_files,
    }


def _swap_bytes() -> tuple[int, int]:
    values: dict[str, int] = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith(("SwapTotal:", "SwapFree:")):
            key, amount, _ = line.split()
            values[key[:-1]] = int(amount) * 1024
    return values.get("SwapTotal", 0), values.get("SwapFree", 0)


def resource_preflight(stage: str, durable_root: Path, scratch_root: Path,
                       request: ResourceRequest, *, write_probe: bool = True) -> dict[str, Any]:
    request.validate()
    durable_root.mkdir(parents=True, exist_ok=True)
    scratch_root.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {"stage": stage, "checked_at_utc": utc_now(), "request": asdict(request)}
    for label, root in (("durable", durable_root), ("scratch", scratch_root)):
        mount = _mount_record(root)
        fs = _fs_record(root)
        if write_probe:
            probe = root / f".compat-write-probe-{os.getpid()}-{time.time_ns()}"
            try:
                fd = os.open(probe, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
                os.write(fd, b"probe\n")
                os.fsync(fd)
                os.close(fd)
            finally:
                probe.unlink(missing_ok=True)
        result[label] = {"path": str(root), "mount": mount, **fs, "owner_uid": root.stat().st_uid,
                         "write_probe": PASS if write_probe else "NOT_RUN"}
    if result["durable"]["free_bytes"] < 2_000_000_000_000 or result["durable"]["free_inodes"] < 1_000_000:
        raise GateError("durable hard floor failed (<2 TB or <1,000,000 inodes)")
    if result["scratch"]["free_bytes"] < 4_000_000_000_000 or result["scratch"]["free_inodes"] < 5_000_000:
        raise GateError("scratch live preflight failed (<4 TB or <5,000,000 inodes)")
    if result["durable"]["free_bytes"] < 2 * request.unfinished_write_bytes:
        raise GateError("durable physical free space lacks 2x unfinished writes")
    if result["scratch"]["free_bytes"] < 2 * request.unfinished_write_bytes:
        raise GateError("scratch physical free space lacks 2x unfinished writes")
    result["swap_total_bytes"], result["swap_free_bytes"] = _swap_bytes()
    result["verdict"] = PASS
    return result


def percent_encode_identifier(raw: bytes) -> str:
    """PanSN field encoding: [A-Za-z0-9._-] passes, all else uppercase %HH."""
    safe = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-"
    return "".join(chr(b) if b in safe else f"%{b:02X}" for b in raw)


def percent_decode_identifier(encoded: str) -> bytes:
    escapes = re.findall(r"%(..)", encoded)
    if re.search(r"%(?![0-9A-F]{2})", encoded) or any(
        not re.fullmatch(r"[0-9A-F]{2}", escape) for escape in escapes
    ):
        raise GateError(f"noncanonical percent escape: {encoded}")
    decoded = unquote_to_bytes(encoded)
    if percent_encode_identifier(decoded) != encoded:
        raise GateError(f"identifier is not canonical one-layer encoding: {encoded}")
    return decoded


def gff_lexical_encode(semantic: str) -> str:
    """Encode a semantic GFF seqid; '%' and '#' necessarily become %25/%23."""
    return quote_from_bytes(semantic.encode("utf-8"), safe="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:^*$@!+_?-|")


def gff_lexical_decode(lexical: str) -> str:
    try:
        value = unquote_to_bytes(lexical).decode("utf-8")
    except (UnicodeDecodeError, ValueError) as exc:
        raise GateError(f"invalid GFF lexical seqid: {lexical}") from exc
    if gff_lexical_encode(value) != lexical:
        raise GateError(f"noncanonical GFF lexical seqid: {lexical}")
    return value


def stage_gff_semantic_alias(source: Path, destination: Path, allowed_semantic_ids: set[str],
                             mapping_path: Path) -> list[dict[str, str]]:
    """Make a raw-ID parser view by strict single-layer GFF seqid decoding.

    The view is reversible because every lexical/semantic pair is recorded and
    duplicate semantic aliases are rejected.  Coordinates and all other bytes
    remain unchanged.
    """
    mapping: dict[str, str] = {}
    output: list[str] = []
    for line_no, line in enumerate(source.read_text().splitlines(keepends=True), 1):
        if not line or line.startswith("#"):
            output.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise GateError(f"GFF row {line_no} does not have 9 columns")
        lexical = fields[0]
        semantic = gff_lexical_decode(lexical)
        if semantic not in allowed_semantic_ids:
            raise GateError(f"GFF semantic seqid is not in FASTA: {semantic}")
        previous = mapping.get(semantic)
        if previous is not None and previous != lexical:
            raise GateError(f"ambiguous GFF lexical aliases for {semantic}")
        mapping[semantic] = lexical
        fields[0] = semantic
        output.append("\t".join(fields) + ("\n" if line.endswith("\n") else ""))
    destination.write_text("".join(output))
    rows = [{"lexical": lexical, "semantic": semantic} for semantic, lexical in sorted(mapping.items())]
    write_json(mapping_path, {"schema_version": "gff-semantic-alias-v1", "rows": rows, "reversible": True})
    return rows


def mash_lower_to_full_phylip(source: Path, destination: Path) -> dict[str, Any]:
    """Validate Mash relaxed lower PHYLIP and materialize RapidNJ full PHYLIP.

    Mash 2.3 emits row i with exactly i off-diagonal values and no diagonal.
    RapidNJ 2.3.2's ``pd`` reader requires N values on every row.  Labels are
    copied byte-for-byte as whitespace-free tokens and the matrix is expanded
    symmetrically with an exact zero diagonal.
    """
    lines = source.read_text().splitlines()
    if not lines:
        raise GateError("empty Mash lower PHYLIP")
    try:
        n = int(lines[0].strip())
    except ValueError as exc:
        raise GateError("invalid Mash lower PHYLIP cardinality") from exc
    if n <= 0 or len(lines[1:]) != n:
        raise GateError("Mash lower PHYLIP row cardinality mismatch")
    names: list[str] = []
    lower: list[list[float]] = []
    for i, line in enumerate(lines[1:]):
        fields = line.split()
        if len(fields) != i + 1:
            raise GateError(f"Mash lower PHYLIP row {i} has {len(fields)-1} values, expected {i}")
        name = fields[0]
        if not name or any(c.isspace() for c in name) or name in names:
            raise GateError(f"invalid/duplicate Mash label: {name!r}")
        try:
            values = [float(x) for x in fields[1:]]
        except ValueError as exc:
            raise GateError(f"nonnumeric Mash distance on row {i}") from exc
        if any(not (0.0 <= value <= 1.0) for value in values):
            raise GateError(f"out-of-range Mash distance on row {i}")
        names.append(name)
        lower.append(values)
    matrix = [[0.0] * n for _ in range(n)]
    for i, values in enumerate(lower):
        for j, value in enumerate(values):
            matrix[i][j] = matrix[j][i] = value
    with destination.open("w") as fh:
        fh.write(f"{n}\n")
        for name, row in zip(names, matrix):
            fh.write(name + "\t" + "\t".join(f"{value:.9g}" for value in row) + "\n")
    return {"rows": n, "off_diagonal_pairs": n * (n - 1) // 2,
            "full_values": n * n, "labels": names, "symmetric": True,
            "diagonal_zero": True}


def parse_fasta(path: Path, *, gzipped: bool | None = None) -> Iterator[tuple[str, bytes]]:
    if gzipped is None:
        gzipped = path.suffix == ".gz"
    opener = gzip.open if gzipped else open
    with opener(path, "rb") as fh:  # type: ignore[arg-type]
        name: str | None = None
        seq = bytearray()
        for raw in fh:
            if raw.startswith(b">"):
                if name is not None:
                    yield name, bytes(seq)
                token = raw[1:].strip().split(None, 1)[0]
                name = token.decode("utf-8")
                seq.clear()
            else:
                if name is None:
                    raise GateError(f"sequence before FASTA header in {path}")
                seq.extend(raw.strip().upper())
        if name is not None:
            yield name, bytes(seq)


def write_fasta(path: Path, records: Iterable[tuple[str, bytes]], width: int = 60) -> None:
    with path.open("wb") as fh:
        for name, seq in records:
            if not name or any(c.isspace() for c in name):
                raise GateError(f"unsafe FASTA primary token: {name!r}")
            fh.write(b">" + name.encode() + b"\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + b"\n")


def tree_size(root: Path) -> tuple[int, int]:
    total = files = 0
    if not root.exists():
        return 0, 0
    for p in root.rglob("*"):
        if p.is_file() and not p.is_symlink():
            files += 1
            total += p.stat().st_size
    return total, files


def fsync_tree(root: Path) -> None:
    for p in sorted(root.rglob("*"), reverse=True):
        if p.is_file():
            with p.open("rb") as fh:
                os.fsync(fh.fileno())
        elif p.is_dir():
            fd = os.open(p, os.O_RDONLY | os.O_DIRECTORY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_promote(staging: Path, final: Path) -> None:
    if final.exists():
        raise GateError(f"refusing overwrite of release: {final}")
    if staging.stat().st_dev != final.parent.stat().st_dev:
        raise GateError("release promotion is not same-filesystem")
    if not (staging / "COMPLETE").is_file():
        raise GateError("refusing promotion without COMPLETE")
    fsync_tree(staging)
    os.rename(staging, final)
    fd = os.open(final.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def require_cleaned(path: Path) -> None:
    if path.exists():
        raise GateError(f"required scratch cleanup did not occur: {path}")
