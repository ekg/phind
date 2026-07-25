#!/usr/bin/env python3
"""Independent, read-only audit of canonical-cohort-250-v1.

This module intentionally does not import the cohort runner or its validator.
It re-derives canonical-only gates from immutable bytes and local tools.  It
never opens a network connection and never writes below a canonical release.
"""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import math
import os
import re
import stat
import struct
import subprocess
import sys
import urllib.parse
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any, BinaryIO, Iterable

RELEASE_ID = "canonical-cohort-250-v1-a6184d7d6ee08bda"
EXTERNAL_RELEASE = Path(
    "/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-250/"
    + RELEASE_ID
)
N100_ID = "canonical-cohort-100-v1-6be4c0dde65f31d0"
N10_ID = "canonical-cohort-010-v1-e71484de9994fc28"
SELECTION_ID = "pilot-cohorts-v1-8afc0ea03d9e50dc"
COMPATIBILITY_ID = "consumer-compatibility-v1-78d7e93f19fa3d87"
SEMANTICS_ID = "prophage-semantics-v1-f5619e221ff272ae"
ROOT_HASHES = {
    "26k_ecoli_accession.txt": "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5",
    "26k_prophage1.csv": "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996",
}
PINS = {
    "release_json": "dcf2b887afa51e4e0e739ae2fef9b5a9d72fb8bc9a4d698a161a99673aaf504a",
    "external_sha256sums": "45fd42b76bf1c1ace3a2e882fe6a9a8f6af2457c0f5d4bc28011cb99b521c5b7",
    "complete": "026ba58d2865915284df8e32abb05ecb1ace862650f5f784144e751dc4b9bce0",
    "external_tree": "5a11927a00486d327d0a33f42ec6e4361b1716118a4e623e1da99401038cad79",
    "tracked_tree": "37e595ef8d745f0b0c9074f8f9a88aab73b117b671e0fde35e97dff37ff6853d",
    "tracked_sha256sums": "4100e202fe4db19dc6e9ea3cf79fac53c00badc0fafa814ad5b2141d38b0f696",
    "cohort_250": "ba2cf2909ccf62a0c1944a76b522edc5600953511ec355479117b4a419acbc9f",
    "cohort_100": "13e203961a9fcec18a8a09e690582652d8085b2a386811e6c6a03184b9489182",
    "n100_release": "3b91b24e23323ef971a13f22825e512a233bb592ed641ea9b270a2f1fd683795",
    "n10_release": "4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23",
    "selection_release": "d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99",
    "compatibility_release": "021719ddadd7bb7fa2932d2ef9cb25da9c666ebe0389988691283011ee12f4c7",
}
ACCESSION_RE = re.compile(r"^GC[AF]_[0-9]{9}\.[1-9][0-9]*$")
SEQ_ACCESSION_RE = re.compile(r"^([A-Za-z]{1,6}_?[0-9]+)\.([1-9][0-9]*)$")
DNA = frozenset(b"ACGTRYSWKMBDHVN")
SAFE_CONTIG = frozenset(b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-")
GFF_SAFE = frozenset(b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:^*$@!+_?-|")
ROLE_NAMES = {
    "source_package", "source_manifest", "canonical_bgzf", "fai", "gzi",
    "contig_crosswalk", "annotation_aliases", "canonical_manifest",
}


class AuditError(RuntimeError):
    """A fail-closed canonical audit error."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_snapshot(root: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    files = 0
    size = 0
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file() and not candidate.is_symlink()):
        relative = str(path.relative_to(root))
        file_digest = sha_file(path)
        digest.update(relative.encode() + b"\0" + file_digest.encode() + b"\n")
        files += 1
        size += path.stat().st_size
    return {"root": str(root), "files": files, "bytes": size, "tree_sha256": digest.hexdigest()}


def parse_inventory(root: Path, *, exact: bool = True, complete: bool = True) -> dict[str, str]:
    sums_path = root / "SHA256SUMS"
    require(sums_path.is_file() and not sums_path.is_symlink(), f"missing/unsafe inventory: {sums_path}")
    rows: dict[str, str] = {}
    for number, line in enumerate(sums_path.read_text().splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, f"malformed inventory line: {sums_path}:{number}")
        digest, relative = match.groups()
        require(relative not in rows and not relative.startswith("/") and ".." not in Path(relative).parts,
                f"unsafe/duplicate inventory path: {root}:{relative}")
        path = root / relative
        require(path.is_file() and not path.is_symlink(), f"inventory target missing/symlink: {path}")
        require(sha_file(path) == digest, f"inventory digest mismatch: {path}")
        rows[relative] = digest
    if exact:
        actual = {
            str(path.relative_to(root))
            for path in root.rglob("*")
            if path.is_file() and str(path.relative_to(root)) not in {"SHA256SUMS", "COMPLETE"}
        }
        require(actual == set(rows), f"inventory coverage mismatch: {root}")
    if complete:
        marker = root / "COMPLETE"
        require(marker.is_file() and not marker.is_symlink(), f"missing/unsafe COMPLETE: {root}")
        text = marker.read_text().strip()
        expected = sha_file(sums_path)
        if text.startswith("{"):
            obj = json.loads(text)
            require(obj.get("sha256sums_sha256") == expected and obj.get("verdict") == "PASS",
                    f"COMPLETE JSON mismatch: {root}")
        else:
            require(text.split()[:1] == [expected], f"COMPLETE inventory pin mismatch: {root}")
    require(not any(path.is_symlink() for path in root.rglob("*")), f"symlink below immutable root: {root}")
    return rows


def row_hash(row: dict[str, Any], fields: list[str]) -> str:
    payload = "\t".join(str(row.get(field, ".")) for field in fields if field != "row_sha256") + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def read_tsv(path: Path, *, compressed: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    require(fields and fields[-1] == "row_sha256", f"row-hash schema missing: {path}")
    for number, row in enumerate(rows, 2):
        require(row.get("row_sha256") == row_hash(row, fields), f"row digest mismatch: {path}:{number}")
    return fields, rows


def read_plain_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    require(fields, f"empty TSV schema: {path}")
    return fields, rows


def require_release(repo: Path, tracked_dir: Path, expected_id: str, expected_sha: str) -> tuple[dict[str, Any], Path]:
    parse_inventory(tracked_dir, complete=False)
    release_path = tracked_dir / "release.json"
    require(sha_file(release_path) == expected_sha, f"tracked release pin mismatch: {expected_id}")
    release = json.loads(release_path.read_text())
    require(release.get("release_id") == expected_id and release.get("verdict") == "PASS" and release.get("immutable") is True,
            f"release identity/verdict/immutability mismatch: {expected_id}")
    external_raw = release.get("external_release_path", release.get("external_path"))
    require(isinstance(external_raw, str) and external_raw, f"release lacks external path: {expected_id}")
    external = Path(external_raw)
    parse_inventory(external)
    require((external / "release.json").read_bytes() == release_path.read_bytes(), f"tracked/external release mismatch: {expected_id}")
    return release, external


def nested_rows(larger: list[dict[str, str]], smaller: list[dict[str, str]]) -> bool:
    # Rung size and design weights are intentionally recomputed per nested rung;
    # every identity, order, stratum, occurrence, and frozen collection digest
    # must remain identical.
    ignored = {"rung_n", "inclusion_probability", "inference_weight", "row_sha256"}
    for left, right in zip(larger, smaller):
        keys = set(left) | set(right)
        if any(left.get(key) != right.get(key) for key in keys - ignored):
            return False
    return len(larger) >= len(smaller)


def all_strings(value: Any, key_hint: str = "") -> Iterable[tuple[str, str]]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield from all_strings(child, str(key))
    elif isinstance(value, list):
        for child in value:
            yield from all_strings(child, key_hint)
    elif isinstance(value, str):
        yield key_hint, value


def catalog_accessions(value: Any) -> set[str]:
    return {
        text for key, text in all_strings(value)
        if key.lower() in {"accession", "assemblyaccession", "assembly_accession"}
        and ACCESSION_RE.fullmatch(text)
    }


def safe_zip_path(raw: str) -> str:
    result = raw[2:] if raw.startswith("./") else raw
    require(not result.startswith("/") and ".." not in Path(result).parts and "\\" not in result,
            f"unsafe ZIP path: {raw}")
    return result


def hash_member(archive: zipfile.ZipFile, member: str) -> tuple[str, str, int]:
    md5 = hashlib.md5(usedforsecurity=False)
    sha = hashlib.sha256()
    size = 0
    with archive.open(member) as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(block)
            sha.update(block)
            size += len(block)
    return md5.hexdigest(), sha.hexdigest(), size


def validate_archive(package: Path, accession: str) -> tuple[dict[str, Any], zipfile.ZipFile]:
    require(zipfile.is_zipfile(package), f"not a ZIP archive: {accession}")
    archive = zipfile.ZipFile(package)
    infos = archive.infolist()
    require(0 < len(infos) <= 10_000, f"invalid ZIP entry count: {accession}")
    names: list[str] = []
    total = 0
    for info in infos:
        name = safe_zip_path(info.filename)
        names.append(name)
        total += info.file_size
        require(not stat.S_ISLNK(info.external_attr >> 16), f"symlink in archive: {accession}:{name}")
    require(total <= 2_000_000_000, f"unbounded ZIP expansion: {accession}")
    require(archive.testzip() is None, f"ZIP CRC failure: {accession}")
    files = {name for name in names if not name.endswith("/")}
    md5_members = [name for name in files if Path(name).name == "md5sum.txt"]
    require(len(md5_members) == 1, f"md5sum cardinality mismatch: {accession}")
    md5_member = md5_members[0]
    md5_rows: dict[str, str] = {}
    for number, line in enumerate(archive.read(md5_member).decode().splitlines(), 1):
        match = re.fullmatch(r"([0-9a-fA-F]{32})\s+\*?(.+)", line)
        require(match is not None, f"malformed upstream MD5: {accession}:{number}")
        relative = safe_zip_path(match.group(2))
        require(relative not in md5_rows, f"duplicate upstream MD5 path: {accession}:{relative}")
        md5_rows[relative] = match.group(1).lower()
    payload = files - {md5_member}
    require(set(md5_rows) - payload == set() and (payload - set(md5_rows)).issubset({"README.md"}),
            f"upstream MD5 coverage mismatch: {accession}")
    member_stats: dict[str, tuple[str, int]] = {}
    for relative, expected in md5_rows.items():
        md5, sha, size = hash_member(archive, relative)
        require(md5 == expected, f"upstream MD5 mismatch: {accession}:{relative}")
        member_stats[relative] = (sha, size)
    assembly_dirs = {
        parts[2] for name in payload
        if len(parts := Path(name).parts) >= 4 and parts[:2] == ("ncbi_dataset", "data")
        and ACCESSION_RE.fullmatch(parts[2])
    }
    require(assembly_dirs == {accession}, f"archive exact accession mismatch: {accession}:{sorted(assembly_dirs)}")
    fasta = [name for name in payload if name.startswith(f"ncbi_dataset/data/{accession}/") and name.endswith("_genomic.fna")]
    gff = [name for name in payload if name.startswith(f"ncbi_dataset/data/{accession}/") and name.endswith(("genomic.gff", "_genomic.gff"))]
    sequence_report = [name for name in payload if name == f"ncbi_dataset/data/{accession}/sequence_report.jsonl"]
    catalogs = [name for name in payload if Path(name).name == "dataset_catalog.json"]
    reports = [name for name in payload if Path(name).name == "assembly_data_report.jsonl"]
    require(len(fasta) == 1 and len(gff) <= 1 and len(sequence_report) <= 1 and len(catalogs) == 1,
            f"archive required member cardinality mismatch: {accession}")
    catalog_ids = catalog_accessions(json.loads(archive.read(catalogs[0])))
    require(not catalog_ids or catalog_ids == {accession}, f"catalog accession mismatch: {accession}")
    report_ids: set[str] = set()
    for report in reports:
        for line in archive.read(report).splitlines():
            if line.strip():
                report_ids.update(catalog_accessions(json.loads(line)))
    require(not report_ids or accession in report_ids, f"assembly report accession mismatch: {accession}")
    fasta_sha, fasta_bytes = member_stats[fasta[0]]
    gff_name = gff[0] if gff else "."
    gff_sha, gff_bytes = member_stats[gff_name] if gff else (".", 0)
    seq_name = sequence_report[0] if sequence_report else "."
    seq_sha = member_stats[seq_name][0] if sequence_report else "."
    return {
        "package_sha256": sha_file(package), "package_bytes": package.stat().st_size,
        "zip_entries": len(infos), "uncompressed_member_bytes": total, "md5_entries": len(md5_rows),
        "fasta_member": fasta[0], "fasta_sha256": fasta_sha, "fasta_bytes": fasta_bytes,
        "gff_member": gff_name, "gff_sha256": gff_sha, "gff_bytes": gff_bytes,
        "sequence_report_member": seq_name, "sequence_report_sha256": seq_sha,
        "dataset_catalog_member": catalogs[0], "dataset_catalog_sha256": member_stats[catalogs[0]][0],
        "assembly_report_accessions": sorted(report_ids),
    }, archive


def encode_contig(token: bytes, sample: str) -> tuple[str, str]:
    require(token, "empty FASTA token")
    encoded = "".join(chr(byte) if byte in SAFE_CONTIG else f"%{byte:02X}" for byte in token)
    encoding = "IDENTITY_V1" if encoded.encode() == token else "PERCENT_UTF8_BYTES_V1"
    if len(f"{sample}#1#{encoded}".encode()) > 240:
        encoded = "CTGSHA256_" + hashlib.sha256(token).hexdigest()
        encoding = "SHA256_ALIAS_V1"
    require(len(f"{sample}#1#{encoded}".encode()) <= 240, "PanSN name exceeds 240 bytes")
    return encoded, encoding


def parse_source_fasta(data: bytes, accession: str) -> list[dict[str, Any]]:
    require(b"\r" not in data, f"CR byte in source FASTA: {accession}")
    records: list[dict[str, Any]] = []
    header: bytes | None = None
    token: bytes | None = None
    digest: hashlib._Hash | None = None
    length = 0
    prefix = bytearray()

    def finish() -> None:
        nonlocal header, token, digest, length, prefix
        if header is None:
            return
        require(token is not None and digest is not None and length > 0, f"empty source FASTA record: {accession}")
        encoded, encoding = encode_contig(token, accession)
        records.append({
            "header": header, "token": token, "display": token.decode("utf-8", "replace"),
            "encoding": encoding, "pansn_contig": encoded, "name": f"{accession}#1#{encoded}",
            "length": length, "sha256": digest.hexdigest(), "prefix": bytes(prefix),
        })
        header = token = digest = None
        length = 0
        prefix = bytearray()

    seen: set[bytes] = set()
    for raw in data.splitlines():
        if raw.startswith(b">"):
            finish()
            header = raw[1:]
            parts = header.split(None, 1)
            require(parts, f"empty source FASTA header: {accession}")
            token = parts[0]
            require(token not in seen, f"duplicate source FASTA token: {accession}")
            seen.add(token)
            digest = hashlib.sha256()
        else:
            require(header is not None and digest is not None and raw and all(byte in DNA for byte in raw),
                    f"invalid source FASTA sequence: {accession}")
            digest.update(raw)
            length += len(raw)
            if len(prefix) < 60:
                prefix.extend(raw[:60 - len(prefix)])
    finish()
    require(records, f"no source FASTA records: {accession}")
    return records


def parse_canonical_bgzf(path: Path) -> tuple[list[dict[str, Any]], str]:
    content = hashlib.sha256()
    records: list[dict[str, Any]] = []
    current: str | None = None
    digest: hashlib._Hash | None = None
    length = 0
    prefix = bytearray()

    def finish() -> None:
        nonlocal current, digest, length, prefix
        if current is not None:
            require(digest is not None and length > 0, f"empty canonical FASTA record: {path}")
            records.append({"name": current, "length": length, "sha256": digest.hexdigest(), "prefix": bytes(prefix)})
        current = digest = None
        length = 0
        prefix = bytearray()

    seen: set[str] = set()
    with gzip.open(path, "rb") as handle:
        for raw in handle:
            content.update(raw)
            require(b"\r" not in raw and raw.endswith(b"\n"), f"non-LF canonical FASTA: {path}")
            line = raw[:-1]
            if line.startswith(b">"):
                finish()
                current = line[1:].decode("ascii")
                fields = current.split("#")
                require(len(fields) == 3 and ACCESSION_RE.fullmatch(fields[0]) and fields[1] == "1" and fields[2],
                        f"invalid PanSN grammar: {current}")
                require(current not in seen and len(current.encode()) <= 240, f"duplicate/overlength PanSN name: {current}")
                seen.add(current)
                digest = hashlib.sha256()
            else:
                require(current is not None and digest is not None and 0 < len(line) <= 60 and all(byte in DNA for byte in line),
                        f"invalid canonical sequence line: {path}")
                digest.update(line)
                length += len(line)
                if len(prefix) < 60:
                    prefix.extend(line[:60 - len(prefix)])
    finish()
    return records, content.hexdigest()


def validate_fai(path: Path, records: list[dict[str, Any]]) -> None:
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        columns = line.split("\t")
        require(len(columns) >= 5, f"malformed FAI: {path}:{number}")
        rows.append((columns[0], int(columns[1]), int(columns[2]), int(columns[3]), int(columns[4])))
    require([(row[0], row[1]) for row in rows] == [(row["name"], row["length"]) for row in records],
            f"FAI names/length/order mismatch: {path}")
    require(all(row[2] >= 0 and row[3] > 0 and row[4] >= row[3] for row in rows), f"invalid FAI offsets/layout: {path}")


def validate_gzi(path: Path, compressed_bytes: int) -> None:
    data = path.read_bytes()
    require(len(data) >= 8 and (len(data) - 8) % 16 == 0, f"malformed GZI: {path}")
    count = struct.unpack("<Q", data[:8])[0]
    require(len(data) == 8 + 16 * count, f"GZI count mismatch: {path}")
    pairs = [struct.unpack("<QQ", data[8 + 16*i:24 + 16*i]) for i in range(count)]
    require(all(0 < compressed < compressed_bytes and uncompressed > 0 for compressed, uncompressed in pairs),
            f"GZI offsets out of range: {path}")
    require(pairs == sorted(pairs), f"non-monotonic GZI offsets: {path}")


def validate_faidx_regions(samtools: str, bgzf: Path, records: list[dict[str, Any]]) -> int:
    checked = 0
    for offset in range(0, len(records), 500):
        chunk = records[offset:offset + 500]
        regions = [f"{row['name']}:1-{min(60, row['length'])}" for row in chunk]
        result = subprocess.run([samtools, "faidx", str(bgzf), *regions], capture_output=True)
        require(result.returncode == 0, f"samtools faidx failed: {bgzf}: {result.stderr.decode(errors='replace')}")
        observed: list[bytes] = []
        sequence = bytearray()
        for line in result.stdout.splitlines():
            if line.startswith(b">"):
                if sequence:
                    observed.append(bytes(sequence))
                    sequence.clear()
            else:
                sequence.extend(line)
        if sequence:
            observed.append(bytes(sequence))
        expected = [row["prefix"][:min(60, row["length"])] for row in chunk]
        require(observed == expected, f"samtools indexed region mismatch: {bgzf}")
        checked += len(chunk)
    return checked


def load_sequence_reports(archive: zipfile.ZipFile, member: str) -> list[dict[str, Any]]:
    if member == ".":
        return []
    return [json.loads(line) for line in archive.read(member).splitlines() if line.strip()]


def report_for_token(reports: list[dict[str, Any]], token: str) -> dict[str, Any] | None:
    matches = [report for report in reports if token in {value for _, value in all_strings(report)}]
    require(len(matches) <= 1, f"ambiguous sequence report token: {token}")
    return matches[0] if matches else None


def report_aliases(report: dict[str, Any] | None) -> set[str]:
    if report is None:
        return set()
    return {
        value for key, value in all_strings(report)
        if key.lower().replace("_", "") in {"refseqaccession", "genbankaccession", "sequencename"}
        and SEQ_ACCESSION_RE.fullmatch(value)
    }


def gff_escape(value: str) -> str:
    return "".join(chr(byte) if byte in GFF_SAFE and byte != ord("%") else f"%{byte:02X}" for byte in value.encode())


def recompute_annotation(
    archive: zipfile.ZipFile, gff_member: str, source_records: list[dict[str, Any]],
    reports: list[dict[str, Any]], accession: str,
) -> tuple[list[dict[str, Any]], dict[str, Any], str | None]:
    if gff_member == ".":
        return [], {"status": "NOT_AVAILABLE", "coordinate_policy": "NOT_APPLICABLE_NO_GFF", "feature_rows": 0}, None
    aliases: dict[str, dict[str, Any]] = {}
    try:
        for record in source_records:
            report = report_for_token(reports, record["display"])
            for candidate in {record["display"], *report_aliases(report)}:
                prior = aliases.get(candidate)
                if prior is not None and prior["name"] != record["name"]:
                    raise AuditError(f"sequence-report alias collision: {accession}:{candidate}")
                aliases[candidate] = record
        used: dict[str, dict[str, Any]] = {}
        feature_rows = 0
        with archive.open(gff_member) as handle:
            for line_number, raw in enumerate(handle, 1):
                if b"\r" in raw:
                    raise AuditError(f"CR byte in GFF: {accession}:{line_number}")
                if raw.startswith(b"##FASTA"):
                    raise AuditError(f"unexpected embedded FASTA in source GFF: {accession}")
                if not raw or raw.startswith(b"#") or not raw.strip():
                    continue
                columns = raw.rstrip(b"\n").split(b"\t")
                if len(columns) != 9:
                    raise AuditError(f"GFF row does not have 9 columns: {accession}:{line_number}")
                lexical = columns[0].decode("ascii")
                decoded_bytes = urllib.parse.unquote_to_bytes(lexical)
                decoded = decoded_bytes.decode("utf-8")
                record = aliases.get(decoded)
                if record is None:
                    raise AuditError(f"GFF seqid does not resolve to FASTA/sequence-report alias: {accession}:{decoded}")
                try:
                    start, end = int(columns[3]), int(columns[4])
                except ValueError as exc:
                    raise AuditError(f"non-integer GFF coordinate: {accession}:{line_number}") from exc
                if not (1 <= start <= end <= record["length"]):
                    raise AuditError(f"out-of-range GFF coordinate: {accession}:{line_number}:{start}-{end}")
                feature_rows += 1
                used.setdefault(lexical, {
                    "accession": accession, "source_gff_member": gff_member,
                    "source_gff_seqid_lexical_b64": base64.b64encode(columns[0]).decode(),
                    "source_gff_seqid_decoded_b64": base64.b64encode(decoded_bytes).decode(),
                    "source_fasta_id_token_b64": base64.b64encode(record["token"]).decode(),
                    "pansn_sequence_name": record["name"],
                    "canonical_gff_seqid_lexical": gff_escape(record["name"]),
                    "canonical_gff_seqid_decoded": record["name"],
                    "source_coordinate_convention": "GFF3_1_BASED_CLOSED",
                })
        if feature_rows == 0:
            raise AuditError(f"present GFF contains no feature rows: {accession}")
        rows = [used[key] | {"alias_order": number} for number, key in enumerate(sorted(used), 1)]
        summary = {"status": "ALIASES_VALIDATED_NO_TRANSFORMED_GFF", "coordinate_policy": "GFF3_1_BASED_CLOSED_PASS",
                   "feature_rows": feature_rows, "distinct_seqids": len(rows)}
        return rows, summary, None
    except (AuditError, UnicodeDecodeError, ValueError) as exc:
        return [], {}, str(exc)


def compare_alias_rows(observed: list[dict[str, str]], expected: list[dict[str, Any]], fields: list[str]) -> bool:
    comparable = [field for field in fields if field != "row_sha256"]
    normalized = [{field: str(row.get(field, ".")) for field in comparable} for row in expected]
    return [{field: row[field] for field in comparable} for row in observed] == normalized


def graph_contract(repo: Path) -> dict[str, Any]:
    def wg_json(*args: str) -> dict[str, Any]:
        result = subprocess.run(["wg", "--json", *args], cwd=repo, capture_output=True, text=True)
        require(result.returncode == 0, f"WG query failed: {' '.join(args)}: {result.stderr}")
        return json.loads(result.stdout)

    prepare = wg_json("show", "prepare-canonical-cohort-250")
    stale = "integrated N=100 automatic `GO_250` verdict"
    require(stale in prepare.get("description", ""), "stale integrated-GO clause was not found verbatim")
    source = wg_json("show", "resolve-prophage-source")
    require(source.get("status") == "failed" and "EXTRACTION_BLOCKED" in source.get("failure_reason", ""),
            "resolve-prophage-source is not failed/EXTRACTION_BLOCKED")
    integrated: dict[str, Any] = {}
    for task_id in ("run-integrated-100-genome", "run-integrated-250-genome", "run-integrated-500-genome", "run-integrated-syng"):
        blocked = wg_json("why-blocked", task_id)
        encoded = json.dumps(blocked, sort_keys=True)
        require(blocked.get("is_blocked") is True and "resolve-prophage-source" in encoded and "EXTRACTION_BLOCKED" in encoded,
                f"integrated task is not blocked by source semantics: {task_id}")
        integrated[task_id] = {"task_status": blocked["task"]["status"], "is_blocked": True,
                               "blocked_by_failed_source_semantics": True}
    semantics_ref = json.loads((repo / "artifacts/prophage_semantics/release_reference.json").read_text())
    consumer = json.loads((repo / "artifacts/prophage_semantics/consumer_gate.json").read_text())
    parse_inventory(repo / "artifacts/prophage_semantics", exact=True, complete=False)
    require(semantics_ref.get("release_id") == SEMANTICS_ID and semantics_ref.get("verdict") == "EXTRACTION_BLOCKED"
            and semantics_ref.get("consumer_action") == "REJECT", "semantics reference mismatch")
    sem_external = Path(semantics_ref["external_path"])
    parse_inventory(sem_external, complete=False)
    sem_complete = json.loads((sem_external / "COMPLETE").read_text())
    require(sha_file(sem_external / "SHA256SUMS") == semantics_ref.get("sha256sums_sha256")
            and sha_file(sem_external / "COMPLETE") == semantics_ref.get("complete_sha256")
            and sem_complete.get("sha256sums_sha256") == semantics_ref.get("sha256sums_sha256")
            and sem_complete.get("verdict") == "EXTRACTION_BLOCKED",
            "blocked semantics COMPLETE/inventory pin mismatch")
    sem_release = json.loads((sem_external / "release.json").read_text())
    require(sem_release.get("release_id") == SEMANTICS_ID and sem_release.get("verdict") == "EXTRACTION_BLOCKED",
            "external semantics release mismatch")
    require(consumer.get("release_id") == SEMANTICS_ID and consumer.get("release_verdict") == "EXTRACTION_BLOCKED"
            and consumer.get("consumer_action") == "REJECT" and consumer.get("strict_validator_exit_code") == 2,
            "strict extraction consumer did not reject")
    return {
        "stale_clause": stale,
        "stale_clause_location": "WG task prepare-canonical-cohort-250 / Automatic dispatch and immutable inputs",
        "disposition": "SUPERSEDED_FOR_CANONICAL_BRANCH_WITHOUT_HISTORY_REWRITE",
        "corrected_contract": "CANONICAL_ONLY_ACQUISITION_AND_CANONICALIZATION_GATES",
        "source_task": {"task_id": "resolve-prophage-source", "status": "failed",
                        "extraction_eligibility": "EXTRACTION_BLOCKED", "consumer_action": "REJECT",
                        "release_id": SEMANTICS_ID},
        "integrated_tasks": integrated,
    }


def audit(repo: Path, external: Path, samtools: str, bgzip: str) -> dict[str, Any]:
    repo = repo.resolve()
    external = external.resolve()
    require(external == EXTERNAL_RELEASE, f"unexpected release path: {external}")
    roots_start = {name: sha_file(repo / name) for name in ROOT_HASHES}
    require(roots_start == ROOT_HASHES, "immutable root input mismatch at audit start")
    external_start = tree_snapshot(external)
    tracked_root = repo / "manifests/canonical-cohort-250-v1"
    tracked_start = tree_snapshot(tracked_root)
    require(external_start["tree_sha256"] == PINS["external_tree"] and external_start["files"] == 2269
            and external_start["bytes"] == 581_575_254, "immutable external release snapshot mismatch")
    require(tracked_start["tree_sha256"] == PINS["tracked_tree"] and tracked_start["files"] == 10
            and tracked_start["bytes"] == 6_892_308, "immutable tracked release snapshot mismatch")

    tracked_inventory = parse_inventory(tracked_root, complete=False)
    external_inventory = parse_inventory(external)
    require(sha_file(tracked_root / "release.json") == PINS["release_json"], "tracked N=250 release digest mismatch")
    require(sha_file(external / "release.json") == PINS["release_json"], "external N=250 release digest mismatch")
    require(sha_file(external / "SHA256SUMS") == PINS["external_sha256sums"], "external inventory digest mismatch")
    require(sha_file(external / "COMPLETE") == PINS["complete"], "external COMPLETE digest mismatch")
    require(sha_file(tracked_root / "SHA256SUMS") == PINS["tracked_sha256sums"], "tracked inventory digest mismatch")
    require((tracked_root / "external_SHA256SUMS").read_bytes() == (external / "SHA256SUMS").read_bytes(),
            "tracked/external SHA inventory bytes differ")
    mapping = {
        "cohort-0250.tsv": "manifests/cohort-0250.tsv", "assemblies.tsv": "manifests/assemblies.tsv",
        "checksums.tsv": "manifests/checksums.tsv", "state.tsv": "manifests/state.tsv",
        "object_refs.tsv": "manifests/object_refs.tsv", "batch_metrics.tsv": "manifests/batch_metrics.tsv",
        "release.json": "release.json",
    }
    for tracked_name, external_name in mapping.items():
        require((tracked_root / tracked_name).read_bytes() == (external / external_name).read_bytes(),
                f"tracked/external manifest bytes differ: {tracked_name}")
    require(gzip.open(tracked_root / "contigs.tsv.gz", "rb").read() == (external / "manifests/contigs.tsv").read_bytes(),
            "tracked deterministic contig view differs from external manifest")

    release = json.loads((external / "release.json").read_text())
    require(release.get("release_id") == RELEASE_ID and release.get("verdict") == "PASS" and release.get("immutable") is True,
            "N=250 release identity/verdict/immutability mismatch")
    require(release.get("external_release_path") == str(external), "external release path contract mismatch")
    require(release.get("input_cohort_0250_sha256") == PINS["cohort_250"], "release cohort pin mismatch")
    require(release.get("predecessor_release_id") == N100_ID and release.get("predecessor_release_json_sha256") == PINS["n100_release"],
            "release N=100 predecessor pin mismatch")
    require(release.get("selection_release_id") == SELECTION_ID and release.get("selection_release_json_sha256") == PINS["selection_release"],
            "release selection pin mismatch")
    require(release.get("compatibility_release_id") == COMPATIBILITY_ID
            and release.get("compatibility_release_json_sha256") == PINS["compatibility_release"],
            "release compatibility pin mismatch")
    release_gates = release.get("applicable_gates", {})
    for name, verdict in release_gates.items():
        if name == "integrated_n100_go_250":
            require(verdict == "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY", "integrated gate was fabricated/bypassed")
        elif name == "scale_trend":
            require(verdict == "NOT_APPLICABLE_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS", "legacy scale applicability mismatch")
        else:
            require(verdict == "PASS", f"applicable release gate is not PASS: {name}={verdict}")

    selection_release, selection_external = require_release(
        repo, repo / "manifests/pilot-cohorts-v1", SELECTION_ID, PINS["selection_release"]
    )
    require(selection_release["manifests"]["cohort-0250.tsv"]["sha256"] == PINS["cohort_250"], "selection N=250 pin mismatch")
    n100_release, n100_external = require_release(
        repo, repo / "manifests/canonical-cohort-100-v1", N100_ID, PINS["n100_release"]
    )
    n10_release, n10_external = require_release(
        repo, repo / "manifests/canonical-cohort-010-v1", N10_ID, PINS["n10_release"]
    )
    compatibility, compatibility_external = require_release(
        repo, repo / "manifests/consumer-compatibility-v1", COMPATIBILITY_ID, PINS["compatibility_release"]
    )
    require(all(value == "PASS" or str(value).startswith("NOT_APPLICABLE")
                for value in compatibility.get("applicable_gates", {}).values()),
            "pinned compatibility release gate failure")
    _, consumers = read_plain_tsv(compatibility_external / "consumers.tsv")
    require(len(consumers) == 19 and all(row.get("verdict") == "PASS" for row in consumers),
            "pinned consumer compatibility row failure")

    _, cohort250 = read_tsv(repo / "manifests/pilot-cohorts-v1/cohort-0250.tsv")
    _, cohort100 = read_tsv(repo / "manifests/pilot-cohorts-v1/cohort-0100.tsv")
    _, cohort10 = read_tsv(repo / "manifests/pilot-cohorts-v1/cohort-0010.tsv")
    require(len(cohort250) == 250 and len(cohort100) == 100 and len(cohort10) == 10, "cohort cardinality mismatch")
    require(nested_rows(cohort250, cohort100) and nested_rows(cohort100, cohort10), "recursive N=250/N=100/N=10 nesting mismatch")
    accessions = [row["exact_assembly_accession_version"] for row in cohort250]
    require([int(row["cohort_order"]) for row in cohort250] == list(range(1, 251)), "N=250 cohort order mismatch")
    require(len(set(accessions)) == 250 and all(ACCESSION_RE.fullmatch(value) for value in accessions),
            "N=250 exact-version identity/uniqueness mismatch")
    require(release.get("sequence_bearing_assembly_accessions") == accessions, "release accession order mismatch")
    require((external / "manifests/cohort-0250.tsv").read_bytes() == (repo / "manifests/pilot-cohorts-v1/cohort-0250.tsv").read_bytes(),
            "external exact cohort bytes mismatch")

    _, assemblies = read_tsv(external / "manifests/assemblies.tsv")
    _, states = read_tsv(external / "manifests/state.tsv")
    _, refs = read_tsv(external / "manifests/object_refs.tsv")
    _, checksums = read_tsv(external / "manifests/checksums.tsv")
    _, batch_rows = read_tsv(external / "manifests/batch_metrics.tsv")
    contig_fields, cohort_contigs = read_tsv(external / "manifests/contigs.tsv")
    require((len(assemblies), len(states), len(refs), len(checksums), len(batch_rows), len(cohort_contigs))
            == (250, 250, 250, 2000, 25, 41050), "release manifest cardinality mismatch")
    for manifest_name, rows in {
        "cohort-0250.tsv": cohort250, "assemblies.tsv": assemblies, "state.tsv": states,
        "object_refs.tsv": refs, "checksums.tsv": checksums, "batch_metrics.tsv": batch_rows,
        "contigs.tsv": cohort_contigs,
    }.items():
        path = external / "manifests" / manifest_name
        record = release["manifests"][manifest_name]
        require(record["rows"] == len(rows) and record["bytes"] == path.stat().st_size and record["sha256"] == sha_file(path),
                f"release manifest accounting mismatch: {manifest_name}")
    require([row["accession"] for row in assemblies] == accessions
            and [row["accession"] for row in states] == accessions
            and [row["accession"] for row in refs] == accessions, "ordered row accounting mismatch")
    require(all(row["source_state"] == row["canonical_state"] == "COMPLETE" and row["terminal_state"] == "VALIDATED" for row in states),
            "non-terminal release state")
    require({(row["accession"], row["artifact_role"]) for row in checksums}
            == {(accession, role) for accession in accessions for role in ROLE_NAMES}, "artifact role inventory mismatch")

    _, refs100 = read_tsv(n100_external / "manifests/object_refs.tsv")
    require(len(refs100) == 100, "N=100 object reference cardinality mismatch")
    for index in range(10):
        require(refs[index] == refs100[index], f"N=10 recursive reference drift: {accessions[index]}")
    for index in range(10, 100):
        current, predecessor = refs[index], refs100[index]
        require(current["storage_release_id"] == N100_ID and Path(current["storage_root"]) == n100_external,
                f"N=100 recursive storage identity mismatch: {current['accession']}")
        for field in ("accession", "source_object_relpath", "source_inventory_sha256", "canonical_object_relpath", "canonical_inventory_sha256"):
            require(current[field] == predecessor[field], f"N=100 recursive object digest drift: {current['accession']}:{field}")
    for index in range(10):
        require(refs100[index]["storage_release_id"] == N10_ID and Path(refs100[index]["storage_root"]) == n10_external,
                f"N=100 did not recursively resolve N=10: {refs100[index]['accession']}")
    for index in range(100, 250):
        require(refs[index]["storage_release_id"] == "SELF" and refs[index]["storage_root"] == "."
                and refs[index]["reuse_status"] == "CREATED_OR_CHECKSUM_RESUMED",
                f"new N=250 object storage mismatch: {refs[index]['accession']}")

    checksum_map = {(row["accession"], row["artifact_role"]): row for row in checksums}
    contigs_by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cohort_contigs:
        contigs_by_accession[row["accession"]].append(row)
    old_assembly_map = {row["accession"]: row for row in read_tsv(n100_external / "manifests/assemblies.tsv")[1]}
    digest_fields = (
        "source_package_sha256", "source_decompressed_sha256", "source_gff_sha256", "canonical_bgzf_sha256",
        "fai_sha256", "gzi_sha256", "crosswalk_sha256", "annotation_aliases_sha256",
    )
    object_inventory_count = archive_count = bgzf_count = faidx_regions = upstream_md5_rows = 0
    annotation_valid = annotation_quarantined = annotation_features = total_bases = 0
    all_names: set[str] = set()
    recursive_roots: set[str] = set()
    for order, (assembly, ref) in enumerate(zip(assemblies, refs), 1):
        accession = assembly["accession"]
        root = external if ref["storage_release_id"] == "SELF" else Path(ref["storage_root"])
        recursive_roots.add(ref["storage_release_id"])
        source = root / ref["source_object_relpath"]
        canonical = root / ref["canonical_object_relpath"]
        source_inventory = parse_inventory(source)
        canonical_inventory = parse_inventory(canonical)
        object_inventory_count += 2
        require(sha_file(source / "SHA256SUMS") == ref["source_inventory_sha256"]
                and sha_file(canonical / "SHA256SUMS") == ref["canonical_inventory_sha256"],
                f"recursive object inventory pin mismatch: {accession}")
        source_manifest = json.loads((source / "manifest.json").read_text())
        canonical_manifest = json.loads((canonical / "manifest.json").read_text())
        require(source_manifest.get("accession") == canonical_manifest.get("accession") == accession
                and source_manifest.get("state") == canonical_manifest.get("state") == "COMPLETE",
                f"object manifest identity/state mismatch: {accession}")
        if order <= 100:
            old = old_assembly_map.get(accession)
            require(old is not None and all(assembly[field] == old[field] for field in digest_fields),
                    f"predecessor digest reuse mismatch: {accession}")
            require(ref["reuse_status"] == "REUSED_PREDECESSOR_BY_DIGEST" and ref["predecessor_digest_match"] == "PASS",
                    f"predecessor reuse status mismatch: {accession}")
        for role in ROLE_NAMES:
            row = checksum_map[(accession, role)]
            require(row["storage_release_id"] == ref["storage_release_id"] and row["storage_root"] == ref["storage_root"],
                    f"artifact storage reference mismatch: {accession}:{role}")
            path = root / row["relative_path"]
            require(path.is_file() and not path.is_symlink() and path.stat().st_size == int(row["bytes"])
                    and sha_file(path) == row["sha256"], f"artifact role checksum mismatch: {accession}:{role}")
        package = source / "package.zip"
        archive_stats, archive = validate_archive(package, accession)
        archive_count += 1
        upstream_md5_rows += archive_stats["md5_entries"]
        validation = source_manifest["validation"]
        for left, right in (
            (archive_stats["package_sha256"], assembly["source_package_sha256"]),
            (archive_stats["fasta_sha256"], assembly["source_decompressed_sha256"]),
            (archive_stats["gff_sha256"], assembly["source_gff_sha256"]),
            (archive_stats["fasta_member"], assembly["source_fasta_member"]),
            (archive_stats["gff_member"], assembly["source_gff_member"]),
        ):
            require(str(left) == str(right), f"archive/assembly manifest mismatch: {accession}")
        for field in ("package_sha256", "package_bytes", "fasta_member", "fasta_sha256", "fasta_bytes", "gff_member", "gff_sha256",
                      "gff_bytes", "sequence_report_member", "sequence_report_sha256", "dataset_catalog_member",
                      "dataset_catalog_sha256", "assembly_report_accessions", "zip_entries", "uncompressed_member_bytes", "md5_entries"):
            require(archive_stats[field] == validation[field], f"archive validation replay mismatch: {accession}:{field}")
        fasta_data = archive.read(archive_stats["fasta_member"])
        require(hashlib.sha256(fasta_data).hexdigest() == archive_stats["fasta_sha256"], f"FASTA member digest drift: {accession}")
        source_records = parse_source_fasta(fasta_data, accession)
        object_fields, object_contigs = read_tsv(canonical / "contigs.tsv")
        require(object_fields == contig_fields and object_contigs == contigs_by_accession[accession],
                f"object/cohort crosswalk mismatch: {accession}")
        require(len(object_contigs) == len(source_records) == int(assembly["contig_count"]),
                f"crosswalk contig count mismatch: {accession}")
        for index, (record, row) in enumerate(zip(source_records, object_contigs), 1):
            require(int(row["stage_b_order"]) == order and int(row["contig_order"]) == index
                    and row["accession"] == row["assembly_id"] == accession,
                    f"crosswalk identity/order mismatch: {accession}:{index}")
            require(row["source_fasta_member"] == archive_stats["fasta_member"]
                    and base64.b64decode(row["source_fasta_header_b64"]) == b">" + record["header"]
                    and base64.b64decode(row["source_fasta_id_token_b64"]) == record["token"],
                    f"crosswalk source identity mismatch: {accession}:{index}")
            require(row["contig_id_encoding"] == record["encoding"] and row["pansn_sample"] == accession
                    and row["pansn_haplotype"] == "1" and row["pansn_contig"] == record["pansn_contig"]
                    and row["pansn_sequence_name"] == row["fasta_seqid"] == record["name"],
                    f"PanSN reversibility mismatch: {accession}:{index}")
            require(int(row["contig_length"]) == record["length"] and row["contig_sequence_sha256"] == record["sha256"],
                    f"rename-only source crosswalk mismatch: {accession}:{index}")
            require(record["name"] not in all_names, f"global PanSN collision: {record['name']}")
            all_names.add(record["name"])
        bgzf = root / canonical_manifest["canonical_bgzf_relpath"]
        fai = root / canonical_manifest["fai_relpath"]
        gzi = root / canonical_manifest["gzi_relpath"]
        require(subprocess.run([bgzip, "-t", str(bgzf)], capture_output=True).returncode == 0,
                f"BGZF integrity failure: {accession}")
        canonical_records, content_sha = parse_canonical_bgzf(bgzf)
        require([(row["name"], row["length"], row["sha256"]) for row in canonical_records]
                == [(row["name"], row["length"], row["sha256"]) for row in source_records],
                f"rename-only canonical sequence mismatch: {accession}")
        require(content_sha == canonical_manifest["canonical_fasta_content_sha256"]
                and sha_file(bgzf) == canonical_manifest["canonical_bgzf_sha256"]
                and sha_file(fai) == canonical_manifest["fai_sha256"]
                and sha_file(gzi) == canonical_manifest["gzi_sha256"],
                f"canonical digest mismatch: {accession}")
        validate_fai(fai, canonical_records)
        validate_gzi(gzi, bgzf.stat().st_size)
        faidx_regions += validate_faidx_regions(samtools, bgzf, canonical_records)
        bgzf_count += 1
        reports = load_sequence_reports(archive, archive_stats["sequence_report_member"])
        expected_aliases, annotation_summary, annotation_error = recompute_annotation(
            archive, archive_stats["gff_member"], source_records, reports, accession
        )
        alias_fields, observed_aliases = read_tsv(canonical / "annotation_aliases.tsv")
        published = canonical_manifest["annotation"]
        if published["status"] == "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW":
            annotation_quarantined += 1
            require(not observed_aliases and published.get("failure_reason") == annotation_error,
                    f"annotation quarantine not reproducible: {accession}")
        else:
            annotation_valid += 1
            require(annotation_error is None and published == annotation_summary
                    and compare_alias_rows(observed_aliases, expected_aliases, alias_fields),
                    f"annotation alias/coordinate recomputation mismatch: {accession}")
            annotation_features += annotation_summary["feature_rows"]
        total_bases += sum(row["length"] for row in source_records)
        archive.close()
    require(object_inventory_count == 500 and archive_count == bgzf_count == 250 and faidx_regions == 41050,
            "deep object/archive/BGZF/index validation count mismatch")
    require(len(all_names) == len(cohort_contigs) == 41050 and total_bases == 1_276_442_466,
            "global PanSN/base accounting mismatch")
    require(annotation_valid == 237 and annotation_quarantined == 13,
            "annotation published/quarantine accounting mismatch")
    require(release["counts"]["annotations_alias_validated"] == annotation_valid
            and release["counts"]["annotation_views_quarantined"] == annotation_quarantined,
            "release annotation count mismatch")

    restart = json.loads((external / "restart_evidence.json").read_text())
    require(all(value is True for key, value in restart.items() if key != "schema"), "restart evidence contains failure")
    state_lines = (external / "state.jsonl").read_text().splitlines()
    state_events = [json.loads(line) for line in state_lines if line.strip()]
    names = [row.get("event") for row in state_events]
    for event in ("INJECTED_ACQUISITION_SIGKILL", "ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE",
                  "INJECTED_CONVERSION_SIGKILL", "INTERRUPTED_CONVERSION_STAGE_DISCARDED", "READY_TO_PROMOTE"):
        require(event in names, f"restart/atomicity event missing: {event}")
    require(names.index("INJECTED_ACQUISITION_SIGKILL") < names.index("ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE")
            and names.index("INJECTED_CONVERSION_SIGKILL") < names.index("INTERRUPTED_CONVERSION_STAGE_DISCARDED")
            and names[-1] == "READY_TO_PROMOTE", "restart event order/final promotion event mismatch")
    requested = {row.get("accession") for row in state_events if row.get("event") == "ACQUISITION_REQUEST_STARTED"}
    require(requested == set(accessions[100:]), "predecessor redownload or missing new acquisition evidence")
    require(not any(part.startswith(".stage.") or part.endswith(".partial") for path in external.rglob("*") for part in path.relative_to(external).parts),
            "partial/staging path survived atomic release")
    require(external.stat().st_dev == external.parent.stat().st_dev, "release is not on promotion parent's filesystem")

    resources = [json.loads(line) for line in (external / "resources.jsonl").read_text().splitlines() if line.strip()]
    require(len(resources) == 399, "resource preflight record count mismatch")
    for row in resources:
        allocations = row.get("allocations", {})
        require(row.get("verdict") == "PASS" and row.get("write_probes") == "PASS"
                and all(row.get("checks", {}).values()) and allocations and all(int(value) > 0 for value in allocations.values()),
                f"resource preflight failure/blank allocation: {row.get('stage')}")
        require(int(row["durable_free_bytes"]) >= 2_000_000_000_000 and int(row["durable_free_inodes"]) >= 1_000_000,
                f"durable floor failure: {row.get('stage')}")
        require(int(row["scratch_free_bytes"]) >= 4_000_000_000_000 and int(row["scratch_free_inodes"]) >= 5_000_000,
                f"scratch floor failure: {row.get('stage')}")
    resource = json.loads((external / "resource_summary.json").read_text())
    require(resource.get("verdict") == "PASS" and all(resource.get("checks", {}).values()), "resource summary gate failure")
    allocations = resource["allocations"]
    projection = resource["disk_projection"]
    require(resource["peak_rss_bytes"] * 100 <= allocations["assigned_ram_bytes"] * 70,
            "RAM 70% gate failure")
    require(resource["process_swap_events"] == resource["system_swap_growth_bytes"] == 0, "swap/OOM proxy gate failure")
    require(resource["measured_release_stage_peak_bytes"] * 100 <= allocations["durable_allocation_bytes"] * 70,
            "measured disk allocation gate failure")
    require(resource["measured_release_stage_peak_files"] * 2 <= allocations["inode_allocation"], "measured inode gate failure")
    require(projection["modeled_upper95_peak_bytes"] <= projection["configured_upper95_peak_bytes"]
            and resource["measured_release_stage_peak_bytes"] <= projection["modeled_upper95_peak_bytes"],
            "N=250 upper-95 projection gate failure")
    require(allocations["durable_allocation_bytes"] >= 2 * allocations["unfinished_write_bytes"]
            and allocations["scratch_allocation_bytes"] >= 2 * allocations["unfinished_write_bytes"],
            "unfinished-write allocation reservation failure")
    partial = [int(row["stage_partial_bytes_finish"]) for row in batch_rows]
    require(partial[-1] == 0 and all(value == 0 for value in partial[10:])
            and max(partial) <= allocations["unfinished_write_bytes"], "partial-write batch accounting failure")

    n100_assemblies = read_tsv(n100_external / "manifests/assemblies.tsv")[1]
    n100_resource = json.loads((n100_external / "resource_summary.json").read_text())
    wall100 = 4 * 60 + 26.78
    wall250 = 6 * 60 + 54.48
    new100 = n100_assemblies[10:]
    new250 = assemblies[100:]
    bases100 = sum(int(row["total_bases"]) for row in new100)
    bases250 = sum(int(row["total_bases"]) for row in new250)
    bytes100 = sum(int(row["source_package_bytes"]) for row in new100)
    bytes250 = sum(int(row["source_package_bytes"]) for row in new250)
    slopes100 = {
        "wall_seconds_per_new_object": wall100 / 90,
        "wall_seconds_per_new_base": wall100 / bases100,
        "wall_seconds_per_new_source_byte": wall100 / bytes100,
        "stage_bytes_per_new_object": n100_resource["measured_release_stage_peak_bytes"] / 90,
        "stage_files_per_new_object": n100_resource["measured_release_stage_peak_files"] / 90,
    }
    slopes250 = {
        "wall_seconds_per_new_object": wall250 / 150,
        "wall_seconds_per_new_base": wall250 / bases250,
        "wall_seconds_per_new_source_byte": wall250 / bytes250,
        "stage_bytes_per_new_object": resource["measured_release_stage_peak_bytes"] / 150,
        "stage_files_per_new_object": resource["measured_release_stage_peak_files"] / 150,
    }
    slope_changes = {key: slopes250[key] / slopes100[key] - 1 for key in slopes100}
    exponent = math.log(wall250 / wall100) / math.log(150 / 90)
    require(exponent <= 1.3 and all(abs(value) <= 0.25 for value in slope_changes.values()),
            "canonical preparation scale exponent/per-base slope gate failure")
    projected_new500 = 250
    project_ratio = projected_new500 / 150
    projected_500 = {
        "incremental_objects": projected_new500,
        "method": "N=250 measured incremental workload linearly scaled to 250 new objects plus 25% upper-95 allowance",
        "stage_upper95_bytes": math.ceil(resource["measured_release_stage_peak_bytes"] * project_ratio * 1.25),
        "files_upper95": math.ceil(resource["measured_release_stage_peak_files"] * project_ratio * 1.25),
        "rss_upper95_bytes": math.ceil(resource["peak_rss_bytes"] * project_ratio * 1.25),
        "wall_upper95_seconds": wall250 * project_ratio * 1.25,
        "comparison_durable_allocation_bytes": allocations["durable_allocation_bytes"],
        "comparison_inode_allocation": allocations["inode_allocation"],
        "comparison_assigned_ram_bytes": allocations["assigned_ram_bytes"],
    }
    require(projected_500["stage_upper95_bytes"] * 100 <= allocations["durable_allocation_bytes"] * 70
            and projected_500["files_upper95"] * 2 <= allocations["inode_allocation"]
            and projected_500["rss_upper95_bytes"] * 100 <= allocations["assigned_ram_bytes"] * 70,
            "conservative N=500 canonical resource projection exceeds comparison allocations")

    rerun = json.loads((repo / "artifacts/canonical_cohort_250/deterministic_rerun.json").read_text())
    require(rerun.get("verdict") == "PASS" and rerun.get("network_requests") == 0
            and rerun.get("objects_downloaded") == rerun.get("objects_recompressed") == 0,
            "published deterministic rerun side-effect failure")
    require(rerun["external_tree_sha256_before"] == rerun["external_tree_sha256_after"] == external_start["tree_sha256"]
            and rerun["tracked_tree_sha256_before"] == rerun["tracked_tree_sha256_after"] == tracked_start["tree_sha256"],
            "published deterministic rerun tree mismatch")
    require(rerun["state_sha256_before"] == rerun["state_sha256_after"] == sha_file(external / "state.jsonl"),
            "published deterministic rerun state mismatch")

    frozen1000 = {row["exact_assembly_accession_version"] for row in read_tsv(repo / "manifests/pilot-cohorts-v1/cohort-1000.tsv")[1]}
    union: set[str] = set()
    scanned = 0
    releases_root = external.parents[2]
    for release_path in sorted(releases_root.rglob("release.json")):
        candidate = json.loads(release_path.read_text())
        values = candidate.get("sequence_bearing_assembly_accessions", [])
        require(isinstance(values, list) and all(isinstance(value, str) and ACCESSION_RE.fullmatch(value) for value in values),
                f"invalid sequence-bearing release inventory: {release_path}")
        union.update(values)
        scanned += 1
    require(union == set(accessions) and union.issubset(frozen1000) and len(union) == 250 <= 1000,
            "global exact union/frozen cap mismatch")

    graph = graph_contract(repo)
    artifact_inventory = parse_inventory(repo / "artifacts/canonical_cohort_250", exact=True, complete=False)
    owned_roots = [tracked_root, repo / "artifacts/canonical_cohort_250", repo / "reports/canonical_cohort_250.md"]
    owned_files = [
        path for root in owned_roots
        for path in ([root] if root.is_file() else list(root.rglob("*")))
        if path.is_file()
    ]
    forbidden_suffixes = (".fa.gz", ".fna.gz", ".fasta.gz", ".fai", ".gzi", ".zip", ".gff")
    require(not any(path.stat().st_size > 10 * 1024 * 1024 for path in owned_files), "task-owned Git evidence exceeds 10 MiB")
    require(not any(path.name.endswith(forbidden_suffixes) for path in owned_files), "sequence/index/archive payload in task-owned Git paths")
    require(not any(path.suffix.lower() in {".fa", ".fna", ".fasta"} for path in external.rglob("*") if path.is_file()),
            "routine plain FASTA survived external publication")

    roots_finish = {name: sha_file(repo / name) for name in ROOT_HASHES}
    external_finish = tree_snapshot(external)
    tracked_finish = tree_snapshot(tracked_root)
    require(roots_finish == roots_start == ROOT_HASHES, "root inputs changed during audit")
    require(external_finish == external_start and tracked_finish == tracked_start, "immutable release bytes changed during audit")

    applicable_gates = {
        "exact_n250_identity_order_and_recursive_nesting": "PASS",
        "recursive_n100_n10_release_and_object_references": "PASS",
        "root_input_immutability": "PASS",
        "external_tracked_artifact_and_object_sha_inventories": "PASS",
        "archive_crc_upstream_md5_local_sha_and_exact_accession": "PASS",
        "bgzf_pansn_rename_only_fai_gzi_and_all_contig_faidx": "PASS",
        "crosswalk_and_source_gff_annotation_policy": "PASS",
        "deterministic_read_only_rerun": "PASS",
        "kill_resume_atomicity_and_no_predecessor_redownload": "PASS",
        "resource_floors_canonical_scale_and_n500_projection": "PASS",
        "global_exact_union_and_cap": "PASS",
        "compact_git_and_external_payload_policy": "PASS",
        "canonical_extraction_graph_separation": "PASS",
    }
    not_applicable_gates = {
        "prophage_extraction_source_coordinate_policy": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED",
        "integrated_extraction_go": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED",
        "integrated_biological_scale_trend": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS",
    }
    return {
        "schema": "canonical-cohort-250-contract-audit-v1",
        "audit_verdict": "PASS",
        "scope": "CANONICAL_ACQUISITION_AND_CANONICALIZATION_ONLY",
        "release_id": RELEASE_ID,
        "release_immutable": True,
        "applicable_gates": applicable_gates,
        "not_applicable_gates": not_applicable_gates,
        "counts": {
            "assemblies": len(assemblies), "contigs": len(cohort_contigs), "total_bases": total_bases,
            "artifact_role_checksums": len(checksums), "recursive_object_references": len(refs),
            "object_sha_inventories": object_inventory_count, "external_inventory_rows": len(external_inventory),
            "tracked_inventory_rows": len(tracked_inventory), "artifact_inventory_rows": len(artifact_inventory),
            "archives_revalidated": archive_count, "upstream_md5_rows_revalidated": upstream_md5_rows,
            "bgzf_revalidated": bgzf_count, "faidx_regions_revalidated": faidx_regions,
            "annotation_views_validated": annotation_valid, "annotation_views_quarantined": annotation_quarantined,
            "annotation_feature_rows_revalidated_in_published_views": annotation_features,
            "resource_preflight_rows": len(resources), "batch_metric_rows": len(batch_rows),
            "global_release_json_scanned": scanned, "global_distinct_exact_assembly_revisions": len(union),
        },
        "recursive_storage_release_ids": sorted(recursive_roots),
        "release_pins": {
            "tracked_release_json_sha256": PINS["release_json"],
            "external_release_json_sha256": PINS["release_json"],
            "external_sha256sums_sha256": PINS["external_sha256sums"],
            "complete_sha256": PINS["complete"],
            "external_tree_sha256": external_start["tree_sha256"],
            "tracked_tree_sha256": tracked_start["tree_sha256"],
            "cohort_0250_sha256": PINS["cohort_250"],
            "predecessor_n100_release_json_sha256": PINS["n100_release"],
            "recursive_n10_release_json_sha256": PINS["n10_release"],
        },
        "resource": {
            "canonical_scale_exponent": exponent, "canonical_slope_changes": slope_changes,
            "n100_slopes": slopes100, "n250_slopes": slopes250, "n500_projection": projected_500,
            "durable_floor_bytes": 2_000_000_000_000, "durable_inode_floor": 1_000_000,
            "scratch_preflight_bytes": 4_000_000_000_000, "scratch_inode_floor": 5_000_000,
        },
        "graph_contract": graph,
        "side_effects": {
            "network_requests": 0, "sequence_downloads": 0, "objects_recomputed": 0,
            "release_writes": 0, "root_inputs_unchanged": True,
            "external_snapshot_start_equals_finish": True, "tracked_snapshot_start_equals_finish": True,
        },
        "snapshots": {"external": external_finish, "tracked": tracked_finish, "root_inputs": roots_finish},
    }


def build_verdict(result: dict[str, Any]) -> dict[str, Any]:
    audit_sha = hashlib.sha256(canonical_json(result)).hexdigest()
    return {
        "schema": "canonical-scale-verdict-v1",
        "verdict": "CANONICAL_GO_500",
        "consumer_task_id": "prepare-canonical-cohort-500",
        "scope": "CANONICAL_ACQUISITION_AND_CANONICALIZATION_ONLY",
        "source_release_id": RELEASE_ID,
        "source_release_json_sha256": PINS["release_json"],
        "source_external_sha256sums_sha256": PINS["external_sha256sums"],
        "source_complete_sha256": PINS["complete"],
        "source_external_tree_sha256": PINS["external_tree"],
        "audit_result_sha256": audit_sha,
        "applicable_gate_rule": "every applicable_gates value must equal PASS",
        "applicable_gates": result["applicable_gates"],
        "not_applicable_gates": result["not_applicable_gates"],
        "extraction_branch": {
            "task_id": "resolve-prophage-source", "task_status": "failed",
            "verdict": "EXTRACTION_BLOCKED", "consumer_action": "REJECT",
            "integrated_go_claimed": False,
        },
        "conditions_for_consumer": [
            "verify this file against artifacts/canonical_cohort_250_contract_audit/SHA256SUMS",
            "verify audit_result_sha256 against audit.json",
            "verify all source release pins before any N=500 work",
            "repeat live resource/root/global-cap preflight; this verdict does not waive N=500 gates",
            "do not use this verdict for extraction, source-coordinate interpretation, or integrated analysis",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--external-release", default=str(EXTERNAL_RELEASE))
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bgzip", default="bgzip")
    parser.add_argument("--audit-output")
    parser.add_argument("--verdict-output")
    args = parser.parse_args()
    try:
        result = audit(Path(args.repo_root), Path(args.external_release), args.samtools, args.bgzip)
        verdict = build_verdict(result)
        if args.audit_output:
            Path(args.audit_output).write_bytes(canonical_json(result))
        if args.verdict_output:
            Path(args.verdict_output).write_bytes(canonical_json(verdict))
        print(canonical_json(result).decode(), end="")
        return 0
    except (AuditError, OSError, ValueError, KeyError, json.JSONDecodeError, zipfile.BadZipFile) as exc:
        no_go = {
            "schema": "canonical-scale-verdict-v1", "verdict": "CANONICAL_NO_GO_500",
            "scope": "CANONICAL_ACQUISITION_AND_CANONICALIZATION_ONLY", "reason": str(exc),
            "integrated_go_claimed": False,
        }
        if args.verdict_output:
            Path(args.verdict_output).write_bytes(canonical_json(no_go))
        print(f"CANONICAL_NO_GO_500: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
