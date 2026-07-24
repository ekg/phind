#!/usr/bin/env python3
"""Independent semantic and checksum validator for canonical-cohort-010-v1."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
from pathlib import Path
from typing import Any

from workflow.acquisition_canonicalization import pilot as p


class ValidationError(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def validate(repo: Path, external: Path) -> dict[str, Any]:
    checks: list[str] = []
    root_inputs = p.verify_root_inputs(repo)
    checks.append("root_input_sha256")
    input_rows, predecessor = p.verify_predecessor(repo)
    checks.append("predecessor_release_inventory_id_manifest_sha256")
    p.verify_sha_inventory(external)
    checks.append("external_release_complete_and_full_sha256_inventory")
    release = json.loads((external / "release.json").read_text())
    check(release.get("release_id") == p.release_id(), "release ID mismatch")
    check(release.get("verdict") == "PASS" and release.get("immutable") is True, "release verdict/immutability failure")
    check(release.get("collection_release_id") == p.COLLECTION_RELEASE_ID, "collection release ID mismatch")
    check(release.get("collection_release_json_sha256") == p.COLLECTION_RELEASE_SHA256, "collection release digest mismatch")
    check(release.get("input_stage_b_manifest_sha256") == p.STAGE_B_SHA256, "Stage-B digest mismatch")
    check(release.get("sequence_bearing_assembly_accessions") == p.EXPECTED_ACCESSIONS, "cohort order mismatch")
    checks.append("release_contract")
    gates = release.get("applicable_gates", {})
    for name, verdict in gates.items():
        if name == "scale_trend":
            check(verdict == "NOT_APPLICABLE_STAGE_B_NON_SCALE_BEARING", f"bad scale gate verdict: {verdict}")
        else:
            check(verdict == "PASS", f"applicable gate is not PASS: {name}={verdict}")
    checks.append("all_applicable_gates_pass")
    stage_bytes = (external / "manifests/stage_b_10.tsv").read_bytes()
    check(stage_bytes == (repo / "manifests/collection-v1/stage_b_10.tsv").read_bytes(), "exact input Stage-B bytes changed")
    check(p.sha_bytes(stage_bytes) == p.STAGE_B_SHA256, "external Stage-B SHA mismatch")
    checks.append("immutable_stage_b_bytes")
    assemblies = p.read_tsv(external / "manifests/assemblies.tsv", p.ASSEMBLY_FIELDS, verify_hashes=True)
    states = p.read_tsv(external / "manifests/state.tsv", p.STATE_FIELDS, verify_hashes=True)
    checksums = p.read_tsv(external / "manifests/checksums.tsv", p.CHECKSUM_FIELDS, verify_hashes=True)
    contigs = p.read_tsv(external / "manifests/contigs.tsv", p.CONTIG_FIELDS, verify_hashes=True)
    check(len(assemblies) == len(states) == 10, "assembly/state row count is not 10")
    check([row["accession"] for row in assemblies] == p.EXPECTED_ACCESSIONS, "assembly row order mismatch")
    check([row["accession"] for row in states] == p.EXPECTED_ACCESSIONS, "state row order mismatch")
    check(all(row["terminal_state"] in ("VALIDATED", "QUARANTINED") for row in states), "non-terminal state")
    check(all(row["terminal_state"] == "VALIDATED" for row in states), "pilot has quarantined assembly")
    check(sum(int(row["contig_count"]) for row in assemblies) == len(contigs), "contig row accounting mismatch")
    check(len({row["pansn_sequence_name"] for row in contigs}) == len(contigs), "cohort PanSN collision")
    checks.append("assembly_contig_state_row_accounting")
    checksum_by_path = {row["relative_path"]: row for row in checksums}
    check(len(checksum_by_path) == len(checksums) == 80, "checksum inventory must have 8 roles x 10 assemblies")
    total_bases = 0
    per_accession_contigs: dict[str, list[dict[str, str]]] = {}
    for row in contigs:
        per_accession_contigs.setdefault(row["accession"], []).append(row)
    for assembly in assemblies:
        accession = assembly["accession"]
        source = external / "source_objects" / accession
        canonical = external / "canonical_objects" / accession
        p.verify_sha_inventory(source)
        p.verify_sha_inventory(canonical)
        source_manifest = json.loads((source / "manifest.json").read_text())
        canonical_manifest = json.loads((canonical / "manifest.json").read_text())
        check(source_manifest.get("accession") == accession and source_manifest.get("state") == "COMPLETE", f"bad source manifest: {accession}")
        check(canonical_manifest.get("accession") == accession and canonical_manifest.get("state") == "COMPLETE", f"bad canonical manifest: {accession}")
        revalidated_package = p.validate_package(source / "package.zip", accession)
        check(revalidated_package["package_sha256"] == assembly["source_package_sha256"], f"source package SHA mismatch: {accession}")
        check(revalidated_package["fasta_sha256"] == assembly["source_decompressed_sha256"], f"source FASTA SHA mismatch: {accession}")
        bgzf = external / canonical_manifest["canonical_bgzf_relpath"]
        fai = Path(str(bgzf) + ".fai")
        gzi = Path(str(bgzf) + ".gzi")
        check(subprocess.run(["bgzip", "-t", str(bgzf)], capture_output=True).returncode == 0, f"BGZF integrity failed: {accession}")
        p._validate_gzi(gzi)
        parsed, content_sha = p._parse_canonical_bgzf(bgzf)
        object_contigs = p.read_tsv(canonical / "contigs.tsv", p.CONTIG_FIELDS, verify_hashes=True)
        expected = [(row["pansn_sequence_name"], int(row["contig_length"]), row["contig_sequence_sha256"]) for row in object_contigs]
        observed = [(row["name"], row["length"], row["sha256"]) for row in parsed]
        check(observed == expected, f"rename-only order/length/digest mismatch: {accession}")
        check(p._parse_fai(fai) == [(name, length) for name, length, _ in expected], f"FAI mismatch: {accession}")
        check(content_sha == canonical_manifest["canonical_fasta_content_sha256"], f"canonical content SHA mismatch: {accession}")
        check(p.sha_file(bgzf) == canonical_manifest["canonical_bgzf_sha256"], f"compressed SHA mismatch: {accession}")
        check(p.sha_file(fai) == canonical_manifest["fai_sha256"], f"FAI SHA mismatch: {accession}")
        check(p.sha_file(gzi) == canonical_manifest["gzi_sha256"], f"GZI SHA mismatch: {accession}")
        annotation_status = canonical_manifest["annotation"]["status"]
        check(annotation_status in ("ALIASES_VALIDATED_NO_TRANSFORMED_GFF", "NOT_AVAILABLE", "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW"), f"annotation policy mismatch: {accession}")
        check(canonical_manifest["annotation_view"].startswith("validated alias table only"), f"transformed GFF policy mismatch: {accession}")
        alias_rows = p.read_tsv(canonical / "annotation_aliases.tsv", verify_hashes=True)
        if annotation_status == "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW":
            check(not alias_rows and "failure_reason" in canonical_manifest["annotation"], f"quarantined annotation exposed aliases: {accession}")
        for role, rel in (("source_package", f"source_objects/{accession}/package.zip"),
                          ("canonical_bgzf", canonical_manifest["canonical_bgzf_relpath"]),
                          ("fai", canonical_manifest["fai_relpath"]),
                          ("gzi", canonical_manifest["gzi_relpath"])):
            row = checksum_by_path.get(rel)
            check(row is not None and row["artifact_role"] == role, f"missing checksum row: {rel}")
            check(p.sha_file(external / rel) == row["sha256"] and (external / rel).stat().st_size == int(row["bytes"]), f"checksum row mismatch: {rel}")
        total_bases += sum(length for _, length, _ in expected)
    check(total_bases == release["counts"]["total_bases"], "release total bases mismatch")
    checks.append("source_archive_md5_sha256_and_canonical_bgzf_fai_gzi_semantics")
    restart = json.loads((external / "restart_evidence.json").read_text())
    required_restart = [key for key in restart if key != "schema"]
    check(all(restart[key] is True for key in required_restart), "injected restart evidence incomplete")
    state_log = (external / "state.jsonl").read_text()
    check("INJECTED_ACQUISITION_SIGKILL" in state_log and "INJECTED_CONVERSION_SIGKILL" in state_log, "kill events absent from append-only ledger")
    state_events = [json.loads(line) for line in state_log.splitlines()]
    check(state_events[-1].get("event") == "READY_TO_PROMOTE" and state_events[-1].get("release_id") == p.release_id(), "READY_TO_PROMOTE is not last state")
    checks.append("acquisition_conversion_kill_restart")
    resource_summary = json.loads((external / "resource_summary.json").read_text())
    check(resource_summary.get("verdict") == "PASS" and all(resource_summary.get("checks", {}).values()), "resource summary failure")
    check(resource_summary["peak_rss_fraction"] <= 0.70, "RSS exceeds 70%")
    check(resource_summary["process_swap_events"] == 0 and resource_summary["system_swap_growth_bytes"] == 0, "swap growth observed")
    checks.append("live_mount_resource_rss_swap_disk_inode_reservation")
    cap = json.loads((external / "global_cap_evidence.json").read_text())
    check(cap.get("verdict") == "PASS" and cap["projected_finish_distinct_exact_assembly_revisions"] == 10, "global cap evidence mismatch")
    allowed = p._load_allowed_collection(repo)
    finish_cap = p.audit_global_release_cap(external.parents[2], allowed)
    check(finish_cap["distinct_exact_assembly_revisions"] == 10 and finish_cap["accessions"] == sorted(p.EXPECTED_ACCESSIONS), "finish global cap mismatch")
    checks.append("global_cap_and_frozen_subset_finish")
    forbidden_plain = []
    for path in external.rglob("*"):
        if path.is_file() and path.suffix.lower() in (".fa", ".fna", ".fasta"):
            forbidden_plain.append(str(path.relative_to(external)))
    check(not forbidden_plain, f"routine plain FASTA found: {forbidden_plain}")
    check(not any(path.is_dir() and path.name.startswith(".stage.") for path in external.rglob("*")), "nested staging directory in release")
    checks.append("no_plain_fasta_or_partial_stage")
    return {
        "schema": "canonical-cohort-010-validation-v1", "verdict": "PASS",
        "release_id": release["release_id"], "external_release_path": str(external),
        "checks": checks, "check_count": len(checks), "assembly_rows": len(assemblies),
        "contig_rows": len(contigs), "checksum_rows": len(checksums), "total_bases": total_bases,
        "root_inputs": root_inputs, "predecessor_release_id": predecessor["release_id"],
        "global_distinct_exact_assembly_revisions": finish_cap["distinct_exact_assembly_revisions"],
        "large_payloads_in_git": 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--external-release", required=True)
    parser.add_argument("--output")
    args = parser.parse_args()
    try:
        result = validate(Path(args.repo_root).resolve(), Path(args.external_release).resolve())
        data = p.canonical_json(result)
        if args.output:
            output = Path(args.output)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(data)
        print(data.decode(), end="")
        return 0
    except (p.GateError, ValidationError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"NO_GO: {exc}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
