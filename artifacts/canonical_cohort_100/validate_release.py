#!/usr/bin/env python3
"""Independent checksum and semantic validator for canonical-cohort-100-v1."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from workflow.acquisition_canonicalization import pilot as p
from artifacts.canonical_cohort_100 import runner as r


class ValidationError(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def _root_for_ref(external: Path, ref: dict[str, str], predecessor: Path) -> Path:
    if ref["storage_release_id"] == r.PREDECESSOR_RELEASE_ID:
        check(ref["storage_root"] == str(predecessor), "predecessor storage root mismatch")
        return predecessor
    check(ref["storage_release_id"] == "SELF" and ref["storage_root"] == ".", "self storage contract mismatch")
    return external


def _validate_region_roundtrip(bgzf: Path, parsed: list[dict[str, Any]], samtools: str) -> None:
    for row in parsed:
        end = min(60, int(row["length"]))
        result = subprocess.run([samtools, "faidx", str(bgzf), f"{row['name']}:1-{end}"], capture_output=True)
        check(result.returncode == 0, f"samtools region lookup failed: {bgzf}:{row['name']}")
        bases = b"".join(result.stdout.splitlines()[1:])
        check(bases == row["prefix"][:end], f"samtools region bases mismatch: {bgzf}:{row['name']}")


def validate(repo: Path, external: Path, samtools: str = "samtools", bgzip: str = "bgzip") -> dict[str, Any]:
    groups: list[str] = []
    inputs = r.verify_inputs(repo)
    r.configure_pinned_primitives(inputs["accessions"])
    groups.append("immutable_roots_and_pinned_predecessors")
    p.verify_sha_inventory(external)
    release = json.loads((external / "release.json").read_text())
    check(release.get("release_id") == r.release_id(inputs["accessions"]), "release ID mismatch")
    check(release.get("verdict") == "PASS" and release.get("immutable") is True, "release verdict/immutability mismatch")
    check(release.get("selection_release_id") == r.SELECTION_RELEASE_ID, "selection release ID mismatch")
    check(release.get("selection_release_json_sha256") == r.SELECTION_RELEASE_JSON_SHA256, "selection release SHA mismatch")
    check(release.get("input_cohort_0100_sha256") == r.COHORT_SHA256, "cohort SHA mismatch")
    check(release.get("predecessor_release_id") == r.PREDECESSOR_RELEASE_ID, "predecessor release mismatch")
    check(release.get("predecessor_release_json_sha256") == r.PREDECESSOR_RELEASE_JSON_SHA256, "predecessor release SHA mismatch")
    check(release.get("compatibility_release_id") == r.COMPATIBILITY_RELEASE_ID, "compatibility release mismatch")
    check(release.get("compatibility_release_json_sha256") == r.COMPATIBILITY_RELEASE_JSON_SHA256, "compatibility release SHA mismatch")
    check(release.get("sequence_bearing_assembly_accessions") == inputs["accessions"], "cohort order mismatch")
    for name, verdict in release.get("applicable_gates", {}).items():
        if name == "scale_trend":
            check(verdict == "NOT_APPLICABLE_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS", "scale gate applicability mismatch")
        else:
            check(verdict == "PASS", f"release gate is not PASS: {name}={verdict}")
    groups.append("release_identity_and_all_applicable_gates")

    cohort_external = external / "manifests/cohort-0100.tsv"
    cohort_tracked = repo / "manifests/pilot-cohorts-v1/cohort-0100.tsv"
    check(cohort_external.read_bytes() == cohort_tracked.read_bytes(), "external exact cohort bytes changed")
    check(r.sha_file(cohort_external) == r.COHORT_SHA256, "external cohort SHA mismatch")
    groups.append("exact_input_manifest_bytes_order_rows")

    assemblies = p.read_tsv(external / "manifests/assemblies.tsv", p.ASSEMBLY_FIELDS, verify_hashes=True)
    states = p.read_tsv(external / "manifests/state.tsv", p.STATE_FIELDS, verify_hashes=True)
    checksums = p.read_tsv(external / "manifests/checksums.tsv", r.CHECKSUM_FIELDS, verify_hashes=True)
    contigs = p.read_tsv(external / "manifests/contigs.tsv", p.CONTIG_FIELDS, verify_hashes=True)
    refs = p.read_tsv(external / "manifests/object_refs.tsv", r.OBJECT_REF_FIELDS, verify_hashes=True)
    batches = p.read_tsv(external / "manifests/batch_metrics.tsv", r.BATCH_FIELDS, verify_hashes=True)
    accessions = inputs["accessions"]
    check(len(assemblies) == len(states) == len(refs) == 100, "assembly/state/ref cardinality mismatch")
    check(len(checksums) == 800, "checksum row cardinality mismatch")
    check(len(batches) == 10, "bounded batch row cardinality mismatch")
    check([row["accession"] for row in assemblies] == accessions, "assembly order mismatch")
    check([row["accession"] for row in states] == accessions, "state order mismatch")
    check([row["accession"] for row in refs] == accessions, "reference order mismatch")
    check(all(row["terminal_state"] in ("VALIDATED", "QUARANTINED") for row in states), "non-terminal assembly state")
    check(all(row["terminal_state"] == "VALIDATED" for row in states), "genome object quarantine present")
    check(sum(int(row["contig_count"]) for row in assemblies) == len(contigs), "contig accounting mismatch")
    check(len({row["pansn_sequence_name"] for row in contigs}) == len(contigs), "PanSN collision")
    groups.append("complete_terminal_row_accounting")

    checksum_by_key = {(row["accession"], row["artifact_role"]): row for row in checksums}
    check(len(checksum_by_key) == 800, "duplicate checksum role rows")
    contigs_by_accession: dict[str, list[dict[str, str]]] = {}
    for row in contigs:
        contigs_by_accession.setdefault(row["accession"], []).append(row)
    predecessor = inputs["predecessor_external"]
    old_assemblies = {row["accession"]: row for row in p.read_tsv(predecessor / "manifests/assemblies.tsv", p.ASSEMBLY_FIELDS, verify_hashes=True)}
    reuse_digest_fields = ["source_package_sha256", "source_decompressed_sha256", "source_gff_sha256", "canonical_bgzf_sha256", "fai_sha256", "gzi_sha256", "crosswalk_sha256", "annotation_aliases_sha256"]
    total_bases = 0
    annotation_quarantines = 0
    for order, (assembly, ref) in enumerate(zip(assemblies, refs), 1):
        accession = assembly["accession"]
        check(int(assembly["stage_b_order"]) == order and int(ref["cohort_order"]) == order, f"cohort order mismatch: {accession}")
        root = _root_for_ref(external, ref, predecessor)
        source = root / ref["source_object_relpath"]
        canonical = root / ref["canonical_object_relpath"]
        p.verify_sha_inventory(source)
        p.verify_sha_inventory(canonical)
        check(r.sha_file(source / "SHA256SUMS") == ref["source_inventory_sha256"], f"source object inventory reference mismatch: {accession}")
        check(r.sha_file(canonical / "SHA256SUMS") == ref["canonical_inventory_sha256"], f"canonical object inventory reference mismatch: {accession}")
        if order <= 10:
            check(ref["reuse_status"] == "REUSED_PREDECESSOR_BY_DIGEST" and ref["predecessor_digest_match"] == "PASS", f"predecessor reuse status mismatch: {accession}")
            old = old_assemblies.get(accession)
            check(old is not None and all(assembly[field] == old[field] for field in reuse_digest_fields), f"N=10 checksum reuse mismatch: {accession}")
        else:
            check(root == external and ref["reuse_status"] == "CREATED_OR_CHECKSUM_RESUMED", f"new object storage mismatch: {accession}")
        source_manifest = json.loads((source / "manifest.json").read_text())
        canonical_manifest = json.loads((canonical / "manifest.json").read_text())
        check(source_manifest.get("accession") == accession and source_manifest.get("state") == "COMPLETE", f"source manifest mismatch: {accession}")
        check(canonical_manifest.get("accession") == accession and canonical_manifest.get("state") == "COMPLETE", f"canonical manifest mismatch: {accession}")
        package_validation = p.validate_package(source / "package.zip", accession)
        check(package_validation["package_sha256"] == assembly["source_package_sha256"], f"source package SHA mismatch: {accession}")
        check(package_validation["fasta_sha256"] == assembly["source_decompressed_sha256"], f"source FASTA SHA mismatch: {accession}")
        bgzf_path = root / canonical_manifest["canonical_bgzf_relpath"]
        fai = root / canonical_manifest["fai_relpath"]
        gzi = root / canonical_manifest["gzi_relpath"]
        check(subprocess.run([bgzip, "-t", str(bgzf_path)], capture_output=True).returncode == 0, f"BGZF integrity failure: {accession}")
        p._validate_gzi(gzi)
        parsed, content_sha = p._parse_canonical_bgzf(bgzf_path)
        object_contigs = p.read_tsv(canonical / "contigs.tsv", p.CONTIG_FIELDS, verify_hashes=True)
        expected = [(row["pansn_sequence_name"], int(row["contig_length"]), row["contig_sequence_sha256"]) for row in object_contigs]
        observed = [(row["name"], row["length"], row["sha256"]) for row in parsed]
        check(observed == expected, f"rename-only sequence/order mismatch: {accession}")
        check(p._parse_fai(fai) == [(name, length) for name, length, _ in expected], f"FAI names/length mismatch: {accession}")
        check(content_sha == canonical_manifest["canonical_fasta_content_sha256"], f"canonical content SHA mismatch: {accession}")
        check(r.sha_file(bgzf_path) == canonical_manifest["canonical_bgzf_sha256"], f"BGZF SHA mismatch: {accession}")
        check(r.sha_file(fai) == canonical_manifest["fai_sha256"], f"FAI SHA mismatch: {accession}")
        check(r.sha_file(gzi) == canonical_manifest["gzi_sha256"], f"GZI SHA mismatch: {accession}")
        _validate_region_roundtrip(bgzf_path, parsed, samtools)
        expected_cohort_contigs = contigs_by_accession.get(accession, [])
        check([(row["pansn_sequence_name"], row["contig_sequence_sha256"]) for row in object_contigs] == [(row["pansn_sequence_name"], row["contig_sequence_sha256"]) for row in expected_cohort_contigs], f"cohort/object crosswalk mismatch: {accession}")
        annotation = canonical_manifest["annotation"]
        check(annotation["status"] in ("ALIASES_VALIDATED_NO_TRANSFORMED_GFF", "NOT_AVAILABLE", "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW"), f"annotation state invalid: {accession}")
        aliases = p.read_tsv(canonical / "annotation_aliases.tsv", verify_hashes=True)
        if annotation["status"] == "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW":
            annotation_quarantines += 1
            check(not aliases and annotation.get("failure_reason"), f"quarantined annotation exposed aliases: {accession}")
        for role in ("source_package", "source_manifest", "canonical_bgzf", "fai", "gzi", "contig_crosswalk", "annotation_aliases", "canonical_manifest"):
            checksum = checksum_by_key.get((accession, role))
            check(checksum is not None, f"missing checksum role: {accession}:{role}")
            check(checksum["storage_release_id"] == ref["storage_release_id"] and checksum["storage_root"] == ref["storage_root"], f"checksum storage mismatch: {accession}:{role}")
            path = root / checksum["relative_path"]
            check(path.is_file() and not path.is_symlink(), f"checksum artifact missing/symlink: {path}")
            check(path.stat().st_size == int(checksum["bytes"]) and r.sha_file(path) == checksum["sha256"], f"checksum artifact mismatch: {path}")
        total_bases += sum(length for _, length, _ in expected)
    groups.append("all_100_source_bgzf_index_crosswalk_annotation_semantics")
    check(total_bases == release["counts"]["total_bases"], "total base accounting mismatch")
    check(annotation_quarantines == release["counts"]["annotation_views_quarantined"], "annotation quarantine accounting mismatch")

    state_text = (external / "state.jsonl").read_text()
    restart = json.loads((external / "restart_evidence.json").read_text())
    check(all(value is True for key, value in restart.items() if key != "schema"), "restart evidence contains failure")
    check("INJECTED_ACQUISITION_SIGKILL" in state_text, "acquisition SIGKILL event absent")
    check("INJECTED_CONVERSION_SIGKILL" in state_text, "conversion SIGKILL event absent")
    state_events = [json.loads(line) for line in state_text.splitlines()]
    check(state_events[-1].get("event") == "READY_TO_PROMOTE", "READY_TO_PROMOTE not final state event")
    request_accessions = {event.get("accession") for event in state_events if event.get("event") == "ACQUISITION_REQUEST_STARTED"}
    check(not request_accessions.intersection(accessions[:10]), "validated predecessor object was re-downloaded")
    check(request_accessions.issubset(set(accessions[10:])), "download outside frozen new 90 rows")
    groups.append("bounded_retry_kill_restart_no_predecessor_redownload")

    resources = [json.loads(line) for line in (external / "resources.jsonl").read_text().splitlines()]
    check(resources and all(row.get("verdict") == "PASS" and all(row.get("checks", {}).values()) for row in resources), "resource preflight failure")
    summary = json.loads((external / "resource_summary.json").read_text())
    check(summary.get("verdict") == "PASS" and all(summary.get("checks", {}).values()), "resource summary failure")
    check(summary["peak_rss_fraction"] <= 0.70, "RAM gate failure")
    check(summary["process_swap_events"] == 0 and summary["system_swap_growth_bytes"] == 0, "swap gate failure")
    check(summary["scale_trend"] == "NOT_APPLICABLE_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS", "scale applicability mismatch")
    groups.append("live_resources_metrics_and_allocation_thresholds")

    root_finish = r.verify_root_inputs(repo)
    check(root_finish == release["root_inputs_start"] == release["root_inputs_finish"], "root hashes changed")
    allowed = p._load_allowed_collection(repo)
    global_finish = r.audit_global_release_cap(external.parents[2], allowed)
    check(global_finish["distinct_exact_assembly_revisions"] == 100, "finish union cardinality is not 100")
    check(set(global_finish["accessions"]) == set(accessions), "finish union differs from exact frozen cohort")
    groups.append("root_hashes_and_global_exact_cap_finish")

    forbidden = [str(path.relative_to(external)) for path in external.rglob("*") if path.is_file() and path.suffix.lower() in (".fa", ".fna", ".fasta")]
    check(not forbidden, f"plain FASTA in release: {forbidden}")
    check(not any(path.name.startswith(".stage.") for path in external.rglob("*")), "nested staging path in final release")
    check(not any(path.is_symlink() for path in external.rglob("*")), "symlink in final release")
    groups.append("atomic_inventory_no_partial_plain_fasta_or_symlink")

    return {
        "schema": "canonical-cohort-100-validation-v1", "verdict": "PASS",
        "release_id": release["release_id"], "external_release_path": str(external),
        "check_groups": groups, "check_group_count": len(groups), "assembly_rows": len(assemblies),
        "contig_rows": len(contigs), "checksum_rows": len(checksums), "object_reference_rows": len(refs),
        "batch_rows": len(batches), "total_bases": total_bases, "annotation_view_quarantines": annotation_quarantines,
        "reused_predecessor_assemblies": 10, "new_assemblies": 90,
        "global_distinct_exact_assembly_revisions": global_finish["distinct_exact_assembly_revisions"],
        "root_inputs": root_finish, "large_payloads_in_git": 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--external-release", required=True)
    parser.add_argument("--output")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bgzip", default="bgzip")
    args = parser.parse_args()
    try:
        result = validate(Path(args.repo_root).resolve(), Path(args.external_release).resolve(), args.samtools, args.bgzip)
        data = r.canonical_json(result)
        if args.output:
            output = Path(args.output)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(data)
        print(data.decode(), end="")
        return 0
    except (p.GateError, r.GateError, ValidationError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
