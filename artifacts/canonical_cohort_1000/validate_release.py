#!/usr/bin/env python3
"""Independent checksum and semantic validator for canonical-cohort-1000-v1."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any

from workflow.acquisition_canonicalization import pilot as p
from artifacts.canonical_cohort_1000 import runner as r


class ValidationError(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def _root_for_ref(external: Path, ref: dict[str, str]) -> Path:
    if ref["storage_release_id"] == "SELF":
        check(ref["storage_root"] == ".", "self storage contract mismatch")
        return external
    root = Path(ref["storage_root"])
    check(root.is_absolute(), "predecessor storage root is not absolute")
    return root


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
    check(release.get("input_cohort_1000_sha256") == r.COHORT_SHA256, "cohort SHA mismatch")
    check(release.get("predecessor_release_id") == r.PREDECESSOR_RELEASE_ID, "predecessor release mismatch")
    check(release.get("predecessor_release_json_sha256") == r.PREDECESSOR_RELEASE_JSON_SHA256, "predecessor release SHA mismatch")
    check(release.get("predecessor_scale_trend_verdict") == "PASS", "canonical N=500 scale trend verdict missing")
    check(release.get("predecessor_scale_trend_sha256") == r.PREDECESSOR_SCALE_TREND_SHA256, "canonical N=500 scale trend checksum mismatch")
    check(release.get("predecessor_n1000_projection_verdict") == "PASS", "canonical N=1,000 projection verdict missing")
    check(release.get("compatibility_release_id") == r.COMPATIBILITY_RELEASE_ID, "compatibility release mismatch")
    check(release.get("compatibility_release_json_sha256") == r.COMPATIBILITY_RELEASE_JSON_SHA256, "compatibility release SHA mismatch")
    check(release.get("sequence_bearing_assembly_accessions") == inputs["accessions"], "cohort order mismatch")
    expected_na = {
        "prophage_source_coordinate_policy": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED",
        "integrated_extraction_verdict": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED",
        "integrated_biological_scale_trend": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS",
    }
    gates = release.get("applicable_gates", {})
    check(gates.get("canonical_n500_scale_trend") == "PASS", "canonical N=500 scale gate missing")
    check(gates.get("canonical_n1000_projection") == "PASS", "canonical N=1,000 projection gate missing")
    check(gates.get("canonical_scale_trend") == "PASS", "final canonical scale gate missing")
    for name, verdict in gates.items():
        if name in expected_na:
            check(verdict == expected_na[name], f"not-applicable gate contract mismatch: {name}")
        else:
            check(verdict == "PASS", f"release gate is not PASS: {name}={verdict}")
    check(all(name in gates for name in expected_na), "required canonical-only NOT_APPLICABLE gates missing")
    groups.append("release_identity_and_all_applicable_gates")

    tracked = repo / "manifests/canonical-cohort-1000-v1"
    tracked_sums = r.verify_tracked_inventory(tracked)
    tracked_actual = {path.name for path in tracked.iterdir() if path.is_file()}
    check(tracked_actual == set(tracked_sums) | {"SHA256SUMS"}, "tracked inventory has missing/unlisted files")
    tracked_mapping = {
        "cohort-1000.tsv": "manifests/cohort-1000.tsv",
        "assemblies.tsv": "manifests/assemblies.tsv",
        "checksums.tsv": "manifests/checksums.tsv",
        "state.tsv": "manifests/state.tsv",
        "object_refs.tsv": "manifests/object_refs.tsv",
        "batch_metrics.tsv": "manifests/batch_metrics.tsv",
        "release.json": "release.json",
        "external_SHA256SUMS": "SHA256SUMS",
    }
    for tracked_name, external_name in tracked_mapping.items():
        check((tracked / tracked_name).read_bytes() == (external / external_name).read_bytes(), f"tracked/external handoff mismatch: {tracked_name}")
    contig_parts = sorted(tracked.glob("contigs.part-*.tsv.gz"))
    check(contig_parts, "tracked compressed contig parts missing")
    digest = hashlib.sha256()
    row_count = 0
    with (external / "manifests/contigs.tsv").open("rb") as handle:
        expected_header = handle.readline()
    digest.update(expected_header)
    for part in contig_parts:
        with gzip.open(part, "rb") as handle:
            check(handle.readline() == expected_header, f"tracked contig part header mismatch: {part.name}")
            for line in handle:
                digest.update(line)
                row_count += 1
    check(row_count == release["manifests"]["contigs.tsv"]["rows"], "tracked contig part row accounting mismatch")
    check(digest.hexdigest() == release["manifests"]["contigs.tsv"]["sha256"], "tracked contig parts do not reconstruct external manifest")
    groups.append("tracked_split_compact_manifest_inventory")

    cohort_external = external / "manifests/cohort-1000.tsv"
    cohort_tracked = repo / "manifests/pilot-cohorts-v1/cohort-1000.tsv"
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
    check(len(assemblies) == len(states) == len(refs) == r.COHORT_ROWS, "assembly/state/ref cardinality mismatch")
    check(len(checksums) == r.COHORT_ROWS * 8, "checksum row cardinality mismatch")
    check(len(batches) == 100, "bounded batch row cardinality mismatch")
    expected_manifest_rows = {
        "cohort-1000.tsv": r.COHORT_ROWS,
        "assemblies.tsv": r.COHORT_ROWS,
        "contigs.tsv": len(contigs),
        "checksums.tsv": r.COHORT_ROWS * 8,
        "state.tsv": r.COHORT_ROWS,
        "object_refs.tsv": r.COHORT_ROWS,
        "batch_metrics.tsv": len(batches),
    }
    for name, rows in expected_manifest_rows.items():
        path = external / "manifests" / name
        recorded = release.get("manifests", {}).get(name, {})
        check(recorded.get("rows") == rows, f"release manifest row count mismatch: {name}")
        check(recorded.get("bytes") == path.stat().st_size, f"release manifest byte count mismatch: {name}")
        check(recorded.get("sha256") == r.sha_file(path), f"release manifest SHA mismatch: {name}")
    check(release.get("counts", {}).get("reused_predecessor_objects") == r.PREDECESSOR_ROWS, "predecessor reuse count mismatch")
    check(release.get("counts", {}).get("new_objects") == r.NEW_ROWS, "new object count mismatch")
    check(release.get("counts", {}).get("validated") == r.COHORT_ROWS, "validated count mismatch")
    check([row["accession"] for row in assemblies] == accessions, "assembly order mismatch")
    check([row["accession"] for row in states] == accessions, "state order mismatch")
    check([row["accession"] for row in refs] == accessions, "reference order mismatch")
    check(all(row["terminal_state"] in ("VALIDATED", "QUARANTINED") for row in states), "non-terminal assembly state")
    check(all(row["terminal_state"] == "VALIDATED" for row in states), "genome object quarantine present")
    check(sum(int(row["contig_count"]) for row in assemblies) == len(contigs), "contig accounting mismatch")
    check(len({row["pansn_sequence_name"] for row in contigs}) == len(contigs), "PanSN collision")
    groups.append("complete_terminal_row_accounting")

    checksum_by_key = {(row["accession"], row["artifact_role"]): row for row in checksums}
    check(len(checksum_by_key) == r.COHORT_ROWS * 8, "duplicate checksum role rows")
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
        root = _root_for_ref(external, ref)
        expected_root, expected_storage_id, expected_storage_root, _ = r._storage_for(
            order, accession, external, inputs
        )
        check(
            root == expected_root
            and ref["storage_release_id"] == expected_storage_id
            and ref["storage_root"] == expected_storage_root,
            f"immutable object storage reference mismatch: {accession}",
        )
        source = root / ref["source_object_relpath"]
        canonical = root / ref["canonical_object_relpath"]
        p.verify_sha_inventory(source)
        p.verify_sha_inventory(canonical)
        check(r.sha_file(source / "SHA256SUMS") == ref["source_inventory_sha256"], f"source object inventory reference mismatch: {accession}")
        check(r.sha_file(canonical / "SHA256SUMS") == ref["canonical_inventory_sha256"], f"canonical object inventory reference mismatch: {accession}")
        if order <= r.PREDECESSOR_ROWS:
            check(ref["reuse_status"] == "REUSED_PREDECESSOR_BY_DIGEST" and ref["predecessor_digest_match"] == "PASS", f"predecessor reuse status mismatch: {accession}")
            old = old_assemblies.get(accession)
            check(old is not None and all(assembly[field] == old[field] for field in reuse_digest_fields), f"N=500 checksum reuse mismatch: {accession}")
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
        with zipfile.ZipFile(source / "package.zip") as archive:
            reports = p.load_sequence_reports(archive, package_validation["sequence_report_member"])
            recompute_contigs = []
            for row in object_contigs:
                report = p.report_for_token(reports, row["source_contig_id_display"])
                recompute_contigs.append(dict(row, sequence_report_aliases=sorted(p.report_aliases(report))))
            try:
                recomputed_aliases, recomputed_annotation = p.validate_gff_aliases(
                    archive, package_validation["gff_member"], recompute_contigs, accession
                )
                annotation_failure = None
            except p.AnnotationValidationError as exc:
                recomputed_aliases, recomputed_annotation, annotation_failure = [], {}, str(exc)
        if annotation["status"] == "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW":
            annotation_quarantines += 1
            check(not aliases and annotation.get("failure_reason"), f"quarantined annotation exposed aliases: {accession}")
            check(annotation_failure == annotation["failure_reason"], f"annotation quarantine is not reproducible: {accession}")
        else:
            check(annotation_failure is None, f"published annotation alias view fails recomputation: {accession}")
            comparable_fields = [field for field in aliases[0] if field != "row_sha256"] if aliases else []
            observed_aliases = [{field: row[field] for field in comparable_fields} for row in aliases]
            expected_aliases = [{field: str(row[field]) for field in comparable_fields} for row in recomputed_aliases]
            check(observed_aliases == expected_aliases, f"annotation alias recomputation mismatch: {accession}")
            check(annotation == recomputed_annotation, f"annotation coordinate-policy summary mismatch: {accession}")
        for role in ("source_package", "source_manifest", "canonical_bgzf", "fai", "gzi", "contig_crosswalk", "annotation_aliases", "canonical_manifest"):
            checksum = checksum_by_key.get((accession, role))
            check(checksum is not None, f"missing checksum role: {accession}:{role}")
            check(checksum["storage_release_id"] == ref["storage_release_id"] and checksum["storage_root"] == ref["storage_root"], f"checksum storage mismatch: {accession}:{role}")
            path = root / checksum["relative_path"]
            check(path.is_file() and not path.is_symlink(), f"checksum artifact missing/symlink: {path}")
            check(path.stat().st_size == int(checksum["bytes"]) and r.sha_file(path) == checksum["sha256"], f"checksum artifact mismatch: {path}")
        total_bases += sum(length for _, length, _ in expected)
    groups.append("all_1000_source_bgzf_index_crosswalk_annotation_semantics")
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
    check(not request_accessions.intersection(accessions[:r.PREDECESSOR_ROWS]), "validated predecessor object was re-downloaded")
    check(request_accessions.issubset(set(accessions[r.PREDECESSOR_ROWS:])), "download outside frozen new 500 rows")
    groups.append("bounded_retry_kill_restart_no_predecessor_redownload")

    resources = [json.loads(line) for line in (external / "resources.jsonl").read_text().splitlines()]
    check(resources and all(row.get("verdict") == "PASS" and all(row.get("checks", {}).values()) for row in resources), "resource preflight failure")
    summary = json.loads((external / "resource_summary.json").read_text())
    check(summary.get("verdict") == "PASS" and all(summary.get("checks", {}).values()), "resource summary failure")
    check(summary["peak_rss_fraction"] <= 0.70, "RAM gate failure")
    check(summary["process_swap_events"] == 0 and summary["system_swap_growth_bytes"] == 0, "swap gate failure")
    check(int(summary["finish"]["swap_free_bytes"]) >= int(summary["start"]["swap_free_bytes"]), "successful retry sampled swap growth")
    resource_failures = [json.loads(line) for line in (external / "failures.jsonl").read_text().splitlines() if line.strip()]
    check(all(row.get("event") == "RUN_FAILED" for row in resource_failures), "unexpected failure-ledger event class")
    check(all("system_swap_growth_zero" in row.get("message", "") for row in resource_failures), "failed attempt was not a recorded resource-only refusal")
    scale = json.loads((external / "scale_trend.json").read_text())
    check(scale.get("verdict") == "PASS" and all(scale.get("checks", {}).values()), "canonical scale trend failure")
    check(summary.get("scale_trend") == scale, "resource summary/scale trend mismatch")
    check(release.get("scale_trend_sha256") == r.sha_file(external / "scale_trend.json"), "scale trend checksum mismatch")
    scale_observation = summary.get("scale_observation", {})
    check(scale_observation.get("source") in ("CURRENT_INVOCATION", "CHECKSUM_VALIDATED_PRIOR_ATTEMPT"), "scale observation provenance missing")
    recomputed_scale = r.compute_scale_trend(
        repo, inputs, assemblies, float(scale_observation["wall_seconds"]), int(scale_observation["peak_rss_bytes"]),
        int(scale_observation["stage_bytes"]), int(scale_observation["stage_files"]),
        p.Allocations(**summary["allocations"]), float(scale["pinned_n1000_projection"]["time_allocation_seconds"]),
    )
    check(recomputed_scale == scale, "canonical final scale trend recomputation mismatch")
    check(scale["authorized_ceiling"] == 1000, "final scale authorization ceiling mismatch")
    check(scale["projection_beyond_ceiling"] == "NOT_AUTHORIZED_NOT_COMPUTED", "unapproved beyond-1,000 projection present")
    check(scale["time_exponent"]["empirical_upper_bound"] <= 1.3, "canonical time exponent upper bound failure")
    check(all(abs(value) <= 0.25 for value in scale["last_two_rung_per_base_slopes"]["relative_changes"].values()), "canonical per-base slope change failure")
    check(all(scale["pinned_n1000_projection"]["checks"].values()), "consumed N=1,000 upper-95 projection failure")
    projection = summary.get("disk_projection", {})
    check(projection.get("modeled_upper95_peak_bytes", 0) <= projection.get("configured_upper95_peak_bytes", -1), "upper-95 projection was underallocated")
    check(summary["measured_release_stage_peak_bytes"] <= projection.get("modeled_upper95_peak_bytes", -1), "unexpected disk growth exceeded upper-95 projection")
    if scale_observation["source"] == "CHECKSUM_VALIDATED_PRIOR_ATTEMPT":
        failures = [json.loads(line) for line in (external / "failures.jsonl").read_text().splitlines() if line.strip()]
        check(any(row.get("event") == "RUN_FAILED" and "system_swap_growth_zero" in row.get("message", "") for row in failures), "preserved scale observation lacks recorded resource-only failure")
    partial_stage_bytes = [int(row["stage_partial_bytes_finish"]) for row in batches]
    check(partial_stage_bytes[-1] == 0, "final batch ended with unfinished partial bytes")
    predecessor_only_batches = r.PREDECESSOR_ROWS // 10
    check(all(value == 0 for value in partial_stage_bytes[predecessor_only_batches:]), "partial stage survived its deterministic resume batch")
    check(max(partial_stage_bytes) <= summary["allocations"]["unfinished_write_bytes"], "partial stage exceeded unfinished-write reservation")
    check(all(int(row["durable_free_bytes"]) >= 2_000_000_000_000 for row in batches), "batch durable floor failure")
    check(all(int(row["durable_free_inodes"]) >= 1_000_000 for row in batches), "batch durable inode floor failure")
    check(all(int(row["scratch_free_bytes"]) >= 4_000_000_000_000 for row in batches), "batch scratch preflight floor failure")
    check(all(int(row["scratch_free_inodes"]) >= 5_000_000 for row in batches), "batch scratch inode floor failure")
    check([int(row["cumulative_transfer_bytes"]) for row in batches] == sorted(int(row["cumulative_transfer_bytes"]) for row in batches), "non-monotonic cumulative transfer metrics")
    check(int(batches[-1]["cumulative_transfer_bytes"]) == sum(int(row["source_package_bytes"]) for row in assemblies[r.PREDECESSOR_ROWS:]), "cumulative transfer accounting mismatch")
    check(int(batches[-1]["cumulative_canonical_bgzf_bytes"]) == sum(int(row["canonical_bgzf_bytes"]) for row in assemblies[r.PREDECESSOR_ROWS:]), "cumulative canonical byte accounting mismatch")
    check(sum(int(row["restart_events"]) for row in batches) >= 4, "batch restart metrics incomplete")
    check(sum(int(row["partial_bytes_observed"]) for row in batches) >= 262_144, "partial-byte accounting incomplete")
    check(all(float(row["wall_seconds"]) >= 0 and float(row["cpu_seconds"]) >= 0 for row in batches), "negative batch wall/CPU metric")
    check(all(int(row["peak_rss_bytes"]) > 0 and int(row["stage_final_bytes_finish"]) >= 0 for row in batches), "batch RAM/final-byte metric missing")
    groups.append("live_resources_metrics_and_allocation_thresholds")

    root_finish = r.verify_root_inputs(repo)
    check(root_finish == release["root_inputs_start"] == release["root_inputs_finish"], "root hashes changed")
    cap_evidence = json.loads((external / "global_cap_evidence.json").read_text())
    check(cap_evidence.get("verdict") == "PASS", "global cap evidence verdict mismatch")
    check(cap_evidence.get("start", {}).get("distinct_exact_assembly_revisions") == r.PREDECESSOR_ROWS, "start union was not exact N=500")
    check(cap_evidence.get("projected_finish_distinct_exact_assembly_revisions") == r.COHORT_ROWS, "projected finish union mismatch")
    check(cap_evidence.get("projected_finish_subset_of_frozen_collection") is True, "projected union not frozen subset")
    allowed = p._load_allowed_collection(repo)
    global_finish = r.audit_global_release_cap(external.parents[2], allowed)
    check(global_finish["distinct_exact_assembly_revisions"] == r.COHORT_ROWS, "finish union cardinality is not 1,000")
    check(set(global_finish["accessions"]) == set(accessions), "finish union differs from exact frozen cohort")
    groups.append("root_hashes_and_global_exact_cap_finish")

    forbidden = [str(path.relative_to(external)) for path in external.rglob("*") if path.is_file() and path.suffix.lower() in (".fa", ".fna", ".fasta")]
    check(not forbidden, f"plain FASTA in release: {forbidden}")
    check(not any(path.name.startswith(".stage.") for path in external.rglob("*")), "nested staging path in final release")
    check(not any(path.is_symlink() for path in external.rglob("*")), "symlink in final release")
    inventory_paths = {
        line.split("  ", 1)[1] for line in (external / "SHA256SUMS").read_text().splitlines()
    }
    actual_paths = {
        str(path.relative_to(external)) for path in external.rglob("*") if path.is_file()
    } - {"SHA256SUMS", "COMPLETE"}
    check(inventory_paths == actual_paths, "external inventory has missing or unlisted files")
    owned_roots = [
        repo / "manifests/canonical-cohort-1000-v1",
        repo / "artifacts/canonical_cohort_1000",
        repo / "reports/canonical_cohort_1000.md",
    ]
    owned_files = [
        path for root in owned_roots
        for path in ([root] if root.is_file() else (list(root.rglob("*")) if root.exists() else []))
        if path.is_file()
    ]
    large_git_payloads = [str(path.relative_to(repo)) for path in owned_files if path.stat().st_size > 10 * 1024 * 1024]
    forbidden_git_payloads = [
        str(path.relative_to(repo)) for path in owned_files
        if path.name.endswith((".fa.gz", ".fna.gz", ".fasta.gz", ".fai", ".gzi", ".zip", ".gff"))
    ]
    check(not large_git_payloads, f"task-owned Git artifact exceeds 10 MiB: {large_git_payloads}")
    check(not forbidden_git_payloads, f"sequence/index payload in task-owned Git paths: {forbidden_git_payloads}")
    groups.append("atomic_inventory_no_partial_plain_fasta_or_symlink")

    return {
        "schema": "canonical-cohort-1000-validation-v1", "verdict": "PASS",
        "release_id": release["release_id"], "external_release_path": str(external),
        "check_groups": groups, "check_group_count": len(groups), "assembly_rows": len(assemblies),
        "contig_rows": len(contigs), "checksum_rows": len(checksums), "object_reference_rows": len(refs),
        "batch_rows": len(batches), "total_bases": total_bases, "annotation_view_quarantines": annotation_quarantines,
        "reused_predecessor_assemblies": r.PREDECESSOR_ROWS, "new_assemblies": r.NEW_ROWS,
        "global_distinct_exact_assembly_revisions": global_finish["distinct_exact_assembly_revisions"],
        "root_inputs": root_finish, "large_payloads_in_git": len(large_git_payloads),
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
