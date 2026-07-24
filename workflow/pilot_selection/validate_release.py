#!/usr/bin/env python3
"""Independent semantic/checksum validator for pilot-cohorts-v1."""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.pilot_selection import selection as s


class ValidationError(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def validate(repo: Path, external: Path) -> dict[str, Any]:
    checks: list[str] = []
    assemblies, stage_b, input_gate = s.verify_inputs(repo, s.DEFAULT_RELEASES_ROOT)
    checks.append("immutable_inputs_predecessor_release_ids_inventories_and_applicable_gates")
    s.verify_inventory(external)
    checks.append("external_COMPLETE_and_full_SHA256SUMS")

    release = json.loads((external / "release.json").read_text(encoding="utf-8"))
    expected_release_id = s.release_id(s.sha_file(Path(s.__file__)))
    check(release.get("release_id") == expected_release_id, "release ID is not code/input/policy pinned")
    check(release.get("verdict") == "PASS" and release.get("immutable") is True, "release verdict/immutability failure")
    check(release.get("external_release_path") == str(external), "external release path mismatch")
    check(release.get("collection_release_id") == s.COLLECTION_RELEASE_ID, "collection release ID mismatch")
    check(release.get("validated_stage_b_release_id") == s.ACQUISITION_RELEASE_ID, "Stage-B acquisition release ID mismatch")
    check(release.get("engineering_control_source_semantics_release_id") == s.SOURCE_SEMANTICS_RELEASE_ID, "source-semantics pointer mismatch")
    check(release.get("source_semantics_extraction_eligibility") == "EXTRACTION_BLOCKED", "blocked extraction status was not preserved")
    check(release.get("source_semantics_consumer_action") == "REJECT", "source semantics is not fail-closed")
    allowed_na = {
        "bgzf_index_name_roundtrip": "NOT_APPLICABLE_METADATA_ONLY_SELECTION",
        "source_coordinate_policy": "NOT_APPLICABLE_HOST_SELECTION_EXTRACTION_BLOCKED",
        "scale_trend": "NOT_APPLICABLE_METADATA_ONLY_SELECTION",
    }
    for name, verdict in release.get("applicable_gates", {}).items():
        check(verdict == "PASS" or allowed_na.get(name) == verdict, f"applicable release gate is not PASS/declared NA: {name}={verdict}")
    checks.append("release_contract_and_applicable_gate_verdicts")

    input_manifest = json.loads((external / "input_manifest.json").read_text(encoding="utf-8"))
    policy = json.loads((external / "selection_policy.json").read_text(encoding="utf-8"))
    check(s.sha_file(external / "input_manifest.json") == release["input_manifest_sha256"], "input manifest digest mismatch")
    check(s.sha_file(external / "selection_policy.json") == release["selection_policy_sha256"], "selection policy digest mismatch")
    check(input_manifest == input_gate, "published exact input manifest differs from fresh gate")
    check(policy.get("phage_blind") is True and policy.get("policy_frozen") is True, "selection policy is not frozen/phage blind")
    check(policy.get("allowed_selection_fields") == list(s.ALLOWED_SELECTION_FIELDS), "allowed selection field list drift")
    check(policy.get("forbidden_selection_fields") == list(s.FORBIDDEN_SELECTION_FIELDS), "forbidden selection field list drift")
    check(policy.get("seed_literal") == s.MAIN_SEED and policy.get("algorithm") == s.SELECTION_ALGORITHM, "seed/algorithm drift")
    checks.append("frozen_phage_blind_policy_seed_and_exact_input_manifest")

    check(s.sha_file(external / "inputs/collection-assemblies.tsv.gz") == s.COLLECTION_ASSEMBLIES_GZ_SHA256, "external exact input assemblies changed")
    check(s.sha_file(external / "inputs/stage_b_10.tsv") == s.STAGE_B_SHA256, "external exact Stage-B input changed")
    check((external / "inputs/stage_b_10.tsv").read_bytes() == (repo / "manifests/collection-v1/stage_b_10.tsv").read_bytes(), "external Stage-B bytes differ")
    check(s.sha_file(external / "inputs/collection-release.json") == s.COLLECTION_RELEASE_JSON_SHA256, "external collection release input changed")
    check(s.sha_file(external / "inputs/canonical-cohort-010-release.json") == s.ACQUISITION_RELEASE_JSON_SHA256, "external acquisition release input changed")
    check(s.sha_file(external / "inputs/source-semantics-release.json") == s.SOURCE_SEMANTICS_RELEASE_JSON_SHA256, "external source semantics input changed")
    checks.append("exact_digest_pinned_input_copies")

    frame = s.read_tsv(external / "manifests/frame.tsv", s.FRAME_FIELDS, verify_hashes=True)
    controls = s.read_tsv(external / "manifests/engineering-controls.tsv", s.CONTROL_FIELDS, verify_hashes=True)
    rungs = {n: s.read_tsv(external / f"manifests/cohort-{n:04d}.tsv", s.COHORT_FIELDS, verify_hashes=True) for n in s.RUNGS}
    check(len(frame) == 26077, "full frame row count is not 26,077")
    check(len({row["assembly_id"] for row in frame}) == len({row["exact_assembly_accession_version"] for row in frame}) == 26077, "frame exact identity duplication")
    frame_counts = Counter(row["selection_stratum"] for row in frame)
    check(frame_counts == Counter({"main_phage_blind_srs": 25353, "stage_b_certainty": 10, "terminal_suppressed_ineligible": 714}), "frame stratum accounting mismatch")
    check(all(row["frame_disposition"] == "INELIGIBLE_TERMINAL_SUPPRESSED" for row in frame if row["assembly_status"] == "suppressed"), "suppressed row eligibility failure")
    checks.append("full_26077_mapping_exact_versions_and_terminal_status_accounting")

    expected = s.select_cohort(assemblies, stage_b)
    expected_order = [row["assembly_id"] for row in expected.cohort]
    observed_order = [row["assembly_id"] for row in rungs[1000]]
    check(observed_order == expected_order, "published order differs from independent deterministic selection")
    for n in s.RUNGS:
        check(len(rungs[n]) == n, f"rung {n} cardinality mismatch")
        check([row["assembly_id"] for row in rungs[n]] == observed_order[:n], f"rung {n} is not exact prefix")
        check([int(row["cohort_order"]) for row in rungs[n]] == list(range(1, n + 1)), f"rung {n} order numbering mismatch")
        check(all(int(row["rung_n"]) == n for row in rungs[n]), f"rung {n} label mismatch")
    check([row["assembly_id"] for row in rungs[10]] == [row["assembly_id"] for row in stage_b], "N=10 is not exact Stage-B order")
    check(len(set(observed_order)) == 1000, "duplicate selected assembly ID")
    checks.append("exact_nested_10_100_250_500_1000_and_stage_b_reuse")

    main_count = 25353
    for n in s.RUNGS:
        for row in rungs[n]:
            if row["selection_stratum"] == "stage_b_certainty":
                check(row["inclusion_probability"] == "1/1" and row["inference_weight"] == "1/1", "certainty design weight mismatch")
            else:
                check(row["inclusion_probability"] == f"{n - 10}/{main_count}", "main-stratum inclusion probability mismatch")
                check(row["inference_weight"] == f"{main_count}/{n - 10}", "main-stratum inverse probability mismatch")
    for row in frame:
        for n in s.RUNGS:
            expected_probability = (
                "1/1" if row["selection_stratum"] == "stage_b_certainty" else
                "0/1" if row["selection_stratum"] == "terminal_suppressed_ineligible" else
                f"{n - 10}/{main_count}"
            )
            check(row[f"inclusion_probability_n{n}"] == expected_probability, "frame inclusion probability mismatch")
    checks.append("exact_strata_inclusion_probabilities_and_inference_weights")

    check(controls, "engineering control manifest empty")
    check(all(row["selection_effect"] == "NONE_POST_SELECTION_LABEL" for row in controls), "engineering control affected selection")
    check(any(row["control_class"] == "accession_exact_revision_ambiguity" for row in controls), "accession ambiguity control absent")
    check(any(row["control_class"].startswith("assembly_size") for row in controls), "assembly size control absent")
    check(any(row["control_class"].startswith("assembly_contiguity") for row in controls), "assembly contiguity control absent")
    check(any(row["control_class"] == "unsafe_source_contig_id" for row in controls), "unsafe-ID fixture absent")
    check(any(row["control_class"] == "source_contig_role" for row in controls), "source-contig role control absent")
    check(any(row["control_class"] == "prophage_interval_edge_extreme" for row in controls), "future interval edge/extreme fixture absent")
    check(all(row["inference_disposition"].startswith(("EXCLUDE_", "USE_ONCE_WITH_FROZEN_DESIGN_WEIGHT")) for row in controls), "control inference disposition absent")
    check(all(row["activation_status"] == "BLOCKED_EXTRACTION_SEMANTICS" for row in controls if row["control_class"].startswith("prophage_interval")), "interval control activated despite blocked semantics")
    checks.append("separate_bounded_engineering_controls_and_inference_disposition")

    for name, metadata in release["manifests"].items():
        path = external / "manifests" / name
        check(path.is_file() and s.sha_file(path) == metadata["sha256"], f"release manifest checksum mismatch: {name}")
        check(path.stat().st_size == metadata["bytes"], f"release manifest byte count mismatch: {name}")
        check(sum(1 for _ in path.open(encoding="utf-8")) - 1 == metadata["rows"], f"release manifest row count mismatch: {name}")
    checks.append("manifest_row_byte_sha256_contract")

    restart = json.loads((external / "restart_evidence.json").read_text(encoding="utf-8"))
    check(restart.get("verdict") == "PASS", "injected restart verdict failure")
    check(all(restart.get(key) is True for key in (
        "injected_kill_recorded", "interrupted_stage_discarded",
        "checksum_validated_selection_unit_reused", "no_partial_final_before_restart",
    )), "injected restart evidence incomplete")
    state = (external / "state.jsonl").read_text(encoding="utf-8")
    check("INJECTED_KILL_BEFORE_COMPLETE" in state and "INTERRUPTED_PUBLICATION_STAGE_DISCARDED" in state and "SELECTION_UNIT_REUSED_CHECKSUM_VALIDATED" in state, "restart state events absent")
    checks.append("injected_kill_checksum_validated_resume_and_no_partial_publication")

    deterministic = json.loads((external / "deterministic_rerun.json").read_text(encoding="utf-8"))
    check(deterministic.get("verdict") == "PASS", "deterministic rerun verdict failure")
    for relative, digest in deterministic.get("manifest_sha256", {}).items():
        check(s.sha_file(external / relative) == digest, f"deterministic manifest digest mismatch: {relative}")
    checks.append("byte_deterministic_manifest_rerun")

    resources = json.loads((external / "resource_summary.json").read_text(encoding="utf-8"))
    check(resources.get("verdict") == "PASS" and all(resources.get("checks", {}).values()), "resource summary failure")
    check(resources["peak_rss_fraction"] <= 0.70 and resources["swap_growth_bytes"] == 0, "RSS/swap threshold failure")
    check(resources["actual_preseal_stage_bytes"] <= resources["upper_95_predicted_durable_peak_bytes"], "actual stage exceeds predicted upper-95 peak")
    check(resources["resource_record_count"] >= 2, "per-stage resource preflights absent")
    checks.append("mount_owner_write_free_space_inode_ram_swap_and_reservation_gates")

    cap = json.loads((external / "global_cap_evidence.json").read_text(encoding="utf-8"))
    check(cap.get("verdict") == "PASS" and cap.get("distinct_exact_assembly_revisions") == 10, "global union count mismatch")
    check(set(cap["accessions"]).issubset({row["exact_assembly_accession_version"] for row in rungs[1000]}), "global sequence union outside frozen cohort")
    check(release["counts"]["distinct_sequence_bearing_assemblies_created"] == 0 and release["sequence_downloads"] == 0 and release["biological_analyses"] == 0, "selection release claims payload/analysis")
    checks.append("global_distinct_assembly_cap_frozen_subset_and_zero_new_payload")

    forbidden_suffixes = (".fa", ".fna", ".fasta", ".gff", ".gff3", ".bgz", ".fai", ".gzi")
    forbidden = [str(path.relative_to(external)) for path in external.rglob("*") if path.is_file() and path.name.lower().endswith(forbidden_suffixes)]
    check(not forbidden, f"sequence/index-bearing output in metadata release: {forbidden}")
    checks.append("metadata_only_file_inventory")

    tracked = repo / "manifests/pilot-cohorts-v1"
    s._read_checksum_map(tracked)
    check((tracked / "release.json").read_bytes() == (external / "release.json").read_bytes(), "tracked release JSON differs")
    for n in s.RUNGS:
        check((tracked / f"cohort-{n:04d}.tsv").read_bytes() == (external / f"manifests/cohort-{n:04d}.tsv").read_bytes(), f"tracked rung {n} bytes differ")
    with gzip.open(tracked / "frame.tsv.gz", "rb") as handle:
        check(handle.read() == (external / "manifests/frame.tsv").read_bytes(), "tracked frame does not decompress to external frame")
    check(all(path.stat().st_size <= 10 * 1024 * 1024 for path in tracked.rglob("*") if path.is_file()), "tracked artifact exceeds 10 MiB")
    checks.append("tracked_handoff_checksums_exact_bytes_and_compactness")

    roots_finish = s.verify_root_inputs(repo)
    check(release["root_inputs_start"] == roots_finish == release["root_inputs_finish"], "root input start/finish checksum mismatch")
    checks.append("root_input_sha256_start_and_finish")

    return {
        "schema": "pilot-cohorts-v1-validation", "verdict": "PASS",
        "release_id": release["release_id"], "external_release_path": str(external),
        "checks": checks, "check_count": len(checks),
        "frame_rows": len(frame), "eligible_rows": 25363,
        "rung_rows": {str(n): len(rungs[n]) for n in s.RUNGS},
        "engineering_control_rows": len(controls),
        "global_distinct_sequence_bearing_assemblies": cap["distinct_exact_assembly_revisions"],
        "sequence_downloads": 0, "biological_analyses": 0,
        "root_inputs": roots_finish,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--external-release", required=True)
    parser.add_argument("--output")
    args = parser.parse_args()
    try:
        result = validate(Path(args.repo_root).resolve(), Path(args.external_release).resolve())
        data = s.canonical_json(result)
        if args.output:
            output = Path(args.output)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(data)
        print(data.decode(), end="")
        return 0
    except (s.GateError, ValidationError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
