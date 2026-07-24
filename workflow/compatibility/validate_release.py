#!/usr/bin/env python3
"""Independent semantic/checksum validator for consumer compatibility v1."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Sequence

from .compatibility import (
    EXPECTED_ASSEMBLIES, GateError, PASS, PREDECESSOR_RELEASE_ID,
    PREDECESSOR_RELEASE_JSON_SHA256, require_cleaned, sha256_file,
    verify_inventory, verify_predecessor, verify_root_hashes,
)

REQUIRED_GATES = {
    "bgzip", "samtools-faidx", "impg-syng-build", "impg-query", "impg-map",
    "mash", "rapidnj", "skani", "quast", "gffread", "pggb", "seqwish",
    "smoothxg", "odgi", "vg", "prodigal", "mmseqs2", "hmmer", "mcl",
}


def validate_release(repo: Path, release: Path) -> dict[str, Any]:
    groups: dict[str, str] = {}
    if not release.is_dir() or not (release / "COMPLETE").is_file():
        raise GateError("compatibility release absent or incomplete")
    inventory_rows = verify_inventory(release, release / "SHA256SUMS", exact=True)
    complete = json.loads((release / "COMPLETE").read_text())
    if complete.get("verdict") != PASS or complete.get("sha256sums_sha256") != sha256_file(release / "SHA256SUMS"):
        raise GateError("COMPLETE does not seal SHA256SUMS with PASS")
    groups["inventory_complete_atomic"] = PASS

    doc = json.loads((release / "release.json").read_text())
    if doc.get("schema_version") != "consumer-compatibility-release-v1" or doc.get("verdict") != PASS:
        raise GateError("release schema/verdict mismatch")
    if doc.get("release_id") != release.name or doc.get("external_path") != str(release):
        raise GateError("release identity/path mismatch")
    if doc.get("predecessor_release_id") != PREDECESSOR_RELEASE_ID:
        raise GateError("predecessor release ID changed")
    if doc.get("predecessor_release_json_sha256") != PREDECESSOR_RELEASE_JSON_SHA256:
        raise GateError("predecessor manifest digest changed")
    counts = doc.get("counts", {})
    if counts.get("assemblies") != 10 or counts.get("contigs") != 1223 or counts.get("bases") != 51731662:
        raise GateError("release row/base/contig accounting mismatch")
    if counts.get("distinct_sequence_bearing_assemblies") != 10 or counts.get("global_distinct_assembly_cap") != 1000:
        raise GateError("global distinct-assembly union/cap mismatch")
    if doc.get("cohort_order") != EXPECTED_ASSEMBLIES:
        raise GateError("immutable cohort order mismatch")
    for name, verdict in doc.get("applicable_gates", {}).items():
        if verdict != PASS and not (name == "scale_trend" and verdict == "NOT_APPLICABLE_NON_SCALE_BEARING_COMPATIBILITY"):
            raise GateError(f"applicable release gate {name} is {verdict}")
    groups["release_identity_rows_cap"] = PASS

    input_doc = json.loads((release / "input_manifest.json").read_text())
    if input_doc.get("assembly_order") != EXPECTED_ASSEMBLIES or input_doc.get("assembly_count") != 10:
        raise GateError("input manifest order/count mismatch")
    if input_doc.get("predecessor_manifest_bytes_sha256") != PREDECESSOR_RELEASE_JSON_SHA256:
        raise GateError("input manifest predecessor checksum mismatch")
    _, _, predecessor = verify_predecessor(repo)
    if input_doc.get("predecessor", {}).get("external_inventory_sha256") != predecessor["external_inventory_sha256"]:
        raise GateError("input manifest external predecessor checksum mismatch")
    groups["predecessor_all_objects_indexes"] = PASS

    gates: dict[str, dict[str, Any]] = {}
    for path in sorted((release / "gates").glob("*.json")):
        gate = json.loads(path.read_text())
        gate_id = gate.get("gate_id")
        if path.name != f"{gate_id}.json" or gate_id in gates:
            raise GateError(f"duplicate/misnamed machine gate: {path}")
        if gate.get("verdict") != PASS:
            raise GateError(f"required consumer {gate_id} is not unqualified PASS")
        if not all(gate.get(k) for k in ("tool", "input_form", "view_contract", "invocation", "output_name_behavior")):
            raise GateError(f"consumer gate lacks exact contract fields: {gate_id}")
        if any(v not in {PASS, True, "PASS"} for v in gate.get("checks", {}).values()):
            raise GateError(f"consumer gate has failed/conditional check: {gate_id}")
        gates[gate_id] = gate
    if set(gates) != REQUIRED_GATES:
        raise GateError(f"required consumer gate set mismatch missing={sorted(REQUIRED_GATES-set(gates))} extra={sorted(set(gates)-REQUIRED_GATES)}")
    if counts.get("required_consumer_gates") != len(REQUIRED_GATES) or counts.get("pass_consumer_gates") != len(REQUIRED_GATES):
        raise GateError("consumer gate count mismatch")
    groups["all_required_consumers_pass"] = PASS

    with (release / "consumers.tsv").open(newline="") as fh:
        consumer_rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(consumer_rows) != len(REQUIRED_GATES) or {r["consumer_id"] for r in consumer_rows} != REQUIRED_GATES:
        raise GateError("consumers.tsv does not account for required gates")
    if any(r["verdict"] != PASS for r in consumer_rows):
        raise GateError("consumers.tsv contains a non-PASS required consumer")
    groups["consumer_table_contracts"] = PASS

    tools = json.loads((release / "tool_versions.json").read_text())
    if tools.get("host_environment_lock", {}).get("sha256") != doc.get("host_environment_lock_sha256"):
        raise GateError("host environment lock pin mismatch")
    if tools.get("host_environment_lock", {}).get("package_sha256_inventory_sha256") != doc.get("host_package_sha256_inventory_sha256"):
        raise GateError("host package SHA-256 inventory pin mismatch")
    if tools.get("graph_environment_lock", {}).get("sha256") != doc.get("graph_environment_lock_sha256"):
        raise GateError("graph environment lock pin mismatch")
    if tools.get("graph_environment_lock", {}).get("package_sha256_inventory_sha256") != doc.get("graph_package_sha256_inventory_sha256"):
        raise GateError("graph package SHA-256 inventory pin mismatch")
    if set(tools.get("tools", {})) != {
        "bgzip", "samtools", "impg", "mash", "rapidnj", "skani", "quast", "gffread",
        "prodigal", "mmseqs", "hmmbuild", "hmmsearch", "mcl", "pggb", "wfmash",
        "seqwish", "smoothxg", "odgi", "vg",
    }:
        raise GateError("tool identity set mismatch")
    for name, identity in tools["tools"].items():
        if (len(identity.get("sha256", "")) != 64 or not identity.get("version_output")
                or not identity.get("help_output") or "help_argv" not in identity):
            raise GateError(f"tool lacks exact digest/version/help capture: {name}")
    groups["tool_version_digest_lock"] = PASS

    views = json.loads((release / "view_contracts.json").read_text())
    if views.get("verdict") != PASS or not views.get("cleanup_required"):
        raise GateError("staged-view contract is not fail-closed PASS")
    by_id = {v.get("view_id", v.get("type", "")): v for v in views.get("views", [])}
    required_view_fragments = ("combined-bgzf-all10-v1", "graph-fixture-v1")
    text = json.dumps(views, sort_keys=True)
    for fragment in required_view_fragments:
        if fragment not in text:
            raise GateError(f"missing staged view contract {fragment}")
    for view in views.get("views", []):
        if view.get("reversible") is not True and view.get("reversible_name_mapping") is not True and view.get("reversible_name_map") != "identity":
            raise GateError(f"view is not declared reversible: {view}")
        if not view.get("quota_bytes") or "cleanup" not in view:
            raise GateError("view lacks explicit quota/cleanup")
    scratch = Path(views["scratch_namespace"])
    if scratch.exists():
        raise GateError(f"external sequence-bearing scratch view was not cleaned: {scratch}")
    groups["staged_views_reversible_quota_cleaned"] = PASS

    resource = json.loads((release / "resource_summary.json").read_text())
    if resource.get("verdict") != PASS or not resource.get("ram_le_70_percent"):
        raise GateError("resource summary is not PASS")
    if resource.get("peak_rss_fraction", 2) > .70 or resource.get("swap_growth_bytes") != 0:
        raise GateError("RAM/swap threshold failed")
    if resource.get("scratch_upper95_fraction", 2) > .70 or resource.get("durable_upper95_fraction", 2) > .70:
        raise GateError("disk upper-95 threshold failed")
    if resource.get("projected_files_fraction", 2) > .50 or resource.get("actual_peak_files_fraction", 2) > .50:
        raise GateError("inode threshold failed")
    if min(resource.get("unfinished_write_reserve_factor_durable", 0), resource.get("unfinished_write_reserve_factor_scratch", 0)) < 2:
        raise GateError("unfinished-write reserve failed")
    groups["resource_thresholds"] = PASS

    restart = json.loads((release / "restart_evidence.json").read_text())
    if restart.get("verdict") != PASS or not restart.get("partial_never_published") or not restart.get("partial_cleaned"):
        raise GateError("kill/restart fail-closed evidence mismatch")
    if not restart.get("clean_restart_promoted_only_after_complete"):
        raise GateError("clean restart promotion evidence mismatch")
    groups["injected_kill_restart"] = PASS

    deterministic = json.loads((release / "deterministic_rerun.json").read_text())
    required_determinism = {
        "combined_bgzf_byte_identical", "impg_six_file_byte_identical",
        "impg_query_byte_identical", "impg_map_byte_identical", "mash_triangle_byte_identical",
    }
    if deterministic.get("verdict") != PASS or any(deterministic.get(key) != PASS for key in required_determinism):
        raise GateError("deterministic rerun evidence is incomplete/non-PASS")
    groups["deterministic_rerun"] = PASS

    roots = json.loads((release / "root_hashes.json").read_text())
    current = verify_root_hashes(repo)
    if roots.get("start") != current or roots.get("finish") != current or roots.get("verdict") != PASS:
        raise GateError("root source immutability start/finish mismatch")
    groups["root_inputs_start_finish"] = PASS

    if (release / "failures.jsonl").read_text().strip():
        raise GateError("immutable successful release contains failure ledger entries")
    states = [json.loads(x) for x in (release / "state.jsonl").read_text().splitlines() if x]
    if not states or states[0].get("status") != "START" or states[-1].get("status") != "VALIDATED":
        raise GateError("append-only state ledger does not end VALIDATED")
    commands = [json.loads(x) for x in (release / "commands.jsonl").read_text().splitlines() if x]
    if not commands or any("exit_code" not in x or "argv" not in x for x in commands):
        raise GateError("command/argv provenance ledger incomplete")
    groups["state_failure_command_provenance"] = PASS

    # Git/task ownership: no sequence or index-bearing files and no >10 MiB task output.
    owned = [repo / "workflow/compatibility", repo / "artifacts/consumer_compatibility",
             repo / "manifests/consumer-compatibility-v1"]
    prohibited_suffixes = (".fa", ".fasta", ".fa.gz", ".fna", ".faa", ".gfa", ".og", ".vg", ".msh", ".impg", ".1gbwt")
    for root in owned:
        if not root.exists(): continue
        for path in root.rglob("*"):
            if path.is_file():
                if path.name.endswith(prohibited_suffixes):
                    raise GateError(f"sequence/index-bearing payload committed under task path: {path}")
                if path.stat().st_size > 10 * 1024 * 1024:
                    raise GateError(f"task-owned git artifact exceeds 10 MiB: {path}")
    groups["compact_git_ownership"] = PASS

    result = {
        "schema_version": "consumer-compatibility-validation-v1", "verdict": PASS,
        "release_id": release.name, "release_sha256": sha256_file(release / "release.json"),
        "sha256_inventory_rows": inventory_rows, "consumer_gate_count": len(gates),
        "groups": groups, "validated_at_utc": __import__("time").strftime("%Y-%m-%dT%H:%M:%SZ", __import__("time").gmtime()),
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--external-release", type=Path, required=True)
    p.add_argument("--output", type=Path)
    args = p.parse_args(argv)
    result = validate_release(args.repo.resolve(), args.external_release)
    data = json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(data)
    print(data, end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except GateError as exc:
        print(f"NO_GO: {exc}", file=__import__("sys").stderr)
        raise SystemExit(2)
