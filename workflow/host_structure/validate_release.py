#!/usr/bin/env python3
"""Independent checksum and semantic validation of host-structure-1000 v1."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Sequence

from . import host_structure as h
from . import runner


def check(value: bool, message: str) -> None:
    if not value:
        raise h.GateError(message)


def plain_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate(repo: Path, external: Path) -> dict[str, Any]:
    inputs = h.verify_inputs(repo)
    inventory = h.verify_external_inventory(external)
    release = json.loads((external / "release.json").read_text())
    expected_id = runner.release_id()
    check(release.get("release_id") == expected_id, "release identity mismatch")
    check(release.get("schema_version") == h.SCHEMA, "release schema mismatch")
    check(release.get("verdict") == h.PASS and release.get("immutable") is True, "release is not immutable PASS")
    h.assert_host_only_manifest(release)
    gates = release.get("applicable_gates", {})
    h.require_pass_or_explicit_na(gates, "host structure release")
    check(gates.get("prophage_source_coordinate_policy") == h.NA_HOST, "prophage coordinate gate applicability is wrong")
    check(gates.get("prophage_extraction_semantics") == h.NA_HOST, "prophage extraction gate applicability is wrong")
    check(gates.get("integrated_rung_scale_criteria") == h.NA_RUNG, "scale applicability is wrong")
    check(release.get("root_inputs_start") == h.ROOT_HASHES == release.get("root_inputs_finish"), "root hash start/finish mismatch")
    check(release.get("host_only_input_allowlist") == list(h.HOST_ONLY_INPUT_ALLOWLIST), "host-only allow-list mismatch")

    input_manifest = h.read_hashed_tsv(external / "input_manifest.tsv", runner.INPUT_MANIFEST_FIELDS)
    check(len(input_manifest) == h.EXPECTED_N, "input manifest is not exact N=1,000")
    check([x["accession"] for x in input_manifest] == inputs["accessions"], "input manifest order/accessions mismatch")
    check(len({x["canonical_fasta_content_sha256"] for x in input_manifest}) <= h.EXPECTED_N, "input content union exceeds cap")
    for row, assembly, ref in zip(input_manifest, inputs["assemblies"], inputs["refs"]):
        check(row["canonical_bgzf_sha256"] == assembly["canonical_bgzf_sha256"], "input BGZF checksum mismatch")
        check(row["canonical_bgzf_path"] == str(h.canonical_object_path(ref, row["accession"])), "input object path mismatch")

    outputs = external / "outputs"
    qc = h.read_hashed_tsv(outputs / "host_qc.tsv", runner.QC_FIELDS)
    membership = h.read_hashed_tsv(outputs / "host_membership.tsv", runner.MEMBERSHIP_FIELDS)
    medoids = h.read_hashed_tsv(outputs / "medoids_and_cases.tsv", runner.MEDOID_FIELDS)
    clades = h.read_hashed_tsv(outputs / "host_clades.tsv", runner.CLADE_FIELDS)
    alternatives = h.read_hashed_tsv(outputs / "alternative_partitions.tsv", runner.ALT_FIELDS)
    check(len(qc) == len(membership) == h.EXPECTED_N, "QC/membership all-host accounting failed")
    check([r["tip_id"] for r in qc] == inputs["accessions"], "QC tip order mismatch")
    check([r["tip_id"] for r in membership] == inputs["accessions"], "membership tip order mismatch")
    check(all(r["eligible"] == h.PASS and r["pansn_name_roundtrip"] == h.PASS for r in qc), "QC eligibility/name gate failed")
    check(all(r["representative_tip"] in set(inputs["accessions"]) for r in membership), "unknown representative tip")
    check(all(r["placement_status"] in {"SUPPORTED_FIXED", "AMBIGUOUS_UNROOTED_BACKBONE_OR_UNSUPPORTED"} for r in membership), "forced/unknown placement status")
    check(len(alternatives) == 3 * h.EXPECTED_N, "alternative partition accounting failed")
    # Three complete N=1,000 alternatives, regardless of each partition's k.
    by_partition: dict[str, list[dict[str, str]]] = {}
    for row in alternatives:
        by_partition.setdefault(row["partition_id"], []).append(row)
    check(set(by_partition) == {"host_genetic_medoid_k12", "host_genetic_medoid_k16", "host_genetic_medoid_k20"}, "alternative partition IDs mismatch")
    check(all(len(rows) == h.EXPECTED_N for rows in by_partition.values()), "alternative partition lost hosts")
    check(all(row["host_only"] == h.PASS for row in alternatives), "alternative partition is not host-only")

    clade_ids = {r["clade_id"] for r in clades}
    frozen_ids = {r["frozen_clade_id"] for r in membership if r["frozen_clade_id"]}
    check(frozen_ids <= clade_ids, "membership references unknown clade")
    for clade in clades:
        check(clade["membership_status"] == "SUPPORTED_FIXED", "unsupported clade was frozen")
        check(clade["tree_release_state"] == "FROZEN_BEFORE_PHAGE_ASSOCIATION", "clade publication order wrong")
        check(float(clade["mash_parameter_support"]) >= 1.0, "clade lacks parameter stability")
        check(float(clade["mash_resampling_support"]) >= .95, "clade lacks resampling stability")
        count = sum(r["frozen_clade_id"] == clade["clade_id"] for r in membership)
        check(count == int(clade["host_count"]), "clade membership count mismatch")

    mash_metrics = json.loads((outputs / "mash_metrics.json").read_text())
    check(mash_metrics["exact_pair_validation"]["unordered_off_diagonal_pairs"] == h.EXPECTED_PAIRS, "Mash pair count mismatch")
    check(mash_metrics["exact_pair_validation"]["directed_records"] == h.EXPECTED_N ** 2, "Mash directed count mismatch")
    check(mash_metrics["exact_pair_validation"]["symmetry"] == h.PASS, "Mash symmetry failed")
    check(len(mash_metrics["sensitivity"]) == 4 and mash_metrics["resampling_replicates"] == 6, "Mash sensitivity/resampling run count mismatch")
    baseline_root = outputs / "mash/k21_s10000_seed42"
    labels, matrix = h.parse_mash_triangle(baseline_root / "triangle.phylip", inputs["accessions"])
    exact = h.validate_directed_mash(baseline_root / "directed_all_pairs.tsv", labels, matrix)
    check(exact == mash_metrics["exact_pair_validation"], "independent Mash exact validation differs")
    for cfg in mash_metrics["configs"]:
        root = outputs / "mash" / cfg["name"]
        labels2, _ = h.parse_mash_triangle(root / "triangle.phylip", inputs["accessions"])
        tree = h.parse_newick((root / "rapidnj.unrooted.nwk").read_text())
        check(set(h.leaf_names(tree)) == set(labels2) and len(h.leaf_names(tree)) == h.EXPECTED_N, "Mash tree exact tips failed")
    overview = h.parse_newick((outputs / "trees/all_host_mash_supported.unrooted.nwk").read_text())
    check(set(h.leaf_names(overview)) == set(inputs["accessions"]), "supported overview lost tips")

    high = json.loads((outputs / "high_fidelity_metrics.json").read_text())
    check(high.get("verdict") == h.PASS and len(high.get("lineages", [])) == 16, "high-fidelity ensemble release is not PASS/16 accounted")
    check(high.get("pass_lineages", 0) >= 1 and high.get("failed_lineages_used_for_clade_inference") == 0,
          "high-fidelity failures were bypassed or no lineage passed")
    check(high.get("pass_lineages", 0) + high.get("blocked_ambiguous_lineages", 0) == 16,
          "high-fidelity PASS/ambiguous lineage accounting mismatch")
    check(high.get("rooting") == runner.PARAMETERS["root_policy"], "root policy mismatch")
    check(high.get("outgroup_acquisition") == "BLOCKED_BY_FROZEN_COHORT", "outgroup acquisition contract mismatch")
    selected_tips = {r["selected_tip"] for r in medoids}
    check(selected_tips <= set(inputs["accessions"]) and len(medoids) >= 64, "high-fidelity selection is missing/noncohort")
    lineage_verdict = {row["lineage_id"]: row["verdict"] for row in high["lineages"]}
    check(all(not row["frozen_clade_id"] or lineage_verdict[row["sampling_partition_k16"]] == h.PASS for row in membership),
          "membership from blocked high-fidelity lineage entered a frozen clade")
    for lineage in high["lineages"]:
        if lineage["verdict"] != h.PASS:
            check(lineage["verdict"].startswith("AMBIGUOUS_") and "BLOCKED_FROM_CLADE_INFERENCE" in lineage["verdict"],
                  "failed high-fidelity lineage was not explicitly ambiguous/blocked")
            continue
        check(lineage["reference_split_concordance_gate"] == h.PASS, "reference bias gate failed")
        for role in ("primary_reference", "alternative_reference"):
            root = outputs / "core" / lineage["lineage_id"] / role
            metric = json.loads((root / "metrics.json").read_text())
            check(all(x == h.PASS for x in metric["scientific_gates"].values()), "core scientific gate failed")
            check(metric["nonrecombinant_parsimony_informative_sites"] >= 100, "informative site gate failed")
            check(metric["recombination_candidate_mask"]["method"].startswith("host-core"), "recombination handling missing")
            support_rows = plain_tsv(root / "split_support.tsv")
            check(sum(float(row["bootstrap_support"]) >= .95 for row in support_rows) == metric["supported_splits_ge_95pct"],
                  "branch-level bootstrap support accounting mismatch")
            check(all(set(row["smaller_side_tips"].split(",")) <= set(metric["samples"]) for row in support_rows),
                  "branch support references unknown tip")
            tree = h.parse_newick((root / "core_snp.support_collapsed.unrooted.nwk").read_text())
            check(set(h.leaf_names(tree)) == set(metric["samples"]), "core tree tip mismatch")
            # Sequence-bearing alignments remain external but are structurally checked.
            for name in ("core_alignment.fa.gz", "core_alignment.recombination_candidate_masked.fa.gz"):
                count = 0
                with gzip.open(root / name, "rt") as handle:
                    for line in handle:
                        if line.startswith(">"):
                            count += 1
                check(count == len(metric["samples"]), "core alignment sample count mismatch")

    resources = json.loads((external / "resource_summary.json").read_text())
    check(resources.get("verdict") == h.PASS and all(resources.get("checks", {}).values()), "resource summary not PASS")
    check(resources["peak_rss_fraction"] <= .70, "peak RSS threshold exceeded")
    state = [json.loads(line) for line in (external / "state.jsonl").read_text().splitlines()]
    check(any(x.get("event") == "INJECTED_MATERIALIZATION_SIGKILL" for x in state), "injected kill evidence missing")
    check(any(x.get("event") == "READY_TO_PROMOTE" for x in state), "ready-to-promote event missing")
    check(not any("DOWNLOAD" in x.get("event", "") or "ACQUISITION" in x.get("event", "") for x in state), "genome acquisition event present")

    compact_semantics = {
        "release_id": release["release_id"], "input_manifest_sha256": h.sha256_file(external / "input_manifest.tsv"),
        "membership_sha256": h.sha256_file(outputs / "host_membership.tsv"),
        "clades_sha256": h.sha256_file(outputs / "host_clades.tsv"),
        "alternative_partitions_sha256": h.sha256_file(outputs / "alternative_partitions.tsv"),
        "overview_tree_sha256": h.sha256_file(outputs / "trees/all_host_mash_supported.unrooted.nwk"),
        "mash_exact": exact, "high_fidelity_lineages": len(high["lineages"]),
        "all_host_rows": len(membership), "root_hashes": h.verify_root_hashes(repo),
    }
    semantic_sha = hashlib.sha256(h.canonical_json(compact_semantics)).hexdigest()
    return {
        "schema": "host-structure-1000-independent-validation-v1", "verdict": h.PASS,
        "release_id": release["release_id"], "inventory_files": len(inventory),
        "assembly_rows": len(input_manifest), "qc_rows": len(qc), "membership_rows": len(membership),
        "clades": len(clades), "ambiguous_memberships": sum(not r["frozen_clade_id"] for r in membership),
        "alternative_partition_rows": len(alternatives), "mash_pairs": exact["unordered_off_diagonal_pairs"],
        "mash_directed_records": exact["directed_records"], "high_fidelity_lineages": len(high["lineages"]),
        "semantic_sha256": semantic_sha, "semantic": compact_semantics,
    }


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--repo", type=Path, default=Path("."))
    p.add_argument("--external", type=Path, required=True)
    p.add_argument("--output", type=Path)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = validate(args.repo.resolve(), args.external.resolve())
    except h.GateError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    data = h.canonical_json(result)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True); args.output.write_bytes(data)
    sys.stdout.buffer.write(data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
