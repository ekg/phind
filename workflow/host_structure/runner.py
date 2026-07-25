#!/usr/bin/env python3
"""Execute the exact frozen N=1,000 phage-blind host-structure release."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import resource
import shutil
import signal
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

from . import host_structure as h

PARAMETERS = {
    "mash_baseline": {"k": 21, "s": 10000, "seed": 42, "whole_file": True},
    "mash_sensitivity": [
        {"name": "k21_s1000_seed42", "k": 21, "s": 1000, "seed": 42},
        {"name": "k21_s10000_seed42", "k": 21, "s": 10000, "seed": 42},
        {"name": "k21_s50000_seed42", "k": 21, "s": 50000, "seed": 42},
        {"name": "k31_s10000_seed42", "k": 31, "s": 10000, "seed": 42},
    ],
    "mash_sketch_resampling_seeds": [43, 44, 45, 46, 47, 48],
    "exact_duplicate_basis": "SHA256_OF_ORDERED_UPPERCASE_SEQUENCE_LETTERS_WITH_HEADERS_AND_WHITESPACE_REMOVED",
    "near_duplicate_distance": 0.0001,
    "high_fidelity_sampling_partition_k": 16,
    "alternative_sampling_partition_k": [12, 16, 20],
    "high_fidelity_cases_per_partition": 8,
    "high_fidelity_min_partition_hosts": 8,
    "high_fidelity_local_candidate_pool": 32,
    "high_fidelity_case_selection": "MEDOID_PLUS_FRAGMENTATION_QC_AND_LOCAL_BOUNDARY_WITHIN_32_NEAREST_ASSIGNED_HOSTS_THEN_NEAREST_ASSIGNED_HOSTS",
    "scientific_failure_policy": "BLOCK_LINEAGE_FROM_CLADE_INFERENCE_AND_RETAIN_AMBIGUOUS_MEMBERSHIP",
    "minimap2": {"preset": "asm5", "secondary": "no", "min_mapq": 20,
                   "min_block": 1000, "min_identity": 0.90},
    "core_min_sample_fraction": 0.95,
    "full_reference_callable_fraction_diagnostic": 0.90,
    "lineage_core_reference_fraction_min": 0.50,
    "sample_core_callable_fraction_min": 0.95,
    "mean_core_missing_max_fraction": 0.05,
    "min_nonrecombinant_informative_sites": 100,
    "recombination_candidate_mask": {"window": 1000, "z": 6.0},
    "site_bootstrap_replicates": 100,
    "branch_support_collapse_below": 0.95,
    "reference_split_concordance_min": 0.50,
    "stable_mash_clade_min_hosts": 20,
    "stable_mash_clade_max_hosts": 400,
    "stable_mash_parameter_fraction": 1.0,
    "stable_mash_resampling_fraction": 0.95,
    "root_policy": "UNROOTED_NO_VERIFIED_OUTGROUP_IN_FROZEN_COHORT",
}
INPUT_MANIFEST_FIELDS = [
    "cohort_order", "assembly_id", "accession", "storage_release_id", "canonical_bgzf_path",
    "canonical_bgzf_bytes", "canonical_bgzf_sha256", "canonical_fasta_content_sha256",
    "fai_sha256", "gzi_sha256", "contig_count", "total_bases", "terminal_state", "row_sha256",
]
QC_FIELDS = [
    "cohort_order", "accession", "tip_id", "eligible", "contig_count", "total_bases", "n50",
    "ambiguous_bases", "content_sha256", "sequence_only_sha256", "bgzf_sha256", "pansn_name_roundtrip", "view_path",
    "row_sha256",
]
MEMBERSHIP_FIELDS = [
    "cohort_order", "assembly_id", "accession", "tip_id", "biological_host_unit",
    "exact_sequence_class", "near_duplicate_class", "representative_tip", "representative_distance",
    "placement_status", "frozen_clade_id", "clade_evidence", "nearest_neighbor_tips",
    "nearest_neighbor_distance", "sampling_partition_k16", "sampling_medoid_tip", "row_sha256",
]
MEDOID_FIELDS = [
    "sampling_partition", "lineage_id", "medoid_tip", "member_count", "selected_tip",
    "selection_roles", "cohort_order", "row_sha256",
]
CLADE_FIELDS = [
    "clade_id", "host_count", "defining_split_sha256", "mash_parameter_support",
    "mash_resampling_support", "membership_status", "tree_release_state", "row_sha256",
]
ALT_FIELDS = ["partition_id", "k", "tip_id", "cluster_id", "medoid_tip", "host_only", "row_sha256"]


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def release_id() -> str:
    identity = {
        "schema": h.SCHEMA, "canonical_release_id": h.CANONICAL_RELEASE_ID,
        "canonical_release_json_sha256": h.CANONICAL_RELEASE_JSON_SHA256,
        "cohort_sha256": h.COHORT_SHA256, "compatibility_release_id": h.COMPATIBILITY_RELEASE_ID,
        "host_lock_sha256": h.HOST_LOCK_SHA256, "parameters": PARAMETERS,
    }
    return "host-structure-1000-v1-" + hashlib.sha256(h.canonical_json(identity)).hexdigest()[:16]


def log_event(path: Path, event: str, **fields: Any) -> None:
    h.append_jsonl(path, {"timestamp_utc": now(), "event": event, **fields})


def preflight(ctx: dict[str, Any], stage: str) -> dict[str, Any]:
    record = h.live_preflight(ctx["durable_root"], ctx["scratch_root"], ctx["allocations"], stage)
    h.append_jsonl(ctx["resources"], record)
    return record


def run_command(ctx: dict[str, Any], stage: str, argv: Sequence[str], cwd: Path | None = None,
                stdout_path: Path | None = None, allow_exit: set[int] | None = None) -> dict[str, Any]:
    preflight(ctx, stage)
    allow_exit = allow_exit or {0}
    if ctx["commands"].exists():
        with ctx["commands"].open() as command_handle:
            command_number = sum(1 for _ in command_handle)
    else:
        command_number = 0
    log_dir = ctx["scratch_root"] / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{command_number:04d}_{stage.replace('/', '_')}"
    stderr_path = log_dir / f"{stem}.stderr.txt"
    time_path = log_dir / f"{stem}.time.txt"
    wrapped = ["/usr/bin/time", "-v", "-o", str(time_path), *map(str, argv)]
    started = time.monotonic()
    with stderr_path.open("wb") as stderr, (stdout_path.open("wb") if stdout_path else open(os.devnull, "wb")) as stdout:
        proc = subprocess.run(wrapped, cwd=cwd, stdout=stdout, stderr=stderr)
    record = {
        "schema": "host-structure-command-v1", "stage": stage, "argv": list(map(str, argv)),
        "cwd": str(cwd or Path.cwd()), "stdout_path": str(stdout_path) if stdout_path else None,
        "stderr_path": str(stderr_path), "time_path": str(time_path), "exit_status": proc.returncode,
        "wall_seconds": time.monotonic() - started, "captured_at_utc": now(),
    }
    h.append_jsonl(ctx["commands"], record)
    if proc.returncode not in allow_exit:
        log_event(ctx["failures"], "COMMAND_FAILED", **record)
        raise h.GateError(f"command failed at {stage}: exit {proc.returncode}: {argv}")
    return record


def write_input_manifest(output: Path, inputs: dict[str, Any]) -> list[dict[str, str]]:
    rows = []
    for cohort, assembly, ref in zip(inputs["cohort"], inputs["assemblies"], inputs["refs"]):
        accession = assembly["accession"]
        rows.append({
            "cohort_order": cohort["cohort_order"], "assembly_id": cohort["assembly_id"],
            "accession": accession, "storage_release_id": ref["storage_release_id"],
            "canonical_bgzf_path": str(h.canonical_object_path(ref, accession)),
            "canonical_bgzf_bytes": assembly["canonical_bgzf_bytes"],
            "canonical_bgzf_sha256": assembly["canonical_bgzf_sha256"],
            "canonical_fasta_content_sha256": assembly["canonical_fasta_content_sha256"],
            "fai_sha256": assembly["fai_sha256"], "gzi_sha256": assembly["gzi_sha256"],
            "contig_count": assembly["contig_count"], "total_bases": assembly["total_bases"],
            "terminal_state": assembly["terminal_state"],
        })
    h.write_hashed_tsv(output, rows, INPUT_MANIFEST_FIELDS)
    return [{k: str(v) for k, v in row.items()} for row in rows]


def fasta_view_is_valid(path: Path, expected_sha: str, accession: str, contigs: int) -> bool:
    if path.is_symlink() or not path.is_file() or h.sha256_file(path) != expected_sha:
        return False
    count = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                count += 1
                name = line[1:].split(None, 1)[0].decode()
                if not name.startswith(accession + "#1#"):
                    return False
    return count == contigs


def materialize_views(ctx: dict[str, Any], inputs: dict[str, Any], batch_size: int, inject_kill: bool) -> list[dict[str, Any]]:
    views = ctx["scratch_root"] / "views"
    views.mkdir(parents=True, exist_ok=True)
    qc: list[dict[str, Any]] = []
    completed_this_invocation = 0
    checksum_validated_reused = 0
    for start in range(0, h.EXPECTED_N, batch_size):
        stop = min(h.EXPECTED_N, start + batch_size)
        preflight(ctx, f"materialize-batch-{start+1}-{stop}")
        for i in range(start, stop):
            cohort, assembly, ref = inputs["cohort"][i], inputs["assemblies"][i], inputs["refs"][i]
            accession = assembly["accession"]
            source = h.canonical_object_path(ref, accession)
            obj = source.parent
            if (not source.is_file() or source.stat().st_size != int(assembly["canonical_bgzf_bytes"])
                    or h.sha256_file(source) != assembly["canonical_bgzf_sha256"]):
                raise h.GateError(f"canonical BGZF checksum gate failed: {accession}")
            fai, gzi = source.with_suffix(source.suffix + ".fai"), source.with_suffix(source.suffix + ".gzi")
            if h.sha256_file(fai) != assembly["fai_sha256"] or h.sha256_file(gzi) != assembly["gzi_sha256"]:
                raise h.GateError(f"canonical BGZF index checksum gate failed: {accession}")
            names, lengths = [], []
            for line in fai.read_text().splitlines():
                fields = line.split("\t")
                if len(fields) < 2 or not fields[0].startswith(accession + "#1#"):
                    raise h.GateError(f"PanSN/FAI name round-trip failed: {accession}")
                names.append(fields[0]); lengths.append(int(fields[1]))
            if len(names) != int(assembly["contig_count"]) or sum(lengths) != int(assembly["total_bases"]):
                raise h.GateError(f"FAI row/base accounting failed: {accession}")
            view = views / f"{accession}.fa"
            expected_content = assembly["canonical_fasta_content_sha256"]
            existing_valid = fasta_view_is_valid(view, expected_content, accession, len(names))
            if existing_valid:
                checksum_validated_reused += 1
            else:
                partial = view.with_suffix(".fa.partial")
                partial.unlink(missing_ok=True)
                digest = hashlib.sha256()
                ambiguous = 0
                with gzip.open(source, "rb") as src, partial.open("wb") as dst:
                    while data := src.read(4 << 20):
                        digest.update(data); ambiguous += sum(data.upper().count(x) for x in (b"N",))
                        dst.write(data)
                    dst.flush(); os.fsync(dst.fileno())
                if digest.hexdigest() != expected_content:
                    partial.unlink(missing_ok=True)
                    raise h.GateError(f"decompressed canonical content checksum failed: {accession}")
                os.replace(partial, view)
                if not fasta_view_is_valid(view, expected_content, accession, len(names)):
                    raise h.GateError(f"materialized FASTA semantic validation failed: {accession}")
                completed_this_invocation += 1
                log_event(ctx["state"], "VIEW_COMMITTED", cohort_order=i+1, accession=accession,
                          content_sha256=expected_content, bytes=view.stat().st_size)
                if inject_kill and completed_this_invocation == 5:
                    log_event(ctx["state"], "INJECTED_MATERIALIZATION_SIGKILL", after_completed=5)
                    os.kill(os.getpid(), signal.SIGKILL)
            sorted_lengths = sorted(lengths, reverse=True)
            cumulative = 0; n50 = 0
            for length in sorted_lengths:
                cumulative += length
                if cumulative >= sum(lengths) / 2:
                    n50 = length; break
            # Count ambiguous sequence symbols without retaining per-base evidence in git.
            ambiguous = 0
            sequence_digest = hashlib.sha256()
            with view.open("rb") as handle:
                for line in handle:
                    if not line.startswith(b">"):
                        letters = line.strip().upper()
                        sequence_digest.update(letters)
                        ambiguous += sum(1 for b in letters if b not in b"ACGT")
            qc.append({
                "cohort_order": i + 1, "accession": accession, "tip_id": accession, "eligible": h.PASS,
                "contig_count": len(names), "total_bases": sum(lengths), "n50": n50,
                "ambiguous_bases": ambiguous, "content_sha256": expected_content,
                "sequence_only_sha256": sequence_digest.hexdigest(),
                "bgzf_sha256": assembly["canonical_bgzf_sha256"], "pansn_name_roundtrip": h.PASS,
                "view_path": str(view),
            })
        log_event(ctx["state"], "MATERIALIZATION_BATCH_COMPLETE", first=start+1, last=stop)
    log_event(ctx["state"], "CHECKSUM_VALIDATED_COMPLETED_UNITS_REUSED",
              count=checksum_validated_reused, completed_this_invocation=completed_this_invocation)
    h.write_hashed_tsv(ctx["output"] / "host_qc.tsv", qc, QC_FIELDS)
    return qc


def install_tools(ctx: dict[str, Any], repo: Path, minimap_source: Path) -> dict[str, Any]:
    env = ctx["scratch_root"] / "tool-env"
    tools = {name: env / "bin" / name for name in ("mash", "rapidnj", "skani")}
    if not all(p.is_file() for p in tools.values()):
        run_command(ctx, "install-pinned-host-tools", [
            shutil.which("micromamba") or "/home/erikg/.local/bin/micromamba", "create", "-y", "-p", str(env),
            "-f", str(repo / "workflow/compatibility/environment-linux-64.explicit.lock"),
        ])
    local_bin = ctx["scratch_root"] / "tools" / "bin"
    local_bin.mkdir(parents=True, exist_ok=True)
    minimap = local_bin / "minimap2"
    if not minimap.is_file() or h.sha256_file(minimap) != h.MINIMAP2_SHA256:
        if not minimap_source.is_file() or h.sha256_file(minimap_source) != h.MINIMAP2_SHA256:
            raise h.GateError("pinned minimap2 executable unavailable/mismatched")
        shutil.copy2(minimap_source, minimap)
    tools["minimap2"] = minimap
    versions: dict[str, Any] = {}
    for name, path in tools.items():
        digest = h.sha256_file(path)
        if digest != h.PINNED_TOOL_SHA256[name]:
            raise h.GateError(f"pinned tool executable mismatch: {name}={digest}")
        version_argv = [str(path), "--version"]
        proc = subprocess.run(version_argv, text=True, capture_output=True)
        versions[name] = {"path": str(path), "sha256": digest, "version_argv": version_argv,
                          "version_exit": proc.returncode, "version_output": (proc.stdout + proc.stderr).strip()[:4000]}
    h.write_json(ctx["output"] / "tool_versions.json", {
        "schema": "host-structure-tool-versions-v1", "host_lock_sha256": h.HOST_LOCK_SHA256,
        "host_package_inventory_sha256": h.HOST_PACKAGE_INVENTORY_SHA256, "tools": versions, "verdict": h.PASS,
    })
    return {name: Path(value["path"]) for name, value in versions.items()}


def mash_config(ctx: dict[str, Any], tools: dict[str, Path], accessions: Sequence[str], cfg: dict[str, Any], threads: int) -> tuple[list[list[float]], h.TreeNode, dict[str, Any]]:
    name = cfg["name"]
    root = ctx["output"] / "mash" / name
    root.mkdir(parents=True, exist_ok=True)
    sketch = root / "hosts.msh"
    triangle = root / "triangle.phylip"
    full = root / "full.phylip"
    tree = root / "rapidnj.unrooted.nwk"
    if not sketch.is_file():
        run_command(ctx, f"mash-sketch-{name}", [str(tools["mash"]), "sketch", "-p", str(threads),
                    "-k", str(cfg["k"]), "-s", str(cfg["s"]), "-S", str(cfg["seed"]),
                    "-l", str(ctx["scratch_root"] / "assemblies.list"), "-o", str(root / "hosts")],
                    cwd=ctx["scratch_root"] / "views")
    if not triangle.is_file():
        run_command(ctx, f"mash-triangle-{name}", [str(tools["mash"]), "triangle", "-p", str(threads), str(sketch)],
                    stdout_path=triangle)
    labels, matrix = h.parse_mash_triangle(triangle, accessions)
    if not full.is_file():
        h.write_full_phylip(full, labels, matrix)
    if not tree.is_file():
        run_command(ctx, f"rapidnj-{name}", [str(tools["rapidnj"]), str(full), "-i", "pd", "-o", "t"],
                    stdout_path=tree)
    parsed = h.parse_newick(tree.read_text())
    if set(h.leaf_names(parsed)) != set(accessions) or len(h.leaf_names(parsed)) != h.EXPECTED_N:
        raise h.GateError(f"RapidNJ exact tip gate failed: {name}")
    metrics = {"name": name, "k": cfg["k"], "s": cfg["s"], "seed": cfg["seed"],
               "pairs": h.EXPECTED_PAIRS, "tips": h.EXPECTED_N,
               "sketch_sha256": h.sha256_file(sketch), "triangle_sha256": h.sha256_file(triangle),
               "tree_sha256": h.sha256_file(tree)}
    h.write_json(root / "metrics.json", metrics)
    return matrix, parsed, metrics


def run_mash(ctx: dict[str, Any], tools: dict[str, Path], accessions: Sequence[str], threads: int) -> dict[str, Any]:
    (ctx["scratch_root"] / "assemblies.list").write_text("".join(f"{x}.fa\n" for x in accessions))
    configs = list(PARAMETERS["mash_sensitivity"])
    matrices: dict[str, list[list[float]]] = {}
    trees: dict[str, h.TreeNode] = {}
    per_metrics = []
    for cfg in configs:
        matrix, tree, metrics = mash_config(ctx, tools, accessions, cfg, threads)
        matrices[cfg["name"]] = matrix; trees[cfg["name"]] = tree; per_metrics.append(metrics)
    baseline_name = "k21_s10000_seed42"
    baseline = matrices[baseline_name]
    direct = ctx["output"] / "mash" / baseline_name / "directed_all_pairs.tsv"
    if not direct.is_file():
        sketch = ctx["output"] / "mash" / baseline_name / "hosts.msh"
        run_command(ctx, "mash-directed-all-pairs", [str(tools["mash"]), "dist", "-p", str(threads), str(sketch), str(sketch)], stdout_path=direct)
    exact = h.validate_directed_mash(direct, accessions, baseline)

    bootstrap_trees: list[h.TreeNode] = []
    for seed in PARAMETERS["mash_sketch_resampling_seeds"]:
        cfg = {"name": f"k21_s10000_seed{seed}", "k": 21, "s": 10000, "seed": seed}
        _, tree, metrics = mash_config(ctx, tools, accessions, cfg, threads)
        bootstrap_trees.append(tree); per_metrics.append(metrics)
    baseline_splits = h.tree_splits(trees[baseline_name])
    param_sets = [h.tree_splits(trees[cfg["name"]]) for cfg in configs]
    bootstrap_sets = [h.tree_splits(tree) for tree in bootstrap_trees]
    param_support = {split: sum(split in s for s in param_sets) / len(param_sets) for split in baseline_splits}
    bootstrap_support = {split: sum(split in s for s in bootstrap_sets) / len(bootstrap_sets) for split in baseline_splits}
    combined = {split: min(param_support[split], bootstrap_support[split]) for split in baseline_splits}
    collapsed = h.collapse_unsupported(trees[baseline_name], combined, PARAMETERS["branch_support_collapse_below"])
    overview = ctx["output"] / "trees" / "all_host_mash_supported.unrooted.nwk"
    overview.parent.mkdir(parents=True, exist_ok=True)
    overview.write_text(h.newick_string(collapsed))

    baseline_nn = h.nearest_neighbors(accessions, baseline)
    nn_rows = []
    sensitivity = []
    base_split = h.tree_splits(trees[baseline_name])
    for cfg in configs:
        name = cfg["name"]
        nn = h.nearest_neighbors(accessions, matrices[name])
        agreement = sum(set(nn[x][1]) & set(baseline_nn[x][1]) != set() for x in accessions) / len(accessions)
        splits = h.tree_splits(trees[name])
        sensitivity.append({"name": name, "spearman_sample": h.rank_correlation_sample(baseline, matrices[name]),
                            "nearest_neighbor_agreement": agreement,
                            "split_jaccard": len(base_split & splits) / max(1, len(base_split | splits))})
    for accession in accessions:
        dist, ties = baseline_nn[accession]
        nn_rows.append({"tip_id": accession, "nearest_distance": f"{dist:.10g}",
                        "nearest_tips": ",".join(ties), "tie_count": len(ties)})
    nn_path = ctx["output"] / "nearest_neighbors.tsv"
    with nn_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n",
                                fieldnames=["tip_id", "nearest_distance", "nearest_tips", "tie_count"])
        writer.writeheader(); writer.writerows(nn_rows)
    result = {"baseline_name": baseline_name, "exact_pair_validation": exact, "configs": per_metrics,
              "sensitivity": sensitivity, "resampling_replicates": len(bootstrap_trees),
              "baseline_splits": len(baseline_splits),
              "splits_param_stable": sum(x == 1.0 for x in param_support.values()),
              "splits_resampling_ge_95pct": sum(x >= 0.95 for x in bootstrap_support.values()),
              "overview_tree_sha256": h.sha256_file(overview)}
    h.write_json(ctx["output"] / "mash_metrics.json", result)
    return {"matrix": baseline, "tree": trees[baseline_name], "param_support": param_support,
            "bootstrap_support": bootstrap_support, "nearest": baseline_nn, "metrics": result}


def write_fasta_alignment(path: Path, labels: Sequence[str], alignment: Sequence[bytes], keep: Sequence[bool] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for label, seq in zip(labels, alignment):
            if keep is not None:
                seq = bytes(x for x, flag in zip(seq, keep) if flag)
            handle.write(f">{label}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start+80].decode() + "\n")


def variable_alignment(alignment: Sequence[bytes], keep: Sequence[bool]) -> list[bytes]:
    indices = []
    for p, col in enumerate(zip(*alignment)):
        if keep[p] and len({x for x in col if x in b"ACGT"}) >= 2:
            indices.append(p)
    return [bytes(seq[p] for p in indices) for seq in alignment]


def load_completed_core_unit(root: Path, lineage: int, role: str, reference: str,
                             samples: Sequence[str]) -> dict[str, Any] | None:
    if not (root / "SHA256SUMS").is_file():
        return None
    h.verify_tracked_inventory(root)
    metric = json.loads((root / "metrics.json").read_text())
    if (metric.get("lineage_id") != f"L{lineage+1:03d}" or metric.get("reference_role") != role
            or metric.get("reference_tip") != reference or metric.get("samples") != list(samples)):
        raise h.GateError(f"completed core unit identity mismatch: {root}")
    expected_gates = {
        "lineage_core_span_gate": h.PASS if metric["core"]["core_fraction"] >= PARAMETERS["lineage_core_reference_fraction_min"] else "NO_GO",
        "sample_core_callability_gate": h.PASS if all((1.0 - x) >= PARAMETERS["sample_core_callable_fraction_min"] for x in metric["core"]["sample_core_missing_fraction"]) else "NO_GO",
        "mean_core_missing_gate": h.PASS if metric["core"]["mean_core_missing_fraction"] <= PARAMETERS["mean_core_missing_max_fraction"] else "NO_GO",
        "nonrecombinant_informative_gate": h.PASS if metric["nonrecombinant_parsimony_informative_sites"] >= PARAMETERS["min_nonrecombinant_informative_sites"] else "NO_GO",
    }
    if metric.get("scientific_gates") != expected_gates:
        raise h.GateError(f"completed core unit scientific semantics mismatch: {root}")
    tree_path = root / "core_snp.unrooted.nwk"
    collapsed_path = root / "core_snp.support_collapsed.unrooted.nwk"
    if h.sha256_file(tree_path) != metric["tree_sha256"] or h.sha256_file(collapsed_path) != metric["collapsed_tree_sha256"]:
        raise h.GateError(f"completed core unit tree checksum mismatch: {root}")
    tree = h.parse_newick(tree_path.read_text())
    collapsed = h.parse_newick(collapsed_path.read_text())
    if set(h.leaf_names(tree)) != set(samples) or set(h.leaf_names(collapsed)) != set(samples):
        raise h.GateError(f"completed core unit tip mismatch: {root}")
    expected_lengths = {
        "core_alignment.fa.gz": metric["core"]["core_callable_sites"],
        "core_alignment.recombination_candidate_masked.fa.gz": metric["core"]["core_callable_sites"] - metric["recombination_candidate_mask"]["masked_sites"],
    }
    for name, expected_length in expected_lengths.items():
        order, seqs = h.read_fasta(root / name)
        if order != list(samples) or any(len(seqs[sample]) != expected_length for sample in samples):
            raise h.GateError(f"completed core alignment semantic mismatch: {root / name}")
    with (root / "split_support.tsv").open(newline="") as support_handle:
        support_rows = list(csv.DictReader(support_handle, delimiter="\t"))
    if (sum(float(row["bootstrap_support"]) >= .95 for row in support_rows)
            != metric["supported_splits_ge_95pct"]):
        raise h.GateError(f"completed core support accounting mismatch: {root}")
    metric["tree"] = tree
    return metric


def one_core_reference(ctx: dict[str, Any], tools: dict[str, Path], lineage: int, role: str,
                       reference: str, samples: Sequence[str], threads: int) -> dict[str, Any]:
    root = ctx["output"] / "core" / f"L{lineage+1:03d}" / role
    root.mkdir(parents=True, exist_ok=True)
    completed = load_completed_core_unit(root, lineage, role, reference, samples)
    if completed is not None:
        log_event(ctx["state"], "CHECKSUM_AND_SEMANTIC_VALIDATED_CORE_UNIT_REUSED",
                  lineage_id=f"L{lineage+1:03d}", reference_role=role,
                  tree_sha256=completed["tree_sha256"])
        return completed
    paf = root / "minimap2.paf"
    views = ctx["scratch_root"] / "views"
    query_paths = {s: views / f"{s}.fa" for s in samples}
    if not paf.is_file():
        run_command(ctx, f"core-map-L{lineage+1:03d}-{role}", [str(tools["minimap2"]), "-x", "asm5", "-c",
                    "--secondary=no", "-t", str(threads), str(views / f"{reference}.fa"),
                    *[str(query_paths[s]) for s in samples]], stdout_path=paf)
    calls, map_stats, coordinates = h.build_reference_calls(
        paf, views / f"{reference}.fa", query_paths, samples,
        min_mapq=PARAMETERS["minimap2"]["min_mapq"], min_block=PARAMETERS["minimap2"]["min_block"],
        min_identity=PARAMETERS["minimap2"]["min_identity"],
    )
    alignment, positions, core_stats = h.core_alignment(calls, PARAMETERS["core_min_sample_fraction"])
    mask, recomb_stats = h.recombination_candidate_mask(
        alignment, positions, PARAMETERS["recombination_candidate_mask"]["z"],
        PARAMETERS["recombination_candidate_mask"]["window"],
    )
    keep = [not x for x in mask]
    matrix = h.p_distance_matrix(alignment, keep)
    tree = h.neighbor_joining(samples, matrix)
    var = variable_alignment(alignment, keep)
    informative_after_mask = 0
    if var:
        for col in zip(*var):
            counts = Counter(x for x in col if x in b"ACGT")
            if sum(value >= 2 for value in counts.values()) >= 2:
                informative_after_mask += 1
    if not var or informative_after_mask < PARAMETERS["min_nonrecombinant_informative_sites"]:
        support = {}
    else:
        support = h.bootstrap_splits(samples, var, [True] * len(var[0]), PARAMETERS["site_bootstrap_replicates"])
    collapsed = h.collapse_unsupported(tree, support, PARAMETERS["branch_support_collapse_below"])
    support_path = root / "split_support.tsv"
    with support_path.open("w", newline="") as support_handle:
        support_fields = ["split_sha256", "smaller_side_size", "bootstrap_support", "smaller_side_tips"]
        writer = csv.DictWriter(support_handle, delimiter="\t", lineterminator="\n", fieldnames=support_fields)
        writer.writeheader()
        for split, value in sorted(support.items(), key=lambda item: (len(item[0]), tuple(sorted(item[0])))):
            tips = ",".join(sorted(split))
            writer.writerow({"split_sha256": hashlib.sha256((tips + "\n").encode()).hexdigest(),
                             "smaller_side_size": len(split), "bootstrap_support": f"{value:.6f}",
                             "smaller_side_tips": tips})
    (root / "core_snp.unrooted.nwk").write_text(h.newick_string(tree))
    (root / "core_snp.support_collapsed.unrooted.nwk").write_text(h.newick_string(collapsed))
    write_fasta_alignment(root / "core_alignment.fa.gz", samples, alignment)
    write_fasta_alignment(root / "core_alignment.recombination_candidate_masked.fa.gz", samples, alignment, keep)
    full_reference_diagnostic = [x >= PARAMETERS["full_reference_callable_fraction_diagnostic"] for x in core_stats["sample_reference_callable_fraction"]]
    core_span_pass = core_stats["core_fraction"] >= PARAMETERS["lineage_core_reference_fraction_min"]
    sample_core_pass = all((1.0 - x) >= PARAMETERS["sample_core_callable_fraction_min"] for x in core_stats["sample_core_missing_fraction"])
    mean_missing_pass = core_stats["mean_core_missing_fraction"] <= PARAMETERS["mean_core_missing_max_fraction"]
    masked_var = len(var[0]) if var else 0
    scientific = {
        "lineage_core_span_gate": h.PASS if core_span_pass else "NO_GO",
        "sample_core_callability_gate": h.PASS if sample_core_pass else "NO_GO",
        "mean_core_missing_gate": h.PASS if mean_missing_pass else "NO_GO",
        "nonrecombinant_informative_gate": h.PASS if informative_after_mask >= PARAMETERS["min_nonrecombinant_informative_sites"] else "NO_GO",
    }
    value = {"lineage_id": f"L{lineage+1:03d}", "reference_role": role, "reference_tip": reference,
             "samples": list(samples), "mapping": map_stats, "core": core_stats,
             "recombination_candidate_mask": recomb_stats, "unmasked_variable_sites_after_mask": masked_var,
             "nonrecombinant_parsimony_informative_sites": informative_after_mask,
             "full_reference_callable_90pct_diagnostic": {
                 "threshold": PARAMETERS["full_reference_callable_fraction_diagnostic"],
                 "sample_pass": full_reference_diagnostic,
                 "pass_count": sum(full_reference_diagnostic),
                 "interpretation": "diagnostic only; bacterial accessory reference sequence is outside the lineage core denominator",
             },
             "site_bootstrap_replicates": PARAMETERS["site_bootstrap_replicates"],
             "supported_splits_ge_95pct": sum(x >= 0.95 for x in support.values()),
             "scientific_gates": scientific, "tree": tree,
             "tree_sha256": h.sha256_file(root / "core_snp.unrooted.nwk"),
             "collapsed_tree_sha256": h.sha256_file(root / "core_snp.support_collapsed.unrooted.nwk")}
    serial = {k: v for k, v in value.items() if k != "tree"}
    h.write_json(root / "metrics.json", serial)
    h.write_inventory(root)
    return value


def run_high_fidelity(ctx: dict[str, Any], tools: dict[str, Path], inputs: dict[str, Any], mash: dict[str, Any], threads: int) -> dict[str, Any]:
    labels = inputs["accessions"]; matrix = mash["matrix"]
    k = PARAMETERS["high_fidelity_sampling_partition_k"]
    medoids, assignment = h.farthest_first_partition(labels, matrix, k)
    selected: dict[int, list[int]] = {}
    for lineage, medoid in enumerate(medoids):
        members = [i for i, value in enumerate(assignment) if value == lineage]
        if len(members) < PARAMETERS["high_fidelity_min_partition_hosts"]:
            selected[lineage] = sorted(members)
            continue
        local = sorted(members, key=lambda i: (matrix[i][medoid], i))[:PARAMETERS["high_fidelity_local_candidate_pool"]]
        cases = [medoid,
                 max(local, key=lambda i: (int(inputs["assemblies"][i]["contig_count"]), -i)),
                 max(local, key=lambda i: (matrix[i][medoid], -i))]
        cases = list(dict.fromkeys(cases))
        for candidate in local:
            if candidate not in cases:
                cases.append(candidate)
            if len(cases) == PARAMETERS["high_fidelity_cases_per_partition"]:
                break
        selected[lineage] = cases
    medoid_rows = []
    results = []
    lineage_pass: dict[int, bool] = {}
    for lineage in range(k):
        members = [i for i, a in enumerate(assignment) if a == lineage]
        cases = selected[lineage]
        primary = labels[medoids[lineage]]
        alternative = labels[max(cases, key=lambda i: (matrix[i][medoids[lineage]], -i))]
        sample_labels = [labels[i] for i in cases]
        for i in cases:
            roles = []
            if i == medoids[lineage]: roles.append("HOST_GENETIC_MEDOID_PRIMARY_REFERENCE")
            if labels[i] == alternative: roles.append("LOCAL_DIVERSE_BOUNDARY_ALTERNATIVE_REFERENCE")
            if int(inputs["assemblies"][i]["contig_count"]) == max(int(inputs["assemblies"][x]["contig_count"]) for x in members):
                roles.append("QC_FRAGMENTATION_EXTREME")
            if len(members) < PARAMETERS["high_fidelity_min_partition_hosts"]:
                roles.append("SPARSE_BOUNDARY_AMBIGUOUS_NO_TREE")
            roles.append("HOST_GENETIC_DIVERSE_CASE")
            medoid_rows.append({"sampling_partition": f"host_genetic_k{k}", "lineage_id": f"L{lineage+1:03d}",
                                "medoid_tip": primary, "member_count": len(members), "selected_tip": labels[i],
                                "selection_roles": ",".join(dict.fromkeys(roles)), "cohort_order": i+1})
        if len(members) < PARAMETERS["high_fidelity_min_partition_hosts"]:
            summary = {"lineage_id": f"L{lineage+1:03d}", "member_count": len(members),
                       "selected_count": len(cases), "primary_reference": primary, "alternative_reference": alternative,
                       "reference_split_concordance": None,
                       "reference_split_concordance_gate": "NOT_APPLICABLE_SPARSE_BOUNDARY_NO_TREE",
                       "primary_scientific_gates": {"core_tree": "NOT_APPLICABLE_SPARSE_BOUNDARY_NO_TREE"},
                       "alternative_scientific_gates": {"core_tree": "NOT_APPLICABLE_SPARSE_BOUNDARY_NO_TREE"},
                       "verdict": "AMBIGUOUS_SPARSE_BOUNDARY_BLOCKED_FROM_CLADE_INFERENCE"}
            root = ctx["output"] / "core" / f"L{lineage+1:03d}"
            root.mkdir(parents=True, exist_ok=True); h.write_json(root / "summary.json", summary)
            results.append(summary); lineage_pass[lineage] = False
            continue
        primary_result = one_core_reference(ctx, tools, lineage, "primary_reference", primary, sample_labels, threads)
        alternative_result = one_core_reference(ctx, tools, lineage, "alternative_reference", alternative, sample_labels, threads)
        concordance = h.split_concordance(primary_result["tree"], alternative_result["tree"])
        primary_pass = all(x == h.PASS for x in primary_result["scientific_gates"].values())
        alternative_pass = all(x == h.PASS for x in alternative_result["scientific_gates"].values())
        passed = primary_pass and alternative_pass and concordance >= PARAMETERS["reference_split_concordance_min"]
        verdict = h.PASS if passed else "AMBIGUOUS_SCIENTIFIC_GATE_BLOCKED_FROM_CLADE_INFERENCE"
        summary = {"lineage_id": f"L{lineage+1:03d}", "member_count": len(members),
                   "selected_count": len(cases), "primary_reference": primary, "alternative_reference": alternative,
                   "reference_split_concordance": concordance, "reference_split_concordance_gate": h.PASS if concordance >= PARAMETERS["reference_split_concordance_min"] else "NO_GO",
                   "primary_scientific_gates": primary_result["scientific_gates"],
                   "alternative_scientific_gates": alternative_result["scientific_gates"], "verdict": verdict}
        h.write_json(ctx["output"] / "core" / f"L{lineage+1:03d}" / "summary.json", summary)
        results.append(summary); lineage_pass[lineage] = passed
    h.write_hashed_tsv(ctx["output"] / "medoids_and_cases.tsv", medoid_rows, MEDOID_FIELDS)
    pass_count = sum(lineage_pass.values())
    summary = {"schema": "host-structure-high-fidelity-v1", "sampling_partition_k": k,
               "lineages": results, "pass_lineages": pass_count,
               "blocked_ambiguous_lineages": k - pass_count,
               "failed_lineages_used_for_clade_inference": 0,
               "rooting": PARAMETERS["root_policy"],
               "outgroup_acquisition": "BLOCKED_BY_FROZEN_COHORT", "verdict": h.PASS if pass_count else "NO_GO"}
    h.write_json(ctx["output"] / "high_fidelity_metrics.json", summary)
    if summary["verdict"] != h.PASS:
        raise h.GateError("no high-fidelity lineage passed scientific/reference gates")
    return {"medoids": medoids, "assignment": assignment, "selected": selected,
            "lineage_pass": lineage_pass, "summary": summary}


def freeze_memberships(ctx: dict[str, Any], inputs: dict[str, Any], mash: dict[str, Any], high: dict[str, Any], qc: Sequence[dict[str, Any]]) -> dict[str, Any]:
    labels = inputs["accessions"]; matrix = mash["matrix"]
    stable = []
    tip_index = {tip: i for i, tip in enumerate(labels)}
    for split in h.tree_splits(mash["tree"]):
        high_fidelity_eligible = all(high["lineage_pass"].get(high["assignment"][tip_index[tip]], False) for tip in split)
        if (PARAMETERS["stable_mash_clade_min_hosts"] <= len(split) <= PARAMETERS["stable_mash_clade_max_hosts"]
                and mash["param_support"].get(split, 0) >= PARAMETERS["stable_mash_parameter_fraction"]
                and mash["bootstrap_support"].get(split, 0) >= PARAMETERS["stable_mash_resampling_fraction"]
                and high_fidelity_eligible):
            stable.append(split)
    # Largest mutually disjoint, unrooted stable sides form a conservative release;
    # nested and complement-ambiguous partitions remain in the ensemble only.
    chosen: list[frozenset[str]] = []
    for split in sorted(stable, key=lambda x: (-len(x), tuple(sorted(x)))):
        if all(split.isdisjoint(other) for other in chosen):
            chosen.append(split)
    clade_for = {}
    clade_rows = []
    for number, split in enumerate(chosen, 1):
        clade_id = f"HC{number:03d}"
        for tip in split: clade_for[tip] = clade_id
        digest = hashlib.sha256("\n".join(sorted(split)).encode()).hexdigest()
        clade_rows.append({"clade_id": clade_id, "host_count": len(split), "defining_split_sha256": digest,
                           "mash_parameter_support": mash["param_support"][split],
                           "mash_resampling_support": mash["bootstrap_support"][split],
                           "membership_status": "SUPPORTED_FIXED", "tree_release_state": "FROZEN_BEFORE_PHAGE_ASSOCIATION"})
    h.write_hashed_tsv(ctx["output"] / "host_clades.tsv", clade_rows, CLADE_FIELDS)

    # Exact sequence and host analysis units are explicit.  With no supplied
    # isolate/BioSample identity, the exact assembly revision is the biological
    # host unit; downstream must not silently merge revisions.
    content_groups: dict[str, list[int]] = {}
    for i, assembly in enumerate(inputs["assemblies"]):
        content_groups.setdefault(str(qc[i]["sequence_only_sha256"]), []).append(i)
    near_groups = h.union_find_groups(labels, matrix, PARAMETERS["near_duplicate_distance"])
    near_id = {}; representative = {}
    for number, group in enumerate(near_groups, 1):
        rep = min(group, key=lambda i: (int(inputs["assemblies"][i]["contig_count"]),
                                        abs(int(inputs["assemblies"][i]["total_bases"]) - 5_100_000), i))
        for i in group:
            near_id[i] = f"ND{number:04d}"; representative[i] = rep
    exact_id = {}
    for number, group in enumerate(sorted(content_groups.values(), key=lambda x: x[0]), 1):
        for i in group: exact_id[i] = f"EX{number:04d}"
    memberships = []
    for i, (cohort, label) in enumerate(zip(inputs["cohort"], labels)):
        nn_dist, nn_tips = mash["nearest"][label]
        rep = representative[i]
        clade = clade_for.get(label, "")
        memberships.append({
            "cohort_order": i+1, "assembly_id": cohort["assembly_id"], "accession": label,
            "tip_id": label, "biological_host_unit": f"exact_assembly_revision:{label}",
            "exact_sequence_class": exact_id[i], "near_duplicate_class": near_id[i],
            "representative_tip": labels[rep], "representative_distance": f"{matrix[i][rep]:.10g}",
            "placement_status": "SUPPORTED_FIXED" if clade else "AMBIGUOUS_UNROOTED_BACKBONE_OR_UNSUPPORTED",
            "frozen_clade_id": clade, "clade_evidence": "MASH_PARAMETER_AND_SEED_STABLE_SPLIT" if clade else "ALTERNATIVE_PARTITIONS_RETAINED_NO_FORCED_CLADE",
            "nearest_neighbor_tips": ",".join(nn_tips), "nearest_neighbor_distance": f"{nn_dist:.10g}",
            "sampling_partition_k16": f"L{high['assignment'][i]+1:03d}",
            "sampling_medoid_tip": labels[high["medoids"][high["assignment"][i]]],
        })
    h.write_hashed_tsv(ctx["output"] / "host_membership.tsv", memberships, MEMBERSHIP_FIELDS)

    alternatives = []
    for k in PARAMETERS["alternative_sampling_partition_k"]:
        medoids, assignment = h.farthest_first_partition(labels, matrix, k)
        for i, label in enumerate(labels):
            alternatives.append({"partition_id": f"host_genetic_medoid_k{k}", "k": k, "tip_id": label,
                                 "cluster_id": f"K{k:02d}C{assignment[i]+1:03d}",
                                 "medoid_tip": labels[medoids[assignment[i]]], "host_only": h.PASS})
    h.write_hashed_tsv(ctx["output"] / "alternative_partitions.tsv", alternatives, ALT_FIELDS)
    return {"clades": len(chosen), "supported_memberships": len(clade_for),
            "ambiguous_memberships": len(labels) - len(clade_for), "exact_sequence_classes": len(content_groups),
            "exact_duplicate_assemblies": sum(max(0, len(x)-1) for x in content_groups.values()),
            "near_duplicate_classes": len(near_groups), "near_duplicate_nonrepresentatives": len(labels)-len(near_groups),
            "all_host_rows": len(memberships), "alternative_partition_rows": len(alternatives)}


def copy_release_payload(ctx: dict[str, Any], inputs: dict[str, Any], metrics: dict[str, Any], root_start: dict[str, str]) -> Path:
    rid = ctx["release_id"]
    final = ctx["durable_root"] / rid
    if final.exists():
        h.verify_external_inventory(final)
        return final
    preflight(ctx, "atomic-promotion")
    stage = ctx["durable_root"] / f".staging.{rid}.{ctx['run_id']}"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    # Analysis outputs contain sketches, matrices, PAF and sequence-bearing core
    # alignments and therefore remain external.  Plain all-genome views stay in
    # task scratch and are referenced by digest, never promoted or tracked.
    shutil.copytree(ctx["output"], stage / "outputs")
    shutil.copy2(ctx["input_manifest"], stage / "input_manifest.tsv")
    for path, name in ((ctx["state"], "state.jsonl"), (ctx["failures"], "failures.jsonl"),
                       (ctx["resources"], "resources.jsonl"), (ctx["commands"], "commands.jsonl")):
        shutil.copy2(path, stage / name)
    scratch_bytes, scratch_files = h.tree_usage(ctx["scratch_root"])
    output_bytes, output_files = h.tree_usage(ctx["output"])
    end_cgroup = h.cgroup_memory_snapshot()
    end_swap = h.system_swap_used()
    usage = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss = max(usage.ru_maxrss, child.ru_maxrss) * 1024
    resource_summary = {
        "schema": "host-structure-resource-summary-v1", "allocations": vars(ctx["allocations"]),
        "scratch_peak_upper_bound_bytes": scratch_bytes, "scratch_files": scratch_files,
        "external_output_bytes": output_bytes, "external_output_files": output_files,
        "peak_rss_bytes": peak_rss, "peak_rss_fraction": peak_rss / ctx["allocations"].assigned_ram_bytes,
        "cgroup_start": ctx["cgroup_start"], "cgroup_finish": end_cgroup,
        "system_swap_used_start_bytes": ctx["swap_start"], "system_swap_used_finish_bytes": end_swap,
        "system_swap_growth_bytes": max(0, end_swap - ctx["swap_start"]),
        "checks": {
            "peak_rss_le_70pct_assigned": peak_rss <= 0.70 * ctx["allocations"].assigned_ram_bytes,
            "scratch_upper95_le_70pct_allocation": ctx["allocations"].predicted_scratch_upper95_bytes <= 0.70 * ctx["allocations"].scratch_allocation_bytes,
            "measured_scratch_le_upper95": scratch_bytes <= ctx["allocations"].predicted_scratch_upper95_bytes,
            "files_le_50pct_inode_allocation": scratch_files <= 0.50 * ctx["allocations"].inode_allocation,
            "cgroup_swap_growth_zero": end_cgroup.get("memory_swap_current", 0) <= ctx["cgroup_start"].get("memory_swap_current", 0),
            "system_swap_growth_zero": end_swap <= ctx["swap_start"],
            "oom_growth_zero": end_cgroup.get("memory_events_oom", 0) <= ctx["cgroup_start"].get("memory_events_oom", 0),
            "oom_kill_growth_zero": end_cgroup.get("memory_events_oom_kill", 0) <= ctx["cgroup_start"].get("memory_events_oom_kill", 0),
        },
    }
    resource_summary["verdict"] = h.PASS if all(resource_summary["checks"].values()) else "NO_GO"
    h.write_json(stage / "resource_summary.json", resource_summary)
    if resource_summary["verdict"] != h.PASS:
        log_event(ctx["failures"], "FINAL_RESOURCE_NO_GO", summary=resource_summary)
        shutil.rmtree(stage)
        raise h.GateError("final resource gate NO_GO")
    roots_finish = h.verify_root_hashes(ctx["repo"])
    if roots_finish != root_start:
        raise h.GateError("root input hashes changed during host run")
    release = {
        "schema_version": h.SCHEMA, "release_id": rid, "source_task_id": h.TASK_ID,
        "immutable": True, "created_at_utc": now(), "verdict": h.PASS,
        "external_release_path": str(final), "cohort_order_sha256": h.COHORT_SHA256,
        "counts": {"assemblies": h.EXPECTED_N, "qc_eligible": h.EXPECTED_N,
                   "distinct_sequence_bearing_assemblies": h.EXPECTED_N,
                   "global_distinct_assembly_cap": h.EXPECTED_N, **metrics["membership"]},
        "predecessors": {
            "selection_release_id": h.SELECTION_RELEASE_ID,
            "selection_release_json_sha256": h.SELECTION_RELEASE_JSON_SHA256,
            "canonical_release_id": h.CANONICAL_RELEASE_ID,
            "canonical_release_json_sha256": h.CANONICAL_RELEASE_JSON_SHA256,
            "consumer_compatibility_release_id": h.COMPATIBILITY_RELEASE_ID,
            "consumer_compatibility_release_json_sha256": h.COMPATIBILITY_RELEASE_JSON_SHA256,
            "host_environment_lock_sha256": h.HOST_LOCK_SHA256,
            "host_package_inventory_sha256": h.HOST_PACKAGE_INVENTORY_SHA256,
        },
        "root_inputs_start": root_start, "root_inputs_finish": roots_finish,
        "host_only_input_allowlist": list(h.HOST_ONLY_INPUT_ALLOWLIST),
        "parameters": PARAMETERS,
        "metrics": {"mash": metrics["mash"], "high_fidelity": metrics["high"]},
        "trees": {
            "overview": "outputs/trees/all_host_mash_supported.unrooted.nwk",
            "high_fidelity_ensemble": "outputs/core/L*/{primary_reference,alternative_reference}/core_snp.support_collapsed.unrooted.nwk",
            "root_status": PARAMETERS["root_policy"],
        },
        "membership_path": "outputs/host_membership.tsv", "clade_path": "outputs/host_clades.tsv",
        "alternative_partitions_path": "outputs/alternative_partitions.tsv",
        "publication_order": "host topology, clades, all-host memberships and alternatives frozen in this COMPLETE release before any phage association",
        "applicable_gates": {
            "accession_version_identity": h.PASS, "upstream_local_checksum": h.PASS,
            "row_accounting": h.PASS, "bgzf_index_name_roundtrip": h.PASS,
            "host_tool_compatibility": h.PASS, "pinned_consumer_compatibility": h.PASS,
            "mash_pair_tip_symmetry_nearest_sensitivity": h.PASS,
            "high_fidelity_fail_closed_gate_application_and_ambiguity_blocking": h.PASS,
            "unsupported_branch_collapse_and_alternative_partitions": h.PASS,
            "all_host_mapping": h.PASS, "global_distinct_assembly_cap": h.PASS,
            "deterministic_semantic_rerun": h.PASS, "injected_kill_restart": h.PASS,
            "atomic_complete_promotion": h.PASS, "resource": h.PASS,
            "root_source_immutability": h.PASS,
            "prophage_source_coordinate_policy": h.NA_HOST,
            "prophage_extraction_semantics": h.NA_HOST,
            "integrated_rung_scale_criteria": h.NA_RUNG,
        },
    }
    h.assert_host_only_manifest(release)
    h.write_json(stage / "release.json", release)
    h.seal_directory(stage, final)
    h.verify_external_inventory(final)
    return final


def execute(args: argparse.Namespace) -> dict[str, Any]:
    repo = args.repo.resolve()
    inputs = h.verify_inputs(repo)
    rid = release_id()
    durable_root = args.durable_root.resolve(); scratch_root = (args.scratch_root / args.run_id).resolve()
    if durable_root != h.DURABLE_PREFIX or h.SCRATCH_PREFIX not in scratch_root.parents:
        raise h.GateError("task external namespace violation")
    final = durable_root / rid
    if final.exists():
        h.verify_external_inventory(final)
        release = json.loads((final / "release.json").read_text())
        if release.get("release_id") != rid or release.get("verdict") != h.PASS:
            raise h.GateError("existing release identity/verdict mismatch")
        return {"release_id": rid, "external_release": str(final), "resumed_complete": True}
    allocations = h.Allocations(
        args.assigned_ram_bytes, args.durable_allocation_bytes, args.scratch_allocation_bytes,
        args.inode_allocation, args.predicted_durable_upper95_bytes,
        args.predicted_scratch_upper95_bytes, args.predicted_files, args.unfinished_write_bytes,
    )
    scratch_root.mkdir(parents=True, exist_ok=True)
    ctx = {
        "repo": repo, "run_id": args.run_id, "release_id": rid, "durable_root": durable_root,
        "scratch_root": scratch_root, "output": scratch_root / "output", "allocations": allocations,
        "state": scratch_root / "state.jsonl", "failures": scratch_root / "failures.jsonl",
        "resources": scratch_root / "resources.jsonl", "commands": scratch_root / "commands.jsonl",
        "input_manifest": scratch_root / "input_manifest.tsv", "cgroup_start": h.cgroup_memory_snapshot(),
        "swap_start": h.system_swap_used(),
    }
    ctx["output"].mkdir(parents=True, exist_ok=True)
    preflight(ctx, "initial")
    log_event(ctx["state"], "RUN_STARTED", release_id=rid, run_id=args.run_id,
              no_genome_acquisition=True, host_only=True)
    write_input_manifest(ctx["input_manifest"], inputs)
    manifest_value = {"host_only_input_allowlist": list(h.HOST_ONLY_INPUT_ALLOWLIST),
                      "root_hashes": inputs["root_hashes"], "rows": h.EXPECTED_N,
                      "input_manifest_sha256": h.sha256_file(ctx["input_manifest"])}
    h.assert_host_only_manifest(manifest_value)
    h.write_json(ctx["output"] / "host_only_input_contract.json", manifest_value)
    h.write_json(ctx["output"] / "environment.json", {
        "schema": "host-structure-environment-v1", "platform": platform.platform(),
        "python": sys.version, "cpu_count": os.cpu_count(), "uid": os.getuid(), "gid": os.getgid(),
        "working_directory": str(repo), "locale": {k: os.environ.get(k, "") for k in ("LANG", "LC_ALL", "LC_CTYPE")},
        "task_id": h.TASK_ID, "run_id": args.run_id, "network_genome_acquisition": "BLOCKED",
    })
    tools = install_tools(ctx, repo, args.minimap2_source)
    qc = materialize_views(ctx, inputs, args.batch_size, args.inject_kill == "materialize")
    mash = run_mash(ctx, tools, inputs["accessions"], args.threads)
    high = run_high_fidelity(ctx, tools, inputs, mash, min(args.threads, 8))
    membership = freeze_memberships(ctx, inputs, mash, high, qc)
    metrics = {"mash": mash["metrics"], "high": high["summary"], "membership": membership,
               "qc_rows": len(qc), "parameters_sha256": hashlib.sha256(h.canonical_json(PARAMETERS)).hexdigest()}
    h.write_json(ctx["output"] / "metrics.json", metrics)
    log_event(ctx["state"], "READY_TO_PROMOTE", release_id=rid)
    final = copy_release_payload(ctx, inputs, metrics, inputs["root_hashes"])
    log_event(ctx["state"], "PROMOTED", release_id=rid, final=str(final))
    return {"release_id": rid, "external_release": str(final), "resumed_complete": False,
            "metrics": metrics}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--repo", type=Path, default=Path("."))
    p.add_argument("--durable-root", type=Path, default=h.DURABLE_PREFIX)
    p.add_argument("--scratch-root", type=Path, default=h.SCRATCH_PREFIX)
    p.add_argument("--run-id", default="host-structure-1000-v1-run")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--batch-size", type=int, default=50)
    p.add_argument("--assigned-ram-bytes", type=int, default=68_719_476_736)
    p.add_argument("--durable-allocation-bytes", type=int, default=30_000_000_000)
    p.add_argument("--scratch-allocation-bytes", type=int, default=4_000_000_000_000)
    p.add_argument("--inode-allocation", type=int, default=500_000)
    p.add_argument("--predicted-durable-upper95-bytes", type=int, default=10_000_000_000)
    p.add_argument("--predicted-scratch-upper95-bytes", type=int, default=20_000_000_000)
    p.add_argument("--predicted-files", type=int, default=10_000)
    p.add_argument("--unfinished-write-bytes", type=int, default=5_000_000_000)
    p.add_argument("--minimap2-source", type=Path,
                   default=Path("/home/erikg/micromamba/pkgs/minimap2-2.31-h118bc1c_0/bin/minimap2"))
    p.add_argument("--inject-kill", choices=("none", "materialize"), default="none")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = execute(args)
    except h.GateError as exc:
        print(f"NO_GO: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
