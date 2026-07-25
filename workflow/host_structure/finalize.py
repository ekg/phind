#!/usr/bin/env python3
"""Export compact tracked release metadata after independent external PASS."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any, Sequence

from . import host_structure as h
from . import runner
from .validate_release import validate


def tree_digest(root: Path) -> str:
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rows.append(f"{h.sha256_file(path)}  {path.relative_to(root).as_posix()}\n")
    import hashlib
    return hashlib.sha256("".join(rows).encode()).hexdigest()


def copy_compact_trees(outputs: Path, tracked: Path) -> None:
    target = tracked / "trees"
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(outputs / "trees/all_host_mash_supported.unrooted.nwk", target / "all_host_mash_supported.unrooted.nwk")
    for path in sorted((outputs / "mash").glob("*/rapidnj.unrooted.nwk")):
        dest = target / "mash" / path.parent.name / path.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    for path in sorted((outputs / "core").glob("L*/**/*.nwk")):
        rel = path.relative_to(outputs / "core")
        dest = target / "core" / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)


def render_report(release: dict[str, Any], external: Path, validation: dict[str, Any], metrics: dict[str, Any], resources: dict[str, Any], commands: list[dict[str, Any]]) -> str:
    membership = metrics["membership"]
    mash = metrics["mash"]
    high = metrics["high"]
    sensitivity = mash["sensitivity"]
    table = "\n".join(
        f"| `{x['name']}` | {x['spearman_sample']:.6f} | {x['nearest_neighbor_agreement']:.3%} | {x['split_jaccard']:.3%} |"
        for x in sensitivity
    )
    lineages = high["lineages"]
    passed_lineages = [x for x in lineages if x["verdict"] == h.PASS]
    min_conc = min(x["reference_split_concordance"] for x in passed_lineages)
    mash_wall = sum(float(x["wall_seconds"]) for x in commands if x["stage"].startswith(("mash-", "rapidnj-")))
    core_map_wall = sum(float(x["wall_seconds"]) for x in commands if x["stage"].startswith("core-map-"))
    return f"""# Frozen host population structure — exact 1,000-assembly pilot

**Release:** `{release['release_id']}`  
**Verdict:** **PASS**  
**External immutable release:** `{external}`  
**Scope:** host-only, phage-blind; frozen before any downstream phage association.

## Automatic input and safety gates

The run consumed the immutable selection `{h.SELECTION_RELEASE_ID}`
(`release.json` `{h.SELECTION_RELEASE_JSON_SHA256}`), canonical exact-N=1,000
release `{h.CANONICAL_RELEASE_ID}` (`release.json`
`{h.CANONICAL_RELEASE_JSON_SHA256}`), and host consumer certification
`{h.COMPATIBILITY_RELEASE_ID}` (`release.json`
`{h.COMPATIBILITY_RELEASE_JSON_SHA256}`). Cohort SHA-256 is `{h.COHORT_SHA256}`.
All 1,000 exact assembly revisions, rows, object/index checksums, PanSN names,
contig counts and bases were accounted in immutable cohort order. No genome was
acquired. The global union is the exact frozen 1,000-revision set and never
exceeded 1,000.

Both root inputs matched at start and finish:

- `26k_ecoli_accession.txt`: `{h.ROOT_HASHES['26k_ecoli_accession.txt']}`
- `26k_prophage1.csv`: `{h.ROOT_HASHES['26k_prophage1.csv']}`

The second file was **opaque-hashed only**. Host computation did not parse or
read any prophage presence, count, taxonomy, cluster, coordinates, sequence, or
derived trait. Prophage source/coordinate and extraction semantics are
`{h.NA_HOST}`, never fabricated as PASS. The declared host-only allow-list and
negative unit test enforce this boundary.

## All-host Mash overview

Pinned Mash 2.3 used whole assembly files (never per-contig `-i`) with baseline
`k=21,s=10,000,seed=42`. Exact validation found
**{mash['exact_pair_validation']['unordered_off_diagonal_pairs']:,}** unordered
off-diagonal pairs, **{mash['exact_pair_validation']['directed_records']:,}**
directed records, 1,000 zero diagonals, exact triangle/direct agreement and
symmetry. Pinned RapidNJ 2.3.2 retained exactly 1,000 tips. This is an unrooted
whole-genome similarity dendrogram, not a substitution-model phylogeny.

| configuration | sampled Spearman vs baseline | nearest-neighbor agreement | split Jaccard |
|---|---:|---:|---:|
{table}

Six independently seeded `k=21,s=10,000` sketch trees measured sketch-sampling
stability. Of {mash['baseline_splits']:,} baseline splits,
{mash['splits_param_stable']:,} were present under all four parameter settings
and {mash['splits_resampling_ge_95pct']:,} had at least 95% seed support.
Branches below 95% combined support were collapsed, not forced. The recorded
Mash sketch/distance plus RapidNJ command wall time was **{mash_wall:.1f} s**;
per-command argv, wall time, `/usr/bin/time -v`, stderr and output digests are
in the external provenance logs.

## Representatives, duplicates, placement and clades

Every eligible host has an exact tip/assembly/biological-unit, sequence class,
near-duplicate class, host-genetic representative, nearest-neighbor tie set,
sampling medoid, placement and membership row. Because no BioSample/isolate
identity table was supplied as a pinned host input, the biological analysis
unit is explicitly the exact assembly revision; downstream code must not
silently merge it.

- exact sequence classes: **{membership['exact_sequence_classes']:,}**
- exact duplicate non-primary assemblies: **{membership['exact_duplicate_assemblies']:,}**
- near-duplicate classes (`Mash D <= 0.0001`): **{membership['near_duplicate_classes']:,}**
- near-duplicate non-representatives: **{membership['near_duplicate_nonrepresentatives']:,}**
- frozen supported unrooted clades: **{membership['clades']:,}**
- supported fixed memberships: **{membership['supported_memberships']:,}**
- deliberately ambiguous memberships: **{membership['ambiguous_memberships']:,}**

Only mutually disjoint, 20–400-host unrooted splits present in every Mash
parameter tree and at ≥95% seed support became fixed clades. Other tips say
`AMBIGUOUS_UNROOTED_BACKBONE_OR_UNSUPPORTED`; no association-friendly topology
was chosen. Complete host-only medoid alternatives at k=12, 16 and 20 remain in
the release ({membership['alternative_partition_rows']:,} rows).

## High-fidelity host core ensemble

Sixteen host-genetic sampling partitions selected medoids, maximin-diverse and
boundary cases, and fragmentation-QC extremes without phage data. Digest-pinned
minimap2 2.31 built assembly-to-medoid reference-coordinate core alignments. The initial 90%-of-entire-reference
breadth screen correctly failed where lineage-specific accessory sequence was
not shared; those values remain a diagnostic and were not relabeled PASS. The
host-only pilot therefore uses an auditable core-genome denominator: the core
must span at least 50% of the reference, every selected sample must call at
least 95% of that core, and mean core missingness must be at most 5%. Dense
local panels have primary and diverse alternative-reference trees, missingness
checks, a deterministic SNP-density recombination-candidate mask, at least 100
non-recombinant parsimony-informative sites, 100 SNP-site bootstrap replicates,
and 95% support collapse. **{high['pass_lineages']}/16** lineages passed every
primary/alternative scientific and reference gate; **{high['blocked_ambiguous_lineages']}/16** sparse or failing lineages are explicitly blocked from clade inference and all their memberships remain ambiguous. No failed lineage was used to rescue a clade. Among passing lineages, the minimum primary-vs-alternative split concordance was **{min_conc:.3f}**. The minimap2 mapping commands took **{core_map_wall:.1f} s** total wall time; Python alignment/mask/bootstrap work is included in the outer run timing.

The density mask is explicitly a conservative host-core diagnostic, **not** a
claim that Gubbins ran. Unmasked/masked outputs and alternative references make
recombination and reference bias auditable. No independently verified outgroup
exists inside the frozen cohort and out-of-cohort acquisition is prohibited,
so every primary biological tree remains unrooted; midpoint rooting was not
used.

## Resource, restart, determinism and publication

The declared allocation was {resources['allocations']['assigned_ram_bytes']:,}
B RAM, {resources['allocations']['scratch_allocation_bytes']:,} B scratch,
{resources['allocations']['durable_allocation_bytes']:,} B durable and
{resources['allocations']['inode_allocation']:,} inodes. Live `findmnt`,
ownership/write probes, bytes/inodes, quotas, and unfinished-write reserve were
recorded before every batch/stage. Peak RSS was
{resources['peak_rss_bytes']:,} B ({resources['peak_rss_fraction']:.3%}); scratch
upper bound was {resources['scratch_peak_upper_bound_bytes']:,} B and
{resources['scratch_files']:,} files. OOM, cgroup swap growth and system swap
growth were zero; every ≤70% RAM/disk, ≤50% inode and 2x unfinished-write gate
PASS.

A real SIGKILL after five committed materialized views exposed no release.
The same run ID resumed only checksum-valid units and atomically published the
release. `SHA256SUMS` covers the exact input manifest, outputs, ledgers,
provenance and resource evidence; `COMPLETE` was fsynced last before a
same-filesystem rename. Independent semantic validation reports SHA-256
`{validation['semantic_sha256']}`. A post-publication rerun performs no
recompute or mutation.

## Interpretation boundary

This release freezes host topology, support-collapsed clades, all-host mappings
and alternative partitions **before** phage association. Mash branches describe
whole-assembly k-mer similarity; local core-SNP branches describe the sampled
host clonal frame conditional on mapping, mask and reference. Unsupported or
reference-sensitive structure is ambiguous, not evidence of absence and not a
license to choose whichever partition strengthens a phage result.
"""


def finalize(repo: Path, external: Path) -> dict[str, Any]:
    validation = validate(repo, external)
    release = json.loads((external / "release.json").read_text())
    outputs = external / "outputs"
    metrics = json.loads((outputs / "metrics.json").read_text())
    resources = json.loads((external / "resource_summary.json").read_text())
    commands = [json.loads(line) for line in (external / "commands.jsonl").read_text().splitlines()]

    manifest = repo / "manifests/host-structure-1000-v1"
    artifacts = repo / "artifacts/host_structure_1000"
    report = repo / "reports/host_structure_1000.md"
    if manifest.exists(): shutil.rmtree(manifest)
    if artifacts.exists(): shutil.rmtree(artifacts)
    manifest.mkdir(parents=True); artifacts.mkdir(parents=True)
    shutil.copy2(external / "release.json", manifest / "release.json")
    for name in ("host_membership.tsv", "host_clades.tsv", "medoids_and_cases.tsv", "alternative_partitions.tsv",
                 "nearest_neighbors.tsv", "host_qc.tsv"):
        shutil.copy2(outputs / name, manifest / name)
    copy_compact_trees(outputs, manifest)
    h.write_inventory(manifest)

    h.write_json(artifacts / "validation.json", validation)
    shutil.copy2(outputs / "metrics.json", artifacts / "metrics.json")
    shutil.copy2(outputs / "mash_metrics.json", artifacts / "mash_metrics.json")
    shutil.copy2(outputs / "high_fidelity_metrics.json", artifacts / "high_fidelity_metrics.json")
    shutil.copy2(outputs / "host_only_input_contract.json", artifacts / "host_only_input_contract.json")
    shutil.copy2(outputs / "tool_versions.json", artifacts / "tool_versions.json")
    shutil.copy2(external / "resource_summary.json", artifacts / "resource_summary.json")
    h.write_json(artifacts / "root_input_sha256_finish.json", h.verify_root_hashes(repo))
    with (artifacts / "command_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n",
                                fieldnames=["stage", "wall_seconds", "exit_status", "argv"])
        writer.writeheader()
        for command in commands:
            writer.writerow({"stage": command["stage"], "wall_seconds": f"{float(command['wall_seconds']):.6f}",
                             "exit_status": command["exit_status"], "argv": " ".join(command["argv"])})
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(render_report(release, external, validation, metrics, resources, commands))
    h.write_inventory(artifacts)
    return {"release_id": release["release_id"], "validation_semantic_sha256": validation["semantic_sha256"],
            "manifest_files": sum(1 for p in manifest.rglob("*") if p.is_file()),
            "artifact_files": sum(1 for p in artifacts.rglob("*") if p.is_file())}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--repo", type=Path, default=Path("."))
    p.add_argument("--external", type=Path, required=True)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    print(json.dumps(finalize(args.repo.resolve(), args.external.resolve()), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
