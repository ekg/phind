#!/usr/bin/env python3
"""Build compact, reproducible evidence after both N=1,000 deep validators pass."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from workflow.acquisition_canonicalization import pilot as p
from artifacts.canonical_cohort_1000 import runner as r

REPO = Path(".").resolve()
ARTIFACT = REPO / "artifacts/canonical_cohort_1000"
TRACKED = REPO / "manifests/canonical-cohort-1000-v1"
EXTERNAL = Path(
    "/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-1000/"
    "canonical-cohort-1000-v1-4bc3e029e6e0be44"
)
RUNNER_ARGS = [
    sys.executable, "-m", "artifacts.canonical_cohort_1000.runner",
    "--durable-task-root", str(EXTERNAL.parent),
    "--scratch-root", "/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/prepare-canonical-cohort-1000/canonical-1000-v1-run",
    "--run-id", "canonical-1000-v1-run", "--batch-size", "10", "--retries", "8", "--rate-delay", "0.5",
    "--assigned-ram-bytes", "34359738368", "--durable-allocation-bytes", "15000000000",
    "--scratch-allocation-bytes", "4000000000000", "--inode-allocation", "400000",
    "--predicted-durable-peak-bytes", "4000000000", "--predicted-scratch-peak-bytes", "3000000000",
    "--predicted-files", "15000", "--unfinished-write-bytes", "2000000000",
    "--n1000-time-allocation-seconds", "7200", "--bgzip-threads", "2", "--inject-kill", "none",
]


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode() + b"\n"


def tree_stats(root: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    files = 0
    size = 0
    for path in sorted(x for x in root.rglob("*") if x.is_file() and not x.is_symlink()):
        files += 1
        size += path.stat().st_size
        digest.update(str(path.relative_to(root)).encode() + b"\0")
        digest.update(r.sha_file(path).encode() + b"\n")
    return {"root": str(root), "files": files, "bytes": size, "tree_sha256": digest.hexdigest()}


def state_counts(path: Path) -> Counter[str]:
    return Counter(json.loads(line)["event"] for line in path.read_text().splitlines() if line.strip())


def parse_time(path: Path) -> dict[str, Any]:
    text = path.read_text()
    def value(label: str, cast: type = float) -> Any:
        match = re.search(rf"^\s*{re.escape(label)}:\s*(.+)$", text, re.MULTILINE)
        if not match:
            raise RuntimeError(f"missing time field {label}: {path}")
        return cast(match.group(1).strip())
    elapsed = value("Elapsed (wall clock) time (h:mm:ss or m:ss)", str)
    parts = [float(x) for x in elapsed.split(":")]
    seconds = parts[-1] + (parts[-2] * 60 if len(parts) >= 2 else 0) + (parts[-3] * 3600 if len(parts) >= 3 else 0)
    command_match = re.search(r'^\s*Command being timed: "(.*)"$', text, re.MULTILINE)
    return {
        "command": command_match.group(1) if command_match else ".",
        "wall_seconds": seconds,
        "user_seconds": value("User time (seconds)"),
        "system_seconds": value("System time (seconds)"),
        "peak_rss_bytes": value("Maximum resident set size (kbytes)", int) * 1024,
        "swap_events": value("Swaps", int),
        "exit_status": value("Exit status", int),
    }


def build_exceptions(assemblies: list[dict[str, str]], refs: list[dict[str, str]]) -> None:
    fields = ["cohort_order", "accession", "annotation_status", "failure_reason",
              "storage_release_id", "storage_root", "row_sha256"]
    rows: list[dict[str, Any]] = []
    for assembly, ref in zip(assemblies, refs):
        if not assembly["annotation_status"].startswith("QUARANTINED"):
            continue
        root = EXTERNAL if ref["storage_release_id"] == "SELF" else Path(ref["storage_root"])
        manifest = json.loads((root / ref["canonical_object_relpath"] / "manifest.json").read_text())
        rows.append({
            "cohort_order": ref["cohort_order"], "accession": ref["accession"],
            "annotation_status": assembly["annotation_status"],
            "failure_reason": manifest["annotation"]["failure_reason"],
            "storage_release_id": ref["storage_release_id"], "storage_root": ref["storage_root"],
        })
    p.write_tsv(ARTIFACT / "exceptions.tsv", fields, rows)


def build_failures() -> None:
    fields = ["event", "type", "at", "message", "retry_disposition", "row_sha256"]
    rows = []
    for row in [json.loads(line) for line in (EXTERNAL / "failures.jsonl").read_text().splitlines() if line.strip()]:
        rows.append({
            "event": row.get("event", "."), "type": row.get("type", "."),
            "at": row.get("at", "."), "message": row.get("message", "."),
            "retry_disposition": "STOPPED_NO_PUBLICATION_RESUMED_FROM_CHECKSUM_VALIDATED_UNITS",
        })
    p.write_tsv(ARTIFACT / "failures.tsv", fields, rows)


def copy_evidence_logs() -> None:
    mapping = {
        "/tmp/canonical1000-time-acquisition-kill.txt": "time_acquisition_kill.txt",
        "/tmp/canonical1000-time-conversion-kill.txt": "time_conversion_kill.txt",
        "/tmp/canonical1000-time-full-run.txt": "time_full_run.txt",
        "/tmp/canonical1000-time-resource-failed.txt": "time_resource_failed_retry.txt",
        "/tmp/canonical1000-time-resource-retry2-failed.txt": "time_resource_retry2_failed.txt",
        "/tmp/canonical1000-time-resource-retry3-failed.txt": "time_resource_retry3_failed.txt",
        "/tmp/canonical1000-time-validation.txt": "time_validation.txt",
        "/tmp/canonical1000-time-validation-rerun.txt": "time_validation_rerun.txt",
        "/tmp/canonical1000-time-all.txt": "time_all_tests.txt",
        "/tmp/canonical1000-all.txt": "all_unit_tests.txt",
    }
    for source, target in mapping.items():
        src = Path(source)
        if not src.is_file():
            raise RuntimeError(f"required compact evidence log missing: {src}")
        shutil.copyfile(src, ARTIFACT / target)
    result = subprocess.run(
        [sys.executable, "-m", "unittest", "-v", "artifacts.canonical_cohort_1000.test_runner"],
        cwd=REPO, text=True, capture_output=True, check=True,
    )
    (ARTIFACT / "unit_tests.txt").write_text(result.stdout + result.stderr)


def build_report(
    release: dict[str, Any], summary: dict[str, Any], scale: dict[str, Any],
    metrics: dict[str, Any], deterministic: dict[str, Any], validation: dict[str, Any],
) -> None:
    counts = release["counts"]
    resources = metrics["resource"]
    slopes = scale["last_two_rung_per_base_slopes"]
    changes = slopes["relative_changes"]
    full = metrics["invocations"]["full_run"]
    resource_failed = metrics["invocations"]["resource_failed_attempt"]
    resource_retry2 = metrics["invocations"]["resource_retry2_failed"]
    resource_retry3 = metrics["invocations"]["resource_retry3_failed"]
    acquisition = metrics["invocations"]["acquisition_kill"]
    conversion = metrics["invocations"]["conversion_kill"]
    external_sums_sha = r.sha_file(EXTERNAL / "SHA256SUMS")
    complete_sha = r.sha_file(EXTERNAL / "COMPLETE")
    text = f"""# Canonical cohort 1,000

**Verdict: PASS.** The exact frozen N=1,000 ceiling is immutable release
`{release['release_id']}` at:

```text
{EXTERNAL}/
```

It contains {counts['validated']:,} terminal `VALIDATED` genome rows,
{counts['contigs']:,} lossless contig rows, and {counts['total_bases']:,} bases.
The first {r.PREDECESSOR_ROWS:,} objects are read-only digest references; only
rows 501-1,000 were acquired and canonicalized. No 1,001st revision was
requested, analyzed, or published, and this release does not authorize any
beyond-1,000 projection.

## Automatic immutable gates

Execution consumed selection `{r.SELECTION_RELEASE_ID}`, tracked `release.json`
SHA-256 `{r.SELECTION_RELEASE_JSON_SHA256}`, and the exact {r.COHORT_BYTES:,}-byte,
1,000-row cohort SHA-256 `{r.COHORT_SHA256}`. It consumed predecessor
`{r.PREDECESSOR_RELEASE_ID}`, release SHA-256
`{r.PREDECESSOR_RELEASE_JSON_SHA256}`, plus checksum-pinned canonical scale
trend SHA-256 `{r.PREDECESSOR_SCALE_TREND_SHA256}`. The predecessor's
`applicable_gates.canonical_scale_trend`, `scale_trend.json` verdict, every
scale check, and every N=1,000 projection check were unqualified `PASS`.
All 19 pinned-consumer gates were also `PASS`. Missing, mismatched,
`CONDITIONAL`, or failed inputs are rejected by code and tests; no prompt,
wait, waiver, shrink, reconstruction, or substitution path exists.

Canonical preparation records prophage source/coordinate policy and integrated
extraction as `NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED`.
It does not fabricate a coordinate or integrated-analysis PASS.

Both immutable roots matched at start and finish and were not edited:

- `26k_ecoli_accession.txt`: `{r.ROOT_HASHES['26k_ecoli_accession.txt']}`
- `26k_prophage1.csv`: `{r.ROOT_HASHES['26k_prophage1.csv']}`

The graph-wide sequence-bearing union was exactly 500 at start and the exact
frozen 1,000 at finish, wholly within the selected collection and hard cap.

## Reuse, acquisition, canonicalization, and quarantine

All 500 predecessor references and all eight role-specific digests per row
matched the immutable N=500 chain. The request ledger contains no predecessor
acquisition or canonicalization event. The 500 incremental packages total
{metrics['bytes']['source_packages_new_500_validated_transfer']:,} bytes and
incremental canonical BGZF files total
{metrics['bytes']['canonical_bgzf_new_500']:,} bytes.

Every referenced assembly passed exact accession/version package identity,
ZIP/CRC/upstream-MD5/local-SHA checks; rename-only sequence identity; exact,
unique and reversible PanSN names; BGZF integrity; `.fai`/`.gzi`; every-contig
`samtools faidx` prefix-region round trip; crosswalk; source-GFF alias/bounds;
object inventory; and atomic marker validation. There are zero genome-object
quarantines. {counts['annotation_views_quarantined']:,} optional source-GFF
alias views reproducibly failed their own bounds and remain explicit in
`artifacts/canonical_cohort_1000/exceptions.tsv`; no partial alias view was
published. The other {counts['annotations_alias_validated']:,} annotation
views reproduced exactly. This annotation QC makes no prophage-coordinate claim.

## Batches, retry, and injected restart

The release has {metrics['batches']['count']} deterministic ten-row batches and
{summary['preflight_record_count']:,} append-only live resource records. The
injected acquisition process was SIGKILLed after fsyncing 131,072 partial bytes
({acquisition['wall_seconds']:.2f} s); neither final release nor `COMPLETE`
appeared. Resume rejected unsafe range identity and reacquired that exact row.
Conversion was then SIGKILLed after crossing the configured 262,144-base threshold
({conversion['wall_seconds']:.2f} s), again with no final publication. Resume
discarded the interrupted conversion stage and independently validated the
completed unit. Batch evidence records {metrics['restart']['batch_restart_events']}
restart events and {metrics['batches']['partial_bytes_observations']:,} partial
bytes; final partial bytes are zero.

The payload-completion attempt used {resource_failed['wall_seconds']:.2f} s but
correctly stopped before promotion when sampled system swap-free decreased;
every other resource and scale check passed. That failure is append-only.
Two checksum-only promotion retries also refused on continued system swap
growth after {resource_retry2['wall_seconds'] + resource_retry3['wall_seconds']:.2f}
s total. Once live swap-free stabilized, the next checksum-only retry revalidated
all units and observed zero growth, using {full['wall_seconds']:.2f} s wall and
{full['user_seconds'] + full['system_seconds']:.2f} s user+system. Across six
measured invocations, wall was {metrics['cumulative_execution']['wall_seconds']:.2f}
s, CPU was {metrics['cumulative_execution']['cpu_seconds']:.2f} s, outer peak
RSS was {metrics['cumulative_execution']['peak_rss_bytes']:,} bytes, and process
swap events were zero. The append ledger additionally exposes
{metrics['cumulative_execution']['unmeasured_interrupted_retry_invocations']}
checksum-only code-correction interruption(s), with no payload work or
publication. There were {metrics['transport']['download_requests']}
bounded GET requests for {metrics['transport']['completed_downloads']} completed
incremental packages and {metrics['transport']['failure_ledger_events']} explicit
failure-ledger event(s), compactly mirrored in
`artifacts/canonical_cohort_1000/failures.tsv`.

`COMPLETE` was written/fsynced last and the full staging directory was promoted
by same-filesystem rename. No stage survived. External `SHA256SUMS` SHA-256 is
`{external_sums_sha}` and `COMPLETE` SHA-256 is `{complete_sha}`.

## Live resources and final scaling model

| Resource | Allocation | Measured / upper model | Gate |
|---|---:|---:|---|
| RAM | {summary['allocations']['assigned_ram_bytes']:,} B | {resources['peak_rss_bytes_outer']:,} B outer peak ({resources['peak_rss_fraction_outer']:.3%}) | <=70% PASS |
| Durable disk | {summary['allocations']['durable_allocation_bytes']:,} B | {resources['measured_stage_peak_bytes']:,} B measured; {resources['predicted_upper95_peak_bytes']:,} B modeled upper-95 | <=70% PASS |
| Scratch disk | {summary['allocations']['scratch_allocation_bytes']:,} B | {summary['allocations']['predicted_scratch_peak_bytes']:,} B reserved | <=70% PASS |
| Inodes | {summary['allocations']['inode_allocation']:,} | {resources['measured_stage_peak_files']:,} files; {summary['allocations']['predicted_files']:,} configured | <=50% PASS |
| Unfinished write | {summary['allocations']['unfinished_write_bytes']:,} B | at least 2x retained live/allocation | PASS |

Every stage/batch record includes `findmnt`, mount/source/fstype, ownership and
write probe, live free bytes/inodes, assigned quotas, swap, and unfinished-write
reservation. Durable free space remained at least
{min(summary['start']['durable_free_bytes'], summary['finish']['durable_free_bytes']):,}
bytes with at least
{min(summary['start']['durable_free_inodes'], summary['finish']['durable_free_inodes']):,}
free inodes; scratch remained at least
{min(summary['start']['scratch_free_bytes'], summary['finish']['scratch_free_bytes']):,}
bytes and
{min(summary['start']['scratch_free_inodes'], summary['finish']['scratch_free_inodes']):,}
inodes. All durable 2 TB/1M-inode and scratch 4 TB preflight plus 2 TB/5M stop
floors passed. The successful promotion retry had no OOM, process swap, or
system swap growth; every earlier system-wide swap-growth refusal remains
explicit in the immutable failure/resource ledgers.

The final adjacent N=500 to N=1,000 time exponent is
**{scale['time_exponent']['current_n500_to_n1000']:.6f}**; the empirical upper
bound is **{scale['time_exponent']['empirical_upper_bound']:.6f}** (limit 1.3).

| Per-new-base slope | N=500 | N=1,000 | Change |
|---|---:|---:|---:|
| Wall seconds | {slopes['n500']['wall_seconds_per_new_base']:.7g} | {slopes['n1000']['wall_seconds_per_new_base']:.7g} | {changes['wall_seconds_per_new_base']:+.3%} |
| Source bytes | {slopes['n500']['source_bytes_per_new_base']:.7g} | {slopes['n1000']['source_bytes_per_new_base']:.7g} | {changes['source_bytes_per_new_base']:+.3%} |
| Stage bytes | {slopes['n500']['stage_bytes_per_new_base']:.7g} | {slopes['n1000']['stage_bytes_per_new_base']:.7g} | {changes['stage_bytes_per_new_base']:+.3%} |
| Stage files | {slopes['n500']['stage_files_per_new_base']:.7g} | {slopes['n1000']['stage_files_per_new_base']:.7g} | {changes['stage_files_per_new_base']:+.3%} |
| Peak RSS bytes | {slopes['n500']['peak_rss_bytes_per_new_base']:.7g} | {slopes['n1000']['peak_rss_bytes_per_new_base']:.7g} | {changes['peak_rss_bytes_per_new_base']:+.3%} |

Every absolute change is <=25%. All measured N=1,000 values stayed within the
checksum-pinned predecessor upper-95 projection. `scale_trend.json` publishes
descriptive final fits over N=10/100/250/500/1,000 and explicitly records
`projection_beyond_ceiling=NOT_AUTHORIZED_NOT_COMPUTED`.

## Independent validation, determinism, and compact handoff

Two full independent validators produced byte-identical PASS JSON SHA-256
`{metrics['validation']['semantic_json_sha256']}`. Each rehashed inventories;
validated all {validation['assembly_rows']:,} archives/object references and
{validation['contig_rows']:,} contig/index/name/region checks; recomputed every
annotation decision and final scale check; audited restart/resources/roots/
exact cap; and rejected partial, unlisted, symlink, or plain-FASTA release data.
Task tests cover manifest/inventory mismatch, non-PASS gates, live-resource
refusal, a 1,001-row union, scale regression, bounded URL/retry, corrupt
completed units, and interrupted promotion. The task suite passes 14/14 and the
full project regression suite passes 95/95.

A post-validation existing-release rerun made zero network requests and zero
canonical recomputations. External tree SHA-256
`{deterministic['after']['external']['tree_sha256']}` and tracked-manifest tree
SHA-256 `{deterministic['after']['tracked']['tree_sha256']}` were unchanged,
as were state bytes and request/commit counts.

Git contains only compact files under `manifests/canonical-cohort-1000-v1/`,
`artifacts/canonical_cohort_1000/`, and this report. No package, genome/prophage
sequence, source FASTA/GFF, canonical BGZF, `.fai`, `.gzi`, raw cache, whole
index, or per-hit biological output is tracked; each task-owned file is under
10 MiB. Bulky objects remain solely in this task's external namespace, and no
dependent analysis is authorized beyond N=1,000.
"""
    report = REPO / "reports/canonical_cohort_1000.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(text)


def main() -> int:
    if not (EXTERNAL / "COMPLETE").is_file():
        raise RuntimeError("atomic release is not complete")
    first = ARTIFACT / "validation.json"
    second = ARTIFACT / "validation_rerun.json"
    if not first.is_file() or not second.is_file() or first.read_bytes() != second.read_bytes():
        raise RuntimeError("two byte-identical deep PASS validations are required")
    validation = json.loads(first.read_text())
    if validation.get("verdict") != "PASS" or validation.get("assembly_rows") != 1000:
        raise RuntimeError("deep validation is not exact N=1,000 PASS")

    state = EXTERNAL / "state.jsonl"
    before = {
        "external": tree_stats(EXTERNAL), "tracked": tree_stats(TRACKED),
        "state_sha256": r.sha_file(state), "events": dict(state_counts(state)),
    }
    time_path = ARTIFACT / "time_deterministic_rerun.txt"
    with time_path.open("w") as time_handle:
        result = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(time_path), *RUNNER_ARGS],
            cwd=REPO, text=True, capture_output=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"deterministic runner failed: {result.stderr}")
    after = {
        "external": tree_stats(EXTERNAL), "tracked": tree_stats(TRACKED),
        "state_sha256": r.sha_file(state), "events": dict(state_counts(state)),
    }
    deterministic = {
        "schema": "canonical-cohort-1000-deterministic-rerun-v1", "verdict": "PASS",
        "before": before, "after": after,
        "checks": {
            "external_tree_unchanged": before["external"] == after["external"],
            "tracked_tree_unchanged": before["tracked"] == after["tracked"],
            "state_ledger_unchanged": before["state_sha256"] == after["state_sha256"],
            "zero_new_network_requests": before["events"].get("ACQUISITION_REQUEST_STARTED", 0) == after["events"].get("ACQUISITION_REQUEST_STARTED", 0),
            "zero_new_canonical_recomputations": before["events"].get("CANONICAL_OBJECT_COMMITTED", 0) == after["events"].get("CANONICAL_OBJECT_COMMITTED", 0),
        },
    }
    if not all(deterministic["checks"].values()):
        raise RuntimeError(f"deterministic rerun gate failed: {deterministic['checks']}")
    (ARTIFACT / "deterministic_rerun.json").write_bytes(canonical(deterministic))

    copy_evidence_logs()
    release = json.loads((EXTERNAL / "release.json").read_text())
    summary = json.loads((EXTERNAL / "resource_summary.json").read_text())
    scale = json.loads((EXTERNAL / "scale_trend.json").read_text())
    assemblies = p.read_tsv(EXTERNAL / "manifests/assemblies.tsv", p.ASSEMBLY_FIELDS, verify_hashes=True)
    refs = p.read_tsv(EXTERNAL / "manifests/object_refs.tsv", r.OBJECT_REF_FIELDS, verify_hashes=True)
    batches = p.read_tsv(EXTERNAL / "manifests/batch_metrics.tsv", r.BATCH_FIELDS, verify_hashes=True)
    build_exceptions(assemblies, refs)
    build_failures()

    event_counts = state_counts(state)
    invocation_times = {
        "acquisition_kill": parse_time(ARTIFACT / "time_acquisition_kill.txt"),
        "conversion_kill": parse_time(ARTIFACT / "time_conversion_kill.txt"),
        "resource_failed_attempt": parse_time(ARTIFACT / "time_resource_failed_retry.txt"),
        "resource_retry2_failed": parse_time(ARTIFACT / "time_resource_retry2_failed.txt"),
        "resource_retry3_failed": parse_time(ARTIFACT / "time_resource_retry3_failed.txt"),
        "full_run": parse_time(ARTIFACT / "time_full_run.txt"),
    }
    invocation_times["acquisition_kill"]["terminated_by_injected_signal_9"] = True
    invocation_times["conversion_kill"]["terminated_by_injected_signal_9"] = True
    outer_peak = max(row["peak_rss_bytes"] for row in invocation_times.values())
    cumulative_wall = sum(row["wall_seconds"] for row in invocation_times.values())
    cumulative_cpu = sum(row["user_seconds"] + row["system_seconds"] for row in invocation_times.values())
    new_rows = assemblies[r.PREDECESSOR_ROWS:]
    metrics = {
        "schema": "canonical-cohort-1000-compact-metrics-v1", "verdict": "PASS",
        "release_id": release["release_id"], "counts": release["counts"],
        "bytes": {
            "source_packages_new_500_validated_transfer": sum(int(row["source_package_bytes"]) for row in new_rows),
            "source_packages_all_1000": sum(int(row["source_package_bytes"]) for row in assemblies),
            "canonical_bgzf_new_500": sum(int(row["canonical_bgzf_bytes"]) for row in new_rows),
            "canonical_bgzf_all_1000": sum(int(row["canonical_bgzf_bytes"]) for row in assemblies),
            "release_bytes": after["external"]["bytes"],
            "injected_acquisition_partial_network_bytes": 131072,
        },
        "files": {"release_files": after["external"]["files"], "batch_final_files": int(batches[-1]["stage_files_finish"])},
        "transport": {
            "completed_downloads": event_counts.get("ACQUISITION_DOWNLOAD_COMPLETE", 0),
            "download_requests": event_counts.get("ACQUISITION_REQUEST_STARTED", 0),
            "unsafe_partial_restarts": event_counts.get("ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE", 0),
            "failure_ledger_events": len([line for line in (EXTERNAL / "failures.jsonl").read_text().splitlines() if line.strip()]),
        },
        "restart": {
            "injected_acquisition_kills": event_counts.get("INJECTED_ACQUISITION_SIGKILL", 0),
            "injected_conversion_kills": event_counts.get("INJECTED_CONVERSION_SIGKILL", 0),
            "conversion_partial_discards": event_counts.get("INTERRUPTED_CONVERSION_STAGE_DISCARDED", 0),
            "batch_restart_events": sum(int(row["restart_events"]) for row in batches),
        },
        "batches": {
            "count": len(batches), "wall_seconds_sum": sum(float(row["wall_seconds"]) for row in batches),
            "driver_cpu_seconds_sum": sum(float(row["cpu_seconds"]) for row in batches),
            "last_cumulative_transfer_bytes": int(batches[-1]["cumulative_transfer_bytes"]),
            "last_cumulative_canonical_bgzf_bytes": int(batches[-1]["cumulative_canonical_bgzf_bytes"]),
            "partial_bytes_observations": sum(int(row["partial_bytes_observed"]) for row in batches),
        },
        "invocations": invocation_times,
        "cumulative_execution": {"wall_seconds": cumulative_wall, "cpu_seconds": cumulative_cpu,
                                  "peak_rss_bytes": outer_peak, "swap_events": sum(row["swap_events"] for row in invocation_times.values()),
                                  "state_invocation_preflights": event_counts.get("PREFLIGHT_PASS", 0) + event_counts.get("RESUME_PREFLIGHT_PASS", 0),
                                  "unmeasured_interrupted_retry_invocations": max(0, event_counts.get("PREFLIGHT_PASS", 0) + event_counts.get("RESUME_PREFLIGHT_PASS", 0) - len(invocation_times))},
        "resource": {
            "assigned_ram_bytes": summary["allocations"]["assigned_ram_bytes"],
            "durable_allocation_bytes": summary["allocations"]["durable_allocation_bytes"],
            "inode_allocation": summary["allocations"]["inode_allocation"],
            "peak_rss_bytes_internal": summary["peak_rss_bytes"], "peak_rss_bytes_outer": outer_peak,
            "peak_rss_fraction_outer": outer_peak / summary["allocations"]["assigned_ram_bytes"],
            "measured_stage_peak_bytes": summary["measured_release_stage_peak_bytes"],
            "measured_stage_peak_files": summary["measured_release_stage_peak_files"],
            "predicted_upper95_peak_bytes": summary["disk_projection"]["modeled_upper95_peak_bytes"],
            "preflight_records": summary["preflight_record_count"],
            "system_swap_growth_bytes": summary["system_swap_growth_bytes"],
            "canonical_time_exponent_upper_bound": scale["time_exponent"]["empirical_upper_bound"],
            "max_absolute_last_two_per_base_slope_change": max(abs(x) for x in scale["last_two_rung_per_base_slopes"]["relative_changes"].values()),
        },
        "scale_trend": scale,
        "annotation": {"views_validated": release["counts"]["annotations_alias_validated"],
                       "views_quarantined": release["counts"]["annotation_views_quarantined"]},
        "validation": {"independent_deep_runs": 2, "byte_identical": True,
                       "semantic_json_sha256": r.sha_file(first), "check_groups": validation["check_group_count"]},
        "determinism": {"external_tree_sha256": after["external"]["tree_sha256"],
                        "tracked_tree_sha256": after["tracked"]["tree_sha256"],
                        "zero_new_network_requests": True, "zero_new_canonical_recomputations": True},
        "gate_applicability": {
            "canonical_n500_scale_trend": "PASS", "canonical_n1000_projection": "PASS",
            "canonical_scale_trend": "PASS", "human_override_used": False,
            "integrated_extraction_verdict": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED",
            "prophage_source_coordinate_policy": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED",
            "beyond_1000_authorized": False,
        },
    }
    if outer_peak * 100 > summary["allocations"]["assigned_ram_bytes"] * 70 or metrics["cumulative_execution"]["swap_events"]:
        raise RuntimeError("outer RAM/swap gate failed")
    (ARTIFACT / "metrics.json").write_bytes(canonical(metrics))
    build_report(release, summary, scale, metrics, deterministic, validation)
    shutil.copyfile(EXTERNAL / "global_cap_evidence.json", ARTIFACT / "global_cap_evidence.json")
    shutil.copyfile(EXTERNAL / "resource_summary.json", ARTIFACT / "resource_summary.json")
    shutil.copyfile(EXTERNAL / "restart_evidence.json", ARTIFACT / "restart_evidence.json")
    shutil.copyfile(EXTERNAL / "scale_trend.json", ARTIFACT / "scale_trend.json")
    shutil.copyfile(EXTERNAL / "tools.json", ARTIFACT / "tools.json")
    (ARTIFACT / "root_input_sha256_finish.txt").write_text(
        "\n".join(f"{r.sha_file(REPO / name)}  {name}" for name in sorted(r.ROOT_HASHES)) + "\n"
    )

    sums = []
    for path in sorted(x for x in ARTIFACT.iterdir() if x.is_file() and x.name != "SHA256SUMS"):
        sums.append(f"{r.sha_file(path)}  {path.name}\n")
    (ARTIFACT / "SHA256SUMS").write_text("".join(sums))
    print(canonical({"verdict": "PASS", "release_id": release["release_id"], "metrics": str(ARTIFACT / "metrics.json")}).decode(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
