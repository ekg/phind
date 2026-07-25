"""Read-only, fail-closed certification of the Phigaro v2.4.0 N=10 release."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

EXPECTED_ACCESSIONS = [
    "GCF_000005845.2",
    "GCF_000812325.1",
    "GCF_002302315.1",
    "GCF_004664255.1",
    "GCF_015644385.1",
    "GCF_020829045.1",
    "GCF_921380995.1",
    "GCF_000167895.3",
    "GCF_001881595.4",
    "GCF_000498835.2",
]
DEFAULT_REFERENCE = "artifacts/phigaro_v2_4_pilot/release_reference.json"
PRODUCER_VALIDATOR = "workflow/phigaro_v2_4_pilot/validate_release.py"
DEFAULT_RELEASE_PARENT = Path(
    "/home/erikg/phind-data/ecoli26k/v1/releases/complete-phigaro-v2-4-n10-pilot"
)
REQUIRED_RELEASE_FILES = [
    "COMPLETE",
    "SHA256SUMS",
    "input_manifest.tsv",
    "tool_versions.json",
    "database.json",
    "config.json",
    "environment.json",
    "processes.jsonl",
    "outcomes.tsv",
    "resources.json",
    "deterministic_rerun.json",
    "restart_evidence.json",
    "root_input_sha256_start.json",
    "root_input_sha256_finish.json",
]
DEEP_GATES = [
    "release_complete_and_inventory_digests",
    "exact_ordered_n10_release_manifest",
    "phigaro_2_4_0_dependency_database_config_pins",
    "real_argv_and_process_logs_for_ten_genomes",
    "ten_outcomes_and_call_count_reconciliation",
    "native_tsv_bed_gff_saved_fasta_outputs",
    "interval_and_base_round_trips",
    "resource_evidence",
    "deterministic_and_restart_evidence",
    "published_root_digest_evidence",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_file_sha256(repo: Path, relpath: str) -> str | None:
    proc = subprocess.run(
        ["git", "-C", str(repo), "show", f"HEAD:{relpath}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return hashlib.sha256(proc.stdout).hexdigest() if proc.returncode == 0 else None


def git_tracked(repo: Path, relpath: str) -> bool:
    return subprocess.run(
        ["git", "-C", str(repo), "ls-files", "--error-unmatch", "--", relpath],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    ).returncode == 0


def ordered_cohort(manifest: Path) -> list[str]:
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if [int(row["cohort_order"]) for row in rows] != list(range(1, 11)):
        raise ValueError("cohort_order is not exactly 1..10")
    if any(row["rung_n"] != "10" for row in rows):
        raise ValueError("not every canonical cohort row has rung_n=10")
    return [row["exact_assembly_accession_version"] for row in rows]


def parse_sha256sums(path: Path) -> tuple[dict[str, str], list[str]]:
    entries: dict[str, str] = {}
    errors: list[str] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw:
            continue
        fields = raw.split(maxsplit=1)
        if (
            len(fields) != 2
            or len(fields[0]) != 64
            or any(char not in "0123456789abcdefABCDEF" for char in fields[0])
        ):
            errors.append(f"malformed SHA256SUMS line {number}")
            continue
        relpath = fields[1].lstrip("*")
        if relpath.startswith("/") or ".." in Path(relpath).parts or relpath in entries:
            errors.append(f"unsafe or duplicate SHA256SUMS path on line {number}")
            continue
        entries[relpath] = fields[0].lower()
    return entries, errors


def verify_inventory(release: Path) -> tuple[bool, dict[str, Any]]:
    inventory = release / "SHA256SUMS"
    entries, errors = parse_sha256sums(inventory)
    actual_files = {
        path.relative_to(release).as_posix()
        for path in release.rglob("*")
        if path.is_file() and path.name not in {"SHA256SUMS", "COMPLETE"}
    }
    if set(entries) != actual_files:
        errors.append("SHA256SUMS coverage does not exactly equal release files excluding markers")
    mismatches: list[str] = []
    for relpath, expected in sorted(entries.items()):
        candidate = release / relpath
        if candidate.is_symlink() or not candidate.is_file() or sha256(candidate) != expected:
            mismatches.append(relpath)
    return not errors and not mismatches, {
        "inventory_rows": len(entries),
        "actual_inventory_files": len(actual_files),
        "inventory_errors": errors,
        "digest_mismatches": mismatches,
        "sha256sums_sha256": sha256(inventory),
    }


def add_failure(result: dict[str, Any], code: str, detail: str) -> None:
    result["failures"].append({"code": code, "detail": detail})


def certify(repo: Path, reference_relpath: str, release_parent: Path) -> dict[str, Any]:
    repo = repo.resolve()
    canonical_relpath = "manifests/pilot-cohorts-v1/cohort-0010.tsv"
    canonical_path = repo / canonical_relpath
    root_relpaths = ["26k_ecoli_accession.txt", "26k_prophage1.csv"]
    root_start = {name: sha256(repo / name) for name in root_relpaths}
    canonical_accessions = ordered_cohort(canonical_path)

    producer_validator = repo / PRODUCER_VALIDATOR
    releases_root = release_parent.parent
    external_candidates = sorted(
        str(path)
        for parent in ([path for path in releases_root.iterdir() if path.is_dir()] if releases_root.is_dir() else [])
        for path in [parent, *[child for child in parent.iterdir() if child.is_dir()]]
        if "phigaro" in path.name.lower()
    )
    result: dict[str, Any] = {
        "schema": "phigaro-v2.4-n10-independent-certification-v1",
        "verdict": "REJECTED",
        "qualified_for_scaling": False,
        "audit_mode": "READ_ONLY_FAIL_CLOSED",
        "source_task_id": "complete-real-phigaro",
        "reference": {
            "relpath": reference_relpath,
            "exists": (repo / reference_relpath).is_file(),
            "git_tracked": git_tracked(repo, reference_relpath),
        },
        "expected_release_parent": str(release_parent),
        "bounded_discovery": {
            "producer_validator_relpath": PRODUCER_VALIDATOR,
            "producer_validator_exists": producer_validator.is_file(),
            "producer_validator_git_tracked": git_tracked(repo, PRODUCER_VALIDATOR),
            "expected_release_parent_exists": release_parent.is_dir(),
            "phigaro_named_release_candidates_max_depth_2": external_candidates,
        },
        "canonical_n10": {
            "manifest_relpath": canonical_relpath,
            "manifest_sha256": sha256(canonical_path),
            "row_count": len(canonical_accessions),
            "ordered_accessions": canonical_accessions,
            "matches_frozen_identity_and_order": canonical_accessions == EXPECTED_ACCESSIONS,
        },
        "checks": {},
        "failures": [],
        "side_effects": {
            "release_writes": 0,
            "host_assembly_downloads": 0,
            "larger_cohorts_executed": 0,
            "historical_csv_semantics_interpreted": False,
        },
    }
    checks = result["checks"]
    producer_validator_ok = producer_validator.is_file() and git_tracked(repo, PRODUCER_VALIDATOR)
    checks["producer_strict_validator_available"] = "PASS" if producer_validator_ok else "FAIL"
    if not producer_validator_ok:
        add_failure(result, "MISSING_PRODUCER_STRICT_VALIDATOR", f"missing tracked {PRODUCER_VALIDATOR}")
    release_parent_ok = release_parent.is_dir()
    checks["expected_external_release_parent_present"] = "PASS" if release_parent_ok else "FAIL"
    if not release_parent_ok:
        add_failure(result, "MISSING_EXPECTED_RELEASE_PARENT", str(release_parent))

    canonical_ok = canonical_accessions == EXPECTED_ACCESSIONS
    checks["independent_canonical_n10_basis"] = "PASS" if canonical_ok else "FAIL"
    if not canonical_ok:
        add_failure(result, "CANONICAL_N10_MISMATCH", "tracked frozen N=10 identity/order differs from the audit constant")

    reference_path = repo / reference_relpath
    reference_ok = reference_path.is_file() and git_tracked(repo, reference_relpath)
    checks["tracked_release_reference"] = "PASS" if reference_ok else "FAIL"
    if not reference_path.is_file():
        add_failure(result, "MISSING_TRACKED_RELEASE_REFERENCE", f"missing {reference_relpath}")
    elif not git_tracked(repo, reference_relpath):
        add_failure(result, "UNTRACKED_RELEASE_REFERENCE", f"{reference_relpath} is not in HEAD")

    reference: dict[str, Any] | None = None
    release: Path | None = None
    if reference_ok:
        try:
            reference = json.loads(reference_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            add_failure(result, "INVALID_RELEASE_REFERENCE", str(exc))
        if isinstance(reference, dict):
            release_value = reference.get("external_release_path")
            if isinstance(release_value, str):
                release = Path(release_value)
            else:
                add_failure(result, "MISSING_EXTERNAL_RELEASE_PATH", "reference lacks string external_release_path")

    release_ok = False
    if release is not None:
        parent_abs = Path(os.path.abspath(release_parent))
        release_abs = Path(os.path.abspath(release))
        path_ok = release.is_absolute() and release_abs.parent == parent_abs and release_abs.name != ""
        checks["absolute_external_release_path"] = "PASS" if path_ok else "FAIL"
        if not path_ok:
            add_failure(result, "INVALID_EXTERNAL_RELEASE_PATH", f"release must be a direct child of {parent_abs}")
        elif release.is_symlink() or not release.is_dir():
            checks["external_release_directory"] = "FAIL"
            add_failure(result, "MISSING_EXTERNAL_RELEASE", str(release))
        else:
            checks["external_release_directory"] = "PASS"
            missing = [name for name in REQUIRED_RELEASE_FILES if not (release / name).is_file()]
            symlinks = [path.relative_to(release).as_posix() for path in release.rglob("*") if path.is_symlink()]
            if missing:
                add_failure(result, "MISSING_REQUIRED_RELEASE_FILES", ", ".join(missing))
            if symlinks:
                add_failure(result, "RELEASE_SYMLINKS_FORBIDDEN", ", ".join(symlinks))
            if not missing and not symlinks:
                inventory_ok, inventory = verify_inventory(release)
                result["release_inventory"] = inventory
                checks["release_complete_and_inventory_digests"] = "PASS" if inventory_ok else "FAIL"
                if not inventory_ok:
                    add_failure(result, "RELEASE_INVENTORY_INVALID", "SHA256SUMS coverage or digest mismatch")
                release_ok = inventory_ok
    else:
        checks["absolute_external_release_path"] = "NOT_EVALUATED_NO_TRACKED_REFERENCE"
        checks["external_release_directory"] = "NOT_EVALUATED_NO_TRACKED_REFERENCE"

    # Deep semantic checks may only run after a tracked reference and complete digest-valid release.
    for gate in DEEP_GATES:
        checks.setdefault(gate, "NOT_EVALUATED_NO_DIGEST_VALID_RELEASE")
    if release_ok:
        # A release reaching here still requires a schema-specific validator from its producer.
        # This independent certifier never infers missing historical or native-output semantics.
        for gate in DEEP_GATES[1:]:
            checks[gate] = "FAIL_NO_PRODUCER_STRICT_VALIDATOR_CONTRACT"
        add_failure(
            result,
            "STRICT_SEMANTIC_CONTRACT_UNAVAILABLE",
            "digest-valid bytes alone cannot certify process/native-output/coordinate semantics",
        )

    root_finish = {name: sha256(repo / name) for name in root_relpaths}
    root_head = {name: git_file_sha256(repo, name) for name in root_relpaths}
    root_ok = root_start == root_finish == root_head
    checks["root_inputs_unchanged_during_independent_audit"] = "PASS" if root_ok else "FAIL"
    result["root_inputs"] = {
        "start_sha256": root_start,
        "finish_sha256": root_finish,
        "head_sha256": root_head,
    }
    if not root_ok:
        add_failure(result, "ROOT_INPUT_IMMUTABILITY_FAILURE", "start, finish, and HEAD content hashes differ")

    all_required_pass = reference_ok and release_ok and canonical_ok and root_ok and all(
        checks.get(gate) == "PASS" for gate in DEEP_GATES
    )
    if all_required_pass and not result["failures"]:
        result["verdict"] = "CERTIFIED_GO"
        result["qualified_for_scaling"] = True
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--reference", default=DEFAULT_REFERENCE)
    parser.add_argument("--release-parent", type=Path, default=DEFAULT_RELEASE_PARENT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    result = certify(args.repo_root, args.reference, args.release_parent)
    payload = json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    sys.stdout.write(payload)
    return 0 if result["verdict"] == "CERTIFIED_GO" else 2


if __name__ == "__main__":
    raise SystemExit(main())
