#!/usr/bin/env python3
"""Validation script for integrated pilot 100-genome release."""

import argparse
import hashlib
import json
import sys
from pathlib import Path


SCHEMA = "integrated-pilot-100-release-v1"
PASS_OR_NA = {"PASS", "NOT_APPLICABLE"}


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_sum_file(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line:
            continue
        try:
            digest, rel = line.split("  ", 1)
        except ValueError as exc:
            raise RuntimeError(f"invalid SHA256SUMS row {path}:{line_number}") from exc
        if len(digest) != 64 or Path(rel).is_absolute() or ".." in Path(rel).parts:
            raise RuntimeError(f"unsafe SHA256SUMS row {path}:{line_number}")
        rows.append((digest, rel))
    return rows


def verify_inventory(root: Path, sums: Path) -> int:
    if not sums.is_file():
        raise RuntimeError(f"missing checksum inventory: {sums}")
    rows = parse_sum_file(sums)
    for digest, rel in rows:
        path = root / rel
        if not path.is_file() or sha_file(path) != digest:
            raise RuntimeError(f"inventory mismatch: {path}")
    return len(rows)


def validate_release(path: Path) -> dict:
    if not (path / "COMPLETE").is_file():
        raise RuntimeError("release is incomplete (COMPLETE absent)")
    
    inventory_rows = verify_inventory(path, path / "SHA256SUMS")
    release = json.loads((path / "release.json").read_text())
    
    if release.get("release_id") != path.name:
        raise RuntimeError("release identity mismatch")
    if release.get("verdict") != "PASS":
        raise RuntimeError("release verdict is not PASS")
    if release.get("consumer_action") != "ALLOW":
        raise RuntimeError("release consumer action is not ALLOW")
    
    required_pass_gates = {
        "root_input_sha256", "integrated_plan_sha256", "canonical_cohort_identity",
        "prophage_semantics_consumer_gate", "assembly_qc_reconciliation",
        "host_sketches_engineering", "syng_prefix_integrity",
        "prophage_joins_lossless", "extraction_controls",
        "impg_query_correctness", "clustering_preliminary", "matrix_states",
        "phage_blind_host", "deterministic_rerun", "injected_kill_restart",
        "resource", "atomic_promotion", "global_distinct_assembly_cap",
    }
    gates = release.get("gates", {})
    if any(gates.get(name) != "PASS" for name in required_pass_gates):
        bad = {name: gates.get(name) for name in required_pass_gates if gates.get(name) != "PASS"}
        raise RuntimeError(f"applicable release gates not PASS: {bad}")
    
    counts = release.get("counts", {})
    if counts.get("assemblies") != 100:
        raise RuntimeError(f"expected 100 assemblies, got {counts.get('assemblies')}")
    if counts.get("global_distinct_exact_assembly_revisions") != 100:
        raise RuntimeError("global distinct assembly revisions must be 100")
    if counts.get("global_cap") != 1000:
        raise RuntimeError("global cap must be 1000")
    if counts.get("new_assembly_downloads") != 0:
        raise RuntimeError("no new assembly downloads allowed")
    
    # Verify predecessor releases
    canonical_cohort_id = release.get("canonical_cohort_release_id")
    prophage_semantics_id = release.get("prophage_semantics_release_id")
    
    if canonical_cohort_id != "canonical-cohort-100-v1-6be4c0dde65f31d0":
        raise RuntimeError("canonical cohort release ID mismatch")
    if prophage_semantics_id != "prophage-semantics-v2-7dc695b85e5fd229":
        raise RuntimeError("prophage semantics release ID mismatch")
    
    return {
        "schema": "integrated-pilot-validation-v1", "verdict": "PASS",
        "release_id": path.name, "release_verdict": release["verdict"],
        "inventory_rows": inventory_rows, "assembly_count": counts["assemblies"],
        "global_distinct_exact_assembly_revisions": counts["global_distinct_exact_assembly_revisions"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release", type=Path)
    args = parser.parse_args()
    
    try:
        result = validate_release(args.release)
        print(json.dumps(result, sort_keys=True, separators=(",", ":")))
        return 0
    except RuntimeError as exc:
        print(f"VALIDATION_FAILED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())