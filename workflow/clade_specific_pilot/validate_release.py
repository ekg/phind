#!/usr/bin/env python3
"""Validation script for clade-specific prophage pan-genome pilot release."""

import argparse
import hashlib
import json
import sys
from pathlib import Path


SCHEMA = "clade-specific-prophage-pilot-v1"
PASS_OR_NA = {"PASS", "NOT_APPLICABLE", "NOT_APPLICABLE_NON_SCALE_BEARING"}

# Predecessor releases
HOST_STRUCTURE_RELEASE_ID = "host-structure-1000-v1-3e16e725f70d0fdd"
PROPHAGE_SEMANTICS_RELEASE_ID = "prophage-semantics-v2-7dc695b85e5fd229"
PILOT_MIN_ASSEMBLIES_PER_CLADE = 5


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
        "root_input_sha256", "integrated_plan_sha256", "host_structure_consumer_gate",
        "prophage_semantics_consumer_gate", "clade_selection", "prophage_extraction",
        "syng_build_integrity", "syng_validation", "ancestral_inference",
        "pairwise_similarity", "phylogeny_construction", "phage_blind_validation",
        "deterministic_rerun", "injected_kill_restart", "resource", "atomic_promotion",
        "global_distinct_assembly_cap",
    }
    gates = release.get("gates", {})
    if any(gates.get(name) != "PASS" for name in required_pass_gates):
        bad = {name: gates.get(name) for name in required_pass_gates if gates.get(name) != "PASS"}
        raise RuntimeError(f"applicable release gates not PASS: {bad}")
    
    counts = release.get("counts", {})
    if counts.get("pilot_clades") < 2:
        raise RuntimeError(f"expected at least 2 pilot clades, got {counts.get('pilot_clades')}")
    if counts.get("assemblies_per_clade_min", 0) < PILOT_MIN_ASSEMBLIES_PER_CLADE:
        raise RuntimeError(f"minimum assemblies per clade not met: {counts.get('assemblies_per_clade_min')}")
    if counts.get("global_cap") != 1000:
        raise RuntimeError("global cap must be 1000")
    if counts.get("new_assembly_downloads") != 0:
        raise RuntimeError("no new assembly downloads allowed")
    
    # Verify predecessor releases
    host_structure_id = release.get("host_structure_release_id")
    prophage_semantics_id = release.get("prophage_semantics_release_id")
    
    if host_structure_id != HOST_STRUCTURE_RELEASE_ID:
        raise RuntimeError(f"host structure release ID mismatch: {host_structure_id}")
    if prophage_semantics_id != PROPHAGE_SEMANTICS_RELEASE_ID:
        raise RuntimeError(f"prophage semantics release ID mismatch: {prophage_semantics_id}")
    
    # Verify SYNG parameters
    syng_params = release.get("syng_parameters", {})
    if syng_params.get("k") != 24:
        raise RuntimeError(f"SYNG k parameter must be 24, got {syng_params.get('k')}")
    if syng_params.get("w") != 8:
        raise RuntimeError(f"SYNG w parameter must be 8, got {syng_params.get('w')}")
    
    return {
        "schema": "clade-pilot-validation-v1", "verdict": "PASS",
        "release_id": path.name, "release_verdict": release["verdict"],
        "inventory_rows": inventory_rows, "pilot_clades": counts["pilot_clades"],
        "total_assemblies": counts["total_assemblies"],
        "total_prophage_loci": counts["total_prophage_loci"],
        "syng_k": syng_params["k"], "syng_w": syng_params["w"],
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