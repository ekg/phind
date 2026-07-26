#!/usr/bin/env python3
"""Independent validation script for the integrated 250-genome scale-bearing release.

This mirrors the in-workflow ``validate_release`` but is standalone so a consumer
(rung N=500) can verify the release by digest without importing the builder.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

SCHEMA = "integrated-pilot-250-release-v1"

# Immutable predecessor/contract digests (must match the builder exactly)
CANONICAL_COHORT_RELEASE_ID = "canonical-cohort-250-v1-a6184d7d6ee08bda"
CANONICAL_COHORT_RELEASE_JSON_SHA = "dcf2b887afa51e4e0e739ae2fef9b5a9d72fb8bc9a4d698a161a99673aaf504a"
PRIOR_INTEGRATED_RELEASE_ID = "integrated-pilot-100-v1-0a11eda244a9def8"
PRIOR_INTEGRATED_RELEASE_JSON_SHA = "6816c4e24f6511e45196d91112da96ab7f56082732c4c363b74ae7010a80e273"
PROPHAGE_SEMANTICS_RELEASE_ID = "prophage-semantics-v2-7dc695b85e5fd229"
PROPHAGE_SEMANTICS_RELEASE_JSON_SHA = "5d8403eb070d8a62140adfe7260b7fde6897598f72ac1c536879e78e8ea2b992"
ACCESSIONS_SHA = "1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5"
SOURCE_SHA = "6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996"
GLOBAL_CAP = 1000
RUNG = 250


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


def validate_release(path: Path, repo: Path | None = None) -> dict:
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
    if release.get("go_500") != "GO_500":
        raise RuntimeError("release go_500 is not GO_500 (scale-up not authorized)")

    required_pass_gates = {
        "root_input_sha256", "integrated_plan_sha256", "canonical_cohort_identity",
        "prior_integrated_reuse", "prophage_semantics_consumer_gate",
        "assembly_qc_reconciliation", "host_sketches_engineering", "syng_prefix_integrity",
        "prophage_joins_lossless", "extraction_controls", "impg_query_correctness",
        "clustering_preliminary", "matrix_states", "phage_blind_host",
        "scale_trend", "deterministic_rerun", "injected_kill_restart",
        "resource", "atomic_promotion", "global_distinct_assembly_cap",
    }
    gates = release.get("gates", {})
    bad = {name: gates.get(name) for name in required_pass_gates if gates.get(name) != "PASS"}
    if bad:
        raise RuntimeError(f"applicable release gates not PASS: {bad}")

    counts = release.get("counts", {})
    if counts.get("assemblies") != RUNG:
        raise RuntimeError(f"expected {RUNG} assemblies, got {counts.get('assemblies')}")
    if counts.get("global_distinct_exact_assembly_revisions") != RUNG:
        raise RuntimeError("global distinct assembly revisions must be the rung size")
    if counts.get("global_cap") != GLOBAL_CAP:
        raise RuntimeError("global cap must be 1000")
    if counts.get("new_assembly_downloads") != 0:
        raise RuntimeError("no new assembly downloads allowed")

    if release.get("canonical_cohort_release_id") != CANONICAL_COHORT_RELEASE_ID:
        raise RuntimeError("canonical cohort release ID mismatch")
    if release.get("prior_integrated_release_id") != PRIOR_INTEGRATED_RELEASE_ID:
        raise RuntimeError("prior integrated release ID mismatch")
    if release.get("prophage_semantics_release_id") != PROPHAGE_SEMANTICS_RELEASE_ID:
        raise RuntimeError("prophage semantics release ID mismatch")

    # Scale-trend checks (hard gate on runtime gating metric; matches builder semantics)
    scale = json.loads((path / "scale_trend.json").read_text())
    if scale.get("verdict") != "PASS" or scale.get("go_500") != "GO_500":
        raise RuntimeError("scale_trend is not PASS/GO_500")
    if scale["time_exponent"]["empirical_upper_bound"] > 1.3:
        raise RuntimeError("time exponent upper bound exceeds 1.3")
    slopes = scale["last_two_rung_per_base_slopes"]
    gating = slopes.get("gating_metric", "wall_seconds_per_new_base")
    gating_change = slopes["relative_changes"][gating]
    if abs(gating_change) > 0.25:
        raise RuntimeError(f"runtime per-base slope change exceeds 25%: {gating}={gating_change}")
    # all applicable machine/projection gates must be PASS
    for name, ok in scale.get("checks", {}).items():
        if not ok:
            raise RuntimeError(f"scale_trend check not satisfied: {name}")

    # Root-input immutability (when repo is supplied)
    if repo is not None:
        for fname, expected in (("26k_ecoli_accession.txt", ACCESSIONS_SHA),
                                ("26k_prophage1.csv", SOURCE_SHA)):
            got = sha_file(repo / fname)
            if got != expected:
                raise RuntimeError(f"root input {fname} digest mismatch: {got} != {expected}")

    return {
        "schema": "integrated-pilot-250-validation-v1", "verdict": "PASS",
        "release_id": path.name, "release_verdict": release["verdict"],
        "go_500": release["go_500"], "inventory_rows": inventory_rows,
        "assemblies": counts["assemblies"],
        "time_exponent_upper_bound": scale["time_exponent"]["empirical_upper_bound"],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("release", type=Path)
    ap.add_argument("--repo", type=Path, default=None)
    args = ap.parse_args()
    try:
        result = validate_release(args.release, args.repo)
        print(json.dumps(result, sort_keys=True, separators=(",", ":")))
        return 0
    except RuntimeError as exc:
        print(f"VALIDATION_FAILED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
