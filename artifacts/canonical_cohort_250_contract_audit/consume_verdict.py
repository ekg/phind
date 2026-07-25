#!/usr/bin/env python3
"""Fail-closed consumer for the checksum-pinned CANONICAL_GO_500 verdict."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


class ConsumerError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ConsumerError(message)


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def consume(repo: Path) -> dict[str, Any]:
    audit_root = repo.resolve() / "artifacts/canonical_cohort_250_contract_audit"
    sums_path = audit_root / "SHA256SUMS"
    require(sums_path.is_file(), "missing audit SHA256SUMS")
    inventory: dict[str, str] = {}
    for number, line in enumerate(sums_path.read_text().splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  ([^/]+)", line)
        require(match is not None, f"malformed audit inventory line {number}")
        digest, name = match.groups()
        require(name not in inventory, f"duplicate audit inventory path: {name}")
        path = audit_root / name
        require(path.is_file() and not path.is_symlink() and sha_file(path) == digest,
                f"audit inventory mismatch: {name}")
        inventory[name] = digest
    actual = {path.name for path in audit_root.iterdir() if path.is_file() and path.name != "SHA256SUMS"}
    require(actual == set(inventory), "audit inventory coverage mismatch")
    require("audit.json" in inventory and "verdict.json" in inventory, "missing machine audit/verdict")
    audit = json.loads((audit_root / "audit.json").read_text())
    verdict = json.loads((audit_root / "verdict.json").read_text())
    require(verdict.get("verdict") == "CANONICAL_GO_500", "canonical GO verdict absent")
    require(verdict.get("scope") == "CANONICAL_ACQUISITION_AND_CANONICALIZATION_ONLY",
            "verdict scope is not canonical-only")
    require(verdict.get("consumer_task_id") == "prepare-canonical-cohort-500", "wrong verdict consumer")
    require(verdict.get("audit_result_sha256") == inventory["audit.json"], "audit result pin mismatch")
    require(audit.get("audit_verdict") == "PASS" and audit.get("release_immutable") is True,
            "audit did not pass immutable release")
    require(audit.get("applicable_gates") == verdict.get("applicable_gates")
            and all(value == "PASS" for value in verdict["applicable_gates"].values()),
            "applicable canonical gate mismatch")
    require(all(str(value).startswith("NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY")
                for value in verdict.get("not_applicable_gates", {}).values()),
            "non-canonical gates were not explicitly NOT_APPLICABLE")
    extraction = verdict.get("extraction_branch", {})
    require(extraction.get("task_status") == "failed" and extraction.get("verdict") == "EXTRACTION_BLOCKED"
            and extraction.get("consumer_action") == "REJECT" and extraction.get("integrated_go_claimed") is False,
            "extraction branch separation mismatch")
    release_path = repo / "manifests/canonical-cohort-250-v1/release.json"
    release = json.loads(release_path.read_text())
    external = Path(release["external_release_path"])
    require(release.get("release_id") == verdict.get("source_release_id"), "source release ID drift")
    require(sha_file(release_path) == verdict.get("source_release_json_sha256")
            and sha_file(external / "release.json") == verdict.get("source_release_json_sha256"),
            "source release JSON pin drift")
    require(sha_file(external / "SHA256SUMS") == verdict.get("source_external_sha256sums_sha256"),
            "source external inventory pin drift")
    require(sha_file(external / "COMPLETE") == verdict.get("source_complete_sha256"),
            "source COMPLETE pin drift")
    roots = audit["snapshots"]["root_inputs"]
    require(all(sha_file(repo / name) == digest for name, digest in roots.items()), "root input pin drift")
    return {
        "schema": "canonical-go-500-consumption-v1", "verdict": "PASS",
        "consumed_verdict": "CANONICAL_GO_500", "release_id": release["release_id"],
        "scope": verdict["scope"], "audit_json_sha256": inventory["audit.json"],
        "verdict_json_sha256": inventory["verdict.json"],
        "extraction_eligibility": "EXTRACTION_BLOCKED",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    args = parser.parse_args()
    try:
        print(canonical_json(consume(Path(args.repo_root))).decode(), end="")
        return 0
    except (ConsumerError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"CANONICAL_NO_GO_500: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
