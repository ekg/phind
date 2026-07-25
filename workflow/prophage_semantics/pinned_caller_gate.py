#!/usr/bin/env python3
"""Independently verify the immutable pinned-Phigaro comparison release and
derive the historical-CSV extraction verdict.

This gate does NOT trust the predecessor's own comparison code or its verdict
strings.  It independently verifies the predecessor external release integrity,
the exact N=10 input identity/order, the version-pinned tools/database/config,
the two SEPARATE machine verdicts, and then RE-DERIVES the historical
attribution evidence with its own code
(``independent_rerun_verification``).  Only if the predecessor declares
``historical_csv_attribution == DECISIVE`` AND the independent re-derivation is
sound is the historical CSV extraction changed to ``EXTRACTION_GO``.

``modern_v2_4_pilot`` is carried through as a strictly separate verdict: a modern
GO alone can never authorize historical CSV extraction (see the task contract).

This is a metadata/checksum/bounded-known-base gate only.  It reads the
predecessor's released native prophage TSV/saved-FASTA (tiny, already extracted
from the 10 validated assemblies) and the immutable CSV header+rows.  It never
reads genome FASTA, never downloads anything, and never exceeds the frozen N=10
pilot cohort.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.prophage_semantics.release import (
    ACCESSIONS_SHA,
    COHORT_ORDER,
    COHORT_RELEASE_ID,
    COHORT_RELEASE_JSON_SHA,
    GateError,
    SOURCE_SHA,
    atomic_write,
    canonical_bytes,
    parse_sum_file,
    sha_bytes,
    sha_file,
    utc_now,
    verify_exact,
)
from workflow.prophage_semantics.independent_rerun_verification import verify as independently_verify

SCHEMA = "pinned-phigaro-consumption-gate-v1"
RELEASE_SCHEMA = "phigaro-version-comparison-release-v1"
TOOLS_SCHEMA = "phigaro-version-comparison-tools-v1"
DEFAULT_EXTERNAL_ROOT = Path(
    "/home/erikg/phind-data/ecoli26k/v1/releases/rerun-phigaro-version-comparison"
)
# Immutable predecessor identity (the external release this task consumes).
PREDECESSOR_RELEASE_ID = "phigaro-version-comparison-v1-e7cfa43b9231aee5"
# Pin the producer versions independently of the predecessor's own declaration.
PHIGARO_V240_COMMIT = "1ff5f85cee31e418bce24e4cd51c7528c43bc968"
PHIGARO_V240_SDIST_SHA256 = "fd764d792a37984bcabaea0da39185dc6c864b8ecfbd8a806553a68ac876d800"
PHIGARO_V230_COMMIT = "aea9469d09cdbfbb528998ebc43232ee9f44decd"


def _hex(value: Any, length: int, label: str) -> str:
    if not isinstance(value, str) or len(value) != length:
        raise GateError(f"{label} is not a {length}-character digest/commit")
    try:
        int(value, 16)
    except ValueError as exc:
        raise GateError(f"{label} is not hexadecimal") from exc
    return value


def _load_json(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise GateError(f"required pinned-caller input missing: {label}: {path}")
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise GateError(f"invalid {label}: {path}") from exc
    if not isinstance(value, dict):
        raise GateError(f"invalid {label}: expected JSON object")
    return value


def _within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def verify_release_inventory(external: Path) -> list[tuple[str, str]]:
    sums = external / "SHA256SUMS"
    if not sums.is_file():
        raise GateError("pinned-caller release lacks SHA256SUMS")
    rows = parse_sum_file(sums)
    by_rel = {rel: digest for digest, rel in rows}
    for rel, digest in by_rel.items():
        path = external / rel
        if not path.is_file():
            raise GateError(f"pinned-caller inventory file missing: {rel}")
        if sha_file(path) != digest:
            raise GateError(f"pinned-caller inventory checksum mismatch: {rel}")
    # every non-excluded file in the tree must be inventoried
    excluded = {"COMPLETE", "SHA256SUMS", "release.json"}
    extra = sorted(
        p.relative_to(external).as_posix()
        for p in external.rglob("*")
        if p.is_file() and p.relative_to(external).as_posix() not in excluded
        and p.relative_to(external).as_posix() not in by_rel
    )
    if extra:
        raise GateError("pinned-caller release has uninventoried files: " + ", ".join(extra[:8]))
    return rows


def decide_extraction(
    historical_csv_attribution: str,
    decisive_evidence_independently_sound: bool,
    modern_v2_4_pilot: str,
) -> tuple[str, str]:
    """Pure decision rule. Historical extraction is GO only when the historical
    attribution is DECISIVE AND independently re-verified as sound. The modern
    v2.4 pilot verdict is NEVER allowed to authorize historical extraction.
    """
    if historical_csv_attribution == "DECISIVE" and decisive_evidence_independently_sound:
        return "EXTRACTION_GO", "ALLOW"
    return "EXTRACTION_BLOCKED", "REJECT"


def validate_pinned_caller_release(repo: Path, external_root: Path) -> dict[str, Any]:
    repo = repo.resolve()
    external_root = external_root.resolve()

    # Recheck global/root immutability before touching a predecessor object.
    verify_exact(repo / "26k_ecoli_accession.txt", ACCESSIONS_SHA, "root accession input")
    verify_exact(repo / "26k_prophage1.csv", SOURCE_SHA, "root prophage input")
    verify_exact(
        repo / "manifests/canonical-cohort-010-v1/release.json",
        COHORT_RELEASE_JSON_SHA,
        "canonical cohort manifest",
    )

    external = external_root / PREDECESSOR_RELEASE_ID
    if not external.is_dir():
        raise GateError("pinned-caller external release directory is absent")
    if not _within(external, external_root):
        raise GateError("pinned-caller external path/namespace mismatch")
    complete = external / "COMPLETE"
    if not complete.is_file():
        raise GateError("pinned-caller release lacks COMPLETE")

    # integrity: every declared file matches, no uninventoried files
    inventory_rows = verify_release_inventory(external)

    release = _load_json(external / "release.json", "pinned-caller release manifest")
    if release.get("schema") != RELEASE_SCHEMA:
        raise GateError("pinned-caller release schema mismatch")
    if release.get("release_id") != PREDECESSOR_RELEASE_ID:
        raise GateError("pinned-caller release identity mismatch")
    if release.get("source_task_id") != "rerun-pinned-phigaro":
        raise GateError("pinned-caller source task mismatch")
    if release.get("predecessor_release_id") != COHORT_RELEASE_ID:
        raise GateError("pinned-caller canonical-cohort predecessor mismatch")
    if release.get("immutable") is not True:
        raise GateError("pinned-caller release is not marked immutable")

    # root input gates (predecessor's own immutable-input checks)
    root_gates = release.get("root_input_gates", {})
    for name, expected in (
        ("26k_ecoli_accession.txt", ACCESSIONS_SHA),
        ("26k_prophage1.csv", SOURCE_SHA),
        ("canonical-cohort-010-v1 assemblies.tsv", "7133058093e3f08c132248b3cf4453c7076b6550e838f2ab9f39a6b5b7b8fcbd"),
    ):
        g = root_gates.get(name, {})
        if g.get("pass") is not True or g.get("required_sha256") != expected or g.get("actual_sha256") != expected:
            raise GateError(f"pinned-caller root input gate failed: {name}")

    # exact N=10 input identity/order + 56-row cardinality
    counts = release.get("counts", {})
    if counts.get("assemblies") != 10 or counts.get("csv_rows_for_n10") != 56:
        raise GateError("pinned-caller N=10/row cardinality mismatch")
    if counts.get("global_distinct_assembly_cap") != 1000 or counts.get("new_assembly_downloads") != 0:
        raise GateError("pinned-caller global cap/download contract mismatch")
    per_asm = release.get("comparison", {}).get("csv_per_assembly", {})
    # The predecessor's per-assembly dict is keyed alphanumerically; the frozen
    # canonical acquisition ORDER is independently enforced by the cohort
    # manifest above.  Here we assert exact set identity + 56-row accounting.
    if set(per_asm) != set(COHORT_ORDER) or sum(per_asm.values()) != 56:
        raise GateError("pinned-caller cohort identity/row accounting mismatch")

    # fixture gate (both versions must pass the official fixture)
    fixture = release.get("fixture_gate", {})
    for ver in ("v2.3.0", "v2.4.0"):
        fg = fixture.get(ver, {})
        if fg.get("pass") is not True or fg.get("prophage_count") != fg.get("expected"):
            raise GateError(f"pinned-caller official fixture gate failed: {ver}")

    # tool/database/config pins (independently checked against required values)
    tools = _load_json(external / "tools.json", "pinned-caller tools manifest")
    if tools.get("schema") != TOOLS_SCHEMA:
        raise GateError("pinned-caller tools schema mismatch")
    pins_v24 = tools.get("phigaro", {}).get("v2.4.0", {})
    if pins_v24.get("version") != "2.4.0" or pins_v24.get("git_tag_commit") != PHIGARO_V240_COMMIT:
        raise GateError("Phigaro v2.4.0 version/commit pin mismatch")
    if pins_v24.get("pypi_sdist_sha256") != PHIGARO_V240_SDIST_SHA256:
        raise GateError("Phigaro v2.4.0 PyPI sdist pin mismatch")
    pins_v23 = tools.get("phigaro", {}).get("v2.3.0", {})
    if pins_v23.get("version") != "2.3.0" or pins_v23.get("git_tag_commit") != PHIGARO_V230_COMMIT:
        raise GateError("Phigaro v2.3.0 version/commit pin mismatch")
    _hex(pins_v23.get("pypi_sdist_sha256"), 64, "Phigaro v2.3.0 PyPI sdist SHA-256")
    # native TSV coordinate conventions are documented and code-backed
    if "1-based inclusive" not in pins_v23.get("tsv_coordinate_convention", ""):
        raise GateError("v2.3.0 TSV coordinate convention not documented as 1-based inclusive")
    if "0-based" not in pins_v24.get("tsv_coordinate_convention", ""):
        raise GateError("v2.4.0 TSV coordinate convention not documented as 0-based")
    # Prodigal/HMMER/pVOG must all be pinned
    if tools.get("prodigal", {}).get("version") != "2.6.3":
        raise GateError("Prodigal version pin mismatch")
    _hex(tools.get("prodigal", {}).get("binary_sha256"), 64, "Prodigal binary SHA-256")
    if tools.get("hmmer", {}).get("version") != "3.3.2":
        raise GateError("HMMER version pin mismatch")
    _hex(tools.get("hmmer", {}).get("hmmsearch_binary_sha256"), 64, "HMMER hmmsearch binary SHA-256")
    pvogs = tools.get("pvogs_database", {}).get("files", {})
    if "allpvoghmms" not in pvogs or _hex(pvogs["allpvoghmms"], 64, "pVOG allpvoghmms SHA-256") is None:
        raise GateError("pVOG database main HMM pin missing")

    # engineering gates (fail-closed on the integrity-bearing booleans)
    eng = release.get("engineering_gates", {})
    required_eng = {
        "canonical_input_digests", "fixture_gate_both_versions", "global_union_le_1000",
        "immutable_root_inputs", "no_new_assembly_download", "out_of_range_intervals_zero",
        "strand_never_inferred", "topology_unknown",
    }
    bad = {k: eng.get(k) for k in required_eng if eng.get(k) is not True}
    if bad:
        raise GateError("pinned-caller engineering gate not satisfied: " + json.dumps(bad, sort_keys=True))

    # the two SEPARATE machine verdicts
    verdicts = release.get("machine_verdicts", {})
    historical = verdicts.get("historical_csv_attribution", {})
    modern = verdicts.get("modern_v2_4_pilot", {})
    if historical.get("verdict") not in {"DECISIVE", "NON_DECISIVE"}:
        raise GateError("historical_csv_attribution machine verdict missing/invalid")
    if modern.get("verdict") not in {"GO", "NO_GO"}:
        raise GateError("modern_v2_4_pilot machine verdict missing/invalid")
    if historical.get("version_identified") != "2.3.0" or historical.get("caller_hypothesis") != "Phigaro":
        raise GateError("historical attribution caller/version claim mismatch")

    # INDEPENDENT re-derivation of the historical attribution evidence.
    # This re-reads the checksum-validated native outputs and the immutable CSV
    # and re-derives whether all 56 rows are reproduced exactly, the version
    # signature, and the coordinate convention — without trusting the
    # predecessor's comparison_pairs.json or its DECISIVE string.
    independent = independently_verify(repo, external, list(COHORT_ORDER))
    if not independent["decisive_evidence_independently_sound"]:
        raise GateError(
            "independent re-derivation of historical attribution is not sound: "
            + json.dumps(
                {
                    "all_fields_exact_count": independent["all_fields_exact_count"],
                    "csv_rows_for_cohort": independent["csv_rows_for_cohort"],
                    "csv_unmatched_count": independent["csv_unmatched_count"],
                    "native_surplus_v23_count": independent["native_surplus_v23_count"],
                    "boundary_signature": independent["boundary_signature"],
                },
                sort_keys=True,
            )
        )

    # Decision: historical CSV extraction is GO only if the predecessor verdict
    # is DECISIVE AND the independent re-derivation is sound. modern_v2_4_pilot
    # is kept strictly separate and cannot authorize historical extraction.
    extraction, consumer_action = decide_extraction(historical["verdict"], True, modern["verdict"])

    return {
        "schema": SCHEMA,
        "captured_at_utc": utc_now(),
        "verdict": "PASS",
        "release_id": PREDECESSOR_RELEASE_ID,
        "external_path": str(external),
        "complete_sha256": sha_file(complete),
        "sha256sums_sha256": sha_file(external / "SHA256SUMS"),
        "release_json_sha256": sha_file(external / "release.json"),
        "comparison_pairs_sha256": sha_file(external / "comparison_pairs.json"),
        "tools_json_sha256": sha_file(external / "tools.json"),
        "inventory_rows": len(inventory_rows),
        "cohort_order": list(COHORT_ORDER),
        "distinct_assemblies": 10,
        "historical_rows": 56,
        "new_assembly_downloads": 0,
        "global_distinct_exact_assembly_revisions": 10,
        "global_cap": 1000,
        "historical_csv_attribution": historical["verdict"],
        "modern_v2_4_pilot": modern["verdict"],
        "modern_v2_4_pilot_separate": True,
        "decisive_evidence_independently_sound": True,
        "independent_reverification": {
            "all_fields_exact_count": independent["all_fields_exact_count"],
            "csv_rows_for_cohort": independent["csv_rows_for_cohort"],
            "exact_coordinate_count": independent["exact_coordinate_count"],
            "boundary_signature": independent["boundary_signature"],
            "coordinate_convention_from_saved_fasta": independent["coordinate_convention_from_saved_fasta"],
        },
        "historical_csv_extraction": extraction,
        "consumer_action": consumer_action,
        "root_input_observed_sha256": {
            "26k_ecoli_accession.txt": ACCESSIONS_SHA,
            "26k_prophage1.csv": SOURCE_SHA,
        },
    }


def failure_result(repo: Path, external_root: Path, exc: Exception) -> dict[str, Any]:
    def digest_or_missing(path: Path) -> str | None:
        return sha_file(path) if path.is_file() else None

    candidates = []
    if external_root.is_dir():
        candidates = sorted(p.name for p in external_root.iterdir() if p.is_dir() and (p / "COMPLETE").is_file())
    return {
        "schema": SCHEMA,
        "captured_at_utc": utc_now(),
        "verdict": "NO_GO",
        "hard_stop": True,
        "reason": str(exc),
        "required_predecessor_task": "rerun-pinned-phigaro",
        "required_release_id": PREDECESSOR_RELEASE_ID,
        "external_root": str(external_root),
        "external_root_exists": external_root.is_dir(),
        "complete_release_candidates": candidates,
        "historical_csv_attribution": "UNAVAILABLE",
        "modern_v2_4_pilot": "UNAVAILABLE",
        "decisive_evidence_independently_sound": False,
        "historical_csv_extraction": "EXTRACTION_BLOCKED",
        "consumer_action": "REJECT",
        "root_input_observed_sha256": {
            "26k_ecoli_accession.txt": digest_or_missing(repo / "26k_ecoli_accession.txt"),
            "26k_prophage1.csv": digest_or_missing(repo / "26k_prophage1.csv"),
        },
        "new_assembly_downloads": 0,
        "sequence_bases_read": 0,
        "global_distinct_exact_assembly_revisions": 10,
        "global_cap": 1000,
    }


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", default=".")
    p.add_argument("--external-root", default=str(DEFAULT_EXTERNAL_ROOT))
    p.add_argument("--output")
    return p


def main() -> int:
    args = parser().parse_args()
    repo = Path(args.repo).resolve()
    external_root = Path(args.external_root).resolve()
    try:
        result = validate_pinned_caller_release(repo, external_root)
        rc = 0
    except (GateError, OSError, ValueError, KeyError) as exc:
        result = failure_result(repo, external_root, exc)
        rc = 2
    if args.output:
        output = Path(args.output)
        if not output.is_absolute():
            output = repo / output
        atomic_write(output, canonical_bytes(result))
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
