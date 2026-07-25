#!/usr/bin/env python3
"""Failure-first unit tests for the independent canonical contract auditor."""
from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
import urllib.parse
from pathlib import Path

import audit


class AuditContractTests(unittest.TestCase):
    def test_inventory_rejects_digest_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "payload").write_text("original\n")
            digest = audit.sha_file(root / "payload")
            (root / "SHA256SUMS").write_text(f"{digest}  payload\n")
            (root / "COMPLETE").write_text(f"{audit.sha_file(root / 'SHA256SUMS')}  SHA256SUMS\n")
            (root / "payload").write_text("mutated\n")
            with self.assertRaisesRegex(audit.AuditError, "digest mismatch"):
                audit.parse_inventory(root)

    def test_inventory_rejects_uninventoried_extra_file(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "payload").write_text("ok\n")
            (root / "extra").write_text("not inventoried\n")
            digest = audit.sha_file(root / "payload")
            (root / "SHA256SUMS").write_text(f"{digest}  payload\n")
            (root / "COMPLETE").write_text(f"{audit.sha_file(root / 'SHA256SUMS')}  SHA256SUMS\n")
            with self.assertRaisesRegex(audit.AuditError, "coverage mismatch"):
                audit.parse_inventory(root)

    def test_row_hash_is_ordered_and_fail_closed(self) -> None:
        fields = ["cohort_order", "accession", "row_sha256"]
        row = {"cohort_order": "1", "accession": "GCF_000005845.2"}
        expected = hashlib.sha256(b"1\tGCF_000005845.2\n").hexdigest()
        self.assertEqual(audit.row_hash(row, fields), expected)
        row["cohort_order"] = "2"
        self.assertNotEqual(audit.row_hash(row, fields), expected)

    def test_recursive_nesting_ignores_only_rung_and_row_digest(self) -> None:
        small = [{"cohort_order": "1", "rung_n": "100", "exact_assembly_accession_version": "GCF_000005845.2", "row_sha256": "a"}]
        large = [{"cohort_order": "1", "rung_n": "250", "exact_assembly_accession_version": "GCF_000005845.2", "row_sha256": "b"}]
        self.assertTrue(audit.nested_rows(large, small))
        large[0]["exact_assembly_accession_version"] = "GCF_000005845.3"
        self.assertFalse(audit.nested_rows(large, small))

    def test_pansn_percent_encoding_is_reversible(self) -> None:
        encoded, policy = audit.encode_contig(b"contig:name#1", "GCF_000005845.2")
        self.assertEqual(encoded, "contig%3Aname%231")
        self.assertEqual(policy, "PERCENT_UTF8_BYTES_V1")
        self.assertEqual(urllib.parse.unquote_to_bytes(encoded), b"contig:name#1")

    def test_machine_verdict_cannot_be_misread_as_integrated_go(self) -> None:
        result = {
            "applicable_gates": {"identity": "PASS"},
            "not_applicable_gates": {
                "prophage_extraction_source_coordinate_policy": "NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED"
            },
        }
        verdict = audit.build_verdict(result)
        self.assertEqual(verdict["verdict"], "CANONICAL_GO_500")
        self.assertEqual(verdict["scope"], "CANONICAL_ACQUISITION_AND_CANONICALIZATION_ONLY")
        self.assertFalse(verdict["extraction_branch"]["integrated_go_claimed"])
        self.assertEqual(verdict["extraction_branch"]["verdict"], "EXTRACTION_BLOCKED")


if __name__ == "__main__":
    unittest.main(verbosity=2)
