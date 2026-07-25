import http.client
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from workflow.acquisition_canonicalization import pilot as p
from artifacts.canonical_cohort_100 import runner as r


class CanonicalCohort100Tests(unittest.TestCase):
    def setUp(self):
        self.old_accessions = list(p.EXPECTED_ACCESSIONS)
        self.old_task = p.TASK_ID
        self.old_agent = p.USER_AGENT

    def tearDown(self):
        p.EXPECTED_ACCESSIONS = self.old_accessions
        p.TASK_ID = self.old_task
        p.USER_AGENT = self.old_agent

    def test_pinned_real_inputs_are_exact_nested_100_and_all_required_gates_pass(self):
        inputs = r.verify_inputs(Path(".").resolve())
        self.assertEqual(len(inputs["rows"]), 100)
        self.assertEqual(len(set(inputs["accessions"])), 100)
        self.assertEqual(inputs["accessions"][:10], inputs["predecessor"]["sequence_bearing_assembly_accessions"])
        self.assertEqual(inputs["compatibility"]["counts"]["pass_consumer_gates"], 19)

    def test_manifest_row_checksum_mismatch_fails_closed(self):
        original = Path("manifests/pilot-cohorts-v1/cohort-0100.tsv").read_text()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "cohort.tsv"
            lines = original.splitlines()
            columns = lines[1].split("\t")
            columns[3] = "GCF_999999999.1"
            lines[1] = "\t".join(columns)
            path.write_text("\n".join(lines) + "\n")
            with self.assertRaisesRegex(p.GateError, "row checksum mismatch"):
                r.read_hashed_tsv(path)

    def test_tracked_inventory_checksum_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "payload").write_text("changed")
            (root / "SHA256SUMS").write_text("0" * 64 + "  payload\n")
            with self.assertRaisesRegex(p.GateError, "tracked checksum mismatch"):
                r.verify_tracked_inventory(root)

    def test_resource_blank_disk_and_inode_refusal(self):
        with self.assertRaisesRegex(p.GateError, "blank"):
            p.Allocations(0, 1, 1, 1, 1, 1, 1, 1).validate()
        with self.assertRaisesRegex(p.GateError, "durable"):
            p.Allocations(100, 100, 100, 100, 71, 1, 1, 1).validate()
        with self.assertRaisesRegex(p.GateError, "inode"):
            p.Allocations(100, 100, 100, 4, 1, 1, 3, 1).validate()

    def test_url_constructor_is_bounded_to_exact_frozen_100(self):
        inputs = r.verify_inputs(Path(".").resolve())
        r.configure_pinned_primitives(inputs["accessions"])
        self.assertIn(inputs["accessions"][99], p.download_url(inputs["accessions"][99]))
        with self.assertRaisesRegex(p.GateError, "outside immutable"):
            p.download_url("GCF_999999999.1")
        with self.assertRaisesRegex(p.GateError, "outside immutable"):
            p.download_url(inputs["accessions"][10].rsplit(".", 1)[0])

    def test_chunked_http_incomplete_read_is_bounded_and_reenters_pinned_resume(self):
        expected = (Path("source"), {"state": "COMPLETE"})
        with tempfile.TemporaryDirectory() as td, mock.patch.object(
            p, "commit_source_object", side_effect=[http.client.IncompleteRead(b"partial"), expected]
        ) as commit:
            root = Path(td)
            with mock.patch.object(r.time, "sleep"):
                observed = r.commit_source_with_transport_retry(
                    "GCF_000005845.2", root / "objects", root / "state.jsonl",
                    root / "failures.jsonl", 2, 0, False, 10,
                )
            self.assertEqual(observed, expected)
            self.assertEqual(commit.call_count, 2)
            self.assertIn("OUTER_HTTP_TRANSPORT_RETRY", (root / "failures.jsonl").read_text())

    def test_invalid_checksum_complete_new_object_is_discarded_for_recompute(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            obj = root / "object"
            obj.mkdir()
            (obj / "manifest.json").write_text(json.dumps({"accession": "GCF_000005845.2", "state": "COMPLETE"}))
            (obj / "SHA256SUMS").write_text("0" * 64 + "  manifest.json\n")
            state, failures = root / "state.jsonl", root / "failures.jsonl"
            r._discard_invalid_completed(obj, "GCF_000005845.2", "source", state, failures)
            self.assertFalse(obj.exists())
            self.assertIn("INVALID_COMPLETED_OBJECT_DISCARDED", failures.read_text())
            self.assertIn("RECOMPUTE_ALLOWED", state.read_text())

    def test_interrupted_promotion_never_exposes_partial_release(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            stage, final = root / ".stage.release", root / "release"
            stage.mkdir()
            (stage / "payload").write_text("valid")
            self.assertFalse(final.exists())
            self.assertFalse((stage / "COMPLETE").exists())
            p.seal_directory(stage, final)
            self.assertFalse(stage.exists())
            self.assertTrue((final / "COMPLETE").exists())
            p.verify_sha_inventory(final)

    def test_release_identity_is_order_and_pin_deterministic(self):
        inputs = r.verify_inputs(Path(".").resolve())
        first = r.release_id(inputs["accessions"])
        second = r.release_id(list(inputs["accessions"]))
        changed = list(inputs["accessions"])
        changed[-1], changed[-2] = changed[-2], changed[-1]
        self.assertEqual(first, second)
        self.assertNotEqual(first, r.release_id(changed))
        self.assertTrue(first.startswith("canonical-cohort-100-v1-"))


if __name__ == "__main__":
    unittest.main()
