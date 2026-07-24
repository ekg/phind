import json
import tempfile
import unittest
from pathlib import Path

from workflow.prophage_semantics import release


ROOT = Path(__file__).resolve().parents[2]


class ProphageSemanticsTests(unittest.TestCase):
    def test_manifest_checksum_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "manifest.json"
            path.write_text("wrong\n")
            with self.assertRaisesRegex(release.GateError, "checksum mismatch"):
                release.verify_exact(path, "0" * 64, "fixture manifest")

    def test_blank_and_overcommitted_resources_refused(self):
        blank = release.Allocations(0, 1, 1, 1, 1, 1, 1, 1)
        with self.assertRaisesRegex(release.GateError, "blank or non-positive"):
            blank.validate()
        disk_over = release.Allocations(100, 100, 100, 100, 71, 1, 1, 1)
        with self.assertRaisesRegex(release.GateError, "durable upper-95"):
            disk_over.validate()
        inode_over = release.Allocations(100, 100, 100, 3, 1, 1, 2, 1)
        with self.assertRaisesRegex(release.GateError, "50% inode"):
            inode_over.validate()

    def test_interrupted_unit_resume_validates_and_mismatch_refuses(self):
        with tempfile.TemporaryDirectory() as td:
            stage = Path(td)
            (stage / "state.jsonl").write_text("")
            release.write_static_unit(stage, "unit.json", b"{\"x\":1}\n")
            release.write_static_unit(stage, "unit.json", b"{\"x\":1}\n")
            events = [json.loads(x)["event"] for x in (stage / "state.jsonl").read_text().splitlines()]
            self.assertEqual(events, ["UNIT_COMMITTED", "RESUME_UNIT_VALIDATED"])
            with self.assertRaisesRegex(release.GateError, "refusing mixed resume"):
                release.write_static_unit(stage, "unit.json", b"{\"x\":2}\n")

    def test_partial_publication_without_complete_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "release"
            path.mkdir()
            (path / "release.json").write_text("{}\n")
            with self.assertRaisesRegex(release.GateError, "COMPLETE absent"):
                release.validate_release(path)

    def test_policy_is_complete_lossless_and_blocked(self):
        policy = json.loads((ROOT / "artifacts/prophage_semantics/semantics_policy_v1.json").read_text())
        release.validate_policy(policy)
        self.assertEqual(policy["extraction_gate"]["verdict"], "EXTRACTION_BLOCKED")
        self.assertEqual(policy["extraction_gate"]["consumer_action"], "REJECT")
        blocked = set(policy["extraction_gate"]["blocking_dimensions"])
        unresolved_critical = {
            d["name"] for d in policy["semantic_dimensions"]
            if d["extraction_critical"] and d["status"] != "RESOLVED"
        }
        self.assertEqual(blocked, unresolved_critical)
        self.assertFalse(any(s["generic_tagged_alias"] for s in policy["scopes"]))
        self.assertEqual({c["id"] for c in policy["coordinate_candidates"]}, {
            "C1_RAW_1_BASED_CLOSED", "C2_RAW_0_BASED_INCLUSIVE"
        })

    def test_root_profile_has_lossless_row_and_scope_accounting(self):
        release.verify_exact(ROOT / "26k_ecoli_accession.txt", release.ACCESSIONS_SHA, "accessions")
        release.verify_exact(ROOT / "26k_prophage1.csv", release.SOURCE_SHA, "prophage")
        profile, rows = release.source_profile(ROOT / "26k_prophage1.csv")
        self.assertEqual(len(rows), 132404)
        self.assertEqual(profile["scope_rows"], {
            "all_records": 132404,
            "transposable_flag_positive": 7695,
            "taxonomy_assigned": 115442,
        })
        self.assertEqual(profile["duplicate_locus_groups"], 0)
        self.assertTrue(profile["raw_rows_preserved_in_place"])
        self.assertFalse(profile["normalization_materialized"])

    def test_exact_predecessor_and_bounded_diagnostic(self):
        predecessor, external, assemblies, inventory_rows = release.validate_predecessor(ROOT)
        self.assertEqual(predecessor["release_id"], release.COHORT_RELEASE_ID)
        self.assertEqual(len(assemblies), 10)
        self.assertGreater(inventory_rows, 0)
        _, rows = release.source_profile(ROOT / "26k_prophage1.csv")
        sentinel = release.annotation_boundary_diagnostic(ROOT, external, assemblies, rows)
        self.assertEqual(sentinel["distinct_assemblies"], 10)
        self.assertEqual(sentinel["source_rows"], 56)
        self.assertEqual(sentinel["new_assembly_downloads"], 0)
        self.assertEqual(sentinel["sequence_bases_read"], 0)
        self.assertEqual(
            sentinel["boundary_matches_after_adding_delta_to_both_raw_coordinates"]["0"]["both_boundaries"],
            45,
        )
        self.assertEqual(
            sentinel["boundary_matches_after_adding_delta_to_both_raw_coordinates"]["1"]["both_boundaries"],
            0,
        )
        self.assertEqual(sentinel["verdict"], "NON_DECISIVE")

    def test_schema_and_evidence_inventory_are_pinned(self):
        schema = json.loads((ROOT / "workflow/prophage_semantics/semantics-policy-v1.schema.json").read_text())
        evidence = json.loads((ROOT / "artifacts/prophage_semantics/evidence_inventory.json").read_text())
        self.assertEqual(schema["properties"]["schema_version"]["const"], "prophage-source-semantics-v1")
        self.assertEqual(evidence["source_provenance_result"]["status"], "MISSING")
        self.assertFalse(evidence["known_base_sentinels_used"])
        self.assertIn("not source attribution", evidence["items"][3]["claim"])


if __name__ == "__main__":
    unittest.main()
