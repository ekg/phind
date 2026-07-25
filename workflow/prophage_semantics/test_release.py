import json
import tempfile
import unittest
from pathlib import Path

from workflow.prophage_semantics import release


ROOT = Path(__file__).resolve().parents[2]


def _current_release_path() -> Path | None:
    ref = ROOT / "artifacts/prophage_semantics/release_reference.json"
    if not ref.is_file():
        return None
    import json
    external = json.loads(ref.read_text()).get("external_path")
    p = Path(external) if external else None
    return p if p and p.is_dir() else None


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

    def test_v1_policy_is_the_blocked_historical_record(self):
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

    def test_v2_policy_is_decisive_go_with_resolved_coordinates(self):
        policy = json.loads((ROOT / release.POLICY_FILE).read_text())
        release.validate_policy(policy)
        gate = policy["extraction_gate"]
        self.assertEqual(gate["verdict"], "EXTRACTION_GO")
        self.assertEqual(gate["consumer_action"], "ALLOW")
        self.assertEqual(gate["blocking_dimensions"], [])
        self.assertEqual(gate["selected_coordinate_candidate"], "C1_RAW_1_BASED_CLOSED")
        # 'tagged' must remain an explicit unresolved user term, never silently relabeled
        tagged = {d["name"]: d for d in policy["semantic_dimensions"]}["tagged"]
        self.assertEqual(tagged["status"], "UNRESOLVED")
        self.assertFalse(tagged["extraction_critical"])
        self.assertFalse(any(s["generic_tagged_alias"] for s in policy["scopes"]))
        # the three scopes remain distinct
        self.assertEqual([s["id"] for s in policy["scopes"]],
                         ["all_records", "transposable_flag_positive", "taxonomy_assigned"])
        # both extraction-critical dimensions are resolved
        for dim in policy["semantic_dimensions"]:
            if dim["extraction_critical"]:
                self.assertTrue(dim["status"].startswith("RESOLVED"), dim["name"])
        # modern v2.4 pilot separation is documented
        self.assertIn("NOT the basis", gate["modern_v2_4_pilot_separation"])

    def test_derive_extraction_fails_closed_without_decisive_gate(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(release.GateError, "fail closed"):
                release.derive_extraction_from_gate(Path(td))

    def test_derive_extraction_fails_closed_on_non_decisive_gate(self):
        with tempfile.TemporaryDirectory() as td:
            repo = Path(td)
            gate = {
                "schema": "pinned-phigaro-consumption-gate-v1",
                "verdict": "PASS",
                "release_id": release.PREDECESSOR_PHIGARO_RELEASE_ID,
                "modern_v2_4_pilot_separate": True,
                "historical_csv_attribution": "NON_DECISIVE",
                "decisive_evidence_independently_sound": True,
                "historical_csv_extraction": "EXTRACTION_BLOCKED",
            }
            (repo / release.PINNED_CALLER_GATE_FILE).parent.mkdir(parents=True, exist_ok=True)
            (repo / release.PINNED_CALLER_GATE_FILE).write_text(json.dumps(gate))
            with self.assertRaisesRegex(release.GateError, "not DECISIVE"):
                release.derive_extraction_from_gate(repo)

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
        # the source-GFF diagnostic is deliberately NON_DECISIVE provenance; the
        # decisive evidence now comes from the pinned-caller rerun, not this probe
        self.assertEqual(sentinel["verdict"], "NON_DECISIVE")

    def test_schema_and_evidence_inventory_are_pinned(self):
        schema = json.loads((ROOT / "workflow/prophage_semantics/semantics-policy-v1.schema.json").read_text())
        evidence = json.loads((ROOT / "artifacts/prophage_semantics/evidence_inventory.json").read_text())
        self.assertEqual(schema["properties"]["schema_version"]["const"], "prophage-source-semantics-v1")
        self.assertFalse(evidence["known_base_sentinels_used"])

    def test_v2_release_validates_as_extraction_go(self):
        release_path = _current_release_path()
        if release_path is None:
            self.skipTest("v2 release not published in this environment")
        validation = release.validate_release(release_path, require_go=True)
        self.assertEqual(validation["verdict"], "PASS")
        self.assertEqual(validation["release_verdict"], "EXTRACTION_GO")
        self.assertEqual(validation["extraction_consumer_gate"], "ALLOW")
        self.assertEqual(validation["policy_version"], "v2")
        self.assertEqual(validation["release_id"], release_path.name)


if __name__ == "__main__":
    unittest.main()
