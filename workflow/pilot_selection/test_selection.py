from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from workflow.pilot_selection import selection as s


def fake_row(number: int, status: str = "current", version: int = 1) -> dict[str, str]:
    accession = f"GCF_{number:09d}.{version}"
    return {
        "input_line_number": str(number),
        "input_occurrence_id": f"occ-{number}",
        "assembly_id": f"asm-{number}",
        "requested_assembly_accession_version": accession,
        "resolved_assembly_accession_version": accession,
        "resolution_status": "EXACT_VERSION_RESOLVED",
        "assembly_status": status,
        "assembly_version": str(version),
        "row_sha256": f"source-{number}",
        # Deliberately forbidden traits must have no effect on selection.
        "taxonomy": "forbidden-a",
        "phage_count": "999",
    }


def stage_row(row: dict[str, str], order: int) -> dict[str, str]:
    return {
        "stage_b_order": str(order),
        "assembly_id": row["assembly_id"],
        "resolved_assembly_accession_version": row["resolved_assembly_accession_version"],
    }


class SelectionTests(unittest.TestCase):
    def test_exact_nested_prefixes_and_suppressed_full_frame_mapping(self):
        rows = [fake_row(i) for i in range(1, 1101)]
        rows[-1]["assembly_status"] = "suppressed"
        stage = [stage_row(row, i) for i, row in enumerate(rows[:10], 1)]
        result = s.select_cohort(rows, stage)
        self.assertEqual(len(result.cohort), 1000)
        self.assertEqual(len(result.frame), 1100)
        self.assertEqual(result.frame_counts["terminal_suppressed_ineligible"], 1)
        for n in s.RUNGS:
            self.assertEqual(len(result.rungs[n]), n)
            self.assertEqual(
                [row["assembly_id"] for row in result.rungs[n]],
                [row["assembly_id"] for row in result.cohort[:n]],
            )
        self.assertEqual(
            [row["assembly_id"] for row in result.cohort[:10]],
            [row["assembly_id"] for row in rows[:10]],
        )
        self.assertEqual(len({row["assembly_id"] for row in result.cohort}), 1000)

    def test_selection_is_phage_blind_and_seed_deterministic(self):
        rows = [fake_row(i) for i in range(1, 1021)]
        stage = [stage_row(row, i) for i, row in enumerate(rows[:10], 1)]
        first = s.select_cohort(rows, stage)
        for row in rows:
            row["taxonomy"] = "changed"
            row["phage_count"] = str(int(row["phage_count"]) + 1)
            row["transposable"] = "1.0"
        second = s.select_cohort(list(reversed(rows)), stage)
        self.assertEqual(
            [row["assembly_id"] for row in first.cohort],
            [row["assembly_id"] for row in second.cohort],
        )
        self.assertTrue(all("taxonomy" not in row and "phage_count" not in row for row in first.frame))

    def test_inclusion_probabilities_are_exact_design_fractions(self):
        rows = [fake_row(i) for i in range(1, 1201)]
        stage = [stage_row(row, i) for i, row in enumerate(rows[:10], 1)]
        result = s.select_cohort(rows, stage)
        certainty = result.rungs[100][0]
        random_unit = result.rungs[100][10]
        self.assertEqual(certainty["inclusion_probability"], "1/1")
        self.assertEqual(certainty["inference_weight"], "1/1")
        self.assertEqual(random_unit["inclusion_probability"], "90/1190")
        self.assertEqual(random_unit["inference_weight"], "1190/90")

    def test_exact_version_mismatch_and_stage_b_drift_fail_closed(self):
        rows = [fake_row(i) for i in range(1, 1021)]
        rows[20]["resolved_assembly_accession_version"] = "GCF_999999999.1"
        stage = [stage_row(row, i) for i, row in enumerate(rows[:10], 1)]
        with self.assertRaisesRegex(s.GateError, "exact-version"):
            s.select_cohort(rows, stage)
        rows[20] = fake_row(21)
        stage[0]["assembly_id"] = "wrong"
        with self.assertRaisesRegex(s.GateError, "Stage-B"):
            s.select_cohort(rows, stage)

    def test_manifest_checksum_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "payload.tsv").write_text("x\n")
            (root / "SHA256SUMS").write_text(f"{'0' * 64}  payload.tsv\n")
            (root / "COMPLETE").write_text(f"{s.sha_file(root / 'SHA256SUMS')}  SHA256SUMS\n")
            with self.assertRaisesRegex(s.GateError, "checksum mismatch"):
                s.verify_inventory(root)

    def test_resource_blank_disk_inode_and_unfinished_write_refusal(self):
        with self.assertRaisesRegex(s.GateError, "blank"):
            s.Allocations(0, 1, 1, 1, 1, 1, 1, 1).validate()
        with self.assertRaisesRegex(s.GateError, "durable"):
            s.Allocations(100, 100, 100, 100, 71, 1, 1, 1).validate()
        with self.assertRaisesRegex(s.GateError, "inode"):
            s.Allocations(100, 100, 100, 4, 1, 1, 3, 1).validate()
        with self.assertRaisesRegex(s.GateError, "unfinished"):
            s.Allocations(100, 100, 100, 100, 1, 1, 1, 60).validate()

    def test_interrupted_promotion_is_never_visible_and_restart_discards_stage(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            stage, final, state = root / ".stage.release", root / "release", root / "state.jsonl"
            stage.mkdir()
            (stage / "manifest.tsv").write_text("id\n1\n")
            with self.assertRaises(s.InjectedInterruption):
                s.seal_and_promote(stage, final, inject=lambda: (_ for _ in ()).throw(s.InjectedInterruption()))
            self.assertFalse(final.exists())
            self.assertFalse((stage / "COMPLETE").exists())
            self.assertTrue((stage / "SHA256SUMS").exists())
            s.discard_interrupted_stage(stage, state)
            self.assertFalse(stage.exists())
            self.assertIn("INTERRUPTED_PUBLICATION_STAGE_DISCARDED", state.read_text())
            stage.mkdir()
            (stage / "manifest.tsv").write_text("id\n1\n")
            s.seal_and_promote(stage, final)
            s.verify_inventory(final)

    def test_engineering_controls_do_not_add_or_reorder_assemblies(self):
        rows = [fake_row(i) for i in range(1, 1021)]
        stage = [stage_row(row, i) for i, row in enumerate(rows[:10], 1)]
        result = s.select_cohort(rows, stage)
        controls = s.synthetic_engineering_controls(stage, "blocked-release")
        self.assertTrue(any(row["control_class"] == "unsafe_source_contig_id" for row in controls))
        self.assertTrue(any(row["activation_status"] == "BLOCKED_EXTRACTION_SEMANTICS" for row in controls))
        self.assertTrue(all(row["selection_effect"] == "NONE_POST_SELECTION_LABEL" for row in controls))
        self.assertEqual(len(result.cohort), 1000)


if __name__ == "__main__":
    unittest.main()
