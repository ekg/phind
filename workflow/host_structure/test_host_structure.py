import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from workflow.host_structure import host_structure as h
from workflow.host_structure import runner


class HostStructureContractTests(unittest.TestCase):
    def test_real_frozen_inputs_are_exact_1000_and_host_tools_pass(self):
        inputs = h.verify_inputs(Path(".").resolve())
        self.assertEqual(len(inputs["accessions"]), 1000)
        self.assertEqual(len(set(inputs["accessions"])), 1000)
        self.assertEqual(inputs["canonical"]["release_id"], h.CANONICAL_RELEASE_ID)
        self.assertEqual(inputs["compatibility"]["verdict"], h.PASS)
        self.assertEqual(h.canonical_object_path(inputs["refs"][500], inputs["accessions"][500]).parent.parent.parent,
                         h.CANONICAL_EXTERNAL_PATH)

    def test_manifest_checksum_mismatch_fails_closed(self):
        original = Path("manifests/canonical-cohort-1000-v1/cohort-1000.tsv").read_text()
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "cohort.tsv"
            lines = original.splitlines()
            fields = lines[1].split("\t")
            fields[3] = "GCF_999999999.1"
            lines[1] = "\t".join(fields)
            p.write_text("\n".join(lines) + "\n")
            with self.assertRaisesRegex(h.GateError, "row checksum mismatch"):
                h.read_hashed_tsv(p)

    def test_host_only_contract_rejects_phage_trait_or_artifact(self):
        for value in (
            {"prophage_count": 1}, {"input": "/tmp/phage_clusters.tsv"},
            {"taxonomy": "virus"}, {"coordinates": [1, 2]},
        ):
            with self.subTest(value=value), self.assertRaisesRegex(h.GateError, "phage"):
                h.assert_host_only_manifest(value)
        h.assert_host_only_manifest({
            "root": {"26k_prophage1.csv": h.ROOT_HASHES["26k_prophage1.csv"]},
            "prophage_source_coordinate_policy": h.NA_HOST,
            "prophage_extraction_semantics": h.NA_HOST,
            "host_input": "manifests/canonical-cohort-1000-v1/cohort-1000.tsv",
        })

    def test_resource_blank_overallocation_and_low_live_floor_refuse(self):
        with self.assertRaisesRegex(h.GateError, "blank"):
            h.Allocations(0, 1, 1, 1, 1, 1, 1, 1).validate()
        with self.assertRaisesRegex(h.GateError, "70%"):
            h.Allocations(100, 100, 100, 100, 71, 1, 1, 1).validate()
        allocation = h.Allocations(100, 1000, 1000, 100, 10, 10, 10, 1)
        low = mock.Mock(f_bavail=1, f_frsize=1, f_favail=10)
        with tempfile.TemporaryDirectory() as td, mock.patch.object(h.os, "statvfs", return_value=low), \
                mock.patch.object(h, "_mount_record", return_value={"target": "/", "source": "x", "fstype": "x"}):
            with self.assertRaisesRegex(h.GateError, "resource gate NO_GO"):
                h.live_preflight(Path(td) / "d", Path(td) / "s", allocation, "test")

    def test_interrupted_promotion_has_no_complete_or_final(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            stage, final = root / ".staging", root / "final"
            stage.mkdir(); (stage / "payload").write_text("valid")
            self.assertFalse(final.exists())
            self.assertFalse((stage / "COMPLETE").exists())
            h.seal_directory(stage, final)
            self.assertFalse(stage.exists())
            self.assertTrue((final / "COMPLETE").is_file())
            h.verify_external_inventory(final)

    def test_release_id_is_order_parameter_and_pin_deterministic(self):
        first = runner.release_id(); second = runner.release_id()
        self.assertEqual(first, second)
        self.assertTrue(first.startswith("host-structure-1000-v1-"))
        old = runner.PARAMETERS["near_duplicate_distance"]
        try:
            runner.PARAMETERS["near_duplicate_distance"] = old + 1e-8
            self.assertNotEqual(first, runner.release_id())
        finally:
            runner.PARAMETERS["near_duplicate_distance"] = old

    def test_newick_tip_splits_support_and_collapse(self):
        tree = h.parse_newick("((A:1,B:1):1,(C:1,D:1):1);")
        self.assertEqual(set(h.leaf_names(tree)), {"A", "B", "C", "D"})
        splits = h.tree_splits(tree)
        self.assertEqual(len(splits), 1)
        collapsed = h.collapse_unsupported(tree, {}, 0.95)
        self.assertEqual(set(h.leaf_names(collapsed)), {"A", "B", "C", "D"})
        self.assertIn("A", h.newick_string(collapsed))

    def test_directed_distance_exact_symmetry_and_mismatch(self):
        labels = ["A", "B", "C"]
        matrix = [[0, .1, .2], [.1, 0, .3], [.2, .3, 0]]
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "dist.tsv"
            with p.open("w") as out:
                for i, a in enumerate(labels):
                    for j, b in enumerate(labels):
                        out.write(f"{a}.fa\t{b}.fa\t{matrix[i][j]}\t0\t1/1\n")
            result = h.validate_directed_mash(p, labels, matrix)
            self.assertEqual(result["directed_records"], 9)
            text = p.read_text().replace("A.fa\tB.fa\t0.1", "A.fa\tB.fa\t0.2")
            p.write_text(text)
            with self.assertRaisesRegex(h.GateError, "mismatch"):
                h.validate_directed_mash(p, labels, matrix)

    def test_neighbor_joining_bootstrap_core_and_recombination_mask(self):
        labels = ["A", "B", "C", "D"]
        alignment = [b"AAAACCCCGGGG", b"AAAACCCCGGGA", b"TTTTCCCCGGGG", b"TTTTCCCCGGGA"]
        matrix = h.p_distance_matrix(alignment)
        tree = h.neighbor_joining(labels, matrix)
        self.assertEqual(set(h.leaf_names(tree)), set(labels))
        support = h.bootstrap_splits(labels, alignment, [True] * 12, 10, seed=7)
        self.assertTrue(all(0 <= x <= 1 for x in support.values()))
        calls = [bytearray(x) for x in alignment]
        core, positions, stats = h.core_alignment(calls, .75)
        self.assertEqual(stats["core_callable_sites"], 12)
        mask, mstats = h.recombination_candidate_mask(core, positions, z=1, window=4)
        self.assertEqual(len(mask), 12)
        self.assertIn("masked_fraction", mstats)

    def test_paf_reference_call_roundtrip_forward_and_reverse(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            ref = root / "R.fa"; q = root / "Q.fa"; paf = root / "x.paf"
            ref.write_text(">R#1#c\nACGTACGT\n")
            q.write_text(">Q#1#c\nACGTACGT\n")
            paf.write_text("Q#1#c\t8\t0\t8\t+\tR#1#c\t8\t0\t8\t8\t8\t60\tcg:Z:8M\n")
            calls, stats, coords = h.build_reference_calls(paf, ref, {"Q": q, "R": ref}, ["R", "Q"], min_block=1)
            self.assertEqual(bytes(calls[0]), b"ACGTACGT")
            self.assertEqual(bytes(calls[1]), b"ACGTACGT")
            self.assertEqual(stats["paf_records_accepted"], 1)
            self.assertEqual(coords[0][0], "R#1#c")

    def test_alternative_partitions_and_near_duplicate_units_are_deterministic(self):
        labels = ["A", "B", "C", "D", "E"]
        matrix = [[abs(i-j)/10 for j in range(5)] for i in range(5)]
        medoids1, assign1 = h.farthest_first_partition(labels, matrix, 2)
        medoids2, assign2 = h.farthest_first_partition(labels, matrix, 2)
        self.assertEqual((medoids1, assign1), (medoids2, assign2))
        groups = h.union_find_groups(labels, matrix, .11)
        self.assertEqual(sum(map(len, groups)), 5)


if __name__ == "__main__":
    unittest.main()
