import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from workflow.compatibility.compatibility import (
    GateError, ResourceRequest, atomic_promote, gff_lexical_decode,
    gff_lexical_encode, mash_lower_to_full_phylip, percent_decode_identifier, percent_encode_identifier,
    sha256_file, stage_gff_semantic_alias, verify_inventory, verify_predecessor,
    write_inventory,
)


class NamingTests(unittest.TestCase):
    def test_unsafe_pansn_and_gff_layers_round_trip(self):
        raw = "ctg#A/plasmid|β and space%23".encode()
        encoded = percent_encode_identifier(raw)
        self.assertEqual(percent_decode_identifier(encoded), raw)
        self.assertIn("%23", encoded)
        self.assertIn("%2F", encoded)
        self.assertIn("%20", encoded)
        semantic = "GCF_000005845.2#1#" + encoded
        lexical = gff_lexical_encode(semantic)
        self.assertEqual(gff_lexical_decode(lexical), semantic)
        self.assertIn("%23", lexical)       # semantic PanSN delimiter
        self.assertIn("%2523", lexical)     # encoded raw delimiter, escaped again

    def test_noncanonical_lowercase_escape_rejected(self):
        with self.assertRaises(GateError):
            percent_decode_identifier("bad%2fslash")
        with self.assertRaises(GateError):
            gff_lexical_decode("sample%231%23bad%2fslash")

    def test_gff_semantic_alias_preserves_coordinates_and_rejects_orphan(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            semantic = "GCF_000005845.2#1#NC_000913.3"
            lexical = gff_lexical_encode(semantic)
            source = root / "source.gff3"
            source.write_text(f"##gff-version 3\n{lexical}\tx\tgene\t11\t20\t.\t-\t.\tID=g1\n")
            out = root / "alias.gff3"
            mapping = root / "map.json"
            rows = stage_gff_semantic_alias(source, out, {semantic}, mapping)
            self.assertEqual(rows, [{"lexical": lexical, "semantic": semantic}])
            fields = out.read_text().splitlines()[1].split("\t")
            self.assertEqual(fields[0], semantic)
            self.assertEqual(fields[3:5], ["11", "20"])
            with self.assertRaises(GateError):
                stage_gff_semantic_alias(source, root / "bad.gff", {"other"}, root / "bad.json")


class InventoryTests(unittest.TestCase):
    def test_package_locks_have_matching_sha256_inventories(self):
        repo = Path(__file__).resolve().parents[2]
        for stem in ("environment", "graph"):
            lock = repo / f"workflow/compatibility/{stem}-linux-64.explicit.lock"
            sha = repo / f"workflow/compatibility/{stem}-package-sha256.tsv"
            locked = {line.split("#", 1)[0]: line.rsplit("#", 1)[1]
                      for line in lock.read_text().splitlines() if line.startswith("https://")}
            rows = [line.split("\t") for line in sha.read_text().splitlines()[1:] if line]
            self.assertEqual(set(locked), {row[0] for row in rows})
            self.assertTrue(all(len(row[1]) == 64 and len(row[2]) == 32 for row in rows))
            self.assertEqual(locked, {row[0]: row[2] for row in rows})

    def test_mash_lower_to_rapidnj_full_matrix_and_malformed_refusal(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            lower, full = root / "lower", root / "full"
            lower.write_text("3\na\nb 0.1\nc 0.2 0.3\n")
            result = mash_lower_to_full_phylip(lower, full)
            self.assertEqual(result["off_diagonal_pairs"], 3)
            rows = full.read_text().splitlines()
            self.assertEqual(rows[0], "3")
            self.assertEqual(rows[1].split(), ["a", "0", "0.1", "0.2"])
            self.assertEqual(rows[2].split(), ["b", "0.1", "0", "0.3"])
            lower.write_text("3\na 0\nb 0.1\nc 0.2 0.3\n")
            with self.assertRaises(GateError):
                mash_lower_to_full_phylip(lower, full)

    def test_checksum_mismatch_and_path_traversal_fail_closed(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "x").write_text("ok\n")
            write_inventory(root)
            verify_inventory(root, root / "SHA256SUMS", exact=True)
            (root / "x").write_text("changed\n")
            with self.assertRaises(GateError):
                verify_inventory(root, root / "SHA256SUMS")
            (root / "SHA256SUMS").write_text("0" * 64 + "  ../escape\n")
            with self.assertRaises(GateError):
                verify_inventory(root, root / "SHA256SUMS")

    def test_manifest_checksum_mutation_rejected_before_consumption(self):
        repo = Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory(dir=repo) as td:
            fake = Path(td)
            for name in ("26k_ecoli_accession.txt", "26k_prophage1.csv"):
                os.link(repo / name, fake / name)
            shutil.copytree(repo / "manifests/canonical-cohort-010-v1", fake / "manifests/canonical-cohort-010-v1")
            release = fake / "manifests/canonical-cohort-010-v1/release.json"
            doc = json.loads(release.read_text())
            doc["release_id"] = "substituted-release"
            release.write_text(json.dumps(doc) + "\n")
            with self.assertRaises(GateError):
                verify_predecessor(fake)


class ResourceAndPromotionTests(unittest.TestCase):
    def good_request(self):
        return ResourceRequest(
            assigned_ram_bytes=8_000_000_000,
            durable_allocation_bytes=5_000_000_000,
            scratch_allocation_bytes=4_000_000_000_000,
            inode_allocation=100_000,
            predicted_durable_peak_bytes=1_000_000_000,
            predicted_scratch_peak_bytes=3_000_000_000,
            predicted_files=5_000,
            unfinished_write_bytes=500_000_000,
        )

    def test_blank_and_overallocation_refused(self):
        r = self.good_request()
        r.validate()
        for bad in (
            ResourceRequest(0, r.durable_allocation_bytes, r.scratch_allocation_bytes, r.inode_allocation,
                            r.predicted_durable_peak_bytes, r.predicted_scratch_peak_bytes, r.predicted_files,
                            r.unfinished_write_bytes),
            ResourceRequest(r.assigned_ram_bytes, r.durable_allocation_bytes, r.scratch_allocation_bytes,
                            r.inode_allocation, 4_000_000_000, r.predicted_scratch_peak_bytes,
                            r.predicted_files, r.unfinished_write_bytes),
            ResourceRequest(r.assigned_ram_bytes, r.durable_allocation_bytes, r.scratch_allocation_bytes,
                            r.inode_allocation, r.predicted_durable_peak_bytes, r.predicted_scratch_peak_bytes,
                            60_000, r.unfinished_write_bytes),
        ):
            with self.assertRaises(GateError):
                bad.validate()

    def test_interrupted_promotion_never_publishes_and_clean_restart_does(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            staging = root / ".release.staging"
            final = root / "release"
            staging.mkdir()
            (staging / "payload").write_text("partial")
            with self.assertRaises(GateError):
                atomic_promote(staging, final)
            self.assertFalse(final.exists())
            shutil.rmtree(staging)
            staging.mkdir()
            (staging / "payload").write_text("complete")
            write_inventory(staging)
            (staging / "COMPLETE").write_text('{"verdict":"PASS"}\n')
            atomic_promote(staging, final)
            self.assertTrue((final / "COMPLETE").is_file())
            verify_inventory(final, final / "SHA256SUMS", exact=True)

    def test_actual_sigkill_worker_leaves_only_unpublished_partial(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            staging, final = root / ".staging", root / "final"
            proc = subprocess.run([sys.executable, "-m", "workflow.compatibility.pilot", "_promotion-worker",
                                   "--staging", str(staging), "--final", str(final)])
            self.assertEqual(proc.returncode, -signal.SIGKILL)
            self.assertTrue(staging.is_dir())
            self.assertFalse((staging / "COMPLETE").exists())
            self.assertFalse(final.exists())


class PredecessorRegressionTests(unittest.TestCase):
    def test_exact_predecessor_all_objects_and_rows_pass(self):
        repo = Path(__file__).resolve().parents[2]
        release, assemblies, evidence = verify_predecessor(repo)
        self.assertEqual(release["release_id"], "canonical-cohort-010-v1-e71484de9994fc28")
        self.assertEqual(len(assemblies), 10)
        self.assertEqual(evidence["contigs"], 1223)
        self.assertEqual(evidence["bases"], 51731662)


if __name__ == "__main__":
    unittest.main()
