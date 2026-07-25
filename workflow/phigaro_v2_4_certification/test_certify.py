from __future__ import annotations

import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

from workflow.phigaro_v2_4_certification.certify import (
    EXPECTED_ACCESSIONS,
    certify,
    parse_sha256sums,
)


class CertificationTests(unittest.TestCase):
    def make_repo(self, root: Path) -> None:
        (root / "manifests/pilot-cohorts-v1").mkdir(parents=True)
        with (root / "manifests/pilot-cohorts-v1/cohort-0010.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["cohort_order", "rung_n", "exact_assembly_accession_version"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for order, accession in enumerate(EXPECTED_ACCESSIONS, 1):
                writer.writerow(
                    {"cohort_order": order, "rung_n": 10, "exact_assembly_accession_version": accession}
                )
        (root / "26k_ecoli_accession.txt").write_text("opaque-root-one\n")
        (root / "26k_prophage1.csv").write_text("opaque-root-two\n")
        subprocess.run(["git", "init", "-q", str(root)], check=True)
        subprocess.run(["git", "-C", str(root), "config", "user.email", "audit@example.invalid"], check=True)
        subprocess.run(["git", "-C", str(root), "config", "user.name", "Audit Test"], check=True)
        subprocess.run(["git", "-C", str(root), "add", "manifests", "26k_ecoli_accession.txt", "26k_prophage1.csv"], check=True)
        subprocess.run(["git", "-C", str(root), "commit", "-qm", "fixture"], check=True)

    def test_missing_reference_rejects_without_release_writes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            repo.mkdir()
            self.make_repo(repo)
            release_parent = Path(tmp) / "external-releases"
            result = certify(repo, "artifacts/phigaro/release_reference.json", release_parent)
            self.assertEqual("REJECTED", result["verdict"])
            self.assertFalse(result["qualified_for_scaling"])
            self.assertEqual("FAIL", result["checks"]["tracked_release_reference"])
            self.assertEqual(0, result["side_effects"]["release_writes"])
            self.assertFalse(release_parent.exists())

    def test_tracked_reference_to_missing_release_rejects(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            repo.mkdir()
            self.make_repo(repo)
            relref = "artifacts/phigaro/release_reference.json"
            ref = repo / relref
            ref.parent.mkdir(parents=True)
            release_parent = Path(tmp) / "external-releases"
            release = release_parent / "release-1"
            ref.write_text(json.dumps({"external_release_path": str(release)}))
            subprocess.run(["git", "-C", str(repo), "add", relref], check=True)
            subprocess.run(["git", "-C", str(repo), "commit", "-qm", "reference"], check=True)
            result = certify(repo, relref, release_parent)
            self.assertEqual("PASS", result["checks"]["tracked_release_reference"])
            self.assertEqual("FAIL", result["checks"]["external_release_directory"])
            self.assertEqual("REJECTED", result["verdict"])

    def test_sha256_inventory_parser_rejects_escape_and_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sums = Path(tmp) / "SHA256SUMS"
            digest = "0" * 64
            sums.write_text(
                f"{digest}  ok.txt\n{digest}  ../escape\n{digest}  ok.txt\n{'z' * 64}  bad.txt\n"
            )
            entries, errors = parse_sha256sums(sums)
            self.assertEqual({"ok.txt": digest}, entries)
            self.assertEqual(3, len(errors))


if __name__ == "__main__":
    unittest.main()
