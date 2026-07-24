import argparse
import hashlib
import io
import json
import shutil
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

from workflow.acquisition_canonicalization import pilot as p


class FakeHeaders(dict):
    def get(self, key, default=None):
        return super().get(key.lower(), default)


class FakeResponse:
    def __init__(self, data: bytes, status: int, url: str, headers: dict[str, str]):
        self._data = io.BytesIO(data)
        self.status = status
        self._url = url
        self.headers = FakeHeaders({key.lower(): value for key, value in headers.items()})

    def read(self, size=-1):
        return self._data.read(size)

    def geturl(self):
        return self._url

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def package_bytes(accession=p.EXPECTED_ACCESSIONS[0], unsafe=False, gff_end=19):
    token = "ctg#A" if unsafe else "NC_000001.1"
    fasta = f">{token} fixture\nACGTRYSWKMBDHVNACGT\n".encode()
    gff_seqid = "ctg%23A" if unsafe else token
    gff = f"##gff-version 3\n{gff_seqid}\tRefSeq\tregion\t1\t{gff_end}\t.\t+\t.\tID=x\n".encode()
    report = (json.dumps({"refseqAccession": token, "role": "assembled-molecule", "topology": "linear"}) + "\n").encode()
    catalog = (json.dumps({"assemblies": [{"accession": accession}]}) + "\n").encode()
    assembly_report = (json.dumps({"accession": accession}) + "\n").encode()
    payload = {
        "README.md": b"fixture\n",
        "ncbi_dataset/data/dataset_catalog.json": catalog,
        "ncbi_dataset/data/assembly_data_report.jsonl": assembly_report,
        f"ncbi_dataset/data/{accession}/{accession}_fixture_genomic.fna": fasta,
        f"ncbi_dataset/data/{accession}/genomic.gff": gff,
        f"ncbi_dataset/data/{accession}/sequence_report.jsonl": report,
    }
    md5 = "".join(f"{hashlib.md5(data, usedforsecurity=False).hexdigest()}  {name}\n" for name, data in payload.items()).encode()
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in payload.items():
            archive.writestr(name, data)
        archive.writestr("md5sum.txt", md5)
    return output.getvalue()


class PilotTests(unittest.TestCase):
    def test_manifest_checksum_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            tracked = root / "manifests/collection-v1"
            tracked.mkdir(parents=True)
            (tracked / "release.json").write_text("{}\n")
            (tracked / "stage_b_10.tsv").write_text("wrong\n")
            (tracked / "SHA256SUMS").write_text(
                f"{p.COLLECTION_RELEASE_SHA256}  release.json\n{p.STAGE_B_SHA256}  stage_b_10.tsv\n"
            )
            with self.assertRaisesRegex(p.GateError, "tracked checksum mismatch"):
                p.verify_predecessor(root)

    def test_resource_blank_disk_and_inode_refusal(self):
        blank = p.Allocations(0, 1, 1, 1, 1, 1, 1, 1)
        with self.assertRaisesRegex(p.GateError, "blank"):
            blank.validate()
        disk = p.Allocations(100, 100, 100, 100, 71, 1, 1, 1)
        with self.assertRaisesRegex(p.GateError, "durable"):
            disk.validate()
        inodes = p.Allocations(100, 100, 100, 4, 1, 1, 3, 1)
        with self.assertRaisesRegex(p.GateError, "inode"):
            inodes.validate()

    def test_payload_bound_rejects_eleventh_or_versionless_identity(self):
        with self.assertRaisesRegex(p.GateError, "outside immutable"):
            p.download_url("GCF_999999999.1")
        with self.assertRaisesRegex(p.GateError, "outside immutable"):
            p.download_url("GCF_000005845")
        self.assertIn("GCF_000005845.2", p.download_url("GCF_000005845.2"))

    def test_package_archive_upstream_md5_and_exact_identity(self):
        accession = p.EXPECTED_ACCESSIONS[0]
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "package.zip"
            path.write_bytes(package_bytes(accession))
            result = p.validate_package(path, accession)
            self.assertEqual(result["archive_validation"], "PASS")
            self.assertEqual(result["upstream_md5_validation"], "PASS")
            with self.assertRaisesRegex(p.GateError, "accession directory identity"):
                p.validate_package(path, p.EXPECTED_ACCESSIONS[1])
            # Corrupt one member while leaving upstream md5 unchanged.
            corrupt = Path(td) / "corrupt.zip"
            with zipfile.ZipFile(path) as source, zipfile.ZipFile(corrupt, "w") as target:
                for info in source.infolist():
                    data = source.read(info.filename)
                    if info.filename.endswith("genomic.gff"):
                        data += b"x"
                    target.writestr(info.filename, data)
            with self.assertRaisesRegex(p.GateError, "upstream MD5 mismatch"):
                p.validate_package(corrupt, accession)

    def test_safe_range_resume_requires_unchanged_strong_identity(self):
        accession = p.EXPECTED_ACCESSIONS[0]
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            state, failures = root / "state.jsonl", root / "failures.jsonl"
            partial = root / "package.zip.partial"
            partial.write_bytes(b"abc")
            url = p.download_url(accession)
            identity = {
                "requested_url": url,
                "final_url": url,
                "status": 200,
                "headers": {"etag": '"fixed"', "last_modified": ".", "accept_ranges": "bytes"},
                "total_bytes": 6,
                "identity_strength": "STRONG_ETAG",
            }
            (root / "remote_identity.json").write_bytes(p.canonical_json(identity))
            response = FakeResponse(
                b"def", 206, url,
                {"ETag": '"fixed"', "Accept-Ranges": "bytes", "Content-Range": "bytes 3-5/6", "Content-Length": "3"},
            )
            with mock.patch.object(p, "_head_identity", return_value=identity), mock.patch.object(p.urllib.request, "urlopen", return_value=response):
                package, receipt = p.acquire_package(accession, root, state, failures, 1, 0, False, 10)
            self.assertEqual(package.read_bytes(), b"abcdef")
            self.assertEqual(receipt["range_resumes_this_invocation"], 1)
            self.assertIn("ACQUISITION_RANGE_RESUME", state.read_text())

    def test_unsafe_partial_identity_restarts_without_append(self):
        accession = p.EXPECTED_ACCESSIONS[0]
        body = package_bytes(accession)
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "package.zip.partial").write_bytes(b"torn-not-prefix")
            weak = {"requested_url": p.download_url(accession), "final_url": p.download_url(accession), "status": 200,
                    "headers": {"etag": ".", "last_modified": ".", "accept_ranges": "."}, "total_bytes": None,
                    "identity_strength": "INSUFFICIENT_FOR_RANGE_RESUME"}
            (root / "remote_identity.json").write_bytes(p.canonical_json(weak))
            response = FakeResponse(body, 200, p.download_url(accession), {"Content-Length": str(len(body)), "Content-Type": "application/zip"})
            with mock.patch.object(p, "_head_identity", return_value=None), mock.patch.object(p.urllib.request, "urlopen", return_value=response):
                package, _ = p.acquire_package(accession, root, root / "state.jsonl", root / "failures.jsonl", 1, 0, False, 10)
            self.assertEqual(package.read_bytes(), body)
            self.assertIn("PARTIAL_RESTART_IDENTITY_UNSAFE", (root / "state.jsonl").read_text())

    @unittest.skipUnless(shutil.which("bgzip") and shutil.which("samtools"), "requires bgzip and samtools")
    def test_streamed_rename_bgzf_indexes_and_alias_roundtrip(self):
        accession = p.EXPECTED_ACCESSIONS[0]
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / "source"
            source.mkdir()
            package = source / "package.zip"
            package.write_bytes(package_bytes(accession, unsafe=True))
            validated = p.validate_package(package, accession)
            manifest = {"validation": validated}
            final, result, rows = p.canonicalize_object(
                accession, 1, source, manifest, root / "canonical", root / "state.jsonl",
                shutil.which("bgzip"), shutil.which("samtools"), 1, 6, False, 10,
            )
            self.assertTrue((final / "COMPLETE").exists())
            self.assertEqual(result["rename_only_validation"], "PASS")
            self.assertEqual(rows[0]["pansn_sequence_name"], f"{accession}#1#ctg%23A")
            alias = p.read_tsv(final / "annotation_aliases.tsv", verify_hashes=True)[0]
            self.assertEqual(alias["canonical_gff_seqid_lexical"], f"{accession}%231%23ctg%2523A")
            self.assertTrue(Path(str(final / f"{accession}.pansn.fa.gz") + ".gzi").exists())

    @unittest.skipUnless(shutil.which("bgzip") and shutil.which("samtools"), "requires bgzip and samtools")
    def test_invalid_source_gff_is_quarantined_without_alias_view(self):
        accession = p.EXPECTED_ACCESSIONS[0]
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / "source"
            source.mkdir()
            package = source / "package.zip"
            package.write_bytes(package_bytes(accession, gff_end=20))
            final, result, _ = p.canonicalize_object(
                accession, 1, source, {"validation": p.validate_package(package, accession)}, root / "canonical",
                root / "state.jsonl", shutil.which("bgzip"), shutil.which("samtools"), 1, 6, False, 10,
            )
            self.assertEqual(result["annotation"]["status"], "QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW")
            self.assertEqual(p.read_tsv(final / "annotation_aliases.tsv", verify_hashes=True), [])
            self.assertIn("ANNOTATION_ALIAS_VIEW_QUARANTINED", (root / "state.jsonl").read_text())

    @unittest.skipUnless(shutil.which("bgzip") and shutil.which("samtools"), "requires bgzip and samtools")
    def test_interrupted_conversion_stage_is_discarded_then_committed(self):
        accession = p.EXPECTED_ACCESSIONS[0]
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / "source"
            source.mkdir()
            package = source / "package.zip"
            package.write_bytes(package_bytes(accession))
            canonical = root / "canonical"
            stale = canonical / f".stage.{accession}"
            stale.mkdir(parents=True)
            (stale / "partial.fa.gz").write_bytes(b"torn")
            final, _, _ = p.canonicalize_object(
                accession, 1, source, {"validation": p.validate_package(package, accession)}, canonical,
                root / "state.jsonl", shutil.which("bgzip"), shutil.which("samtools"), 1, 6, False, 10,
            )
            self.assertFalse(stale.exists())
            self.assertTrue((final / "COMPLETE").exists())
            self.assertIn("INTERRUPTED_CONVERSION_STAGE_DISCARDED", (root / "state.jsonl").read_text())

    def test_interrupted_release_promotion_never_exposes_partial(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            stage, final = root / ".stage.release", root / "release"
            stage.mkdir()
            (stage / "payload").write_text("valid")
            # Simulated kill: final does not exist and staging has no COMPLETE.
            self.assertFalse(final.exists())
            self.assertFalse((stage / "COMPLETE").exists())
            p.seal_directory(stage, final)
            self.assertFalse(stage.exists())
            self.assertTrue((final / "COMPLETE").exists())
            p.verify_sha_inventory(final)

    def test_nested_object_inventories_are_covered_by_release_inventory(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            stage, final = root / ".stage", root / "final"
            nested = stage / "object"
            nested.mkdir(parents=True)
            (nested / "payload").write_text("x")
            p.seal_directory(nested, stage / "sealed-object")
            p.seal_directory(stage, final)
            inventory = (final / "SHA256SUMS").read_text()
            self.assertIn("sealed-object/COMPLETE", inventory)
            self.assertIn("sealed-object/SHA256SUMS", inventory)


if __name__ == "__main__":
    unittest.main()
