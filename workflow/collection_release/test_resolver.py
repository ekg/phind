import argparse, gzip, hashlib, json, tempfile, unittest
from pathlib import Path
from unittest import mock

from workflow.collection_release import resolver as r

class ResolverTests(unittest.TestCase):
    def test_checksum_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"x"; p.write_bytes(b"wrong\n")
            with self.assertRaisesRegex(r.GateError,"checksum mismatch"):
                r.read_input(p,"0"*64)

    def test_external_inventory_mismatch_refused(self):
        with tempfile.TemporaryDirectory() as td:
            d=Path(td); (d/"payload").write_text("good"); (d/"SHA256SUMS").write_text(f"{r.sha_file(d/'payload')}  payload\n"); (d/"COMPLETE").write_text(f"{r.sha_file(d/'SHA256SUMS')}  SHA256SUMS\n")
            (d/"payload").write_text("corrupt")
            with self.assertRaisesRegex(r.GateError,"checksum mismatch"):
                r.verify_external(d)

    def test_resource_blank_and_overallocation_refused(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            blank=r.Resources(0,1,1,1,1,1,1)
            with self.assertRaisesRegex(r.GateError,"blank"):
                r.preflight(p,p,blank,False)
            over=r.Resources(1,100,100,2,80,2,1)
            with self.assertRaisesRegex(r.GateError,"NO_GO"):
                r.preflight(p,p,over,False)

    def test_exact_version_identity_rejects_unrequested_report(self):
        with tempfile.TemporaryDirectory() as td:
            d=Path(td)
            def fake(endpoint, accs, req, resp, receipt, retries, sleep):
                req.write_text("{}")
                obj={"valid_assemblies":accs} if endpoint==r.CHECK_URL else {"reports":[{"accession":"GCF_999999999.1"}]}
                resp.write_text(json.dumps(obj)); out={"response_sha256":r.sha_file(resp),"request_sha256":r.sha_file(req),"completed_at_utc":"x"}; receipt.write_text(json.dumps(out)); return out
            with self.assertRaisesRegex(r.GateError,"version identity violation"):
                self._tiny_freeze(d,fake)

    def test_lossless_freeze_and_atomic_resume(self):
        with tempfile.TemporaryDirectory() as td:
            d=Path(td)
            def fake(endpoint, accs, req, resp, receipt, retries, sleep):
                req.write_text(json.dumps({"accessions":accs}))
                if endpoint==r.CHECK_URL: obj={"valid_assemblies":accs}
                else: obj={"reports":[{"accession":a,"current_accession":a,"paired_accession":"GCA_"+a[4:],"organism":{"organism_name":"Escherichia coli","infraspecific_names":{"strain":"K-12"}},"assembly_info":{"assembly_status":"current","biosample":{"accession":"SAMN1"}},"annotation_info":{"name":"RefSeq"}} for a in accs]}
                resp.write_text(json.dumps(obj)); out={"response_sha256":r.sha_file(resp),"request_sha256":r.sha_file(req),"completed_at_utc":"x"}; receipt.write_text(json.dumps(out)); return out
            result=self._tiny_freeze(d,fake)
            final=Path(result["external_release_path"])
            self.assertTrue((final/"COMPLETE").is_file())
            self.assertEqual(result["counts"]["physical_occurrences"],3)
            rows=gzip.decompress((d/"tracked"/"occurrences.tsv.gz").read_bytes()).decode().splitlines()
            self.assertEqual(len(rows),4); self.assertIn("REJECTED_MALFORMED_TERMINAL_TOKEN",rows[-1])
            # Complete releases are immutable and accepted on resume without new requests.
            with mock.patch.object(r,"request_batch",side_effect=AssertionError("network called")):
                resumed=self._tiny_freeze(d,None)
            self.assertEqual(resumed["release_id"],result["release_id"])

    def test_interrupted_promotion_stage_is_cleanly_restarted(self):
        with tempfile.TemporaryDirectory() as td:
            d=Path(td)
            # A stale incomplete stage must not become a release input.
            (d/"artifacts").mkdir(); (d/"26k_ecoli_accession.txt").write_bytes(b"GCF_000000001.1\nGCF_000000002.2\ngenome\n"); (d/"26k_prophage1.csv").write_bytes(b"x\n"); (d/"artifacts/project_manifest_template.tsv").write_bytes(b"schema\n")
            seed=self._seed(d)
            rid="collection-v1-"+r.sha_bytes(r.canonical(seed))[:16]
            stale=d/"durable"/(".stage."+rid+".test-run"); stale.mkdir(parents=True); (stale/"partial").write_text("torn")
            def fake(endpoint, accs, req, resp, receipt, retries, sleep):
                req.write_text("{}")
                obj={"valid_assemblies":accs} if endpoint==r.CHECK_URL else {"reports":[]}
                resp.write_text(json.dumps(obj)); out={"response_sha256":r.sha_file(resp),"request_sha256":r.sha_file(req),"completed_at_utc":"x"}; receipt.write_text(json.dumps(out)); return out
            result=self._tiny_freeze(d,fake)
            self.assertFalse(stale.exists()); self.assertTrue((Path(result["external_release_path"])/"COMPLETE").exists())

    def _seed(self,d):
        access=(d/"26k_ecoli_accession.txt").read_bytes(); prop=(d/"26k_prophage1.csv").read_bytes(); templ=(d/"artifacts/project_manifest_template.tsv").read_bytes()
        return {"schema":r.SCHEMA,"accession_sha256":hashlib.sha256(access).hexdigest(),"prophage_sha256":hashlib.sha256(prop).hexdigest(),"template_sha256":hashlib.sha256(templ).hexdigest(),"exact_version_policy":"requested-version-is-identity","candidate_count":2}

    def _tiny_freeze(self,d,fake):
        if not (d/"26k_ecoli_accession.txt").exists():
            (d/"artifacts").mkdir(); (d/"26k_ecoli_accession.txt").write_bytes(b"GCF_000000001.1\nGCF_000000002.2\ngenome\n"); (d/"26k_prophage1.csv").write_bytes(b"x\n"); (d/"artifacts/project_manifest_template.tsv").write_bytes(b"schema\n")
        args=argparse.Namespace(repo_root=str(d),tracked_root=str(d/"tracked"),durable_root=str(d/"durable"),scratch_root=str(d/"scratch"),run_id="test-run",batch_size=2,retries=1,sleep_seconds=0,assigned_ram_bytes=1_000_000_000,durable_allocation_bytes=1000,scratch_allocation_bytes=1000,inode_allocation=100,predicted_peak_disk_bytes=10,predicted_files=10,unfinished_write_bytes=1,test_relaxed_floors=True)
        patches=[mock.patch.object(r,"EXPECTED_ACCESSION_SHA",r.sha_file(d/"26k_ecoli_accession.txt")),mock.patch.object(r,"EXPECTED_PROPHAGE_SHA",r.sha_file(d/"26k_prophage1.csv")),mock.patch.object(r,"EXPECTED_LINES",3),mock.patch.object(r,"EXPECTED_CANDIDATES",2),mock.patch.object(r,"STAGE_B",[("a",1),("b",2)])]
        for p in patches:p.start()
        try:
            if fake is None: return r.freeze(args)
            with mock.patch.object(r,"request_batch",side_effect=fake): return r.freeze(args)
        finally:
            for p in reversed(patches):p.stop()

if __name__=="__main__": unittest.main()
