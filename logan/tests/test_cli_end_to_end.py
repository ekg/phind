#!/usr/bin/env python3
"""End-to-end CLI integration test: prep -> record -> fetch -> status.

Spins up a local HTTP server that mimics the documented kmviz download
endpoint (GET /api/download/<session> -> ZIP; unknown session -> 400 with the
kmviz error string). The fetch CLI is pointed at the local server, so this
test performs NO network access beyond loopback and never touches the real
service.
"""

from __future__ import annotations

import http.server
import io
import json
import os
import socketserver
import subprocess
import sys
import threading
import zipfile
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
LOGAN = HERE.parent
SEQ = ("GTGTCAGCTTTCGTGGTGTGCAGCTGGCGTCAGATGACAACATGCTGCCAGACAGCCTGAAAGGGTTTGC"
       "GCCTGTGGTGCGTGGTATCGCCAAAAGCAATGCCCAGATAACGATTAAGCAAAATGGTTACACCATTTAC"
       "CAAACTTATG")
SESSION = "kmviz-12345678-1234-1234-1234-123456789abc"


def make_result_zip() -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        cols = "acc\tkmer_coverage\tANI_estimation\torganism\n"
        zf.writestr("ecoli_smoke_t50.tsv", cols +
                    "SRR111\t0.987\t0.9995\tEscherichia coli\n"
                    "SRR110\t0.611\t0.9831\tEscherichia coli\n")
        zf.writestr("session.json", json.dumps({SESSION: {}}))
        zf.writestr("README.md", "# results")
    return buf.getvalue()


class Handler(http.server.BaseHTTPRequestHandler):
    zip_body = make_result_zip()

    def do_GET(self):  # noqa: N802
        if self.path.startswith("/api/download/"):
            sid = self.path.rsplit("/", 1)[1]
            if sid == SESSION:
                body = self.zip_body
                self.send_response(200)
                self.send_header("Content-Type", "application/zip")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            body = b"An error occured while processing your request"
            self.send_response(400)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        self.send_response(404)
        self.end_headers()

    def log_message(self, *a):  # silence
        pass


@pytest.fixture()
def server():
    srv = socketserver.TCPServer(("127.0.0.1", 0), Handler)
    port = srv.server_address[1]
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    yield f"http://127.0.0.1:{port}"
    srv.shutdown()
    srv.server_close()


def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(LOGAN / "run_pilot.py"), *map(str, args)],
        capture_output=True, text=True, cwd=str(LOGAN))


def test_cli_end_to_end(tmp_path, server):
    fasta = tmp_path / "baits.fa"
    fasta.write_text(f">ecoli_smoke\n{SEQ}\n")
    run_dir = tmp_path / "run"

    # dryrun: no network, validates and plans
    r = run_cli("dryrun", "--fasta", fasta, "--group", "GenBank_RefSeq",
                "--thresholds", "0.5")
    assert r.returncode == 0, r.stderr
    report = json.loads(r.stdout)
    assert report["planned_download_requests"] == 1
    assert report["dry_run"] is True

    # prep: refuses over-limit input before building anything
    (tmp_path / "toolong.fa").write_text(">x\n" + "ACGT" * 626 + "\n")
    r = run_cli("prep", "--run-dir", tmp_path / "bad", "--fasta",
                tmp_path / "toolong.fa")
    assert r.returncode == 2 and "too long" in r.stderr

    r = run_cli("prep", "--run-dir", run_dir, "--run-id", "it1", "--fasta",
                fasta, "--group", "GenBank_RefSeq", "--thresholds", "0.5")
    assert r.returncode == 0, r.stderr
    assert (run_dir / "manifest.json").exists()
    assert (run_dir / "INSTRUCTIONS.md").exists()
    sub = (run_dir / "submissions" / "ecoli_smoke_t50.fa")
    assert sub.exists() and sub.read_text().startswith(">ecoli_smoke_t50\n")

    # record: rejects bad session ids, accepts documented format (offline)
    r = run_cli("record", "--run-dir", run_dir, "--bait", "ecoli_smoke_t50",
                "--session", "bogus")
    assert r.returncode == 2 and "kmviz-" in r.stderr
    r = run_cli("record", "--run-dir", run_dir, "--bait", "ecoli_smoke_t50",
                "--session", SESSION)
    assert r.returncode == 0, r.stderr

    # status shows submitted before fetch
    r = run_cli("status", "--run-dir", run_dir)
    assert json.loads(r.stdout)["status_counts"] == {"submitted": 1}

    # fetch: unknown session 400 first (soft pending), then real fetch by
    # pointing the row at the good session; here we fetch the good one
    r = run_cli("fetch", "--run-dir", run_dir, "--base-url", server,
                "--min-interval", "0.05", "--retries", "1",
                "--max-requests", "3")
    assert r.returncode == 0, r.stderr
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["rows"][0]["status"] == "complete"
    assert manifest["rows"][0]["n_hits"] == 2
    norm = run_dir / "normalized" / "ecoli_smoke_t50.normalized.tsv"
    content = norm.read_text()
    assert content.splitlines()[0] == "acc\tani_estimation\tkmer_coverage\torganism"
    assert content.splitlines()[1].startswith("SRR111")

    # ledger: submission + request + checksums recorded; fetch is idempotent
    ledger = [json.loads(l) for l in
              (run_dir / "ledger.jsonl").read_text().splitlines() if l.strip()]
    events = [e["event"] for e in ledger]
    assert "submission_recorded" in events and "http_request" in events
    req = [e for e in ledger if e["event"] == "http_request"][0]
    assert len(req["sha256"]) == 64 and req["status_code"] == 200
    r = run_cli("fetch", "--run-dir", run_dir, "--base-url", server)
    assert r.returncode == 0
    assert "already complete" in r.stdout

    # cache integrity: payload tampering is detected
    cache_zip = run_dir / "cache" / f"{SESSION}.zip"
    cache_zip.write_bytes(cache_zip.read_bytes() + b"x")
    from logan_search_client import ResponseCache, CacheIntegrityError
    with pytest.raises(CacheIntegrityError):
        ResponseCache(str(run_dir / "cache")).get(SESSION)


def test_cli_fetch_pending_stays_submitted(tmp_path, server):
    fasta = tmp_path / "b.fa"
    fasta.write_text(f">b1\n{SEQ}\n")
    run_dir = tmp_path / "run"
    run_cli("prep", "--run-dir", run_dir, "--fasta", fasta,
            "--thresholds", "0.5")
    # point at a session the fake server does not know (400 => pending)
    run_cli("record", "--run-dir", run_dir, "--bait", "b1_t50",
            "--session", SESSION.replace("12345678", "87654321"))
    r = run_cli("fetch", "--run-dir", run_dir, "--base-url", server,
                "--poll", "2", "--min-interval", "0.05",
                "--max-requests", "3")
    assert r.returncode == 0
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["rows"][0]["status"] == "submitted"  # resume later
    ledger = [json.loads(l) for l in
              (run_dir / "ledger.jsonl").read_text().splitlines() if l.strip()]
    assert any(e["event"] == "poll_wait" for e in ledger)
    statuses = [e["status_code"] for e in ledger if e["event"] == "http_request"]
    assert statuses == [400, 400]
