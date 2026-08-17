#!/usr/bin/env python3
"""Unit tests for the bounded Logan Search pilot client.

All network-facing tests use fake transports and synthetic/fixture responses;
no test touches the live service. Run:  python3 -m pytest tests/ -q
"""

from __future__ import annotations

import io
import json
import os
import sys
import zipfile
from typing import List, Optional

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from logan_search_client import (  # noqa: E402
    CacheIntegrityError, ClientConfig, HttpTransport, Ledger,
    LoganRequestError, LoganSearchClient, LoganValidationError,
    MalformedResponseError, Plan, RateLimiter, ResponseCache, BaitQuery,
    build_manifest, dry_run_report, load_manifest, normalize_sequence,
    parse_results_zip, plan_from_manifest, save_manifest,
    write_normalized_tsv,
)

SEQ_150 = ("GTGTCAGCTTTCGTGGTGTGCAGCTGGCGTCAGATGACAACATGCTGCCAGACAGCCTGAAAGGGTTTGC"
           "GCCTGTGGTGCGTGGTATCGCCAAAAGCAATGCCCAGATAACGATTAAGCAAAATGGTTACACCATTTAC"
           "CAAACTTATG")


# ---------------------------------------------------------------------------
# fixtures / helpers
# ---------------------------------------------------------------------------


def make_zip(rows_by_query, session_json=None, include_readme=True) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for qname, rows in rows_by_query.items():
            cols = ["acc", "kmer_coverage", "ANI_estimation", "organism"]
            out = io.StringIO()
            out.write("\t".join(cols) + "\n")
            for r in rows:
                out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
            zf.writestr(f"{qname}.tsv", out.getvalue())
        if session_json is not None:
            zf.writestr("session.json", json.dumps(session_json))
        if include_readme:
            zf.writestr("README.md", "# results\n")
    return buf.getvalue()


class FakeResponse:
    def __init__(self, status_code: int, content: bytes):
        self.status_code = status_code
        self.content = content


class FakeTransport(HttpTransport):
    """Scripted transport: pops scripted responses; counts calls."""

    def __init__(self, scripted: List[FakeResponse],
                 exceptions: Optional[List[Exception]] = None):
        super().__init__(base_url="https://fake.example")
        self.scripted = list(scripted)
        self.exceptions = list(exceptions or [])
        self.calls: List[str] = []

    def get(self, url: str):
        self.calls.append(url)
        if self.exceptions:
            raise self.exceptions.pop(0)
        return self.scripted.pop(0)


class FakeClock:
    def __init__(self):
        self.now = 1_000_000.0

    def time(self):
        return self.now

    def sleep(self, s: float):
        self.now += s


def make_client(tmp_path, transport, clock=None, cfg=None) -> LoganSearchClient:
    cfg = cfg or ClientConfig(min_request_interval_s=0.001,
                              max_retries=2, backoff_base_s=0.001,
                              timeout_s=1.0)
    limiter = RateLimiter(cfg.min_request_interval_s,
                          state_path=str(tmp_path / "rl.json"),
                          sleep=(clock.sleep if clock else (lambda s: None)))
    if clock:
        import logan_search_client as L
        L.time.sleep = clock.sleep
        L.time.time = clock.time
    return LoganSearchClient(
        cfg, transport,
        ResponseCache(str(tmp_path / "cache")),
        Ledger(str(tmp_path / "ledger.jsonl")),
        rate_limiter=limiter,
        sleep=(clock.sleep if clock else (lambda s: None)),
    )


def read_ledger(tmp_path) -> List[dict]:
    with open(tmp_path / "ledger.jsonl") as fh:
        return [json.loads(line) for line in fh if line.strip()]


@pytest.fixture(autouse=True)
def restore_time():
    import logan_search_client as L
    real_time, real_sleep = L.time.time, L.time.sleep
    yield
    L.time.time, L.time.sleep = real_time, real_sleep


# ---------------------------------------------------------------------------
# query validation (documented limits enforced, over-limit refused)
# ---------------------------------------------------------------------------


class TestValidation:
    def test_accepts_documented_query(self):
        q = BaitQuery(bait_id="b1", sequence=SEQ_150,
                      group="GenBank_RefSeq", threshold=0.5).validate()
        assert len(q.sequence) == 150
        assert q.kmers() == {"len": 150, "kmers_total": 120,
                             "kmers_without_n": 120}

    def test_refuses_over_limit_sequence(self):
        with pytest.raises(LoganValidationError, match="too long"):
            BaitQuery(bait_id="b1", sequence="ACGT" * 626,  # 2504 nt > 2500
                      group="All", threshold=0.5).validate()

    def test_refuses_short_sequence(self):
        with pytest.raises(LoganValidationError, match="too short"):
            BaitQuery(bait_id="b1", sequence="ACGT" * 7,  # 28 nt < 31
                      group="All", threshold=0.5).validate()

    def test_refuses_multi_sequence_submission(self):
        with pytest.raises(LoganValidationError, match="single bare sequence"):
            BaitQuery(bait_id="b1", sequence=f"{SEQ_150}\n>other\nACGT",
                      group="All", threshold=0.5).validate()

    def test_refuses_non_dna_characters(self):
        with pytest.raises(LoganValidationError, match="non-ACGTN"):
            BaitQuery(bait_id="b1", sequence=SEQ_150 + "U",
                      group="All", threshold=0.5).validate()

    def test_refuses_threshold_out_of_bounds(self):
        for bad in (0.1, 0.24, 1.01, 2.0):
            with pytest.raises(LoganValidationError, match="threshold"):
                BaitQuery(bait_id="b1", sequence=SEQ_150,
                          group="All", threshold=bad).validate()

    def test_refuses_unknown_group(self):
        with pytest.raises(LoganValidationError, match="unknown group"):
            BaitQuery(bait_id="b1", sequence=SEQ_150,
                      group="Everything", threshold=0.5).validate()

    def test_refuses_oversized_plan_batch(self):
        queries = [q(f"b{i:03d}") for i in range(41)]
        plan = Plan(queries=queries, max_requests_per_run=40)
        with pytest.raises(LoganValidationError, match="hard cap"):
            plan.validate()

    def test_normalization_is_deterministic(self):
        raw = f">hdr\n{SEQ_150[:70]}\n{SEQ_150[70:]}\n"
        a = normalize_sequence(raw)
        b = normalize_sequence(SEQ_150.lower())
        assert a == b == SEQ_150


# ---------------------------------------------------------------------------
# response parsing / normalization
# ---------------------------------------------------------------------------


class TestParsing:
    def test_parses_and_sorts_hits(self):
        z = make_zip({"b1": [
            {"acc": "SRR2", "kmer_coverage": 0.42, "organism": "x"},
            {"acc": "SRR1", "kmer_coverage": 0.99, "organism": "y"},
            {"acc": "SRR3", "kmer_coverage": 0.42, "organism": "z"},
        ]})
        parsed = parse_results_zip(z, session="kmviz-x")
        rows = parsed.hits_by_query["b1"]
        assert [r["acc"] for r in rows] == ["SRR1", "SRR2", "SRR3"]  # cov desc, acc asc
        assert rows[0]["kmer_coverage"] == 0.99

    def test_empty_hits_is_valid(self):
        z = make_zip({"b1": []})
        parsed = parse_results_zip(z)
        assert parsed.hits_by_query["b1"] == []

    def test_missing_required_columns_is_malformed(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr("b1.tsv", "acc\torganism\nSRR1\tx\n")
        with pytest.raises(MalformedResponseError, match="kmer_coverage"):
            parse_results_zip(buf.getvalue())

    def test_non_zip_is_malformed(self):
        with pytest.raises(MalformedResponseError, match="not a ZIP"):
            parse_results_zip(b"<html>error page</html>")

    def test_zip_without_tsv_is_malformed(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr("README.md", "hi")
        with pytest.raises(MalformedResponseError, match="no per-query"):
            parse_results_zip(buf.getvalue())

    def test_normalized_tsv_is_deterministic(self, tmp_path):
        rows = [
            {"acc": "SRR2", "kmer_coverage": 0.42, "organism": "x"},
            {"acc": "SRR1", "kmer_coverage": 0.99, "organism": "y"},
        ]
        d1, d2 = tmp_path / "a", tmp_path / "b"
        p1 = write_normalized_tsv(parse_results_zip(make_zip({"b1": rows})), str(d1))
        p2 = write_normalized_tsv(parse_results_zip(make_zip({"b1": rows})), str(d2))
        assert open(p1["b1"]).read() == open(p2["b1"]).read()
        content = open(p1["b1"]).read()
        assert content.splitlines()[0] == "acc\tani_estimation\tkmer_coverage\torganism"
        assert "SRR1" in content.splitlines()[1]

    def test_session_json_captured(self):
        z = make_zip({"b1": []}, session_json={"kmviz-x": {}})
        parsed = parse_results_zip(z)
        assert parsed.session_json is not None
        assert "README.md" in parsed.files


# ---------------------------------------------------------------------------
# rate limiting
# ---------------------------------------------------------------------------


class TestRateLimiter:
    def test_enforces_spacing(self, tmp_path):
        clock = FakeClock()
        rl = RateLimiter(30.0, state_path=str(tmp_path / "rl.json"),
                         sleep=clock.sleep)
        import logan_search_client as L
        real_time = L.time.time
        L.time.time = clock.time
        try:
            rl.wait()                      # first request: no wait
            waited = rl.wait()             # second: must space 30 s
        finally:
            L.time.time = real_time
        assert waited == pytest.approx(30.0)

    def test_state_persists_across_instances(self, tmp_path):
        clock = FakeClock()
        import logan_search_client as L
        real_time = L.time.time
        L.time.time = clock.time
        try:
            rl1 = RateLimiter(10.0, state_path=str(tmp_path / "rl.json"),
                              sleep=clock.sleep)
            rl1.wait()
            clock.now += 1.0
            rl2 = RateLimiter(10.0, state_path=str(tmp_path / "rl.json"),
                              sleep=clock.sleep)
            waited = rl2.wait()
        finally:
            L.time.time = real_time
        assert waited == pytest.approx(9.0)


# ---------------------------------------------------------------------------
# fetch: retries, backoff, caching, ledger
# ---------------------------------------------------------------------------


class TestFetch:
    def test_retries_then_succeeds_and_records_backoff(self, tmp_path):
        z = make_zip({"b1": [{"acc": "SRR1", "kmer_coverage": 0.9}]})
        t = FakeTransport([
            FakeResponse(500, b"boom"),
            FakeResponse(500, b"boom"),
            FakeResponse(200, z),
        ])
        c = make_client(tmp_path, t)
        parsed, meta = c.fetch_session("kmviz-aaaa1111-bbbb-cccc-dddd-eeeeeeeeeeee")
        assert t.calls  # three attempts
        assert parsed.hits_by_query["b1"][0]["acc"] == "SRR1"
        events = [e["event"] for e in read_ledger(tmp_path)]
        assert events.count("http_request") == 3
        assert events.count("backoff") == 2

    def test_retry_exhaustion_raises_and_logs(self, tmp_path):
        t = FakeTransport([], exceptions=[__import__("requests").ConnectionError(
            "refused")] * 3)
        c = make_client(tmp_path, t)
        with pytest.raises(LoganRequestError, match="after 3 attempts"):
            c.fetch_session("kmviz-aaaa1111-bbbb-cccc-dddd-eeeeeeeeeeee")
        events = [e["event"] for e in read_ledger(tmp_path)]
        assert events.count("http_error") == 3

    def test_cache_prevents_second_download(self, tmp_path):
        z = make_zip({"b1": []})
        t = FakeTransport([FakeResponse(200, z)])
        c = make_client(tmp_path, t)
        sid = "kmviz-aaaa1111-bbbb-cccc-dddd-eeeeeeeeeeee"
        c.fetch_session(sid)
        n_calls = len(t.calls)
        c.fetch_session(sid)  # served from cache
        assert len(t.calls) == n_calls
        assert any(e["event"] == "cache_hit" for e in read_ledger(tmp_path))

    def test_cache_integrity_checked(self, tmp_path):
        z = make_zip({"b1": []})
        t = FakeTransport([FakeResponse(200, z)])
        c = make_client(tmp_path, t)
        sid = "kmviz-aaaa1111-bbbb-cccc-dddd-eeeeeeeeeeee"
        c.fetch_session(sid)
        payload = tmp_path / "cache" / f"{sid}.zip"
        payload.write_bytes(payload.read_bytes() + b"tampered")
        with pytest.raises(CacheIntegrityError):
            c.cache.get(sid)

    def test_hard_request_cap_refused(self, tmp_path):
        z = make_zip({"b1": []})
        t = FakeTransport([FakeResponse(200, z)] * 10)
        cfg = ClientConfig(min_request_interval_s=0.001, max_retries=0,
                           backoff_base_s=0.001, max_requests_per_run=2)
        c = make_client(tmp_path, t, cfg=cfg)
        for i in range(2):
            sid = f"kmviz-aaaa111{i}-bbbb-cccc-dddd-eeeeeeeeeeee"
            c.fetch_session(sid)
        with pytest.raises(Exception, match="hard request cap"):
            c.fetch_session("kmviz-aaaa1119-bbbb-cccc-dddd-eeeeeeeeeeee")

    def test_pending_400_is_soft_miss(self, tmp_path):
        t = FakeTransport([FakeResponse(400, b"An error occured")])
        c = make_client(tmp_path, t)
        parsed, meta = c.fetch_session(
            "kmviz-aaaa1111-bbbb-cccc-dddd-eeeeeeeeeeee", expect_pending=True)
        assert parsed is None and meta is None

    def test_ledger_records_checksum_and_params(self, tmp_path):
        z = make_zip({"b1": [{"acc": "SRR9", "kmer_coverage": 1.0}]})
        t = FakeTransport([FakeResponse(200, z)])
        c = make_client(tmp_path, t)
        sid = "kmviz-aaaa1111-bbbb-cccc-dddd-eeeeeeeeeeee"
        c.fetch_session(sid, params={"bait_id": "b1", "threshold": 0.5})
        entry = [e for e in read_ledger(tmp_path)
                 if e["event"] == "http_request"][0]
        assert entry["sha256"] and len(entry["sha256"]) == 64
        assert entry["params"]["threshold"] == 0.5
        assert entry["status_code"] == 200


# ---------------------------------------------------------------------------
# manifest / resume behavior
# ---------------------------------------------------------------------------


def q(bait_id="b1", seq=SEQ_150, group="GenBank_RefSeq", thr=0.5):
    return BaitQuery(bait_id=bait_id, sequence=seq, group=group, threshold=thr)


class TestResume:
    def test_manifest_roundtrip_and_plan_rebuild(self, tmp_path):
        plan = Plan(queries=[q("b1"), q("b2")])
        m = build_manifest(plan, "run-x")
        for row, query in zip(m["rows"], plan.queries):
            row["sequence"] = query.sequence
        path = str(tmp_path / "manifest.json")
        save_manifest(m, path)
        m2 = load_manifest(path)
        plan2 = plan_from_manifest(m2)
        assert [x.bait_id for x in plan2.queries] == ["b1", "b2"]

    def test_fetch_skips_complete_and_pending_rows(self, tmp_path):
        m = build_manifest(Plan(queries=[q("b1"), q("b2"), q("b3")]), "run-y")
        z = make_zip({"b2": [{"acc": "SRR1", "kmer_coverage": 0.8}]})
        t = FakeTransport([FakeResponse(200, z)])
        c = make_client(tmp_path, t)
        for row in m["rows"]:
            row["sequence"] = SEQ_150
            row["session_id"] = f"kmviz-aaaa111{row['bait_id'][-1]}-bbbb-cccc-dddd-eeeeeeeeeeee"
        m["rows"][0]["status"] = "complete"   # already done -> skip network
        m["rows"][1]["status"] = "submitted"  # will be fetched
        m["rows"][2]["status"] = "pending"    # no session -> untouched
        path = str(tmp_path / "manifest.json")
        save_manifest(m, path)
        # emulate cmd_fetch loop on the RELOADED manifest (resume semantics)
        m2 = load_manifest(path)
        for row in m2["rows"]:
            if row["status"] != "submitted" or not row.get("session_id"):
                continue
            parsed, meta = c.fetch_session(row["session_id"],
                                           params={"bait_id": row["bait_id"]})
            write_normalized_tsv(parsed, str(tmp_path / "norm"))
            row["status"] = "complete"
        save_manifest(m2, path)
        final = load_manifest(path)
        assert [r["status"] for r in final["rows"]] == [
            "complete", "complete", "pending"]
        assert len(t.calls) == 1  # only the submitted row hit the network

    def test_failed_row_is_marked_not_lost(self, tmp_path):
        t = FakeTransport([], exceptions=[__import__("requests").Timeout(
            "t")] * 3)
        c = make_client(tmp_path, t)
        m = build_manifest(Plan(queries=[q("b1")]), "run-z")
        m["rows"][0]["sequence"] = SEQ_150
        m["rows"][0]["session_id"] = "kmviz-aaaa1111-bbbb-cccc-dddd-eeeeeeeeeeee"
        m["rows"][0]["status"] = "submitted"
        try:
            c.fetch_session(m["rows"][0]["session_id"])
            m["rows"][0]["status"] = "complete"
        except LoganRequestError as exc:
            m["rows"][0]["status"] = "failed"
            m["rows"][0]["error"] = str(exc)
        assert m["rows"][0]["status"] == "failed"
        assert "after 3 attempts" in m["rows"][0]["error"]


# ---------------------------------------------------------------------------
# dry run
# ---------------------------------------------------------------------------


class TestDryRun:
    def test_dry_run_reports_plan_without_network(self, tmp_path):
        plan = Plan(queries=[q("b1"), q("b1_t70"[:3] and "b1_t70", thr=0.7)])
        report = dry_run_report(plan)
        assert report["dry_run"] is True
        assert report["planned_download_requests"] == 2
        assert report["threshold_bounds"] == [0.25, 1.0]
        assert report["max_query_len_nt"] == 2500
        assert not os.path.exists(str(tmp_path / "ledger.jsonl"))


# ---------------------------------------------------------------------------
# real smoke fixture (present only after the bounded live smoke run; skipped
# otherwise so tests never depend on network)
# ---------------------------------------------------------------------------

REAL_FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "fixtures", "smoke_ecoli_k12.zip")


class TestRealSmokeFixture:
    @pytest.mark.skipif(not os.path.exists(REAL_FIXTURE),
                        reason="live smoke response not yet cached")
    def test_real_zip_parses_and_normalizes(self, tmp_path):
        data = open(REAL_FIXTURE, "rb").read()
        parsed = parse_results_zip(data, session="smoke-fixture")
        assert parsed.hits_by_query, "real zip must contain at least one query tsv"
        for qname, rows in parsed.hits_by_query.items():
            covs = [r["kmer_coverage"] for r in rows]
            assert covs == sorted(covs, reverse=True)
            for r in rows:
                assert r["acc"]
        paths = write_normalized_tsv(parsed, str(tmp_path))
        again = write_normalized_tsv(parse_results_zip(data), str(tmp_path / "2"))
        for q in paths:
            assert open(paths[q]).read() == open(again[q]).read()


# ---------------------------------------------------------------------------
# session id format
# ---------------------------------------------------------------------------


class TestSubmissionRecording:
    def test_rejects_malformed_session_id(self, tmp_path):
        t = FakeTransport([])
        c = make_client(tmp_path, t)
        with pytest.raises(LoganValidationError, match="kmviz-"):
            c.record_submission(q("b1"), "not-a-session-id")

    def test_records_valid_session(self, tmp_path):
        t = FakeTransport([])
        c = make_client(tmp_path, t)
        sid = "kmviz-08a28c6c-9691-4eca-a8aa-ef62f098f62c"
        c.record_submission(q("b1"), sid)
        e = [e for e in read_ledger(tmp_path)
             if e["event"] == "submission_recorded"][0]
        assert e["session"] == sid
        assert e["threshold"] == 0.5
        assert e["seq_sha256"]
