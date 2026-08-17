#!/usr/bin/env python3
"""Safe, bounded, cached, rate-limited client for Logan Search (IndexThePlanet).

Scope of this module (verified against official documentation on 2026-08-17,
see logan/README.md for full citations and retrieval dates):

* Logan Search (https://logan-search.org) is a kmviz (https://tlemane.github.io/kmviz)
  front-end over a kmindex k=31 index of Logan unitigs (SRA assemblies).
* Query submission is DASHBOARD-ONLY on the public instance: the kmviz REST
  query routes (POST /api/query, POST /api/query/<db>) are not registered
  (verified live 2026-08-17: POST /api/query -> 405 from the Dash catch-all).
* The kmviz results-download route IS registered on the public instance:
  GET <base>/api/download/<session> returns a ZIP with one TSV per query
  sequence plus session.json and README.md (verified live 2026-08-17 against
  the kmviz source at github.com/tlemane/kmviz, kmviz/api.py).
* Documented query constraints (docs.logan-search.org):
    - ONE FASTA sequence per submission, max length 2.5 kb
    - search group in {All, All_No_viral_human, Fast, Fast_No_human, GenBank_RefSeq}
    - threshold = minimum proportion of query k-mers shared, in [0.25, 1.0]
    - results retained one month; query IDs look like kmviz-<uuid>

This client therefore:
  - validates every query against the documented limits BEFORE any request;
  - refuses over-limit sequences / batches (hard bound, configurable);
  - rate-limits every HTTP request with a persistent, cross-run spacer;
  - caches raw responses content-addressed (sha256) and never re-downloads;
  - records every request/response (params, timestamps, status, checksum,
    retries, backoff) in an append-only JSONL ledger;
  - parses + deterministically normalizes the kmviz download ZIP;
  - supports dry-run, resume, and a threshold sensitivity ladder;
  - performs NO unbounded work: a hard cap on requests per run is enforced.

It does NOT scrape the dashboard: submissions are prepared as operator-ready
packages (see run_pilot.py prep) or performed by an explicit browser-assisted
step, and only the documented download endpoint is used programmatically.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import time
import zipfile
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, Optional, Tuple

import requests

# ---------------------------------------------------------------------------
# Documented service facts (see logan/README.md; retrieved 2026-08-17)
# ---------------------------------------------------------------------------

DEFAULT_BASE_URL = "https://logan-search.org"
DOWNLOAD_ROUTE_FMT = "/api/download/{session}"

K = 31  # k-mer size used by the Logan Search index (docs.logan-search.org)

#: Maximum query length: one FASTA sequence, max 2.5 kb (docs.logan-search.org).
MAX_QUERY_LEN = 2500
#: Minimum sensible query length: at least one k-mer (k = 31).
MIN_QUERY_LEN = K

#: Threshold bounds: minimum proportion of query k-mers, in [0.25, 1.0]
#: (docs.logan-search.org, "Threshold (optional)").
THRESHOLD_MIN = 0.25
THRESHOLD_MAX = 1.0

#: Candidate sensitivity ladder for the NTM bait pilot (all within bounds).
THRESHOLD_LADDER = (0.5, 0.7, 0.9)

#: Documented search groups with coverage notes (docs.logan-search.org,
#: retrieved 2026-08-17). Percentages are of the 23.4M-sample search index.
#: NOTE: the live dashboard combobox (observed 2026-08-17) additionally
#: exposes Fast_No_RefSeq, Transcriptomic, Metatranscriptomic, Metagenomic;
#: these are marked documented=False until the docs catch up.
GROUPS: Dict[str, Dict[str, Any]] = {
    "All": {
        "coverage": "all Logan unitigs, SRA up until 2023 - 23.4M total samples",
        "pct_samples": "100%",
        "documented": True,
    },
    "All_No_viral_human": {
        "coverage": "excludes viral and human samples",
        "pct_samples": "65.49%",
        "documented": True,
    },
    "Fast": {
        "coverage": "optimized subset for faster queries, excludes viral samples; "
                    "~99.5% of non-viral samples",
        "pct_samples": "77.59%",
        "documented": True,
    },
    "Fast_No_human": {
        "coverage": "same as Fast but excluding human samples",
        "pct_samples": "64.99%",
        "documented": True,
    },
    "GenBank_RefSeq": {
        "coverage": "only GenBank and RefSeq reference genomes",
        "pct_samples": "~45k samples",
        "documented": True,
    },
    "Fast_No_RefSeq": {
        "coverage": "live dashboard option (2026-08-17); not in docs",
        "pct_samples": "unknown",
        "documented": False,
    },
    "Transcriptomic": {
        "coverage": "live dashboard option (2026-08-17); not in docs",
        "pct_samples": "unknown",
        "documented": False,
    },
    "Metatranscriptomic": {
        "coverage": "live dashboard option (2026-08-17); not in docs",
        "pct_samples": "unknown",
        "documented": False,
    },
    "Metagenomic": {
        "coverage": "live dashboard option (2026-08-17); not in docs",
        "pct_samples": "unknown",
        "documented": False,
    },
}

#: Index snapshot of record for the search service (docs.logan-search.org,
#: retrieved 2026-08-17). NOTE: the Logan S3 dataset itself has newer releases
#: (v1.2 = December 2025 SRA freeze); the *search index* documented above
#: still describes SRA up until 2023. Record both; do not conflate.
INDEX_SNAPSHOT = {
    "search_index": "kmindex k=31 index over Logan unitigs (Bloom-filter based), "
                    "SRA up until 2023, 23.4M accessions + GenBank/RefSeq references",
    "dataset_latest": "Logan v1.2 (SRA freeze 2025-12-31): 38.1M accessions with "
                      "unitigs, 37.3M with contigs (IndexThePlanet/Logan Stats-v1.2.md)",
    "verified": "2026-08-17",
    "sources": [
        "https://docs.logan-search.org/",
        "https://github.com/IndexThePlanet/Logan",
        "https://github.com/IndexThePlanet/LoganSearch",
    ],
}

RESULTS_RETENTION = "one month (docs.logan-search.org)"

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class LoganError(Exception):
    """Base error for the Logan pilot client."""


class LoganValidationError(LoganError):
    """A query/plan violates documented limits; no request was made."""


class MalformedResponseError(LoganError):
    """A response could not be parsed (missing columns, corrupt zip, ...)."""


class CacheIntegrityError(LoganError):
    """Cached payload does not match its recorded checksum."""


class LoganRequestError(LoganError):
    """HTTP request failed after the configured number of retries."""


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def normalize_sequence(raw: str) -> str:
    """Deterministically normalize a raw FASTA sequence body.

    Uppercases, removes whitespace, '*' and '-' padding, and FASTA header
    lines. Rejects (LoganValidationError) anything outside A/C/G/T/N after
    normalization, because the kmindex k-mer alphabet is DNA.
    """
    seq = "".join(
        line.strip() for line in raw.splitlines() if not line.startswith(">")
    )
    seq = re.sub(r"[\s\*\-]", "", seq).upper()
    bad = sorted(set(seq) - set("ACGTN"))
    if bad:
        raise LoganValidationError(
            f"sequence contains non-ACGTN characters after normalization: {bad}"
        )
    return seq


def count_kmers(seq: str, k: int = K) -> Dict[str, int]:
    """Count (31-)mers and the subset without N, canonical-form agnostic.

    The engine counts canonical k-mers (kmtricks/kmindex); we only report
    counts for validation/provenance, not to reimplement the engine.
    """
    n = len(seq)
    total = max(0, n - k + 1)
    no_n = sum(1 for i in range(total) if "N" not in seq[i : i + k])
    return {"len": n, "kmers_total": total, "kmers_without_n": no_n}


# ---------------------------------------------------------------------------
# Rate limiter (persistent, cross-run)
# ---------------------------------------------------------------------------


class RateLimiter:
    """Enforce a minimum wall-clock interval between HTTP requests.

    State (last request wall time) is persisted so spacing is enforced across
    process restarts, protecting the public service during resumed pilots.
    """

    def __init__(self, min_interval_s: float, state_path: Optional[str] = None,
                 sleep: Callable[[float], None] = time.sleep):
        if min_interval_s <= 0:
            raise ValueError("min_interval_s must be > 0")
        self.min_interval_s = float(min_interval_s)
        self.state_path = state_path
        self._sleep = sleep

    def _last(self) -> Optional[float]:
        if self.state_path and os.path.exists(self.state_path):
            try:
                with open(self.state_path) as fh:
                    return float(json.load(fh)["last_request_ts"])
            except (ValueError, KeyError, OSError):
                return None
        return None

    def wait(self, now: Optional[float] = None) -> float:
        """Block until the next request may be issued; return the wait time."""
        now = time.time() if now is None else now
        waited = 0.0
        last = self._last()
        if last is not None:
            delta = now - last
            if delta < self.min_interval_s:
                waited = self.min_interval_s - delta
                self._sleep(waited)
        self.record()
        return waited

    def record(self) -> None:
        ts = time.time()
        if self.state_path:
            os.makedirs(os.path.dirname(self.state_path) or ".", exist_ok=True)
            tmp = self.state_path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump({"last_request_ts": ts}, fh)
            os.replace(tmp, self.state_path)


# ---------------------------------------------------------------------------
# Ledger (append-only JSONL, every request/attempt/failure)
# ---------------------------------------------------------------------------


class Ledger:
    def __init__(self, path: str):
        self.path = path
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    def append(self, event: str, **fields: Any) -> Dict[str, Any]:
        entry = {"ts": utc_now(), "event": event}
        entry.update(fields)
        with open(self.path, "a") as fh:
            fh.write(json.dumps(entry, sort_keys=True) + "\n")
        return entry


# ---------------------------------------------------------------------------
# Cache (content-addressed raw responses)
# ---------------------------------------------------------------------------


@dataclass
class CachedResponse:
    session: str
    url: str
    sha256: str
    size_bytes: int
    fetched_at: str
    params: Dict[str, Any]


class ResponseCache:
    """Content-addressed store of raw (bytes) responses keyed by session id.

    A cached entry is only trusted if the stored payload still hashes to the
    recorded sha256; otherwise CacheIntegrityError.
    """

    def __init__(self, root: str):
        self.root = root
        os.makedirs(root, exist_ok=True)

    def _payload_path(self, session: str) -> str:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", session)
        return os.path.join(self.root, f"{safe}.zip")

    def _meta_path(self, session: str) -> str:
        return self._payload_path(session) + ".meta.json"

    def has(self, session: str) -> bool:
        return os.path.exists(self._payload_path(session)) and os.path.exists(
            self._meta_path(session)
        )

    def get(self, session: str) -> Tuple[bytes, CachedResponse]:
        with open(self._payload_path(session), "rb") as fh:
            data = fh.read()
        with open(self._meta_path(session)) as fh:
            meta = json.load(fh)
        digest = sha256_bytes(data)
        if digest != meta["sha256"]:
            raise CacheIntegrityError(
                f"cached payload for {session} hashes to {digest}, "
                f"meta records {meta['sha256']}"
            )
        return data, CachedResponse(
            session=session,
            url=meta["url"],
            sha256=digest,
            size_bytes=len(data),
            fetched_at=meta["fetched_at"],
            params=meta.get("params", {}),
        )

    def put(self, session: str, data: bytes, url: str, params: Dict[str, Any]) -> CachedResponse:
        with open(self._payload_path(session), "wb") as fh:
            fh.write(data)
        meta = CachedResponse(
            session=session,
            url=url,
            sha256=sha256_bytes(data),
            size_bytes=len(data),
            fetched_at=utc_now(),
            params=params,
        )
        with open(self._meta_path(session), "w") as fh:
            json.dump(asdict(meta), fh, sort_keys=True, indent=2)
        return meta


# ---------------------------------------------------------------------------
# Query records + validation
# ---------------------------------------------------------------------------


@dataclass
class BaitQuery:
    """One planned Logan Search submission (one sequence, one threshold)."""

    bait_id: str
    sequence: str
    group: str
    threshold: float
    note: str = ""

    def validate(self) -> "BaitQuery":
        if not self.bait_id or not re.fullmatch(r"[A-Za-z0-9_.\-]+", self.bait_id):
            raise LoganValidationError(f"invalid bait_id: {self.bait_id!r}")
        if self.group not in GROUPS:
            raise LoganValidationError(
                f"unknown group {self.group!r}; documented groups: {sorted(GROUPS)}"
            )
        if not (THRESHOLD_MIN <= self.threshold <= THRESHOLD_MAX):
            raise LoganValidationError(
                f"threshold {self.threshold} outside documented range "
                f"[{THRESHOLD_MIN}, {THRESHOLD_MAX}]"
            )
        seq = normalize_sequence(self.sequence)
        if "\n" in self.sequence or ">" in self.sequence:
            raise LoganValidationError(
                "submission must be a single bare sequence (no FASTA header, "
                "one sequence per submission per docs.logan-search.org)"
            )
        if len(seq) < MIN_QUERY_LEN:
            raise LoganValidationError(
                f"sequence too short: {len(seq)} nt < k={K}"
            )
        if len(seq) > MAX_QUERY_LEN:
            raise LoganValidationError(
                f"sequence too long: {len(seq)} nt > documented max "
                f"{MAX_QUERY_LEN} nt (2.5 kb, one sequence per submission)"
            )
        self.sequence = seq
        self.threshold = round(float(self.threshold), 2)
        return self

    def kmers(self) -> Dict[str, int]:
        return count_kmers(self.sequence)

    def request_params(self) -> Dict[str, Any]:
        """Everything that determines/records this request (for the ledger)."""
        km = self.kmers()
        return {
            "bait_id": self.bait_id,
            "group": self.group,
            "threshold": self.threshold,
            "seq_len": km["len"],
            "kmers_total": km["kmers_total"],
            "kmers_without_n": km["kmers_without_n"],
            "seq_sha256": sha256_bytes(self.sequence.encode()),
            "k": K,
            "engine": "kmindex/kmviz dashboard submission + /api/download retrieval",
        }


@dataclass
class Plan:
    """A bounded set of BaitQueries with a hard request cap."""

    queries: List[BaitQuery] = field(default_factory=list)
    max_requests_per_run: int = 40

    def validate(self) -> "Plan":
        seen = set()
        for q in self.queries:
            q.validate()
            key = (q.bait_id, q.group, q.threshold)
            if key in seen:
                raise LoganValidationError(f"duplicate query key: {key}")
            seen.add(key)
        if len(self.queries) > self.max_requests_per_run:
            raise LoganValidationError(
                f"plan has {len(self.queries)} queries > hard cap "
                f"{self.max_requests_per_run} per run (bounded pilot)"
            )
        return self


# ---------------------------------------------------------------------------
# Response parsing / deterministic normalization
# ---------------------------------------------------------------------------

#: Columns expected in every Logan-Search per-query result TSV. The metadata
#: schema is documented at docs.logan-search.org (Table view); `acc` and
#: `kmer_coverage` are the minimum required for a well-formed hit table.
REQUIRED_TSV_COLUMNS = ("acc", "kmer_coverage")

#: Normalization order (hit rows sorted by these keys, descending coverage).
_SORT_KEYS = (("kmer_coverage", True), ("acc", False))

_FLOAT_COLS = {
    "kmer_coverage", "ani_estimation", "p_value", "e_value",
    "contigs_n50", "contigs_nbseq", "contigs_maxlen", "contigs_sumlen",
    "unitigs_n50", "unitigs_nbseq", "unitigs_maxlen", "unitigs_sumlen",
    "mbytes", "avgspotlen", "mbases",
}


def _coerce(col: str, value: str) -> Any:
    if value == "" or value is None:
        return None
    if col in _FLOAT_COLS:
        try:
            return float(value)
        except ValueError:
            return value
    return value


@dataclass
class ParsedSession:
    session: str
    hits_by_query: Dict[str, List[Dict[str, Any]]]
    session_json: Optional[Dict[str, Any]] = None
    files: List[str] = field(default_factory=list)


def parse_results_zip(data: bytes, session: str = "unknown") -> ParsedSession:
    """Parse a kmviz download ZIP into deterministic, normalized hit tables.

    Rules (deterministic):
      * one TSV per query sequence; columns lower-cased; numeric-ish columns
        coerced to float; rows sorted by kmer_coverage DESC then acc ASC;
      * a TSV with a header and zero rows is a valid EMPTY result;
      * a TSV missing `acc` or `kmer_coverage`, or any non-ZIP payload, is
        malformed (MalformedResponseError) -- never silently invent hits.
    """
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise MalformedResponseError("payload is not a ZIP archive")
    hits: Dict[str, List[Dict[str, Any]]] = {}
    session_json: Optional[Dict[str, Any]] = None
    files: List[str] = []
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        for name in sorted(zf.namelist()):
            files.append(name)
            if name == "README.md":
                continue
            if name == "session.json":
                try:
                    session_json = json.loads(zf.read(name).decode("utf-8"))
                except (ValueError, UnicodeDecodeError) as exc:
                    raise MalformedResponseError(f"bad session.json: {exc}")
                continue
            if not name.endswith(".tsv"):
                continue
            query_name = name[: -len(".tsv")]
            text = zf.read(name).decode("utf-8", errors="replace")
            reader = csv.DictReader(io.StringIO(text), delimiter="\t")
            raw_cols = reader.fieldnames or []
            cols = [c.strip().lower() for c in raw_cols]
            missing = [c for c in REQUIRED_TSV_COLUMNS if c not in cols]
            if missing:
                raise MalformedResponseError(
                    f"{name}: missing required column(s) {missing}; "
                    f"found {cols}"
                )
            rows = []
            for row in reader:
                norm = {
                    (k or "").strip().lower(): _coerce((k or "").strip().lower(), v)
                    for k, v in row.items()
                    if k is not None
                }
                if norm.get("acc") in (None, ""):
                    continue  # blank padding row, not a hit
                rows.append(norm)
            rows.sort(
                key=lambda r: (
                    -(r.get("kmer_coverage") or 0.0)
                    if isinstance(r.get("kmer_coverage"), (int, float))
                    else 0.0,
                    str(r.get("acc") or ""),
                )
            )
            hits[query_name] = rows
    if not hits:
        raise MalformedResponseError(
            f"ZIP contains no per-query .tsv result files (files: {files})"
        )
    return ParsedSession(session=session, hits_by_query=hits,
                         session_json=session_json, files=files)


def write_normalized_tsv(parsed: ParsedSession, out_dir: str) -> Dict[str, str]:
    """Write one deterministic normalized TSV per query; returns paths."""
    os.makedirs(out_dir, exist_ok=True)
    out: Dict[str, str] = {}
    for query_name, rows in sorted(parsed.hits_by_query.items()):
        cols: List[str] = []
        for r in rows:
            for c in r:
                if c not in cols:
                    cols.append(c)
        cols = sorted(cols)  # stable schema independent of row order
        path = os.path.join(out_dir, f"{query_name}.normalized.tsv")
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(cols)
            for r in rows:
                w.writerow([r.get(c, "") for c in cols])
        out[query_name] = path
    return out


# ---------------------------------------------------------------------------
# Transport + client
# ---------------------------------------------------------------------------


class HttpTransport:
    """Real HTTP transport for the documented download endpoint only."""

    def __init__(self, base_url: str = DEFAULT_BASE_URL, timeout_s: float = 60.0,
                 session: Optional[requests.Session] = None,
                 user_agent: str = "phind-logan-pilot/0.1 (bounded; contact: erikg)"):
        self.base_url = base_url.rstrip("/")
        self.timeout_s = timeout_s
        self._session = session or requests.Session()
        self._session.headers["User-Agent"] = user_agent

    def download_url(self, session_id: str) -> str:
        return self.base_url + DOWNLOAD_ROUTE_FMT.format(session=session_id)

    def get(self, url: str) -> requests.Response:
        return self._session.get(url, timeout=self.timeout_s)


@dataclass
class ClientConfig:
    base_url: str = DEFAULT_BASE_URL
    min_request_interval_s: float = 30.0   # >= 2 req/min to the public service
    max_retries: int = 4
    backoff_base_s: float = 5.0            # 5, 10, 20, 40 (deterministic)
    backoff_cap_s: float = 300.0
    timeout_s: float = 60.0
    max_requests_per_run: int = 40
    max_sessions_per_run: int = 40


class LoganSearchClient:
    """Bounded fetch/parse/cache driver for Logan Search session results.

    NOTE on submission: the public Logan Search instance exposes no REST
    query endpoint (verified 2026-08-17, see README). Submissions happen via
    the dashboard (human or browser-assisted) and are recorded here as
    `record_submission(...)`; this client then fetches, checksums, caches,
    parses and normalizes results from the documented download endpoint.
    """

    def __init__(self, cfg: ClientConfig, transport: HttpTransport,
                 cache: ResponseCache, ledger: Ledger,
                 rate_limiter: Optional[RateLimiter] = None,
                 sleep: Callable[[float], None] = time.sleep):
        self.cfg = cfg
        self.transport = transport
        self.cache = cache
        self.ledger = ledger
        self.limiter = rate_limiter or RateLimiter(cfg.min_request_interval_s)
        self._sleep = sleep
        self._requests_this_run = 0

    # -- submission bookkeeping (no network) --------------------------------

    def record_submission(self, query: BaitQuery, session_id: str,
                          submitted_at: Optional[str] = None) -> None:
        query.validate()
        if not re.fullmatch(r"kmviz-[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-"
                            r"[0-9a-f]{4}-[0-9a-f]{12}", session_id):
            raise LoganValidationError(
                f"session id {session_id!r} does not match documented "
                "kmviz-<uuid> format"
            )
        self.ledger.append("submission_recorded", session=session_id,
                           submitted_at=submitted_at or utc_now(),
                           **query.request_params())

    # -- bounded fetch with retries/backoff ----------------------------------

    def _bound_check(self) -> None:
        if self._requests_this_run >= self.cfg.max_requests_per_run:
            raise LoganError(
                f"hard request cap for this run reached "
                f"({self.cfg.max_requests_per_run}); refusing to continue "
                "(bounded pilot)"
            )

    def fetch_session(self, session_id: str, params: Optional[Dict[str, Any]] = None,
                      expect_pending: bool = False) -> Tuple[ParsedSession, CachedResponse]:
        """Fetch (or load from cache) one session's results.

        Retries transient failures with deterministic exponential backoff and
        records every attempt in the ledger. `expect_pending=True` treats the
        documented 400 "session not ready/unknown" response as a soft miss
        (returns None) so callers can poll gently.
        """
        params = dict(params or {})
        if self.cache.has(session_id):
            data, meta = self.cache.get(session_id)
            self.ledger.append("cache_hit", session=session_id,
                               sha256=meta.sha256, size_bytes=meta.size_bytes,
                               params=params)
            return parse_results_zip(data, session=session_id), meta

        url = self.transport.download_url(session_id)
        last_error: Optional[str] = None
        for attempt in range(self.cfg.max_retries + 1):
            self._bound_check()
            self._requests_this_run += 1
            waited = self.limiter.wait()
            t0 = time.time()
            try:
                resp = self.transport.get(url)
                dt_ms = int((time.time() - t0) * 1000)
                body = resp.content
                self.ledger.append(
                    "http_request",
                    session=session_id, url=url, attempt=attempt,
                    status_code=resp.status_code,
                    size_bytes=len(body), sha256=sha256_bytes(body),
                    duration_ms=dt_ms, rate_limit_wait_s=round(waited, 3),
                    params=params,
                )
                if resp.status_code == 200:
                    meta = self.cache.put(session_id, body, url, params)
                    return parse_results_zip(body, session=session_id), meta
                if resp.status_code == 400 and expect_pending:
                    # documented kmviz API error body for unknown/not-ready
                    # sessions; treat as soft "not ready yet"
                    return None, None  # type: ignore[return-value]
                last_error = f"HTTP {resp.status_code}: {body[:200]!r}"
            except requests.RequestException as exc:
                dt_ms = int((time.time() - t0) * 1000)
                last_error = f"{type(exc).__name__}: {exc}"
                self.ledger.append(
                    "http_error", session=session_id, url=url, attempt=attempt,
                    error=last_error, duration_ms=dt_ms, params=params,
                )
            if attempt < self.cfg.max_retries:
                backoff = min(self.cfg.backoff_base_s * (2 ** attempt),
                              self.cfg.backoff_cap_s)
                self.ledger.append("backoff", session=session_id,
                                   attempt=attempt, backoff_s=backoff)
                self._sleep(backoff)
        raise LoganRequestError(
            f"session {session_id}: failed after {self.cfg.max_retries + 1} "
            f"attempts; last error: {last_error}"
        )

    # -- dry run --------------------------------------------------------------

def dry_run_report(plan: Plan, cfg: Optional["ClientConfig"] = None) -> Dict[str, Any]:
    """Validate a plan and describe exactly what WOULD be requested (pure)."""
    plan.validate()
    submissions = [q.request_params() for q in plan.queries]
    return {
        "dry_run": True,
        "generated_at": utc_now(),
        "base_url": (cfg.base_url if cfg else DEFAULT_BASE_URL),
        "download_route": DOWNLOAD_ROUTE_FMT,
        "index_snapshot": INDEX_SNAPSHOT,
        "groups_available": {g: GROUPS[g] for g in sorted(GROUPS)},
        "threshold_bounds": [THRESHOLD_MIN, THRESHOLD_MAX],
        "max_query_len_nt": MAX_QUERY_LEN,
        "k": K,
        "results_retention": RESULTS_RETENTION,
        "planned_submissions": submissions,
        "planned_download_requests": len(plan.queries),
        "max_requests_per_run": plan.max_requests_per_run,
        "min_request_interval_s": (cfg.min_request_interval_s if cfg else 30.0),
        "note": "submission is dashboard-assisted; only /api/download is "
                "called programmatically, one request per session id",
    }


# ---------------------------------------------------------------------------
# Resume manifest
# ---------------------------------------------------------------------------

#: Manifest row statuses.
STATUS_PENDING = "pending"            # validated, not yet submitted
STATUS_SUBMITTED = "submitted"        # dashboard session id recorded
STATUS_COMPLETE = "complete"          # results fetched, parsed, normalized
STATUS_FAILED = "failed"              # retries exhausted; see ledger


def load_manifest(path: str) -> Dict[str, Any]:
    with open(path) as fh:
        return json.load(fh)


def save_manifest(manifest: Dict[str, Any], path: str) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(manifest, fh, sort_keys=True, indent=2)
    os.replace(tmp, path)


def manifest_row(query: BaitQuery, status: str = STATUS_PENDING,
                 session_id: Optional[str] = None,
                 error: Optional[str] = None) -> Dict[str, Any]:
    row = dict(query.request_params())
    row.update({"status": status, "session_id": session_id, "error": error})
    return row


def build_manifest(plan: Plan, run_id: str) -> Dict[str, Any]:
    plan.validate()
    return {
        "run_id": run_id,
        "created_at": utc_now(),
        "base_url": DEFAULT_BASE_URL,
        "index_snapshot": INDEX_SNAPSHOT,
        "results_retention": RESULTS_RETENTION,
        "rows": [manifest_row(q) for q in plan.queries],
    }


def plan_from_manifest(manifest: Dict[str, Any]) -> Plan:
    """Rebuild a validated Plan from a manifest (sequence via seq_sha256 check
    happens for rows that still carry one; fetch-only rows skip it)."""
    queries = []
    for row in manifest["rows"]:
        if "sequence" not in row:
            continue  # fetch-only manifest (e.g., recorded sessions)
        queries.append(
            BaitQuery(
                bait_id=row["bait_id"],
                sequence=row["sequence"],
                group=row["group"],
                threshold=row["threshold"],
                note=row.get("note", ""),
            ).validate()
        )
    return Plan(queries=queries,
                max_requests_per_run=manifest.get("max_requests_per_run", 40))
