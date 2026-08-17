#!/usr/bin/env python3
"""CLI driver for the bounded Logan Search pilot.

Subcommands
-----------
dryrun  Validate a bait FASTA + options; print exactly what would be requested.
        No network access at all.
prep    Build a resumable run directory: manifest.json, per-bait submission
        FASTA files (one sequence each, <= 2.5 kb), and INSTRUCTIONS.md for
        dashboard-assisted submission. No network access.
record  Record a dashboard session id (kmviz-<uuid>) for one bait row.
fetch   Fetch results for every `submitted` row through the documented
        GET /api/download/<session> endpoint: rate-limited, retried with
        backoff, checksummed, cached, parsed, normalized. Resumable.
status  Print manifest summary (pending/submitted/complete/failed counts).

Boundedness: every network-touching command enforces
  * hard cap on HTTP requests per invocation (--max-requests),
  * a minimum interval between requests (--min-interval, default 30 s),
  * one sequence per submission, <= 2.5 kb, threshold in [0.25, 1.0]
    (documented limits; see logan/README.md).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import List, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from logan_search_client import (  # noqa: E402
    DEFAULT_BASE_URL, GROUPS, INDEX_SNAPSHOT, K, MAX_QUERY_LEN,
    RESULTS_RETENTION, THRESHOLD_LADDER, THRESHOLD_MAX, THRESHOLD_MIN,
    BaitQuery, ClientConfig, HttpTransport, Ledger, LoganError,
    LoganSearchClient, Plan, RateLimiter, ResponseCache, build_manifest,
    dry_run_report, load_manifest, save_manifest, write_normalized_tsv,
)

INSTRUCTIONS_TEMPLATE = """# Logan Search dashboard submission instructions

Run: {run_id}
Generated: {generated_at}

For each file in submissions/ (in this order, one at a time):

1. Open https://logan-search.org/dashboard
2. "Submit a query": upload (or paste) the single-sequence FASTA.
3. Group: {group}   (recorded in the manifest row; do not change)
4. Threshold: {threshold}   (recorded in the manifest row; do not change)
5. Email: leave empty (no notifications; results retained {retention}).
6. Submit. The result page URL ends with the session id, e.g.
   https://logan-search.org/dashboard/kmviz-<uuid>  ->  session id
   kmviz-<uuid>. Copy it.
7. Record it WITHOUT any network request from this machine:
     python3 run_pilot.py record --run-dir . --bait <bait_id> \\
         --session kmviz-<uuid>
8. Only after all submissions are recorded (or in bounded batches):
     python3 run_pilot.py fetch --run-dir .

Rate-limit policy: at most one dashboard submission per {min_interval} s and
at most {max_requests} HTTP requests per fetch invocation. Do not run
multiple fetch processes concurrently against the same run directory.
"""


def read_bait_fasta(path: str, group: str, thresholds: List[float]) -> List[BaitQuery]:
    """Each FASTA record becomes len(thresholds) BaitQueries (sensitivity ladder)."""
    queries: List[BaitQuery] = []
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    _extend(queries, name, "".join(chunks), group, thresholds)
                name = line[1:].split()[0] if len(line) > 1 else "bait"
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        _extend(queries, name, "".join(chunks), group, thresholds)
    if not queries:
        raise LoganError(f"no FASTA records found in {path}")
    return queries


def _extend(queries: List[BaitQuery], name: str, seq: str,
            group: str, thresholds: List[float]) -> None:
    for t in thresholds:
        queries.append(
            BaitQuery(bait_id=f"{name}_t{int(round(t * 100))}",
                      sequence=seq, group=group, threshold=t,
                      note="sensitivity ladder")
        )


def cmd_dryrun(args: argparse.Namespace) -> int:
    queries = read_bait_fasta(args.fasta, args.group,
                              _ladder(args))
    plan = Plan(queries=queries, max_requests_per_run=args.max_requests)
    report = dry_run_report(plan, ClientConfig(
        base_url=args.base_url,
        min_request_interval_s=args.min_interval,
        max_requests_per_run=args.max_requests))
    json.dump(report, sys.stdout, indent=2, sort_keys=True)
    print()
    return 0


def cmd_prep(args: argparse.Namespace) -> int:
    queries = read_bait_fasta(args.fasta, args.group, _ladder(args))
    plan = Plan(queries=queries, max_requests_per_run=args.max_requests)
    plan.validate()  # refuses over-limit sizes/batches before anything is built
    os.makedirs(args.run_dir, exist_ok=True)
    subs = os.path.join(args.run_dir, "submissions")
    os.makedirs(subs, exist_ok=True)
    manifest = build_manifest(plan, args.run_id)
    for q in plan.queries:
        # manifest rows carry the sequence so `fetch`/audit can re-validate
        for row in manifest["rows"]:
            if row["bait_id"] == q.bait_id:
                row["sequence"] = q.sequence
        with open(os.path.join(subs, f"{q.bait_id}.fa"), "w") as fh:
            fh.write(f">{q.bait_id}\n")
            s = q.sequence
            for i in range(0, len(s), 70):
                fh.write(s[i : i + 70] + "\n")
    manifest["max_requests_per_run"] = args.max_requests
    manifest["group"] = args.group
    save_manifest(manifest, os.path.join(args.run_dir, "manifest.json"))
    with open(os.path.join(args.run_dir, "INSTRUCTIONS.md"), "w") as fh:
        fh.write(INSTRUCTIONS_TEMPLATE.format(
            run_id=args.run_id, generated_at=manifest["created_at"],
            group=args.group,
            threshold=" / ".join(f"{t:.2f}" for t in _ladder(args)),
            retention=RESULTS_RETENTION,
            min_interval=args.min_interval, max_requests=args.max_requests,
        ))
    print(f"prepared {len(plan.queries)} submissions in {args.run_dir}")
    print(f"next: submit each submissions/*.fa via the dashboard, then "
          f"`record` session ids, then `fetch`")
    return 0


def cmd_record(args: argparse.Namespace) -> int:
    mpath = os.path.join(args.run_dir, "manifest.json")
    manifest = load_manifest(mpath)
    hit = False
    for row in manifest["rows"]:
        if row["bait_id"] == args.bait:
            hit = True
            q = BaitQuery(bait_id=row["bait_id"], sequence=row["sequence"],
                          group=row["group"], threshold=row["threshold"])
            client = _offline_client(run_dir=args.run_dir)
            client.record_submission(q, args.session)
            row["session_id"] = args.session
            row["status"] = "submitted"
    if not hit:
        raise LoganError(f"bait {args.bait!r} not in manifest")
    save_manifest(manifest, mpath)
    print(f"recorded {args.bait} -> {args.session}")
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    mpath = os.path.join(args.run_dir, "manifest.json")
    manifest = load_manifest(mpath)
    client = _online_client(args)
    norm_dir = os.path.join(args.run_dir, "normalized")
    fetched = skipped = failed = 0
    for row in manifest["rows"]:
        if row["status"] != "submitted" or not row.get("session_id"):
            if row["status"] == "complete":
                skipped += 1
            continue
        params = {k: row.get(k) for k in (
            "bait_id", "group", "threshold", "seq_len", "kmers_total",
            "kmers_without_n", "seq_sha256")}
        attempts = args.poll or 1
        parsed = meta = None
        try:
            for attempt in range(attempts):
                parsed, meta = client.fetch_session(
                    row["session_id"], params=params,
                    expect_pending=bool(args.poll))
                if parsed is not None:
                    break
                if attempt < attempts - 1:
                    client.ledger.append(
                        "poll_wait", session=row["session_id"],
                        bait_id=row["bait_id"], attempt=attempt,
                        wait_s=args.min_interval)
                    time.sleep(args.min_interval)
            if parsed is None:
                # not ready after bounded polls; stays submitted for resume
                print(f"pending: {row['bait_id']} ({row['session_id']})")
                continue
            paths = write_normalized_tsv(parsed, norm_dir)
            summary = {
                "bait_id": row["bait_id"],
                "session": row["session_id"],
                "n_hits": {q: len(h) for q, h in parsed.hits_by_query.items()},
                "normalized_tsv": paths,
                "zip_sha256": meta.sha256,
                "zip_size_bytes": meta.size_bytes,
            }
            with open(os.path.join(norm_dir,
                                   f"{row['bait_id']}.summary.json"), "w") as fh:
                json.dump(summary, fh, indent=2, sort_keys=True)
            row["status"] = "complete"
            row["zip_sha256"] = meta.sha256
            row["n_hits"] = max(
                (len(h) for h in parsed.hits_by_query.values()), default=0)
            fetched += 1
            print(f"complete: {row['bait_id']} hits={row['n_hits']} "
                  f"sha256={meta.sha256[:12]}")
        except LoganError as exc:
            row["status"] = "failed"
            row["error"] = str(exc)
            failed += 1
            print(f"FAILED: {row['bait_id']}: {exc}", file=sys.stderr)
        save_manifest(manifest, mpath)
    print(f"fetch done: {fetched} complete, {skipped} already complete, "
          f"{failed} failed")
    return 1 if failed else 0


def cmd_status(args: argparse.Namespace) -> int:
    manifest = load_manifest(os.path.join(args.run_dir, "manifest.json"))
    counts: dict = {}
    for row in manifest["rows"]:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print(json.dumps({
        "run_id": manifest["run_id"],
        "created_at": manifest["created_at"],
        "index_snapshot": manifest["index_snapshot"],
        "status_counts": counts,
        "rows": [
            {"bait_id": r["bait_id"], "status": r["status"],
             "session_id": r.get("session_id"),
             "n_hits": r.get("n_hits"), "error": r.get("error")}
            for r in manifest["rows"]],
    }, indent=2, sort_keys=True))
    return 0


# ---------------------------------------------------------------------------


def _ladder(args: argparse.Namespace) -> List[float]:
    if args.thresholds:
        return [float(t) for t in args.thresholds]
    return list(THRESHOLD_LADDER)


def _offline_client(run_dir: str = ".") -> LoganSearchClient:
    cfg = ClientConfig()
    return LoganSearchClient(
        cfg, HttpTransport(base_url=cfg.base_url),
        ResponseCache(os.path.join(run_dir, "cache")),
        Ledger(os.path.join(run_dir, "ledger.jsonl")),
    )


def _online_client(args: argparse.Namespace) -> LoganSearchClient:
    cfg = ClientConfig(
        base_url=args.base_url,
        min_request_interval_s=args.min_interval,
        max_retries=args.retries,
        max_requests_per_run=args.max_requests,
    )
    return LoganSearchClient(
        cfg,
        HttpTransport(base_url=cfg.base_url, timeout_s=args.timeout),
        ResponseCache(os.path.join(args.run_dir, "cache")),
        Ledger(os.path.join(args.run_dir, "ledger.jsonl")),
        rate_limiter=RateLimiter(
            cfg.min_request_interval_s,
            state_path=os.path.join(args.run_dir, "rate_limit_state.json")),
    )


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--fasta", help="bait FASTA (one record per bait)")
    common.add_argument("--group", default="GenBank_RefSeq",
                        choices=sorted(GROUPS),
                        help="search group (default: GenBank_RefSeq, smallest)")
    common.add_argument("--thresholds", nargs="*", type=float,
                        default=None,
                        help="sensitivity ladder (default: 0.5 0.7 0.9; "
                             "bounds [0.25, 1.0] per docs)")
    common.add_argument("--max-requests", type=int, default=40,
                        help="hard cap on HTTP requests per invocation")
    common.add_argument("--min-interval", type=float, default=30.0,
                        help="min seconds between HTTP requests")
    common.add_argument("--base-url", default=DEFAULT_BASE_URL)
    common.add_argument("--timeout", type=float, default=60.0)
    common.add_argument("--retries", type=int, default=4)

    d = sub.add_parser("dryrun", parents=[common],
                       help="validate + plan; no network")
    d.set_defaults(func=cmd_dryrun)

    pr = sub.add_parser("prep", parents=[common],
                        help="build resumable run dir + submission packages")
    pr.add_argument("--run-dir", required=True)
    pr.add_argument("--run-id", default=None)
    pr.set_defaults(func=cmd_prep)

    rec = sub.add_parser("record", help="record a dashboard session id (offline)")
    rec.add_argument("--run-dir", required=True)
    rec.add_argument("--bait", required=True)
    rec.add_argument("--session", required=True,
                     help="kmviz-<uuid> from the dashboard result URL")
    rec.set_defaults(func=cmd_record)

    f = sub.add_parser("fetch", help="fetch submitted sessions (network)")
    f.add_argument("--run-dir", required=True)
    f.add_argument("--poll", type=int, default=None, metavar="MAX_TRIES",
                   help="treat HTTP 400 as 'session not ready' and allow up "
                        "to MAX_TRIES bounded polls spaced --min-interval")
    f.add_argument("--max-requests", type=int, default=40)
    f.add_argument("--min-interval", type=float, default=30.0)
    f.add_argument("--base-url", default=DEFAULT_BASE_URL)
    f.add_argument("--timeout", type=float, default=60.0)
    f.add_argument("--retries", type=int, default=4)
    f.set_defaults(func=cmd_fetch)

    st = sub.add_parser("status", help="manifest summary")
    st.add_argument("--run-dir", required=True)
    st.set_defaults(func=cmd_status)

    args = p.parse_args(argv)
    try:
        return args.func(args)
    except LoganError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
