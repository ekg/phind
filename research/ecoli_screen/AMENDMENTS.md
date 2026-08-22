# AMENDMENTS — local-e-coli run

Amendments are logged per RUN_PLAN §0. None of these change thresholds,
budgets, bait panel, tier definitions, or the hit definition.

## A1 — 2026-08-22T00:5xZ — esummary switched to HTTP POST after HTTP 414

**Observed:** the metadata sweep (chunk=400 uids, GET query string) crashed
at uid ~9,971,xxx (row 141,076) with `HTTP Error 414: Request-URI Too Long`
— once uids reached 7–8 digits, 400 ids exceeded the GET URL limit. The
same chunk size worked in the NTM run only because its uid range was
shorter (≤7 digits, ~21k uids total).

**Remediation (transport-only, semantics unchanged):** `esummary.fcgi` is
now called with HTTP POST (body `db=sra&id=<csv>&retmode=json`) — the
documented EUtils mechanism for large id lists. Response bodies, parsing,
rate limit (≤ 1 req/s), retry/backoff, and the per-response sha256 receipts
are byte-identical in intent to the GET form. The sweep resumes from the
last persisted row (script is resumable by design; `done_uids` from the
output TSV).

No query string, filter, threshold, tier rule, or budget is affected. This
amendment is logged before the affected (re)start.
