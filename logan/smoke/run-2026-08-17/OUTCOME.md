# Smoke run outcome — `smoke-ecoli-k12-150bp` (2026-08-17)

**Status: retrieval blocked.** No live results were obtained, and none are
inferred. This file records the precise blocker so `execute-bounded-ntm` can
re-run the smoke gate (§"Handoff") before any pilot submission.

## Submission (succeeded)

* Bait `ecoli_k12_mg1655_1000000_1000149_t50`: 150 nt of *E. coli* K-12
  MG1655 `NC_000913.3:1,000,000–1,000,149` (public BSL-1 reference;
  non-sensitive). Sequence sha256
  `ac58a6648da0ca4ac7f14728e50f3e110bc3e1c25cdc52b544650ec2b4c801ed`;
  120/120 k-mers (k=31) without N.
* Group `GenBank_RefSeq` (smallest sub-index, ~45k samples), threshold 0.5
  (dashboard default), no email — exactly one dashboard submission at
  ~2026-08-17T17:08Z, human-driven per `INSTRUCTIONS.md`.
* Session id `kmviz-9b877a85-62a8-46da-94d7-b056c3fda36d` copied by a human
  from the dashboard result-page URL and recorded offline via
  `run_pilot.py record` (ledger event `submission_recorded`).

## Retrieval (blocked)

`GET https://logan-search.org/api/download/kmviz-9b877a85-62a8-46da-94d7-b056c3fda36d`
— **14 bounded requests between 2026-08-17T17:08:23Z and 2026-08-17T20:17:37Z**
(13 during the original bounded polling windows, spaced 75–120 s; 1 final
check on task resume ~3 h later). Every response:

* HTTP **400**, body **46 bytes**, sha256
  `980998f54f89f6f8d5b1ac0f2cedbe79a1c494d3abb529e288f96e6dbd31b04a` —
  byte-identical across all 14 responses.
* 46 bytes is exactly the length of the kmviz generic API error string
  `An error occured while processing your request` (kmviz `api.py`; the same
  signature observed 2026-08-17 with a deliberately bogus session id — see
  `logan/README.md` §1). The kmviz download API uses this one body both for
  *unknown* and *not-yet-ready* sessions, so the response alone cannot
  distinguish the two.

## Interpretation (what is and is not established)

Established:

* The download endpoint is registered on the public instance and responding
  (kmviz-specific error; not the Dash 405 catch-all that unregistered routes
  return).
* The rate limiter, ledger, checksums, poll bounds, and resume semantics all
  worked as designed under live conditions (enforced waits 48.4 s and 78.5 s
  are recorded in `ledger.jsonl`; the resume check honored the persistent
  `rate_limit_state.json`).

NOT established — do not treat the session as expired-and-lost or the API as
broken. Docs (docs.logan-search.org, retrieved 2026-08-17) promise results in
"a few minutes" and retention of one month; a 150-nt query against the
smallest group should not take > 3 h. A 400 with the unknown-session
signature persisting > 3 h is most consistent with the **recorded session id
not being recognized by the API**. Plausible causes, in order:

1. **Mis-transcription** of the session id when copied from the dashboard URL
   (human step; single point of failure, no machine verification possible
   offline).
2. Server-side failure of the submission after the result page was shown,
   with the session purged from the API store.
3. The download API resolving sessions only in some state this session never
   reached (not documented; cannot be distinguished from (1)/(2) remotely).

## Actions deliberately NOT taken (bounds respected)

* No duplicate/re-submission of the smoke query (task + retry constraints).
* No full-panel or multi-threshold search; no threshold fishing.
* No automated retry storms: max 6 requests per fetch invocation, ≥ 30 s
  spacing, bounded backoff — all visible in `ledger.jsonl`.

## What the smoke did and did not demonstrate

| Requirement | Demonstrated? | Evidence |
|---|---|---|
| Documented interfaces verified w/ retrieval dates | yes | `logan/README.md` §1–3 (2026-08-17) |
| Bounded, rate-limited, checksummed, logged requests (live) | yes | `ledger.jsonl` (14 requests, params + sha256 + waits) |
| Cross-process/persistent rate limiting | yes | waits 48.4 s / 78.5 s enforced from `rate_limit_state.json` |
| Resumable run state | yes | manifest row still `submitted`; resume check at 20:17Z touched network only for that row |
| Response parsing / caching / deterministic normalization | **offline only** | 33 unit + CLI end-to-end tests vs local fake kmviz server (fixtures; no network) |
| Live 200 result parse | **no — blocked** | this file |
| Over-limit refusal | offline (unit tests) | `tests/test_logan_client.py` |

## Handoff to `execute-bounded-ntm` — Stage 0 (smoke gate, re-run)

The PILOT_PLAN smoke gate is **not** satisfied by this run. Before any
Stage-1 submission:

1. Submit ONE fresh non-sensitive smoke query (same E. coli slice, group
   `GenBank_RefSeq`, threshold 0.5) via the dashboard.
2. **Wait until the dashboard result page actually shows results** (table
   rendered), THEN copy the session id from the URL — verify it
   character-by-character against the URL bar before recording.
3. `run_pilot.py record`, then `fetch --poll 6 --min-interval 60` (≤ 7
   requests). Gate: HTTP 200, parseable ZIP, cache + normalized TSV written,
   ledger reconciles.
4. Only on gate pass proceed to Stage 1. If the fresh smoke also 400s with
   the same signature after a dashboard-verified result page, escalate the
   "download API vs dashboard state" discrepancy to the maintainers
   (logan-search.org contact) — do not spend pilot budget diagnosing.
