# Escalation draft — Logan Search session not retrievable (Stage-0 smoke)

Status: DRAFT — to be posted as a GitHub issue on
`IndexThePlanet/LoganSearch` **only if** the preregistered >=24 h re-check
(2026-08-18T21:00Z) still returns HTTP 400 for the session below.
Channel chosen because the repo README ("Submit an issue") asks users to
report query problems there *including their session-id*.

---

**Title:** Session `kmviz-c112ba44-…` never retrievable via `/api/download` ("query still running" after 24+ h)

**Body:**

Hi — reporting a query whose results never became downloadable, per the
README's request to include the session id.

- **Session id:** `kmviz-c112ba44-680f-4682-b2f5-696c9264e934`
- **Submitted:** 2026-08-17 20:55:30 UTC via the dashboard
  (https://logan-search.org/dashboard), single 150-nt query
  (*E. coli* K-12 MG1655 `NC_000913.3:1,000,000–1,000,149`), group
  `GenBank_RefSeq`, threshold 0.5 (default), no email supplied.
- The session id was captured from the app's own
  `localStorage['user-sessions']` a few seconds after submit, so it is not a
  transcription error.
- **Download API:** `GET https://logan-search.org/api/download/kmviz-c112ba44-680f-4682-b2f5-696c9264e934`
  returned **HTTP 400** with a 46-byte error body
  (sha256 `980998f54f89f6f8d5b1ac0f2cedbe79a1c494d3abb529e288f96e6dbd31b04a`)
  on every bounded poll from 2026-08-17 20:56 UTC through 21:31 UTC
  (21 polls, >= 60 s apart), and again on
  **[RECHECK_TS — fill from ledger final_status entry]** after the
  recommended wait.
- **Dashboard session loader** (Load session by id), checked 2026-08-17
  21:03 and 21:24 UTC: *"Session not found. Invalid session id, query still
  running, or results erased."*
- An earlier human-submitted session the same day
  (`kmviz-9b877a85-…`, 14 polls over ~3 h) showed the identical signature.

Questions:

1. Is `GET /api/download/<session>` expected to serve dashboard-submitted
   sessions on the public instance, or is the dashboard result page the only
   retrieval surface?
2. If the query is genuinely still queued after 24+ h on the smallest group
   (`GenBank_RefSeq`, 150-nt query), is there a way to check queue position
   or expected latency beyond "a few minutes"?
3. Does "results retained one month" start at submission or at completion?

Happy to run any additional read-only check you suggest. Thanks!

---

## Poster checklist (do not post before the recheck)

- [ ] Re-check executed at >= 2026-08-18T21:00Z (`logan/resume_stage0_recheck.sh`)
- [ ] Fill `[RECHECK_TS]` and the recheck outcome line from
      `logan/runs/stage0-smoke/ledger.jsonl` (`final_status` entry)
- [ ] If the recheck returned HTTP 200 instead: do NOT post; proceed to
      Stage-1 per `ntm/v2/pilot/REPORT.md` section 4.
- [ ] Post to https://github.com/IndexThePlanet/LoganSearch/issues (new issue),
      then record the issue URL in `ntm/v2/pilot/REPORT.md` and in the WG
      task log. No further pilot spend until maintainers respond.
