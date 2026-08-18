# Pre-wake verification (attempt-1-2, 2026-08-18T08:10Z)

This file was written by the pre-wake incarnation of `resume-logan-stage-0`
(agent-67, attempt-1-2). The task is parked by WG on a timer until
**2026-08-18T21:00:29Z** (the preregistered >=24 h NO-GO pause from
PILOT_PLAN §4 / REPORT.md §6), so no network request was made in this
incarnation and none should be made before the timer fires.

## What this incarnation verified (all offline, no bait spend)

| Check | Result |
|---|---|
| `logan/tests/test_evaluate_gates.py` | 6/6 passed (2026-08-18T08:09Z) |
| `run_pilot.py status --run-dir runs/stage0-smoke` | OK; smoke row still `submitted` (idempotent resume state preserved) |
| `evaluate_gates.py --run-dir ../logan/runs/stage1-gbrefseq` (repo root) | exit 0; verdict INCOMPLETE (expected pre-run); G4 ledger/manifest reconcile PASS; gates.json regenerated |
| Stage-1 submissions | 24/24 `.fa` files present in `logan/runs/stage1-gbrefseq/submissions/`, matching manifest rows (threshold 0.7, group GenBank_RefSeq) |
| `gen_submit_js.py <fa> GenBank_RefSeq 0.7` | generates 4328-byte JS payload offline (base64 mega-eval) |
| `confirm_hits.py --help` | CLI intact (argparse, budgets) |
| `resume_stage0_recheck.sh` | sleep math targets 2026-08-18T21:00:30Z; at wake (~21:00:29Z) it fires the single bounded fetch immediately |

## Reminders for the woken agent (post-2026-08-18T21:00:29Z)

1. Exactly ONE bounded Stage-0 recheck — `cd logan && python3 run_pilot.py
   fetch --run-dir runs/stage0-smoke --poll 6 --min-interval 60`
   (or run `./resume_stage0_recheck.sh`, which does the same after its sleep).
2. HTTP 200 + parseable ZIP -> Stage-0 gate PASSES -> submit the 24 frozen
   rows per `logan/runs/stage1-gbrefseq/INSTRUCTIONS.md` (>=2 min apart,
   record every session id via `run_pilot.py record`), fetch within 24 h,
   then `ntm/v2/pilot/evaluate_gates.py`. Honor G1/G2/G3 verbatim; any
   failure = documented NO-GO, no threshold fishing.
3. Still 400 -> do NOT diagnose; post `ntm/v2/pilot/ESCALATION_DRAFT.md`
   content to the logan-search.org maintainers with session id
   `kmviz-c112ba44-680f-4682-b2f5-696c9264e934` + timestamps, and record the
   escalation in REPORT.md.
4. Stage-1 hits -> `ntm/v2/pilot/confirm_hits.py` (<=5 accessions/bait,
   20 GB / 300-acc stops); report per the 5-level ladder in REPORT.md §3;
   no analogue claims from scores alone.

## Provenance note

Graph writes (`wg log`/`msg`/`done`) were correctly fenced off for this
incarnation after the dispatcher parked the attempt (`attempt-parked`,
`wait:...:timer:46433s`, revision 11). This file + the accompanying commit are
the filesystem-level breadcrumbs; no daemon-internal state was modified.
