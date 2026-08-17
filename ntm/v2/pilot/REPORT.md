# Execute bounded NTM Logan pilot — run report (2026-08-17)

**Outcome: documented NO-GO / preregistered 24-h pause. Stage-0 smoke gate NOT
passed; zero Stage-1 bait submissions spent.** No threshold fishing; no
results invented; all state resumable.

## 1. What was executed (all live actions, full provenance in ledgers)

### Stage 0 — smoke gate (mandatory before Stage 1)

Two dashboard submissions of the identical non-sensitive smoke query
(*E. coli* K-12 MG1655 `NC_000913.3:1,000,000–1,000,149`, 150 nt, seq sha256
`ac58a6648da0ca4ac7f14728e50f3e110bc3e1c25cdc52b544650ec2b4c801ed`,
group `GenBank_RefSeq`, threshold 0.5 = dashboard default):

| # | Time (UTC) | Submission path | Session id | Fate |
|---|---|---|---|---|
| 1 | ~20:42 | browser dashboard (mega-eval, Submit clicked, form verified: 1 query loaded, GenBank_RefSeq pill, threshold 0.5) | — | **id lost**: automation wrapper reset the tab to `about:blank` before the app recorded the id (localStorage `user-sessions` empty on next load). Counted as a lost-session submission. |
| 2 | 20:55:30 | same path | `kmviz-c112ba44-680f-4682-b2f5-696c9264e934` | **machine-captured** from the app's own `localStorage['user-sessions']` 3 s post-submit (no human transcription). Recorded via `run_pilot.py record` at 20:55:57Z. |

Retrieval of session 2 (`runs/stage0-smoke/ledger.jsonl`, 22 bounded,
rate-limited, checksummed requests; 20 via client + 1 manual final re-check + 1 dashboard session-load cross-check):

* 8 polls 20:56–21:04 (`--poll 8 --min-interval 60`), then 12 polls
  21:04–21:23 (`--poll 12 --min-interval 90`), then 1 manual re-check at
  ~21:31. **All HTTP 400**, body 46 bytes, sha256
  `980998f54f89f6f8d5b1ac0f2cedbe79a1c494d3abb529e288f96e6dbd31b04a`
  (byte-identical kmviz generic error = unknown-or-not-ready).
* Dashboard cross-check (session mode → Load, 21:03 and 21:24): app reports
  *"Session not found. Invalid session id, query still running, or results
  erased."* — i.e. the service itself says the query is still running at
  ~30 min, for a 150-nt query on the smallest sub-index.

### Stage 1 — prepared, NOT submitted (gate holds)

`ntm/v2/pilot/` holds the preregistered 24-submission panel
(gen/select_stage1.py, selection fixed before any Stage-1 result could be
seen):

* **A** — one representative `interior_module` bait per A-tier clade (15;
  clade 0_0021 → its `member_interior` bait), 1200 bp each.
* **B** — top-5 junction baits spanning **purely reconstructed** joins
  (`observed_adjacency_count == 0`), ranked co_occurrence DESC:
  0_0024_J2 (14), 0_0073_J1 (9), 0_0000_J2 (8), 0_0032_J2 (8), 0_0032_J3 (8).
* **C** — positive controls: mid-genome 1200 bp slices of complete public
  references **D29** (`AF022214.2`) and **L5** (`NC_001335.1`) from the
  curated external panel (checksums in `stage1_selection.tsv`).
* **D** — negative controls: `NTMBAIT_SHUF_CONTROL_01/02` (Markov
  composition-matched; expected zero 31-mer hits).

Group `GenBank_RefSeq`, threshold **0.7**, run dir
`logan/runs/stage1-gbrefseq/` (manifest + per-bait submission FASTAs +
INSTRUCTIONS.md, all rows `pending`). **Zero submissions made** — Stage-1 is
gated on the Stage-0 pass (PILOT_PLAN §4), which did not occur.

### Confirmation toolchain (built + validated, no live hits yet)

`ntm/v2/pilot/confirm_hits.py`: bounded S3 retrieval (HEAD size/ETag/Last-
Modified + sha256 ledger, ≤5 accessions/bait, 20 GB / 300-accession hard
stops), `back_to_sequences --output-mapping-positions` k-mer localization,
`minimap2 -c --eqx` (preset ladder sr → asm5 → asm20, cg:Z: CIGAR parsing,
query coverage + identity). Validated end-to-end on public probe
`DRR000016` contigs (189,571 B zst; 1,851 contigs): 150-nt probe sliced from
contig `DRR000016_60` → b2s matched exactly that contig with 120/120 k-mers
at positions 2000–2119; minimap2 `-x sr` = 150/150 aligned, identity 1.0.
minimap2 2.31-r1302 + back_to_sequences 0.8.4 installed at `/tmp/mmenv`
(micromamba, bioconda).

`ntm/v2/pilot/gen_submit_js.py`: generates the dashboard-submission
mega-eval (2-bit-packed ACGT payload; textarea → Load → group MultiSelect
toggle-safe select → threshold slider via synthetic Arrow keys (verified
0.5→0.7→restore) → Submit → in-page session-id capture from
`localStorage['user-sessions']`). Slider/group/pill state all verified live
without extra submissions.

## 2. Why this is a NO-GO pause, not a failure of the method

Preregistered rule (PILOT_PLAN §4, "Service degradation: > 30 min/query mean
latency … → pause ≥ 24 h, resume from manifest"): the smallest-group 150-nt
smoke query had not produced a retrievable result after 35 min and 21
checksummed polls. Both the documented download API and the dashboard's own
session loader agree the query is "still running". Because the session id
this time was **machine-captured from the app's own storage** (not
human-copied), the earlier mis-transcription hypothesis for the 2026-08-17
morning smoke (logan/smoke/run-2026-08-17/OUTCOME.md) is **not supported**
for this run: the id is exactly what the app recorded. This is a
service-latency / retrievability issue on the public instance.

Actions taken per preregistration:

* **No threshold fishing** — the smoke stayed at the default 0.5 throughout.
* **No Stage-1 spend** — bait budget preserved (0 of 64 ceiling used;
  2 smoke submissions total, both E. coli, non-sensitive).
* **No retry storms** — 21 polls total, ≥60 s spacing, all logged.
* Pause ≥ 24 h; results live one month; `runs/stage0-smoke/manifest.json`
  row stays `submitted` → `run_pilot.py fetch` resumes idempotently.
* Recommended escalation for the service maintainers (logan-search.org
  contact, docs "Submitting a query"): a GenBank_RefSeq-group 150-nt query
  (session `kmviz-c112ba44-680f-4682-b2f5-696c9264e934`, submitted
  2026-08-17T20:55:30Z) never became retrievable via
  `GET /api/download/<session>` (HTTP 400, byte-identical 46-byte body) or
  the dashboard session loader ("query still running") after 35 min,
  contradicting the documented "a few minutes" latency. The identical
  signature was seen for a human-submitted session earlier the same day
  (`kmviz-9b877a85-…`, 14 polls over 3 h).

## 3. Classification-ladder status (no analogue claims)

| Bait | Status |
|---|---|
| ecoli smoke | submitted; result not retrievable within preregistered window |
| all 24 Stage-1 rows (15 clade reps, 5 reconstructed junctions, D29/L5 positives, 2 shuffled negatives) | prepared, pending — **no screen result exists** |

No bait has any Logan-derived claim. Ladder for when results exist:
`no hit` → `accession-level screen hit` → `module-only homology` →
`supported junction/synteny` → `near-complete analogue` (PILOT_PLAN §5.5),
with same-accession (preferably same-contig) multi-bait/junction support
required before any "module/junction present" statement.

## 4. Resume procedure (after the ≥ 24-h pause)

1. `cd logan && python3 run_pilot.py fetch --run-dir runs/stage0-smoke --poll 6 --min-interval 60`
   — one bounded re-check of `kmviz-c112ba44-…`.
2. If HTTP 200 + parseable ZIP: gate passes → submit the 24 Stage-1 rows via
   the dashboard (`ntm/v2/pilot/gen_submit_js.py` per row,
   `logan/runs/stage1-gbrefseq/submissions/*.fa`, ≥2 min apart), record each
   session id, then `fetch` within 24 h.
3. If still 400 after the pause: escalate to maintainers before any further
   spend (do not diagnose on pilot budget).
4. On Stage-1 hits: `ntm/v2/pilot/confirm_hits.py --hits <normalized tsv>
   --baits-fa ntm/v2/pilot/stage1_panel.fa --out-dir confirm/` (budgets and
   provenance enforced by the script).

## 5. Files

| Path | Purpose |
|---|---|
| `logan/runs/stage0-smoke/` | manifest (row `submitted`), ledger.jsonl (21 requests + poll waits), rate-limit state |
| `logan/runs/stage1-gbrefseq/` | prepared Stage-1 run dir (24 pending rows, submissions/, INSTRUCTIONS.md) |
| `ntm/v2/pilot/select_stage1.py` | preregistered deterministic Stage-1 selection (reruns identical) |
| `ntm/v2/pilot/stage1_panel.fa` + `stage1_selection.tsv` + `stage1_selection_meta.json` | frozen panel with per-bait sha256, roles, sources |
| `ntm/v2/pilot/gen_submit_js.py` | dashboard submission driver (base64 mega-eval generator) |
| `ntm/v2/pilot/confirm_hits.py` | bounded S3 + b2s + minimap2 confirmation pipeline (validated) |
