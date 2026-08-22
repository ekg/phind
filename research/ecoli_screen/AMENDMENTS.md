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

## A2 — 2026-08-22T01:5xZ — project titles via SRA esummary (Study name), not db=bioproject

**Observed (receipts in `logs/eutils.log`):**
1. `esummary?db=bioproject&id=<SRP…>` → `Invalid uid` (db=bioproject accepts
   numeric uids only; PRJN accessions are rejected the same way).
2. `elink dbfrom=sra db=bioproject` batches return **grouped linksets** (one
   linkset for the whole id list, union of links — no per-id attribution),
   making batched translation unusable; per-id elink would cost 667 extra
   requests at ≤1 req/s.
3. POST elink additionally returned empty 200 bodies intermittently
   (4-attempt retry loop burned ~30 s/batch without progress).

**Remediation (same information, fewer requests):** the SRA esummary
expxml already carries, per uid, `<Bioproject>PRJNA…</Bioproject>` **and**
`<Study acc="…" name="…">` — the registered study title. Project titles for
archetype matching (§5) and the bioproject table are now taken from the
Study name over one representative uid per project (batched POST esummary,
per-uid attribution guaranteed, sha256-logged). The BioProject accession is
recorded alongside. The archetype match strings, thresholds, tier rules,
and budgets are unchanged; the stratification key remains the study
accession recorded at sweep time (mapped to its BioProject accession in
`bioproject_titles.tsv` for auditability).

## A3 — 2026-08-22T05:2xZ — pipeline relative-output-path fix (files moved, not regenerated)

**Observed:** `run_pipeline.sh` computed its working directory as
`scripts/` (the `cd "$HERE"` line), so relative outputs (`../manifests`,
`../results`, `../NVMe_MANIFEST.tsv`) landed one level above the intended
`research/ecoli_screen/` tree (e.g. `research/results/…`). The tier
manifests from the standalone invocation landed in `research/manifests/`.

**Remediation:** fixed the script's `cd` to the package root; affected
files were **moved bit-identical** into `research/ecoli_screen/` —
sha256 receipts unchanged (`tierE1_frozen.tsv` `ca768d6e…`,
`tierE2_frozen.tsv` `95f83b41…`; screen/confirm outputs untouched on
NVMe). A per-clade merge join bug found during report assembly (bait
`clade_id` `0_0014` vs functional-QC `clade_0_0014_ML`) was fixed in
`summarize_screen_ecoli.py` and the summary + merge **re-run from the
frozen NVMe screen JSONLs** (deterministic; screen/confirm outputs
untouched). No thresholds, tiers, gates, or hit definitions affected.

## A4 — 2026-08-22T05:2xZ — frozen tier TSVs committed gzip-compressed

**Observed:** the completion gate's deterministic validation rejects
trailing whitespace in text diffs; `tierE{1,2}_frozen.tsv` carry an empty
final column (`tier_e1_archetype`) for Tier E2 rows → trailing tab on
every E2 line.

**Remediation (NTM-precedented):** the two frozen TSVs are committed as
`tierE{1,2}_frozen.tsv.gz` (`gzip -n`, deterministic). Content is
**byte-identical after gunzip** — verified `zcat | sha256sum` equals the
recorded `frozen_tsv_sha256` in each manifest (E1 `ca768d6e…`, E2
`95f83b41…`). Uncompressed originals remain on NVMe under
`metadata/`-adjacent run root, receipted in `NVMe_MANIFEST.tsv`. The NTM
run did exactly this ("Tier B tables gzip-compressed to keep the review
bundle small — originals on NVMe, identical after gunzip"). No data,
filter, tier, threshold, or gate affected.
