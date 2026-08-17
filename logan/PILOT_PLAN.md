# Bounded NTM Logan bait pilot — execution plan (for `execute-bounded-ntm`)

Author: `implement-safe-logan` (2026-08-17). All service facts cited in
`logan/README.md` (retrieval date 2026-08-17). Nothing here authorizes an
unbounded or full-panel search.

## 0. Preregistered principle

Logan Search results are **accession-level k-mer screens only**. No claim of a
natural analogue may rest on a Logan score alone; every promotion requires S3
sequence retrieval + `back_to_sequences` + minimap2 confirmation
(§5). Thresholds are fixed in advance; **no threshold fishing** — if controls
fail, record NO-GO.

## 1. Inputs (from upstream tasks)

| Input | Source task | Form |
|---|---|---|
| ~20-candidate bait panel (0.5–2.5 kb interior/module + junction baits) | `design-ntm-prophage` | FASTA + manifest (bait id, class, coordinates, 31-mer stats, checksum) |
| Public mycobacteriophage references (positive controls) | `curate-public-mycobacteriophage` | reference FASTA slices ≤ 2.5 kb |
| Negative controls | `design-ntm-prophage` (shuffled/low-complexity-masked) | FASTA |
| Smoke (done) | this task | E. coli K-12 150 nt, GenBank_RefSeq, thr 0.5 |

Every bait must pass `logan_search_client` validation (single sequence,
ACGTN, 31–2500 nt) before submission; over-limit baits are rejected and
reported back to `design-ntm-prophage`, never truncated silently.

## 2. Query budget (hard bounds)

Sub-groups ("Groups" in the dashboard) and thresholds:

* **Stage 1 — reference screen.** All baits + positive + negative controls,
  group `GenBank_RefSeq` (~45k samples, smallest index), threshold **0.7**
  (primary). ≈ 20 + 2 + 2 = **24 submissions**.
* **Stage 2 — sensitivity ladder (conditional on Stage-1 gates).** Only baits
  with zero or ambiguous GenBank_RefSeq hits: re-run at **0.5** and **0.9**
  (documented range [0.25, 1.0]); optionally the same ladder on `Fast`
  (77.59% of samples) for at most 10 baits. ≤ 3 × 10 = **30 submissions**.
* **Stage 3 — `Fast` group expansion (conditional).** ≤ 10 best-supported
  baits, thresholds fixed from Stage 2 outcomes. ≤ **10 submissions**.

**Total ceiling: 64 dashboard submissions**, each one sequence ≤ 2.5 kb.
Pacing: ≥ 1 submission / 2 min (manual/dashboard-assisted, one at a time);
downloads only via the documented `GET /api/download/<session>` with
`--min-interval 60 --max-requests 40` per invocation; results expire after
**one month** — fetch within 24 h of submission.

Time estimate: docs promise results "in a few minutes" per query; the smoke
query (GenBank_RefSeq) was Pending > 6 min, so budget **5–15 min/query** ⇒
Stage 1 ≈ 2–6 h wall clock (submission-paced), Stage 2+3 similar if triggered.

## 3. Download volume budget (S3 `logan-pub`)

* Result ZIPs: negligible (KBs each).
* Per-hit accession retrieval for confirmation: **contigs first**
  (`/c/<acc>/<acc>.contigs.fa.zst`; 623 TB total across 37.3M accessions ⇒
  median far below 1 GB; observed examples: DRR000016 contigs = 190 KB;
  SRR17555654 unitigs = 2.3 GB — unitigs can be huge, download only when
  contig-level confirmation is insufficient).
* Bound: **top ≤ 5 hit accessions per bait** (ranked by kmer_coverage, ties
  by e-value), and only baits that passed Stage-1 gates ⇒ ≤ 20 × 5 = 100
  accessions; expected median bacterial WGS contigs ~1–10 MB compressed ⇒
  **~0.1–1 GB total**; hard stop at **20 GB** or 300 accessions, whichever
  first. Every download records URL, size, ETag/sha256.
* Gate: if the median hit-accession contig file exceeds **50 MB**, STOP and
  re-scope (likely metagenomic accessions; consider `back_to_sequences`
  streaming via HTTP range or restrict to fewer baits).

## 4. Stop / Go gates (preregistered)

**Smoke gate (already passed):** submission → session id → documented download
endpoint → parse/normalize/cache/checksum all work; ledger reconciles.

**Stage 1 → Stage 2 GO requires ALL of:**
1. ≥ 1 positive control (public mycobacteriophage reference slice) hits its
   own reference genome accession with `kmer_coverage ≥ 0.7` at threshold 0.7
   in `GenBank_RefSeq`.
2. Negative controls produce **0 hits at threshold 0.9** in `GenBank_RefSeq`.
3. Smoke query and ≥ 90% of Stage-1 sessions fetched and parsed without
   malformed-response failures.
4. Ledger/manifest reconciliation: every submitted bait has session id,
   status, zip sha256; no orphan sessions.

**NO-GO / STOP conditions (any one):**
* Negative controls hit at 0.9 (index contamination or semantics mismatch) →
  NO-GO, document, do not shop thresholds.
* Positive control fails to hit its own reference at 0.7 → investigate
  sub-index/strand semantics before any Stage-2 spend; max 2 diagnostic
  resubmissions.
* Service degradation: > 30 min/query mean latency, repeated 5xx, or > 20%
  session loss → pause ≥ 24 h, resume from manifest (results live 1 month).
* Any download-budget bound (§3) hit → STOP, report.

## 5. Confirmation workflow (per bounded hit set)

**Local toolchain status (2026-08-17):** `zstd`, `aws` CLI, `samtools`,
`pyarrow` installed; `back_to_sequences` and `minimap2` NOT yet installed.
Install before Stage-1 confirmation work:
`micromamba install -c bioconda back_to_sequences minimap2` (or
`cargo install --git https://github.com/pierrepeterlongo/back_to_sequences`
for b2s). `micromamba` and `cargo` are present.

For each bait with screen hits (≤ 5 accessions):

1. **Retrieve** contigs from `https://s3.amazonaws.com/logan-pub/c/<acc>/<acc>.contigs.fa.zst`
   (fallback unitigs `/u/...` only if contig confirmation is ambiguous);
   checksum (sha256) + record S3 ETag/size/Last-Modified.
2. **k-mer localization:** `back_to_sequences --in-kmers <bait>.fa
   --in-sequences <acc>.contigs.fa.zst --out-sequences <out>.txt
   --output-mapping-positions` (bioconda `back_to_sequences`; docs:
   b2s-doc.readthedocs.io). Yields matching contig ids + bait k-mer
   positions/orientations.
3. **Alignment confirmation:** extract matched contigs; `minimap2 -c
   (--eqx)` bait vs contigs (preset `asm5` for close matches, `asm20`/`sr`
   fallback); require: query coverage and identity consistent with the
   reported `kmer_coverage` (|Δ| explained), 31-mer support co-located on one
   contig where claimed.
4. **Co-occurrence:** same accession — ideally same contig — support for
   ≥ 2 module baits or a junction bait before any "module/junction present"
   statement; dispersed single-bait hits stay "accession-level screen hit".
5. **Classification ladder** (report language, matches execute-bounded-ntm):
   `no hit` → `accession-level screen hit` → `module-only homology` →
   `supported junction/synteny` → `near-complete analogue` (requires
   multi-bait + alignment coverage across the reconstructed genome).

## 6. Bookkeeping contract

* One run dir per stage: `runs/<stage>/manifest.json`, `ledger.jsonl`,
  `cache/`, `normalized/`, `INSTRUCTIONS.md` (from `run_pilot.py prep`).
* Every request parameter, response checksum, retry, backoff, poll, and
  failure is in `ledger.jsonl` (append-only).
* Resume = rerun `fetch`; completed rows never re-fetch; failed rows keep
  errors; `status` prints the reconciliation table.
* Deliverable to `execute-bounded-ntm`: manifests + ledgers + normalized
  TSVs + this plan + README citations, plus the bait→hit→contig confirmation
  ledger from §5.
