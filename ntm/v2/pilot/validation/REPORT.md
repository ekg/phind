# Validation of the bounded NTM Logan pilot + external-evidence report (2026-08-17/18)

**Task:** `validate-logan-candidates` (independent validation; performed in worktree
`agent-66`, no service requests spent).

**Headline:** the pilot's documented NO-GO state is **independently reproduced and
confirmed**: zero Logan candidates exist (Stage-0 smoke never retrievable; 0 of 24
Stage-1 submissions spent). Therefore **no genome can be promoted above
`no hit (screen not executed)` on the Logan axis** — and this is explicitly *not*
evidence of biological absence. On the independent **public-reference axis** (6,266
curated actinobacteriophage records), all 17 tested NTM generated genomes now have
calibrated external-evidence categories, module/junction alignments, and functional-QC
cross-links. Two modules are fully present in public genomes; **zero of five
reconstructed junctions are observed in any public genome**.

## 1. Independent re-verification of the pilot state (receipt-backed)

Every check below was re-executed from raw artifacts (not read off REPORT.md):

| Claim in `ntm/v2/pilot/REPORT.md` | Re-verification result |
|---|---|
| Stage-0 session `kmviz-c112ba44-…` submitted 20:55:30Z, machine-captured | `logan/runs/stage0-smoke/ledger.jsonl` entry `submission_recorded` at 20:55:57Z with identical seq sha256 `ac58a664…`; manifest row `submitted` |
| 21 download-API polls, all HTTP 400 / 46 B / byte-identical body | 21 `http_request` entries (20 client polls 20:56–21:21 + 1 manual re-check 21:31:00), all `status_code=400`, `size_bytes=46`, sha256 `980998f5…`; min spacing 60 s, max 595 s |
| NO-GO decision per PILOT_PLAN §4 | `final_status` entry 21:31:12Z, `elapsed_since_submit_s=2136`, decision text matches preregistered rule |
| 0 of 24 Stage-1 submissions spent | `logan/runs/stage1-gbrefseq/manifest.json`: 24 rows, all `status=pending`, **0 rows carry a session_id**; no result ZIPs/`normalized/`/`cache/` anywhere under `logan/runs/` |
| Frozen panel integrity | all 24 `stage1_panel.fa` sequences re-hashed → sha256 match `stage1_selection.tsv` and manifest rows exactly; all 24 `submissions/*.fa` sequences byte-match their panel entries (headers carry `_t70` tag) |
| Deterministic preregistered selection | `select_stage1.py` re-run → byte-identical outputs, repo tree unchanged (`git status` clean); A=15 B=5 C=2 D=2 in preregistered order |

Counting note (honest discrepancy, not affecting the decision): REPORT.md says
"22 checksummed requests (20 client + 1 manual + 1 dashboard cross-check)". The
ledger holds 21 checksummed download-API requests; the **two** dashboard
session-loader cross-checks (21:03, 21:24, browser-side, both "query still running")
are described in prose but not ledger rows. The NO-GO trigger (35.6 min, 21 polls) is
unaffected.

## 2. Tested NTM generated genomes (what "tested" means here)

The frozen 24-row panel tests **17 unique ML-reconstructed clade-representative
genomes** (14 `interior_module` baits + 1 `member_interior` + 5 `junction` baits;
junction baits from clades 0_0000, 0_0024, 0_0073, 0_0032×2). Controls: D29/L5
positive slices + 2 shuffled negatives. All sequences extracted to
`tested_genomes.fa` with per-genome sha256 (`tested_genomes.sha256`).

## 3. External-evidence categories (public-reference axis)

Full table: `genome_external_evidence.tsv` (repo copy) / external root. Evidence:
mash 2.3 (k=21, s=20000, `-i`) vs 6,266-record panel; minimap2 2.31-r1302
`-x asm20` genome-level and `sr`+`asm20` bait-level alignments; categories
preregistered *here* as: module-complete (bait cov ≥0.9 & id ≥0.9), partial module
relative (cov ≥0.4 & id ≥0.85), distant mosaic (mash d ≤0.20 but no qualifying
module), no public relative (d >0.20 and zero alignments).

| Clade | Member species (n) | Category | Best public relative | Module bait | Junction bait |
|---|---|---|---|---|---|
| 0_0000 | M. abscessus ssp. massiliense (100) | partial module | MW584184.1 d=0.034 | 0.71 cov @0.99 | 0.48 @0.98 ✗ |
| 0_0006 | M. abscessus ssp. abscessus (9) | partial module | MW353181.1 d=0.079 | 0.78 @0.97 | – |
| 0_0019 | M. chelonae (6) | partial module | OQ417963.1 d=0.103 | 0.41 @0.89 | – |
| 0_0021 | M. abscessus ssp. abscessus (13) | partial module (borderline) | MW584166.1 d=0.064 | 0.97 @0.885 | – |
| 0_0023 | M. chelonae (2) | distant mosaic | MW584188.1 d=0.153 | 0.24 @0.94 | – |
| 0_0024 | M. immunogenum (19) | distant mosaic | MW570842.1 d=0.097 | – | 0.35 @0.87 ✗ |
| 0_0032 | M. abscessus (9) | distant mosaic | OQ417973.1 d=0.053 | – | 0.78 @0.93 ✗ |
| 0_0038 | M. franklinii (4) | partial module | MW584190.1 d=0.134 | 0.41 @0.91 | – |
| 0_0049 | M. abscessus (6) | **module-complete** | MW584165.1 d=0.035 | **1.00 @0.97 ✓** | – |
| 0_0063 | M. abscessus ssp. massiliense (5) | **module-complete** | MW584187.1 d=0.091 | **0.99 @0.92 ✓** (genome qcov 0.31) | – |
| 0_0066 | M. wolinskyi (2) | **no public relative** | d=0.301 | none | – |
| 0_0073 | M. immunogenum (9) | partial module | MW314859.1 d=0.118 | 0.50 @0.89 | 0.49 @0.90 ✗ |
| 0_0077 | Mycobacteroides sp. (2) | distant mosaic | OQ417963.1 d=0.151 | none | – |
| 0_0085 | M. abscessus ssp. massiliense (5) | **no public relative** | MW584167.1 d=0.264 | none | – |
| 0_0169 | M. chelonae (4) | **no public relative** | MW584160.1 d=0.268 | none | – |
| 0_0218 | M. paragordonae/pumcae (2) | **no public relative** | MW584168.1 d=0.406 | none | – |
| 0_0219 | M. trivialis (4) | **no public relative** | BabeRuth d=0.340 | none | – |

Key findings:

* **Junction support: zero.** None of the 5 reconstructed junction baits
  (all `observed_adjacency_count == 0` internally, co-occurrence 8–14) is observed
  contiguously (≥0.9 cov @ ≥0.9 id) in **any** of the 3,064 public
  mycobacteriophage genomes. Best partial: 0_0032_J3 0.78 cov @0.93 id (OQ417987-
  class relatives); flanks have relatives, the **join itself is unobserved anywhere
  public**. This is exactly the hypothesis the Logan screen was designed to test.
* **Mash overstates relatedness for these mosaics.** E.g. 0_0000 vs MW584184.1:
  d=0.034 (Jaccard-inflated by thousands of short shared tracts) yet contiguous
  nucleotide alignment covers only 9% of the reconstruction at 89% id. Alignment
  confirmation (the preregistered §5 workflow) is indispensable — k-mer screen hits
  alone would overclaim.
* Two A-tier clades contain a module that exists complete in a public genome
  (0_0049 portal-module in MW584165.1; 0_0063 terminase-module in MW584187.1) —
  accession-level Logan hits are *likely* for these if any SRA assembly carries
  the module; promotion to "module present" still requires same-accession
  co-occurrence + alignment per PILOT_PLAN §5.4.
* 5/17 clades have **no close relative in curated public space** (0_0066, 0_0085,
  0_0169, 0_0218, 0_0219) — the highest-information Logan targets; public absence
  is an absence *in the panel*, not in nature.

## 4. Functional-QC cross-link (Pharokka/CheckV; version match)

The established NTM v2 run (Pharokka v1.10.1, DB 9Aug2025; CheckV 1.1.1, DB v1.5;
`annotation/run_state/tool_provenance.json`) already includes all 17 tested genomes
(cohort `ntm2_ml_reconstructed`), so **tool versions/databases are identical by
construction** — no re-run was needed or performed (no divergence to explain).

* CheckV: 8/17 High-quality (completeness 91–100%), 9/17 Genome-fragment
  (58–90% complete), **0% contamination everywhere**, `provirus=No` for all
  (no host carry-over in the reconstructions). See `genome_external_evidence.tsv`.
* Key-gene content: all 17 carry ≥4 of 6 key phage protein categories
  (terminase/portal/capsid/tail; `n_key_proteins` 4–5).
* **Integrated vs free/circular context:** the member prophages are *integrated*
  (extracted with host flanks from NTM assemblies; source accessions in
  `pilot_selection.tsv`); 3/17 reconstructions carry an integrase gene (0_0024,
  0_0032, 0_0218) → temperate/integration-competent marker present; the public
  panel relatives are free/circular isolate genomes (PhagesDB/NCBI complete
  records). No circularization evidence exists for any reconstruction (they are
  in-silico traversals, several length-capped at ~150 kb).
* geNomad was not part of the NTM v2 annotation workflow (prophage calling used
  the established v2 caller); stated explicitly rather than silently omitted.

## 5. Control performance (executed offline where possible)

| Control | Logan screen | Offline alignment leg (executed now) |
|---|---|---|
| D29 slice (`AF022214.2`) | **not run** (gated by Stage-0) | recovers own reference in panel: PHPUB-000587 = AF022214.2, **100% cov** ✓ |
| L5 slice (`NC_001335.1`) | **not run** | recovers own reference: PHPUB-001403 = NC_001335.1, **100% cov** ✓ |
| Shuffled negatives ×2 | **not run** | **0 hits** vs 3,064 public mycobacteriophage genomes ✓ (expected-zero) |
| E. coli smoke | submitted ×2, never retrievable | n/a |

The preregistered Stage-1 control gates (positive ≥0.7 self-hit; negatives 0 hits at
thr 0.9) therefore remain **unevaluated on the Logan side and must run before any
candidate validation**; the panel rows are frozen and ready.

## 6. Threshold sensitivity, false-positive review, bias, limitations

* **Threshold sensitivity:** preregistered thr 0.7 fixed; ledgers show **no threshold
  fishing** (smoke stayed 0.5; no Stage-2/3 spend). Sensitivity expectation is now
  evidence-informed: k=31 exact-matching (Logan/kmindex) tolerates ≤~10–15% divergence;
  the module-complete clades (≥0.92 id) are within range; mosaic/distant clades may
  need the preregistered 0.5 rung — which is *pre-authorized* by PILOT_PLAN §2, not
  fishing. Note k-mer *containment* semantics: reconstructed-genome baits query
  1,170 31-mers; partial-module public relatives sharing ~50% of a bait's kmers sit
  near the 0.5 rung, consistent with observed alignment fractions.
* **False-positive review:** all 20 NTM baits passed the design-stage filters
  (`n_ambiguous=0`, `n_low_complexity=0`, `n_internal_dup=0`,
  `cross_bait_max_shared_frac=0`); host-like k-mers were filtered at design
  (max 20/1170 for one bait → 1,150 usable; see `rejection_ledger.tsv` for the
  rejected baits, e.g. host_like_frac 0.40–0.43). Residual FP risks documented for
  the future screen: (i) junction baits can match **assembly chimeras/contig breaks**
  in SRA unitigs (Logan results are assembly-derived, not read-derived); (ii) shared
  mosaic tracts can produce accession-level hits without module presence — the
  §5.4 same-accession co-occurrence rule exists for exactly this; (iii) integrated
  prophage loci in RefSeq assemblies may hit interior baits without representing a
  free natural analogue (context must be checked per contig).
* **Sampling bias:** public panel is dominated by *M. smegmatis*-host SEA-PHAGES
  isolates; clinical-NTM (abscessus/chelonae/immunogenum) phages are rare (the
  MW5841xx batch and a handful of OQ/PQ accessions are the closest relatives found).
  Logan's index is SRA assemblies (freeze 2025-12-31; kmindex sub-index SRA ≤2023)
  + GenBank/RefSeq references; NTM clinical WGS in SRA is sparse and geographically
  biased. "No public relative" ≠ absent in nature; it means absent from *this* panel.
* **Logan assembly limitations (documented, feed into GO/NO-GO):** unitig/contig
  artifacts (chimeras, breaks at repeats) can both create and destroy junction
  evidence; contigs—not reads—are retrieved, so read-level junction confirmation is
  not directly possible from Logan artifacts (ENA raw-read backfill is feasible for
  ~79% of run accessions per `ntm/v2/run_assemblies/ena_backfill_report.md` and is
  the designated follow-up for any junction-supported candidate); one-month result
  retention; and the observed service retrievability failure (2 sessions, morning
  human-run + machine-captured, both unretrievable >30 min).

## 7. Traceability & artifacts

* Every tested genome/bait traces to: bait id → `bait_manifest.tsv` (genome_id,
  coordinates, module, checksums) → frozen panel sha256 → pending manifest row.
  Zero candidates exist, so there are no raw Logan responses, accessions, contigs or
  read retrievals to trace — the absence itself is receipt-backed (§1).
* Bulky artifacts (mash sketches incl. 1.0 GB panel sketch, PAFs, per-genome ref
  FASTAs) live in the external evidence root
  `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/pilot_validation/` with
  `MANIFEST.sha256` (63 entries). Only small tables + this report are committed.
* Reproducibility: mash 2.3 (`-k21 -s20000 -i`), minimap2 2.31-r1302 (`asm20`,
  `sr`), both locally installed; commands recorded in MANIFEST header notes and
  this report.

## 8. GO/NO-GO recommendation for scaling beyond the bounded pilot

**NO-GO on scaling Logan spend today; GO on everything else.** Justification:

1. **Service reliability is the sole blocker and is unremediated:** two independent
   submissions (one human, one machine-captured) failed retrieval within the
   preregistered window; the pause rule (≥24 h) is active and correct. Scaling
   submissions before a documented Stage-0 pass + maintainer response
   (logan-search.org) would burn budget on an unreliable retrieval path.
   `resume-logan-stage-0` already encodes exactly this gate — it must complete first.
2. **Method GO:** bait panel frozen + verified deterministic; confirmation pipeline
   (S3 → b2s → minimap2) validated end-to-end on DRR000016; offline control leg now
   also passes (D29/L5 self-recovery, negatives zero).
3. **Science GO with recalibrated expectations:** module-level relatives are
   confirmed for 12/17 clades (2 module-complete) → accession-level hits are
   plausible for those; 5 clades have no public relative at all → highest value;
   **junction baits remain wholly unvalidated in any public genome** — the
   reconstructed-junction hypothesis is genuinely open, and only Logan (SRA space)
   can test it at scale.
4. **Trigger to revisit:** Stage-0 retrieval pass + control gates honored → release
   the frozen 24 submissions; any junction-supported hit → immediate ENA read-level
   confirmation before promotion (per §6 FP review).
