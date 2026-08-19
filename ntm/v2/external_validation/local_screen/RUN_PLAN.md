# RUN_PLAN — Local NTM Logan-subset screen (Tier A+B) via S3 assemblies

**Task:** `local-ntm-logan` · **Preregistered:** 2026-08-19T14:08:52Z (UTC)
**Supersedes spend on:** nothing — the public Logan Search queue is stalled
(both smoke sessions `Pending 0/1` 37+ h; `ntm/v2/pilot/REPORT.md` §7–8), so
this run spends **zero Logan Search budget** (public sessions stay
`submitted`; recheck only if the queue visibly drains).

Everything below is fixed **before any tier download**. Amendments require a
logged entry in `AMENDMENTS.md` (none at preregistration time).

## 1. Objective

Screen the frozen 24-bait NTM prophage panel against Logan's freely
downloadable per-accession assemblies (S3 `logan-pub`, contigs) for:

- **Tier A** — SRA runs from *M. smegmatis* (and other NTM where present)
  mitomycin-C / prophage-induction studies. Highest prior for NTM prophage
  baits. Run first; report before any Tier B spend.
- **Tier B** — all runs for *M. abscessus*, *M. avium*, *M. smegmatis*
  (≈21k runs per EUtils 2026-08-19; reproduced at preregistration:
  12,772 + 5,874 + 2,541 = 21,187).

## 2. Inputs (frozen — no redesign)

| Input | Path | Integrity |
|---|---|---|
| Bait panel (24 rows) | `ntm/v2/pilot/stage1_panel.fa` + `stage1_selection.tsv` | per-bait sha256 in TSV (15 interior-module, 1 member-interior, 5 junction, D29/L5 positive controls, 2 shuffled negatives) |
| Confirmation toolchain | `ntm/v2/pilot/confirm_hits.py` (validated on DRR000016) | extended into a batch driver here; ladder unchanged |
| Preregistered thresholds | `logan/PILOT_PLAN.md` §0,§2,§4 | 0.7 primary; no threshold fishing |

## 3. Screen definition

- **Bait k-mer set:** all distinct canonical (strand-min) 31-mers of the bait
  sequence (both orientations count as one k-mer), computed with the same
  canonicalization as `back_to_sequences` (verified on DRR000016 probe +
  reverse-complement probe before Tier A — see `logs/toolchain_check.log`).
- **Per-(bait, accession) kmer_coverage** = |bait 31-mers present anywhere in
  the accession's contigs (either strand)| / |distinct bait 31-mers|.
- **Screen hit** ⇔ `kmer_coverage ≥ 0.7` (mirrors PILOT_PLAN §2 Stage-1
  threshold 0.7). Report 0.5/0.9 sub-threshold counts for context only; they
  never promote a hit.
- Screening target: `https://s3.amazonaws.com/logan-pub/c/<acc>/<acc>.contigs.fa.zst`
  (contigs only; unitigs never downloaded in this run).

## 4. Control gates (evaluated FIRST; any failure → STOP + document)

1. **Positive controls:** `POSCON_D29_AF022214_2_MID1200` must reach
   `kmer_coverage ≥ 0.7` against full AF022214.2 (D29) **and** minimap2
   identity ≥ 0.7 on the matched contig; same for
   `POSCON_L5_NC_001335_1_MID1200` vs full NC_001335.1 (L5). Reference
   genomes fetched from NCBI EFetch, sha256-logged.
2. **Negative controls:** both shuffled baits must produce **0 screen hits**
   (`kmer_coverage ≥ 0.7`) against the two control references **and every
   Tier A accession actually screened**. Checked continuously through Tier A;
   any negative-control hit anywhere ⇒ STOP before Tier B.
3. No threshold adjustment is permitted between gate evaluation and reporting.

## 5. Tier A manifest (preregistered queries — exact strings)

EUtils `esearch.fcgi?db=sra`, `retmode=json`, executed 2026-08-19 (UTC),
uid lists + query strings + sha256 recorded in `manifests/tierA_manifest.json`
**before any Tier A download**. Flat queries only (nested-paren queries
return 0 on db=sra — calibration note, `logs/eutils_calibration.log`).

| id | term | count@prereg |
|---|---|---|
| A1 | `"Mycobacterium smegmatis"[Organism] AND mitomycin[All Fields]` | 46 |
| A2 | `"Mycobacterium smegmatis"[Organism] AND prophage[All Fields]` | 1,590 |
| A3 | `"Mycobacterium smegmatis"[Organism] AND temperate[All Fields]` | 22 |
| A4 | `"Mycobacterium"[Organism] AND mitomycin[All Fields]` | 9 (post-filter to NTM organisms via runinfo) |

Tier A set = union of uid lists, deduplicated, with per-uid query-provenance
(union size + final NTM-filtered size recorded post-fetch). Expected order:
hundreds to low thousands of accessions.

## 6. Tier B manifest (preregistered queries — exact strings)

Same execution protocol, recorded in `manifests/tierB_manifest.json`
**before any Tier B download**:

| id | term | count@prereg |
|---|---|---|
| B1 | `"Mycobacterium abscessus"[Organism]` | 12,772 |
| B2 | `"Mycobacterium avium"[Organism]` | 5,874 |
| B3 | `"Mycobacterium smegmatis"[Organism]` | 2,541 |

Tier B set = union (21,187 expected). Accessions already screened in Tier A
are not re-downloaded; tier attribution is kept per accession.

## 7. Budgets (hard stops — checked before every tier and continuously)

| Budget | Value |
|---|---|
| Downloaded bytes (zst, cumulative) | 400 GB |
| Distinct accessions downloaded | 25,000 |
| Wall clock from RUN_START (2026-08-19T14:08:52Z) | 48 h |
| S3 concurrent connections | ≤ 4 |
| Retries | ≤ 3 per file, exponential backoff (5/20/60 s) |
| EUtils request rate | ≤ 1 req/s (no API key) |
| Logan Search spend | 0 |

Disk preflight: S3 `Content-Length` sampled per tier → projected total vs
NVMe free space (4.0 TB at preregistration). Bulky files
(`downloads/`, per-accession screen scratch) live only under
`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation/local_screen/`;
the repo mirrors only code, configs, manifests, ledgers, hit tables,
summaries, and the report (see `NVMe_MANIFEST.tsv`).

Every downloaded file records: URL, HTTP status, size, ETag, Last-Modified,
sha256, timestamp in `downloads/download_ledger.jsonl` (append-only).
Files missing from S3 (404) are recorded as `not_in_logan` (never retried).

## 8. Confirmation (per screen hit; ladder from PILOT_PLAN §5, unchanged)

1. Extract contigs carrying bait k-mers (b2s localization,
   `--output-mapping-positions`).
2. `minimap2 -c --eqx` preset ladder `sr → asm5 → asm20`, cg:Z CIGAR
   parsing → query coverage + identity.
3. Junction/synteny claims require **same-contig** support (contigs-only run;
   the unitig-pair rule of the task does not apply). Co-occurrence of ≥ 2
   module baits or a junction bait on one contig precedes any
   "module/junction present" statement.
4. Classification ladder (frozen): `no hit` → `accession-level screen hit` →
   `module-only homology` → `junction/synteny supported` → `near-complete
   analogue`. "Not detected" ≠ biological absence (subset ≠ all SRA;
   assembly-only evidence; Logan contigs are assemblies, not raw reads).

## 9. Environment (reproducible)

`setup_env.sh` recreates the micromamba env at
`<nvme_root>/env` with pinned `back_to_sequences=0.8.4`, `minimap2=2.31`
(conda-forge/bioconda); `env/VERSIONS.txt` records exact versions. The
previous `/tmp/mmenv` copy (same versions, verified 2.31-r1302 / 0.8.4) did
not survive-by-design; this install is on NVMe.

## 10. Deliverables

1. `manifests/tierA_manifest.json`, `manifests/tierB_manifest.json` (queries,
   date, uid lists, sha256) — preregistered before downloads.
2. `downloads/download_ledger.jsonl` + `downloads/availability.tsv` (S3
   HEAD results incl. 404s).
3. Per-bait hit tables: `results/tier{A,B}_hits.tsv`, screen coverage
   matrices + alignment summaries `results/confirm_{tier}.jsonl`.
4. `results/genome_external_evidence_update.tsv` merged with
   `ntm/v2/pilot/validation/` categories.
5. `REPORT.md`: control-gate results, per-tier hit rates,
   junction-supported / near-complete-analogue candidates with contig-level
   evidence, limitations, GO/NO-GO for 102-bait panel or broader
   Mycobacterium set.
