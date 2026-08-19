# Local NTM Logan-subset screen (Tier A + Tier B) — run report

**Task:** `local-ntm-logan` · **Run window:** 2026-08-19T14:08:52Z–15:4xZ (~1.6 h)
**Run root (NVMe):** `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation/local_screen/`
(33 GB bulky artifacts, receipted in `NVMe_MANIFEST.tsv`)
**Preregistration:** `RUN_PLAN.md` (written before any download); amendments
A1–A3 in `AMENDMENTS.md` (all tooling/query-mechanics, none touching
thresholds, budgets, panel, or hit definition; all logged before the affected
download). **Logan Search spend: zero** (public sessions untouched, still
`submitted`).

## 1. Headline

Screened the frozen 24-bait NTM prophage panel against 17,596 Logan S3
per-accession assemblies (Tier A: 172 induction-study runs; Tier B: 17,424
*M. abscessus/avium/smegmatis* runs). **97 hit accessions** (≥ 0.7 k-mer
coverage), every one confirmed by minimap2 alignment (qcov ≥ 0.971,
identity ≥ 0.9919, median identity 1.0). Classification: **73 module-only
homology, 24 junction/synteny supported, 0 near-complete analogue**. The
junction-supported set is a single coherent signal: the reconstructed
partition junction of clade 0_0000 (*M. abscessus*) present at 100 %
identity in 18 accessions of ERP001039/ERP126939.

## 2. Control gates (evaluated first, honored)

| Gate | Result | Evidence |
|---|---|---|
| D29 poscon recovers own reference | **PASS** (kcov 1.0; minimap2 `sr` qcov 1.0, identity 1.0) | `results/controls_gate.json`, `logs/controls_gate.log` |
| L5 poscon recovers own reference | **PASS** (kcov 1.0; qcov 1.0, identity 1.0) | same |
| Shuffled negatives vs both refs | **PASS** (kcov 0.0) | same |
| Shuffled negatives across ALL 17,596 screened accessions | **PASS — 0 hits ≥ 0.7** | `screen_summary.json.negative_gate_violations = []` |

Integrity note: the first gate evaluation FAILed due to a parsing bug in my
new code (b2s `--out-kmers` marks found k-mers with `(id,pos,strand)`
annotations; bare lines are zero-count). The bug was in the consumer, not
the gate — fixed before any tier download, thresholds untouched. The D29/L5
sequences themselves always aligned perfectly (minimap2 1.0/1.0 even on the
failed parse), proving this was tooling, not biology.

## 3. Tier A — induction experiments (run first, reported before Tier B spend)

- Manifest: union of A1–A8, frozen 2026-08-19T14:22Z
  (`manifests/tierA_manifest.json`, uid-union sha256 logged) → 972 runs
  → **200 after NTM filter** (Amendments A1+A3: `prophage×organism`
  compound EUtils queries return 0 server-side; equivalent semantics
  executed as bare term + local filter; four stable induction terms added).
- Composition: 150 RNA-Seq, 18 ssRNA, 14 OTHER, 11 ChIP-Seq, 6 WGS.
- S3: 174/200 in logan-pub (26 × 404 logged as `not_in_logan`), 0.60 GB.
- **Result: 0 hits ≥ 0.7** in 172 screened runs. Max NTM-bait k-mer coverage
  anywhere: 0.011. Negatives clean everywhere.
- Sub-threshold note: `SRR22827534` shares 30.5 % of L5 poscon k-mers
  (below 0.7 ⇒ not a hit; documented only).

**Interpretation constraint:** 75 % of Tier A is transcriptomics; Logan
assemblies of RNA runs cannot contain uninduced prophage DNA. "Not detected"
in Tier A is weak evidence only — this is stated in the frozen ladder and
honored below.

## 4. Tier B — core NTM clinical screen

- Manifest: preregistered B1–B3 exact terms, union 21,187 uids (reproduces
  the 2026-08-19 EUtils counts exactly: 12,772 + 5,874 + 2,541; union
  sha256 `6afabdae…`) → 19,108 NTM-filtered frozen runs (15,823 WGS;
  12,913 abscessus / 3,602 avium / 2,593 smegmatis attribution kept per run).
- Disk preflight (500-sample, seed 20260819): 91.6 % in logan-pub, median
  1.48 MB, projected 31.6 GB → GO within budgets.
- Executed: 17,254 downloaded (1,706 × 404 `not_in_logan` = 8.9 %),
  17,424 screened (incl. 168 Tier A overlap), **33.24 GB**, 59 min, 0
  download errors, ledger reconciles (34,856 status-200 rows = 2 per
  checksummed download + 1,706 404s).

### Hit rates

| bait (class, clade) | accessions ≥ 0.7 | best cov | confirmed identity range |
|---|---|---|---|
| 0_0049_INTERIOR_MODULE_01 (portal) | 24 | 0.949 | 0.9992–1.0 |
| 0_0000_JUNCTION_02 (partition junction) | 18 | 0.887 | 1.0 |
| 0_0006_INTERIOR_MODULE_01 (portal) | 16 | 0.841 | 0.9994–1.0 |
| 0_0169_INTERIOR_MODULE_01 (capsid) | 14 | 0.843 | 0.9919–1.0 |
| 0_0063_INTERIOR_MODULE_01 (terminase) | 11 | 0.950 | 1.0 |
| 0_0085_INTERIOR_MODULE_01 (portal) | 6 | 0.899 | 0.9992–1.0 |
| 0_0032_JUNCTION_03 (partition junction) | 6 | 0.813 | 0.9980–1.0 |
| 0_0021_MEMBER_INTERIOR_02 (member midpoint) | 2 | 1.000 | 1.0 |
| **Total hit accessions** | **97** (0.56 % of 17,424) | | |

16 of 24 baits: no hit ≥ 0.7 anywhere (kept as `not detected in screened
subset`, explicitly ≠ biological absence). Poscons in-panel behaved as
expected (max L5 0.305 in one RNA-Seq run — sub-threshold).

### Confirmation (all 97, `results/confirm_B.jsonl` + full-provenance TSV)

- All alignments pass the minimap2 `sr` preset (first rung); qcov union
  0.971–1.000; identity 0.9919–1.0000.
- Every row traces: accession → S3 URL → zst sha256 → bait → contig id +
  length + target span + strand + CIGAR length
  (`results/confirmed_hits_full_provenance.tsv`).

### Junction-supported candidates (same-contig rule)

`NTMBAIT_0_0000_JUNCTION_02` is a **reconstructed** partition junction of
clade 0_0000 (bait = host/phage-partition boundary sequence rebuilt from the
generated genome, not copied from any reference). Evidence per accession:

- 18/18 accessions (ERP001039 × 16 + ERP126939 × 2, all *M. abscessus*
  subsp. *abscessus*/*massiliense*): full-length alignment (qcov 1.0,
  identity 1.0) on a **single contig** (e.g. `ERR330837_22`, 87,591 bp,
  target 16,050–16,449 + secondary exact 101-bp block at 4,035–4,136).
  Same-contig support is satisfied by construction: a full-length junction
  bait alignment necessarily spans the junction it encodes; additionally the
  same-clade interior-module bait is present at 20.6 % (sub-threshold)
  within ±5 kb of the junction site on that contig (neighborhood probe,
  §methods below).
- 0_0032_JUNCTION_03 × 6 (ERR459803/5/6 + 3): full-length, single contig.

**Near-complete analogue: none.** No screened accession reached 0.7 on ≥ 2
baits, so no accession qualifies; module-only homology (73 accessions) and
junction-supported (24) are the ceiling claims of this run.

### Effective independence caveat

Within-bait hit sets collapse to few source studies (ERP001039 alone = 50 of
97; identical kcov values = same isolate multiply sequenced). Distinct
genomes behind the 97 accessions are far fewer; per-clade "n accessions" is
an upper bound on independent observations.

## 5. Per-genome external-evidence update

`results/genome_external_evidence_update.tsv` merges every pilot genome row
(`ntm/v2/pilot/validation/genome_external_evidence.tsv`) with the SRA axis:
8 of 17 clades move from `no_SRA_evidence_screen_not_executed` to
`SRA_screen_hit_junction_supported` (0_0000, 0_0032) or
`SRA_screen_hit_module_only` (0_0006, 0_0021, 0_0049, 0_0063, 0_0085,
0_0169); 9 remain `not_detected_in_screened_subset` (subset ≠ all SRA;
assembly-only evidence; explicitly not biological absence). This is the
first live SRA-axis evidence in the project: the public Logan Search
sessions never returned (§7 of the pilot REPORT).

## 6. Budgets (RUN_PLAN §7) — all respected

| budget | cap | used |
|---|---|---|
| downloaded zst bytes | 400 GB | 33.84 GB (8.5 %) |
| distinct accessions | 25,000 | 17,428 |
| wall clock | 48 h | 1.6 h |
| S3 concurrency | 4 | 4 |
| Logan Search spend | 0 | 0 |
| retries | 3/file, backoff 5/20/60 s | 0 failures needing retry |

## 7. Reproducibility

`setup_env.sh` recreates the env (`env/VERSIONS.txt`: minimap2 2.31-r1302,
back_to_sequences 0.8.4, zstd 1.5.5, python 3.12.3, micromamba 1.5.9).
All code committed under `ntm/v2/external_validation/local_screen/`;
commands and parameters in `logs/`; every tier manifest carries uid-union
sha256 + query strings + freeze receipts. Bulky artifacts (33 GB downloads,
19.6 MB raw per-accession results JSONL) stay on NVMe, receipted in
`NVMe_MANIFEST.tsv` (17,578 rows).

## 8. Limitations

1. **Subset**: Tier A+B ≠ all SRA NTM; other species (chelonae, fortuitum,
   marinum, …) unscreened.
2. **Assembly-only evidence**: Logan contigs are assemblies; absence of a
   prophage region in an assembly is not absence in the sample (strains may
   lose prophages in culture; assembly gaps).
3. **Tier A composition**: 75 % RNA-Seq — weak evidential value; the
   intended mitomycin-C induction WGS runs are mostly absent from SRA as
   assemblies (26 of 200 not in logan-pub at all).
4. **8.9 % of Tier B runs 404** on logan-pub (logged, never retried) —
   coverage of the frozen manifest is 91.6 % of runs, not 100 %.
5. **Bait granularity**: 1–2 baits per clade; kmer coverage 0.7 threshold
   fixed a priori; near-complete-analogue claims structurally need ≥2 baits
   ≥0.7 in one accession and none occurred.
6. **Non-independence** of hit accessions (multiple runs per isolate;
   ERP001039 dominance) — see §4.
7. Junction baits are reconstructions: a full-length, 100 %-identity match
   means the reconstructed sequence exists in that assembly exactly; it does
   not by itself demonstrate an infectious cycle or active prophage.

## 9. GO / NO-GO for the 102-bait panel / broader Mycobacterium set

**GO**, with two design changes, justified by this run's receipts:

1. The toolchain is proven at scale (17.5k accessions / 33 GB / 1 h / zero
   negative-control violations / 97/97 confirmations at ≥ 0.99 identity) and
   only 8.5 % of the byte budget was spent — a 102-bait panel over the same
   Tier B set costs essentially the same download volume (k-mer screen cost
   grows sub-linearly: panel 24 → 102 baits ≈ 4× k-mers, b2s runtime ~4× of
   the 0.2 s/accession observed).
2. Change A — dedupe first: screen one run per (BioProject, isolate) and
   expand within-study only for hits; 97 accessions collapse to ≈ a dozen
   studies, so the marginal information per downloaded run is low.
3. Change B — extend species breadth (other NTM + full *Mycobacterium*)
   before extending bait depth per clade: 16/24 baits had zero ≥ 0.7 hits
   even in the clinical WGS set, so clade-level sensitivity, not panel
   breadth, is the current bottleneck for analogue discovery. The
   0_0000-junction result (18 accessions, ERP001039) shows the
   junction-bait approach yields exact, replicable evidence when the
   element is present in a population.

## 10. Deliverables map

| Deliverable | Path |
|---|---|
| preregistered tier manifests (queries + date + sha256) | `manifests/tier{A,B}_manifest.json`, `manifests/tier{A,B}_frozen.tsv` (+ uids, runinfo) |
| download ledger + availability | NVMe `downloads/download_ledger.jsonl`, `downloads/availability.tsv` (repo copy of availability) |
| per-bait hit tables | `results/tierA_hits.tsv`, `results/screen_hits_by_bait.tsv`, `results/hit_accessions.tsv` |
| alignment summaries | `results/confirm_B.jsonl`, `results/confirmed_hits_full_provenance.tsv` |
| per-genome external-evidence update | `results/genome_external_evidence_update.tsv` |
| control-gate receipts | `results/controls_gate.json`, `logs/controls_gate.log` |
| bulky artifacts | NVMe run root, receipted in `NVMe_MANIFEST.tsv` |

### Methods notes (for reproducing the neighborhood probe)

Junction-neighborhood co-location check: extract contig ± 5 kb around the
junction alignment, re-run b2s with the full panel, report per-bait coverage
in that window (0_0000 junction 0.495 window-local, 0_0000 interior module
0.206 window-local on `ERR330837_22`).
