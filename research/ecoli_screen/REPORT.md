# Local E. coli Logan-subset screen (stratified Tier E1+E2) — run report

**Task:** `local-e-coli` · **Run window:** 2026-08-22T00:25:00Z (preregistration) →
2026-08-22T05:15Z (pipeline complete; ~4.8 h wall)
**Run root (NVMe):** `/mnt/nvme3n1/erikg/phind-genome-work/ecoli_screen/`
(86 GB bulky artifacts incl. env + downloads, receipted in `NVMe_MANIFEST.tsv`,
25,509 rows)
**Preregistration:** `RUN_PLAN.md` written before any assembly download;
amendments A1–A3 in `AMENDMENTS.md` (all tooling/transport, none touching
thresholds, budgets, panel, tier rules, or hit definition; all logged before
the affected execution). **Logan Search spend: zero.**

## 1. Headline

Screened the frozen 113-bait E. coli prophage panel
(`research/ecoli_bait/panel/baits.fa`, run id `2c0055272b71cc48`) against
**15,105 Logan S3 per-accession assemblies** — a preregistered stratified
subset of the 538,347-run E. coli WGS universe (Tier E1: 4,266 screened of
11,975 frozen — 62.8 % of E1 runs are not in logan-pub; Tier E2: 10,839
screened of 12,000 frozen, 88.4 % availability). **10,389 hit accessions**
(≥ 0.7 k-mer coverage on ≥ 1 non-control bait), every one
alignment-confirmed by the minimap2 ladder (22,625 confirmed (bait,
accession) pairs; qcov 0.746–1.0 median 1.0; identity 0.719–1.0 median
0.995). Classification (frozen ladder): **3,731 module-only homology, 322
junction/synteny supported, 6,336 near-complete analogue** (task rule: ≥ 2
non-control baits ≥ 0.7 in one accession). 73 of 107 non-control baits hit
≥ 0.7 somewhere; 39 of 59 junction baits.

The NTM-calibration predictions are confirmed in absolute terms: junction
confirmations 24 → 322 accessions (13×), near-complete analogues 0 → 6,336
(∞), known-phage novelty 92.3 % (§6).

## 2. Stratification (frozen before any download — receipts)

- **Universe query (frozen):** `"Escherichia coli"[Organism] AND wgs[strategy]`
  → 538,919 experiments (calibration in `logs/eutils_calibration.log`;
  `strategy wgs` bare-token form misparses to 1,050 — rejected). Uid list
  sha256 `94a2a5f8…`.
- **Metadata sweep:** esummary over all uids (Amendment A1: HTTP POST after
  HTTP 414 at 7-digit uids) → 538,784 rows → frozen WGS universe
  `library_strategy==WGS AND scientific_name startswith "Escherichia coli"`
  → **538,347 runs**, sha256 `af56e85b…` (freeze receipt in
  `ecoli_universe_manifest.json`).
- **Tier E1 (depth, seed 20260822):** BioProjects ≥ 100 frozen WGS runs
  whose study title matches an archetype by preregistered substring rules →
  18 projects / **11,975 runs** (in the 8–12k window): 16 ST131 studies
  (largest SRP362296 551, SRP068615 417, ERP165651 363, SRP254876 344 —
  population dynamics, transmission, carbapenemase, food+human) + 2
  global-diversity studies (ERP179860 "Global distribution of O-serotypes
  and antibiotic resistance" 7,498; SRP456113 "Global serotyping of
  invasive E. coli" 882). Longitudinal archetype matched **0** projects
  under the frozen title rule (documented here; several clearly
  longitudinal studies exist, e.g. ERP166924 "Longitudinal veal study",
  but their titles matched no preregistered string — no post-hoc scope
  growth). Manifest: `manifests/tierE1_manifest.json` (per-project
  justification, bioproject accessions, sha256 `ca768d6e…`).
- **Tier E2 (breadth, seed 20260822):** frame = frozen universe minus all
  E1 BioProjects; stratified by BioProject with **cap 60 runs/project**
  (anti-clone-collapse; 95 projects sit at the cap, incl. the 85k-run
  SRP046387 and 38k-run SRP071789 surveillance monsters which would
  otherwise dominate exactly as ERP001039 did in NTM); per-project RNG
  `random.Random("20260822|<project>")`; accumulate to exactly 12,000 →
  **842 projects**. Manifest: `manifests/tierE2_manifest.json`
  (strata sizes + exclusion rules + sha256 `95f83b41…`).
- **E1∩E2 run overlap: 0** (enforced and verified).

## 3. Control gates (evaluated FIRST, honored)

| Gate | Result | Evidence |
|---|---|---|
| λ poscon (`ECBAIT_POSCON_LAMBDA_01`) vs NC_001416.1 | **PASS** (kcov 1.0; minimap2 `sr` qcov 1.0, identity 1.0) | `results/controls_gate.json`, `logs/controls_gate.log` |
| HK97 poscon (`ECBAIT_POSCON_HK97_02`) vs NC_002167.1 | **PASS** (kcov 1.0; qcov 1.0, identity 1.0) | same |
| Shuffled negatives vs both references | **PASS** (kcov 0.0) | same |
| Shuffled negatives across ALL 15,105 screened accessions | **PASS — 0 hits ≥ 0.7** | `screen_summary.json.negative_gate_violations = []` |
| Host controls (diagnostic, panel-expected `host_like`) | host-like as designed: HOST_01 ≥ 0.7 in 4,252/15,105 (28 %), > 0 in 10,711; HOST_02 ≥ 0.7 in 1,757, > 0 in 14,636 | excluded from gates and classification; confirms both the host mask and that the screened set is genuinely E. coli |

Gates were evaluated 2026-08-22T00:5xZ, before the first assembly download
(02:00Z). No threshold adjustment at any point.

## 4. Screen + confirmation

- **Screen:** canonical 31-mers, hit ⇔ k-mer coverage ≥ 0.7 (0.5/0.9
  context only). E1: 4,266 screened → 3,990 hit (93.5 %); E2: 10,839 →
  6,399 (59.0 %). 8,870 frozen runs were 404 (`not_in_logan`, never
  retried); 0 download failures; ledger reconciles
  (`downloads/download_ledger.jsonl` + `availability.tsv`).
- **Confirmation:** every non-control hit pair b2s-localized then
  minimap2 `-c --eqx` `sr → asm5 → asm20`; 22,625/22,625 aligned
  (qcov 0.746–1.0, identity 0.719–1.0). Every row of
  `results/confirmed_hits_full_provenance.tsv` traces accession → S3 URL →
  zst sha256 → bait → contig id + length + target span + strand + CIGAR +
  identity + classification (22,625 rows).
- **Junction/synteny (322 accessions):** same-contig rule satisfied by
  construction — a full-length junction-bait alignment spans the
  host/phage partition junction it encodes; junction-bait alignments:
  qcov 0.764–1.0 (median 0.995), identity 0.807–1.0 (median 0.996).
  Example: `ERR2205963_248` (22,463 bp contig), ECBAIT_0_0000_JUNCTION_02
  at 12,676–13,172 (strand −), identity 0.980.
- **Near-complete analogue (6,336 accessions):** task-definition rule
  (≥ 2 non-control baits ≥ 0.7 in one accession, each alignment-confirmed)
  — the preregistered primary statistic. The NTM run's stricter secondary
  rule (≥ 1 junction + ≥ 2 modules same contig + combined span ≥ 20 kb)
  fires 0 times here **and is structurally near-unreachable with this
  panel**: baits are 1.0–1.2 kb, so even 4 aligned baits span ≈ 4.8 kb ≪
  20 kb. Reported for calibration transparency, not promoted.
- **Baits with ≥ 1 confirmed hit:** 73/107 non-control (68 %): 39/59
  junction, 34/48 interior (module+member). Top: ECBAIT_10_0057_MEMBER
  (8,271 accessions), ECBAIT_10_0002_MEMBER (3,328), ECBAIT_10_0172_MEMBER
  (1,537), ECBAIT_0_0000_JUNCTION_02 (1,505 — a reconstructed junction,
  63× the NTM best junction count).
- **Effective-independence audit:** the 10,389 hit accessions span **510
  distinct BioProjects**; top contributors SRP456113 (748), SRP362296
  (471), SRP068615 (414), ERP165651 (361), SRP254876 (344) — no
  ERP001039-style single-study collapse (largest share 7.2 % of hits).
  Identical-download-sha256 groups: 0 (every assembly file unique).
  Within-study re-sequencing of related isolates certainly remains
  (ST131 depth tier especially), so per-clade "n accessions" stays an
  upper bound on independent observations.

## 5. Per-clade external-evidence update

`results/genome_external_evidence_update.tsv` merges the screen axis into
the 585-row per-genome functional report
(`research/phage_annotation/per_genome_annotation_qc.tsv`, source
`ecoli_ml`): **24/585 clades move to `SRA_screen_hit`** (561 remain
`no hit (not detected in screened subset)` — subset ≠ all SRA;
assembly-only evidence; explicitly not biological absence). Top clades by
hit volume: 10_0057 (8,318 hits, 3 baits, best cov 1.0), 10_0002 (5,332,
5 baits), 0_0000 (3,220, 4 baits), 10_0172 (1,592), 10_0044 (674), 1_0004
(595). This is the first live SRA-axis evidence for the E. coli ML clades.

## 6. Calibration vs NTM (`ntm/v2/external_validation/local_screen/REPORT.md`)

| metric | NTM (17,424 acc) | E. coli (15,105 acc) | prediction | outcome |
|---|---|---|---|---|
| hit accessions | 97 (0.56 %) | 10,389 (68.8 %) | higher | **confirmed, 123× rate** |
| junction/synteny supported | 24 acc (2 junction baits of 5) | **322 acc (39 junction baits of 59)** | more junction confirmations from dense phage sampling + dense SRA | **confirmed in absolute terms (13× accessions, 66 % of junction baits fire); as a share of hits it is lower (3.1 % vs 24.7 %) because the near-complete class absorbs multi-bait accessions** |
| near-complete analogue | 0 | **6,336 (41.9 % of screened)** | expected higher in E. coli | **confirmed, dramatically** |
| baits with ≥ 1 hit | 8/24 | 73/107 | higher | confirmed |
| identity of confirmations | 0.9919–1.0 | 0.719–1.0 (median 0.995) | — | comparable at the median; wider floor reflects genuine E. coli prophage diversity |
| clone-collapse | ERP001039 = 50/97 hits | max share 7.2 % (510 projects) | cap prevents collapse | **confirmed (Tier E2 cap works)** |
| Logan availability of manifest | 91.6 % | E1 37.2 % / E2 88.4 % | — | E1's curated ST131 studies are largely absent from logan-pub assemblies (documented limitation) |

**Known-phage calibration (novelty estimate):** of 22,625 confirmed
(bait, accession) pairs, aligning each hit bait against λ NC_001416.1,
HK97 NC_002167.1, P2 NC_001895.1 (rule: qcov ≥ 0.5 AND identity ≥ 0.7):
**1,748 known-phage-matched (7.7 %), 20,877 novel (92.3 %)**
(`results/known_phage_calibration.tsv/.json`). The matched fraction
concentrates in the phage-family-conserved module baits (terminase/portal);
reconstructed junction baits are novel by construction (host-prophage
boundaries).

## 7. Budgets (RUN_PLAN §8) — all respected

| budget | cap | used |
|---|---|---|
| downloaded zst bytes | 400 GB | **25.35 GB (6.3 %)** |
| distinct accessions | 30,000 | **15,105** |
| wall clock | 48 h | **≈ 3.2 h** pipeline (≈ 4.8 h incl. metadata sweep from preregistration) |
| S3 concurrency | ≤ 4 | 4 |
| retries | 3/file, 5/20/60 s | 0 failures needing retry |
| EUtils rate | ≤ 1 req/s | held (0.05 s HEAD spacing; 1.1 s EUtils spacing) |
| Logan Search spend | 0 | **0** |

Disk preflight per tier (500-sample, seed 20260822): E1 37.2 % in
logan-pub, projected 6.8 GB; E2 88.4 %, projected 17.1 GB → GO. 4.0 TB
free at preflight.

## 8. Limitations

1. **Subset**: 15,105 of 538,347 frozen WGS runs (2.8 %), and only runs
   with Logan assemblies — 8,870/23,975 frozen runs 404. E1's curated
   ST131 studies are 62.8 % absent from logan-pub.
2. **Assembly-only evidence**: absence of a prophage region in an assembly
   is not absence in the sample (culture loss, assembly gaps).
3. **Near-complete rule granularity**: the task-definition ≥ 2 baits rule
   counts baits from the same reconstructed element (module + member +
   junctions of one clade) — it evidences a largely-complete element of
   that clade, not a complete novel phage genome. The NTM 20 kb span
   clause cannot fire with 1.0–1.2 kb baits (§4).
4. **Non-independence**: 510 BioProjects is a floor on study diversity,
   not on genome diversity; ST131-tier runs are related isolates by
   design.
5. **Host controls behave as host-like** (28 %/11.6 % ≥ 0.7) — expected;
   they validate the host mask but also show host DNA survives in baits'
  k-mer sets at sub-threshold levels for most accessions.
6. Junction baits are reconstructions: exact matches mean the
   reconstructed boundary sequence exists in that assembly; they do not
   demonstrate an active prophage or infectious cycle.

## 9. GO / NO-GO

**GO — extend to the 102-bait full panel over the same stratified frame.**
Justification from this run's receipts:

1. The toolchain is proven at 15.1k accessions / 25 GB / 3.2 h / zero
   negative-gate violations / 22,625/22,625 confirmations — with 93.7 % of
   the byte budget and 49.6 % of the accession budget unused.
2. Full-coverage budget math: 538,347 frozen runs × ~88 % logan
   availability × ~1.6 MB ≈ **760 GB** — 1.9× the byte cap. Full coverage
   therefore needs either a raised byte budget (2× is enough at current
   prices/time: extrapolated ≈ 33 h at 4-way concurrency) or the same
   stratified frame retained (recommended: it demonstrably prevents
   surveillance-project collapse).
3. The 102-bait full panel (100-genome selection, frozen in
   `research/ecoli_bait/panel/full_panel_selection.tsv`) costs the same
   download volume (k-mer screen grows sub-linearly; panel 113 → ~470
   baits ≈ 4× b2s runtime of the observed ~0.2 s/accession → ≈ 4 h).
4. Highest-value extensions, in order: (a) E1-replacement — a
   logan-enriched depth tier (the curated ST131 studies were mostly
   absent; pick depth projects by *observed* logan availability);
   (b) the remaining 482 clades currently `no hit` (561/585 clades
   unsampled by hits — breadth, not depth, is the bottleneck, same
   conclusion as NTM); (c) dedupe-first screening for full coverage
   (one run per BioProject × k-mer-cluster, expand only on hits).

## 10. Deliverables map

| Deliverable | Path |
|---|---|
| preregistered plan + amendments | `RUN_PLAN.md`, `AMENDMENTS.md` |
| universe + tier manifests (queries, seed, strata, sha256, justification) | `manifests/tierE1_manifest.json`, `tierE2_manifest.json`, frozen TSVs (NVMe originals + `ecoli_universe_manifest.json`) |
| control-gate receipts | `results/controls_gate.json`, `logs/controls_gate.log` |
| disk preflights | `results/E1_preflight.json`, `E2_preflight.json` |
| download ledger + availability | NVMe `downloads/download_ledger.jsonl`, `downloads/availability.tsv` |
| per-bait hit tables | `results/screen_hits_by_bait.tsv`, `results/hit_accessions.tsv` |
| alignment summaries (full provenance) | `results/confirmed_hits_full_provenance.tsv` (22,625 rows), NVMe `confirm/confirm_E1E2.jsonl` |
| per-clade external-evidence update | `results/genome_external_evidence_update.tsv` |
| known-phage calibration | `results/known_phage_calibration.tsv`, `.json` |
| run summary + independence audit | `results/screen_summary.json` |
| NVMe artifact receipt | `NVMe_MANIFEST.tsv` (25,509 rows; `.gz` in-repo) |
| screen logs | `logs/` (eutils, esummary, controls, screen, confirm) |
