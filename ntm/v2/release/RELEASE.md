# PHIND-NTM v2 — ML & ancestral phage genomes (QC-passed collaborator cohort)

**Deliverable.** Maximum-likelihood (and ancestral) phage genomes mined from
prophages called on the collaborator's **QC-passed NTM cohort**, one ML genome
per prophage tight clade, **labeled by v2 host clade**.

**Date.** 2026-08-16 · tasks: ntm-v2-download → ntm-v2-host → ntm-v2-clades →
ntm-v2-ml → **ntm-v2-release** (this file). Replaces the v1 release
(`ntm/RELEASE.md`, 913 genomes, geNomad calls); see
`v1_v2_comparison.md` for the redo's effect.

## Cohort composition

The input is the collaborator's **QC-passed accession list** (33,082 rows:
`NTM_QC_passed_accession_list.tsv`), not the v1 NCBI-family query. Composition:

| stratum | genomes | prophage calls | status in this release |
|---|--:|--:|---|
| NCBI subset (canonical, deduped GCA-preferred + 332 GCF-only) | **13,122** | 8,502 used, on 3,527 accessions (3,466 canonical + 61 GCF twins of canonical genomes) | **used** (downloaded/linked, host-clade-joined) |
| Run-assemblies (ERR/SRR/DRR) | **9,543** | 14,741 calls on **5,907** prophage-bearing runs | prophage calls **shipped** by collaborator; run FASTAs **pending** collaborator (task ntm-v2-run) — 0 ingested, excluded here |

Notes:

- The canonical NCBI stratum is 13,122 assemblies = 12,790 GCA + 332
  GCF-only (GCA/GCF twins deduped by assembly numeric;
  `host_clades.tsv` holds 13,122 data rows / 13,123 lines incl. header —
  the chat-agent update's "13,123" is the line count).
- **TB-complex kept verbatim, decision pending.** The QC-passed list contains
  **13,227 TB-complex-labeled accessions** (13,216 exact *M. tuberculosis*,
  9 *M. tuberculosis* subsp. *tuberculosis*, 2 *M. [tuberculosis]*; 7,056 GCA /
  6,171 GCF, all NCBI) despite being labelled "NTM". By ANI-to-closest-type-strain
  13,220/13,227 (99.95%) are closest to ***M. decipiens***, i.e. the label
  "M. tuberculosis" is a collaborator species call, not an MTC membership call.
  7,244 of them survive canonical dedup into the cohort (host_clade_0001,
  7,465 members with near-neighbours). Per the chat-agent update (2026-08-16)
  they are **kept verbatim pending user decision**; if the user later rules
  them out-of-scope, the 23 ML genomes joined to host_clade_0001 (see
  `per_ntm_clade_summary.tsv`) and downstream counts should be re-derived.
  The update's "~6.6k unique" figure matches no exact cut of the data
  (closest measured: 6,171 GCF-only TB rows excluded by canonical dedup);
  the measured values above are authoritative.

## Pipeline (end-to-end on the v2 cohort)

1. **Inputs** — collaborator QC-passed accession list (33,082) + prophage
   master manifest (51,004 rows; 30,165 prophage rows);
   `ntm/v2/scripts/validate_inputs.py`.
2. **Download** — 13,122/13,122 canonical genomes resolved as PanSN bgzip
   (11,972 downloaded + 1,150 linked from v1, 0 failures);
   `ntm/v2/scripts/download_ntm_genomes_v2.py`.
3. **Prophage extraction** — collaborator calls on canonical cohort →
   `full_prophages.fa`, **8,502 prophages** (run-assembly calls pending).
4. **Prophage MASH + tight clades** (thr 0.25 / max 100, community 0) →
   **2,388 clades** (1,251 alignable ≥2 members + 1,137 singletons; median
   internal mash distance 0.0060).
5. **Per-clade allwave + partition** (k-nearest `tree:5:0:0.0`, window 500) →
   ≤1 kb partitions.
6. **Two-level ML + ancestral traversal** (`--n-samples 5 --seed 42`, max
   length 150 kb; singletons pass through) → 2,388 ML + 1,251 ancestral
   genomes (task ntm-v2-ml; 0 traversal failures, 39 warnings in
   `ml/warnings.tsv`).
7. **Host clades** — full re-sketch of the 13,122-genome cohort
   (`mash sketch -k 21 -s 10000`, cluster at dist ≤ 0.05) → 1,044 host clades
   (`host_clade_0001`–`1044`); task ntm-v2-host.
8. **Release join** (this task, `ntm/v2/scripts/build_v2_release.py`) —
   prophage → source accession → host clade; FASTA headers and tables
   annotated with host clades + species.

## Files (`ntm/v2/release/`)

| File | Contents |
|---|---|
| `all_ntm_ml_phage_genomes.fa` | **2,388 ML phage genomes** (1,251 ML + 1,137 singletons) |
| `all_ntm_ancestral_phage_genomes.fa` | 1,251 ancestral-state genomes |
| `release_manifest.tsv` | per-clade: members, status, ancestral, length, median member length, source genomes, host clades, species, join status, ML-stage flags |
| `host_range.tsv` | prophage clade → host clade ids + species spanned |
| `per_ntm_clade_summary.tsv` | host clade → # ML genomes / # clades |
| `unjoined_genomes.tsv` | member-level join exceptions (all 117 are GCA↔GCF twin joins; 0 unjoined) |
| `qc_table.tsv` / `qc_flagged.tsv` | per-genome QC (`scripts/qc_ml_phage_genomes.py`) |
| `v1_v2_comparison.md` | v1 (913-genome) vs v2 (2,388-genome) redo comparison |

Headers: `>ntm2_<cid>_ML status=ml|singleton n_members=N length=L host_clades=… species=…`

## QC (2,388 ML genomes, ~68 MB; `qc_table.tsv`)

- **Lengths:** min 237 bp, **median 20,268 bp**, mean 28,284 bp, max 150,000 bp
  (traversal cap). Buckets: <10kb **707** · 10–50kb **1,311** · 50–100kb **330**
  · 100–150kb **39** · ≥150kb **1**. Alignable-clade genomes only (n=1,251):
  median 30,277 bp. 707 genomes <10 kb are flagged `too_short` — 467 of them
  are singleton pass-throughs of genuinely short collaborator calls (237–9,999
  bp; the ML stage already warned on 23 genomes <1 kb).
- **GC median 63.81%** (min 34.5%, max 72.7%) — matches the mycobacterial
  host expectation (~64%); alignable-only median 63.72%.
- **N-runs:** median max-N-run **0**; 14 genomes contain any N; **6 genomes
  have an N-run ≥100** (max 5,507 — clade 0_1848's ML genome; see
  `qc_flagged.tsv`).
- **Member recovery (completeness/identity):** each member prophage was mapped
  back to its ML genome (mappy/asm5, union of aligned blocks; ≥50% of member
  length covered counts as recovered). Singletons: 1,137/1,137 at 100%/100%.
  Alignable clades: **1,014/1,251 (81%) below the 50% recovery threshold**
  (median 0.0), i.e. the ML mosaic represents a typical member only partially,
  although aligned blocks are near-identical (overall identity median 1.0;
  147 genomes <0.8 flagged `low_identity`). Mechanism: the partition set is
  built from k-nearest pairwise alignments and can be redundant (e.g. clade
  0_0013: three 100%-identical 27.7-kb members → 55 partitions → 54.8-kb ML
  genome, each member ~21% recovered), so the stitched ML genome is a
  weighted-path mosaic over redundant partitions, not a member-length
  consensus. **v1 shows the same property under the same metric** (461/913
  below 0.5, median 0.25) — this was never measured before because the QC
  script's completeness/identity computation was broken (see below). Treat
  ML genomes of large alignable clades as clade-level mosaics; use the
  ancestral set (`--mode ancestral`) or per-partition intermediates when a
  member-faithful consensus is needed.
- **Flags (1,689/2,388):** too_short 707 · long_n_run 6 · low_completeness
  1,014 · low_identity 147 · too_long 0.
- **QC script fixes (this task):** `scripts/qc_ml_phage_genomes.py` had (a) an
  inverted identity formula — mappy `hit.mlen` is *matched* bases, the script
  computed `1 − matches/length` — and (b) best-single-hit coverage that cannot
  see multi-block mosaics. Both fixed; also accepts `ntm_`/`ntm2_` FASTA id
  prefixes and flat clade-dir layouts. v1 baseline re-run with the fixed
  script lives in `v1_qc_baseline/` for like-for-like comparison.

## Host-clade join (every genome joined or accounted)

- **2,388/2,388 ML genomes joined to ≥1 host clade; 0 unjoined.**
- 117 members (61 unique accessions, all GCF) sit on the GCF side of a
  GCA/GCF twin whose canonical cohort copy is the GCA assembly; joined via
  the same-numeric twin and listed with reason `twin_joined` in
  `unjoined_genomes.tsv` (85 clades touched; no clade relied on twins alone
  for its only host clade after fallback — 0 fully-unjoinable clades).
- Host clades represented: **623 of 1,044** (59.9%).

## ML genomes per NTM clade (top host clades)

| host clade | species (dominant) | ML genomes | clade size |
|---|--:|--:|--:|
| host_clade_0002 | *Mycobacteroides abscessus* subsp. *abscessus* | 1,182 | 2,336 |
| host_clade_0006 | *Mycobacterium intracellulare* | 87 | 152 |
| host_clade_0004 | *Mycobacterium avium* subsp. *hominissuis* | 84 | 507 |
| host_clade_0011 | *Mycobacteroides immunogenum* | 36 | 17 |
| host_clade_0008 | *Mycolicibacterium fortuitum* | 36 | 61 |
| host_clade_0005 | *Mycobacterium kansasii* | 30 | 187 |
| host_clade_0009 | *Mycolicibacterium senegalense* | 30 | 25 |
| host_clade_0023 | *Mycolicibacterium gilvum* | 26 | 6 |
| host_clade_0001 | *M. tuberculosis* label (ANI: *M. decipiens*) — see cohort note | 23 | 7,465 |
| host_clade_0003 | *Mycobacterium ulcerans* (+ *M. marinum*) | 20 | 1,139 |

## Caveats

1. **Run-assemblies excluded.** 9,543 run ids (5,907 prophage-bearing,
   14,741 calls) await collaborator FASTAs; when they land, re-run steps 3–8
   — the schema extends unchanged (see `host_clades/host_clades_report.md`).
2. **TB-complex verbatim** (see cohort composition): 23 ML genomes ride on
   host_clade_0001 pending the user's scope decision.
3. **Collaborator vs geNomad caller disagreement is large** — on 2,867 shared
   genomes only 3,448/8,351 v1 geNomad prophages have a ≥50%-reciprocal
   collaborator counterpart (44% have no overlapping call at all). See
   `v1_v2_comparison.md` before assuming continuity with v1.
4. **ML genomes of alignable clades are mosaics** (QC above); lengths can
   exceed member lengths (up to the 150-kb cap) and per-member recovery is
   partial. v1 behaves the same way under the fixed metric.
5. Mash-based clades/trees are sketch approximations; k-nearest
   sparsification limits pairwise evidence in 100-member clades; rare
   partitions are sampled (seed 42, one draw per clade).
6. Species labels are the collaborator's calls verbatim (v2 accession list
   has no organism column); `species=` header fields list up to 3 distinct
   source-genome labels (full set in `host_range.tsv`).
