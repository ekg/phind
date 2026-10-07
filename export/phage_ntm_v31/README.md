# NTM v3.1 — prophage-derived phage fragments for synthesis

**1,304 ML fragments** (+ **893** ancestral-state equivalents) reconstructed from
**36,857 prophages** across the union NTM cohort (34,163 Mycobacteriaceae genomes,
26,890 NTM + 7,273 MTC).

> **Read this first.** These are **candidate fragments, not finished phage
> genomes.** They are genome-path reconstructions: each follows observed
> adjacencies from a real prophage and terminates at a natural terminus. ~46% are a
> single genome's sequence; the rest join 2+ genome contexts. Every record carries
> the metrics needed to decide how much to trust it. If you want the *previous*
> release, do not use it — `export/phage_ntm_v3/` is a chimera catalog (see its
> banner and `v3_mosaic_audit.tsv`).

## Files

| file | |
|---|---|
| `all_ntm_v31_ml_phage_genomes.fa.gz` (+ `.fai`/`.gzi`) | 1,304 ML fragments, **bgzip + faidx-indexed**; headers `>ntm31_<clade>_ML mode=ml n_partitions=N len=L seed=42` |
| `all_ntm_v31_ancestral_phage_genomes.fa.gz` | 893 ancestral-state equivalents |
| `release_manifest.tsv` | 1,304 rows — **the trust table** (see below) |
| `chimerism.tsv` | single-member vs chimeric, longest single-member run |
| `targets.tsv` | per-clade genome-length target + how it was derived |
| `v3_mosaic_audit.tsv` | why the previous release was withdrawn |

Decompress with `bgzip -d`/`gunzip`, or `samtools faidx <file>.gz ntm31_0_0000_ML`.

## The trust table (`release_manifest.tsv`)

| column | meaning |
|---|---|
| `n_members` | prophages in the clade (1 = singleton pass-through, an observed sequence) |
| `status` | `ml` (reconstructed) or `singleton` (observed, not reconstructed) |
| `target_len`, `target_source` | genome-length target: `checkv_aai_median` (marker-based, 607 clades) or `member_seq_median` (fallback, 697 clades — **biased ~22% low**) |
| `ml_length_bp`, `length_ratio_vs_target` | length and how it compares to the target |
| `checkv_expected_length`, `checkv_completeness`, `checkv_quality` | CheckV on the fragment |
| `n_partitions`, `terminated_by` | path length; `terminus` = natural stop, `budget` = safety cap (should not occur) |
| **`chimerism_verdict`** | **`single_member`** (the path is one genome's partition order — an observed sequence) or **`chimeric`** (contains a junction no single genome has) |
| **`longest_single_member_run_frac`** | fraction of the path that is one contiguous single-genome run |
| **`n_members_to_cover`** | minimum genomes needed to tile the path |

### How to pick fragments

- **Safest (observed sequence, no novel junction):** `chimerism_verdict == single_member` → **408 fragments**.
- **High confidence of being genome-like:** `length_ratio_vs_target` in [0.8, 1.2] (820/1,304) **and** `checkv_completeness >= 90` (**124 fragments**).
- **Novel:** 86.7% have no close nucleotide relative (MASH > 0.30) among 6,060
  curated actinobacteriophages — see `ntm/v3.1/NOVELTY.md` for controls and caveats
  (the reference is biased toward fast-growing-host isolates).
- **Avoid / treat as low confidence:** `target_source == member_seq_median` **and**
  `checkv_completeness < 50` — those are the systematically-short fallback clades.

## How these were made

1. **Cohort** — union of the 2026-09-22 delivery (26,499) and v2 holdings (34,163
   objects); prophage coordinates are the BV-BRC phigaro QC-passed calls.
2. **Prophages** — 36,857 extracted, 0 errors.
3. **Host clades** — MASH over 34,163 genomes → 443 clades.
4. **Prophage tight clades** — MASH thr 0.25, max 100 → 1,304 clades.
5. **Alignment + partition** — allwave (sparsified k-nearest) + impg window 500 →
   726,121 partitions.
6. **Reconstruction** — a walk over the clade's partition graph restricted to
   **observed adjacencies**, terminating at a natural terminus, targeting a
   per-clade genome length (CheckV marker-based where available).
   `--n-samples 25 --seed 42`.
7. **Annotation** — Pharokka + CheckV.

Full method, A/B and analysis: `ntm/v3.1/REPORT.md`, `NOVELTY.md`, `SPEC.md`.

## QC summary (1,304 ML fragments)

| | |
|---|---|
| median length | **14,324 bp** (max 105,960) |
| length ÷ CheckV-expected, median | **1.00×** (820/1,304 within ±20%) |
| CheckV quality | 124 high, 255 medium, 647 low, 278 not-determined |
| CheckV completeness | median **28.4%** |
| chimerism | **45.7% single-member**, median longest single-member run 85.7% |

**Completeness is low and that is expected and honest.** The previous release
scored 99.6% by concatenating every common block — a chimera is "complete" because
it contains every marker gene. At true genome scale the number is ~28% overall and
**73.3% for the 336 largest clades**. High completeness and true genome length are
mutually exclusive for these clades.

## Caveats

1. ~54% of reconstructions contain at least one **junction not observed in any
   single genome** — check `chimerism_verdict` before synthesising.
2. 697/1,304 clades used the **biased-low fallback target** (~22% short).
3. Many source prophage calls are **fragments**, which bounds what stitching can
   recover.
4. Novelty is a **nucleotide** measure; protein-level comparison is still pending.
5. Singletons (411) are **observed sequences passed through**, not reconstructions.
