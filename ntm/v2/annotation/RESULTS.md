# NTM v2 functional-QC run results (Pharokka + CheckV)

Execution record for task `run-ntm-v2`: the reviewed NTM v2 functional-QC
pipeline (see `README.md`, task `build-ntm-v2`) run on all **3,639 current
generated genomes — 2,388 ML + 1,251 ancestral** in the dedicated external
analysis root.

## Run provenance

| Item | Value |
|---|---|
| Analysis root (external) | `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation` |
| Input | `input/all_v2_phage_genomes.fa` — 3,639 records, SHA-256 `fa1403019ec6954dd6f422f72c55b010d29f20edde63912ab0c7e9c94f0799b9` |
| Pharokka | v1.10.1 (`micromamba/envs/pharokka`), `--mmseqs2_only --skip_extra_annotations --skip_mash -g prodigal-gv --locustag NTMV2 -t 64` |
| CheckV | v1.1.1 (`micromamba/envs/phage-annot`), `end_to_end … -t 64` |
| Databases (shared with v1, read-only) | `annotation/pharokka_db` (PHROG, 9 Sep 2025; full inventory in `run_state/*.done.json`), `annotation/checkv_db/checkv-db-v1.5` |
| Threads | capped at 64 |
| Stage markers | `run_state/{prepare,pharokka,checkv,report}.done.json` — all `exit_code=0` (finished 2026-08-17 09:20:13 → 09:31:18 UTC; driver process exited 0) |

An earlier attempt died mid-`mmseqs search` with its worker; the run was
resumed in-place after three driver fixes on this branch (resume-path KeyError,
PATH prepending of the tool env for the Phanotate dependency check, and `-f`
for stale partial Pharokka outdirs). The stale failure marker was overwritten
by the successful completion; no rerun was performed after completion.

## Output cardinality validation

- Per-genome report: **3,639 rows = 3,639 unique genome IDs** — exactly
  **2,388 `ntm2_ml` + 1,251 `ntm2_anc`** (set-equal to the prepared input index
  and to CheckV `quality_summary.tsv`; verified by `diff` of sorted ID lists).
- **Zero-gene rows: 0** — min `gene_count` = 1 (no identifier-mismatch or
  `no_annotation` tier anywhere), 206,249 CDS total (mean 56.7/genome).
- **CheckV coverage: 3,639/3,639 rows, 0 missing** — the builder's
  incomplete-coverage abort did not trigger (no `--missing-checkv-ok` used).
  Grades: 1 Complete, 519 High-quality, 926 Medium-quality, 1,566 Low-quality,
  627 Not-determined. miUVIG categories: 520 High-quality, 3,119 Genome-fragment.
- **Determinism**: rerunning `build_annotation_report_v2.py --root …` is
  byte-identical (SHA-256 of all four report TSVs unchanged on rerun).

Report SHA-256s (external, authoritative copies):

```
301b8addb5ba9d91b373fa2c6314c4a7e88401e2152b363068834850b46efa2b  per_genome_functional_qc.tsv
b46db82cbb2dd71a4a4ef0cc7e219cc86492d1a14de3f2c0acc2238027221cc6  summary_by_cohort.tsv
a561be7a32fe2df67c421918d701856b15e3b20e87f9cace601c8ae69fac0c4e  summary_by_source.tsv
833e0d9e8381dcd17d16a206d868e8d3ac3bb05ff916902b19f8b66c789cc829  candidate_functionality_by_cohort.tsv
```

The three small summary TSVs are committed alongside this file; the bulky
per-genome TSV, Pharokka/CheckV outputs and databases stay external.

## Results

### Candidate-functionality tiers (all 3,639 genomes; definitions in `README.md`)

| Tier | Total | % | ntm2_ml_reconstructed | ntm2_ml_singleton | ntm2_anc |
|---|---|---|---|---|---|
| A_strong_candidate | 748 | 20.6 | 269 | 207 | 272 |
| B_moderate_candidate | 455 | 12.5 | 199 | 61 | 195 |
| C_partial_module_evidence | 1,459 | 40.1 | 437 | 583 | 439 |
| D_modules_not_detected | 977 | 26.8 | 346 | 286 | 345 |
| no_annotation | 0 | 0.0 | 0 | 0 | 0 |
| **n** | **3,639** | | **1,251** | **1,137** | **1,251** |

### Cohort summary (medians / % of cohort)

| Cohort | n | median genes | % terminase | % capsid | % tail | % lysis | % temperate | % low-completeness | % contaminated | % tier A+B |
|---|---|---|---|---|---|---|---|---|---|---|
| ntm2_anc | 1,251 | 55 | 34.5 | 29.2 | 55.5 | 43.3 | 21.3 | 53.9 | 0.8 | 37.3 |
| ntm2_ml_reconstructed | 1,251 | 55 | 34.1 | 28.5 | 55.5 | 43.2 | 21.5 | 53.8 | 0.8 | 37.4 |
| ntm2_ml_singleton | 1,137 | 16 | 31.6 | 30.5 | 43.7 | 30.9 | 31.8 | 74.4 | 3.6 | 23.6 |

Full columns (truncation, hypothetical fraction, per-tier %) in
`summary_by_cohort.tsv` / `summary_by_source.tsv` / `candidate_functionality_by_cohort.tsv`.

## Reading

- **ML reconstructed genomes track the ancestral comparator closely** on every
  module-detection and QC metric (A+B: 37.4% vs 37.3%; median 55 genes; near-
  identical module profiles) — consistent with clade-faithful reconstruction.
- **Singletons are the low-quality tail** (median 16 genes, 74% low-completeness,
  3.6% contamination): short pass-through fragments dominate tiers C/D, as
  pre-registered in the README cohort design; they are reported separately so
  they don't drown the reconstructed-clade signal.
- `D_modules_not_detected` (26.8%) is a **detection** statement under
  mmseqs2-only PHROG matching (novel/divergent ORFans, short fragments), not
  proof of biological absence; caveats in `README.md`.
- Integrase/temperate signal is recorded as supportive context
  (prophage-derived genomes) and deliberately excluded from the tier ladder.
