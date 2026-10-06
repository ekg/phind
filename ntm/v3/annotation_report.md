# NTM v3 annotation report (Pharokka + CheckV)

Genomes: 2197 (1304 ML + 893 ancestral).

Numeric detail: `report/per_genome_functional_qc.tsv`, `report/summary_by_cohort.tsv`, `report/summary_by_source.tsv`.


## Candidate functionality by cohort

| cohort | A_strong_candidate | B_moderate_candidate | C_partial_module_evidence | D_modules_not_detected | no_annotation | total |
|---|---|---|---|---|---|---|
| ntm3_anc | 207 | 91 | 338 | 257 | 0 | 893 |
| ntm3_ml_reconstructed | 208 | 91 | 337 | 257 | 0 | 893 |
| ntm3_ml_singleton | 39 | 15 | 192 | 165 | 0 | 411 |

## Cohort module metrics

| cohort | n | median genes | %terminase | %capsid | %tail | %lysis | %tier A | %tier D |
|---|---|---|---|---|---|---|---|---|
| ntm3_anc | 893 | 108 | 34.5 | 29.2 | 47.4 | 37.5 | 23.2 | 28.8 |
| ntm3_ml_reconstructed | 893 | 108 | 34.9 | 29.1 | 47.4 | 37.4 | 23.3 | 28.8 |
| ntm3_ml_singleton | 411 | 10 | 19.7 | 18.7 | 27.7 | 17.8 | 9.5 | 40.1 |

Tier `D_modules_not_detected` is a detection result under mmseqs2-only PHROG matching, NOT proof of biological absence.

