# E. coli prophage-derived maximum-likelihood (ML) phage genomes

**585 synthetically-amenable phage genomes** mined from prophages found in
26,232 *E. coli* genome assemblies (132,393 prophage elements total).

Companion files:
- `phage_ml_genomes.fa.gz` — one ML phage genome per clade (585 entries, bgzipped + faidx-indexed)
- `phage_ml_genomes.tsv` — metadata table describing each genome

## What these genomes are

Each entry is a **maximum-likelihood (typical) genome** reconstructed from a
coherent clade of related prophages:

1. **MASH distance triangle/tree** over all 132,393 prophage elements clustered
   them by sequence content (not by host genome).
2. **Tight clades** (585) were formed so that every member prophage is within
   ~75%+ ANI of its clade representative (leader clustering, clade size ≤ 100).
   Divergent isolate prophages with no close relatives form singletons (106).
3. Within each clade, **all-wave alignment** (`allwave`, biWFA, k-nearest only),
   then **`impg partition`** into small homologous blocks (≤ 1 kb partitions).
4. A **weighted traversal** samples the most common path through each clade's
   partitions (probability ∝ partition occurrence); the **ML genome** is the
   majority-rule consensus of that typical path, with the option of an
   ancestral-state reconstruction (NJ + Fitch parsimony) as an alternative.

Headers encode provenance:
```
>clade_<id>_ML n_members=N n_partitions=P length=L status=ml|singleton genomes=G1,G2,...
```

## TSV columns

| column | meaning |
|---|---|
| `clade_id` | clade identifier (`<community>_<NNNN>`) |
| `community` | originating community (0-11) |
| `n_members` | prophages in the clade (1 for singletons) |
| `n_partitions` | aligned blocks used to build the genome |
| `length_bp` | genome length |
| `gc_pct` | GC content (%) |
| `status` | `ml` (reconstructed) or `singleton` (observed sequence) |
| `n_source_genomes` | number of distinct *E. coli* genomes carrying this clade |
| `source_genomes` | NCBI assembly accessions (GCF_...) of those genomes |

## QC summary (585 genomes)

| length | count | note |
|---|---|---|
| < 10 kb | 78 | mostly singletons (short fragments) |
| 10–50 kb | 235 | typical small–medium phages |
| 50–100 kb | 198 | typical phages |
| 100–150 kb | 54 | large phages |
| ≥ 150 kb | 20 | large/jumbo prophage clades (max 255 kb) |

Median 45.2 kb, mean 54.0 kb. All non-singleton clades validated:
internal median pairwise MASH distance ≤ 0.25.

## Caveats for synthesis

1. **Sparse pairwise evidence in large clades.** k-nearest alignment aligns
   ~7% of possible pairs in 100-member clades. The ML genomes are statistically
   grounded (weighted sampling) but represent an inferred typical path — treat
   as candidates requiring follow-up validation, especially the 20 genomes
   ≥ 150 kb.
2. **Long constructs.** Genomes ≥ 150 kb may exceed practical single-fragment
   synthesis limits; they can be split at partition boundaries.
3. **Rare modules are sampled, not excluded.** Rare prophage genes appear with
   low probability; each genome is one sample (seed 42) of the clade's typical
   distribution. Minor variants arise across seeds.
4. **Singletons** are observed sequences, not reconstructions.

## Reproducibility

All inputs, scripts, and validation reports live in the source repository
(`research/mash_tree/`, `research/clades/`, `research/ml_phage_genomes/`,
`scripts/`). Core pipeline: `build_tight_clades.py`, `per_clade_alignment_pipeline.py`,
`segment_paf.py`, `traverse_partitions.py` (43 tests).
