# PHIND — ML (maximum-likelihood) phage genome release

**Purpose.** Synthesizable ML phage genomes mined from *E. coli* prophages, for
collaborator evaluation and potential use in a phage therapy trial.

**Date.** 2026-08-06 (clade/partition run completed 2026-08-05)

## Files

| File | Description |
|---|---|
| `all_ml_phage_genomes.fa` | 585 ML phage genomes, one per tight clade |
| `release_manifest.tsv` | Per-genome metadata: clade, community, n_members, n_partitions, length, status, source genomes |
| `clade_genomes.json` | clade_id → list of source *E. coli* genomes carrying the clade's prophages |

## How these were made

1. **MASH triangle/tree** over all 132,393 prophage elements
   (`prophage_homology_survey/full_prophages.fa`, `research/mash_tree/`).
2. **Tight clades** (585) from the 12 communities by leader clustering on MASH
   distance: `--threshold 0.25 --max-size 100` (`scripts/build_tight_clades.py`).
   Every clade member is within 0.25 MASH of its representative; clade size
   <= 100. 106 singletons (divergent isolates) are included as their own
   genome.
3. **All-wave alignment** per clade with `allwave` (biWFA, pangenome/allwave):
   k-nearest only, k-farthest = 0 (no stranger-joining), scores
   `0,5,8,2,24,1`, mash orientation.
4. **Small partitions**: PAF CIGAR-segmented (`scripts/segment_paf.py`,
   window 500, max span 1000 bp) then `impg partition` → per-clade BED + MAF.
5. **ML traversal** (`scripts/traverse_partitions.py`, two-level ML):
   weighted common-path sampling across partitions (weight =
   occurrence^alpha, alpha=1, n_samples=5, seed=42) × majority-rule consensus
   within each partition; per-clade length budget = 1.6 × max member length
   (30 kb floor) so the genome tracks the clade's real phage size.

Headers: `>clade_<cid>_ML n_members=N n_partitions=P length=L status=ml|singleton genomes=G1,G2,...`

## QC summary (585 genomes, 31,613,212 bp total)

| length bucket | count | note |
|---|---|---|
| < 10 kb | 78 | mostly singletons (short/degenerate prophage fragments) |
| 10–50 kb | 235 | typical small-medium phages |
| 50–100 kb | 198 | typical medium phages |
| 100–150 kb | 54 | large phages (incl. 100-member clades) |
| ≥ 150 kb | 20 | jumbo/large-phage clades (max 255 kb; members up to 160 kb) |

- Median 45.2 kb, mean 54.0 kb, q25 20.6 kb, q75 77.8 kb — phage-typical.
- 0 traversal failures; all 585 clades represented.
- All non-singleton clades validated: internal median pairwise MASH ≤ 0.25.

## Caveats for synthesis / trial consideration

1. **Sparse alignment evidence in large clades.** k-nearest alignment
   (k=5–10) aligns ~7% of possible pairs in 100-member clades. ML genomes from
   these clades are built from a k-nearest graph — statistically grounded
   (weighted sampling) but with sparse pairwise evidence. For the ~20 clades
   with genomes ≥ 150 kb, consider treating them as candidates requiring
   follow-up validation, not confirmed genome sequences.
2. **Length outliers.** Genomes ≥ 150 kb come from clades whose members
   genuinely reach 100–160 kb (large/jumbo phages) — real, but at the extreme
   of what is commonly synthesized. Constructs > 100 kb may be impractical to
   synthesize in one piece; the ML genome can be split at partition boundaries
   if needed.
3. **Rare partitions are sampled, not excluded.** The traversal gives rare
   modules a low-but-nonzero chance to appear (per project decision). Each
   genome is one sample from the clade's typical-path distribution (seed 42);
   re-running with other seeds yields minor variants. `--mode ancestral`
   (NJ + Fitch parsimony per partition) is available for an ancestral-state
   alternative genome per clade.
4. **Singletons** (106) are single prophages with no close relatives — their
   "genome" is the observed sequence, not an ML reconstruction.

## Reproducibility

- Commands per clade: `research/clades/<community>/<clade>/commands.log`
- Clade manifests: `research/clades/<community>/<clade>/manifest.json`
- Builder: `scripts/build_tight_clades.py`, `scripts/extract_clade_fasta.py`,
  `scripts/per_clade_alignment_pipeline.py`, `scripts/segment_paf.py`,
  `scripts/validate_clades.py`
- Traversal: `scripts/traverse_partitions.py` (+ 43 tests in
  `scripts/test_traverse_partitions.py`)
- Full upstream data: `research/mash_tree/` (sketch, 35 GB triangle, tree,
  labels), `research/ecor/` (ECOR mapping), `research/artifact_audit.md`.
