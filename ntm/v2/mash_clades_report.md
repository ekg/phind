# NTM v2 — prophage MASH + tight clades report

Generated: 2026-08-16T14:15:41Z

## Prophage set

- `full_prophages.fa`: **8502 prophages** (NCBI-only; v1 was 10,438 — run-assembly prophages pending collaborator FASTAs, task ntm-v2-run; count inside the [5000, 35000] band, no flag)
- Mash sketch: `-i -k 21 -s 10000` (same as v1)

## Clades (threshold 0.25, max-size 100, community 0)

- clade count: **2388**
- alignable (>=2 members): **1251**
- singletons: **1137**
- median internal mash distance per clade: **0.0060** (v1: 0.074)
- assignment check: clade members + singletons == 8502 == 8502 prophages (every prophage assigned exactly once)

## Outputs

- `mash_clades/`: `prophages.msh`, `prophages.dist.tsv`, `prophages_mash.dist` (float32 upper triangle), `ids.txt`, `labels.csv`, `full_prophages.idx.json`, `prophages_tree.nwk`, `tree_stats.json`
- `clades/`: `0/` (tight_clades.json, clade_similarity.json, members.json, distances.npz, commands.log), `tight_clades_summary.json`, `clade_summary.tsv`, `alignable_clades.tsv`, `singletons.tsv`


## Comparison with v1

v1 (all 10,438 geNomad prophages over 7,352 NTM assemblies): 913 clades,
472 alignable, 441 singletons, median internal mash 0.074.

v2 (8,502 collaborator-manifest prophages, NCBI-only, GCA/GCF-deduped):
2,388 clades, 1,251 alignable, 1,137 singletons, median internal mash 0.006.

The v2 median is lower and the singleton fraction higher than v1. This is
expected for the collaborator manifest subset: it is dominated by highly
redundant NCBI assemblies of a few species (e.g. 2,634 M. abscessus
assemblies), so identical/near-identical prophages recur across isolates
(736 clades with internal median 0.0) while the long tail of singletons is
larger than v1's. Triangle values were read back against the mash dist TSV
(50 sampled pairs, exact match) and `mash triangle` exits 0, so the numbers
reflect the data, not a pipeline change (identical recipe: sketch
`-i -k 21 -s 10000`, threshold 0.25, max-size 100, community 0).
