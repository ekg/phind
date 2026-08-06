# Full 26k E. coli host MASH tree

A MASH distance triangle and UPGMA tree over all **26,074 E. coli assembly
accessions** of the full 26k PanSN cohort — the host phylogeny that pairs with
the prophage tree `research/mash_tree/full_prophages_tree.nwk` (132,393
leaves). Leaf names are the assembly accessions (`GCF_...`), exactly the host
prefix of the `ACC_prophage_N` prophage leaf names, so the two trees can be
joined leaf-for-leaf by accession for cophylogeny analysis.

## Files

| File | Size | Contents |
|---|---|---|
| `hosts.msh` | 2.1 GB | MASH sketch of all 26,074 genomes (`-i`, k=21, s=10000, seed 42) |
| `host_mash.dist` | 1.36 GB | float32 upper triangle, n·(n−1)/2 = 339,913,701 values (format below) |
| `host_tree.nwk` | 0.9 MB | UPGMA tree, 26,074 leaves, Newick with branch lengths |
| `ids.txt` | 0.4 MB | assembly accessions, one per line; line index = triangle row/col index |
| `triangle_verify.json` | — | triangle size/spot-check verification |
| `tree_verify.json` | — | tree leaf/ultrametric/accession verification |
| `host_tree_stats.json` | — | tree leaf/internal counts, height, branch stats |
| `data/` (git-ignored) | ~129 GB | chunk fastas, chunk sketches, part files, logs |
| `scripts/` | — | full pipeline scripts |
| `COMMANDS.log` | — | reproducible command log |

Large derived files (`hosts.msh`, `host_mash.dist`, `data/`) are git-ignored
and regenerable with `bash run_full.sh`.

## Triangle format (`host_mash.dist`)

Little-endian **float32**, length n·(n−1)/2 = 339,913,701 values, stored
row-major over the upper triangle (diagonal excluded). For rows/cols indexed
0..n−1 in `ids.txt` order, the offset of pair (a,b) with a<b is:

```
offset(a,b) = a·(2n−a−1)/2 + (b−a−1)
```

## Method

1. **Inputs.** 26,074 uncalled E. coli assemblies at
   `/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/26k/canonical_objects/*/*.pansn.fa.gz`
   (PanSN bgzip; header `ACC#hap#chr`). Each genome is concatenated into one
   record named by the assembly accession. 26,074 unique accessions; the three
   `canonical_objects` subdirs without a `.pansn.fa.gz` file
   (GCF_001291365.1, GCF_003886435.1, GCF_012029685.1) are excluded by the
   glob.

2. **Sketch + pairwise triangle.** Mirror of the prophage `mash_tree` method
   (k=21, seed 42) with **s=10000** per the host-structure-1000-v1
   methodology. 6 chunks of 5000 genomes → 21 chunk-pair jobs, 8 mash
   processes × 32 threads, per-job `part_writer` part files to avoid the
   shared-file write pathology, then merged into `host_mash.dist`.

3. **UPGMA tree.** scipy average linkage on the float64 condensed matrix,
   converted to Newick with branch lengths = half the node-height differences
   so leaf-to-leaf path length equals the UPGMA cophenetic distance
   (= tree height = ultrametric).

## Verification summary (`triangle_verify.json`, `tree_verify.json`)

- **Triangle exact size**: 339,913,701 values × 4 bytes = 1,359,654,804 bytes
  (`size_ok=true`).
- **50-pair spot check** vs direct `mash dist` on re-sketched genomes:
  **0 mismatches**.
- **Tree leaves** = 26,074 == sketched count (`leaf_count_ok=true`); all names
  unique and equal the `ids.txt` set.
- **Ultrametric**: root-to-leaf distance range 0.130229–0.130241 (tol < 1e-4)
  (`ultrametric_tol_1e4=true`).
- **Accession names**: all 26,074 leaves are `GCF_`/`GCA_` assembly accessions;
  the host accessions **exactly equal** the prophage→host prefix set
  (26,074/26,074 found among prophage leaf host prefixes) — leaf-for-leaf
  joinable with the prophage tree.
- **Cophenetic vs source**: R² = 0.96 on 300 full-range pairs, mean abs err
  0.001, max abs err 0.007 (matrix nearly ultrametric).

## Reproducibility

```bash
cd research/host_mash_tree
# edits: set MASH + SRC in run_full.sh
bash run_full.sh
```

Environment: Mash v2.3, Python 3.12 (numpy, scipy, ete3), g++ (part writer);
study machine has 256 cores / ~1 TB RAM / NVMe.

## Notes / limitations

- **No saturated pairs at 1.0** (all genomes share many 21-mers at s=10000) —
  unlike the saturated prophage matrix. Distances cluster tightly near mean
  0.024, median 0.025 (range ~0–0.26).
- **346 zero-distance pairs** in the triangle — likely near/exact duplicate
  strains (they form zero-length clades, kept as leaves).
- UPGMA assumes an approximately ultrametric signal (molecular clock); branch
  lengths are mean-merge distances, not divergence times.
