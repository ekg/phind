# All-prophage MASH tree

A MASH distance triangle and UPGMA tree over **all 132,393 prophage elements**
of `full_prophages.fa` (3.1 GB, headers `GCF_..._prophage_N`, from the prophage
homology survey), with cluster/community/MDS labels merged from the existing
survey CSVs.

## Files

| File | Size | Contents |
|---|---|---|
| `full_prophages.msh` | 1.1 GB | MASH sketch of all 132,393 prophages (`-i`, k=21, s=1000, seed 42) |
| `full_prophages_mash.dist` | 35.1 GB | float32 upper triangle, n·(n−1)/2 values (see format below) |
| `full_prophages_tree.nwk` | ~4 MB | UPGMA tree, 132,393 leaves, Newick with branch lengths |
| `full_prophages_tree_labeled.nwk` | ~5 MB | same tree, leaf names annotated `..._N\|C<community>` |
| `full_prophages_labels.csv` | ~4 MB | per-leaf: sequence, community, cluster, genome, MDS1, MDS2 |
| `data/ids.txt` | 1.4 MB | sequence ids, one per line; line index = triangle row/col index |
| `data/chunks/` | ~1.1 GB | 27 chunk FASTA+sketches (5000 seqs each, last 2393) |
| `data/parts/` | ~37 GB | per-job part files (float32 in triangle order); deleted after merge |
| `scripts/` | — | all pipeline scripts (see below) |
| `tree_stats.json` | — | tree leaf/internal counts, height, branch stats |
| `label_merge_stats.json` | — | label merge rates |
| `triangle_verify.json` | — | triangle verification results |

Large derived files (`full_prophages.msh`, `full_prophages_mash.dist`,
`data/`) are git-ignored and regenerable with the commands below.

### Triangle format (`full_prophages_mash.dist`)

Little-endian **float32**, length n·(n−1)/2 = 8,763,887,028 values, stored
row-major over the upper triangle (diagonal excluded). For rows/cols indexed
0..n−1 in `data/ids.txt` order, the offset of pair (a,b) with a<b is:

```
offset(a,b) = a·(2n−a−1)/2 + (b−a−1)
```

A read-back example is in `scripts/verify_triangle.py` (numpy memmap).

## Method

1. **Sketch.** `mash sketch -i -k 21 -s 1000 -p 256 -o full_prophages.msh
   full_prophages.fa` (Mash v2.3; MurmurHash3_x64_128, seed 42, canonical
   ACGT alphabet). All 132,393 sequences were ≥ 261 bp, so **none were
   dropped** by mash (0 dropped elements; sketch contains 132,393 sketches —
   `mash info` confirms). Default s=1000 was used (mean prophage 24.6 kb →
   Jaccard std.err ≈ 0.012 at the mean distance 0.18).

2. **Pairwise triangle.** The full all-vs-all is 17.5e9 ordered pairs
   (8.76e9 unordered). The file was computed in 27 chunks of 5000 sequences:
   `mash dist -p 32 chunk_j.msh chunk_i.msh` (i<j) emits each chunk-pair block
   with ref=later chunk so the value stream is already in triangle row order;
   diagonal blocks are buffered and transposed by the writer. Eight concurrent
   mash processes feed per-job `part_writer` processes writing **separate**
   part files — essential: we measured ~8× wall-time degradation when ≥8
   writers pwrite concurrently into one 35 GB file (411 s vs 62 s for 8 jobs;
   the shared-file pathology persisted even with `fallocate` preallocation).
   Part files are then merged into `full_prophages_mash.dist` with contiguous
   per-row writes (single process, ~5 min).

   **Validation:** on a 2,000-sequence pilot, all 1,999,000 triangle entries
   matched direct `mash dist` output exactly (0 mismatches), including runs
   with unequal chunk sizes; the full file is verified in
   `triangle_verify.json` (size, row-0 range, zero-distance pairs, sample
   statistics, and a 50-pair spot-check vs direct mash).

3. **Tree.** UPGMA via `scipy.cluster.hierarchy.linkage(method="average")` on
   the float64 condensed matrix, converted to Newick with branch lengths =
   half the node-height differences so that the sum of branch lengths between
   two leaves equals their UPGMA merge height (= the tree cophenetic distance,
   matching the source matrix scale). 132,393 leaves / 132,392 internal nodes.

   **Why UPGMA and not NJ:** neighbor-joining is O(n³) naively (2.3e15 ops at
   n=132,393). RapidNJ (O(n²), heap-based) was benchmarked on the pilot:
   it fitted the source distances *worse* than UPGMA on this data
   (cophenetic R² = −0.09 vs 0.45; mean |err| 0.31 vs 0.20) — the matrix is
   dominated by saturated distances (≈70% of pairs at exactly 1.0, no shared
   k-mer) plus 30,125 exact-duplicate copies (10,366 identical-sequence
   groups), which NJ handles poorly (negative-length artifacts) while UPGMA's
   ultrametric structure correctly collapses duplicate prophages into
   zero-length clades. UPGMA is also far cheaper at this scale. The trade-off
   is the molecular-clock assumption (documented limitation).

4. **Labels.** `full_heatmap_clusters.csv` (all 132,393 prophages) assigns a
   community id per sequence: **12 shared communities (19,638 prophages)** —
   Louvain communities over mash-distance edges (d < 0.5, weight 1−d, from
   the survey's `full_heatmap_generator.py`) — plus **112,755 singleton
   communities** (isolated prophages with no close relative; each has its own
   id). `full_prophage_clusters.csv` (428 clusters) and
   `full_prophage_mds_coords.csv` cover a 5,000-prophage subset. All three
   are merged onto the tree leaves (`merge_labels.py`): per-leaf
   community/cluster/genome/MDS in `full_prophages_labels.csv`, and a
   community-annotated Newick where leaves of the 12 shared communities are
   tagged `..._N|C<id>` and isolates `..._N|iso`.
   Merge rates in `label_merge_stats.json`: **community coverage 100%**
   (≥95% validation met), cluster 3.8%, MDS 3.8%.

   **Tree/community agreement** (`tree_verify.json`): among the 19,638
   non-isolate leaves, sibling-clade community purity is **mean 0.992,
   median 1.000**, and mean within-community tree distance (0.61) is well
   below between-community (0.88) — the tree strongly recovers the shared
   communities. The 112,755 isolates are, by construction, singletons that
   share no close sequence similarity with anything (they attach near the
   tree root at saturated distances).

## Reproducibility

Environment: Mash v2.3, Python 3.12 (numpy, scipy, ete3), g++ (part writer),
RapidNJ 2.3.3 (bioconda, used only for the NJ comparison), 256 cores / ~1 TB
RAM / NVMe (machine is shared; timings below are wall-clock under load).

```bash
cd research/mash_tree
# 1. sketch (all 132,393 seqs, none dropped)
mash sketch -i -k 21 -s 1000 -p 256 -o full_prophages.msh \
    /home/erikg/phind/prophage_homology_survey/full_prophages.fa

# 2. split into 27 chunks, sketch chunks, compute part files, merge triangle
python3 scripts/run_pairwise.py \
    --fasta /home/erikg/phind/prophage_homology_survey/full_prophages.fa \
    --workdir data --out full_prophages_mash.dist --ids data/ids.txt \
    --chunk-size 5000 --procs 8 --threads 32

# 3. verify triangle
python3 scripts/verify_triangle.py

# 4. UPGMA tree -> full_prophages_tree.nwk (+ tree_stats.json)
python3 scripts/build_tree.py

# 5. merge labels -> labels.csv, labeled tree, merge-rate report
python3 scripts/merge_labels.py
```

## Timing / resource report (honest)

| Step | Wall time | Notes |
|---|---|---|
| mash sketch (full) | 1 min 15 s | `-p 256`; 1.08 GB output |
| pairwise (378 chunk-pair jobs, 8 procs × 32 threads) | 45 min 37 s | 17.5e9 ordered pairs computed; 8-way process concurrency (mash `-p` does not scale past ~32 threads); per-job part files |
| part-file merge | 1 min 30 s | single process, contiguous per-row writes, 8.76e9 values |
| UPGMA linkage (scipy) | 23 min 47 s | 132,393 leaves, peak RSS ~137 GB (condensed float64 70 GB + scipy working set); single-threaded |
| Newick emission + stats | < 1 s | iterative stack traversal |
| triangle verify | ~2 min | 50-pair spot check vs direct mash: 0 mismatches |
| tree validate | ~3 min | ete3 parse + checks (see tree_verify.json) |
| label merge | < 1 s | 100% community coverage |

Peak machine usage: 8 mash processes × 32 OpenMP threads + 8 writer processes
(~16 cores under load); tree build peaked at ~137 GB RSS. Input file
3.1 GB; intermediate part files ~37 GB (deleted after merge); outputs
35.1 GB (triangle) + ~20 MB (tree, labels, reports).

## Known limitations

- **UPGMA assumes an approximately ultrametric signal** (molecular clock);
  branch lengths are mean-merge distances, not divergence times. NJ was
  benchmarked and performed worse on this saturated-distance data (see above).
- **Saturated distances:** ~70% of pairs have distance exactly 1.0 (no shared
  21-mer). These are uninformative for fine structure; the tree mainly
  resolves the ~30% of pairs that share k-mers (including all
  exact-duplicate groups, distance 0).
- **Exact duplicates:** 30,125 prophages are exact sequence copies of others
  (10,366 groups, up to 135 members — the same prophage in different E. coli
  strains). They form zero-length clades in the tree; all copies are retained
  as leaves (leaf count == sketched count).
- **Community labels:** only 19,638 prophages belong to the 12 shared
  communities; the other 112,755 are isolates (each its own community id).
  The tree recovers the shared communities almost perfectly (clade purity
  0.99), and isolates attach at saturated distances near the root.
- **Sketch size s=1000** is the mash default; finer resolution would need a
  larger sketch at ~10× pairwise cost.
