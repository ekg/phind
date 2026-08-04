# Tight-clade all-wave alignment + small partitions (per-clade-fastga)

Per-clade MASH-tightened clades, all-wave alignment (biWFA via `allwave`),
CIGAR-aware PAF segmentation, and small-window partitioning for ALL 12 old
communities (0–11). This is the clade definition + alignment input for the
downstream ML-phage-genome release and the partition traversal.

## Why tight clades (measured problem)

The old 12-community labels are internally incoherent (measured from
`research/mash_tree/full_prophages_mash.dist`):

| community | n seqs | % pairs at mash dist 1.0 | median pairwise dist |
|---|---:|---:|---:|
| 0 | 4999 | 42.5% | 0.263 |
| 1 | 2128 | 19.8% | 0.120 |
| 2 | 2826 | 60.4% | 1.000 |
| 3 | 1234 | 21.0% | 0.178 |
| 4 | 3030 | 71.7% | 1.000 |
| 10 | 5208 | 33.4% | 0.220 |

Aligning all pairs within these communities (e.g. with stranger-joining
`tree:5:2:…`) wastes days of compute on unrelated sequences. Therefore:
clades are tightened by MASH distance BEFORE any alignment, and
stranger-joining (k-farthest > 0) is NEVER used.

## Clade definition (step 1)

`scripts/build_tight_clades.py` reads the precomputed all-prophage MASH
triangle (`research/mash_tree/full_prophages_mash.dist`, float32 upper
triangle in `research/mash_tree/data/ids.txt` order) and sub-clusters each of
the 12 communities:

- **greedy leader clustering**: members processed in decreasing
  within-community neighbor count; each member joins the nearest leader
  within `--threshold` (default **0.10** mash ≈ ~90% ANI), else becomes a
  new leader (singleton clade if nothing is close).
- **size cap** `--max-size 800` — clades that fill up split.
- **tightening pass**: any clade whose internal **median pairwise** mash
  distance still exceeds the threshold is split (tighter leader threshold
  `0.6×` with a distance-to-representative half-split fallback), so every
  non-singleton clade has median pairwise mash ≤ 0.10.

Result (`research/clades/<community>/`):

| file | contents |
|---|---|
| `tight_clades.json` | clade_id (`<c>_<NNNN>`) → member list, representative first |
| `clade_similarity.json` | per-clade internal pairwise mash stats (median/min/max, frac ≤ threshold) + community report |
| `members.json` | community members in triangle row order |
| `commands.log` | reproducible commands |

Summary: **1245 clades** over 19,638 community members (469 singletons),
max clade size 800. All 776 non-singleton clades have median pairwise
mash ≤ 0.10 (max observed median 0.0995); every member is within 0.10 of
its clade representative.

## Alignment (step 3)

Per clade, `allwave` (biWFA, k-mer-tree sparsification) with
**k-farthest = 0 everywhere — never stranger-joining**:

| clade size | sparsification | pair count |
|---|---:|---|
| n = 1 | none (no alignment) | 0 |
| 2 ≤ n ≤ 30 | `-p none` (all-pairs) | n(n−1)/2 |
| 31 ≤ n ≤ 200 | `-p tree:5:0:0.0` | ~5n directed |
| n > 200 | `-p tree:10:0:0.0` | ~10n directed (k=10 keeps the k-nearest graph connected at linear cost) |

Alignment scores default `0,5,8,2,24,1` (85–95% ANI). Runtime, pair count,
alignment rate and connectivity (`sequences_in_paf`) are recorded per clade
in `manifest.json`.

## Segmentation (step 4)

`scripts/segment_paf.py` — CIGAR-aware, window **500 bp**, max chunk span
**1000 bp (2×window)** on BOTH query and target, contiguous chunks (no
inter-chunk gaps).

## Partitioning (step 5)

`impg partition` per clade with
`-w 500 -d 0 --min-boundary-distance 0 --min-missing-size 0 -m 1
--no-rehome-singletons`:

- `partitions.bed` — combined 4-col BED (prophage, start, end, partition_id)
- `partitions/partition<N>.maf` — per-partition alignment blocks (MAF)
  (identical partition ids as the BED — rehoming disabled so the two outputs
  are consistent)

`-m 1` (max transitive depth 1) keeps partition intervals at ~500 bp
(chunk-aligned) instead of letting transitive gathering merge long
contiguous stretches. Partition intervals are typically 500 bp; the observed
max is ~1.7–1.8 kb on conserved blocks shared by many sequences (two merged
500 bp chunks) and on long private stretches in the most divergent clades
(e.g. c9, max pairwise 0.167). Those are real homologous/private modules
with no internal split evidence — documented per clade in `manifest.json` →
`partition_bed.interval_stats` (`n_gt_1000`).

## Outputs per clade

```
research/clades/<community>/<clade_id>/
    sequences.fa              extracted clade FASTA (byte-identical to source)
    allwave.paf               all-wave alignment
    allwave.segmented.paf     segmented PAF
    partitions.bed            combined partition assignments
    partitions/partition*.maf per-partition alignments
    manifest.json             members, similarity stats, strategy, params,
                              partition size distribution, alignment rate, runtimes
    commands.log              exact reproducible commands
```

Large derived files (FASTA/PAF/MAF) are git-ignored and regenerable via
`scripts/run_all_clade_pipelines.sh` (commands in each `commands.log`).

## Reproduce

```bash
# 1. clade definitions (reads the 35 GB triangle; ~25 s)
python3 scripts/build_tight_clades.py --communities 0,1,2,3,4,5,6,7,8,9,10,11

# 2. offset index + per-clade pipeline (extract → allwave → segment → partition)
bash scripts/run_all_clade_pipelines.sh
```

## References

- MASH triangle & tree: `research/mash_tree/README.md`
- Partition traversal / ML genome: `scripts/traverse_partitions.py`,
  `research/stitching/stitch_algorithm.py`
- Downstream release: task `ml-phage-genome`
