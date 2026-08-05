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
  within `--threshold` (FINAL user-approved **0.25** mash ≈ ~75%+ ANI), else
  becomes a new leader (singleton clade if nothing is close).
- **size cap** `--max-size 100` (FINAL user-approved) — clades that fill up
  split.
- **tightening pass**: any clade whose internal **median pairwise** mash
  distance still exceeds the threshold is split (tighter leader threshold
  `0.6×` with a distance-to-representative half-split fallback), so every
  non-singleton clade has median pairwise mash ≤ 0.25.

Result (`research/clades/<community>/`):

| file | contents |
|---|---|
| `tight_clades.json` | clade_id (`<c>_<NNNN>`) → member list, representative first |
| `clade_similarity.json` | per-clade internal pairwise mash stats (median/min/max, frac ≤ threshold) + community report |
| `members.json` | community members in triangle row order |
| `commands.log` | reproducible commands |

Summary: **585 clades** over 19,638 community members (106 singletons), max
clade size 100. All 479 non-singleton clades have median pairwise mash
≤ 0.25 (max observed median 0.244); every member is within 0.25 of its
clade representative. (The earlier 0.10/800 run produced 1,245 mostly-tiny
clades — median size 3, 38% singletons — and was rejected by the user as
too small for ML/ancestral reconstruction; the 0.25/100 clades are the final
definition.)

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

## Data availability (large inputs/outputs are NOT in git)

- `full_prophages.fa` (3.3 GB) lives in `prophage_homology_survey/` and is
  symlinked into this checkout (`prophage_homology_survey/full_prophages.fa`,
  ignored via `prophage_homology_survey/.gitignore`).
- The 35 GB MASH triangle `research/mash_tree/full_prophages_mash.dist` +
  `research/mash_tree/data/ids.txt` live in the main repo
  (`/home/erikg/phind/...`); `build_tight_clades.py` reads them by absolute
  path, so it works from any worktree.
- Per-clade derived outputs (`sequences.fa`, `*.paf`, `partitions.bed`,
  `partitions/*.maf`) are git-ignored; the last full run's outputs are on disk
  under `research/clades/<c>/<clade_id>/` (validate with
  `python3 scripts/validate_clades.py`).

## Reproduce

```bash
# 1. clade definitions (reads the 35 GB triangle; ~25 s)
python3 scripts/build_tight_clades.py --threshold 0.25 --max-size 100 \
    --communities 0,1,2,3,4,5,6,7,8,9,10,11

# 2. offset index + per-clade pipeline (extract → allwave → segment → partition)
bash scripts/run_all_clade_pipelines.sh

# 3. validate every clade (manifest presence, internal similarity,
#    sparsification, partition size distribution, alignment rate)
python3 scripts/validate_clades.py
```

## References

- MASH triangle & tree: `research/mash_tree/README.md`
- Partition traversal / ML genome: `scripts/traverse_partitions.py`,
  `research/stitching/stitch_algorithm.py`
- Downstream release: task `ml-phage-genome`
