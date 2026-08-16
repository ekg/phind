# NTM v2 — per-clade allwave + segment + impg partition report

Generated: 2026-08-16T14:47:00Z (task ntm-v2-per)

## Run

Command (via `ntm/scripts/resume_per_clade_partition.py`, which re-invokes
`scripts/per_clade_alignment_pipeline.py` with built-in resume):

```
python3 scripts/per_clade_alignment_pipeline.py \
  --community 0 \
  --outdir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/clades \
  --clades  /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/clades/0/tight_clades.json \
  --fasta   /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/full_prophages.fa \
  --index   /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/mash_clades/full_prophages.idx.json \
  --threads 8 --jobs 16 --all
```

Per clade: extract `sequences.fa` (seek-based, offset index) -> **allwave**
-> `scripts/segment_paf.py` (window 500, gap 0, max-span 1000) -> **impg
partition** (`-w 500 -d 0 --min-boundary-distance 0 --min-missing-size 0
-m 1 --no-rehome-singletons`; BED + per-partition MAF).

Allwave parameters (same as v1 mainline):
- n <= 30: `-p none` (all-pairs)
- 30 < n <= 200: `-p tree:5:0:0.0` (k-nearest 5, k-farthest 0, no
  stranger-joining; v2 has no clade > 200, so `tree:10:0:0.0` never used)
- scores: allwave default `0,5,8,2,24,1`

Runtime: **2008 s** (~33.5 min) wall for all 2388 clades, 16 concurrent
clade workers x 8 threads (128 threads; cap on concurrent allwave
processes, v1 practice). Log: `clades/partition_run.log`; per-clade
`commands.log` + `manifest.json` for reproducibility.

## Coverage

- Total tight clades: **2388**
- Alignable (n>=2): **1251** — all with `allwave.paf`,
  `allwave.segmented.paf`, `partitions.bed`, `partitions/partition<N>.maf`
- Singletons: **1137** — `sequences.fa` + empty `allwave.paf` (v1
  convention; pass-through for the ML task)
- Failures: **0** (no clade skipped; no FAILED.log)

## Partition stats

- Total partitions (distinct blocks): **86,587**
- Total partition intervals: **425,090**
- Median partition length: **500 bp** per-clade (median of per-clade
  medians 500 bp; overall interval median 500 bp) — matches v1 (500 bp)
- Overall interval distribution: mean 389 bp, min 1, max 2134;
  0.13% > 1000 bp (544), 12.2% < 100 bp
- Strategy distribution (alignable): `tree:5:0:0.0` x 23 (n>30),
  `none` (all-pairs) x 1228

## Outputs

`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/clades/`:
- `0_0000` .. `0_2387`: per-clade dirs (sequences.fa, allwave.paf,
  allwave.segmented.paf, partitions.bed, partitions/, manifest.json,
  commands.log)
- `pipeline_results.json` (per-clade manifests, aggregate)
- `partition_run.log`, `partition_report.md`

## Comparison with v1

| | v1 (10,438 prophages) | v2 (8,502 prophages) |
|---|---|---|
| alignable clades | 472 | 1251 |
| singletons | 441 | 1137 |
| total partitions | 292,365 | 86,587 |
| median partition length | 500 bp | 500 bp |
| failures | 0 | 0 |

v2 has 2.6x the alignable clades but far fewer large ones (v2 max clade
size 100; only 6 clades of 100 members vs v1's 74 clades in 51-200), so
total partitions are fewer while the per-partition median stays 500 bp.
