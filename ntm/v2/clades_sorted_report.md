# NTM v2 — corrected clade counts (sorted ids.txt re-derivation)

Generated: 2026-09-13 (task verify-v2-clade-order)

## Context — why the published v2 clade count is inflated

`build_tight_clades.py` sorts community members with `sorted(set(members))`
and its `read_community_matrix` fills `D[i][j]` only for member pairs whose
member-list order agrees with triangle row order. The v2 pipeline
(`ntm/v2/scripts/full_ntm_mash_clades.py`) wrote `ids.txt` in FASTA header
order while the float32 triangle rows follow that same FASTA order, so the
string-sorted member list disagreed with triangle row order for a large
fraction of pairs. Measured on the frozen v2 outputs: **12,730,123 of
36,137,751 within-community pairs (35.23%) were never filled (NaN)** in
`clades/0/distances.npz`. NaN compares `False` against every threshold, so
those pairs silently behaved as maximally distant: missed joins → fragmented
(inflated) clade counts, and an artificially low median internal distance
(clades could only survive tightening if all their internal pairs happened to
be filled). v3 fixed this upstream by writing `ids.txt` pre-sorted
(`ntm/v3/mash_clades_report.md`, "Mechanics note"). The published v2 clade
count is therefore an **upper bound**.

## Method

Re-derive the counterfactual "v2 as if `ids.txt` had been written pre-sorted"
without touching the frozen v2 outputs:

- a bare sorted `ids.txt` would **misindex** the frozen triangle (triangle row
  *k* is the *k*-th FASTA-order id, not the *k*-th sorted id), so the driver
  also **permutes the frozen triangle into sorted-id row order** — a lossless
  reordering of the same float32 values (MASH pairwise distances do not depend
  on sketch input order, so this is exactly the triangle `mash` would have
  produced over a pre-sorted FASTA; v3's approach);
- then runs `scripts/build_tight_clades.py` **unmodified**, same flags as the
  original v2 run (`--threshold 0.25 --max-size 100 --communities 0`);
- everything goes into NEW outdirs (`mash_clades_sorted/`, `clades_sorted/`);
  the frozen `mash_clades/` and `clades/` are only ever read.

Reproduce:

```bash
python3 ntm/scripts/rederive_sorted_clades.py --version v2
```

## Results

| | published (fragmented) | corrected (sorted ids) | delta |
|---|--:|--:|--:|
| clade count | **2,388** | **767** | **−1,621 (−68%)** |
| alignable (≥2 members) | 1,251 | 413 | −838 |
| singletons | 1,137 | 354 | −783 |
| median internal mash (per-clade median of medians) | 0.0060 | 0.0639 | +0.058 |
| max clade size | 100 | 100 | — |
| size distribution (1 / 2-10 / 11-50 / 51-200) | 1,137 / 1,106 / 139 / 6 | 354 / 281 / 78 / 54 | — |

Assignment still total: clade members sum to 8,502 = 8,502 prophages, every
prophage in exactly one clade.

### Why the delta

- With all pairs filled, the degree ordering and nearest-leader assignment see
  the 35.23% of pairs that were previously NaN: **814 of the 1,137 published
  "singletons" join a ≥2-member clade** (323 remain true singletons).
- 270 corrected clades unite ≥2 published clades; on average each corrected
  clade absorbs 3.81 published clades. The biggest merges are dramatic: e.g.
  corrected `0_0025` (100 members, at the size cap) unites **50** published
  clades (16 of them published singletons).
- The median internal distance rises 0.0060 → 0.0639 because previously
  unjoinable (NaN) relatives now merge — the published 0.0060 reflected only
  clades that had survived tightening with a fully-filled pair set, a strongly
  biased subsample.
- 390 published clades straddle >1 corrected clade: the greedy leader order
  changes when degrees are complete, so the corrected partition is a re-cluster
  (same algorithm, complete matrix), not a pure merge of the old one. Both
  partitions satisfy every stated invariant (every member within threshold of
  its representative; every non-singleton median ≤ 0.25; size ≤ 100).

## Validation

- corrected community matrix has **0 NaN** pairs (published matrix: 12,730,123
  of 36,137,751 = 35.23%)
- wherever the published matrix was filled, the corrected matrix is **bitwise
  identical** (proves the permuted triangle reuses the frozen float32 bytes;
  member order verified equal)
- TSV spot check: 5 random pairs re-read from the frozen
  `prophages.dist.tsv` text output vs the permuted triangle — worst |Δ| 4.1e-10
- clade invariants: sum of clade sizes == 8,502; all ids unique; max size 100;
  worst non-singleton internal median 0.2488 ≤ 0.25 (per-clade medians
  recomputed independently from the matrix, not trusted from the run's JSON)
- frozen outputs untouched: only new outdirs written
  (`mash_clades_sorted/`, `clades_sorted/`); md5/mtimes of the frozen
  `ids.txt` / triangle / `tight_clades.json` recorded before and after —
  unchanged
- permutation math covered by a built-in synthetic selftest
  (`rederive_sorted_clades.py --selftest`)
- full machine-readable comparison: `clades_sorted/sorted_vs_published.json`

## Downstream impact (v2 artifacts that consumed the clade count)

These frozen v2 artifacts consumed the **fragmented** counts and therefore
overstate clade/ML-genome counts; they remain accurate as records of what was
run, but their counts should be read as inflated upper bounds (corrected
numbers above, or use v3+ outputs which have the fix):

- `ntm/v2/mash_clades_report.md` — published 2,388 / 1,251 / 1,137, median
  0.0060
- `ntm/v2/ml_report.md`, `ntm/v2/ml_validation_report.md` — ML genomes 2,388
  (= alignable 1,251 + singletons 1,137), ancestral genomes 1,251; with
  corrected clades these would be 767 ML / 413 ancestral genomes (814
  published singletons would instead be represented inside multi-member
  clades)
- `ntm/v2/partition_report.md` — "all 2388 clades" partitioned (86,587
  partitions; 1,251 alignable with `allwave.paf`)
- `ntm/v2/release/RELEASE.md`, `ntm/v2/release/v1_v2_comparison.md` — release
  headline "2,388 clades (1,251 alignable + 1,137 singletons)", ML 2,388 /
  ancestral 1,251, and the v1→v2 "×2.6" clade-count ratio (corrected: 913 →
  767 on different prophage sets; see `ntm/v1/clades_sorted_report.md`)

No v2 downstream artifact is modified by this task; per-prophage release
joins (prophage → accession → host clade) are unaffected because they key on
prophage id, not clade id.

## Outputs

- NVMe `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/mash_clades_sorted/`:
  `ids.txt` (sorted), `labels.csv` (matching order), `prophages_mash_sorted.dist`
  (permuted float32 triangle), `provenance.json` (sha256 of frozen + sorted
  triangles)
- NVMe `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/clades_sorted/`:
  `0/` (tight_clades.json, clade_similarity.json, members.json,
  distances.npz, commands.log), `tight_clades_summary.json`,
  `clade_summary.tsv`, `alignable_clades.tsv`, `singletons.tsv`,
  `sorted_vs_published.json`, `build_tight_clades.stdout.log`
- repo: this report + driver `ntm/scripts/rederive_sorted_clades.py`
