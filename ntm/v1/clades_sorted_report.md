# NTM v1 — corrected clade counts (sorted ids.txt re-derivation)

Generated: 2026-09-13 (task verify-v2-clade-order; same method as
`ntm/v2/clades_sorted_report.md` — see there for the full mechanics)

## Result

Re-running `build_tight_clades.py` (thr 0.25 / max 100 / community 0,
unmodified) over the frozen v1 triangle with sorted `ids.txt` + matching
`labels.csv` + the triangle permuted into sorted-id row order reproduces the
published v1 numbers **exactly**:

| | published | corrected (sorted ids) |
|---|--:|--:|
| clade count | 913 | **913** |
| alignable (≥2) | 472 | 472 |
| singletons | 441 | 441 |
| median internal mash | 0.0740 | 0.0741 (identical) |
| size distribution | 441 / 315 / 83 / 74 | identical |

The partition is member-identical (same set-of-sets; 910/913 clade ids map to
byte-identical member lists, the remaining 3 are renumbered by creation order).

## Why v1 was (essentially) unaffected

v1's `ids.txt` was also written in FASTA header order, but that order was
already almost lexicographic: only **104 of 54,470,703** within-community
pairs (0.0002%) were NaN in `clades/0/distances.npz` — vs **35.23%** in v2,
whose FASTA order was far from sorted. Those 104 pairs did not change any
clustering decision.

**v1's published numbers (913 / 472 / 441, median 0.074) stand.** The v1→v2
"×2.6 more clades" comparison in `ntm/v2/release/v1_v2_comparison.md` was
driven by v2's fragmentation: corrected, v2 has *fewer* clades than v1 (767 vs
913) on a smaller, more redundant prophage set — see
`ntm/v2/clades_sorted_report.md`.

## Validation

- corrected matrix 0 NaN (published matrix: 104/54,470,703 = 0.0002%)
- bitwise-identical to the published matrix wherever it was filled
- TSV spot check 5 pairs vs frozen `prophages.dist.tsv`: worst |Δ| 0.0
- clade invariants: sum of sizes == 10,438; ids unique; max size 100; worst
  non-singleton median 0.2488 ≤ 0.25 (recomputed from the matrix)
- frozen outputs untouched (new outdirs only; md5s recorded before/after)
- machine-readable comparison: `clades_sorted/sorted_vs_published.json`

## Outputs

- NVMe `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/mash_clades_sorted/`
  (sorted `ids.txt`, matching `labels.csv`, permuted triangle,
  `provenance.json`) and `.../ntm/v1/clades_sorted/` (full build output +
  `sorted_vs_published.json`)
- reproduce: `python3 ntm/scripts/rederive_sorted_clades.py --version v1`

Downstream v1 artifacts that consumed 913/472/441 (e.g. `ntm/RELEASE.md`,
`ntm/v2/release/v1_v2_comparison.md` v1 columns) need **no** correction.
