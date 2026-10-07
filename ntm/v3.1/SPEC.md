# NTM v3.1 — method revision: structure-driven clustering + genome-path reconstruction

**Goal.** Fix the two defects measured in v3, prove it on a subset, and only then
decide whether to re-run the full pipeline. v3 is **not** modified; v3.1 builds
alongside so both can be compared.

## Defects being fixed (both measured, see `ntm/v3/*validation*`)

1. **Arbitrary size cap in clustering.** `build_tight_clades.py --max-size 100`
   is a count, not a structure criterion. Uncapped leader clustering at MASH 0.25
   yields **967** clades (largest **4,278**, 81.7% of prophages in clades >100,
   median size 2). The cap splits only the ~50 largest into 1,304, and siblings
   are near-identical (min cross-clade distance **0.0000**) → duplicate ML genomes.
2. **Chimeric reconstruction.** `traverse_partitions.py` treats partitions as
   independent tiles and only *biases* toward adjacency (`1 + beta*adj`); it never
   *requires* an observed adjacency. So the walk accumulates high-support blocks
   until `--max-length 150000` stops it. Result: 57% of ML genomes are >1.5× the
   median member length (51–99-member clades: **16.6×**), CheckV median
   completeness **49.7%**, 532/1,304 low-quality.

## Changes (all backward-compatible; v3 defaults unchanged)

### A. Clustering — `scripts/build_tight_clades.py`
Add `--split-mode {cap,medoids}` (default `cap` = current behaviour) and
`--medoid-radius D` (default 0.10):
- `cap`: unchanged.
- `medoids`: leader clustering at `--threshold` with **no size cap**, then
  recursively split any cluster whose members exceed `D` from its medoid
  (k-medoids, k=2 per split, until every member is within `D`).
  Split points are then structure-driven, not count-driven.
- Emit the same `tight_clades.json` schema + a `split_trace.json`
  (per-split: parent size, radius before/after, iterations).

### B. Reconstruction — `scripts/traverse_partitions.py`
Add `--path-mode {free,observed}` (default `free` = current) and `--close-circle`
(default off):
- `observed`: candidate set restricted to partitions with an **observed
  adjacency** from the current partition (`adj[current][next] > 0`); if none,
  the walk terminates (natural terminus). First partition that is *excluded*; prefer
  `first_counts`. `--close-circle`: stop early when a candidate is observed to be
  adjacent back to the path's start (circular phage).
- Length budget stays as a **safety cap only**; in `observed` mode the natural
  terminus should be reached well before it.
- Emit extra per-clade stats: `n_path_steps`, `path_adj_fraction`, `terminated_by`
  (`terminus` | `closed` | `budget`), `budget_hit` (bool).

### C. Quality scoring
For every v3.1 ML genome: CheckV completeness/quality, plus derived flags
`budget_hit`, `duplicate_of` (nearest other clade within MASH 0.05),
`n_members`. Emit a `quality.tsv` and a `quality_summary` with the fraction
high-quality — the number that decides whether v3.1 wins.

### D. Comparison + validation protocol
**Subset first (mandatory).** Run v3.1 clustering + reconstruction on:
- the largest v3 natural cluster (~4,278 members) — the worst case,
- 2 mid-size clusters (200–1,500 members),
- 2 small v3 clades (11–50),
- 5 singletons.
Then CheckV the resulting genomes and compare against the same region in v3:
| metric | v3 | v3.1 | win condition |
|---|---|---|---|
| ML genomes produced for the subset | | | fewer duplicates at equal coverage |
| % high-quality (CheckV) | | | **higher** |
| median completeness | | | **higher** |
| ML length / member median length | ~5–17× | | **→ ~1×** |
| clades hitting the length budget | high | | **→ ~0** |

**Full re-run only if** the subset shows higher high-quality fraction AND
length ratio near 1. Otherwise stop and report.

## Paths
- repo: `ntm/v3.1/` (specs, reports, comparison)
- NVMe: `$NVME/ntm/v3.1/` ({clades,ml,release}; never touches `$NVME/ntm/v3/`)
- original v3 artifacts are read-only inputs

## Non-goals
- No change to prophage calling, extraction, host clades, or the annotation stack.
- No change to v3 outputs on disk.
