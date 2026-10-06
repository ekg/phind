# NTM v3 — per-clade allwave + segment + impg partition report

Generated: 2026-09-14T00:45:00Z (task ntm-v3-per)

## Run

Command (via `ntm/scripts/resume_per_clade_partition.py`, which re-invokes
`scripts/per_clade_alignment_pipeline.py` with built-in resume; driver run
detached via `setsid nohup` so worker/provider deaths never lose compute):

```
python3 ntm/scripts/resume_per_clade_partition.py \
  --outdir /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades \
  --clades  /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades/0/tight_clades.json \
  --fasta   /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/full_prophages.fa \
  --index   /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/mash_clades/full_prophages.idx.json \
  --threads 8 --jobs 16
```

Per clade: extract `sequences.fa` (seek-based, offset index
`mash_clades/full_prophages.idx.json`) -> **allwave** ->
`scripts/segment_paf.py` (window 500, gap 0, max-span 1000) -> **impg
partition** (`-w 500 -d 0 --min-boundary-distance 0 --min-missing-size 0
-m 1 --no-rehome-singletons`; BED + per-partition MAF).

Allwave parameters (FINAL, user-approved; identical to v2):
- n <= 30: `-p none` (all-pairs)
- 30 < n <= 200: `-p tree:5:0:0.0` (k-nearest 5, k-farthest 0, no
  stranger-joining)
- n > 200: `-p tree:10:0:0.0` (k-nearest 10, k-farthest 0) — **not used**:
  v3 max clade size is 100 (cap from the clades step)
- scores: allwave default `0,5,8,2,24,1` (proven 85-95% ANI regime)

Runtime: **~70 min wall** for the completing driver pass (16 concurrent
clade workers x 8 threads), after ~40 clades had been completed by earlier
passes under killed workers. Per-clade runtime across all attempts:
total 65,970 s, median 8.3 s, max 3,121.4 s (clade 0_0034, n=100, first
pass under full 128-thread contention). Logs: NVMe
`ntm/v3/clades/driver_attempt5.log`, `clade_0_0034_retry.log`;
per-clade `commands.log` + `manifest.json` for reproducibility.

## Resume proof

- Three earlier attempts were killed mid-run by provider stream timeouts
  (not by pipeline failures); they left 24 complete clade manifests on
  NVMe.
- The completing driver invocation reported `789 incomplete clades (of
  813)` at launch — the 24 pre-existing manifests were detected and
  skipped, and partially-aligned clades were recomputed idempotently
  (allwave/segment outputs overwritten, no stale artifacts reused except
  one impg index, see Incident below).
- Worker death therefore cost zero clade-days of compute.

## Coverage

- Total tight clades: **813** (over 9,446 prophages)
- Alignable (n>=2): **456** — all with non-empty `allwave.paf`,
  `allwave.segmented.paf`, `partitions.bed`, `partitions/partition<N>.maf`
  and `manifest.json` (validator: 456/456 full outputs, 0 missing)
- Singletons: **357** — enumerated in `clades/singletons.tsv` (repo) and
  `ntm/v3/clades/singletons.tsv` (NVMe copy from the prophage task);
  per-clade dirs hold `sequences.fa` + empty `allwave.paf` +
  `manifest.json` with `strategy: "none (singleton, no alignment)"`
  (v1/v2 convention; pass-through for the ML task, no alignment needed)
- Failures: **0** (no clade skipped; no FAILED.log)

## Partition stats

- Total partitions (distinct blocks): **151,237**
- Total partition intervals: **500,584**
- Median partition length: **500 bp** per-clade (median of per-clade
  medians 500 bp; highest per-clade median 853 bp) — matches v1/v2 (500 bp)
- Overall interval distribution: mean 402.9 bp, max 2,049 bp;
  **0.135% > 1000 bp (676 intervals, in 137 of 456 clades; worst clade
  0_0000 with 29)**, 11.10% < 100 bp
- Strategy distribution (alignable): `none` (all-pairs) x 367,
  `tree:5:0:0.0` x 89 (31 <= n <= 100); `tree:10:0:0.0` x 0 (no clade
  > 200)
- Alignment rate (pairs aligned / possible pairs): median **1.0000**,
  min 0.0679 (tree:5 on n=100 clades aligns ~10 directed k-nearest pairs
  per member by design; k-farthest=0, no stranger-joining)
- Member-hit rate (fraction of clade members with >= 1 PAF hit):
  min **1.0000** — every prophage in every alignable clade has at least
  one alignment

### Justification: intervals > 1000 bp

The 676 intervals > 1000 bp (0.135%, max 2,049 bp) are boundary
intervals where impg's transitive closure (`-m 1`) merges adjacent 500 bp
chunks across small unaligned gaps at block edges. This is the same
phenomenon and rate as v2 (544 intervals, 0.13%, max 2,134 bp) under
identical FINAL parameters, where it was accepted; the body of every
partition remains ~500 bp (median 500). Splitting these boundary
intervals would change the approved impg parameter set, so they are
reported and accepted as-is.

## Incident: clade 0_0034 stale impg index (fixed)

The 21:35Z pass under a killed worker left a stale
`allwave.segmented.paf.impg` index in `clades/0_0034/`; all subsequent
retries regenerated the PAF but impg read CIGARs through the old index's
byte offsets and panicked (`Invalid CIGAR operation: \t`), crashing the
pipeline 5x until the driver's attempt budget ran out. Root cause was
external worker death mid-run, not bad data. Fix (two-part):

1. `scripts/per_clade_alignment_pipeline.py` now removes any stale
   `*.segmented.paf.impg` before the bed partition step, so a killed run
   can never poison a retry (commit in this task).
2. Clade 0_0034 was re-run cleanly: n=100, `tree:5:0:0.0`,
   28.2 s, 1,533 MAF partitions, manifest complete.

After the fix the independent validator
(`ntm/v3/scripts/validate_v3_partitions.py`) reports:
`clades total=813 alignable=456 singletons=357; alignable with full
outputs: 456; missing: 0 failed: 0 manifest issues: 0; VALIDATION: PASS`
(parameter audit incl. k_farthest=0/no-stranger, window 500/max-span
1000, scores, per-clade strategy-vs-n, repo==NVMe clades sha).
The full executed validator output is checked in as
`ntm/v3/partition_validation.txt` (exit 0, VALIDATION: PASS), including
the per-clade `>1000bp` breakdown backing the justification above
(137/456 clades affected, global max 2,049 bp, worst clades listed with
their per-clade maxima).

## Outputs

`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/clades/`:
- `0_0000` .. `0_0812`: per-clade dirs (`sequences.fa`, `allwave.paf`,
  `allwave.segmented.paf`, `partitions.bed`, `partitions/partition<N>.maf`,
  `manifest.json` — members, internal similarity from the clades step,
  sparsification, window, partition size distribution, alignment rate —
  and `commands.log`)
- `partition_summary.tsv` — per-clade summary (456 rows: n, strategy,
  pairs aligned, alignment rate, member-hit rate, partitions, interval
  median/max, >1000bp count, runtime)
- `driver_attempt5.log`, `clade_0_0034_retry.log` — driver logs
- `0/tight_clades.json` — input clade definitions (sha-identical to repo
  `ntm/v3/clades/tight_clades.json.gz`)

Repo: `ntm/v3/partition_report.md` (this file),
`ntm/v3/partition_validation.txt` (executed validator output, exit 0),
`ntm/v3/scripts/validate_v3_partitions.py` (validator),
`scripts/per_clade_alignment_pipeline.py` (stale-index fix).

## Comparison with v2

| | v2 (8,502 prophages) | v3 (9,446 prophages) |
|---|---|---|
| tight clades | 2,388 | 813 |
| alignable clades | 1,251 | 456 |
| singletons | 1,137 | 357 |
| total partitions | 86,587 | 151,237 |
| intervals | 425,090 | 500,584 |
| median partition length | 500 bp | 500 bp |
| intervals > 1000 bp | 544 (0.13%) | 676 (0.135%) |
| failures | 0 | 0 |

v3 has ~3x fewer but much larger clades than v2 (65 clades of 51-100
members vs v2's 6 of 100), so it yields ~1.75x the partitions from ~1.1x
the prophages while the per-partition median stays 500 bp. The clades
step's thr 0.25 / max 100 / community 0 (v3) concentrates related
prophages more aggressively than v2's settings, and 38 v3 clades have
median internal distance 0 (near-identical members).
