# NTM v3 — per-clade allwave + impg partition report (full cohort)

Generated: 2026-10-06 (full-cohort regeneration on the 2026-09-22 collaborator
delivery). Supersedes the 2026-09-14 public-portion report, archived at
`archive_partition_report.public-partial.md`.

## Run

Per-clade `scripts/per_clade_alignment_pipeline.py` driven to completion by
`ntm/scripts/resume_per_clade_partition.py` (resume-first, `--max-attempts 5`):

```
python3 ntm/scripts/resume_per_clade_partition.py \
  --outdir  $NVME/ntm/v3/clades \
  --clades  $NVME/ntm/v3/clades/0/tight_clades.json \
  --fasta   $NVME/ntm/v3/full_prophages.fa \
  --index   $NVME/ntm/v3/mash_clades/full_prophages.idx.json \
  --threads 8 --jobs 24 --max-attempts 5
```

Allwave sparsified k-nearest (`tree:k:0:0`, k-farthest = 0, no stranger-joining),
scores `0,5,8,2,24,1`; impg partition window 500. Singletons pass through
without alignment. `$NVME = /mnt/nvme3n1/erikg/phind-genome-work`.

**Note (stale-resume hazard, avoided):** the 813 public-portion per-clade dirs
were moved to `$NVME/ntm/v3/clades_public_partial/` before this run, because the
re-derived clade IDs had different membership and the resume logic treats any
`manifest.json` as complete. A `attempt 1: 1304 incomplete clades (of 1304)`
start is the expected evidence that this was handled.

## Result

| metric | value |
|---|---:|
| clades total | 1,304 |
| alignable (n>=2, full PAF+BED+MAF+manifest) | 893 |
| singletons (pass-through) | 411 |
| missing / failed / manifest issues | 0 / 0 / 0 |
| strategies | `tree:5:0:0.0` 374, `none` 519 |
| alignment rate (pairs/possible) | median 1.0000, min 0.0661 |
| member-hit rate (seqs in PAF / n) | median 1.0000, min 1.0000 |
| total partitions | 726,121 |
| total intervals | 1,847,004 |
| per-clade median interval | 500 bp (max 926) |
| intervals > 1000 bp | 1,579 (0.085%) |
| intervals < 100 bp | 169,744 (9.19%) |
| global max interval span | 2,180 bp |
| clades with >1000 bp intervals | 305 of 893 |
| per-clade runtime | total 105,448 s, median 8.6 s, max 4,198.6 s |

Validator output: `partition_validation.txt` (`VALIDATION: PASS`).

## Outputs

`$NVME/ntm/v3/clades/`:
- `0_0000` .. `0_1303`: per-clade dirs (`sequences.fa`, `allwave.paf`,
  `allwave.segmented.paf`, `partitions.bed`, `partitions/partition<N>.maf`,
  `manifest.json`, `commands.log`)
- `partition_summary.tsv` — 893 rows
- `0/tight_clades.json` — input clade definitions (1304 clades)

Repo: this report, `partition_validation.txt`, `ntm/v3/scripts/validate_v3_partitions.py`,
`scripts/per_clade_alignment_pipeline.py`.

## Public-portion → full-cohort comparison

| | public portion (2026-09-14) | **full cohort (2026-10-06)** |
|---|---:|---:|
| prophages | 9,446 | 36,857 |
| tight clades | 813 | 1,304 |
| alignable | 456 | 893 |
| singletons | 357 | 411 |
| total partitions | 151,237 | 726,121 |
| intervals | 500,584 | 1,847,004 |
| median partition length | 500 bp | 500 bp |
| intervals > 1000 bp | 676 (0.135%) | 1,579 (0.085%) |
| failures | 0 | 0 |

## Scope

The cohort is the **union** of the export and v2 holdings, including the MTC
block. Scope is a **label**, not a filter (see `PLAYBOOK.md` scope decision):
391 prophage members derive from MTC host genomes and participate in normal
clade alignment. No population was excluded at this stage.
