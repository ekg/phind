# S0b — PanSN canonical-object build for delivered SRA runs

## Deliverable
Driver: `ntm/v3/scripts/build_v3_pansn_objects.py`

Converts delivered SRA-run FASTAs (bare-name headers) into canonical PanSN bgzip
genome objects. Reuses `contig_name_from_header`, `rename_fasta_to_pansn`, and
`write_bgzip_faidx` imported directly from `ntm/v3/scripts/acquire_v3_genomes.py`
(no re-implementation of naming/bgzip logic).

Layout produced (one object per run, disjoint dirs):
```
$NVME/ntm/v3/genomes/canonical_objects/<run>/<run>.pansn.fa.gz
                                          <run>.pansn.fa.gz.fai
                                          <run>.pansn.fa.gz.gzi
```
`<run>.pansn.fa.gz` is written via temp file + `os.replace` (atomic); temp
`.fa.tmp`/`.pansn.fa.gz.tmp` are removed in a `finally`. Pre-existing object
dirs are never deleted or pruned.

## Scope check (read-only, no build run)
- Master table `NTM_Data/NTM_master_table_complete_20260922.tsv`: `source == ASSEMBLY` → **17,920** unique `genome_id`, 17,920 rows (26,500 total rows). Matches the manifest's `plan == blocked_run` / `source_ns == ASSEMBLY` count (18,055 manifest rows = 17,920 runs after collapsing multi-row keys).
- `--dry-run` against defaults: `runs to consider: 17920 ... reused 0`, `missing input: 12582` (extraction into `delivery_raw/sra_assembled` is still in progress; count grows over time). Resume-first means the later full run will build whatever is present and reuse anything already built.

## Exact full-build invocation (for the parent, run detached via the `process` tool)
```bash
python3 ntm/v3/scripts/build_v3_pansn_objects.py --jobs 16
```
Defaults: `--raw-dir $NVME/ntm/v3/genomes/delivery_raw/sra_assembled`,
`--out-dir $NVME/ntm/v3/genomes/canonical_objects`,
`--manifest NTM_Data/NTM_master_table_complete_20260922.tsv`.
Summary written to `$NVME/ntm/v3/genomes/pansn_build_summary.{json,tsv}`.

It is resumable: re-running skips every run whose `<run>.pansn.fa.gz.fai`
exists and is non-empty. Final stdout line: `BUILD_DONE built=<n> reused=<n> failed=<n>`
(exit 1 iff `failed>0`). `--limit N`, `--run <id>` (repeatable), `--no-resume`,
`--dry-run` are supported. Parallelism is `ProcessPoolExecutor(jobs)` over
independent run dirs — no shared mutable state.

## Test performed (3 runs only; full build NOT run)
Test root `/tmp/s0b_test/{raw,out}`. Members streamed from
`NTM_Data/Genomes_sra_assembled_17920_20260922.zip` via python `zipfile`
(no dependency on the concurrent full extraction).

Commands + observed:
```
# 1. stream 3 members
python3 -c "zipfile ... read DRR015955/56/57.fasta -> /tmp/s0b_test/raw"
# 2. build
python3 ntm/v3/scripts/build_v3_pansn_objects.py --raw-dir /tmp/s0b_test/raw \
    --out-dir /tmp/s0b_test/out --jobs 4 --run DRR015955 --run DRR015956 --run DRR015957
=> BUILD_DONE built=3 reused=0 failed=0
# 3. objects + indexes
ls out/DRR01595{5,6,7}/ => each has .pansn.fa.gz + .pansn.fa.gz.fai + .pansn.fa.gz.gzi
# 4. coordinate/sequence resolution
samtools faidx out/DRR015955/DRR015955.pansn.fa.gz \
    'DRR015955#1#NODE_14_length_150130_cov_107.750611'
=> returns sequence, length = 150130  (tr -d '\n' | wc -c)
# 5. bare names == coordinates-CSV scaffolds for DRR015955
fai bare names (strip 'run#1#') superset of
['NODE_14_length_150130_cov_107.750611','NODE_69_length_20432_cov_121.040192'] -> True
# 6. resume
same build command again => BUILD_DONE built=0 reused=3 failed=0
# 7. integrity
bgzip -t DRR015955.pansn.fa.gz => ok ; .fai col1 is PanSN 'DRR015955#1#...'
```

### Per-object results (from `/tmp/s0b_test/pansn_build_summary.tsv`)
| run | contigs | total_bp | object_bytes |
|---|---|---|---|
| DRR015955 | 229 | 7,563,317 | 2,109,371 |
| DRR015956 | 401 | 6,962,622 | 1,937,982 |
| DRR015957 | 632 | 6,985,838 | 1,958,052 |

Contig counts equal the master-table `n_contigs` for these runs (229/401/632),
and `.fai` first column is the PanSN name (`<run>#1#<contig>`).

Resume proof: second invocation → `built=0 reused=3 failed=0`.

## Residual risk
- Input extraction is concurrent; a file being written mid-copy could yield a
  truncated FASTA. Mitigation: resume skips only when `.fai` is non-empty; a
  corrupt build would still be recorded as `failed` in the summary. Parent
  should re-run once extraction completes to sweep any failures (resume-safe).
- `--no-resume` intentionally rebuilds existing objects (used only if a bad
  object must be regenerated); it still never deletes the object dir.
- No other risk identified; no shared state, no writes outside the run dirs and
  the two summary files.
