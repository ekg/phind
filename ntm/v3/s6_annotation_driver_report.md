# S6 — NTM v3 Pharokka annotation driver

## Deliverables (repo)
- `ntm/v3/scripts/run_annotation_v3.py` — restart-safe driver (prepare / pharokka /
  checkv / split / report), adapted from `ntm/v2/scripts/run_annotation_v2.py`.
- `ntm/v3/scripts/build_annotation_report_v3.py` — faithful port of the v2
  functional-QC report builder to `ntm3_*` IDs (approved option (a)); adds a
  deterministic markdown summary.
- `ntm/v3/scripts/test_annotation_v3.py` — focused tests (prepare round-trip,
  duplicate/count aborts, v1-tree guard, GFF3 split, report cohort labels).

## Inputs / pinned provenance
- Combined input: `$NVME/ntm/v3/ml/all_ntm_v3_ml_phage_genomes.fa` (1304) +
  `all_ntm_v3_ancestral_phage_genomes.fa` (893).
- Default root `$NVME/ntm/v3/annotation` (`$NVME=/mnt/nvme3n1/erikg/phind-genome-work`).
- pharokka v1.10.1 (`/home/erikg/micromamba/envs/pharokka/bin/pharokka`),
  checkv/CheckV v1.1.1 (`/home/erikg/micromamba/envs/phage-annot/bin/checkv`).
- Shared read-only DBs: `$NVME/annotation/pharokka_db`,
  `$NVME/annotation/checkv_db/checkv-db-v1.5` (v1 tree is read-only; driver
  refuses any `--root` inside it and records the DB paths in `run_state/`).

## Full invocation for the parent (long compute — run detached with the process tool)
```
python3 ntm/v3/scripts/run_annotation_v3.py            # all stages, threads=64, root=$NVME/ntm/v3/annotation
# staged / resume variants:
python3 ntm/v3/scripts/run_annotation_v3.py --stages pharokka
python3 ntm/v3/scripts/run_annotation_v3.py --stages split
python3 ntm/v3/scripts/run_annotation_v3.py --force --stages checkv
```
Restart-safe: each stage writes `run_state/<stage>.done.json` (command, tool
versions, DB identity, threads, exit code) and is skipped when the marker plus
outputs exist. Pharokka's `-f` is added automatically if `pharokka_out/` exists
without a done-marker (stale partial).

## Stages
- **prepare** → `input/all_v3_phage_genomes.fa` (bare IDs) + `input/genome_index.tsv`
  (genome_id, source, cohort, clade_id, length, status, host_clades, origin_fa);
  aborts on duplicates / wrong counts / round-trip mismatch.
- **pharokka** → `pharokka run -m --mmseqs2_only --skip_extra_annotations --skip_mash -g prodigal-gv -t 64 --locustag NTMV3`.
- **checkv** → `checkv end_to_end -t 64`.
- **split (NEW)** → `pharokka_out/per_genome_gff/<genome_id>.gff`, one GFF3 per
  genome (`##gff-version 3` + per-genome `##sequence-region`, embedded `##FASTA`
  omitted), plus the committed index `ntm/v3/annotation_gff_index.tsv`
  (columns `genome_id`, `gff_path`, `gene_count`) and an on-NVMe
  `per_genome_gff/manifest.tsv` mirror. Aborts if a feature seqid is absent from
  the prepared input.
- **report** → `build_annotation_report_v3.py --root <root>`; writes the four v2
  TSVs under `ROOT/report/` and `ntm/v3/annotation_report.md`.

## Observed evidence
### prepare (real inputs) — `--dry-run`
```
prepare: wrote .../ntm/v3/annotation/input/all_v3_phage_genomes.fa
         (2197 records, 1304 ML + 893 ancestral; round-trip exact, duplicate-free)
index rows: 2197  cohorts: 893 ntm3_anc | 893 ntm3_ml_reconstructed | 411 ntm3_ml_singleton
input FASTA sha256 = b67cbc1e1ea26ccb…
pharokka=pharokka v1.10.1  checkv=CheckV, version 1.1.1  threads=64 (cap 64)
```
Expected `1304 + 893 = 2197` confirmed.

### split test — 3 genomes in a temp dir
Synthetic pytest case (3 genomes, one CDS-free): 3 GFFs written, regions
preserved, `##FASTA` omitted, index == `genome_id/gff_path/gene_count`
(`…_0000_ML→2`, `…_ANCESTRAL→1`, `…_0001_ML→0`), manifest mirror byte-identical.
Manual run on 3 **real** v3 ML records (`ntm3_0_0000_ML`, `_0001_ML`, `_0002_ML`):
3 GFFs, 6 CDS features, 0 zero-feature, committed index emitted with absolute
GFF paths.

### report builder — synthetic 3-genome root
Produces the four TSVs plus `annotation_report.md`; cohort labels
`ntm3_ml_reconstructed` / `ntm3_ml_singleton` / `ntm3_anc`; tier ladder unchanged
(A for full module set, D for no-detect, `no_annotation` for 0 CDS).

### tests
```
python3 -m pytest ntm/v3/scripts/test_annotation_v3.py -q   → 8 passed
python3 -m py_compile run_annotation_v3.py build_annotation_report_v3.py   → OK
```

## Residual risks
- `split` groups GFF lines by seqid (col 1) and drops everything after the first
  `##FASTA`; pharokka output is CDS-only today. If a future pharokka version
  emits parent/child or non-CDS features the lines are still copied verbatim per
  genome, but gene_count counts all feature lines.
- The committed index stores absolute NVMe GFF paths; running the driver from a
  git worktree writes the index into that worktree's `ntm/v3/`.
- The `report` stage requires CheckV `quality_summary.tsv`; missing genomes abort
  by default (v2 behaviour).
- No full annotation was run here (long compute; parent runs it detached).

## Acceptance report
```acceptance-report
{
  "criteriaSatisfied": [
    {"id": "criterion-1", "status": "satisfied", "evidence": "Added only the three S6 files (driver, report builder, test); stages/flags/roots match the task; no v2/v1 code or outputs modified."},
    {"id": "criterion-2", "status": "satisfied", "evidence": "Real prepare dry-run (1304+893=2197), 3-genome split test (synthetic pytest + real-record manual run), report-builder integration test, and 8 passing pytest cases are reported with commands."}
  ],
  "changedFiles": [
    "ntm/v3/scripts/run_annotation_v3.py",
    "ntm/v3/scripts/build_annotation_report_v3.py",
    "ntm/v3/scripts/test_annotation_v3.py"
  ],
  "testsAddedOrUpdated": [
    "ntm/v3/scripts/test_annotation_v3.py"
  ],
  "commandsRun": [
    {"command": "python3 -m py_compile ntm/v3/scripts/run_annotation_v3.py ntm/v3/scripts/build_annotation_report_v3.py", "result": "passed", "summary": "both modules compile"},
    {"command": "python3 -m pytest ntm/v3/scripts/test_annotation_v3.py -q", "result": "passed", "summary": "8 passed"},
    {"command": "python3 ntm/v3/scripts/run_annotation_v3.py --dry-run", "result": "passed", "summary": "prepare wrote 2197 records (1304 ML + 893 ancestral), round-trip exact"},
    {"command": "python3 - << manual 3-real-genome split_gff >>", "result": "passed", "summary": "3 GFF3 written, committed index genome_id/gff_path/gene_count"},
    {"command": "python3 ntm/v3/scripts/build_annotation_report_v3.py --root <synthetic> --markdown <tmp>", "result": "passed", "summary": "4 TSVs + annotation_report.md produced"}
  ],
  "validationOutput": [
    "prepare: 2197 records, 1304 ML + 893 ancestral; round-trip exact, duplicate-free; cohorts 893 anc / 893 ml_reconstructed / 411 ml_singleton",
    "split: ntm3_0_0000_ML=2, ntm3_0_0000_ANCESTRAL=1, ntm3_0_0001_ML=0; ##FASTA omitted; manifest == committed index",
    "pytest: 8 passed"
  ],
  "residualRisks": [
    "split assumes CDS feature lines grouped by seqid (current pharokka output); future feature types copied verbatim",
    "committed index stores absolute NVMe paths; run from repo checkout to land the index in ntm/v3/",
    "report stage needs CheckV outputs; missing genomes abort by default"
  ],
  "noStagedFiles": true,
  "diffSummary": "New ntm/v3 annotation driver + v3 report builder + tests; no existing files modified.",
  "reviewFindings": ["no blockers"],
  "manualNotes": "Full annotation intentionally not run (long compute); parent should launch detached via the process tool. --dry-run already prepared $NVME/ntm/v3/annotation/input for the real run."
}
```
