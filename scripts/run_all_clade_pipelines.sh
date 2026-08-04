#!/usr/bin/env bash
# run_all_clade_pipelines.sh — run the per-clade alignment pipeline for all
# 12 communities concurrently (one process per community, jobs×threads ≈ 216).
# Log: research/clades/pipeline_driver.log (plus per-community pipeline logs).
set -u
cd "$(dirname "$0")/.."
ROOT=$(pwd)
PY=python3
PIPE="$ROOT/scripts/per_clade_alignment_pipeline.py"
FASTA="$ROOT/prophage_homology_survey/full_prophages.fa"
IDX="$ROOT/research/clades/full_prophages.idx.json"
mkdir -p "$ROOT/research/clades"
[ -f "$IDX" ] || "$PY" "$PIPE" --build-index --fasta "$FASTA" --index "$IDX"

declare -A JOBS=( [0]=6 [1]=2 [2]=2 [3]=2 [4]=3 [5]=1 [6]=1 [7]=1 [8]=1 [9]=1 [10]=6 [11]=1 )
THREADS=8
PIDS=()
for c in 0 1 2 3 4 5 6 7 8 9 10 11; do
  OUT="$ROOT/research/clades/$c"
  CLADES="$OUT/tight_clades.json"
  [ -f "$CLADES" ] || continue
  echo "[$(date -u +%FT%TZ)] community $c: jobs=${JOBS[$c]} threads=$THREADS" \
    >> "$ROOT/research/clades/pipeline_driver.log"
  "$PY" "$PIPE" --community "$c" --outdir "$OUT" --clades "$CLADES" \
      --fasta "$FASTA" --index "$IDX" \
      --threads "$THREADS" --jobs "${JOBS[$c]}" --all \
      > "$OUT/pipeline_run.log" 2>&1 &
  PIDS+=($!)
done
rc=0
for p in "${PIDS[@]}"; do
  wait "$p" || rc=1
done
echo "[$(date -u +%FT%TZ)] all communities done rc=$rc" >> "$ROOT/research/clades/pipeline_driver.log"
exit $rc
