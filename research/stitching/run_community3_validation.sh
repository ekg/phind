#!/usr/bin/env bash
# Reproducible community-3 stitching validation.
#
# Runs the partition-stitching pipeline end-to-end from BED + MAF inputs and
# regenerates the stitched consensus FASTA(s) + results JSON in
# research/stitching/.
#
# Inputs (committed in this repo):
#   inputs/community_3_partitions/  502 partition MAF alignment blocks
#   inputs/3_ancestral.fa           pggb cluster-3 ancestral (113,502 bp; see
#                                   validation_report.md for provenance caveats)
#   community_3_partitions.bed      partition assignments per prophage
#
# Outputs:
#   community_3_stitched_mean.fa        50% threshold (strict core)
#   community_3_stitched_mean_45pct.fa  45% threshold (extended core)
#   stitching_results.json              machine-readable results (50%)
#   stitching_results_45pct.json        machine-readable results (45%)
#
# Committed outputs are byte-identical to the original build at 2363ece
# (verified sha256), so this script is the reproducible entry point.

set -euo pipefail
cd "$(dirname "$0")"

PY=python3
PART_DIR=inputs/community_3_partitions
BED=community_3_partitions.bed
ANC=inputs/3_ancestral.fa

echo "== [1/2] 50% accessory threshold (strict core genome) =="
$PY stitch_algorithm.py \
    --partition-dir "$PART_DIR" --bed "$BED" \
    --output community_3_stitched_mean.fa \
    --ancestral "$ANC" --accessory-threshold 0.5 \
    --coverage-threshold 0.0 --json stitching_results.json

echo
echo "== [2/2] 45% accessory threshold (extended core) =="
$PY stitch_algorithm.py \
    --partition-dir "$PART_DIR" --bed "$BED" \
    --output community_3_stitched_mean_45pct.fa \
    --ancestral "$ANC" --accessory-threshold 0.45 \
    --coverage-threshold 0.0 --json stitching_results_45pct.json

echo
echo "== Outputs =="
for f in community_3_stitched_mean.fa community_3_stitched_mean_45pct.fa \
         stitching_results.json stitching_results_45pct.json; do
    echo "  $f : $(sha256sum "$f" | cut -d' ' -f1)"
done
