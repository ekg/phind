#!/usr/bin/env bash
# NTM v2 run-assembly watcher: fetch collaborator uploads from hypervolu.me,
# ingest into PanSN bgzip layout, verify against the manifest.
#
# Poll the destination; when files appear, rsync them to a local staging
# dir, run the ingest+verify script, and report. Designed to be re-run
# idempotently (rsync --partial resumable; ingest skips already-verified? no —
# ingest re-processes everything present; that is fine at cohort scale but if
# it becomes slow, add a --done marker per accession).
#
# Usage:
#   bash ntm/v2/scripts/fetch_ingest_run_assemblies.sh [--dry-run]
set -euo pipefail

DEST="erik@hypervolu.me:~/www/phage/ntm_v2_run_assemblies/"
STAGING=/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/incoming/run_assemblies
OUT=/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/genomes/run_assemblies
MANIFEST=/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/inputs/NTM_QC_passed_prophage_master_manifest.tsv
ACC_LIST=/home/erikg/phind/ntm/v2/run_assemblies/run_assemblies_needed.txt
REPORT="$OUT/ingest_report.md"
SUMMARY="$OUT/ingest_summary.json"

DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

echo "[watch] checking $DEST"
if ! ssh -o BatchMode=yes -o ConnectTimeout=15 erik@hypervolu.me \
      "test -d ~/www/phage/ntm_v2_run_assemblies && ls ~/www/phage/ntm_v2_run_assemblies/ | head -20"; then
  echo "[watch] destination empty or unreachable — no collaborator upload yet"
  exit 0
fi

mkdir -p "$STAGING" "$OUT"
if [ "$DRY" = 1 ]; then
  echo "[watch] DRY-RUN: would rsync from $DEST"
  exit 0
fi

echo "[watch] fetching from hypervolu.me..."
rsync -avz --partial \
  "$DEST" "$STAGING/"

echo "[watch] ingesting + verifying..."
python3 /home/erikg/phind/ntm/v2/scripts/ingest_run_assemblies.py \
  --manifest "$MANIFEST" \
  --staging "$STAGING" \
  --out "$OUT" \
  --accessions "$ACC_LIST" \
  --checksums "$STAGING/sha256sums.txt" \
  --report "$REPORT" \
  --summary "$SUMMARY"
echo "[watch] done — see $REPORT"
