#!/usr/bin/env bash
# Post-metadata pipeline for local-e-coli (RUN_PLAN §4-§10).
# Usage: bash run_pipeline.sh   (assumes metadata sweep complete)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NVME=/mnt/nvme3n1/erikg/phind-genome-work/ecoli_screen
PANEL="$REPO_ROOT/research/ecoli_bait/panel/baits.fa"
B2S="$NVME/env/bin/back_to_sequences"
MM2="$NVME/env/bin/minimap2"
RUN_START_FILE="$NVME/RUN_START.txt"
cd "$HERE"

if [ ! -f "$RUN_START_FILE" ]; then
  date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_START_FILE"
fi
RUN_START="$(cat "$RUN_START_FILE")"
DEADLINE="$(date -u -d "$RUN_START + 48 hours" +%Y-%m-%dT%H:%M:%SZ)"
echo "RUN_START=$RUN_START deadline=$DEADLINE"

step() { echo; echo "=== [$(date -u +%H:%M:%S)] $* ==="; }

step "1-2. freeze universe + tier manifests (SKIP if frozen manifests exist)"
if [ ! -f ../manifests/tierE1_manifest.json ] || [ ! -f ../manifests/tierE2_manifest.json ]; then
  python3 scripts/freeze_wgs_universe.py \
    --universe-manifest "$NVME/metadata/ecoli_universe_manifest.json" \
    --runinfo "$NVME/metadata/ecoli_universe_runinfo.tsv" \
    --out "$NVME/metadata/ecoli_wgs_frozen.tsv"
  python3 scripts/build_bioproject_table.py \
    --frozen "$NVME/metadata/ecoli_wgs_frozen.tsv" \
    --out-dir "$NVME/metadata" --manifest-dir ../manifests --log logs/eutils.log
else
  echo "frozen tier manifests present — tier design is immutable (RUN_PLAN §5/§6)"
fi

step "3. disk preflight per tier (seed 20260822)"
for T in E1 E2; do
  python3 scripts/preflight_disk.py --frozen ../manifests/tier${T}_frozen.tsv \
    --tier "$T" --sample 500 --seed 20260822 --out-dir ../results || true
done
df -h /mnt/nvme3n1 | tee -a logs/preflight.log

step "4. screen Tier E1 then Tier E2 (budgets: 400GB / 30k acc / deadline)"
for T in E1 E2; do
  python3 scripts/screen_accessions.py \
    --runs-tsv ../manifests/tier${T}_frozen.tsv \
    --baits "$PANEL" --run-root "$NVME" --tier "$T" \
    --b2s "$B2S" --concurrency 4 --max-bytes $((400*1024**3)) \
    --max-accessions 30000 --deadline "$DEADLINE" 2>&1 | tee logs/screen_${T}.log
done

step "5. negative-gate check over ALL screened accessions + summarize"
python3 scripts/summarize_screen_ecoli.py \
  --results "$NVME/screen/E1_results.jsonl" "$NVME/screen/E2_results.jsonl" \
  --tiers E1 E2 \
  --functional-qc-tsv "$REPO_ROOT/research/phage_annotation/per_genome_annotation_qc.tsv" \
  --bait-manifest "$REPO_ROOT/research/ecoli_bait/panel/bait_manifest.tsv" \
  --availability "$NVME/downloads/availability.tsv" \
  --runs-tsv ../manifests/tierE1_frozen.tsv ../manifests/tierE2_frozen.tsv \
  --out-dir ../results

step "6. confirmation of hits (minimap2 ladder, frozen classification)"
: > "$NVME/confirm/confirm_E1E2.jsonl"
for T in E1 E2; do
  python3 scripts/confirm_tier_hits_ecoli.py \
    --results "$NVME/screen/${T}_results.jsonl" \
    --baits "$PANEL" --downloads "$NVME/downloads" \
    --b2s "$B2S" --minimap2 "$MM2" \
    --availability "$NVME/downloads/availability.tsv" \
    --runs-tsv ../manifests/tier${T}_frozen.tsv \
    --out "$NVME/confirm/confirm_E1E2.jsonl" 2>&1 | tee logs/confirm_${T}.log
done

step "7. full-provenance TSV"
python3 scripts/emit_full_provenance.py \
  --confirm-jsonl "$NVME/confirm/confirm_E1E2.jsonl" \
  --availability "$NVME/downloads/availability.tsv" \
  --out ../results/confirmed_hits_full_provenance.tsv

step "8. known-phage calibration (lambda / HK97 / P2)"
python3 scripts/known_phage_calibration.py \
  --confirm-jsonl "$NVME/confirm/confirm_E1E2.jsonl" \
  --baits "$PANEL" \
  --lambda-fa "$REPO_ROOT/research/ecoli_bait/panel/public_refs/NC_001416.1.fa" \
  --hk97-fa "$REPO_ROOT/research/ecoli_bait/panel/public_refs/NC_002167.1.fa" \
  --minimap2 "$MM2" --out-dir ../results --log logs/eutils.log

step "9. copy small tables to repo + NVMe manifest receipt"
cp "$NVME/results/controls_gate.json" ../results/ 2>/dev/null || true
cp "$NVME/results/"*_preflight.json ../results/ 2>/dev/null || true
python3 scripts/emit_nvme_manifest.py --root "$NVME" --out ../NVMe_MANIFEST.tsv
gzip -kf ../NVMe_MANIFEST.tsv
echo "pipeline complete"
