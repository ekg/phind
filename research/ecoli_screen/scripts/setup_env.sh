#!/usr/bin/env bash
# Deterministic environment for the local NTM screen (RUN_PLAN §9).
# Recreates the micromamba env used by ntm/v2/pilot/confirm_hits.py
# (minimap2 2.31 + back_to_sequences 0.8.4, previously at /tmp/mmenv).
set -euo pipefail

NVME_ROOT="${NVME_ROOT:-/mnt/nvme3n1/erikg/phind-genome-work/ecoli_screen}"
MICROMAMBA="${MICROMAMBA:-$HOME/.local/bin/micromamba}"

"$MICROMAMBA" create -y -p "$NVME_ROOT/env" \
  -c conda-forge -c bioconda \
  'back_to_sequences=0.8.4' 'minimap2=2.31'

{
  echo "installed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "micromamba=$("$MICROMAMBA" --version)"
  echo "minimap2=$("$NVME_ROOT/env/bin/minimap2" --version)"
  echo "back_to_sequences=$("$NVME_ROOT/env/bin/back_to_sequences" --version 2>&1 | head -1)"
  echo "zstd=$(zstd --version 2>&1 | head -1)"
  echo "python3=$(python3 --version)"
} > "$NVME_ROOT/env/VERSIONS.txt"

cat "$NVME_ROOT/env/VERSIONS.txt"
