#!/usr/bin/env bash
# Stage 0a — extract collaborator delivery FASTAs to NVMe (raw, verifiable).
# Source: /home/erikg/phind/NTM_Data (42 GB zips) -> NVMe raw tree.
# Expected: 17,920 sra_assembled + 8,579 ncbi_bvbrc = 26,499 files.
# No WG involvement. Re-runnable (unzip -o).
set -euo pipefail
DATA=/home/erikg/phind/NTM_Data
OUT=/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/genomes/delivery_raw
mkdir -p "$OUT/sra_assembled" "$OUT/ncbi_bvbrc"

echo "INGEST_START $(date -u +%FT%TZ)"
unzip -q -o "$DATA/Genomes_sra_assembled_17920_20260922.zip" -d "$OUT/sra_assembled"
echo "SRA_UNZIP_DONE $(date -u +%FT%TZ)"
unzip -q -o "$DATA/Genomes_ncbi_bvbrc_8579_20260922.zip" -d "$OUT/ncbi_bvbrc"
echo "NCBI_UNZIP_DONE $(date -u +%FT%TZ)"

n1=$(find "$OUT/sra_assembled" -maxdepth 1 -name '*.fasta' | wc -l)
n2=$(find "$OUT/ncbi_bvbrc" -maxdepth 1 -name '*.fasta' | wc -l)
echo "INGEST_DONE sra=$n1 ncbi_bvbrc=$n2 total=$((n1+n2)) expected=17920/8579/26499"
