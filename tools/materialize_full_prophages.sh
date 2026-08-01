#!/usr/bin/env bash
set -euo pipefail

repo=$(git rev-parse --show-toplevel)
compressed="$repo/prophage_homology_survey/full_prophages.fa.gz"
output="$repo/prophage_homology_survey/full_prophages.fa"
expected="ed85b2fb549be18bc638d8485f5b5add7c2d394f3822efe66a90ca6d979758d3"

if [[ ! -f "$compressed" ]]; then
  echo "missing LFS artifact: $compressed" >&2
  echo "run: git lfs pull" >&2
  exit 1
fi

if [[ -e "$output" ]]; then
  actual=$(sha256sum "$output" | awk '{print $1}')
  if [[ "$actual" == "$expected" ]]; then
    echo "already materialized and verified: $output"
    exit 0
  fi
  echo "refusing to replace existing file with unexpected digest: $output" >&2
  exit 1
fi

tmp="$output.tmp.$$"
trap 'rm -f "$tmp"' EXIT
bgzip -cd "$compressed" > "$tmp"
printf '%s  %s\n' "$expected" "$tmp" | sha256sum -c -
mv "$tmp" "$output"
trap - EXIT
echo "materialized and verified: $output"
