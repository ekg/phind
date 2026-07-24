#!/usr/bin/env bash
# Recreate the local executable/help and host-environment evidence used by
# reports/phylogeny_design.md.  This script does NOT run a sequence benchmark,
# download data, convert BGZF, or execute any distance/tree workload.
set -uo pipefail

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
mkdir -p "$HERE"

{
  printf 'captured_at=%s\n' "$(date --iso-8601=seconds)"
  printf 'working_directory=%s\n' "$(pwd -P)"
  printf 'shell=%s\n' "${SHELL:-unknown}"
  printf '\n[uname]\n'
  uname -a
  printf '\n[lscpu selected fields]\n'
  lscpu 2>&1 | grep -E '^(Architecture|CPU\(s\)|On-line CPU|Model name|Thread|Core|Socket|NUMA)' || true
  printf '\n[free -h]\n'
  free -h 2>&1 || true
  printf '\n[df -h current filesystem]\n'
  df -h . 2>&1 || true
  printf '\n[locale]\n'
  locale 2>&1 || true
  printf '\n[python]\n'
  command -v python3 2>&1 || true
  python3 --version 2>&1 || true
} > "$HERE/environment.txt"

# These are the relevant executable names, including the user's literal "mesh"
# spellings.  Help is bounded to 120 lines per command.
tools=(
  mesh mesh-distance mesh_triangle mesh-triangle
  mash mashtree MashTree rapidnj quicktree FastTree fasttree FastTreeMP
  iqtree iqtree2 iqtree3 snippy snippy-core run_gubbins.py gubbins
  fastANI skani checkm checkm2 quast.py quast
)

{
  printf 'captured_at=%s\n' "$(date --iso-8601=seconds)"
  printf 'PATH=%s\n' "$PATH"
  printf 'purpose=Installed executable/version/help discovery; no tool was installed or run on sequences.\n'
  for tool in "${tools[@]}"; do
    printf '\n[%s]\n' "$tool"
    if path=$(command -v "$tool" 2>/dev/null); then
      printf 'status=FOUND_ON_PATH\npath=%s\n' "$path"
      if [ -f "$path" ]; then
        sha256sum "$path" 2>&1 | sed 's/^/executable_sha256=/' || true
      fi
      printf '%s\n' '-- version attempt --'
      "$tool" --version 2>&1 | sed -n '1,20p' || true
      printf '%s\n' '-- help attempt --'
      "$tool" --help 2>&1 | sed -n '1,120p' || true
    else
      printf 'status=NOT_FOUND_ON_PATH\n'
    fi
  done

  printf '\n[dpkg-query]\n'
  for package in mash mashtree rapidnj iqtree snippy gubbins fastani skani checkm2 quast; do
    dpkg-query -W -f='${Package}\t${Status}\t${Version}\n' "$package" 2>&1 || true
  done

  printf '\n[micromamba environments and matching package metadata]\n'
  if command -v micromamba >/dev/null 2>&1; then
    micromamba env list 2>&1 | sed 's/[[:space:]]*$//' || true
    while read -r prefix; do
      [ -d "$prefix/conda-meta" ] || continue
      find "$prefix/conda-meta" -maxdepth 1 -type f \
        \( -name 'mash-*.json' -o -name 'mashtree-*.json' \
        -o -name 'rapidnj-*.json' -o -name 'iqtree-*.json' \
        -o -name 'snippy-*.json' -o -name 'gubbins-*.json' \
        -o -name 'fastani-*.json' -o -name 'skani-*.json' \
        -o -name 'checkm2-*.json' -o -name 'quast-*.json' \) \
        -print 2>/dev/null || true
    done < <(micromamba env list 2>/dev/null | awk '$NF ~ /^\// {print $NF}')
  else
    printf 'micromamba=NOT_FOUND_ON_PATH\n'
  fi
} > "$HERE/tool_versions.txt"

printf 'Wrote %s and %s\n' "$HERE/environment.txt" "$HERE/tool_versions.txt"
printf 'No benchmark was run (input_count=0; input_bytes=0; output_bytes=0).\n'
