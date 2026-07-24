#!/usr/bin/env bash
# Bounded, synthetic-only BGZF x literal-# compatibility probe.
# This script deliberately does not read or download any production genome.
set -uo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SYN="$ROOT/synthetic"
export LC_ALL=C
umask 022
rm -rf -- "$SYN"
mkdir -p -- "$SYN/input" "$SYN/canonical" "$SYN/logs" "$SYN/results" "$SYN/impg_index"

record_command() {
    local label=$1
    shift
    {
        printf 'argv:'
        printf ' %q' "$@"
        printf '\n'
    } >"$SYN/logs/${label}.command.txt"
    set +e
    "$@" >"$SYN/logs/${label}.stdout.txt" 2>"$SYN/logs/${label}.stderr.txt"
    local rc=$?
    set -e
    printf '%s\n' "$rc" >"$SYN/logs/${label}.exit_status.txt"
    return "$rc"
}

record_shell() {
    local label=$1
    local script=$2
    printf 'bash -o pipefail -c %q\n' "$script" >"$SYN/logs/${label}.command.txt"
    set +e
    bash -o pipefail -c "$script" >"$SYN/logs/${label}.stdout.txt" 2>"$SYN/logs/${label}.stderr.txt"
    local rc=$?
    set -e
    printf '%s\n' "$rc" >"$SYN/logs/${label}.exit_status.txt"
    return "$rc"
}

# Record the execution environment and exact installed/absent tools.  Binary hashes
# disambiguate locally rebuilt programs that report the same semantic version.
{
    printf 'generated_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'probe_root=%s\n' "$ROOT"
    printf 'kernel='; uname -srmo
    printf 'hostname='; hostname
    printf 'architecture='; uname -m
    printf 'locale=%s\n' "${LC_ALL}"
    printf 'shell=%s\n' "${BASH_VERSION}"
    printf 'path=%s\n' "$PATH"
    if [[ -r /etc/os-release ]]; then
        grep -E '^(PRETTY_NAME|VERSION_ID)=' /etc/os-release
    fi
    printf 'scope=synthetic_only; no network; no production inputs read\n'
} >"$ROOT/environment.txt"

{
    printf '# Exact commands, paths, exit statuses, version output, and executable SHA-256\n'
    printf '# Trailing horizontal whitespace in version output is normalized.\n'
    for tool in bgzip samtools wfmash impg gzip python3; do
        printf '\n[%s]\n' "$tool"
        if path=$(command -v "$tool" 2>/dev/null); then
            printf 'path=%s\n' "$path"
            printf 'executable_sha256=%s\n' "$(sha256sum "$path" | awk '{print $1}')"
            case "$tool" in
                gzip) args=(--version) ;;
                python3) args=(--version) ;;
                *) args=(--version) ;;
            esac
            printf 'command=%s' "$tool"
            printf ' %q' "${args[@]}"
            printf '\n'
            set +e
            output=$("$tool" "${args[@]}" 2>&1)
            rc=$?
            set -e
            printf 'exit_status=%s\n' "$rc"
            printf '%s\n' "$output" | sed 's/[[:blank:]]*$//'
        else
            printf 'status=NOT_INSTALLED\n'
        fi
    done
    printf '\n[optional_or_likely_downstream]\n'
    for tool in mash mashmap skani pggb seqwish odgi vg minigraph minigraph-cactus cactus-pangenome Panaroo panaroo roary pangraph panta; do
        if path=$(command -v "$tool" 2>/dev/null); then
            printf '%s\tINSTALLED\t%s\tsha256=%s\n' "$tool" "$path" "$(sha256sum "$path" | awk '{print $1}')"
        else
            printf '%s\tNOT_INSTALLED\n' "$tool"
        fi
    done
} >"$ROOT/tool_versions.txt"

set -e

# Exactly two fictional assemblies are generated.  There are four retained FASTA
# representations: one plain source, one ordinary-gzip source, and two canonical
# BGZF files.  The names and metadata are visibly synthetic.
record_shell generate_inputs "python3 - '$SYN/input/EXAMPLE_ASM_A_v1.source.fa' '$SYN/input/EXAMPLE_ASM_B_v2.source.fa.gz' <<'PY'
import gzip
import io
import random
import sys
from pathlib import Path

plain_path, gzip_path = map(Path, sys.argv[1:])
rng = random.Random(20260724)
chrom_a = ''.join(rng.choice('ACGT') for _ in range(24000))
chrom_b = list(chrom_a)
for i in range(211, len(chrom_b), 887):
    chrom_b[i] = {'A':'C','C':'G','G':'T','T':'A'}[chrom_b[i]]
chrom_b = ''.join(chrom_b)
plasmid = ''.join(rng.choice('ACGT') for _ in range(4800))
unplaced = ''.join(rng.choice('ACGT') for _ in range(3600))

def fasta(records):
    out = io.StringIO()
    for name, seq in records:
        out.write('>' + name + '\\n')
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + '\\n')
    return out.getvalue().encode('ascii')

a = fasta([
    ('EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1', chrom_a),
    ('EXAMPLE_ASM_A_v1#1#EXAMPLE_PLASMID_v1', plasmid),
])
b = fasta([
    ('EXAMPLE_ASM_B_v2#1#EXAMPLE_CHROMOSOME_v3', chrom_b),
    ('EXAMPLE_ASM_B_v2#1#EXAMPLE_UNPLACED_01', unplaced),
])
plain_path.write_bytes(a)
with gzip.GzipFile(filename='', mode='wb', fileobj=gzip_path.open('wb'), mtime=0) as fh:
    fh.write(b)
PY"

# Stream plain FASTA and ordinary gzip into independent .part files; validate before
# same-filesystem promotion.  No uncompressed copy of assembly B is retained.
A_SRC="$SYN/input/EXAMPLE_ASM_A_v1.source.fa"
B_SRC="$SYN/input/EXAMPLE_ASM_B_v2.source.fa.gz"
A_BGZF="$SYN/canonical/EXAMPLE_ASM_A_v1.pansn.fa.gz"
B_BGZF="$SYN/canonical/EXAMPLE_ASM_B_v2.pansn.fa.gz"
record_command gzip_test_source_b gzip -t "$B_SRC"
record_shell convert_plain_to_bgzf "bgzip -@ 2 -l 6 --binary -c '$A_SRC' > '$A_BGZF.part' && bgzip -t '$A_BGZF.part' && mv -- '$A_BGZF.part' '$A_BGZF'"
record_shell convert_gzip_to_bgzf "gzip -cd -- '$B_SRC' | bgzip -@ 2 -l 6 --binary -c > '$B_BGZF.part' && bgzip -t '$B_BGZF.part' && mv -- '$B_BGZF.part' '$B_BGZF'"
record_command bgzip_test_a bgzip -t "$A_BGZF"
record_command bgzip_test_b bgzip -t "$B_BGZF"

# Check decompressed-byte identity and record both content and compressed-byte hashes.
record_shell checksums "{
  printf 'source_a_content_sha256  '; sha256sum '$A_SRC' | awk '{print \$1}'
  printf 'canonical_a_content_sha256  '; bgzip -cd '$A_BGZF' | sha256sum | awk '{print \$1}'
  printf 'source_b_content_sha256  '; gzip -cd '$B_SRC' | sha256sum | awk '{print \$1}'
  printf 'canonical_b_content_sha256  '; bgzip -cd '$B_BGZF' | sha256sum | awk '{print \$1}'
  printf 'canonical_a_bgzf_sha256  '; sha256sum '$A_BGZF' | awk '{print \$1}'
  printf 'canonical_b_bgzf_sha256  '; sha256sum '$B_BGZF' | awk '{print \$1}'
} > '$SYN/results/checksums.txt'
a=\$(awk '/source_a_content/{print \$2}' '$SYN/results/checksums.txt')
aa=\$(awk '/canonical_a_content/{print \$2}' '$SYN/results/checksums.txt')
b=\$(awk '/source_b_content/{print \$2}' '$SYN/results/checksums.txt')
bb=\$(awk '/canonical_b_content/{print \$2}' '$SYN/results/checksums.txt')
[[ \$a == \$aa && \$b == \$bb ]]"

# samtools must build both .fai and .gzi and retrieve a region whose sequence name
# contains two literal '#'.  The argument is quoted and passed as one argv element.
record_command samtools_faidx_a samtools faidx "$A_BGZF"
record_command samtools_faidx_b samtools faidx "$B_BGZF"
record_command samtools_region_hash samtools faidx "$A_BGZF" 'EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1:11-40'
record_shell verify_samtools "test -s '$A_BGZF.fai' && test -s '$A_BGZF.gzi' && test -s '$B_BGZF.fai' && test -s '$B_BGZF.gzi'
printf '%s\\n' \
 'EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1' \
 'EXAMPLE_ASM_A_v1#1#EXAMPLE_PLASMID_v1' > '$SYN/results/expected_a_names.txt'
printf '%s\\n' \
 'EXAMPLE_ASM_B_v2#1#EXAMPLE_CHROMOSOME_v3' \
 'EXAMPLE_ASM_B_v2#1#EXAMPLE_UNPLACED_01' > '$SYN/results/expected_b_names.txt'
cut -f1 '$A_BGZF.fai' > '$SYN/results/observed_a_fai_names.txt'
cut -f1 '$B_BGZF.fai' > '$SYN/results/observed_b_fai_names.txt'
diff -u '$SYN/results/expected_a_names.txt' '$SYN/results/observed_a_fai_names.txt' > '$SYN/results/samtools_a_names.diff'
diff -u '$SYN/results/expected_b_names.txt' '$SYN/results/observed_b_fai_names.txt' > '$SYN/results/samtools_b_names.diff'
grep -Fx '>EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1:11-40' '$SYN/logs/samtools_region_hash.stdout.txt'"

# wfmash is the installed MashMap-family mapper.  Pass both BGZF inputs and require
# at least one PAF row with exact, unmodified PanSN query and target names.
record_command wfmash_bgzf_hash wfmash "$A_BGZF" "$B_BGZF" -m -p 80 -l 500 -w 100 -s 100 -k 15 -t 1 -f --quiet
cp -- "$SYN/logs/wfmash_bgzf_hash.stdout.txt" "$SYN/results/wfmash_bgzf_hash.paf"
record_shell verify_wfmash "awk -F '\\t' '
BEGIN { ok=0 }
NF>=12 {
  if (\$1 !~ /^EXAMPLE_ASM_B_v2#1#/) exit 20
  if (\$6 !~ /^EXAMPLE_ASM_A_v1#1#/) exit 21
  ok=1
}
END { if (!ok) exit 22 }
' '$SYN/results/wfmash_bgzf_hash.paf'"

# IMPG's syng FASTA reader and map output are tested; this is an index/read boundary
# probe, not a pangenome construction.  Exact source names must survive in names and
# projected PAF.  Use one thread and raw small synthetic inputs only.
IMPG_PREFIX="$SYN/impg_index/example_a"
record_command impg_syng_bgzf_hash impg syng -f "$A_BGZF" -o "$IMPG_PREFIX" -t 1 -v 1
record_command impg_map_bgzf_hash impg map -a "$IMPG_PREFIX" -q "$B_BGZF" -o paf -O "$SYN/results/impg_map_bgzf_hash.paf" -t 1 -v 1
record_shell verify_impg "test -s '$IMPG_PREFIX.syng.names'
grep -F 'EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1' '$IMPG_PREFIX.syng.names' >/dev/null
grep -F 'EXAMPLE_ASM_A_v1#1#EXAMPLE_PLASMID_v1' '$IMPG_PREFIX.syng.names' >/dev/null
test -s '$SYN/results/impg_map_bgzf_hash.paf'
awk -F '\\t' '
BEGIN { ok=0 }
NF>=12 {
  if (\$1 !~ /^EXAMPLE_ASM_B_v2#1#/) exit 30
  if (\$6 !~ /^EXAMPLE_ASM_A_v1#1#/) exit 31
  ok=1
}
END { if (!ok) exit 32 }
' '$SYN/results/impg_map_bgzf_hash.paf'"

# Machine-readable summary of every recorded probe command.  PASS means exit 0;
# content-specific assertions above prevent a mere successful exit from being enough.
{
    printf 'case\texit_status\tresult\n'
    for f in "$SYN"/logs/*.exit_status.txt; do
        label=$(basename "$f" .exit_status.txt)
        rc=$(tr -d '\n' <"$f")
        if [[ $rc == 0 ]]; then result=PASS; else result=FAIL; fi
        printf '%s\t%s\t%s\n' "$label" "$rc" "$result"
    done
} >"$SYN/results/status_matrix.tsv"

# Enforce assignment safety limits on retained artifacts.
fasta_count=$(find "$SYN" -type f \( -name '*.fa' -o -name '*.fa.gz' \) | wc -l)
retained_bytes=0
for _ in 1 2 3; do
    {
        printf 'synthetic_assemblies=2\n'
        printf 'retained_fasta_files=%s\n' "$fasta_count"
        printf 'retained_probe_bytes=%s\n' "$retained_bytes"
        printf 'fasta_limit=5\n'
        printf 'byte_limit=104857600\n'
    } >"$SYN/results/safety_limits.txt"
    measured_bytes=$(du -sb "$SYN" | awk '{print $1}')
    [[ $measured_bytes == "$retained_bytes" ]] && break
    retained_bytes=$measured_bytes
done
[[ $fasta_count -le 5 ]]
[[ $retained_bytes -le 104857600 ]]

printf 'PASS: bounded BGZF x literal-# probe; see %s/results/status_matrix.tsv\n' "$SYN"
