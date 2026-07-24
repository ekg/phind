#!/usr/bin/env bash
# Reproduce the read-only installed-binary inspection and the bounded ordinary
# synthetic FASTA/SYNG smoke probe for reports/impg_syng_assessment.md.
# No production inputs are read. The script recreates only artifacts below
# artifacts/impg_probe/synthetic/ and the three evidence text files beside it.
set -u

ROOT="$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)"
OUT="$ROOT/artifacts/impg_probe"
SYN="$OUT/synthetic"
IMPG=/home/erikg/.cargo/bin/impg
SOURCE=/home/erikg/impg
ACCESS_DATE=2026-07-24
mkdir -p "$OUT"
rm -rf "$SYN"
mkdir -p "$SYN"

record() {
    local destination=$1
    shift
    {
        printf '$'
        printf ' %q' "$@"
        printf '\n'
        "$@" 2>&1
        local rc=$?
        printf '[exit_status=%d]\n' "$rc"
        return "$rc"
    } >>"$destination"
}

: >"$OUT/environment.txt"
{
    printf 'Evidence access date: %s\n' "$ACCESS_DATE"
    printf 'Command script: artifacts/impg_probe/commands.sh\n'
} >>"$OUT/environment.txt"
record "$OUT/environment.txt" date --iso-8601=seconds || true
record "$OUT/environment.txt" uname -a || true
record "$OUT/environment.txt" bash --version || true
record "$OUT/environment.txt" lscpu || true
record "$OUT/environment.txt" free -h || true
record "$OUT/environment.txt" df -h "$ROOT" "$SYN" || true
record "$OUT/environment.txt" rustc --version --verbose || true
record "$OUT/environment.txt" cargo --version --verbose || true
record "$OUT/environment.txt" /usr/bin/time --version || true
# Restrict samtools output to stable version lines; some distro build-path bytes
# in its verbose compilation section are not valid UTF-8.
record "$OUT/environment.txt" bash -c 'samtools --version 2>&1 | head -3' || true
record "$OUT/environment.txt" bgzip --version || true
for tool in FastGA fastga wfmash sweepga seqwish agc gfasort gfaffix odgi; do
    {
        printf '$ command -v %q\n' "$tool"
        command -v "$tool" 2>&1
        rc=$?
        printf '[exit_status=%d]\n' "$rc"
    } >>"$OUT/environment.txt"
done
for tool in FastGA wfmash gfasort gfaffix; do
    record "$OUT/environment.txt" "$tool" --version || true
done

: >"$OUT/version.txt"
{
    printf 'Evidence access date: %s\n' "$ACCESS_DATE"
    printf 'Installed path requested: /home/erikg/.cargo/bin/impg\n'
} >>"$OUT/version.txt"
record "$OUT/version.txt" command -v impg || true
record "$OUT/version.txt" ls -ld "$IMPG" || true
record "$OUT/version.txt" readlink "$IMPG" || true
record "$OUT/version.txt" readlink -f "$IMPG" || true
record "$OUT/version.txt" file "$IMPG" || true
record "$OUT/version.txt" stat "$IMPG" || true
record "$OUT/version.txt" sha256sum "$IMPG" || true
record "$OUT/version.txt" "$IMPG" --version || true
record "$OUT/version.txt" ldd "$IMPG" || true
record "$OUT/version.txt" readelf -p .comment "$IMPG" || true
record "$OUT/version.txt" git -C "$SOURCE" remote -v || true
record "$OUT/version.txt" git -C "$SOURCE" rev-parse HEAD || true
record "$OUT/version.txt" git -C "$SOURCE" describe --tags --always --dirty || true
record "$OUT/version.txt" git -C "$SOURCE" status --short || true
record "$OUT/version.txt" sha256sum "$SOURCE/target/guix/release/impg" || true
{
    printf '$ python3 - (read Cargo installation registry entry)\n'
    python3 - <<'PY'
import json
p = '/home/erikg/.cargo/.crates2.json'
d = json.load(open(p))['installs']
for key, value in d.items():
    if key.startswith('impg '):
        print(key)
        for field in ('bins', 'features', 'profile', 'target', 'rustc'):
            print(f'{field}: {value.get(field)}')
PY
    rc=$?
    printf '[exit_status=%d]\n' "$rc"
} >>"$OUT/version.txt" 2>&1

: >"$OUT/help.txt"
for args in \
    "--help" \
    "index --help" \
    "lace --help" \
    "partition --help" \
    "query --help" \
    "refine --help" \
    "similarity --help" \
    "genotype --help" \
    "genotype cos --help" \
    "project --help" \
    "infer --help" \
    "stats --help" \
    "graph --help" \
    "normalize-self-loops --help" \
    "crush --help" \
    "gfa2vcf --help" \
    "graph-report --help" \
    "render --help" \
    "align --help" \
    "map --help" \
    "read-index --help" \
    "syng --help" \
    "syng2gfa --help" \
    "syng-repair --help"; do
    # Intentional shell word-splitting of this fixed, argument-free vocabulary.
    # shellcheck disable=SC2086
    record "$OUT/help.txt" "$IMPG" $args || true
done
record "$OUT/help.txt" "$IMPG" version || true

# Exactly four deterministic, tiny, ordinary-name synthetic FASTA records:
# three 1,200-bp panel paths and one 360-bp query. The latter is a subsequence
# of genomeA; genomeB has substitutions and genomeC has an insertion/deletion.
python3 - "$SYN" <<'PY'
from pathlib import Path
import random, sys
out = Path(sys.argv[1])
rng = random.Random(20260724)
a = ''.join(rng.choice('ACGT') for _ in range(1200))
b = list(a)
for i in range(37, len(b), 101):
    b[i] = {'A':'C','C':'G','G':'T','T':'A'}[b[i]]
b = ''.join(b)
c = a[:610] + 'GATTACAGATTACA' + a[624:]
assert len(a) == len(b) == len(c) == 1200
q = a[180:540]
with (out/'panel.fa').open('w') as f:
    for name, seq in [('genomeA', a), ('genomeB', b), ('genomeC', c)]:
        f.write(f'>{name} ordinary synthetic record\n')
        for i in range(0, len(seq), 60): f.write(seq[i:i+60]+'\n')
with (out/'probe.fa').open('w') as f:
    f.write('>probeA ordinary synthetic subsequence\n')
    for i in range(0, len(q), 60): f.write(q[i:i+60]+'\n')
PY

: >"$SYN/probe.log"
record "$SYN/probe.log" wc -c "$SYN/panel.fa" "$SYN/probe.fa" || true
record "$SYN/probe.log" /usr/bin/time -v "$IMPG" syng \
    -f "$SYN/panel.fa" -o "$SYN/serial.syng" \
    --syncmer-length 31 --smer-length 5 --position-sample-rate 16 -t 1 || true
record "$SYN/probe.log" /usr/bin/time -v "$IMPG" syng \
    -f "$SYN/panel.fa" -o "$SYN/parallel.syng" \
    --syncmer-length 31 --smer-length 5 --position-sample-rate 16 \
    --parallel-dictionary -t 2 || true
record "$SYN/probe.log" /usr/bin/time -v "$IMPG" syng \
    -f "$SYN/panel.fa" -o "$SYN/parallel_repeat.syng" \
    --syncmer-length 31 --smer-length 5 --position-sample-rate 16 \
    --parallel-dictionary -t 2 || true
record "$SYN/probe.log" bash -c \
    'for s in 1khash 1gbwt names pstep spos meta; do sha256sum "$1/serial.syng.$s" "$1/parallel.syng.$s" "$1/parallel_repeat.syng.$s"; done' \
    _ "$SYN" || true
record "$SYN/probe.log" bash -c \
    'for s in 1khash 1gbwt names pstep spos meta; do cmp -s "$1/parallel.syng.$s" "$1/parallel_repeat.syng.$s" || exit 1; done' \
    _ "$SYN" || true
record "$SYN/probe.log" bash -c \
    '"$1" query -a "$2/parallel.syng" -r genomeA:180-540 -d 50 --syng-raw -o bed -t 2 >"$2/query.bed"' \
    _ "$IMPG" "$SYN" || true
record "$SYN/probe.log" bash -c \
    '"$1" map -a "$2/parallel.syng" -q "$2/probe.fa" -o gaf -t 2 >"$2/map.gaf"' \
    _ "$IMPG" "$SYN" || true
record "$SYN/probe.log" bash -c \
    '"$1" map -a "$2/parallel.syng" -q "$2/probe.fa" -o paf --min-anchors 2 --max-hits 10 -t 2 >"$2/map.paf"' \
    _ "$IMPG" "$SYN" || true
record "$SYN/probe.log" "$IMPG" syng2gfa -a "$SYN/parallel.syng" \
    --sequence-files "$SYN/panel.fa" --gfa-mode raw -o "$SYN/panel.raw.gfa" -t 2 || true
record "$SYN/probe.log" bash -c \
    '"$1" query -a "$2/parallel.syng" -r genomeA:180-540 -d 50 --syng-raw -o fasta --sequence-files "$2/panel.fa" -t 2 >"$2/query.fa"' \
    _ "$IMPG" "$SYN" || true
record "$SYN/probe.log" wc -l "$SYN/query.bed" "$SYN/map.gaf" "$SYN/map.paf" "$SYN/panel.raw.gfa" "$SYN/query.fa" || true
record "$SYN/probe.log" du -ah "$SYN" || true
record "$SYN/probe.log" bash -c 'find "$1" -type f -printf "%s\t%P\n" | sort -n' _ "$SYN" || true
record "$SYN/probe.log" bash -c 'find "$1" -type f -printf "%s\n" | awk "{s+=\$1} END{print s}"' _ "$SYN" || true

printf 'Probe complete. Evidence under %s\n' "$OUT"
