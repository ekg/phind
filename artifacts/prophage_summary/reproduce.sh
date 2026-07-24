#!/bin/sh
# Regenerate the prophage audit offline from the two immutable root inputs.
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=${1:-$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)}
CSV="$ROOT/26k_prophage1.csv"
ACCESSIONS="$ROOT/26k_ecoli_accession.txt"
EXPECTED_CSV_SHA256=6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996
EXPECTED_ACCESSIONS_SHA256=1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5

hash_file() {
    sha256sum "$1" | awk '{print $1}'
}

[ -f "$CSV" ] || { echo "Missing input: $CSV" >&2; exit 1; }
[ -f "$ACCESSIONS" ] || { echo "Missing input: $ACCESSIONS" >&2; exit 1; }

csv_before=$(hash_file "$CSV")
accessions_before=$(hash_file "$ACCESSIONS")
[ "$csv_before" = "$EXPECTED_CSV_SHA256" ] || {
    echo "Unexpected 26k_prophage1.csv SHA-256: $csv_before" >&2; exit 1;
}
[ "$accessions_before" = "$EXPECTED_ACCESSIONS_SHA256" ] || {
    echo "Unexpected 26k_ecoli_accession.txt SHA-256: $accessions_before" >&2; exit 1;
}

export LC_ALL=C
export PYTHONHASHSEED=0
python3 "$SCRIPT_DIR/analyze.py" "$ROOT"

csv_after=$(hash_file "$CSV")
accessions_after=$(hash_file "$ACCESSIONS")
[ "$csv_before" = "$csv_after" ] || {
    echo "26k_prophage1.csv changed during reproduction" >&2; exit 1;
}
[ "$accessions_before" = "$accessions_after" ] || {
    echo "26k_ecoli_accession.txt changed during reproduction" >&2; exit 1;
}
printf '%s\n' "Reproduced prophage audit; input SHA-256 values unchanged."
