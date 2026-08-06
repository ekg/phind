#!/usr/bin/env python3
"""Normalize Phigaro v2.4.0 native TSV -> E. coli-compatible prophage table.

Phigaro v2.4.0 emits 0-based inclusive coordinates (project `C2` convention).
The E. coli normalized table uses 1-based closed coordinates (`C1`), i.e.
[begin, end] inclusive, with canonical transform [begin-1, end) 0-based
half-open. This script therefore converts begin/end to 1-based closed and
assigns prophage_id = {accession}_prophage_N, matching 26k_prophage1.csv.

Columns (matches E. coli schema):
  genome, scaffold, begin, end, transposable, taxonomy, prophage_id, length

Also validates every (genome, scaffold, begin, end) is resolvable against the
PanSN bgzip FASTA via `samtools faidx`, and reports length statistics.
"""
import argparse
import csv
import glob
import os
import subprocess
import sys


HEADER = ["genome", "scaffold", "begin", "end", "transposable",
          "taxonomy", "prophage_id", "length"]


def phigaro_rows(tsv):
    """Read a Phigaro native TSV (skip header)."""
    with open(tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if not r.get("scaffold"):
                continue
            yield r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv-glob", required=True,
                    help="glob of Phigaro *.phigaro.tsv files")
    ap.add_argument("--genomes-dir", required=True,
                    help="dir of {acc}.pansn.fa.gz PanSN FASTAs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    for tsv in sorted(glob.glob(args.tsv_glob)):
        genome = os.path.basename(tsv).rsplit(".phigaro.tsv", 1)[0]
        n = 0
        for r in phigaro_rows(tsv):
            begin_c2 = int(r["begin"])  # 0-based inclusive
            end_c2 = int(r["end"])
            # C2 -> C1 1-based closed
            begin = begin_c2 + 1
            end = end_c2 + 1
            length = end - begin + 1
            n += 1
            rows.append({
                "genome": genome,
                "scaffold": r["scaffold"],
                "begin": begin,
                "end": end,
                "transposable": float(r["transposable"] == "True"),
                "taxonomy": r["taxonomy"],
                "prophage_id": f"{genome}_prophage_{n}",
                "length": length,
            })

    # ---- coordinate resolvability check against PanSN FASTA ----
    resolvable = 0
    failed = []
    for row in rows:
        fa = os.path.join(args.genomes_dir, f"{row['genome']}.pansn.fa.gz")
        if not os.path.exists(fa):
            failed.append((row["prophage_id"], "FASTA missing"))
            continue
        region = f"{row['scaffold']}:{row['begin']}-{row['end']}"
        p = subprocess.run(["samtools", "faidx", fa, region],
                           capture_output=True, text=True)
        if p.returncode == 0 and len(p.stdout.strip()) > 0:
            resolvable += 1
        else:
            failed.append((row["prophage_id"],
                           (p.stderr or "").strip()[:120] or "no seq"))

    # ---- write normalized table ----
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    # ---- stats ----
    lens = sorted(r["length"] for r in rows)
    n = len(lens)
    median = lens[n // 2] if n else 0
    print(f"Genomes with >0 prophages: "
          f"{len(set(r['genome'] for r in rows))}/{len(set(
              os.path.basename(x).split('.')[0]
              for x in glob.glob(os.path.join(args.genomes_dir, '*.pansn.fa.gz'))))}")
    print(f"Total prophages: {n}")
    print(f"Median length: {median} bp ({median/1000:.1f} kb)")
    print(f"Max length: {max(lens):,} bp ({max(lens)/1000:.1f} kb)"
          if lens else "Max: n/a")
    print(f"Min length: {min(lens):,} bp" if lens else "Min: n/a")
    print(f"Lengths (kb): {[f'{l/1000:.1f}' for l in lens]}")
    print(f"Resolution: {resolvable}/{n} resolvable "
          f"(0-based C2 -> 1-based C1 conversion applied)")
    for pid, err in failed:
        print(f"  FAIL {pid}: {err}")

    print(f"\nWrote {args.out} ({n} rows)")


if __name__ == "__main__":
    main()