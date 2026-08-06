#!/usr/bin/env python3
"""Consolidate geNomad provirus.tsv -> E. coli-compatible normalized table +
resolvability check against PanSN bgzip FASTA.

geNomad provirus.tsv columns: seq_name, source_seq, start, end, length,
n_genes, v_vs_c_score, in_seq_edge, integrases.
  - source_seq is the PanSN contig name (matches FASTA header).
  - start/end are 1-based closed (C1), matching the E. coli schema directly
    (prophage_id, genome=accession, scaffold=source_seq, begin=start, end=end).
  - prophage_id = {accession}_prophage_N.
"""
import argparse
import csv
import glob
import os
import subprocess

HEADER = ["genome", "scaffold", "begin", "end", "length",
          "n_genes", "v_vs_c_score", "taxonomy", "prophage_id"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv-glob", required=True)
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--taxonomy-suffix", default="",
                    help="optional extra taxonomy tsv per genome dir")
    args = ap.parse_args()

    rows = []
    for tsv in sorted(glob.glob(args.tsv_glob)):
        base = os.path.basename(tsv).rsplit("_provirus.tsv", 1)[0]
        genome = base
        # taxonomy from virus_summary.tsv if present alongside
        tax_map = {}
        summ = os.path.join(os.path.dirname(tsv), "summary",
                            f"{genome}_virus_summary.tsv")
        if os.path.exists(summ):
            with open(summ) as fh:
                sr = csv.DictReader(fh, delimiter="\t")
                for r in sr:
                    if r.get("seq_name") and r.get("taxonomy"):
                        tax_map[r["seq_name"]] = r["taxonomy"]
        n = 0
        with open(tsv) as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for rec in r:
                if not rec.get("seq_name"):
                    continue
                n += 1
                rows.append({
                    "genome": genome,
                    "scaffold": rec["source_seq"],
                    "begin": int(rec["start"]),   # 1-based closed
                    "end": int(rec["end"]),
                    "length": int(rec["length"]),
                    "n_genes": rec["n_genes"],
                    "v_vs_c_score": rec["v_vs_c_score"],
                    "taxonomy": tax_map.get(rec["seq_name"], ""),
                    "prophage_id": f"{genome}_prophage_{n}",
                })

    # resolvability
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
        if p.returncode == 0 and p.stdout.strip():
            resolvable += 1
        else:
            failed.append((row["prophage_id"], (p.stderr or "").strip()[:120]))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    lens = sorted(r["length"] for r in rows)
    n = len(lens)
    med = lens[n // 2] if n else 0
    print(f"Total prophages: {n}")
    print(f"Median length: {med} bp ({med/1000:.1f} kb)" if n else "Median: n/a")
    print(f"Max length: {max(lens):,} bp ({max(lens)/1000:.1f} kb)" if n else "")
    print(f"Min length: {min(lens):,} bp" if n else "")
    print(f"Lengths (kb): {[f'{l/1000:.1f}' for l in lens]}")
    print(f"Resolvable: {resolvable}/{n}")
    for pid, err in failed:
        print(f"  FAIL {pid}: {err}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()