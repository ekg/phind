#!/usr/bin/env python3
"""extract_clade_fasta.py — extract per-clade FASTA from full_prophages.fa.

Usage:
  extract_clade_fasta.py --clades research/clades/10/tight_clades.json \
      --fasta prophage_homology_survey/full_prophages.fa \
      --outdir research/clades/10 --clade 10_0134 [--clade 10_0000 ...]

If --clade is omitted, all clades in the JSON are extracted. Output:
<outdir>/<clade_id>/sequences.fa.
"""
import argparse
import json
import os
import sys


def extract(fasta_path, members, out_path):
    want = set(members)
    n = 0
    with open(fasta_path) as fin, open(out_path, "w") as fout:
        cur = None
        buf = []
        for line in fin:
            if line.startswith(">"):
                if cur is not None and cur in want:
                    fout.write(">" + cur + "\n" + "".join(buf) + "\n")
                    n += 1
                cur = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line)
        if cur is not None and cur in want:
            fout.write(">" + cur + "\n" + "".join(buf) + "\n")
            n += 1
    assert n == len(members), (n, len(members))
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clades", required=True, help="tight_clades.json path")
    ap.add_argument("--fasta", required=True, help="full_prophages.fa path")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--clade", action="append", help="clade id(s); default all")
    args = ap.parse_args()

    tc = json.load(open(args.clades))
    ids = args.clade or list(tc.keys())
    for cid in ids:
        outdir = os.path.join(args.outdir, cid)
        os.makedirs(outdir, exist_ok=True)
        n = extract(args.fasta, tc[cid], os.path.join(outdir, "sequences.fa"))
        print(f"{cid}: {n} sequences -> {outdir}/sequences.fa", flush=True)


if __name__ == "__main__":
    sys.exit(main())
