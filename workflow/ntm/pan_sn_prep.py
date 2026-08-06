#!/usr/bin/env python3
"""PanSN-rename public NCBI FASTA -> {acc}.pansn.fa.gz + .fai + .gzi.

Mirrors the convention used by download-7-352-ntm-genomes / the E. coli
26k canonical_objects (contig header {acc}#1#{orig_contig}, bgzip + faidx).
Core: {acc} = assembly accession WITHOUT versionified filename differences;
contig token is the unversioned contig id (matches Phigaro scaffold naming).

Usage:
  python workflow/ntm/pan_sn_prep.py \
      --in-dir /tmp/ntm_dl \
      --out-dir ntm/v1/genomes/canonical_objects
"""
import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys

PANS_PREFIX_RE = re.compile(r"^>([^\s]+)")


def iter_records(fna_path):
    """Yield (header_id, seq) from a gz or plain FASTA."""
    opener = gzip.open if fna_path.endswith(".gz") else open
    with opener(fna_path, "rt") as fh:
        header = None
        seq = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None and seq:
                    yield header, "".join(seq)
                m = PANS_PREFIX_RE.match(line)
                header = m.group(1) if m else line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None and seq:
            yield header, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for acc_dir in sorted(os.listdir(args.in_dir)):
        data_dir = os.path.join(args.in_dir, acc_dir, "ncbi_dataset", "data", acc_dir)
        if not os.path.isdir(data_dir):
            print(f"SKIP {acc_dir}: no data dir", file=sys.stderr)
            continue
        fna = [f for f in os.listdir(data_dir) if f.endswith("_genomic.fna")]
        if not fna:
            print(f"SKIP {acc_dir}: no _genomic.fna", file=sys.stderr)
            continue
        fna = os.path.join(data_dir, fna[0])

        out_stem = os.path.join(args.out_dir, f"{acc_dir}.pansn")
        out_fa = out_stem + ".fa"
        n_contigs = 0
        total_bp = 0
        with open(out_fa, "w") as out:
            for orig, seq in iter_records(fna):
                # Unversioned contig token: strip trailing version (.1 etc.)
                token = orig.split(".")[0]
                header = f"{acc_dir}#1#{token}"
                out.write(f">{header}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")
                n_contigs += 1
                total_bp += len(seq)

        # bgzip then remove plain fa
        subprocess.run(["bgzip", "-f", out_fa], check=True)
        gz = out_fa + ".gz"
        subprocess.run(["samtools", "faidx", gz], check=True)
        print(f"OK {acc_dir}: {n_contigs} contigs, {total_bp:,} bp -> {gz}")


if __name__ == "__main__":
    main()