#!/usr/bin/env python3
"""
abscessus_subset_pipeline.py — fast-track end-to-end on the abscessus prophages
already called, to produce a first NTM ML phage genome result while the full
7,303-genome geNomad call continues.

Steps:
  1. Extract every abscessus prophage (geNomad provirus.tsv -> samtools faidx)
     -> subset_abscessus/prophages.fa
  2. mash sketch + dist -> float32 upper-triangle (build_tight_clades format)
     + ids.txt + single-community labels.csv
  3. scripts/build_tight_clades.py -> tight_clades.json (validated clustering)

Outputs under /mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/subset_abscessus/.
"""
import csv
import glob
import os
import struct
import subprocess
import sys

ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1"
MAN = f"{ROOT}/accessions/ntm_accession_manifest.tsv"
PG = f"{ROOT}/prophages/per_genome"
GENOMES = f"{ROOT}/genomes/canonical_objects"
OUT = f"{ROOT}/subset_abscessus"
TARGET_SPECIES = "Mycobacteroides abscessus"
REPO = "/home/erikg/phind"
BUILD_TIGHT = f"{REPO}/scripts/build_tight_clades.py"


def run(cmd):
    print("+", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(f"FAILED: {' '.join(cmd)}")


def main():
    os.makedirs(OUT, exist_ok=True)
    sp = {r["primary_acc"]: r["species"]
          for r in csv.DictReader(open(MAN), delimiter="\t")}

    # 1. extract abscessus prophages
    fa_path = f"{OUT}/prophages.fa"
    ids = []
    with open(fa_path, "w") as fout:
        for pv in sorted(glob.glob(f"{PG}/*_genomad/*_find_proviruses/*_provirus.tsv")):
            acc = os.path.basename(pv).replace(".pansn_provirus.tsv", "")
            if sp.get(acc) != TARGET_SPECIES:
                continue
            gz = f"{GENOMES}/{acc}/{acc}.pansn.fa.gz"
            if not os.path.exists(gz):
                continue
            rows = [r for r in csv.DictReader(open(pv), delimiter="\t")]
            for i, r in enumerate(rows, 1):
                seq = r["source_seq"]
                beg, end = int(r["start"]), int(r["end"])
                if beg > end:
                    beg, end = end, beg
                region = f"{seq}:{beg}-{end}"
                pid = f"{acc}_prophage_{i}"
                res = subprocess.run(
                    ["samtools", "faidx", gz, region],
                    capture_output=True, text=True)
                if res.returncode != 0 or not res.stdout.strip():
                    continue
                # rewrite header to pid, keep sequence
                lines = res.stdout.splitlines()
                seq_body = "".join(l for l in lines[1:] if l)
                if not seq_body:
                    continue
                fout.write(f">{pid}\n{seq_body}\n")
                ids.append(pid)
    print(f"extracted {len(ids)} abscessus prophages -> {fa_path}", flush=True)
    if len(ids) < 20:
        raise SystemExit("too few prophages; wait for more abscessus to be called")

    # 2. mash sketch + dist -> triangle
    run(["mash", "sketch", "-i", "-k", "21", "-s", "10000", "-o", f"{OUT}/prophages",
         fa_path])
    dist_tsv = f"{OUT}/prophages.dist.tsv"
    with open(dist_tsv, "w") as fdist:
        subprocess.run(["mash", "dist", f"{OUT}/prophages.msh", f"{OUT}/prophages.msh"],
                       stdout=fdist)
    # parse -> float32 upper triangle in `ids` order
    idx = {s: i for i, s in enumerate(ids)}
    n = len(ids)
    pair = {}
    for line in open(dist_tsv):
        a, b, d = line.split("\t")[:3]
        ia, ib = idx[a], idx[b]
        if ia == ib:
            continue
        pair[(min(ia, ib), max(ia, ib))] = float(d)

    def off(a, b):
        return a * (2 * n - a - 1) // 2 + (b - a - 1)

    tri_path = f"{OUT}/prophages_mash.dist"
    with open(tri_path, "wb") as f:
        for a in range(n):
            for b in range(a + 1, n):
                f.write(struct.pack("<f", pair.get((a, b), 1.0)))
    with open(f"{OUT}/ids.txt", "w") as f:
        f.write("\n".join(ids) + "\n")
    # single-community labels.csv (build_tight_clades reads a `community` col)
    with open(f"{OUT}/labels.csv", "w") as f:
        f.write("sequence,community\n")
        for s in ids:
            f.write(f"{s},0\n")
    print(f"triangle: {n} seqs, {n*(n-1)//2} pairs -> {tri_path}", flush=True)

    # 3. build_tight_clades (reuse validated E. coli script)
    run(["python3", BUILD_TIGHT,
         "--threshold", "0.25", "--max-size", "100",
         "--communities", "0",
         "--outdir", f"{OUT}/clades",
         "--ids-file", f"{OUT}/ids.txt",
         "--triangle", tri_path,
         "--labels-csv", f"{OUT}/labels.csv"])
    print("DONE abscessus tight clades ->", f"{OUT}/clades", flush=True)


if __name__ == "__main__":
    main()
