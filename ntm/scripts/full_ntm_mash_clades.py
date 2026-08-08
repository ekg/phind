#!/usr/bin/env python3
"""
full_ntm_mash_clades.py — mash sketch + dist -> float32 upper-triangle (in
build_tight_clades format) + ids.txt + single-community labels.csv, then run
scripts/build_tight_clades.py over ALL 10,438 NTM prophages.

Single community (0): leader clustering forms tight clades naturally across
species; cross-species prophage clades (host range) are preserved and resolved
later by the host-clade join.

Outputs under ntm/v1/mash_clades/.
"""
import csv
import os
import struct
import subprocess
import sys

ROOT = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1"
FA = f"{ROOT}/full_prophages.fa"
OUT = f"{ROOT}/mash_clades"
REPO = "/home/erikg/phind"
BUILD_TIGHT = f"{REPO}/scripts/build_tight_clades.py"


def run(cmd):
    print("+", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(f"FAILED: {' '.join(cmd)}")


def main():
    os.makedirs(OUT, exist_ok=True)
    # 1. ids in file order
    ids = []
    with open(FA) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    print(f"{len(ids)} prophages", flush=True)

    # 2. mash sketch (individual) + dist
    run(["mash", "sketch", "-i", "-k", "21", "-s", "10000", "-o", f"{OUT}/prophages", FA])
    dist_tsv = f"{OUT}/prophages.dist.tsv"
    with open(dist_tsv, "w") as fdist:
        subprocess.run(["mash", "dist", f"{OUT}/prophages.msh", f"{OUT}/prophages.msh"],
                       stdout=fdist)
    print("mash dist done -> triangle", flush=True)

    # 3. float32 upper triangle in ids order
    idx = {s: i for i, s in enumerate(ids)}
    n = len(ids)
    pair = {}
    with open(dist_tsv) as f:
        for line in f:
            a, b, d = line.split("\t")[:3]
            ia, ib = idx.get(a), idx.get(b)
            if ia is None or ib is None or ia == ib:
                continue
            pair[(min(ia, ib), max(ia, ib))] = float(d)

    def off(a, b):
        return a * (2 * n - a - 1) // 2 + (b - a - 1)

    tri = f"{OUT}/prophages_mash.dist"
    with open(tri, "wb") as f:
        for a in range(n):
            for b in range(a + 1, n):
                f.write(struct.pack("<f", pair.get((a, b), 1.0)))
    with open(f"{OUT}/ids.txt", "w") as f:
        f.write("\n".join(ids) + "\n")
    with open(f"{OUT}/labels.csv", "w") as f:
        f.write("sequence,community\n")
        for s in ids:
            f.write(f"{s},0\n")
    print(f"triangle {n} seqs {n*(n-1)//2} pairs -> {tri}", flush=True)

    # 4. tight clades
    run(["python3", BUILD_TIGHT,
         "--threshold", "0.25", "--max-size", "100", "--communities", "0",
         "--outdir", f"{OUT}/clades", "--ids-file", f"{OUT}/ids.txt",
         "--triangle", tri, "--labels-csv", f"{OUT}/labels.csv"])
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
