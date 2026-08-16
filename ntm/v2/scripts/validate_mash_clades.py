#!/usr/bin/env python3
"""
validate_mash_clades.py — independent validation for task ntm-v2-prophage.

Checks (mirrors the task's Validation section):
  1. `mash triangle` on the sketch exits 0 (lower-tri matrix to a temp file,
     deleted afterwards).
  2. dist matrix row count == prophage count: ids.txt lines == FASTA header
     count; triangle byte size == n*(n-1)/2 * 4 (float32 upper triangle).
  3. Readback spot-check: sample pairs from the triangle agree with
     `mash dist` TSV (same float32 values, offset formula a*(2n-a-1)/2+(b-a-1)).
  4. Tight clade assignment: every prophage assigned exactly once
     (sum of tight_clades.json member counts == n; ids match, no dupes).
  5. Per-clade internal median mash distribution (v1: 0.074) reported from
     clade_similarity.json; every non-singleton clade median <= threshold 0.25.

Exit 0 if all pass, 1 otherwise.
"""
import argparse
import json
import os
import struct
import subprocess
import sys

import numpy as np

V2 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
FA = f"{V2}/full_prophages.fa"
OUT = f"{V2}/mash_clades"
CLADES_OUT = f"{V2}/clades"
THRESHOLD = 0.25


def load_ids():
    with open(f"{OUT}/ids.txt") as f:
        return [l.strip() for l in f if l.strip()]


def fasta_headers(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def offset(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    args = ap.parse_args()

    errors = []
    info = []

    # 1. mash triangle exit 0
    tri_tmp = f"{OUT}/prophages.triangle.phylip.tmp"
    r = subprocess.run(["mash", "triangle", "-p", "32",
                        f"{OUT}/prophages.msh"],
                       stdout=open(tri_tmp, "w"), stderr=subprocess.DEVNULL)
    if os.path.exists(tri_tmp):
        os.remove(tri_tmp)
    if r.returncode == 0:
        info.append("mash triangle: exit 0")
    else:
        errors.append(f"mash triangle: exit {r.returncode}")

    # 2. row count == prophage count
    ids = load_ids()
    n = len(ids)
    fh = fasta_headers(FA)
    if n != len(fh):
        errors.append(f"ids.txt rows {n} != FASTA headers {len(fh)}")
    else:
        info.append(f"dist matrix rows == prophage count == {n}")
    if ids != fh:
        errors.append("ids.txt order != FASTA header order")
    if len(set(ids)) != n:
        errors.append(f"duplicate ids in ids.txt ({n} rows, {len(set(ids))} unique)")

    npairs = n * (n - 1) // 2
    tri_path = f"{OUT}/prophages_mash.dist"
    if os.path.getsize(tri_path) != npairs * 4:
        errors.append(f"triangle size {os.path.getsize(tri_path)} != {npairs * 4}")
    else:
        info.append(f"triangle float32 {npairs} values ({npairs*4} bytes) OK")

    # 3. readback spot-check vs mash dist TSV
    idx = {s: i for i, s in enumerate(ids)}
    sample = {}
    with open(f"{OUT}/prophages.dist.tsv") as f:
        for i, line in enumerate(f):
            if i >= 200000:
                break
            a, b, d = line.split("\t")[:3]
            ia, ib = idx.get(a), idx.get(b)
            if ia is None or ib is None or ia == ib:
                continue
            key = (min(ia, ib), max(ia, ib))
            if len(sample) < 50:
                sample[key] = float(d)
    mm = np.memmap(tri_path, dtype="<f4", mode="r", shape=(npairs,))
    mismatches = 0
    for (a, b), d in sample.items():
        got = float(mm[offset(a, b, n)])
        if abs(got - d) > 1e-6:
            mismatches += 1
            if mismatches <= 3:
                errors.append(f"readback mismatch {(a, b)}: triangle {got} vs tsv {d}")
    if mismatches == 0:
        info.append(f"readback: {len(sample)} sampled pairs match dist TSV exactly")
    del mm

    # 4. assignment completeness
    with open(f"{CLADES_OUT}/0/tight_clades.json") as f:
        tc = json.load(f)
    all_members = []
    for cid, members in tc.items():
        all_members.extend(members)
    assigned = len(all_members)
    if assigned != n:
        errors.append(f"assignment incomplete: {assigned} != {n}")
    elif len(set(all_members)) != n:
        errors.append(f"duplicate assignments: {len(set(all_members))} unique of {assigned}")
    else:
        info.append(f"assignment complete: {assigned} members across {len(tc)} clades, all unique")

    # 5. per-clade stats distribution
    with open(f"{CLADES_OUT}/0/clade_similarity.json") as f:
        sim = json.load(f)["per_clade"]
    medians = [s["median"] for s in sim.values() if s.get("median") is not None]
    over = [cid for cid, s in sim.items()
            if s.get("n", 1) > 1 and s.get("median", 0) > args.threshold + 1e-6]
    n_clades = len(tc)
    n_singletons = sum(1 for s in sim.values() if s["n"] == 1)
    n_alignable = n_clades - n_singletons
    info.append(f"clades={n_clades} alignable(>=2)={n_alignable} "
                f"singletons={n_singletons}")
    info.append(f"internal median mash: median={np.median(medians):.6f} "
                f"min={np.min(medians):.6f} max={np.max(medians):.6f} "
                f"(v1 median: 0.074)")
    if over:
        errors.append(f"non-singleton clades with median > threshold: {over}")

    print("\n".join(info))
    if errors:
        print("ERRORS:", file=sys.stderr)
        for e in errors:
            print("  -", e, file=sys.stderr)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
