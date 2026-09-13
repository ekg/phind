#!/usr/bin/env python3
"""validate_mash_clades_v3.py — independent validation for task ntm-v3-prophage.

Prior art: ntm/v2/scripts/validate_mash_clades.py (v2). Checks mirror the
task's Validation section:

  1. sketch integrity: `mash info v3_prophages.msh` exits 0 and reports
     Sketches == n == 9,446 (with --mash-triangle: additionally a full
     `mash triangle` re-run exits 0; off by default because the driver's
     spotcheck already re-measured 50 pairs from sequence data).
  2. row count == prophage count: ids.txt lines == FASTA headers, same
     order, unique; v3_prophages.dist byte size == n*(n-1)/2 * 4 (float32
     upper triangle).
  3. readback spot-check: sampled pairs from v3_prophages.triangle.txt agree
     with the float32 triangle (offset a*(2n-a-1)/2+(b-a-1)).
  4. tight clade assignment: every prophage assigned exactly once (sum of
     tight_clades.json member counts == n; ids match ids.txt, no dupes).
  5. per-clade internal median MASH <= 0.25 for every non-singleton clade
     (clade_similarity.json); distribution reported.
  6. repo clade definitions: ntm/v3/clades/tight_clades.json.gz decompresses
     byte-identically to the NVMe tight_clades.json; SHA256SUMS verifies;
     clade_summary.tsv row count == clade count.
  7. prophage length band: lengths recomputed from full_prophages.fa match
     the extract_report references (n 9,446, min 237, max 93,713,
     total 206,284,318, 91 outside [1,000, 100,000]).

Exit 0 if all pass, 1 otherwise.
"""
import argparse
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys

import numpy as np

WORK = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3"
FA = f"{WORK}/full_prophages.fa"
OUT = f"{WORK}/mash_clades"
CLADES = f"{WORK}/clades"
THRESHOLD = 0.25
BAND = (1000, 100000)
EXTRACT_REF = {"n": 9446, "len_min": 237, "len_max": 93713,
               "len_total": 206284318, "outside_band": 91}


def offset(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    ap.add_argument("--mash-triangle", action="store_true",
                    help="also re-run `mash triangle` on the sketch (exit "
                         "code check; ~minutes)")
    args = ap.parse_args()

    errors, info = [], []

    # 1. sketch integrity
    r = subprocess.run(["mash", "info", f"{OUT}/v3_prophages.msh"],
                       capture_output=True, text=True)
    n_sketch = None
    if r.returncode == 0:
        for line in r.stdout.splitlines():
            ls = line.strip()
            if ls.startswith("Sketches:"):
                val = ls[len("Sketches:"):].strip()
                if val.isdigit():
                    n_sketch = int(val)
                    break
        info.append(f"mash info: exit 0, sketches={n_sketch}")
    else:
        errors.append(f"mash info: exit {r.returncode}")
    if args.mash_triangle:
        rt = subprocess.run(["mash", "triangle", "-p", "32",
                             f"{OUT}/v3_prophages.msh"],
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
        (info if rt.returncode == 0 else errors).append(
            f"mash triangle re-run: exit {rt.returncode}")

    # 2. row count == prophage count
    with open(f"{OUT}/ids.txt") as f:
        ids = [l.strip() for l in f if l.strip()]
    n = len(ids)
    fh = []
    lengths = []
    cur = None
    with open(FA) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    lengths.append(cur)
                fh.append(line[1:].strip().split()[0])
                cur = 0
            else:
                cur += len(line.strip())
    if cur is not None:
        lengths.append(cur)
    if n != len(fh):
        errors.append(f"ids.txt rows {n} != FASTA headers {len(fh)}")
    else:
        info.append(f"dist matrix rows == prophage count == {n}")
    if set(ids) != set(fh):
        errors.append("ids.txt ids != FASTA header ids")
    if ids != sorted(ids):
        errors.append("ids.txt is not in sorted() order "
                      "(build_tight_clades row-order contract; v3 "
                      "convention — see report mechanics note)")
    if len(set(ids)) != n:
        errors.append(f"duplicate ids in ids.txt ({n} rows, "
                      f"{len(set(ids))} unique)")
    if n_sketch is not None and n_sketch != n:
        errors.append(f"sketch sketches {n_sketch} != ids.txt rows {n}")

    npairs = n * (n - 1) // 2
    tri_path = f"{OUT}/v3_prophages.dist"
    size = os.path.getsize(tri_path)
    if size != npairs * 4:
        errors.append(f"triangle size {size} != {npairs * 4}")
    else:
        info.append(f"triangle float32 {npairs} values ({size} bytes) == "
                    f"n*(n-1)/2*4 OK")

    # 3. readback spot-check vs the mash triangle text (rows are in mash's
    # own sketch order — map by NAME to ids.txt indices, as the driver does)
    idx = {s: i for i, s in enumerate(ids)}
    rng = random.Random(7)
    want = set()
    n_sample = min(100, n * (n - 1) // 2)
    while len(want) < n_sample:
        i = rng.randrange(1, n)
        j = rng.randrange(0, i)
        want.add((i, j))
    by_row = {}
    for i, j in want:
        by_row.setdefault(i, []).append(j)
    got = {}
    textpos_to_ididx = np.empty(n, dtype=np.int64)
    with open(f"{OUT}/v3_prophages.triangle.txt") as f:
        head = f.readline().rstrip("\n")
        if not (head.startswith("\t") and head[1:].isdigit()
                and int(head[1:]) == n):
            errors.append(f"triangle text header {head!r} != \\t{n}")
        rows = 0
        for line in f:
            toks = line.rstrip("\n").split("\t")
            if toks[0] not in idx:
                errors.append(f"triangle text row {rows} name unknown")
                break
            b = idx[toks[0]]
            textpos_to_ididx[rows] = b
            js = by_row.get(rows)
            if js:
                for j in js:
                    got[(rows, j)] = (b, int(textpos_to_ididx[j]),
                                      float(toks[1 + j]))
            rows += 1
        if rows != n:
            errors.append(f"triangle text rows {rows} != {n}")
    mm = np.memmap(tri_path, dtype="<f4", mode="r", shape=(npairs,))
    mismatches = 0
    for (i, j), (b, a, d) in got.items():
        if a == b:
            continue
        if abs(float(mm[offset(min(a, b), max(a, b), n)]) - d) > 1e-6:
            mismatches += 1
            if mismatches <= 3:
                errors.append(f"readback mismatch {(a, b)}: f32 "
                              f"{float(mm[offset(min(a, b), max(a, b), n)])} "
                              f"vs text {d}")
    if mismatches == 0:
        info.append(f"readback: {len(got)} sampled pairs match triangle text")
    del mm

    # 4. assignment completeness
    with open(f"{CLADES}/0/tight_clades.json") as f:
        tc = json.load(f)
    all_members = [m for members in tc.values() for m in members]
    if len(all_members) != n:
        errors.append(f"assignment incomplete: {len(all_members)} != {n}")
    elif len(set(all_members)) != n:
        errors.append(f"duplicate assignments: {len(set(all_members))} "
                      f"unique of {len(all_members)}")
    elif set(all_members) != set(ids):
        errors.append("clade member ids != ids.txt ids")
    else:
        info.append(f"assignment complete: {len(all_members)} members across "
                    f"{len(tc)} clades, all unique, == ids.txt")

    # 5. per-clade stats distribution
    with open(f"{CLADES}/0/clade_similarity.json") as f:
        sim = json.load(f)["per_clade"]
    medians = [s["median"] for s in sim.values() if s.get("median") is not None]
    over = [cid for cid, s in sim.items()
            if s.get("n", 1) > 1 and s.get("median", 0) > args.threshold + 1e-6]
    n_clades = len(tc)
    n_singletons = sum(1 for s in sim.values() if s["n"] == 1)
    n_alignable = n_clades - n_singletons
    max_size = max(s["n"] for s in sim.values())
    info.append(f"clades={n_clades} alignable(>=2)={n_alignable} "
                f"singletons={n_singletons} max_size={max_size}")
    info.append(f"internal median mash: median={np.median(medians):.6f} "
                f"min={np.min(medians):.6f} max={np.max(medians):.6f} "
                f"(v2 median: 0.0060; v1: 0.074)")
    if over:
        errors.append(f"non-singleton clades with median > threshold: {over}")
    else:
        info.append(f"every non-singleton clade median <= {args.threshold}")
    if max_size > 100:
        errors.append(f"max clade size {max_size} > 100 cap")

    # 6. repo clade definitions
    # <repo>/ntm/v3/scripts/validate_mash_clades_v3.py -> <repo>
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    repo_clades = os.path.join(repo_root, "ntm", "v3", "clades")
    gz = os.path.join(repo_clades, "tight_clades.json.gz")
    if os.path.exists(gz):
        with open(gz, "rb") as f:
            dec = gzip.decompress(f.read())
        with open(f"{CLADES}/0/tight_clades.json", "rb") as f:
            ref = f.read()
        if dec == ref:
            info.append(f"repo ntm/v3/clades/{os.path.basename(gz)} "
                        f"decompresses byte-identical to NVMe "
                        f"tight_clades.json ({len(ref):,} bytes)")
        else:
            errors.append("repo tight_clades.json.gz != NVMe tight_clades.json")
        sums = os.path.join(repo_clades, "SHA256SUMS")
        if os.path.exists(sums):
            bad = []
            for line in open(sums):
                h, name = line.split()
                path = os.path.join(repo_clades, name)
                d = hashlib.sha256(open(path, "rb").read()).hexdigest()
                if d != h:
                    bad.append(name)
            if bad:
                errors.append(f"SHA256SUMS mismatch: {bad}")
            else:
                info.append("repo SHA256SUMS verifies (4 files)")
        n_summary = sum(1 for _ in open(
            os.path.join(repo_clades, "clade_summary.tsv"))) - 1
        if n_summary != n_clades:
            errors.append(f"clade_summary.tsv rows {n_summary} != clades "
                          f"{n_clades}")
        else:
            info.append(f"clade_summary.tsv rows == clade count == {n_clades}")
    else:
        errors.append(f"missing repo clade definitions {gz}")

    # 7. prophage length band (recomputed from the FASTA)
    outside = sum(1 for L in lengths if not (BAND[0] <= L <= BAND[1]))
    lencheck = {
        "n": len(lengths), "len_min": min(lengths), "len_max": max(lengths),
        "len_total": sum(lengths), "outside_band": outside,
    }
    bad = {k: (v, EXTRACT_REF[k]) for k, v in lencheck.items()
           if v != EXTRACT_REF[k]}
    if bad:
        errors.append(f"length band mismatch vs extract_report: {bad}")
    else:
        info.append(f"length band: n={len(lengths)}, min={min(lengths)}, "
                    f"max={max(lengths)}, total={sum(lengths):,}, "
                    f"outside [{BAND[0]},{BAND[1]}] = {outside} "
                    f"(== extract_report)")

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
