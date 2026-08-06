#!/usr/bin/env python3
"""verify_triangle.py — structural + spot-check validation of the host triangle.

Checks:
  1. File size == n*(n-1)/2 * 4 bytes.
  2. Spot-check 50 random pairs (i<j) against direct `mash dist` on the full
     hosts.msh sketch (extract those accessions' sequences and re-sketch).
  3. Zero-count (exact duplicates).
  4. Row-0 sanity (values in [0,1]).
  5. Sample stats.
Writes: triangle_verify.json
"""
import json, os, random, subprocess, sys, argparse, glob, gzip
import numpy as np

def off(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--dist", required=True)
    ap.add_argument("--msh", required=True)
    ap.add_argument("--src", required=True, help="canonical_objects dir for spot-check seqs")
    ap.add_argument("--mash", default="mash")
    ap.add_argument("--out", default="triangle_verify.json")
    args = ap.parse_args()

    ids = [l.strip() for l in open(args.ids) if l.strip()]
    n = len(ids)
    tri_len = n * (n - 1) // 2
    size = os.path.getsize(args.dist)
    ok_size = size == tri_len * 4
    print(f"[verify] n={n} size={size} expected={tri_len*4} match={ok_size}")
    if not ok_size:
        sys.exit(1)

    tri = np.memmap(args.dist, dtype="<f4", mode="r", shape=(tri_len,))
    row0 = tri[:n - 1]
    print(f"[verify] row0 min={row0.min():.4f} max={row0.max():.4f} "
          f"in_range={row0.min() >= 0.0 and row0.max() <= 1.0}")
    nz = int((tri == 0.0).sum())
    print(f"[verify] pairs with distance 0.0: {nz} ({nz/tri_len*100:.3f}%)")

    rng = np.random.default_rng(42)
    offs = rng.integers(0, tri_len, size=1_000_000)
    sample = tri[offs]
    print(f"[verify] sample(1M) mean={sample.mean():.4f} median={np.median(sample):.4f} "
          f"frac==1.0={(sample==1.0).mean():.4f}")

    # spot check 50 random distinct pairs vs direct mash dist on re-sketch
    random.seed(7)
    pairs = set()
    while len(pairs) < 50:
        i, j = sorted(random.sample(range(n), 2))
        pairs.add((i, j))
    pairs = sorted(pairs)
    want = set()
    for i, j in pairs:
        want.add(ids[i]); want.add(ids[j])

    # build a single fasta of the wanted accessions from the source pansn files
    seqs = {}
    for acc in want:
        p = glob.glob(os.path.join(args.src, acc, "*.pansn.fa.gz"))
        if not p:
            print(f"[verify] missing source for {acc}"); sys.exit(1)
        chars = []
        with gzip.open(p[0], "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    continue
                chars.append(line.rstrip())
        seqs[acc] = "".join(chars).upper()
    fa = "/tmp/host_verify_pairs.fa"
    with open(fa, "w") as f:
        for acc in sorted(seqs):
            f.write(f">{acc}\n{seqs[acc]}\n")
    msh = "/tmp/host_verify_pairs.msh"
    subprocess.run([args.mash, "sketch", "-i", "-k", "21", "-s", "10000",
                    "-o", msh, fa], check=True, capture_output=True)
    out = subprocess.run([args.mash, "dist", "-p", "32", msh, msh],
                         capture_output=True, text=True).stdout
    d = {}
    for line in out.splitlines():
        a, b, dist, _, _ = line.split("\t")
        d[(a, b)] = float(dist)
    bad = 0
    for i, j in pairs:
        got = float(tri[off(i, j, n)])
        exp = d.get((ids[i], ids[j])) or d.get((ids[j], ids[i]))
        if exp is None:
            print(f"[verify] pair ({ids[i]},{ids[j]}) missing")
            bad += 1
        elif abs(got - exp) > 1e-6:
            print(f"[verify] MISMATCH {ids[i]} {ids[j]} got={got} exp={exp}")
            bad += 1
    print(f"[verify] spot-checked {len(pairs)} pairs vs direct mash: mismatches={bad}")

    summary = {
        "n": n, "tri_len": tri_len, "size_bytes": size, "size_ok": ok_size,
        "row0_min": float(row0.min()), "row0_max": float(row0.max()),
        "zero_distance_pairs": int(nz),
        "sample_mean": float(sample.mean()),
        "sample_median": float(np.median(sample)),
        "sample_frac_one": float((sample == 1.0).mean()),
        "spot_check_pairs": len(pairs), "spot_check_mismatches": bad,
    }
    with open(args.out, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    sys.exit(0 if (ok_size and bad == 0) else 1)

if __name__ == "__main__":
    main()