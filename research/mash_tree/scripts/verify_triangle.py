#!/usr/bin/env python3
"""verify_triangle.py — structural + spot-check validation of the triangle.

Checks:
  1. File size == n*(n-1)/2 * 4 bytes (every pair written).
  2. Spot-check N random pairs (i<j) against direct `mash dist` on the full
     sketch (allows a small tolerance for the spot-check set).
  3. Zero-count: number of pairs with distance exactly 0 (exact duplicates).
  4. Row 0 sanity: values in [0,1].
  5. Report full-matrix stats: mean/median distance, fraction of pairs at 1.0.
"""
import json, os, random, subprocess, sys
import numpy as np

def off(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)

def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ids = [l.strip() for l in open(os.path.join(base, "data", "ids.txt")) if l.strip()]
    n = len(ids)
    tri_path = os.path.join(base, "full_prophages_mash.dist")
    msh_path = os.path.join(base, "full_prophages.msh")

    tri_len = n * (n - 1) // 2
    size = os.path.getsize(tri_path)
    ok_size = size == tri_len * 4
    print(f"[verify] n={n} size={size} expected={tri_len*4} match={ok_size}")
    if not ok_size:
        sys.exit(1)

    tri = np.memmap(tri_path, dtype="<f4", mode="r", shape=(tri_len,))

    # row-0 sanity
    row0 = tri[:n - 1]
    print(f"[verify] row0 min={row0.min():.4f} max={row0.max():.4f} "
          f"in_range={row0.min() >= 0.0 and row0.max() <= 1.0}")

    # zero-distance pairs (exact duplicates)
    nz = int((tri == 0.0).sum())
    print(f"[verify] pairs with distance 0.0: {nz} "
          f"({nz/tri_len*100:.2f}% of triangle)")

    # overall stats on a uniform random sample of 1M pairs
    rng = np.random.default_rng(42)
    offs = rng.integers(0, tri_len, size=1_000_000)
    sample = tri[offs]
    print(f"[verify] sample(1M) mean={sample.mean():.4f} "
          f"median={np.median(sample):.4f} frac==1.0={(sample == 1.0).mean():.4f}")

    # spot check vs direct mash dist (random pairs, batch of 50)
    random.seed(7)
    pairs = [tuple(sorted(random.sample(range(n), 2))) for _ in range(50)]
    # build query file of pairs for mash? simpler: full mash dist on full sketch
    # is 17.5e9 lines — too big. Instead use mash dist with the sketch and the
    # two specific sequences via a small fasta of the 50 pairs' sequences.
    # Alternative: recompute via chunk sketches (already validated exhaustively
    # on the pilot). Here: check symmetry against recomputation from chunk
    # sketches is covered; do a direct fasta-based check for a few pairs.
    seqs = {}
    with open("/home/erikg/phind/prophage_homology_survey/full_prophages.fa") as f:
        cur = None
        for line in f:
            if line.startswith(">"):
                cur = line.strip()[1:].split()[0]
                if cur in [ids[i] for i, _ in pairs] or cur in [ids[j] for _, j in pairs]:
                    seqs[cur] = []
                else:
                    cur = None
            elif cur is not None:
                seqs[cur].append(line.strip())
    with open("/tmp/verify_pairs.fa", "w") as f:
        for name in seqs:
            f.write(f">{name}\n{''.join(seqs[name])}\n")
    subprocess.run(["mash", "sketch", "-i", "-o", "/tmp/verify_pairs.msh",
                    "/tmp/verify_pairs.fa"], check=True, capture_output=True)
    out = subprocess.run(["mash", "dist", "-p", "32", "/tmp/verify_pairs.msh",
                          "/tmp/verify_pairs.msh"], capture_output=True, text=True).stdout
    d = {}
    for line in out.splitlines():
        a, b, dist, _, _ = line.split("\t")
        d[(a, b)] = float(dist)
    bad = 0
    for i, j in pairs:
        got = float(tri[off(i, j, n)])
        exp = d.get((ids[i], ids[j])) or d.get((ids[j], ids[i]))
        if exp is None:
            print(f"[verify] pair ({ids[i]},{ids[j]}) missing from mash output")
            bad += 1
        elif abs(got - exp) > 1e-6:
            print(f"[verify] MISMATCH {ids[i]} {ids[j]} got={got} exp={exp}")
            bad += 1
    print(f"[verify] spot-checked {len(pairs)} pairs vs direct mash: "
          f"mismatches={bad}")

    # symmetry: recompute a few (j,i) from the chunk sketches is O(1) here —
    # the triangle stores a<b only, so nothing to check beyond offsets.
    summary = {
        "n": n,
        "tri_len": tri_len,
        "size_bytes": size,
        "size_ok": ok_size,
        "row0_min": float(row0.min()),
        "row0_max": float(row0.max()),
        "zero_distance_pairs": int(nz),
        "sample_mean": float(sample.mean()),
        "sample_median": float(np.median(sample)),
        "sample_frac_one": float((sample == 1.0).mean()),
        "spot_check_pairs": len(pairs),
        "spot_check_mismatches": bad,
    }
    with open(os.path.join(base, "triangle_verify.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    sys.exit(0 if (ok_size and bad == 0) else 1)

if __name__ == "__main__":
    main()
