#!/usr/bin/env python3
"""build_tree.py — UPGMA tree from a float32 triangle via scipy average linkage.

Reads:  ids.txt (n ids)
        <dist>  (little-endian float32 condensed upper triangle, n*(n-1)/2)
Writes: <out>.nwk   (Newick, leaves named by id)
        <out>_stats.json
"""
import json, os, sys, time, argparse
import numpy as np
from scipy.cluster.hierarchy import linkage

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--dist", required=True)
    ap.add_argument("--out", required=True, help="output nwk path (no ext)")
    args = ap.parse_args()

    t0 = time.time()
    with open(args.ids) as f:
        ids = [l.strip() for l in f if l.strip()]
    n = len(ids)
    print(f"[load] n={n} ids in {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    tri32 = np.memmap(args.dist, dtype="<f4", mode="r", shape=(n * (n - 1) // 2,))
    print(f"[load] mmap triangle {tri32.nbytes/1e9:.2f} GB in {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    condensed = np.array(tri32, dtype=np.float64, copy=True)
    del tri32
    print(f"[load] condensed float64 {condensed.nbytes/1e9:.2f} GB in {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    print("[tree] running scipy average linkage (UPGMA)...", flush=True)
    Z = linkage(condensed, method="average")
    del condensed
    print(f"[tree] linkage done in {time.time()-t0:.1f}s, Z shape={Z.shape}", flush=True)

    # ---- Newick emission (iterative) ----
    t0 = time.time()
    left = np.zeros(2 * n - 1, dtype=np.int64)
    right = np.zeros(2 * n - 1, dtype=np.int64)
    for k, row in enumerate(Z):
        i, j = int(row[0]), int(row[1])
        node = n + k
        left[node], right[node] = i, j
    root = n + (n - 2)

    parts = [None] * (2 * n - 1)
    stack = [(root, False)]
    while stack:
        node, visited = stack.pop()
        if visited:
            l, r = left[node], right[node]
            h = 0.0 if node < n else float(Z[node - n, 2])
            hl = 0.0 if l < n else float(Z[l - n, 2])
            hr = 0.0 if r < n else float(Z[r - n, 2])
            bl = f"{(h - hl) / 2:.6f}" if node >= n else "0"
            br = f"{(h - hr) / 2:.6f}" if node >= n else "0"
            sl = parts[l] if l >= n else f"{ids[l]}"
            sr = parts[r] if r >= n else f"{ids[r]}"
            parts[node] = f"({sl}:{bl},{sr}:{br})"
        else:
            stack.append((node, True))
            l, r = left[node], right[node]
            if l >= n:
                stack.append((l, False))
            if r >= n:
                stack.append((r, False))
    newick = parts[root] + ";"
    with open(args.out + ".nwk", "w") as f:
        f.write(newick + "\n")
    print(f"[tree] wrote {args.out}.nwk ({len(newick)/1e6:.2f} MB) in {time.time()-t0:.1f}s", flush=True)

    # ---- stats ----
    heights = np.zeros(2 * n - 1)
    for k, row in enumerate(Z):
        heights[n + k] = row[2]
    max_height = float(heights[-1]) if len(heights) else 0.0
    child_dists = []
    for k in range(n - 1):
        node = n + k
        i, j = int(Z[k, 0]), int(Z[k, 1])
        h = float(Z[k, 2])
        hi = 0.0 if i < n else float(Z[i - n, 2])
        hj = 0.0 if j < n else float(Z[j - n, 2])
        child_dists.append(h - hi)
        child_dists.append(h - hj)
    child_dists = np.array(child_dists)
    stats = {
        "n_leaves": n,
        "n_internal": n - 1,
        "method": "UPGMA (scipy average linkage)",
        "tree_height": max_height,
        "branch_length_mean": float(child_dists.mean()),
        "branch_length_median": float(np.median(child_dists)),
        "branch_length_min": float(child_dists.min()),
        "branch_length_max": float(child_dists.max()),
        "branch_length_zero_count": int((child_dists < 1e-9).sum()),
        "newick_bytes": len(newick),
        "elapsed_s": 0.0,
    }
    with open(args.out + "_stats.json", "w") as f:
        json.dump(stats, f, indent=2)
    print(f"[tree] stats: {json.dumps(stats, indent=2)}", flush=True)

if __name__ == "__main__":
    main()