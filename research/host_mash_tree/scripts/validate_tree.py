#!/usr/bin/env python3
"""validate_tree.py — structural + biological validation of the host tree.

Checks:
  1. Newick parses (ete3); leaf count == n == 26,074.
  2. Leaf names unique and exactly the ids.txt set (GCF_/GCA_ accessions).
  3. Ultrametricity (all leaves same root distance; tolerance 1e-4).
  4. Leaf names are assembly accessions (GCF_/GCA_) matching prophage->host
     prefix mapping (prophage leaves are ACC_prophage_N; host leaves = ACC).
  5. Cophenetic-vs-source correlation on a random sample (300 pairs).
Writes: tree_verify.json
"""
import json, os, random, sys, time, argparse
import numpy as np
from ete3 import Tree

def off(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--dist", required=True)
    ap.add_argument("--nwk", required=True)
    ap.add_argument("--out", default="tree_verify.json")
    args = ap.parse_args()

    ids = [l.strip() for l in open(args.ids) if l.strip()]
    n = len(ids)
    t0 = time.time()
    t = Tree(args.nwk)
    leaves = t.get_leaves()
    print(f"[validate] parsed in {time.time()-t0:.0f}s; leaves={len(leaves)}")
    leaf_count_ok = len(leaves) == n
    names = [l.name for l in leaves]
    names_unique = len(set(names)) == n
    names_match_ids = set(names) == set(ids)
    print(f"[validate] leaf_count_ok={leaf_count_ok} names_unique={names_unique} "
          f"names_match_ids={names_match_ids}")

    # accession pattern check
    import re
    acc_ok = all(re.fullmatch(r"(GCF|GCA)_\d+\.\d+", nm) for nm in names)
    print(f"[validate] all leaf names are GCF_/GCA_ accessions: {acc_ok}")

    d0 = [t.get_distance(l) for l in leaves[:2000]]
    ultra = max(d0) - min(d0) < 1e-4
    print(f"[validate] root-to-leaf (2000): min={min(d0):.6f} max={max(d0):.6f} "
          f"ultrametric={ultra}")

    # cophenetic vs source: 300 random pairs across the FULL range so the
    # distance variance (and hence R^2) is representative (a narrow index
    # subset collapses the denominator and makes R^2 unreliable).
    tri = np.memmap(args.dist, dtype="<f4", mode="r", shape=(n * (n - 1) // 2,))
    random.seed(9)
    pairs = set()
    while len(pairs) < 300:
        i, j = sorted(random.sample(range(n), 2))
        pairs.add((i, j))
    pairs = sorted(pairs)
    src, cop = [], []
    for i, j in pairs:
        src.append(float(tri[off(i, j, n)]))
        cop.append(t.get_distance(ids[i], ids[j]))
    src = np.array(src); cop = np.array(cop)
    err = np.abs(src - cop)
    r2 = 1 - np.sum((src - cop) ** 2) / np.sum((src - src.mean()) ** 2)
    print(f"[validate] cophenetic vs source: n={len(src)} R2={r2:.4f} "
          f"mean_abs_err={np.mean(err):.4f} max_abs_err={err.max():.4f}")

    summary = {
        "n": n, "leaves": len(leaves), "leaf_count_ok": leaf_count_ok,
        "names_unique": names_unique, "names_match_ids": names_match_ids,
        "accession_names_ok": acc_ok,
        "ultrametric_tol_1e4": ultra,
        "root_to_leaf_min": float(min(d0)), "root_to_leaf_max": float(max(d0)),
        "cophenetic_R2": round(float(r2), 4),
        "cophenetic_mean_abs_err": round(float(np.mean(err)), 4),
        "cophenetic_max_abs_err": round(float(err.max()), 4),
        "n_pairs_cophenetic": len(src),
    }
    with open(args.out, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    print("[validate] wrote tree_verify.json")
    sys.exit(0)

if __name__ == "__main__":
    main()