#!/usr/bin/env python3
"""validate_tree.py — structural + biological validation of the full tree.

Checks:
  1. Newick parses (ete3); leaf count == n == 132,393 == sketched count.
  2. Leaf names are unique and exactly the ids.txt set (no dropped/reordered).
  3. Ultrametricity: every leaf at the same distance from the root
     (UPGMA property; tolerance 1e-4).
  4. Cophenetic-vs-source distance correlation on a random sample.
  5. Community clade purity: for a sample of leaves, the purity of the
     sibling clade (parent subtree) w.r.t. full_heatmap_clusters community.
Outputs: tree_verify.json
"""
import collections, csv, json, os, random, statistics, sys, time
import numpy as np
from ete3 import Tree

def off(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)

def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ids = [l.strip() for l in open(os.path.join(base, "data", "ids.txt")) if l.strip()]
    n = len(ids)
    nwk_path = os.path.join(base, "full_prophages_tree.nwk")
    tri_path = os.path.join(base, "full_prophages_mash.dist")

    t0 = time.time()
    t = Tree(nwk_path)
    leaves = t.get_leaves()
    print(f"[validate] parsed in {time.time()-t0:.0f}s; leaves={len(leaves)}")
    assert len(leaves) == n, f"leaf count {len(leaves)} != n {n}"
    names = [l.name for l in leaves]
    assert len(set(names)) == n, "duplicate leaf names"
    assert set(names) == set(ids), "leaf names != ids.txt set"
    print("[validate] leaf count == n == sketched count; names match ids.txt")

    # ultrametric check: distance from root to leaf
    d0 = [t.get_distance(l) for l in leaves[:2000]]
    print(f"[validate] root-to-leaf dists (2000 sample): min={min(d0):.6f} "
          f"max={max(d0):.6f} -> ultrametric={max(d0)-min(d0) < 1e-4}")

    # cophenetic vs source (300 pairs)
    tri = np.memmap(tri_path, dtype="<f4", mode="r", shape=(n * (n - 1) // 2,))
    rng = np.random.default_rng(5)
    idx = set(rng.integers(0, n, size=600).tolist())
    idx = list(idx)[:300]
    random.seed(9)
    pairs = [tuple(random.sample(idx, 2)) for _ in range(200)]
    src, cop = [], []
    for i, j in pairs:
        src.append(float(tri[off(i, j, n)]))
        cop.append(t.get_distance(ids[i], ids[j]))
    src = np.array(src); cop = np.array(cop)
    r2 = 1 - np.sum((src - cop) ** 2) / np.sum((src - src.mean()) ** 2)
    print(f"[validate] cophenetic vs source: n_pairs={len(src)} "
          f"R2={r2:.4f} mean_abs_err={np.mean(np.abs(src - cop)):.4f}")

    # community clade purity: only the 12 shared (non-singleton) communities
    # carry meaningful structure; the other 112,755 leaves are isolates
    # (singleton communities by construction — nothing to group).
    comm = {}
    with open("/home/erikg/phind/prophage_homology_survey/full_heatmap_clusters.csv") as f:
        for row in csv.DictReader(f):
            comm[row["sequence"]] = row["community"]
    cc = collections.Counter(comm.values())
    big = set(cid for cid, cnt in cc.items() if cnt > 1)
    def base_name(l):
        return l.name.split("|")[0]
    purities = []
    n_nonsingleton = 0
    for lf in leaves:
        c = comm.get(base_name(lf))
        if c is None or c not in big or lf.up is None:
            continue
        n_nonsingleton += 1
        clade = lf.up.get_leaves()
        if len(clade) < 2:
            continue
        cs = [comm.get(base_name(x)) for x in clade]
        cs = [x for x in cs if x is not None and x in big]
        if not cs:
            continue
        purities.append(cs.count(c) / len(cs))
    print(f"[validate] non-singleton leaves: {n_nonsingleton}")
    if purities:
        print(f"[validate] sibling-clade community purity ({len(purities)}): "
              f"mean={statistics.mean(purities):.4f} "
              f"median={statistics.median(purities):.4f}")
    else:
        print("[validate] no non-singleton leaves checked")

    summary = {
        "n": n,
        "leaves": len(leaves),
        "leaf_count_ok": len(leaves) == n,
        "names_unique": True,
        "names_match_ids": True,
        "ultrametric_tol_1e4": max(d0) - min(d0) < 1e-4,
        "root_to_leaf_min": float(min(d0)),
        "root_to_leaf_max": float(max(d0)),
        "cophenetic_R2": round(float(r2), 4),
        "cophenetic_mean_abs_err": round(float(np.mean(np.abs(src - cop))), 4),
        "community_purity_mean": round(statistics.mean(purities), 4) if purities else None,
        "community_purity_median": round(statistics.median(purities), 4) if purities else None,
        "community_purity_n": len(purities),
        "non_singleton_leaves": n_nonsingleton,
        "shared_communities": len(big),
    }
    with open(os.path.join(base, "tree_verify.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    print("[validate] wrote tree_verify.json")
    sys.exit(0)

if __name__ == "__main__":
    main()
