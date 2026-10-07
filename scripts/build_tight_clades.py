#!/usr/bin/env python3
"""
build_tight_clades.py — MASH-distance tight clades from the 12 old communities.

Context (measured problem, see task per-clade-fastga): the old 12-community
labels are internally incoherent — e.g. c2/c4 have ~55-72% of within-community
pairs at MASH distance 1.0. Aligning all pairs within those communities wastes
compute on unrelated sequences. This script sub-clusters each community by MASH
distance so that EVERY member is within `--threshold` (default 0.10, ~90% ANI)
of its clade representative, caps clade size at `--max-size` (default 800,
splitting further if a clade fills up), and leaves sequences too divergent to
join any clade as singleton clades.

Distances are read from the precomputed all-prophage triangle
research/mash_tree/full_prophages_mash.dist (float32 upper triangle, rows in
research/mash_tree/data/ids.txt order; offset(a,b) = a*(2n-a-1)/2 + (b-a-1)).

Clustering algorithm (greedy leader clustering):
  1. Build the within-community pairwise distance matrix.
  2. Order members by decreasing "degree" = number of members within the
     similarity threshold (most-connected first); they become leaders in that
     order, which favors compact clades anchored at dense centers.
  3. For each member in that order: if it is within threshold of an existing
     leader whose clade is not full, assign it to the NEAREST such leader;
     otherwise it becomes a new leader. Every member ends up within threshold
     of its leader by construction. Clades that would exceed max-size split
     (the member becomes a new leader), so no clade exceeds max-size.

Split modes (`--split-mode`):
  cap (default): the v3 behaviour described above — leader clustering capped
    at `--max-size`, then `tighten_clades` (median-criterion re-anchoring).
    Output files are unchanged by this option when it is left at the default.
  medoids: leader clustering at `--threshold` with NO size cap, then any clade
    whose members are further than `--medoid-radius` (default 0.10) from its
    medoid is recursively split with k-medoids (k=2) until every member is
    within the radius of its clade medoid. Split points are structure-driven
    rather than count-driven. `--max-size` is ignored in this mode.
    Clade ids keep the v3 `<community>_<NNNN>` shape for initial clades and
    add `_s1`/`_s2` suffixes for the children produced by each split, so the
    ids in `split_trace.json` resolve against `tight_clades.json`.
    A missing (NaN) MASH distance cannot be certified within the radius, so
    a clade containing one is split as well — consistent with leader
    clustering, where `D <= threshold` is False for NaN. Its radius is
    reported as `null` in the trace.

Outputs per community C (into --outdir/C/):
  tight_clades.json      clade_id -> [member ids (FASTA headers), ...]
  split_trace.json       (medoids mode only) one record per k-medoids split:
                         parent id, size before/after, radius before/after,
                         iterations
  clade_similarity.json  per-clade internal pairwise MASH stats + global summary
  members.json           community member ids in triangle row order
  distances.npz          within-community float32 matrix + member order
                         (reusable, avoids re-reading the 35 GB triangle)
  commands.log           reproducible commands (appended)

Clade ids are stable: <community>_<NNNN> in creation order (representative
first, then members). The representative is listed first.

Usage:
  build_tight_clades.py [--threshold 0.10] [--max-size 800]
                        [--split-mode cap|medoids] [--medoid-radius 0.10]
                        [--communities 0,1,2,3,4,5,6,7,8,9,10,11]
                        [--outdir research/clades]
"""

import argparse
import csv
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np

REPO = "/home/erikg/phind"
TRIANGLE = os.path.join(REPO, "research/mash_tree/full_prophages_mash.dist")
IDS_FILE = os.path.join(REPO, "research/mash_tree/data/ids.txt")
LABELS_CSV = os.path.join(REPO, "research/mash_tree/full_prophages_labels.csv")


def load_ids(path=IDS_FILE):
    ids = [l.strip() for l in open(path)]
    return ids, {s: i for i, s in enumerate(ids)}


def load_communities(path=LABELS_CSV):
    comm = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f):
            c = row["community"]
            if c.isdigit() and int(c) < 12:
                comm[int(c)].append(row["sequence"])
    return comm


def triangle_offset(a, b, n):
    return a * (2 * n - a - 1) // 2 + (b - a - 1)


def read_community_matrix(ids, idx, members, mm):
    """Read the within-community pairwise distance matrix (float32, n x n).
    members must be sorted by triangle row index."""
    n = len(ids)
    m = len(members)
    D = np.full((m, m), np.nan, dtype=np.float32)
    np.fill_diagonal(D, 0.0)
    row_idx = [idx[s] for s in members]
    pos_of_row = {r: i for i, r in enumerate(row_idx)}   # triangle row -> local idx
    for i, a in enumerate(row_idx):
        if a + 1 >= n:
            continue
        row = mm[triangle_offset(a, a + 1, n): triangle_offset(a, n - 1, n) + 1]
        # columns in the triangle row are a+1 .. n-1; pick member columns > a
        cols = [b for b in row_idx[i + 1:] if b > a]
        if not cols:
            continue
        rel = np.array(cols, dtype=np.int64) - (a + 1)
        D[i, [pos_of_row[c] for c in cols]] = row[rel]
        D[[pos_of_row[c] for c in cols], i] = row[rel]
    return D


def leader_cluster(D, members, threshold, max_size):
    """Greedy leader clustering on the distance matrix D (m x m).
    Returns list of clades: each clade is a list of member ids (0-based
    indices into `members`), representative first."""
    m = len(members)
    close = D <= threshold
    # degree: number of members within threshold (excluding self)
    degree = close.sum(axis=1).astype(int)
    order = sorted(range(m), key=lambda i: (-degree[i], i))

    leaders = []          # index of leader member
    clades = []           # list of member indices per clade
    assigned = [False] * m

    for i in order:
        best = None
        best_dist = np.inf
        for c, lead in enumerate(leaders):
            if len(clades[c]) >= max_size:
                continue
            d = D[i, lead]
            if d <= threshold and d < best_dist:
                best = c
                best_dist = d
        if best is None:
            # new clade; i is the representative
            leaders.append(i)
            clades.append([i])
            assigned[i] = True
        else:
            clades[best].append(i)
            assigned[i] = True

    assert all(assigned)
    # verify invariant: every member within threshold of its leader
    for c in clades:
        lead = c[0]
        for j in c[1:]:
            d = D[j, lead]
            assert d <= threshold + 1e-6, (j, d, threshold)
    return clades, leaders


def tighten_clades(D, members, clades, threshold, max_size):
    """Split clades whose internal MEDIAN pairwise distance exceeds
    `threshold`, so every final clade (size>1) has median <= threshold,
    members within `threshold` of their leader, and size <= max_size.

    Iterative work queue: a clade with median > threshold is re-anchored via
    leader clustering at a tighter threshold (0.6x, floor 0.01). If that
    makes no progress (hub-like clade where all members hug one leader), the
    clade is split in half by distance to its representative and both halves
    are re-anchored at the full threshold. Every split strictly reduces the
    piece size, so the queue terminates; singletons trivially satisfy the
    median criterion."""
    from collections import deque

    def median_of(idx):
        if len(idx) < 2:
            return None
        sub = D[np.ix_(idx, idx)]
        vals = sub[np.triu_indices(len(idx), 1)]
        return float(np.median(vals))

    out = []
    queue = deque(clades)
    guard = 0
    while queue:
        guard += 1
        if guard > 1000000:
            raise RuntimeError("tighten_clades did not terminate")
        c = queue.popleft()
        if len(c) < 2:
            out.append(c)
            continue
        med = median_of(c)
        if med is not None and med <= threshold:
            out.append(c)
            continue
        # try tighter leader clustering on the subset
        t = max(threshold * 0.6, 0.01)
        sub_D = D[np.ix_(c, c)]
        local = np.arange(len(c))
        sub_local, _ = leader_cluster(sub_D, local, t, max_size)
        pieces = [[c[j] for j in p] for p in sub_local]
        if len(pieces) == 1 and set(pieces[0]) == set(c):
            # no progress: split by distance to representative
            rep = c[0]
            ds = D[c, rep]
            order = np.argsort(ds)
            half = max(1, len(c) // 2)
            pieces = [[c[j] for j in order[:half]],
                      [c[j] for j in order[half:]]]
            # re-anchor each half at the full threshold
            anchored = []
            for p in pieces:
                if len(p) < 2:
                    anchored.append(p)
                    continue
                pD = D[np.ix_(p, p)]
                plocal, _ = leader_cluster(pD, np.arange(len(p)), threshold, max_size)
                anchored.extend([[p[j] for j in q] for q in plocal])
            pieces = anchored
        queue.extend(pieces)
    return out


def medoid_of(D, idx):
    """Medoid of the members in `idx`: the member with the smallest total
    distance to the other members of `idx`. Only meaningful for clades with
    no missing distances; `cluster_radius` reports such clades as infinite
    and therefore split before a medoid is ever reported."""
    if len(idx) == 1:
        return idx[0]
    sub = D[np.ix_(idx, idx)]
    return idx[int(np.argmin(sub.sum(axis=1)))]


def cluster_radius(D, idx):
    """Max distance from the clade medoid to any of its members.
    0.0 for clades of size <= 1.

    Returns inf when `idx` contains a missing (NaN) pairwise distance. A
    missing distance cannot be certified to be within the radius, and the
    rest of the script already treats NaN as "not close" (leader clustering's
    `D <= threshold` is False for NaN, v3's median test fails on NaN), so such
    clades are split until they are NaN-free."""
    if len(idx) <= 1:
        return 0.0
    sub = D[np.ix_(idx, idx)]
    if np.isnan(sub).any():
        return float("inf")
    med = medoid_of(D, idx)
    return float(D[idx, med].max())


def kmedoids_split(D, idx, max_iter=100):
    """Split `idx` into two clades with k-medoids (k=2) on the pairwise
    distances in D. Seeds are the clade medoid and the member farthest from
    it; assignment (nearest medoid, ties to the lower index) and medoid
    update alternate until the assignment is stable.

    Returns (left, right, iterations) where both sides are non-empty member
    index lists and `iterations` is the number of assignment/update rounds.
    If the Lloyd loop degenerates to one non-empty side, the clade is split
    in half by distance to its medoid instead."""
    sub = D[np.ix_(idx, idx)]
    k = len(idx)
    c0 = int(np.argmin(sub.sum(axis=1)))
    c1 = int(np.argmax(sub[c0]))
    centers = [c0, c1]
    labels = None
    iters = 0
    for _ in range(max_iter):
        iters += 1
        new = np.argmin(sub[:, centers], axis=1)
        if labels is not None and np.array_equal(new, labels):
            labels = new
            break
        labels = new
        for side in (0, 1):
            memb = np.where(labels == side)[0]
            if len(memb):
                centers[side] = int(memb[np.argmin(
                    sub[np.ix_(memb, memb)].sum(axis=1))])
    left = np.where(labels == 0)[0]
    right = np.where(labels == 1)[0]
    if len(left) == 0 or len(right) == 0:
        order = np.argsort(sub[c0], kind="stable")
        half = max(1, k // 2)
        left = order[:half]
        right = order[half:]
    return ([idx[i] for i in left], [idx[i] for i in right], iters)


def split_until_radius(D, clades, names, radius, max_iter=100):
    """Recursively split every clade whose members exceed `radius` from its
    medoid (k-medoids, k=2 per split) until every clade satisfies the radius
    criterion. Each split strictly reduces size, so this terminates.

    Returns (out, trace):
      out   -- list of (clade_id, member index list), in input order with each
               split parent replaced by its two children (left then right)
      trace -- one record per split: parent id, sizes before/after, radii
               before/after (`null` radius = the clade contains a missing
               MASH distance) and the k-medoids iterations used."""
    out = []
    trace = []

    def radius_or_null(D_, idx_):
        r = cluster_radius(D_, idx_)
        return None if not np.isfinite(r) else round(r, 6)

    def visit(name, idx):
        idx = list(idx)
        r_before = cluster_radius(D, idx)
        if len(idx) <= 1 or r_before <= radius + 1e-9:
            out.append((name, idx))
            return
        left, right, iters = kmedoids_split(D, idx, max_iter=max_iter)
        if not left or not right:
            # cannot make progress; keep the clade rather than loop
            out.append((name, idx))
            return
        lname, rname = f"{name}_s1", f"{name}_s2"
        trace.append({
            "parent_clade_id": name,
            "size_before": len(idx),
            "sizes_after": [len(left), len(right)],
            "radius_before": radius_or_null(D, idx),
            "radii_after": [radius_or_null(D, left),
                            radius_or_null(D, right)],
            "iterations": iters,
        })
        visit(lname, left)
        visit(rname, right)

    for name, idx in zip(names, clades):
        visit(name, idx)
    return out, trace


def clade_stats(D, clades):
    """Per-clade internal similarity stats."""
    out = {}
    for ci, clade in enumerate(clades):
        if len(clade) == 1:
            out[ci] = {
                "n": 1, "pairs": 0, "median": None, "min": None, "max": None,
                "frac_le_threshold": None, "rep_index": 0,
            }
            continue
        idx = np.array(clade)
        sub = D[np.ix_(idx, idx)]
        vals = sub[np.triu_indices(len(clade), 1)]
        out[ci] = {
            "n": len(clade), "pairs": int(len(vals)),
            "median": float(np.median(vals)) if len(vals) else None,
            "min": float(vals.min()) if len(vals) else None,
            "max": float(vals.max()) if len(vals) else None,
            "frac_le_threshold": float((vals <= threshold).mean()) if len(vals) else None,
            "rep_index": 0,
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threshold", type=float, default=0.10,
                    help="max MASH distance to clade representative (~90%% ANI)")
    ap.add_argument("--max-size", type=int, default=800,
                    help="hard cap on clade size (split further if exceeded); "
                         "ignored by --split-mode medoids")
    ap.add_argument("--split-mode", choices=("cap", "medoids"), default="cap",
                    help="cap: leader clustering capped at --max-size "
                         "(v3 behaviour, default); medoids: uncapped leader "
                         "clustering then recursive k-medoids splits until "
                         "every member is within --medoid-radius of its medoid")
    ap.add_argument("--medoid-radius", type=float, default=0.10,
                    help="max MASH distance from a clade medoid allowed in "
                         "--split-mode medoids")
    ap.add_argument("--communities", default="0,1,2,3,4,5,6,7,8,9,10,11",
                    help="comma-separated community ids")
    ap.add_argument("--outdir", default=os.path.join(REPO, "research/clades"))
    ap.add_argument("--ids-file", default=IDS_FILE)
    ap.add_argument("--triangle", default=TRIANGLE)
    ap.add_argument("--labels-csv", default=LABELS_CSV)
    args = ap.parse_args()

    global threshold
    threshold = args.threshold

    communities = [int(x) for x in args.communities.split(",")]
    ids, idx = load_ids(args.ids_file)
    n = len(ids)
    comm = load_communities(args.labels_csv)

    print(f"triangle rows: {n}; communities: {communities}", flush=True)
    mm = np.memmap(args.triangle, dtype="<f4", mode="r",
                   shape=(n * (n - 1) // 2,))

    os.makedirs(args.outdir, exist_ok=True)
    summary_all = {}
    for c in communities:
        members = comm.get(c, [])
        if not members:
            print(f"community {c}: no members, skipping", flush=True)
            continue
        members = sorted(set(members))
        # sanity: all members in ids
        missing = [s for s in members if s not in idx]
        assert not missing, missing[:5]
        t0 = time.time()
        D = read_community_matrix(ids, idx, members, mm)
        t_read = time.time() - t0
        split_trace = None
        if args.split_mode == "cap":
            clades, leaders = leader_cluster(D, members, args.threshold, args.max_size)
            clades = tighten_clades(D, members, clades, args.threshold, args.max_size)
            clade_ids = [f"{c}_{ci:04d}" for ci in range(len(clades))]
        else:
            # re-anchor uncapped (max_size = number of members = no cap), then
            # split on structure until every member is within the medoid radius
            base, _ = leader_cluster(D, members, args.threshold, len(members))
            base_names = [f"{c}_{ci:04d}" for ci in range(len(base))]
            split, split_trace = split_until_radius(
                D, base, base_names, args.medoid_radius)
            clade_ids = [n for n, _ in split]
            clades = [m for _, m in split]
        stats = clade_stats(D, clades)

        cdir = os.path.join(args.outdir, str(c))
        os.makedirs(cdir, exist_ok=True)

        # tight_clades.json: clade_id -> member list (rep first)
        tc = {}
        for ci, clade in enumerate(clades):
            cid = clade_ids[ci]
            tc[cid] = [members[j] for j in clade]
        with open(os.path.join(cdir, "tight_clades.json"), "w") as f:
            json.dump(tc, f, indent=1)

        if split_trace is not None:
            with open(os.path.join(cdir, "split_trace.json"), "w") as f:
                json.dump({
                    "community": c,
                    "split_mode": args.split_mode,
                    "threshold": args.threshold,
                    "medoid_radius": args.medoid_radius,
                    "n_splits": len(split_trace),
                    "splits": split_trace,
                }, f, indent=1)

        sim = {}
        for ci, clade in enumerate(clades):
            cid = clade_ids[ci]
            sim[cid] = stats[ci]
        medians = [s["median"] for s in sim.values() if s["median"] is not None]
        report = {
            "community": c,
            "n_members": len(members),
            "threshold": args.threshold,
            "max_size": args.max_size,
            "n_clades": len(clades),
            "n_singletons": sum(1 for s in sim.values() if s["n"] == 1),
            "max_clade_size": max(s["n"] for s in sim.values()),
            "clade_size_median": float(np.median([s["n"] for s in sim.values()])),
            "clade_size_distribution": {
                "1": sum(1 for s in sim.values() if s["n"] == 1),
                "2-10": sum(1 for s in sim.values() if 2 <= s["n"] <= 10),
                "11-50": sum(1 for s in sim.values() if 11 <= s["n"] <= 50),
                "51-200": sum(1 for s in sim.values() if 51 <= s["n"] <= 200),
                "201-500": sum(1 for s in sim.values() if 201 <= s["n"] <= 500),
                "501-800": sum(1 for s in sim.values() if 501 <= s["n"] <= 800),
                ">800": sum(1 for s in sim.values() if s["n"] > 800),
            },
            "internal_similarity": {
                "median_median": float(np.median(medians)) if medians else None,
                "min_median": float(np.min(medians)) if medians else None,
                "max_median": float(np.max(medians)) if medians else None,
                "clades_with_median_le_threshold": sum(
                    1 for s in sim.values()
                    if s["median"] is not None and s["median"] <= args.threshold
                ),
                "clades_with_max_le_threshold": sum(
                    1 for s in sim.values()
                    if s["max"] is not None and s["max"] <= args.threshold
                ),
            },
            "time_read_triangle_s": round(t_read, 1),
        }
        with open(os.path.join(cdir, "clade_similarity.json"), "w") as f:
            json.dump({"report": report, "per_clade": sim}, f, indent=1)

        with open(os.path.join(cdir, "members.json"), "w") as f:
            json.dump(members, f, indent=0)

        np.savez_compressed(
            os.path.join(cdir, "distances.npz"),
            D=D.astype(np.float32), members=np.array(members, dtype=object),
            threshold=args.threshold,
        )

        with open(os.path.join(cdir, "commands.log"), "a") as f:
            if args.split_mode == "cap":
                f.write(
                    f"# build_tight_clades.py --threshold {args.threshold} "
                    f"--max-size {args.max_size} --communities {c}\n"
                )
            else:
                f.write(
                    f"# build_tight_clades.py --threshold {args.threshold} "
                    f"--split-mode {args.split_mode} "
                    f"--medoid-radius {args.medoid_radius} "
                    f"--communities {c}\n"
                )

        summary_all[c] = report
        print(json.dumps(report, indent=1), flush=True)

    with open(os.path.join(args.outdir, "tight_clades_summary.json"), "w") as f:
        json.dump(summary_all, f, indent=1)
    print("done")


if __name__ == "__main__":
    sys.exit(main())
