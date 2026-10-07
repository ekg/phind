#!/usr/bin/env python3
"""v3.1 chimerism metric — is each reconstruction a single member's path, or a mosaic?

For the chosen ML path P = [p1..pn] (from ml.traversal.json best_sample.path), and
each member's ordered partition list (from partitions.bed), compute:

  longest_run          max k such that some path window of length k appears
                       CONTIGUOUSLY in a single member's partition order
  longest_run_frac     longest_run / n
  n_members_to_cover   greedy minimum number of members whose contiguous runs
                       can tile the whole path
  verdict              "single_member" if longest_run == n, else "chimeric"

This is the synthesis-relevant check: a junction not present in any single genome
is a novel sequence, not an observed one. (`path_adj_fraction=1.0` only proves each
step is an observed adjacency *somewhere*, not that the path came from one genome.)
"""
from __future__ import annotations
import csv, json, os, sys
from collections import defaultdict

NV = "/mnt/nvme3n1/erikg/phind-genome-work"


def member_partitions(bed):
    """member -> list of partition ids ordered by start; and partition interval bp."""
    per = defaultdict(list)
    bp = {}
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            seq, st, en, pid = f[0], int(f[1]), int(f[2]), int(f[3])
            per[seq].append((st, pid))
            bp[pid] = max(bp.get(pid, 0), en - st)
    return {k: [p for _, p in sorted(v)] for k, v in per.items()}, bp


def longest_run(path, per):
    pos = {p: i for i, p in enumerate(path)}
    best = 0
    for m, plist in per.items():
        idxs = [pos[p] for p in plist if p in pos]
        run = 1
        for a, b in zip(idxs, idxs[1:]):
            if b == a + 1:
                run += 1
                best = max(best, run)
            else:
                run = 1
            if best == len(path):
                return best
        if idxs:
            best = max(best, 1)
    return best


def greedy_cover(path, per):
    """minimum number of members' contiguous windows covering the whole path."""
    pos = {p: i for i, p in enumerate(path)}
    windows = []
    for m, plist in per.items():
        idxs = [pos[p] for p in plist if p in pos]
        s = 0
        for e in range(1, len(idxs) + 1):
            if e == len(idxs) or idxs[e] != idxs[e - 1] + 1:
                if e - s >= 1:
                    windows.append((idxs[s], idxs[e - 1]))
                s = e
    windows.sort(key=lambda w: (-(w[1] - w[0]), w[0]))
    covered = [False] * len(path)
    used = 0
    for a, b in windows:
        if any(not covered[i] for i in range(a, b + 1)):
            for i in range(a, b + 1):
                covered[i] = True
            used += 1
        if all(covered):
            break
    return used


def main():
    ids = [l.strip() for l in open(sys.argv[1])] if len(sys.argv) > 1 else \
        sorted(d for d in os.listdir(f"{NV}/ntm/v3.1/ml") if d.startswith("0_"))
    rows = []
    for cid in ids:
        tj = f"{NV}/ntm/v3.1/ml/{cid}/ml.traversal.json"
        bed = f"{NV}/ntm/v3/clades/{cid}/partitions.bed"
        if not (os.path.exists(tj) and os.path.exists(bed)):
            continue
        d = json.load(open(tj))
        path = (d.get("best_sample") or {}).get("path")
        if not path:
            continue
        per, bp = member_partitions(bed)
        lr = longest_run(path, per)
        cov = greedy_cover(path, per)
        rows.append({"clade_id": cid, "n_path_steps": len(path), "n_members": len(per),
                     "longest_run": lr, "longest_run_frac": round(lr / len(path), 4),
                     "n_members_to_cover": cov,
                     "verdict": "single_member" if lr == len(path) else "chimeric"})
    out = f"{NV}/ntm/v3.1/chimerism.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    single = sum(1 for r in rows if r["verdict"] == "single_member")
    import statistics
    fr = [r["longest_run_frac"] for r in rows]
    cov = [r["n_members_to_cover"] for r in rows]
    print(f"clades scored: {n}")
    print(f"  SINGLE-MEMBER (path == one genome's order): {single} ({100*single/n:.1f}%)")
    print(f"  longest_run_frac: median={statistics.median(fr):.3f}  p25={sorted(fr)[n//4]:.3f}  p75={sorted(fr)[3*n//4]:.3f}")
    print(f"  n_members_to_cover: median={statistics.median(cov):.0f}  p90={sorted(cov)[int(n*0.9)]}  max={max(cov)}")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
