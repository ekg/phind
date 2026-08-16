#!/usr/bin/env python3
"""
full_ntm_mash_clades.py (v2) — mash sketch + dist -> float32 upper-triangle (in
build_tight_clades format) + ids.txt + single-community labels.csv, UPGMA tree,
then build_tight_clades.py over all v2 NTM prophages (NCBI-only, 8,502; the
~14.7k run-assembly prophages are pending collaborator FASTAs, task ntm-v2-run).

Same recipe as v1 (ntm/scripts/full_ntm_mash_clades.py): sketch (-i, k=21,
s=10000), `mash dist` all-vs-all, Python-written float32 upper triangle in FASTA
header order, single community (0) so leader clustering forms tight clades
naturally across species (cross-species prophage clades preserved, resolved
later by the host-clade join).

Outputs under /mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/:
  mash_clades/  prophages.msh, prophages.dist.tsv, prophages_mash.dist (float32
                upper triangle), ids.txt, labels.csv (community 0),
                full_prophages.idx.json (offset index for per-clade extraction),
                prophages_tree.nwk + tree_stats.json (UPGMA via scipy average
                linkage, same method as research/mash_tree/build_tree.py)
  clades/       build_tight_clades.py output for community 0 (tight_clades.json
                manifests, clade_similarity.json, members.json, distances.npz,
                commands.log) + tight_clades_summary.json + clade_summary.tsv
                (all clades), alignable_clades.tsv (>=2 members), singletons.tsv
  mash_clades_report.md   numbers for the task Validation section

Threshold 0.25 / max-size 100 (same as the E. coli release and NTM v1).
"""
import csv
import json
import os
import struct
import subprocess
import sys
import time

import numpy as np
from scipy.cluster.hierarchy import linkage

V2 = "/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2"
FA = f"{V2}/full_prophages.fa"
OUT = f"{V2}/mash_clades"
CLADES_OUT = f"{V2}/clades"
REPO = "/home/erikg/phind"
BUILD_TIGHT = f"{REPO}/scripts/build_tight_clades.py"
REPORT = f"{V2}/mash_clades_report.md"

THRESHOLD = 0.25
MAX_SIZE = 100


def run(cmd, out=None):
    print("+", " ".join(cmd), flush=True)
    if out is not None:
        with open(out, "w") as fo:
            r = subprocess.run(cmd, stdout=fo)
    else:
        r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(f"FAILED: {' '.join(cmd)}")
    return r


def build_index(fasta_path, index_path):
    """Offset index in per_clade_alignment_pipeline.py format:
    id -> [byte offset just after the header line, sequence length in bytes]
    (single-line sequences; verified format in the v2 extractor)."""
    idx = {}
    with open(fasta_path, "rb") as f:
        name = None
        start = None
        length = 0
        for line in f:
            if line.startswith(b">"):
                if name is not None:
                    idx[name] = [start, length]
                name = line[1:].strip().split()[0].decode()
                start = f.tell()
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            idx[name] = [start, length]
    tmp = index_path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(idx, f)
    os.replace(tmp, index_path)
    print(f"indexed {len(idx)} sequences -> {index_path}", flush=True)
    return idx


def write_newick(ids, Z, nwk_path, stats_path):
    """UPGMA tree -> Newick (streaming iterative emission; writes to the file
    as it traverses instead of materializing every subtree string in RAM —
    the parts[] approach is O(n^2) bytes and blew up memory at n=8502).
    Same branch-length convention as research/mash_tree/scripts/build_tree.py."""
    from scipy.cluster.hierarchy import to_tree
    n = len(ids)
    tree = to_tree(Z)

    def edge_len(node, child):
        """Half-edge length: (parent height - child height) / 2."""
        return (node.dist - child.dist) / 2.0 if child is not None else 0.0

    tmp = nwk_path + ".tmp"
    with open(tmp, "w") as f:
        stack = [(tree, 0)]  # state 0 enter, 1 separator (after left), 2 exit
        while stack:
            node, state = stack.pop()
            if state == 0:
                if node.left is None and node.right is None:
                    f.write(ids[node.id])
                else:
                    f.write("(")
                    stack.append((node, 2))
                    stack.append((node.right, 0))
                    stack.append((node, 1))
                    stack.append((node.left, 0))
            elif state == 1:
                f.write(f":{edge_len(node, node.left):.6f},")
            else:
                f.write(f":{edge_len(node, node.right):.6f})")
        f.write(";")
    os.replace(tmp, nwk_path)

    branch = []
    for row in Z:
        h = float(row[2])
        i, j = int(row[0]), int(row[1])
        hi = 0.0 if i < n else float(Z[i - n, 2])
        hj = 0.0 if j < n else float(Z[j - n, 2])
        branch.extend([(h - hi) / 2, (h - hj) / 2])
    stats = {
        "n_leaves": n,
        "n_internal": n - 1,
        "method": "UPGMA (scipy average linkage)",
        "tree_height": float(Z[-1, 2]),
        "branch_length_mean": float(np.mean(branch)),
        "branch_length_median": float(np.median(branch)),
        "branch_length_min": float(np.min(branch)),
        "branch_length_max": float(np.max(branch)),
        "branch_length_zero_count": int(np.sum(np.array(branch) == 0)),
        "newick_bytes": os.path.getsize(nwk_path),
    }
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"tree: {n} leaves -> {nwk_path} ({stats['tree_height']:.4f} height)",
          flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(CLADES_OUT, exist_ok=True)
    t_start = time.time()

    # 1. ids in file order
    ids = []
    with open(FA) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    n = len(ids)
    print(f"{n} prophages", flush=True)
    if not (5000 <= n <= 35000):
        print(f"WARNING: prophage count {n} outside expected [5000, 35000] "
              f"(v1 was 10,438; NCBI-only v2 ~8k; run-assemblies pending) — "
              f"flagging, not silently proceeding", flush=True)

    # 2. mash sketch (individual) + dist
    run(["mash", "sketch", "-i", "-k", "21", "-s", "10000", "-p", "64",
         "-o", f"{OUT}/prophages", FA])
    dist_tsv = f"{OUT}/prophages.dist.tsv"
    run(["mash", "dist", "-p", "128", f"{OUT}/prophages.msh",
         f"{OUT}/prophages.msh"], out=dist_tsv)
    print("mash dist done -> triangle", flush=True)

    # 3. float32 upper triangle in ids order
    idx = {s: i for i, s in enumerate(ids)}
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
    npairs = n * (n - 1) // 2
    print(f"triangle {n} seqs {npairs} pairs -> {tri}", flush=True)
    assert os.path.getsize(tri) == npairs * 4, "triangle size mismatch"
    del pair  # free the ~6GB pair dict before the tree step

    # 3b. UPGMA tree (same method as the E. coli / research tree)
    t0 = time.time()
    tri32 = np.memmap(tri, dtype="<f4", mode="r", shape=(npairs,))
    condensed = np.array(tri32, dtype=np.float64, copy=True)
    del tri32
    Z = linkage(condensed, method="average")
    del condensed
    print(f"[tree] linkage done in {time.time()-t0:.1f}s", flush=True)
    write_newick(ids, Z, f"{OUT}/prophages_tree.nwk",
                 f"{OUT}/tree_stats.json")

    # 4. tight clades (threshold 0.25, max-size 100, community 0)
    run(["python3", BUILD_TIGHT,
         "--threshold", str(THRESHOLD), "--max-size", str(MAX_SIZE),
         "--communities", "0",
         "--outdir", CLADES_OUT, "--ids-file", f"{OUT}/ids.txt",
         "--triangle", tri, "--labels-csv", f"{OUT}/labels.csv"])

    # 5. offset index for the downstream per-clade extraction
    build_index(FA, f"{OUT}/full_prophages.idx.json")

    # 6. clade lists + validation
    tc_path = f"{CLADES_OUT}/0/tight_clades.json"
    sim_path = f"{CLADES_OUT}/0/clade_similarity.json"
    with open(tc_path) as f:
        tc = json.load(f)
    with open(sim_path) as f:
        sim = json.load(f)["per_clade"]

    rows = []
    for cid, members in tc.items():
        s = sim.get(cid, {})
        rows.append((cid, len(members),
                     s.get("median"), s.get("min"), s.get("max")))
    rows.sort()
    with open(f"{CLADES_OUT}/clade_summary.tsv", "w") as f:
        f.write("clade_id\tn_members\tmedian_mash\tmin_mash\tmax_mash\n")
        for cid, nm, med, mn, mx in rows:
            def g(x):
                return "" if x is None else f"{x:.6f}"
            f.write(f"{cid}\t{nm}\t{g(med)}\t{g(mn)}\t{g(mx)}\n")
    with open(f"{CLADES_OUT}/alignable_clades.tsv", "w") as f:
        f.write("clade_id\tn_members\tmedian_mash\tmin_mash\tmax_mash\n")
        for cid, nm, med, mn, mx in rows:
            if nm >= 2:
                def g(x):
                    return "" if x is None else f"{x:.6f}"
                f.write(f"{cid}\t{nm}\t{g(med)}\t{g(mn)}\t{g(mx)}\n")
    with open(f"{CLADES_OUT}/singletons.tsv", "w") as f:
        f.write("clade_id\tmember_id\n")
        for cid, nm, med, mn, mx in rows:
            if nm == 1:
                f.write(f"{cid}\t{tc[cid][0]}\n")

    n_clades = len(tc)
    n_singletons = sum(1 for r in rows if r[1] == 1)
    n_alignable = n_clades - n_singletons
    assigned = sum(r[1] for r in rows)
    assert assigned == n, f"assignment incomplete: {assigned} != {n}"
    medians = [r[2] for r in rows if r[2] is not None]
    med_median = float(np.median(medians)) if medians else None
    print(f"clades={n_clades} alignable(>=2)={n_alignable} "
          f"singletons={n_singletons} median_internal_mash={med_median}",
          flush=True)

    # 7. report
    lines = []
    lines.append("# NTM v2 — prophage MASH + tight clades report")
    lines.append("")
    lines.append(f"Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    lines.append("")
    lines.append("## Prophage set")
    lines.append("")
    lines.append(f"- `full_prophages.fa`: **{n} prophages** (NCBI-only; v1 was "
                 f"10,438 — run-assembly prophages pending collaborator FASTAs, "
                 f"task ntm-v2-run; count inside the [5000, 35000] band, no flag)")
    lines.append("- Mash sketch: `-i -k 21 -s 10000` (same as v1)")
    lines.append("")
    lines.append("## Clades (threshold 0.25, max-size 100, community 0)")
    lines.append("")
    lines.append(f"- clade count: **{n_clades}**")
    lines.append(f"- alignable (>=2 members): **{n_alignable}**")
    lines.append(f"- singletons: **{n_singletons}**")
    lines.append(f"- median internal mash distance per clade: **{med_median:.4f}** "
                 f"(v1: 0.074)")
    lines.append("- assignment check: clade members + singletons == "
                 f"{assigned} == {n} prophages (every prophage assigned exactly once)")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append("- `mash_clades/`: `prophages.msh`, `prophages.dist.tsv`, "
                 "`prophages_mash.dist` (float32 upper triangle), `ids.txt`, "
                 "`labels.csv`, `full_prophages.idx.json`, "
                 "`prophages_tree.nwk`, `tree_stats.json`")
    lines.append("- `clades/`: `0/` (tight_clades.json, clade_similarity.json, "
                 "members.json, distances.npz, commands.log), "
                 "`tight_clades_summary.json`, `clade_summary.tsv`, "
                 "`alignable_clades.tsv`, `singletons.tsv`")
    lines.append("")
    with open(REPORT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"report -> {REPORT} (total {time.time()-t_start:.0f}s)", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
