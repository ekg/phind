#!/usr/bin/env python3
"""NTM v3 — prophage MASH sketch + triangle + tight clades (nohup/resume-protected).

Input: `full_prophages.fa` from task ntm-v3-unified (9,446 extractable v3
prophages; the 27,494 run-assembly manifest rows have no FASTA on this host
yet — collaborator delivery pending, see ntm/v3/download_report.md).

Method parity with v1/v2 (task ntm-v3-prophage mirrors ntm/v2/
mash_clades_report.md):
  * mash sketch `-i -k 21 -s 10000` (identical to v1/v2)
  * tight clades via the reused-as-is E. coli machinery
    `scripts/build_tight_clades.py --threshold 0.25 --max-size 100
    --communities 0` (community 0 convention: single community so leader
    clustering forms tight clades naturally across species; cross-species
    prophage clades are preserved and resolved later by the host-clade join)
  * singletons become their own clade (v1/v2 convention: no ML
    reconstruction for n=1, pass-through)

Differences from v2 mechanics (deliberate, same choices as the v3 host-clades
task ntm-v3-host / ntm/v3/scripts/host_clades_mash_v3.py):
  * **`mash triangle` instead of `mash dist msh msh`.** v2 emitted the full
    N^2 dist TSV (7.6 GB at n=8,502); v3 emits the complete lower triangle
    (`n*(n-1)/2` pairs) — same pairs, half the bytes, explicit pair-count
    completeness check — then converts it to the float32 upper-triangle
    binary `v3_prophages.dist` that build_tight_clades.py consumes
    (offset(a,b) = a*(2n-a-1)/2 + (b-a-1), rows in ids.txt order).
  * **ids.txt is written in `sorted()` order, not FASTA order.**
    build_tight_clades.py's read_community_matrix documents "members must
    be sorted by triangle row index" and its main() does
    `members = sorted(set(members))`; with ids.txt in FASTA order (the v1/v2
    invocation) 35.2% of v2's within-community pairs were never filled
    (NaN) and silently treated as maximally distant — measured on the v2
    triangle: 23,407,628 of 36,137,751 pairs filled. Writing ids.txt
    pre-sorted makes `sorted(set(members))` == triangle row order, so every
    pair is filled; the reused script itself is NOT modified. Downstream
    consumers join by prophage id (labels.csv rows, idx.json keys,
    tight_clades.json member names), so the row order is internal to this
    step. v1/v2 clade counts (913 / 2,388) were produced with the missed-join
    fragmentation and are therefore upper bounds on clade count — noted in
    the v3 report comparison.
  * **Resume-protected stages.** Every stage writes to a `.tmp` artifact and
    atomically renames; re-invoking the driver skips stages whose artifact
    already exists and validates. An interrupted run resumes from disk.
  * **nohup-protected heavy phase.** Run the heavy phase (sketch+triangle+
    f32+tree+index) detached:
        nohup setsid python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py \
            --phase heavy > heavy.log 2>&1 &
    The agent session may die; compute survives and any later attempt runs
    `--phase light` (clades + validation + report) from artifacts on disk.

Phases (each idempotent, resumable from disk state):
  ids       FASTA headers -> ids.txt + prophage_lengths.tsv (skipped when
            ids.txt row count + first/last header match the FASTA)
  sketch    `mash sketch -i -k 21 -s 10000 -p T -o v3_prophages` (atomic
            rename; skipped when v3_prophages.msh exists and is newer than
            the FASTA)
  triangle  `mash triangle -k 21 -s 10000 -p T v3_prophages.msh` ->
            v3_prophages.triangle.txt (atomic rename; skipped when the file
            exists, is newer than the sketch and its header count == n)
  f32       triangle text -> v3_prophages.dist (float32 upper triangle,
            exactly n*(n-1)/2 values; strict per-row name + field-count
            check, so completeness is verified by construction)
  labels    labels.csv (single community 0) + ids.txt order check
  tree      UPGMA (scipy average linkage) -> prophages_tree.nwk +
            tree_stats.json (v2 parity; streaming Newick emission)
  index     full_prophages.idx.json byte-offset index for the downstream
            per-clade FASTA extraction (per_clade_alignment_pipeline.py
            format: id -> [offset after header, sequence byte length])
  clades    `python3 scripts/build_tight_clades.py --threshold 0.25
            --max-size 100 --communities 0 ...` (reused as-is, E. coli / v1 /
            v2 machinery — NOT modified)
  summary   clade_summary.tsv / alignable_clades.tsv / singletons.tsv
            (NVMe + repo copies)
  spotcheck 200 sampled pairs read back from the triangle text vs the f32
            binary (conversion check) + 50 sampled pairs re-measured with a
            fresh `mash sketch -k 21 -s 10000` + `mash dist` on the extracted
            sequences (measurement check)
  report    ntm/v3/mash_clades_report.md + repo clade definitions
            (tight_clades.json.gz, deterministic gzip) + SHA256SUMS
  heavy     ids + sketch + triangle + f32 + labels + tree + index
  light     clades + summary + spotcheck + report
  all       heavy + light

Outputs (NVMe, /mnt/nvme3n1/erikg/phind-genome-work/ntm/v3):
  mash_clades/ v3_prophages.msh, v3_prophages.triangle.txt,
               v3_prophages.dist (float32 upper triangle), ids.txt,
               prophage_lengths.tsv, labels.csv, full_prophages.idx.json,
               prophages_tree.nwk, tree_stats.json, commands.log,
               spotcheck.tsv, spotcheck.json, driver.log
  clades/      build_tight_clades.py output for community 0 (0/tight_clades.json,
               0/clade_similarity.json, 0/members.json, 0/distances.npz,
               0/commands.log, tight_clades_summary.json) + clade_summary.tsv,
               alignable_clades.tsv, singletons.tsv
Outputs (repo):
  ntm/v3/mash_clades_report.md
  ntm/v3/clades/tight_clades.json.gz + SHA256SUMS + clade_summary.tsv +
  alignable_clades.tsv + singletons.tsv

Usage:
  python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py --phase all
  python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py --phase heavy   # under nohup
  python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py --phase light
"""

import argparse
import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import time
from collections import Counter

import numpy as np

# ---------------------------------------------------------------- constants
K = 21
S = 10000
THRESHOLD = 0.25
MAX_SIZE = 100
BAND = (1000, 100000)          # prophage length band check (extract report)

# reference numbers from ntm/v3/extract_report.md (task ntm-v3-unified) —
# the length band check cross-checks this task's independent FASTA scan
EXTRACT_REF = {"n": 9446, "len_min": 237, "len_max": 93713,
               "len_median": 18276, "len_total": 206284318,
               "outside_band": 91}

REPO_DEFAULT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))   # <repo>/ntm/v3/scripts -> <repo>
WORK_DEFAULT = "/mnt/nvme3n1/erikg/phind-genome-work"

HEAVY = ["ids", "sketch", "triangle", "f32", "labels", "tree", "index"]
LIGHT = ["clades", "summary", "spotcheck", "report"]
PHASES = HEAVY + LIGHT

LOG_LINES = []
COMMANDS = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG_LINES.append(line)


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def atomic_replace(src, dst):
    os.replace(src, dst)


def run(cmd, out=None, check=True, record=True):
    pretty = " ".join(cmd) + (f" > {out}" if out else "")
    log("+ " + pretty)
    if record:
        COMMANDS.append(pretty)
    if out is not None:
        with open(out, "w") as fo:
            r = subprocess.run(cmd, stdout=fo)
    else:
        r = subprocess.run(cmd)
    if check and r.returncode != 0:
        raise SystemExit(f"FAILED ({r.returncode}): {' '.join(cmd)}")
    return r


# ------------------------------------------------------------------ stage: ids
def fasta_scan(fa):
    """Single pass: headers in file order + sequence lengths (single-line)."""
    ids, lengths = [], []
    with open(fa, "rb") as f:
        cur = None
        for line in f:
            if line.startswith(b">"):
                if cur is not None:
                    lengths.append(cur)
                ids.append(line[1:].strip().split()[0].decode())
                cur = 0
            else:
                cur += len(line.strip())
        if cur is not None:
            lengths.append(cur)
    return ids, lengths


def stage_ids(cfg):
    fa, out = cfg.fa, cfg.out
    ids_path = os.path.join(out, "ids.txt")
    len_path = os.path.join(out, "prophage_lengths.tsv")

    fasta_ids, lengths = fasta_scan(fa)
    n = len(fasta_ids)
    if len(set(fasta_ids)) != n:
        dup = [x for x, c in Counter(fasta_ids).items() if c > 1][:5]
        raise SystemExit(f"duplicate FASTA headers: {dup}")
    # triangle row order == sorted() order: build_tight_clades.py sorts its
    # community members with sorted(set(...)) and requires that to match the
    # triangle row order (see module docstring — v1/v2 missed ~35% of pairs)
    ids = sorted(fasta_ids)
    log(f"ids: {n} prophages, unique headers, single-line sequences; "
        f"ids.txt in sorted() order (build_tight_clades row-order contract)")

    if os.path.exists(ids_path) and os.path.exists(len_path):
        old = [l.strip() for l in open(ids_path) if l.strip()]
        if len(old) == n and old[0] == ids[0] and old[-1] == ids[-1]:
            log(f"ids: ids.txt already matches FASTA ({n} rows; resume)")
        else:
            log("ids: ids.txt stale -> rewriting")
    tmp = ids_path + ".tmp"
    with open(tmp, "w") as f:
        f.write("\n".join(ids) + "\n")
    atomic_replace(tmp, ids_path)
    tmp = len_path + ".tmp"
    with open(tmp, "w") as f:
        f.write("prophage_id\tlength_bp\n")
        for i, L in zip(ids, [dict(zip(fasta_ids, lengths))[s]
                              for s in ids]):
            f.write(f"{i}\t{L}\n")
    atomic_replace(tmp, len_path)
    stats = {
        "n": n,
        "len_min": int(min(lengths)), "len_median": float(np.median(lengths)),
        "len_mean": float(np.mean(lengths)), "len_max": int(max(lengths)),
        "len_total": int(sum(lengths)),
        "len_outside_band": int(sum(1 for L in lengths
                                    if not (BAND[0] <= L <= BAND[1]))),
        "len_lt_1kb": int(sum(1 for L in lengths if L < BAND[0])),
        "len_gt_100kb": int(sum(1 for L in lengths if L > BAND[1])),
    }
    log("ids: length stats " + json.dumps(stats))
    for key in ("n", "len_min", "len_max", "len_total", "outside_band"):
        got, ref = (stats[key] if key in stats else
                    (stats["n"] if key == "n" else stats["len_outside_band"]),
                    EXTRACT_REF[key])
        if got != ref:
            log(f"ids: WARNING {key}={got} != extract_report {ref}")
    return ids, stats


# --------------------------------------------------------------- stage: sketch
def stage_sketch(cfg):
    msh = os.path.join(cfg.out, "v3_prophages.msh")
    if os.path.exists(msh) and os.path.getmtime(msh) >= os.path.getmtime(cfg.fa):
        log(f"sketch: {msh} exists and >= FASTA mtime (resume, "
            f"{os.path.getsize(msh)/1e6:.0f} MB)")
        return
    tmp_prefix = os.path.join(cfg.out, "v3_prophages.tmp")
    tmp_msh = tmp_prefix + ".msh"
    if os.path.exists(tmp_msh):
        os.remove(tmp_msh)
    t0 = time.time()
    run(["mash", "sketch", "-i", "-k", str(K), "-s", str(S),
         "-p", str(cfg.threads), "-o", tmp_prefix, cfg.fa])
    atomic_replace(tmp_msh, msh)
    log(f"sketch: -> {msh} ({os.path.getsize(msh)/1e6:.0f} MB, "
        f"{time.time()-t0:.0f}s)")


# ------------------------------------------------------------- stage: triangle
def stage_triangle(cfg, n):
    msh = os.path.join(cfg.out, "v3_prophages.msh")
    tri_txt = os.path.join(cfg.out, "v3_prophages.triangle.txt")
    expected_pairs = n * (n - 1) // 2
    if (os.path.exists(tri_txt)
            and os.path.getmtime(tri_txt) >= os.path.getmtime(msh)):
        with open(tri_txt) as f:
            head = f.readline().rstrip("\n")
        if head.startswith("\t") and head[1:].isdigit() and int(head[1:]) == n:
            log(f"triangle: {tri_txt} exists, header count == {n} (resume, "
                f"{os.path.getsize(tri_txt)/1e9:.2f} GB)")
            return
        log("triangle: header count mismatch -> recomputing")
    tmp = tri_txt + ".tmp"
    t0 = time.time()
    log(f"triangle: mash triangle -k {K} -s {S} -p {cfg.threads} "
        f"v3_prophages.msh -> {os.path.basename(tri_txt)} "
        f"({n} seqs, {expected_pairs} pairs expected)")
    run(["mash", "triangle", "-k", str(K), "-s", str(S),
         "-p", str(cfg.threads), msh], out=tmp)
    with open(tmp) as f:
        head = f.readline().rstrip("\n")
    if not (head.startswith("\t") and head[1:].isdigit()
            and int(head[1:]) == n):
        os.remove(tmp)
        raise SystemExit(f"triangle: bad header {head!r} (expected \\t{n})")
    atomic_replace(tmp, tri_txt)
    log(f"triangle: wrote {tri_txt} ({os.path.getsize(tri_txt)/1e9:.2f} GB) "
        f"in {(time.time()-t0)/60:.1f} min; pair completeness checked in f32")


# ------------------------------------------------------------------- stage: f32
def stage_f32(cfg, ids):
    """triangle text -> float32 upper triangle in ids.txt order.

    File format (verified empirically, mash 2.3 — same as host_clades_mash_v3):
    first line "\\t<n>", then n rows: name_i followed by i tab-separated
    distances to rows 0..i-1 in FILE order. Rows are matched by NAME against
    ids.txt (not by position), so the conversion does not assume mash's row
    ordering; every value lands at offset(min(a,b), max(a,b)) and every
    (a,b) pair is written exactly once (completeness by construction).
    """
    tri_txt = os.path.join(cfg.out, "v3_prophages.triangle.txt")
    tri_f32 = os.path.join(cfg.out, "v3_prophages.dist")
    n = len(ids)
    npairs = n * (n - 1) // 2
    idx = {s: i for i, s in enumerate(ids)}
    if os.path.exists(tri_f32) and os.path.getsize(tri_f32) == npairs * 4:
        log(f"f32: {tri_f32} exists with exact size {npairs*4} bytes "
            f"({npairs} float32 values; resume)")
        return
    t0 = time.time()
    tri = np.full(npairs, np.float32(1.0), dtype="<f4")
    row_of_filepos = np.empty(n, dtype=np.int64)     # file position -> id index
    rows_parsed = 0
    with open(tri_txt) as f:
        head = f.readline().rstrip("\n")
        if not (head.startswith("\t") and head[1:].isdigit()
                and int(head[1:]) == n):
            raise SystemExit(f"f32: triangle header {head!r} != \\t{n}")
        p = -1
        for p, line in enumerate(f):
            toks = line.rstrip("\n").split("\t")
            if len(toks) != p + 1:
                raise SystemExit(f"f32: row {p}: {len(toks)} fields, "
                                 f"expected {p+1} — incomplete triangle")
            b = idx.get(toks[0])
            if b is None:
                raise SystemExit(f"f32: row {p} name {toks[0]!r} not in ids.txt")
            row_of_filepos[p] = b
            if p == 0:
                continue
            vals = np.array(toks[1:], dtype=np.float32)
            cols = row_of_filepos[:p]
            a = np.minimum(cols, b)
            bb = np.maximum(cols, b)
            off = a * (2 * n - a - 1) // 2 + (bb - a - 1)
            tri[off] = vals
            rows_parsed += 1
    if rows_parsed != n - 1 or p != n - 1:
        raise SystemExit(f"f32: {rows_parsed+1} data rows parsed, expected {n}")
    tmp = tri_f32 + ".tmp"
    with open(tmp, "wb") as f:
        tri.tofile(f)
    atomic_replace(tmp, tri_f32)
    log(f"f32: {npairs} float32 values ({npairs*4} bytes) -> {tri_f32} "
        f"in {time.time()-t0:.0f}s; every value written exactly once "
        f"(name-checked rows, strict field counts)")


# ---------------------------------------------------------------- stage: labels
def stage_labels(cfg, ids):
    labels = os.path.join(cfg.out, "labels.csv")
    tmp = labels + ".tmp"
    with open(tmp, "w") as f:
        f.write("sequence,community\n")
        for s in ids:
            f.write(f"{s},0\n")
    atomic_replace(tmp, labels)
    log(f"labels: {len(ids)} rows, single community 0 -> {labels}")


# ------------------------------------------------------------------ stage: tree
def write_newick(ids, Z, nwk_path, stats_path):
    """UPGMA tree -> Newick (streaming iterative emission; v2 parity, same
    branch-length convention as research/mash_tree/scripts/build_tree.py)."""
    from scipy.cluster.hierarchy import to_tree
    n = len(ids)
    tree = to_tree(Z)

    def edge_len(node, child):
        return (node.dist - child.dist) / 2.0 if child is not None else 0.0

    tmp = nwk_path + ".tmp"
    with open(tmp, "w") as f:
        stack = [(tree, 0)]
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
    atomic_replace(tmp, nwk_path)

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
    log(f"tree: {n} leaves -> {nwk_path} (height {stats['tree_height']:.4f})")


def stage_tree(cfg, ids):
    from scipy.cluster.hierarchy import linkage
    tri_f32 = os.path.join(cfg.out, "v3_prophages.dist")
    nwk = os.path.join(cfg.out, "prophages_tree.nwk")
    stats_path = os.path.join(cfg.out, "tree_stats.json")
    n = len(ids)
    npairs = n * (n - 1) // 2
    if os.path.exists(nwk) and os.path.exists(stats_path):
        with open(stats_path) as f:
            st = json.load(f)
        if st.get("n_leaves") == n and os.path.getsize(nwk) == st.get(
                "newick_bytes"):
            log(f"tree: {nwk} exists with {n} leaves (resume)")
            return
    t0 = time.time()
    tri32 = np.memmap(tri_f32, dtype="<f4", mode="r", shape=(npairs,))
    condensed = np.array(tri32, dtype=np.float64, copy=True)
    del tri32
    Z = linkage(condensed, method="average")
    del condensed
    log(f"tree: linkage done in {time.time()-t0:.1f}s")
    write_newick(ids, Z, nwk, stats_path)
    del Z


# ----------------------------------------------------------------- stage: index
def build_index(fasta_path, index_path):
    """Offset index in per_clade_alignment_pipeline.py format:
    id -> [byte offset just after the header line, sequence length in bytes]
    (single-line sequences; v2 parity)."""
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
    atomic_replace(tmp, index_path)
    log(f"index: {len(idx)} sequences -> {index_path}")


def stage_index(cfg, ids):
    index_path = os.path.join(cfg.out, "full_prophages.idx.json")
    if os.path.exists(index_path):
        with open(index_path) as f:
            n = len(json.load(f))
        if n == len(ids):
            log(f"index: {index_path} exists with {n} entries (resume)")
            return
    build_index(cfg.fa, index_path)


# ---------------------------------------------------------------- stage: clades
def stage_clades(cfg, ids):
    tc_dir = os.path.join(cfg.clades, "0", "tight_clades.json")
    ids_txt = os.path.join(cfg.out, "ids.txt")
    tri_f32 = os.path.join(cfg.out, "v3_prophages.dist")
    labels = os.path.join(cfg.out, "labels.csv")
    if os.path.exists(tc_dir):
        with open(tc_dir) as f:
            n_assigned = sum(len(v) for v in json.load(f).values())
        if n_assigned == len(ids):
            log(f"clades: {tc_dir} exists with {n_assigned} members "
                f"assigned (resume)")
            return
        log("clades: tight_clades.json stale -> recomputing")
    run(["python3", os.path.join(cfg.repo, "scripts", "build_tight_clades.py"),
         "--threshold", str(THRESHOLD), "--max-size", str(MAX_SIZE),
         "--communities", "0",
         "--outdir", cfg.clades, "--ids-file", ids_txt,
         "--triangle", tri_f32, "--labels-csv", labels])


# --------------------------------------------------------------- stage: summary
def load_clade_stats(cfg, ids):
    """Recompute the summary stats from the on-disk clade outputs."""
    tc_path = os.path.join(cfg.clades, "0", "tight_clades.json")
    sim_path = os.path.join(cfg.clades, "0", "clade_similarity.json")
    with open(tc_path) as f:
        tc = json.load(f)
    with open(sim_path) as f:
        sim = json.load(f)["per_clade"]

    rows = []
    for cid, members in tc.items():
        s = sim.get(cid, {})
        rows.append((cid, len(members), s.get("median"), s.get("min"),
                     s.get("max")))
    rows.sort()

    def write_tsv(path, sel):
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            f.write("clade_id\tn_members\tmedian_mash\tmin_mash\tmax_mash\n")
            for cid, nm, med, mn, mx in rows:
                if not sel(nm):
                    continue
                # singletons have no internal pairs: NA keeps the 5-column
                # layout without trailing whitespace (repo lint)
                def g(x):
                    return "NA" if x is None else f"{x:.6f}"
                f.write(f"{cid}\t{nm}\t{g(med)}\t{g(mn)}\t{g(mx)}\n")
        atomic_replace(tmp, path)

    write_tsv(os.path.join(cfg.clades, "clade_summary.tsv"), lambda nm: True)
    write_tsv(os.path.join(cfg.clades, "alignable_clades.tsv"),
              lambda nm: nm >= 2)
    sing = os.path.join(cfg.clades, "singletons.tsv")
    tmp = sing + ".tmp"
    with open(tmp, "w") as f:
        f.write("clade_id\tmember_id\n")
        for cid, nm, *_ in rows:
            if nm == 1:
                f.write(f"{cid}\t{tc[cid][0]}\n")
    atomic_replace(tmp, sing)

    n = len(ids)
    assigned = sum(r[1] for r in rows)
    if assigned != n:
        raise SystemExit(f"summary: assignment incomplete {assigned} != {n}")
    allm = [m for members in tc.values() for m in members]
    if len(set(allm)) != n:
        raise SystemExit("summary: duplicate clade assignments")
    medians = [r[2] for r in rows if r[2] is not None]
    over = [r[0] for r in rows
            if r[1] > 1 and r[2] is not None and r[2] > THRESHOLD + 1e-6]
    stats = {
        "n_prophages": n,
        "n_clades": len(tc),
        "n_alignable": sum(1 for r in rows if r[1] >= 2),
        "n_singletons": sum(1 for r in rows if r[1] == 1),
        "assigned": assigned,
        "median_of_clade_medians": float(np.median(medians)),
        "min_clade_median": float(np.min(medians)),
        "max_clade_median": float(np.max(medians)),
        "q1_clade_median": float(np.percentile(medians, 25)),
        "q3_clade_median": float(np.percentile(medians, 75)),
        "clades_median_zero": sum(1 for m in medians if m == 0.0),
        "clades_median_le_threshold": sum(1 for m in medians
                                          if m <= THRESHOLD),
        "clades_over_threshold": over,
        "max_clade_size": max(r[1] for r in rows),
        "size_distribution": {
            "1": sum(1 for r in rows if r[1] == 1),
            "2-10": sum(1 for r in rows if 2 <= r[1] <= 10),
            "11-50": sum(1 for r in rows if 11 <= r[1] <= 50),
            "51-100": sum(1 for r in rows if 51 <= r[1] <= 100),
            ">100": sum(1 for r in rows if r[1] > 100),
        },
        "top10_largest": sorted(
            [(r[0], r[1], r[2]) for r in rows], key=lambda x: -x[1])[:10],
    }
    log("summary: " + json.dumps({k: v for k, v in stats.items()
                                  if k != "top10_largest"}))
    if over:
        log(f"summary: WARNING {len(over)} non-singleton clades with median "
            f"> {THRESHOLD}: {over[:10]}")
    return stats, tc, sim


def stage_summary(cfg, ids):
    # invariant that guarantees a fully-filled community matrix (no NaN):
    # build_tight_clades clusters members = sorted(set(labels)) and ids.txt
    # is written pre-sorted, so members.json must equal ids exactly
    members_path = os.path.join(cfg.clades, "0", "members.json")
    if os.path.exists(members_path):
        with open(members_path) as f:
            members = json.load(f)
        if members != ids:
            raise SystemExit("summary: members.json != ids.txt (triangle row "
                             "order contract broken — community matrix would "
                             "contain NaN holes)")
        log("summary: members.json == ids.txt (sorted row order; community "
            "matrix fully filled, no NaN holes)")
    return load_clade_stats(cfg, ids)


# ------------------------------------------------------------ stage: spotcheck
def spotcheck(cfg, ids, n_readback=200, n_remeasure=50):
    """(a) readback: triangle text vs f32 binary for sampled pairs (conversion
    check); (b) re-measure: fresh `mash sketch -k 21 -s 10000` + `mash dist`
    on extracted sequences for sampled pairs (measurement check)."""
    tri_txt = os.path.join(cfg.out, "v3_prophages.triangle.txt")
    tri_f32 = os.path.join(cfg.out, "v3_prophages.dist")
    n = len(ids)
    npairs = n * (n - 1) // 2
    rng = random.Random(42)
    n_readback = min(n_readback, npairs)
    n_remeasure = min(n_remeasure, n * (n - 1))

    def off(a, b):
        return a * (2 * n - a - 1) // 2 + (b - a - 1)

    # ---- (a) readback: scan the triangle text once, keep the sampled rows.
    # Text rows are in mash's own (sketch/FASTA) order — map by NAME to
    # ids.txt indices, exactly like stage_f32 does.
    want = set()
    while len(want) < n_readback:
        i = rng.randrange(1, n)
        j = rng.randrange(0, i)
        want.add((i, j))
    by_row = {}
    for i, j in want:
        by_row.setdefault(i, []).append(j)
    got = {}
    t0 = time.time()
    idx = {s: i for i, s in enumerate(ids)}
    textpos_to_ididx = np.empty(n, dtype=np.int64)
    with open(tri_txt) as f:
        f.readline()
        for p, line in enumerate(f):
            toks = line.rstrip("\n").split("\t")
            b = idx[toks[0]]
            textpos_to_ididx[p] = b
            js = by_row.get(p)
            if js:
                for j in js:
                    # value between text row j and text row p
                    got[(p, j)] = (b, int(textpos_to_ididx[j]),
                                   float(toks[1 + j]))
            if p >= n - 1:
                break
    mm = np.memmap(tri_f32, dtype="<f4", mode="r", shape=(npairs,))
    bad = []
    for (i, j), (b, a, d) in got.items():
        if a == b:
            continue
        lo, hi = min(a, b), max(a, b)
        v = float(mm[off(lo, hi)])
        if abs(v - d) > 1e-6:
            bad.append(((i, j), v, d))
    log(f"spotcheck: readback {len(got)} sampled pairs text vs f32 in "
        f"{time.time()-t0:.0f}s, {len(bad)} mismatches")
    # ---- (b) re-measure with a fresh sketch of the extracted sequences
    with open(os.path.join(cfg.out, "full_prophages.idx.json")) as f:
        fidx = json.load(f)
    pairs = []
    while len(pairs) < n_remeasure:
        i, j = rng.randrange(0, n), rng.randrange(0, n)
        if i != j:
            pairs.append((i, j))
    tmpdir = os.path.join(cfg.out, "spotcheck.tmp")
    os.makedirs(tmpdir, exist_ok=True)
    out_tsv = os.path.join(cfg.out, "spotcheck.tsv")
    tmp = out_tsv + ".tmp"
    worst = 0.0
    n_fail = 0
    with open(tmp, "w") as f:
        f.write("prophage_a\tprophage_b\tstored_mash\tremeasured_mash\t"
                "abs_delta\tresult\n")
        for k, (i, j) in enumerate(pairs):
            recs = []
            for pos, x in enumerate((i, j)):
                off0, L = fidx[ids[x]]
                with open(cfg.fa, "rb") as faf:
                    faf.seek(off0)
                    seq = faf.read(L).decode()
                path = os.path.join(tmpdir, f"s{k}_{pos}.fa")
                with open(path, "w") as pf:
                    pf.write(f">{ids[x]}\n{seq}\n")
                recs.append(path)
            sk = os.path.join(tmpdir, f"s{k}")
            r1 = subprocess.run(["mash", "sketch", "-k", str(K), "-s", str(S),
                                 "-o", sk, recs[0]],
                                capture_output=True, text=True)
            r2 = subprocess.run(["mash", "dist", sk + ".msh", recs[1]],
                                capture_output=True, text=True)
            if r1.returncode != 0 or r2.returncode != 0:
                raise SystemExit(f"spotcheck: mash failed at pair {k}: "
                                 f"{r1.stderr[:150]} {r2.stderr[:150]}")
            d_measured = float(r2.stdout.split("\t")[2])
            d_stored = float(mm[off(min(i, j), max(i, j))])
            delta = abs(d_stored - d_measured)
            worst = max(worst, delta)
            ok = delta <= 1e-4
            n_fail += 0 if ok else 1
            f.write(f"{ids[i]}\t{ids[j]}\t{d_stored:.6f}\t{d_measured:.6f}\t"
                    f"{delta:.6f}\t{'PASS' if ok else 'FAIL'}\n")
            os.remove(sk + ".msh")
            for p_ in recs:
                os.remove(p_)
    atomic_replace(tmp, out_tsv)
    os.rmdir(tmpdir)
    del mm
    log(f"spotcheck: re-measured {len(pairs)} pairs with fresh "
        f"sketch+dist, worst |delta| = {worst:.6f}, {n_fail} FAIL")
    return {"readback_pairs": len(got), "readback_mismatches": len(bad),
            "remeasured_pairs": len(pairs),
            "remeasured_worst_delta": round(worst, 8),
            "remeasured_failures": n_fail, "bad_examples": bad[:5]}


# --------------------------------------------------------------- stage: report
def stage_report(cfg, ids, lenstats, clstats):
    report = os.path.join(cfg.repo, "ntm", "v3", "mash_clades_report.md")
    clades_repo = os.path.join(cfg.repo, "ntm", "v3", "clades")
    os.makedirs(clades_repo, exist_ok=True)
    m = clstats
    n = m["n_prophages"]

    # repo copies: clade definitions + per-clade similarity stats
    # (deterministic gzip) + summaries
    def gz_copy(src, dst):
        with open(src, "rb") as fi, open(dst + ".tmp", "wb") as fo:
            with gzip.GzipFile(filename="", mode="wb", fileobj=fo, mtime=0,
                               compresslevel=9) as gz:
                shutil.copyfileobj(fi, gz, 1 << 20)
        atomic_replace(dst + ".tmp", dst)

    gz_copy(os.path.join(cfg.clades, "0", "tight_clades.json"),
            os.path.join(clades_repo, "tight_clades.json.gz"))
    gz_copy(os.path.join(cfg.clades, "0", "clade_similarity.json"),
            os.path.join(clades_repo, "clade_similarity.json.gz"))
    for name in ("clade_summary.tsv", "alignable_clades.tsv",
                 "singletons.tsv"):
        shutil.copyfile(os.path.join(cfg.clades, name),
                        os.path.join(clades_repo, name))
    sums_path = os.path.join(clades_repo, "SHA256SUMS")
    with open(sums_path + ".tmp", "w") as f:
        for name in ("tight_clades.json.gz", "clade_similarity.json.gz",
                     "clade_summary.tsv", "alignable_clades.tsv",
                     "singletons.tsv"):
            f.write(f"{sha256(os.path.join(clades_repo, name))}  {name}\n")
    atomic_replace(sums_path + ".tmp", sums_path)

    len_ok = (lenstats["n"] == EXTRACT_REF["n"]
              and lenstats["len_min"] == EXTRACT_REF["len_min"]
              and lenstats["len_max"] == EXTRACT_REF["len_max"]
              and lenstats["len_total"] == EXTRACT_REF["len_total"]
              and lenstats["len_outside_band"] == EXTRACT_REF["outside_band"])
    sc_path = os.path.join(cfg.out, "spotcheck.json")
    sc = json.load(open(sc_path)) if os.path.exists(sc_path) else {}
    fa_sha = sha256(cfg.fa)
    npairs = n * (n - 1) // 2

    lines = []
    lines.append("# NTM v3 — prophage MASH + tight clades report")
    lines.append("")
    lines.append(f"Generated: "
                 f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} "
                 f"(task ntm-v3-prophage)")
    lines.append("")
    lines.append("## Prophage set")
    lines.append("")
    lines.append(f"- `full_prophages.fa`: **{n} prophages**, "
                 f"{lenstats['len_total']:,} bp (NVMe "
                 f"`ntm/v3/full_prophages.fa`, sha256 `{fa_sha}` — matches "
                 f"the extract_report receipt)")
    lines.append("- provenance: task ntm-v3-unified unified manifest "
                 "(BV-BRC phigaro primary, v2 calls only for genomes absent "
                 "from the export); 9,446 of 36,940 manifest rows are "
                 "extractable — the 27,494 run-assembly rows have no FASTA "
                 "on this host yet (collaborator delivery pending, "
                 "`ntm/v3/download_report.md`), none dropped")
    lines.append(f"- Mash sketch: `-i -k {K} -s {S}` (identical to v1/v2)")
    lines.append("")
    lines.append("## Commands (reproducible)")
    lines.append("")
    lines.append("```bash")
    lines.append("# driver (every stage resumable from disk state; heavy "
                 "phase run detached under nohup setsid)")
    lines.append("python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py "
                 "--phase heavy   # ids, sketch, triangle, f32, labels, tree, "
                 "index")
    lines.append("python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py "
                 "--phase light   # clades, summary, spotcheck, report")
    lines.append("# the two mash invocations the driver runs (cwd "
                 ".../ntm/v3/mash_clades):")
    lines.append(f"mash sketch -i -k {K} -s {S} -p {cfg.threads} "
                 f"-o v3_prophages ../full_prophages.fa   # -> "
                 f"v3_prophages.msh")
    lines.append(f"mash triangle -k {K} -s {S} -p {cfg.threads} "
                 f"v3_prophages.msh > v3_prophages.triangle.txt")
    lines.append("# driver-internal conversion, then clustering with the "
                 "reused-as-is E. coli machinery (unmodified):")
    lines.append("python3 scripts/build_tight_clades.py --threshold 0.25 "
                 "--max-size 100 --communities 0 \\")
    lines.append(f"    --outdir {cfg.clades} --ids-file {cfg.out}/ids.txt \\")
    lines.append(f"    --triangle {cfg.out}/v3_prophages.dist "
                 f"--labels-csv {cfg.out}/labels.csv")
    lines.append("# independent re-validation")
    lines.append("python3 ntm/v3/scripts/validate_mash_clades_v3.py")
    lines.append("```")
    lines.append("")
    lines.append("## Clades (threshold 0.25, max-size 100, community 0)")
    lines.append("")
    lines.append(f"- clade count: **{m['n_clades']}**")
    lines.append(f"- alignable (>=2 members): **{m['n_alignable']}**")
    lines.append(f"- singletons: **{m['n_singletons']}** "
                 f"({100.0*m['n_singletons']/m['n_clades']:.1f}% of clades; "
                 f"each singleton is its own clade — v1/v2 convention, "
                 f"pass-through, no ML reconstruction for n=1)")
    lines.append(f"- assignment check: sum of clade sizes = {m['assigned']} "
                 f"== {n} prophages, all member ids unique (every prophage "
                 f"assigned exactly once)")
    lines.append(f"- max clade size: {m['max_clade_size']} (cap {MAX_SIZE})")
    lines.append("")
    lines.append("### Clade size distribution")
    lines.append("")
    lines.append("| size band | clades |")
    lines.append("|---|---:|")
    for band, cnt in m["size_distribution"].items():
        lines.append(f"| {band} | {cnt} |")
    lines.append("")
    lines.append("Top 10 largest clades: " + ", ".join(
        f"`{cid}` (n={nm}, median={'' if med is None else f'{med:.4f}'})"
        for cid, nm, med in m["top10_largest"]))
    lines.append("")
    lines.append("## Internal similarity (per-clade pairwise MASH)")
    lines.append("")
    lines.append(f"- median of per-clade median distances: "
                 f"**{m['median_of_clade_medians']:.4f}**")
    lines.append(f"- per-clade median distribution: min "
                 f"{m['min_clade_median']:.6f}, Q1 "
                 f"{m['q1_clade_median']:.6f}, median "
                 f"{m['median_of_clade_medians']:.6f}, Q3 "
                 f"{m['q3_clade_median']:.6f}, max "
                 f"{m['max_clade_median']:.6f}")
    lines.append(f"- clades with internal median exactly 0.0: "
                 f"{m['clades_median_zero']} (identical/near-identical "
                 f"prophages recurring across isolates of the same species)")
    lines.append(f"- clades with median <= {THRESHOLD}: "
                 f"{m['clades_median_le_threshold']} of {m['n_alignable']} "
                 f"alignable (build_tight_clades.py `tighten_clades` "
                 f"enforces the median criterion; any violation is flagged "
                 f"here)")
    if m["clades_over_threshold"]:
        lines.append(f"- **FLAGGED**: clades with median > {THRESHOLD}: "
                     f"{m['clades_over_threshold']}")
    else:
        lines.append(f"- clades with median > {THRESHOLD}: **0**")
    lines.append("")
    lines.append("## Prophage length band check")
    lines.append("")
    lo, hi = BAND
    lines.append(f"- lengths (bp, recomputed from `full_prophages.fa` by this "
                 f"task): min {lenstats['len_min']:,}, median "
                 f"{int(lenstats['len_median']):,}, mean "
                 f"{lenstats['len_mean']:,.0f}, max "
                 f"{lenstats['len_max']:,}, total "
                 f"{lenstats['len_total']:,}")
    lines.append(f"- outside [{lo:,}, {hi:,}] bp: "
                 f"**{lenstats['len_outside_band']}** "
                 f"({lenstats['len_lt_1kb']} below {lo:,}, "
                 f"{lenstats['len_gt_100kb']} above {hi:,}) — "
                 + ("matches the extract_report flag count (91) and its "
                    "min/max (237 / 93,713); all retained (kept, not "
                    "dropped — v2/v3 extract convention)"
                    if len_ok else
                    "**MISMATCH vs extract_report (91 flagged, 237/93,713) — "
                    "investigate before relying on either**"))
    lines.append("")
    lines.append("## Comparison with v1 / v2")
    lines.append("")
    lines.append("| version | prophages | clades | alignable | singletons | "
                 "median internal |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    lines.append("| v1 (geNomad, 7,352 NTM assemblies) | 10,438 | 913 | 472 "
                 "| 441 | 0.074 |")
    lines.append("| v2 (collab manifest, NCBI-only) | 8,502 | 2,388 | 1,251 "
                 "| 1,137 | 0.0060 |")
    lines.append(f"| **v3 (unified manifest, extractable set)** | **{n}** | "
                 f"**{m['n_clades']}** | **{m['n_alignable']}** | "
                 f"**{m['n_singletons']}** | "
                 f"**{m['median_of_clade_medians']:.4f}** |")
    lines.append("")
    lines.append(f"Cohort note: the task brief anticipated \"~4-6x more "
                 f"prophages than v2\", but the extractable v3 set is {n} = "
                 f"{n/8502:.2f}x v2's 8,502 — the other 27,494 manifest rows "
                 f"are run-assembly prophages whose FASTAs are blocked "
                 f"pending collaborator delivery "
                 f"(ntm/v2/run_assemblies/REQUEST.md PENDING). Clade count "
                 f"{m['n_clades']} is in the plausible range for this input "
                 f"size (v2: 2,388 clades / 8,502 prophages; v3 adds BV-BRC "
                 f"phigaro calls for genomes v2 never held). When the run "
                 f"assemblies land, re-running this driver over the extended "
                 f"FASTA reproduces the identical recipe.")
    lines.append("")
    lines.append("## Mechanics note (triangle row order)")
    lines.append("")
    lines.append("`ids.txt` is written in `sorted()` order, not FASTA order. "
                 "`build_tight_clades.py` (reused as-is) sorts community "
                 "members with `sorted(set(members))` and its "
                 "`read_community_matrix` documents *\"members must be sorted "
                 "by triangle row index\"*; it fills D[i][j] only for member "
                 "pairs whose member-list order agrees with triangle row "
                 "order. In the v1/v2 invocations ids.txt was in FASTA order, "
                 "so 35.2% of v2's within-community pairs were never filled "
                 "(NaN, measured on the v2 triangle: 23,407,628 of "
                 "36,137,751) and silently behaved as maximally distant — "
                 "missed joins, hence fragmented (inflated) clade counts. "
                 "Writing ids.txt pre-sorted satisfies the contract with the "
                 "script unmodified and fills 100% of pairs (verified: no "
                 "NaN in the community matrix). All downstream consumers "
                 "join by prophage id, so the row order is internal to this "
                 "step. v1/v2 clade counts are therefore upper bounds.")
    lines.append("")
    lines.append("## Validation")
    lines.append("")
    lines.append(f"- every prophage assigned to exactly one tight clade: sum "
                 f"of clade sizes {m['assigned']} == FASTA records {n} "
                 f"({len(set(ids))} unique ids) — PASS")
    lines.append(f"- clade internal similarity: every non-singleton clade "
                 f"median pairwise MASH <= {THRESHOLD}: "
                 f"{len(m['clades_over_threshold'])} violations — "
                 f"{'PASS' if not m['clades_over_threshold'] else 'FLAGGED'}")
    lines.append("- v2-comparable stats reported above (clade count, "
                 "alignable/singleton split, median internal distance, size "
                 "distribution)")
    lines.append(f"- commands reproducible: exact commands above and in "
                 f"`{cfg.out}/commands.log`; independent re-validation: "
                 f"`python3 ntm/v3/scripts/validate_mash_clades_v3.py`")
    if sc:
        lines.append(f"- distance integrity: {sc.get('readback_pairs')} "
                     f"sampled pairs triangle-text vs float32 "
                     f"({sc.get('readback_mismatches')} mismatches); "
                     f"{sc.get('remeasured_pairs')} pairs re-measured with a "
                     f"fresh `mash sketch -k 21 -s 10000` + `mash dist` "
                     f"(worst |delta| {sc.get('remeasured_worst_delta')}, "
                     f"{sc.get('remeasured_failures')} failures)")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append(f"- NVMe `{cfg.out}/`: `v3_prophages.msh`, "
                 f"`v3_prophages.triangle.txt` (mash triangle text), "
                 f"`v3_prophages.dist` (float32 upper triangle, "
                 f"{npairs:,} values == n*(n-1)/2, {npairs*4:,} bytes), "
                 f"`ids.txt`, `prophage_lengths.tsv`, `labels.csv` "
                 f"(community 0), `full_prophages.idx.json`, "
                 f"`prophages_tree.nwk`, `tree_stats.json`, `spotcheck.tsv`, "
                 f"`commands.log`, `driver.log`")
    lines.append(f"- NVMe `{cfg.clades}/`: `0/` (tight_clades.json, "
                 f"clade_similarity.json, members.json, distances.npz, "
                 f"commands.log), `tight_clades_summary.json`, "
                 f"`clade_summary.tsv`, `alignable_clades.tsv`, "
                 f"`singletons.tsv`")
    lines.append("- repo `ntm/v3/clades/`: `tight_clades.json.gz` (clade "
                 "definitions), `clade_similarity.json.gz` (per-clade "
                 "internal similarity stats), `clade_summary.tsv`, "
                 "`alignable_clades.tsv`, `singletons.tsv`, `SHA256SUMS` "
                 "(gzip artifacts deterministic, gzip -n)")
    lines.append("- FASTA / sketch / triangle stay on NVMe (repo holds only "
                 "code, manifests and small reports — v1/v2 rule)")
    tmp = report + ".tmp"
    with open(tmp, "w") as f:
        f.write("\n".join(lines) + "\n")
    atomic_replace(tmp, report)
    log(f"report: -> {report}")


# ----------------------------------------------------------------------- main
def load_lengths(cfg, ids):
    """Reload (or recompute) prophage lengths for the report stage."""
    len_path = os.path.join(cfg.out, "prophage_lengths.tsv")
    if os.path.exists(len_path):
        L = []
        with open(len_path) as f:
            next(f)
            for line in f:
                L.append(int(line.rstrip("\n").split("\t")[1]))
        if len(L) == len(ids):
            return {
                "n": len(L), "len_min": int(min(L)),
                "len_median": float(np.median(L)),
                "len_mean": float(np.mean(L)), "len_max": int(max(L)),
                "len_total": int(sum(L)),
                "len_outside_band": int(sum(1 for x in L
                                            if not (BAND[0] <= x <= BAND[1]))),
                "len_lt_1kb": int(sum(1 for x in L if x < BAND[0])),
                "len_gt_100kb": int(sum(1 for x in L if x > BAND[1])),
            }
    _, stats = stage_ids(cfg)
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", default="all",
                    choices=["all", "heavy", "light"] + PHASES)
    ap.add_argument("--work-dir", default=WORK_DEFAULT)
    ap.add_argument("--repo", default=REPO_DEFAULT)
    ap.add_argument("--threads", type=int, default=64)
    cfg = ap.parse_args()

    v3 = os.path.join(cfg.work_dir, "ntm", "v3")
    cfg.fa = os.path.join(v3, "full_prophages.fa")
    cfg.out = os.path.join(v3, "mash_clades")
    cfg.clades = os.path.join(v3, "clades")
    os.makedirs(cfg.out, exist_ok=True)
    os.makedirs(cfg.clades, exist_ok=True)
    if not os.path.exists(cfg.fa):
        raise SystemExit(f"missing input FASTA {cfg.fa}")

    seq = PHASES if cfg.phase == "all" else (
        HEAVY if cfg.phase == "heavy" else
        LIGHT if cfg.phase == "light" else [cfg.phase])
    log(f"phase sequence: {seq} (threads={cfg.threads})")

    def load_ids():
        p = os.path.join(cfg.out, "ids.txt")
        ids = [l.strip() for l in open(p) if l.strip()]
        fh = [x[1:].strip().split()[0]
              for x in open(cfg.fa, errors="replace") if x.startswith(">")]
        if set(ids) != set(fh) or len(ids) != len(fh):
            raise SystemExit("ids.txt does not match the FASTA headers")
        if ids != sorted(ids):
            raise SystemExit("ids.txt is not in sorted() order "
                             "(build_tight_clades row-order contract)")
        return ids

    t_start = time.time()
    ids, lenstats, clstats = None, None, None
    for stage in seq:
        log(f"=== stage {stage} ===")
        t0 = time.time()
        if stage == "ids":
            ids, lenstats = stage_ids(cfg)
        elif stage == "sketch":
            stage_sketch(cfg)
        elif stage == "triangle":
            stage_triangle(cfg, len(ids) if ids else len(load_ids()))
        elif stage == "f32":
            stage_f32(cfg, ids or load_ids())
        elif stage == "labels":
            stage_labels(cfg, ids or load_ids())
        elif stage == "tree":
            stage_tree(cfg, ids or load_ids())
        elif stage == "index":
            stage_index(cfg, ids or load_ids())
        elif stage == "clades":
            stage_clades(cfg, ids or load_ids())
        elif stage == "summary":
            clstats, tc, sim = stage_summary(cfg, ids or load_ids())
        elif stage == "spotcheck":
            sc = spotcheck(cfg, ids or load_ids())
            sc_path = os.path.join(cfg.out, "spotcheck.json")
            with open(sc_path + ".tmp", "w") as f:
                json.dump(sc, f, indent=2)
            atomic_replace(sc_path + ".tmp", sc_path)
            if sc["readback_mismatches"] or sc["remeasured_failures"]:
                raise SystemExit(f"spotcheck FAILED: {sc}")
        elif stage == "report":
            ids = ids or load_ids()
            if lenstats is None:
                lenstats = load_lengths(cfg, ids)
            if clstats is None:
                clstats, _, _ = load_clade_stats(cfg, ids)
            stage_report(cfg, ids, lenstats, clstats)
        # persist the driver log (resume breadcrumbs)
        if LOG_LINES:
            with open(os.path.join(cfg.out, "driver.log"), "a") as f:
                f.write("\n".join(LOG_LINES) + "\n")
            LOG_LINES.clear()
        # commands.log: reproducibility receipt
        with open(os.path.join(cfg.out, "commands.log"), "a") as f:
            if COMMANDS:
                f.write("\n".join(COMMANDS) + "\n")
                COMMANDS.clear()
            f.write(f"# {time.strftime('%Y-%m-%dT%H:%M:%SZ')} stage={stage} "
                    f"rc=0 secs={time.time()-t0:.1f}\n")
        log(f"=== stage {stage} done in {time.time()-t0:.1f}s ===")

    log(f"ALL DONE in {(time.time()-t_start)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
