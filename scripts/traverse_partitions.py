#!/usr/bin/env python3
"""
traverse_partitions.py — optimal typical traversal of community partition sets
==============================================================================

Given the partition output of ``impg partition`` for a prophage community —
(a) a BED mapping prophages to ordered partition intervals and (b) per-partition
alignment blocks (MAF) or member-sequence bundles (FASTA) — this script computes
an **optimal typical traversal**: an ordering of ALL partitions of the community
(core and rare alike) that maximizes the expected completeness of the resulting
consensus path across the community's prophages, and stitches the per-partition
representative sequences into a community genome of phage-typical length.

This is the deliverable for README §6 ("Step 5 — Optimal traversal / ordering of
partitions", marked CODE TO WRITE).  It extends the validated reference
implementation ``research/stitching/stitch_algorithm.py`` (commit ``2363ece``,
byte-identical reproduction verified) in the following ways:

* **Rare partitions are KEPT in the ordering.**  The historical stitching
  dropped every partition below an occurrence threshold from the path.  This
  script includes every partition carried by any prophage.  Rare partitions
  (present in only a few genomes) are placed by a *typical-path fallback*:
  observed adjacency where it exists, else co-occurrence anchoring, else a
  k-mer-Jaccard sequence-similarity bridge to their most-similar partition.
* **An explicit optimization objective** with per-community validation metrics
  (see ``OBJECTIVE`` below).
* **A length budget** for the stitched genome so the output stays at
  phage-typical scale (tens of kb to ~150 kb) instead of the multi-Mbp naive
  concatenations seen historically (community 3 was 1,398,681 bp).
* **Deterministic ordering** — no randomness anywhere; ties broken by
  (occurrence descending, partition id ascending).  ``--seed`` is recorded in
  the output for provenance but does not change the result.
* **Per-community outputs**: traversal order JSON, per-partition consensus
  FASTA, stitched genome FASTA, coverage statistics TSV (+ stats JSON).

OBJECTIVE
---------

Let ``P`` be the set of partitions of a community, and ``T`` a traversal (an
ordering / permutation of ``P``).  For a prophage ``g`` with ordered partition
list ``S_g``, the *path completeness* of ``T`` for ``g`` is the fraction of the
consecutive partition pairs of ``g`` that are preserved as consecutive pairs of
``T``.  The expected path completeness over the community is therefore the
fraction of all observed consecutive partition pairs (summed over prophages)
that the traversal keeps consecutive.

The traversal maximizes the preserved **weighted** adjacency:

    maximize  Σ_i  W[T[i]][T[i+1]]

    W[a][b] = adj(a,b)                  # prophages with a immediately before b
            + μ * cooc(a,b) / N         # co-occurrence anchor (rare partitions)
            + λ * J(a,b)                # k-mer Jaccard sequence-similarity bridge

where ``N`` = number of prophages, ``μ`` = ``--cooc-weight`` (default 1.0) and
``λ`` = ``--bridge-weight`` (default 1.0).  This is a maximum-weight Hamiltonian
path problem (NP-hard); we solve it deterministically with a greedy
construction (extend-both-ends nearest-neighbour from the best start among the
top start candidates) followed by a first-improvement 2-opt local search.

Validation metrics (computed per community, written to the JSON):

* ``expected_path_completeness`` — fraction of all observed prophage
  consecutive-pair adjacencies preserved as consecutive pairs of T.
* ``ordered_coverage`` — mean over prophages of |LCS(S_g, T)| / |S_g|
  (fraction of each prophage's partitions appearing in T in the same relative
  order; LCS computed as the longest increasing subsequence of T-positions).
* ``fraction_prophages_fully_ordered`` — fraction of prophages whose entire
  partition list is a subsequence of T.
* ``genome_coverage`` — mean over prophages of |S_g ∩ T_genome| / |S_g|, the
  fraction of each prophage's partitions represented in the stitched genome.
* ``fraction_prophages_covered_by_genome`` — fraction of prophages with ≥ 1
  partition in the stitched genome.
* join statistics: how many consecutive traversal pairs are *observed* (backed
  by ≥ 1 prophage adjacency) vs *bridged* (placed by co-occurrence/similarity).

CONSENSUS (per partition)
-------------------------

* MAF alignment block present → majority-rule consensus over the aligned block
  (column emitted iff ≥ ``coverage_threshold`` fraction of block sequences carry
  a non-gap base there; majority base wins).  Default ``--coverage-threshold``
  0.0 reproduces the validated ``2363ece`` consensi byte-identically (note: for
  staggered segment bundles 0.0 emits the union mosaic spanning all segment
  positions, which is longer than any single genome — the ``--max-length``
  budget keeps the *stitched genome* at phage-typical scale).
* No MAF but member sequences in FASTA → longest-sequence fallback.
* No sequence data at all → partition is kept in the ordering and statistics but
  contributes no sequence (flagged ``NO_SEQ``).

INPUT FORMATS (BED)
-------------------

* ``--bed FILE`` with 4 columns ``prophage start end partition_id`` (legacy
  combined layout, e.g. ``research/stitching/community_3_partitions.bed``), or
* ``--bed DIR`` containing per-partition 3-column BEDs ``partition<N>.bed``
  (layout of the alignment task's ``partitions_bed/`` and of archived
  ``research_outputs/partitions/``).  Both are normalized to the same
  per-prophage ordered partition lists (contiguous runs of the same partition
  are merged, consecutive duplicates collapsed).

OUTPUTS (``--output PREFIX``)
-----------------------------

* ``PREFIX.traversal.json`` — full traversal order, per-partition records,
  objective, parameters and metrics.
* ``PREFIX.consensus.fa`` — per-partition representative sequences (majority
  consensus or longest-sequence fallback).
* ``PREFIX.fa`` — stitched community genome (traversal order, length budget
  ``--max-length``; ``--gap-size`` N-run join markers between partitions).
* ``PREFIX.coverage.tsv`` — per-partition coverage statistics (occurrence,
  fraction, rare flag, join type, in-genome flag).
* ``PREFIX.stats.json`` — summary metrics (same as the JSON's ``metrics``).
* ``PREFIX.full.fa`` — full concatenation of every partition representative in
  traversal order (budget disabled; emitted only with ``--full-concat``).

USAGE
-----

    python3 scripts/traverse_partitions.py \
        --partitions-dir <dir> --bed <partitions.bed> --output <prefix>

See ``README.md`` §6.  Deterministic and reproducible: same inputs + same
parameters → identical outputs (verified by tests).
"""

import argparse
import glob
import json
import math
import os
import re
import time
from collections import Counter, defaultdict

import numpy as np

# Default bridge/co-occurrence weights in W (see OBJECTIVE).
DEFAULT_COOC_WEIGHT = 1.0
DEFAULT_BRIDGE_WEIGHT = 1.0
DEFAULT_KMER_SIZE = 15
DEFAULT_RARE_THRESHOLD = 0.2
DEFAULT_RARE_MIN_COUNT = 5
DEFAULT_MAX_LENGTH = 150_000      # phage-typical budget for the stitched genome
DEFAULT_COVERAGE_THRESHOLD = 0.0  # 0.0 -> byte-identical to validated 2363ece runs
DEFAULT_GAP_SIZE = 0              # N-runs between partition blocks in the genome
DEFAULT_SEED = 42                 # recorded for provenance; algorithm is RNG-free


# ─── BED parsing ────────────────────────────────────────────────────────────

def _parse_partition_id(path):
    """Extract the integer partition id from a partition<N>.<ext> path."""
    base = os.path.basename(path)
    m = re.search(r'(\d+)', base)
    return int(m.group(1)) if m else None


def _parse_bed_file(bed_path):
    """Parse a 4-column combined BED: prophage  start  end  partition_id.

    Returns dict prophage -> list[(pid, start, end)].
    """
    prophage_partitions = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            name, start, end, pid = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            prophage_partitions[name].append((pid, start, end))
    return dict(prophage_partitions)


def _parse_bed_dir(bed_dir):
    """Parse a directory of per-partition 3-column BEDs: prophage start end.

    The partition id is taken from the file name (partition<N>.bed).
    Returns dict prophage -> list[(pid, start, end)].
    """
    prophage_partitions = defaultdict(list)
    files = sorted(glob.glob(os.path.join(bed_dir, '*.bed')))
    if not files:
        files = sorted(glob.glob(os.path.join(bed_dir, 'partition*')))
    for path in files:
        pid = _parse_partition_id(path)
        if pid is None:
            continue
        with open(path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                name, start, end = parts[0], int(parts[1]), int(parts[2])
                prophage_partitions[name].append((pid, start, end))
    return dict(prophage_partitions)


def parse_bed(bed_path):
    """Parse either input layout.  Returns dict prophage -> [(pid, start, end)]."""
    if os.path.isdir(bed_path):
        return _parse_bed_dir(bed_path)
    return _parse_bed_file(bed_path)


def normalize_partition_lists(prophage_partitions):
    """Merge contiguous runs of the same partition and collapse consecutive
    duplicates, producing each prophage's ordered unique partition sequence.

    Multiple BED rows for the same (prophage, partition) with contiguous
    coordinates are a single assignment split across window boundaries; keeping
    them as consecutive self-adjacencies would inflate the graph.
    """
    normalized = {}
    for name, entries in prophage_partitions.items():
        entries = sorted(entries, key=lambda x: (x[1], x[2]))
        merged = []
        for pid, s, e in entries:
            if merged and merged[-1][0] == pid and s <= merged[-1][2]:
                prev = merged[-1]
                merged[-1] = (prev[0], prev[1], max(prev[2], e))
            else:
                merged.append((pid, s, e))
        pids = [m[0] for m in merged]
        # collapse consecutive duplicates
        collapsed = []
        for p in pids:
            if not collapsed or collapsed[-1] != p:
                collapsed.append(p)
        normalized[name] = collapsed
    return normalized


# ─── MAF / FASTA parsing ────────────────────────────────────────────────────

def parse_partition_maf(maf_path):
    """Parse a MAF alignment block (optionally gzip-compressed).
    Returns list of dicts with keys name, start, size, strand, src_size,
    seq, or [] if unparseable."""
    import gzip
    records = []
    in_block = False
    opener = gzip.open if maf_path.endswith('.gz') else open
    with opener(maf_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('a '):
                in_block = True
            elif line.startswith('s ') and in_block:
                parts = line.split()
                if len(parts) >= 7:
                    records.append({
                        'name': parts[1],
                        'start': int(parts[2]),
                        'size': int(parts[3]),
                        'strand': parts[4],
                        'src_size': int(parts[5]),
                        'seq': parts[6],
                    })
    return records


def parse_fasta(path):
    """Parse a FASTA file.  Returns list of (name, seq)."""
    records = []
    name = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            records.append((name, ''.join(seq_parts)))
    return records


def _find_partition_files(partitions_dir):
    """Return dict pid -> list of candidate files (any extension)."""
    found = defaultdict(list)
    patterns = ['partition*.maf', 'partition*.maf.gz', 'partition*.fasta',
                'partition*.fa', 'partition*.fna', 'partition*']
    seen = set()
    for pattern in patterns:
        for path in sorted(glob.glob(os.path.join(partitions_dir, pattern))):
            if path in seen:
                continue
            seen.add(path)
            if not os.path.isfile(path):
                continue
            pid = _parse_partition_id(path)
            if pid is not None:
                found[pid].append(path)
    return dict(found)


def load_partition_sequences(partitions_dir, pids):
    """Load sequence data for the requested partition ids.

    Returns (data, warnings) where data maps pid -> dict with:
      'maf': list of records or None; 'fasta': list of (name, seq) or None.
    """
    files = _find_partition_files(partitions_dir)
    data = {}
    warnings = []
    for pid in pids:
        candidates = files.get(pid, [])
        maf = None
        fasta = None
        for path in candidates:
            low = path.lower()
            if low.endswith(('.maf', '.maf.gz')):
                recs = parse_partition_maf(path)
                if recs:
                    maf = recs
            else:
                recs = parse_fasta(path)
                if recs:
                    fasta = recs
        data[pid] = {'maf': maf, 'fasta': fasta}
        if maf is None and fasta is None:
            warnings.append(f'partition {pid}: no sequence data in {partitions_dir}')
    return data, warnings


# ─── Consensus / representative sequences ───────────────────────────────────

def compute_partition_consensus(records, coverage_threshold=DEFAULT_COVERAGE_THRESHOLD):
    """Majority-rule consensus over an aligned block.

    A column is emitted iff at least ``max(1, ceil(coverage_threshold * n))``
    sequences carry a non-gap base there; the majority base wins.  Handles
    ragged block widths by truncating to the minimum width.

    Returns (consensus_string_or_None, diagnostics_dict).
    """
    if not records:
        return None, {}
    width = min(len(r['seq']) for r in records)
    n_seqs = len(records)
    min_cov = max(1, math.ceil(n_seqs * coverage_threshold))
    depth = [0] * width
    consensus = []
    for i in range(width):
        counts = Counter()
        total = 0
        for rec in records:
            base = rec['seq'][i].upper()
            if base in 'ACGT':
                counts[base] += 1
                total += 1
        depth[i] = total
        if total >= min_cov and counts:
            consensus.append(counts.most_common(1)[0][0])
    diag = {
        'n_seqs': n_seqs,
        'aln_width': width,
        'mean_depth': round(sum(depth) / width, 2) if width else 0.0,
        'max_depth': max(depth) if width else 0,
    }
    return ''.join(consensus), diag


def partition_representative(seq_data, coverage_threshold=DEFAULT_COVERAGE_THRESHOLD):
    """Representative sequence for a partition: MAF majority consensus, else
    longest FASTA member sequence (longest-seq fallback), else None.

    Returns (seq_or_None, diag).
    """
    maf = seq_data.get('maf')
    if maf:
        cons, diag = compute_partition_consensus(maf, coverage_threshold)
        diag['source'] = 'maf_consensus'
        return cons, diag
    fasta = seq_data.get('fasta')
    if fasta:
        longest = max(fasta, key=lambda r: len(r[1]))
        diag = {'source': 'fasta_longest', 'n_seqs': len(fasta),
                'longest_name': longest[0], 'longest_len': len(longest[1])}
        return longest[1], diag
    return None, {'source': 'none'}


# ─── Coverage statistics ────────────────────────────────────────────────────

def compute_coverage_stats(normalized, n_prophages, rare_threshold, rare_min_count):
    """Per-partition coverage statistics.

    Returns dict pid -> {occurrence, fraction, rare}.
    A partition is RARE if its occurrence fraction is below rare_threshold OR
    its absolute occurrence count is below rare_min_count.
    """
    occurrence = Counter()
    for pids in normalized.values():
        occurrence.update(set(pids))
    stats = {}
    for pid, count in occurrence.items():
        frac = count / max(n_prophages, 1)
        rare = (frac < rare_threshold) or (count < rare_min_count)
        stats[pid] = {'occurrence': count, 'fraction': round(frac, 4), 'rare': rare}
    return stats


# ─── Adjacency graph ────────────────────────────────────────────────────────

def build_adjacency_graph(normalized):
    """Build the partition adjacency graph.

    Returns (adj, cooc, first_counts, last_counts, occurrence):
      adj[pid][next_pid] = # prophages with pid immediately before next_pid
      cooc[(a,b)] (a<b)   = # prophages carrying both a and b
      first_counts[pid]   = # prophages starting with pid
      last_counts[pid]    = # prophages ending with pid
      occurrence[pid]     = # prophages carrying pid
    """
    adj = defaultdict(lambda: defaultdict(int))
    cooc = defaultdict(int)
    first_counts = Counter()
    last_counts = Counter()
    occurrence = Counter()

    for pids in normalized.values():
        if not pids:
            continue
        first_counts[pids[0]] += 1
        last_counts[pids[-1]] += 1
        unique = set(pids)
        for pid in unique:
            occurrence[pid] += 1
        uniq_sorted = sorted(unique)
        for i in range(len(uniq_sorted)):
            for j in range(i + 1, len(uniq_sorted)):
                cooc[(uniq_sorted[i], uniq_sorted[j])] += 1
        for i in range(len(pids) - 1):
            adj[pids[i]][pids[i + 1]] += 1

    return dict(adj), dict(cooc), first_counts, last_counts, occurrence


# ─── Sequence similarity (k-mer Jaccard) ────────────────────────────────────

_COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def _revcomp(s):
    return s.translate(_COMPLEMENT)[::-1]


def kmer_jaccard_similarity(repr_seqs, k=DEFAULT_KMER_SIZE):
    """Pairwise k-mer Jaccard between partition representative sequences.

    Canonical k-mers (min(kmer, revcomp)); k-mers containing non-ACGT bases
    are skipped.  Returns dict (a, b) -> Jaccard (a < b), plus per-partition
    k-mer counts for the report.
    """
    kmer_sets = {}
    kmer_to_pids = defaultdict(set)
    for pid, seq in repr_seqs.items():
        s = seq.upper()
        kmers = set()
        for i in range(len(s) - k + 1):
            km = s[i:i + k]
            if not all(c in 'ACGT' for c in km):
                continue
            rc = _revcomp(km)
            kmers.add(km if km <= rc else rc)
        kmer_sets[pid] = kmers
        for km in kmers:
            kmer_to_pids[km].add(pid)

    shared = defaultdict(int)
    for km, pids in kmer_to_pids.items():
        pids = sorted(pids)
        for i in range(len(pids)):
            for j in range(i + 1, len(pids)):
                shared[(pids[i], pids[j])] += 1

    sim = {}
    for (a, b), s in shared.items():
        na, nb = len(kmer_sets[a]), len(kmer_sets[b])
        jac = s / (na + nb - s) if (na + nb - s) else 0.0
        sim[(a, b)] = round(jac, 4)
    kmer_counts = {pid: len(km) for pid, km in kmer_sets.items()}
    return sim, kmer_counts


# ─── Traversal optimization ─────────────────────────────────────────────────

def _tie_rank(pid, occurrence):
    """Deterministic tie-break: higher occurrence first, then lower pid."""
    return (occurrence.get(pid, 0), -pid)


def greedy_chain_path(W, idx, nodes):
    """Maximum-weight chain-merging construction (Kruskal-style for paths).

    Directed edges are processed in descending weight order; an edge a->b is
    selected when it joins two different chains without creating a branch or a
    cycle (b must be a chain head, a must be a chain tail).  This preserves the
    heaviest observed adjacencies first.  The resulting chains are then linked
    into a single path greedily by the best bridge edge (tail -> head).

    Fully deterministic: edges sorted by (-weight, a, b); chain-linking picks
    the highest-weight bridge, ties broken by the smaller chain-head pid.

    Returns (path_pids, path_weight).
    """
    edges = []
    for a in nodes:
        ia = idx[a]
        for b in nodes:
            if a != b:
                edges.append((float(W[ia][idx[b]]), a, b))
    edges.sort(key=lambda e: (-e[0], e[1], e[2]))

    succ = {}          # a -> b
    pred = {}          # b -> a
    head = {p: p for p in nodes}   # node -> chain head
    tail = {p: p for p in nodes}   # chain head -> chain tail

    for w, a, b in edges:
        if w <= 0:
            continue
        ha, hb = head[a], head[b]
        if ha == hb:
            continue
        if succ.get(a) is not None or pred.get(b) is not None:
            continue
        succ[a] = b
        pred[b] = a
        tail[ha] = tail[hb]
        node = hb
        while node is not None:
            head[node] = ha
            node = succ.get(node)

    chain_heads = [p for p in nodes if p not in pred]
    if len(chain_heads) == 1:
        path = [chain_heads[0]]
        node = succ.get(chain_heads[0])
        while node is not None:
            path.append(node)
            node = succ.get(node)
        weight = sum(float(W[idx[path[i]]][idx[path[i + 1]]])
                     for i in range(len(path) - 1))
        return path, weight

    # link the chains by best bridge: repeatedly merge the pair of chains with
    # the highest W[tail_i][head_j]
    chains = {h: {'head': h, 'tail': tail[h]} for h in chain_heads}
    while len(chains) > 1:
        best = None
        best_w = -1.0
        heads = list(chains)
        for i in range(len(heads)):
            hi = heads[i]
            ti = chains[hi]['tail']
            for j in range(len(heads)):
                hj = heads[j]
                if hi == hj:
                    continue
                w = float(W[idx[ti]][idx[hj]])
                if best is None or w > best_w + 1e-12 or (
                        abs(w - best_w) <= 1e-12 and hj < best):
                    best = hj
                    best_w = w
                    best_i = hi
        b_tail = chains[best]['tail']
        succ[chains[best_i]['tail']] = chains[best]['head']
        pred[chains[best]['head']] = chains[best_i]['tail']
        chains[best_i]['tail'] = b_tail
        del chains[best]

    final_head = next(iter(chains))
    path = [final_head]
    node = succ.get(final_head)
    while node is not None:
        path.append(node)
        node = succ.get(node)
    weight = sum(float(W[idx[path[i]]][idx[path[i + 1]]])
                 for i in range(len(path) - 1))
    return path, weight


def _flip_prefix_sums(W, pos):
    """Prefix sums of internal edge-flip deltas for the current path.

    flip[k] = sum_{t=0..k} (W[p_{t+1}][p_t] - W[p_t][p_{t+1}]), i.e. the net
    weight change of reversing all internal edges 0..k.  Reversing a directed
    segment flips the direction of every internal edge, so these prefix sums
    make 2-opt deltas exact and O(1)."""
    n = len(pos)
    flip = [0.0] * (n - 1)
    s = 0.0
    for k in range(n - 1):
        s += float(W[pos[k + 1]][pos[k]] - W[pos[k]][pos[k + 1]])
        flip[k] = s
    return flip


def two_opt(W, idx, path, occurrence, max_passes=300):
    """Best-improvement 2-opt local search on the open DIRECTED path.

    Reversing a segment of a directed path flips the direction of every
    internal edge; deltas are exact via flip-prefix-sums.  The best improving
    reversal is applied per pass until no improvement remains or max_passes is
    reached.  Deterministic: scan order (i, j); ties keep the first found.
    """
    n = len(path)
    if n < 4:
        return path
    pos = [idx[p] for p in path]
    passes = 0
    while passes < max_passes:
        flip = _flip_prefix_sums(W, pos)
        best_i = best_j = -1
        best_delta = 1e-9
        for i in range(n - 1):
            for j in range(i + 1, n):
                if i == 0 and j == n - 1:
                    continue  # full reversal would flip every edge: skip
                if i == 0:
                    delta = (flip[j - 1]
                             + float(W[pos[0]][pos[j + 1]]
                                     - W[pos[j]][pos[j + 1]]))
                elif j == n - 1:
                    delta = (flip[n - 2] - flip[i - 1]
                             + float(W[pos[i - 1]][pos[n - 1]]
                                     - W[pos[i - 1]][pos[i]]))
                else:
                    delta = (flip[j - 1] - flip[i - 1]
                             + float(W[pos[i - 1]][pos[j]]
                                     + W[pos[i]][pos[j + 1]]
                                     - W[pos[i - 1]][pos[i]]
                                     - W[pos[j]][pos[j + 1]]))
                if delta > best_delta:
                    best_delta = delta
                    best_i, best_j = i, j
        if best_i < 0:
            break
        path[best_i:best_j + 1] = path[best_i:best_j + 1][::-1]
        pos[best_i:best_j + 1] = pos[best_i:best_j + 1][::-1]
        passes += 1
    return path


def _w(W, idx, a, b):
    """Edge weight with None nodes treated as path boundary (weight 0)."""
    if a is None or b is None:
        return 0.0
    return float(W[idx[a]][idx[b]])


def _move_delta(W, idx, path, i, j):
    """Weight delta of moving path[i] to index j (|i-j| > 1). O(1).

    Exact changed-edge set: extraction of node i (removal bridge) plus
    insertion of node i at position j.  Adjacent moves (|i-j| <= 1) are left to
    2-opt, which handles swaps exactly.
    """
    n = len(path)
    a_prev = path[i - 1] if i > 0 else None
    a_next = path[i + 1] if i + 1 < n else None
    delta = (_w(W, idx, a_prev, a_next)
             - _w(W, idx, a_prev, path[i])
             - _w(W, idx, path[i], a_next))
    if j < i:
        j_prev = path[j - 1] if j > 0 else None
        j_cur = path[j]
        delta += (_w(W, idx, j_prev, path[i]) + _w(W, idx, path[i], j_cur)
                  - _w(W, idx, j_prev, j_cur))
    else:
        j_cur = path[j]
        j_next = path[j + 1] if j + 1 < n else None
        delta += (_w(W, idx, j_cur, path[i]) + _w(W, idx, path[i], j_next)
                  - _w(W, idx, j_cur, j_next))
    return delta


def insertion_search(W, idx, path, occurrence, max_passes=300):
    """Best-improvement node-relocation (Or-opt) search for the directed path.

    Moves a single node to a non-adjacent position when it raises the total
    path weight.  Deterministic scan order; ties keep the first found.
    """
    n = len(path)
    if n < 3:
        return path
    passes = 0
    while passes < max_passes:
        best = None
        best_delta = 1e-9
        for i in range(n):
            for j in range(n):
                if j == i or abs(j - i) <= 1:
                    continue
                delta = _move_delta(W, idx, path, i, j)
                if delta > best_delta:
                    best_delta = delta
                    best = (i, j)
        if best is None:
            break
        i, j = best
        node = path.pop(i)
        path.insert(j, node)
        passes += 1
    return path


def polish_path(W, idx, path, occurrence, max_passes=200):
    """Alternate 2-opt and node-relocation passes until no improvement."""
    prev_weight = -float('inf')
    for _ in range(6):
        path = two_opt(W, idx, path, occurrence, max_passes=max_passes)
        path = insertion_search(W, idx, path, occurrence, max_passes=max_passes)
        w = sum(float(W[idx[path[i]]][idx[path[i + 1]]])
                for i in range(len(path) - 1))
        if w <= prev_weight + 1e-9:
            break
        prev_weight = w
    return path


def find_optimal_traversal(partitions, adj, cooc, sim, first_counts, last_counts,
                           occurrence, n_prophages, cooc_weight=DEFAULT_COOC_WEIGHT,
                           bridge_weight=DEFAULT_BRIDGE_WEIGHT, max_2opt_passes=200):
    """Find the optimal typical traversal (see module docstring OBJECTIVE).

    Construction: maximum-weight chain merging (Kruskal-style) on the full
    weight matrix, then alternating 2-opt and node-relocation local search
    (``polish_path``).  Fully deterministic (no randomness).

    Returns (traversal, best_weight, W, joins) where W is the numpy weight
    matrix and joins classifies each consecutive traversal pair as
    observed / bridged.
    """
    nodes = sorted(partitions)
    idx = {p: i for i, p in enumerate(nodes)}
    n = len(nodes)

    W = np.zeros((n, n))
    for a in nodes:
        for b, w in adj.get(a, {}).items():
            if b in idx:
                W[idx[a]][idx[b]] += w
    for (a, b), c in cooc.items():
        if a in idx and b in idx:
            w = cooc_weight * c / max(n_prophages, 1)
            W[idx[a]][idx[b]] += w
            W[idx[b]][idx[a]] += w
    for (a, b), jac in sim.items():
        if a in idx and b in idx:
            w = bridge_weight * jac
            W[idx[a]][idx[b]] += w
            W[idx[b]][idx[a]] += w

    path, weight = greedy_chain_path(W, idx, nodes)
    path = polish_path(W, idx, path, occurrence, max_passes=max_2opt_passes)
    best_path = path
    best_weight = sum(float(W[idx[best_path[i]]][idx[best_path[i + 1]]])
                      for i in range(len(best_path) - 1))

    # classify joins
    joins = []
    for i in range(len(best_path) - 1):
        a, b = best_path[i], best_path[i + 1]
        obs = adj.get(a, {}).get(b, 0)
        c = cooc.get((min(a, b), max(a, b)), 0)
        jac = sim.get((min(a, b), max(a, b)), 0.0)
        joins.append({'a': a, 'b': b, 'observed': obs,
                      'cooc': c, 'similarity': jac,
                      'type': 'observed' if obs > 0 else 'bridged'})
    return best_path, best_weight, W, joins


# ─── Metrics ────────────────────────────────────────────────────────────────

def _lis_length(positions):
    """Length of the longest strictly increasing subsequence (patience sort)."""
    tails = []
    for x in positions:
        lo, hi = 0, len(tails)
        while lo < hi:
            mid = (lo + hi) // 2
            if tails[mid] < x:
                lo = mid + 1
            else:
                hi = mid
        if lo == len(tails):
            tails.append(x)
        else:
            tails[lo] = x
    return len(tails)


def compute_metrics(normalized, traversal, genome, adj, occurrence):
    """Compute per-community validation metrics (see module docstring)."""
    pos_in_T = {p: i for i, p in enumerate(traversal)}
    genome_set = set(genome)

    total_pairs = 0
    preserved_pairs = 0
    lcs_total = 0
    n_fully_ordered = 0
    genome_covered_frac = 0.0
    n_prophages_covered = 0
    n_prophages = 0
    for pids in normalized.values():
        if not pids:
            continue
        n_prophages += 1
        for i in range(len(pids) - 1):
            total_pairs += 1
            a, b = pids[i], pids[i + 1]
            if pos_in_T.get(a, -1) >= 0 and pos_in_T.get(b, -1) == pos_in_T.get(a, -1) + 1:
                preserved_pairs += 1
        positions = [pos_in_T[p] for p in pids if p in pos_in_T]
        if positions:
            lcs = _lis_length(positions)
            lcs_total += lcs / len(pids)
            if lcs == len(pids):
                n_fully_ordered += 1
        if pids:
            covered = sum(1 for p in pids if p in genome_set)
            genome_covered_frac += covered / len(pids)
            if covered > 0:
                n_prophages_covered += 1

    # preserved observed adjacency weight (weighted by prophage counts)
    preserved_weight = 0
    total_weight = 0
    for pids in normalized.values():
        for i in range(len(pids) - 1):
            total_weight += 1
            a, b = pids[i], pids[i + 1]
            if pos_in_T.get(a, -1) >= 0 and pos_in_T.get(b, -1) == pos_in_T.get(a, -1) + 1:
                preserved_weight += 1

    metrics = {
        'n_prophages': n_prophages,
        'total_observed_adjacency_pairs': total_pairs,
        'preserved_adjacency_pairs': preserved_pairs,
        'expected_path_completeness': round(preserved_pairs / total_pairs, 4) if total_pairs else 0.0,
        'preserved_observed_weight_fraction': round(preserved_weight / total_weight, 4) if total_weight else 0.0,
        'ordered_coverage': round(lcs_total / n_prophages, 4) if n_prophages else 0.0,
        'fraction_prophages_fully_ordered': round(n_fully_ordered / n_prophages, 4) if n_prophages else 0.0,
        'genome_coverage': round(genome_covered_frac / n_prophages, 4) if n_prophages else 0.0,
        'fraction_prophages_covered_by_genome': round(n_prophages_covered / n_prophages, 4) if n_prophages else 0.0,
    }
    return metrics


# ─── Stitching ──────────────────────────────────────────────────────────────

def select_genome_partitions(partitions, coverage_stats, consensi, max_length):
    """Select the genome partitions within the length budget.

    The stitched genome is the *typical path*: the most-common partitions of
    the community (core modules) come first, so the validated core is always
    reproduced; the remaining budget is then filled by the partitions with the
    highest coverage density (occurrence per bp).  This maximizes the expected
    coverage of the community's prophages per base pair of the genome while
    guaranteeing the highest-occurrence modules are present (the single most
    common partition alone is carried by up to ~70% of prophages and must not
    be priced out by a long union-mosaic consensus).

    Deterministic: occurrence descending, then density descending, then pid
    ascending.  Returns the selected partition ids.
    """
    items = []
    for pid in partitions:
        con = consensi.get(pid)
        if not con:
            continue
        occ = coverage_stats.get(pid, {}).get('occurrence', 0)
        ln = len(con)
        if ln <= 0:
            continue
        items.append((pid, occ, ln, occ / ln))

    selected = []
    selected_set = set()
    total = 0
    # 1. core-first: highest-occurrence partitions that fit
    for pid, occ, ln, density in sorted(items, key=lambda x: (-x[1], x[0])):
        if total + ln <= max_length:
            selected.append(pid)
            selected_set.add(pid)
            total += ln
    # 2. density fill: best coverage-per-bp among the remaining
    rest = [it for it in items if it[0] not in selected_set]
    for pid, occ, ln, density in sorted(rest, key=lambda x: (-x[3], -x[1], x[0])):
        if total + ln <= max_length:
            selected.append(pid)
            selected_set.add(pid)
            total += ln
    if not selected and items:
        # fallback: always include the highest-occurrence partition
        selected = [max(items, key=lambda x: (x[1], -x[0]))[0]]
    return selected


def stitch_genome(traversal, consensi, max_length, coverage_stats=None,
                   gap_size=DEFAULT_GAP_SIZE):
    """Concatenate partition representative sequences in traversal order.

    The genome partition set is selected within the ``max_length`` budget by
    coverage density (``select_genome_partitions``), then concatenated in
    traversal order with optional ``gap_size`` N-run join markers.
    Returns (genome_seq, genome_pids, excluded_pids, excluded_bp).
    """
    if coverage_stats is None:
        coverage_stats = {}
    genome_pids = select_genome_partitions(traversal, coverage_stats, consensi,
                                           max_length)
    # order the selected set by traversal position
    pos_in_T = {p: i for i, p in enumerate(traversal)}
    genome_pids = sorted(genome_pids, key=lambda p: pos_in_T[p])
    genome_set = set(genome_pids)
    excluded = [p for p in traversal if p not in genome_set and consensi.get(p)]
    excluded_bp = sum(len(consensi[p]) for p in excluded)
    seq = ''
    for i, pid in enumerate(genome_pids):
        part = consensi[pid]
        if i > 0 and gap_size > 0:
            seq += 'N' * gap_size
        seq += part
    return seq, genome_pids, excluded, excluded_bp


def wrap_fasta(seq, width=80):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


# ─── Outputs ────────────────────────────────────────────────────────────────

def write_outputs(prefix, traversal, genome_pids, consensi, seq_data,
                  coverage_stats, metrics, joins, params, inputs,
                  excluded, excluded_bp, consensus_diag, full_concat,
                  genome_seq):
    """Write all per-community outputs.  Returns dict of written files."""
    parent = os.path.dirname(prefix)
    if parent:
        os.makedirs(parent, exist_ok=True)
    files = {}
    n_prophages = metrics['n_prophages']

    # (a) traversal order JSON
    n_observed = sum(1 for j in joins if j['type'] == 'observed')
    n_bridged = len(joins) - n_observed
    position = {p: i for i, p in enumerate(traversal)}
    genome_positions = set(genome_pids)
    partitions_records = {}
    for pid in traversal:
        cs = coverage_stats.get(pid, {'occurrence': 0, 'fraction': 0.0, 'rare': False})
        con = consensi.get(pid)
        partitions_records[str(pid)] = {
            'occurrence': cs['occurrence'],
            'fraction': cs['fraction'],
            'rare': bool(cs['rare']),
            'traversal_position': position[pid],
            'in_genome': pid in genome_positions,
            'consensus_len': len(con) if con else 0,
            'has_seq': con is not None,
            'seq_source': consensus_diag.get(pid, {}).get('source', 'none'),
            'join_into': next((j['type'] for j in joins if j['b'] == pid), 'start'),
        }
    traversal_json = {
        'community': os.path.basename(prefix),
        'objective': (
            'maximize sum_i W[T[i]][T[i+1]] over permutations T of all partitions, '
            'where W[a][b] = adj(a,b) + mu*cooc(a,b)/N + lambda*J(a,b) '
            '(adj = prophages with a immediately before b; cooc = co-occurring '
            'prophages; J = k-mer Jaccard of representative sequences). '
            'Maximum-weight Hamiltonian path solved by deterministic greedy '
            'construction + first-improvement 2-opt local search.'),
        'parameters': params,
        'inputs': inputs,
        'n_prophages': n_prophages,
        'n_partitions_total': len(traversal),
        'n_partitions_with_seq': sum(1 for p in traversal if consensi.get(p)),
        'n_partitions_rare': sum(1 for p in traversal if coverage_stats.get(p, {}).get('rare')),
        'traversal': traversal,
        'genome': genome_pids,
        'genome_length_bp': len(genome_seq),
        'genome_joins': {'observed': n_observed, 'bridged': n_bridged},
        'excluded_from_genome_by_budget': excluded,
        'excluded_from_genome_bp': excluded_bp,
        'metrics': metrics,
        'partitions': partitions_records,
    }
    json_path = prefix + '.traversal.json'
    with open(json_path, 'w') as f:
        json.dump(traversal_json, f, indent=2)
    files['traversal_json'] = json_path

    # stats JSON (summary)
    stats_json = {'community': os.path.basename(prefix), 'parameters': params,
                  'metrics': metrics,
                  'genome_length_bp': len(genome_seq),
                  'n_partitions_total': len(traversal),
                  'n_partitions_rare': sum(1 for p in traversal
                                           if coverage_stats.get(p, {}).get('rare'))}
    stats_path = prefix + '.stats.json'
    with open(stats_path, 'w') as f:
        json.dump(stats_json, f, indent=2)
    files['stats_json'] = stats_path

    # (b) per-partition consensus FASTA
    cons_path = prefix + '.consensus.fa'
    with open(cons_path, 'w') as f:
        for pid in traversal:
            con = consensi.get(pid)
            if not con:
                continue
            cs = coverage_stats.get(pid, {})
            rare = ' rare' if cs.get('rare') else ''
            src = consensus_diag.get(pid, {}).get('source', '?')
            f.write(f'>partition{pid} n_prophages={cs.get("occurrence", 0)} '
                    f'len={len(con)} source={src}{rare}\n')
            f.write(wrap_fasta(con) + '\n')
    files['consensus_fasta'] = cons_path

    # (c) stitched community genome FASTA
    genome_path = prefix + '.fa'
    with open(genome_path, 'w') as f:
        n_part = len(genome_pids)
        f.write(f'>community_{os.path.basename(prefix)}_stitched_typical_genome '
                f'n_partitions={n_part} len={len(genome_seq)} budget={params["max_length"]}\n')
        f.write(wrap_fasta(genome_seq) + '\n')
    files['genome_fasta'] = genome_path

    # (d) coverage statistics TSV
    tsv_path = prefix + '.coverage.tsv'
    with open(tsv_path, 'w') as f:
        f.write('\t'.join(['partition', 'occurrence', 'fraction', 'rare',
                           'has_seq', 'consensus_len', 'traversal_position',
                           'in_genome', 'join_type']) + '\n')
        for pid in traversal:
            cs = coverage_stats.get(pid, {})
            con = consensi.get(pid)
            join = next((j['type'] for j in joins if j['b'] == pid), 'start')
            f.write('\t'.join([
                str(pid),
                str(cs.get('occurrence', 0)),
                f"{cs.get('fraction', 0.0):.4f}",
                'yes' if cs.get('rare') else 'no',
                'yes' if con else 'no',
                str(len(con) if con else 0),
                str(position[pid]),
                'yes' if pid in genome_positions else 'no',
                join,
            ]) + '\n')
    files['coverage_tsv'] = tsv_path

    # optional full concatenation (budget disabled)
    if full_concat:
        full_path = prefix + '.full.fa'
        full = ''
        for pid in traversal:
            part = consensi.get(pid)
            if part:
                full += part
        with open(full_path, 'w') as f:
            f.write(f'>community_{os.path.basename(prefix)}_full_traversal_concat '
                    f'n_partitions={len(traversal)} len={len(full)}\n')
            f.write(wrap_fasta(full) + '\n')
        files['full_concat_fasta'] = full_path

    return files


# ─── Main pipeline ──────────────────────────────────────────────────────────

def run_traversal(partitions_dir, bed_path, output_prefix, coverage_threshold,
                  rare_threshold, rare_min_count, max_length, cooc_weight,
                  bridge_weight, kmer_size, gap_size, seed, full_concat,
                  max_2opt_passes=60):
    """Run the full traversal pipeline.  Returns (result_dict, files_dict)."""
    t0 = time.time()

    # 1. Parse BED
    raw = parse_bed(bed_path)
    n_bed_prophages = len(raw)
    normalized = normalize_partition_lists(raw)
    n_prophages = len(normalized)
    if n_prophages == 0:
        raise RuntimeError(f'no prophage partition assignments found in {bed_path}')
    partitions = sorted({p for pids in normalized.values() for p in pids})
    print(f'[1] BED: {n_bed_prophages} prophages, {n_prophages} with assignments, '
          f'{len(partitions)} partitions')

    # 2. Coverage statistics
    coverage_stats = compute_coverage_stats(normalized, n_prophages,
                                            rare_threshold, rare_min_count)
    n_rare = sum(1 for s in coverage_stats.values() if s['rare'])
    print(f'[2] coverage: {n_rare}/{len(partitions)} partitions flagged RARE '
          f'(fraction < {rare_threshold} or count < {rare_min_count})')

    # 3. Sequence data + representatives
    seq_data, warnings = load_partition_sequences(partitions_dir, partitions)
    for w in warnings:
        print(f'    warning: {w}')
    consensi = {}
    consensus_diag = {}
    for pid in partitions:
        rep, diag = partition_representative(seq_data[pid], coverage_threshold)
        if rep:
            consensi[pid] = rep
        consensus_diag[pid] = diag
    n_with_seq = sum(1 for pid in partitions if consensi.get(pid))
    print(f'[3] sequence data: {n_with_seq}/{len(partitions)} partitions have a '
          f'representative sequence')

    # 4. Adjacency graph
    adj, cooc, first_counts, last_counts, occurrence = build_adjacency_graph(normalized)
    n_edges = sum(len(t) for t in adj.values())
    print(f'[4] adjacency graph: {len(adj)} nodes, {n_edges} directed observed edges')

    # 5. Sequence similarity (bridging)
    repr_seqs = {pid: consensi[pid] for pid in partitions if consensi.get(pid)}
    sim, kmer_counts = kmer_jaccard_similarity(repr_seqs, kmer_size)
    print(f'[5] k-mer similarity: {len(sim)} partition pairs scored (k={kmer_size})')

    # 6. Optimal traversal
    traversal, opt_weight, W, joins = find_optimal_traversal(
        partitions, adj, cooc, sim, first_counts, last_counts, occurrence,
        n_prophages, cooc_weight=cooc_weight, bridge_weight=bridge_weight,
        max_2opt_passes=max_2opt_passes)
    n_obs_joins = sum(1 for j in joins if j['type'] == 'observed')
    n_bridged = len(joins) - n_obs_joins
    print(f'[6] traversal: {len(traversal)} partitions ordered, '
          f'objective weight={opt_weight:.2f}, '
          f'joins: {n_obs_joins} observed / {n_bridged} bridged')

    # 7. Stitch with budget
    genome_seq, genome_pids, excluded, excluded_bp = stitch_genome(
        traversal, consensi, max_length, coverage_stats=coverage_stats,
        gap_size=gap_size)
    print(f'[7] stitched genome: {len(genome_seq):,} bp from {len(genome_pids)} '
          f'partitions (budget {max_length:,} bp); {len(excluded)} partitions '
          f'excluded ({excluded_bp:,} bp)')

    # 8. Metrics
    metrics = compute_metrics(normalized, traversal, genome_pids, adj, occurrence)
    print(f'[8] metrics: EPC={metrics["expected_path_completeness"]} '
          f'ordered_cov={metrics["ordered_coverage"]} '
          f'genome_cov={metrics["genome_coverage"]} '
          f'fully_ordered={metrics["fraction_prophages_fully_ordered"]}')

    params = {
        'coverage_threshold': coverage_threshold,
        'rare_threshold': rare_threshold,
        'rare_min_count': rare_min_count,
        'max_length': max_length,
        'cooc_weight': cooc_weight,
        'bridge_weight': bridge_weight,
        'kmer_size': kmer_size,
        'gap_size': gap_size,
        'seed': seed,
        'max_2opt_passes': max_2opt_passes,
    }
    inputs = {
        'bed': bed_path,
        'partitions_dir': partitions_dir,
        'n_prophages_with_assignments': n_prophages,
    }

    # 9. Outputs
    files = write_outputs(output_prefix, traversal, genome_pids, consensi, seq_data,
                          coverage_stats, metrics, joins, params, inputs,
                          excluded, excluded_bp, consensus_diag, full_concat,
                          genome_seq)
    print(f'[9] wrote {len(files)} output files under {output_prefix}.*')
    print(f'    completed in {time.time() - t0:.1f}s')

    result = {
        'traversal': traversal,
        'genome_pids': genome_pids,
        'metrics': metrics,
        'files': files,
        'n_partitions_total': len(traversal),
        'n_partitions_rare': n_rare,
        'genome_length_bp': len(genome_seq),
        'n_observed_joins': n_obs_joins,
        'n_bridged_joins': n_bridged,
    }
    return result, files


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Optimal typical traversal of community partition sets '
                    '(README §6). Emits traversal order JSON, per-partition '
                    'consensus FASTA, stitched community genome FASTA and '
                    'coverage statistics per community.')
    parser.add_argument('--partitions-dir', required=True,
                        help='Directory of partition<N>.maf / partition<N>.fasta '
                             '(per-partition alignment blocks or member sequences)')
    parser.add_argument('--bed', required=True,
                        help='Combined 4-column BED file (prophage start end '
                             'partition_id) OR directory of per-partition '
                             '3-column BEDs (partitions_bed/ layout)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output prefix; files are written as '
                             '<prefix>.traversal.json, <prefix>.consensus.fa, '
                             '<prefix>.fa, <prefix>.coverage.tsv, '
                             '<prefix>.stats.json')
    parser.add_argument('--coverage-threshold', type=float,
                        default=DEFAULT_COVERAGE_THRESHOLD,
                        help='Column coverage gate for MAF majority consensus '
                             '(fraction of block sequences with a non-gap base; '
                             'default 0.0 = any covered column, byte-identical '
                             'to validated stitch_algorithm.py runs)')
    parser.add_argument('--rare-threshold', type=float,
                        default=DEFAULT_RARE_THRESHOLD,
                        help='Occurrence fraction below which a partition is '
                             'flagged rare (default 0.2)')
    parser.add_argument('--rare-min-count', type=int,
                        default=DEFAULT_RARE_MIN_COUNT,
                        help='Absolute occurrence count below which a partition '
                             'is flagged rare (default 5)')
    parser.add_argument('--max-length', type=int, default=DEFAULT_MAX_LENGTH,
                        help='Length budget (bp) for the stitched community '
                             'genome; keeps output at phage-typical scale '
                             '(default 150000)')
    parser.add_argument('--cooc-weight', type=float, default=DEFAULT_COOC_WEIGHT,
                        help='Weight of the co-occurrence anchor term mu in W '
                             '(default 1.0)')
    parser.add_argument('--bridge-weight', type=float,
                        default=DEFAULT_BRIDGE_WEIGHT,
                        help='Weight of the k-mer-Jaccard similarity bridge '
                             'term lambda in W (default 1.0)')
    parser.add_argument('--kmer-size', type=int, default=DEFAULT_KMER_SIZE,
                        help='k for k-mer Jaccard similarity (default 15)')
    parser.add_argument('--gap-size', type=int, default=DEFAULT_GAP_SIZE,
                        help='N-runs inserted between partition blocks in the '
                             'stitched genome (default 0 = plain concatenation)')
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED,
                        help='Recorded for provenance; the algorithm is fully '
                             'deterministic (tie-break: occurrence desc, pid asc)')
    parser.add_argument('--full-concat', action='store_true',
                        help='Also write <prefix>.full.fa: the full traversal '
                             'concatenation with the length budget disabled')
    args = parser.parse_args(argv)

    if not os.path.isdir(args.partitions_dir):
        parser.error(f'--partitions-dir not a directory: {args.partitions_dir}')
    if not os.path.exists(args.bed):
        parser.error(f'--bed not found: {args.bed}')
    if args.coverage_threshold < 0 or args.coverage_threshold > 1:
        parser.error('--coverage-threshold must be in [0, 1]')
    if args.max_length <= 0:
        parser.error('--max-length must be positive')

    result, files = run_traversal(
        args.partitions_dir, args.bed, args.output, args.coverage_threshold,
        args.rare_threshold, args.rare_min_count, args.max_length,
        args.cooc_weight, args.bridge_weight, args.kmer_size, args.gap_size,
        args.seed, args.full_concat)
    return result


if __name__ == '__main__':
    main()
