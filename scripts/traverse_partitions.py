#!/usr/bin/env python3
"""
traverse_partitions.py — two-level ML inference: typical (ML) and ancestral
phage genomes per community partition set
==============================================================================

Given the partition output of ``impg partition`` for a prophage community —
(a) a BED mapping prophages to ordered partition intervals and (b) per-partition
alignment blocks (MAF) or member-sequence bundles (FASTA) — this script builds,
for every community:

  * a set of **sampled traversals** of the partition adjacency graph
    (Level 1: ML *across* partitions, by weighted path sampling), and
  * per-partition representative sequences (Level 2: ML *within* each
    partition), in one of two modes:
      - ``--mode ml``        : majority-rule consensus of each aligned block,
      - ``--mode ancestral`` : per-partition tree (neighbor-joining) +
                               Fitch-parsimony ancestral-state reconstruction.
  * the **stitched community genome** per mode: the partition representatives
    concatenated in the order of the highest-support sampled traversal
    (the maximum-likelihood path among the draws), capped by a length budget
    so the output stays at phage-typical scale (tens of kb to ~150 kb) instead
    of the multi-Mbp naive concatenations of all blocks (community 3 alone was
    ~1.29 Mbp).

LEVEL 1 — ML ACROSS PARTITIONS (weighted path sampling)
-------------------------------------------------------
Each partition ``p`` carried by ``occ_p`` of the community's ``N`` prophages
gets a sampling weight

    w_p = (occ_p / N)^alpha                     (alpha, default 1.0)

Every partition with occ_p >= 1 has w_p > 0, so NOTHING is excluded by a hard
threshold and nothing is forced in deterministically: rare partitions are
sampled exactly like common ones, only less likely.

A traversal is sampled as a random walk over the partition set (seeded RNG):

  * the first partition is drawn with probability proportional to w_p;
  * at each step, the next partition ``b`` is drawn from the unvisited
    partitions whose representative fits the remaining genome budget with

        P(b)  proportional to  w_b * (1 + beta * adj[current][b])

    where ``adj[a][b]`` = number of prophages carrying ``a`` immediately
    before ``b`` (observed adjacency support; ``beta``, default 1.0).  The
    ``(1 + beta*adj)`` factor makes the *common paths* (partitions that occur
    in many prophages *and* are observed adjacent) dramatically more likely,
    while the ``w_b`` factor keeps every unvisited partition reachable in
    proportion to its support.

  * the walk stops when no unvisited partition fits the remaining budget
    (``--max-length``) or all partitions have been visited.

``--n-samples`` traversals are drawn; per-partition **sampling frequency**
(share of runs in which the partition is on the sampled path) is reported —
common partitions are sampled in nearly every run, rare ones only
occasionally, at a rate consistent with their weight.  The *stitched genome*
is built from the sampled traversal with the highest path support
``sum(log2 w_p)`` (the ML path among the draws).

LEVEL 2 — ML WITHIN EACH PARTITION (``--mode ml``)
--------------------------------------------------
The representative of a partition is the **majority-rule consensus** of its
aligned block (MAF): each column independently chooses the most likely base
given the observed column frequencies (multinomial MLE = majority base), and
the column is emitted iff at least ``--coverage-threshold`` fraction of block
sequences carry a non-gap base there (default 0.25, the validated
stitch_algorithm.py default; for staggered impg segment bundles this emits
the true conserved locus — 0.0 emits the union mosaic, byte-identical to the
committed 2363ece consensus).  This is the documented "ML within" choice; the
alternative (closest-member-to-consensus) is not used because per-column
majority is the maximum-likelihood per-column estimate and is deterministic.

ANCESTRAL MODE (``--mode ancestral``)
-------------------------------------
Level 2 is replaced by ancestral-state reconstruction per partition:

  1. Neighbor-joining tree (Studier & Keppler 1988) of the aligned block
     members, from pairwise p-distances over columns where BOTH members carry
     a non-gap base (gap = missing data); pairs with no shared columns get
     distance 1.0; negative branch lengths are clamped to 0.  Ties in the NJ
     Q-matrix are broken by row-major order (deterministic).
  2. Unordered Fitch parsimony over {A,C,G,T} (gap = MISSING data: a gap
     leaf is unconstrained, contributing no state and no counts): an internal
     node's set is the intersection of its children's sets if nonempty, else
     their union.
  3. The reconstructed **root** (ancestral) state per column is the
     highest-frequency state in the root's set (ties broken by the fixed
     state order A < C < G < T).  The parsimony intersection rule constrains
     the root whenever the two root subtrees share exactly one possible
     state, so the result is not plain column majority.
  4. Columns are reconstructed only if at least ``--coverage-threshold``
     fraction of block sequences carry a non-gap base there (same column
     gate as the ML consensus, default 0.25); lower-coverage columns are
     insertions / alignment artifacts relative to the ancestor and are
     excluded.  This keeps the ancestral genome on the same scaffold as the
     ML genome (phage-typical after the length budget) instead of collapsing
     staggered impg bundles to tiny cores.

Assumptions: (i) members of a partition block are orthologous segments, so a
single gene tree approximates the locus tree; (ii) a gap in a member means it
does not carry that position (staggered segment bundles are handled exactly
this way); (iii) parsimony is appropriate for closely related prophage
segments (short divergence); (iv) p-distances with gap-missing are a
sufficient dissimilarity for the QC tree.

The per-partition NJ trees are emitted (Newick) for QC in ancestral mode.

DETERMINISM
-----------
All randomness is a single ``random.Random(seed)`` (default seed 42, recorded
in the output).  Same inputs + same parameters + same seed -> identical
samples and byte-identical outputs.

INPUT FORMATS (BED)
-------------------
* ``--bed FILE`` with 4 columns ``prophage start end partition_id`` (legacy
  combined layout, e.g. ``research/stitching/community_3_partitions.bed``), or
* ``--bed DIR`` containing per-partition 3-column BEDs ``partition<N>.bed``.
Both are normalized to per-prophage ordered partition lists (contiguous runs
merged, consecutive duplicates collapsed).

OUTPUTS (``--output PREFIX``)
-----------------------------
* ``PREFIX.traversal.json``      — all sampled traversals with per-partition
                                   weights, best (genome) sample, metrics.
* ``PREFIX.consensus.fa``        — per-partition ML representatives.
* ``PREFIX.ancestral.fa``        — per-partition ancestral representatives
                                   (ancestral mode only).
* ``PREFIX.ml.fa``               — stitched ML community genome (both modes).
* ``PREFIX.ancestral.genome.fa`` — stitched ancestral community genome
                                   (ancestral mode only).
* ``PREFIX.coverage.tsv``        — occurrence, fraction, weight, sampling
                                   frequency per partition.
* ``PREFIX.trees.nwk``           — per-partition NJ trees, one per line
                                   (ancestral mode only).
* ``PREFIX.stats.json``          — summary metrics.

USAGE
-----
    python3 scripts/traverse_partitions.py \
        --partitions-dir <dir> --bed <partitions.bed> --output <prefix> \
        [--mode ml|ancestral] [--n-samples N] [--alpha A] [--seed S]
"""

import argparse
import glob
import json
import math
import os
import random
import re
import time
from collections import Counter, defaultdict

import numpy as np

# ─── defaults ───────────────────────────────────────────────────────────────

DEFAULT_ALPHA = 1.0            # weight = occurrence_fraction ** alpha
DEFAULT_BETA = 1.0             # adjacency modulation in the step distribution
DEFAULT_N_SAMPLES = 10         # weighted traversal draws per community
DEFAULT_MAX_LENGTH = 150_000   # genome budget (phage-typical scale)
# Column coverage gate for MAF consensus AND ancestral reconstruction.
# 0.25 = the validated stitch_algorithm.py default: for staggered impg segment
# bundles it emits the true conserved locus (members overlap only partially),
# not the union mosaic (0.0 spans every segment position and can be 10-20x
# longer than any single member).  The union-mosaic behaviour is available via
# --coverage-threshold 0.
DEFAULT_COVERAGE_THRESHOLD = 0.25
DEFAULT_GAP_SIZE = 0           # N-runs between partition blocks in the genome
DEFAULT_SEED = 42              # seeded RNG -> reproducible samples
DEFAULT_MODE = 'ml'

STATES = 'ACGT-'                       # Fitch state order (tie-break order)
_STATE_IDX = {c: i for i, c in enumerate(STATES)}
_TRANSLATE_ACGT = str.maketrans('ACGTacgt', '01230123')
_NON_ACGT = re.compile(r'[^ACGTacgt]')
_NEWICK_SAFE = re.compile(r'[^\w|.-]')


# ─── BED parsing ────────────────────────────────────────────────────────────

def _parse_partition_id(path):
    """Extract the integer partition id from a partition<N>.<ext> path."""
    base = os.path.basename(path)
    m = re.search(r'(\d+)', base)
    return int(m.group(1)) if m else None


def _parse_bed_file(bed_path):
    """Parse a 4-column combined BED: prophage  start  end  partition_id."""
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

    The partition id is taken from the file name (partition<N>.bed)."""
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
    duplicates, producing each prophage's ordered unique partition sequence."""
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
        collapsed = []
        for p in (m[0] for m in merged):
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


# ─── Level 2 (ML within): consensus / representative sequences ─────────────

def compute_partition_consensus(records, coverage_threshold=DEFAULT_COVERAGE_THRESHOLD):
    """Majority-rule consensus over an aligned block (documented ML-within).

    A column is emitted iff at least ``max(1, ceil(coverage_threshold * n))``
    sequences carry a non-gap base there; the majority base wins (per-column
    multinomial MLE).  Handles ragged block widths by truncating to the
    minimum width.  Returns (consensus_string_or_None, diagnostics_dict)."""
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
    """ML representative for a partition: MAF majority consensus, else
    longest FASTA member sequence (fallback), else None.

    Returns (seq_or_None, diag)."""
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

def compute_coverage_stats(normalized, n_prophages):
    """Per-partition coverage statistics.

    Returns dict pid -> {'occurrence': count, 'fraction': count / n_prophages}.
    """
    occurrence = Counter()
    for pids in normalized.values():
        occurrence.update(set(pids))
    stats = {}
    for pid, count in occurrence.items():
        stats[pid] = {'occurrence': count,
                      'fraction': round(count / max(n_prophages, 1), 6)}
    return stats


# ─── Adjacency graph ────────────────────────────────────────────────────────

def build_adjacency_graph(normalized):
    """Build the partition adjacency graph from consecutive partition pairs.

    Returns (adj, first_counts, last_counts, occurrence):
      adj[pid][next_pid] = # prophages with pid immediately before next_pid
      first_counts[pid]  = # prophages starting with pid
      last_counts[pid]   = # prophages ending with pid
      occurrence[pid]    = # prophages carrying pid
    """
    adj = defaultdict(lambda: defaultdict(int))
    first_counts = Counter()
    last_counts = Counter()
    occurrence = Counter()

    for pids in normalized.values():
        if not pids:
            continue
        first_counts[pids[0]] += 1
        last_counts[pids[-1]] += 1
        for pid in set(pids):
            occurrence[pid] += 1
        for i in range(len(pids) - 1):
            adj[pids[i]][pids[i + 1]] += 1

    return dict(adj), first_counts, last_counts, occurrence


# ─── Level 1 (ML across partitions): weights + weighted path sampling ──────

def compute_sampling_weights(coverage_stats, alpha=DEFAULT_ALPHA):
    """Sampling weight per partition: w_p = occurrence_fraction^alpha.

    Every partition carried by >= 1 prophage has fraction > 0, hence w_p > 0:
    nothing is hard-excluded and nothing is forced in deterministically."""
    return {pid: stats['fraction'] ** alpha
            for pid, stats in coverage_stats.items()}


def sample_traversal(weights, adj, part_lengths, max_length, rng, beta=DEFAULT_BETA):
    """Draw ONE weighted traversal (a path) over the partition set.

    First partition proportional to w; each next partition b among the
    unvisited partitions that fit the remaining budget with
    P(b) ~ w_b * (1 + beta * adj[current][b]).  Stops when nothing fits (or
    all partitions visited).  No hard threshold on partition rarity.

    Returns (path, support, adj_pairs, genome_length_bp):
      support    = sum(log2 w_p) over the path (relative path support),
      adj_pairs  = number of consecutive path pairs that are observed
                   adjacencies in some prophage,
      genome_length_bp = sum of representative lengths on the path.
    """
    pids = list(weights)
    path = []
    visited = set()
    remaining = max_length
    total_len = 0
    current = None
    while True:
        cand = [p for p in pids
                if p not in visited and part_lengths.get(p, 0) <= remaining]
        if not cand:
            break
        if current is None:
            scores = [weights[p] for p in cand]
        else:
            adj_cur = adj.get(current, {})
            scores = [weights[p] * (1.0 + beta * adj_cur.get(p, 0.0))
                      for p in cand]
        pick = rng.choices(cand, weights=scores, k=1)[0]
        path.append(pick)
        visited.add(pick)
        ln = part_lengths.get(pick, 0)
        total_len += ln
        remaining -= ln
        current = pick
    support = sum(math.log2(max(weights[p], 1e-12)) for p in path)
    adj_pairs = sum(1 for i in range(len(path) - 1)
                    if adj.get(path[i], {}).get(path[i + 1], 0) > 0)
    return path, round(support, 4), adj_pairs, total_len


def sample_traversals(weights, adj, part_lengths, max_length, n_samples, seed,
                      beta=DEFAULT_BETA):
    """Draw n_samples traversals with a seeded RNG.

    Same seed -> identical samples.  Returns (samples, rng) where each sample
    is a dict {index, path, support, adj_pairs, genome_length_bp}."""
    rng = random.Random(seed)
    samples = []
    for s in range(n_samples):
        path, support, adj_pairs, length = sample_traversal(
            weights, adj, part_lengths, max_length, rng, beta)
        samples.append({'index': s, 'path': path, 'support': support,
                        'adj_pairs': adj_pairs, 'genome_length_bp': length})
    return samples, rng


# ─── Ancestral mode: NJ tree + Fitch parsimony ─────────────────────────────

def encode_alignment(records, width):
    """Encode aligned sequences as an (n, width) int8 array: 0-3 = A,C,G,T,
    4 = gap or any non-ACGT base."""
    rows = []
    for rec in records:
        s = _NON_ACGT.sub('4', rec['seq'][:width].upper())
        s = s.translate(_TRANSLATE_ACGT)
        rows.append(s)
    arr = np.frombuffer(''.join(rows).encode('ascii'),
                        dtype=np.int8).reshape(len(records), width)
    return arr - 48  # ASCII '0'-'4' -> 0-4


def compute_pairwise_distances(encoded):
    """Pairwise p-distances over the aligned block (gap = missing).

    distance(a,b) = mismatches / shared columns (columns where BOTH carry a
    non-gap base); pairs with no shared column get distance 1.0.  Identical
    columns are grouped (weighted by count) to keep the computation
    O(distinct_columns * n^2) — tractable on large staggered bundles.

    Returns (D, shared_counts): (n, n) float matrices."""
    n, _ = encoded.shape
    if n <= 1:
        return np.zeros((n, n)), np.zeros((n, n))
    patterns, counts = np.unique(encoded.T, axis=0, return_counts=True)
    patterns = patterns.astype(np.float64)          # (n_patterns, n)
    counts = counts.astype(np.float64)              # (n_patterns,)
    nongap = (patterns < 4).astype(np.float64)
    shared = nongap.T @ (nongap * counts[:, None])  # (n, n)
    eq = np.zeros((n, n))
    for v in range(4):
        ev = (patterns == v).astype(np.float64)
        eq += ev.T @ (ev * counts[:, None])
    mism = shared - eq
    with np.errstate(divide='ignore', invalid='ignore'):
        d = np.where(shared > 0, mism / np.maximum(shared, 1e-9), 1.0)
    np.fill_diagonal(d, 0.0)
    d = np.minimum(d, 1.0)
    return d, shared


def _sanitize_name(name):
    """Newick-safe leaf label."""
    return _NEWICK_SAFE.sub('_', str(name)) or 'node'


def neighbor_joining(D, names):
    """Neighbor-joining tree (Studier & Keppler 1988), deterministic.

    Ties for the minimum Q entry are broken by row-major order of the Q
    matrix (i.e. smaller node index first); negative branch lengths are
    clamped to 0.

    Returns (newick, children, root): children maps node id -> [child1,
    child2]; leaves are 0..n-1 in `names` order, internal nodes n..root."""
    n = len(names)
    if n == 1:
        return f'{_sanitize_name(names[0])};', {}, 0
    D = np.array(D, dtype=np.float64)
    active = list(range(n))
    children = {}
    branch = {}
    next_id = n
    while len(active) > 2:
        m = len(active)
        rsum = D.sum(axis=1)
        Q = (m - 2) * D - rsum[:, None] - rsum[None, :]
        np.fill_diagonal(Q, np.inf)
        flat = int(np.argmin(Q))
        i, j = divmod(flat, m)              # i < j by row-major argmin
        a, b = active[i], active[j]
        delta = (rsum[i] - rsum[j]) / (m - 2)
        li = max(0.0, (D[i, j] + delta) / 2.0)
        lj = max(0.0, (D[i, j] - delta) / 2.0)
        u = next_id
        next_id += 1
        children[u] = [a, b]
        branch[a] = li
        branch[b] = lj
        newrow = (D[i, :] + D[j, :] - D[i, j]) / 2.0
        rem = [k for k in range(m) if k not in (i, j)]
        D2 = np.zeros((m - 1, m - 1))
        D2[:m - 2, :m - 2] = D[np.ix_(rem, rem)]
        D2[m - 2, :m - 2] = newrow[rem]
        D2[:m - 2, m - 2] = newrow[rem]
        D = D2
        active = [active[k] for k in rem] + [u]
    a, b = active[0], active[1]
    root = next_id
    children[root] = [a, b]
    branch[a] = max(0.0, D[0, 1] / 2.0)
    branch[b] = max(0.0, D[0, 1] / 2.0)

    def fmt(node):
        if node < n:
            return f'{_sanitize_name(names[node])}:{branch.get(node, 0.0):.6f}'
        c1, c2 = children[node]
        return f'({fmt(c1)},{fmt(c2)}):0.000000'

    return fmt(root) + ';', children, root


def fitch_ancestral(encoded, children, root, n_leaves, coverage_threshold,
                   chunk=512):
    """Fitch-parsimony ancestral (root) sequence, deterministic.

    Gap (state 4) is treated as MISSING DATA: a gap leaf is unconstrained
    (its Fitch set is all four bases, contributing nothing to the state
    frequencies).  Columns are reconstructed only if at least
    ``coverage_threshold`` fraction of block sequences carry a non-gap base
    there (same column gate as the ML consensus); columns below the gate are
    insertions / alignment artifacts relative to the ancestor and are
    excluded.  Bottom-up over the tree: leaf sets = observed base; internal
    sets = intersection of child sets if nonempty else union.  The root state
    per column is the highest-frequency state in the root set (ties broken by
    state order A < C < G < T); the parsimony intersection rule constrains
    the root when the two root subtrees share exactly one possible state.

    Returns (ancestral_seq, n_columns, n_emitted, n_gated_out)."""
    n_nodes = root + 1
    post = []

    def dfs(u):
        chs = children.get(u)
        if chs:
            dfs(chs[0])
            dfs(chs[1])
        post.append(u)

    dfs(root)
    internal = [u for u in post if u in children]
    ch1 = [children[u][0] for u in internal]
    ch2 = [children[u][1] for u in internal]

    n, width = encoded.shape
    min_cov = max(1, math.ceil(n_leaves * coverage_threshold))
    ALL_BASES = 0b1111
    emitted = []
    n_gated = 0
    n_columns = 0
    for c0 in range(0, width, chunk):
        c1 = min(c0 + chunk, width)
        cols = encoded[:, c0:c1].astype(np.int16)   # (n, cw)
        cw = c1 - c0
        nongap = (cols < 4)                          # (n, cw) bool
        cov = nongap.sum(axis=0)                     # (cw,) base count
        n_columns += cw
        n_gated += int((cov < min_cov).sum())
        active = cov >= min_cov                      # columns to reconstruct
        masks = np.zeros((n_nodes, cw), dtype=np.int64)
        freqs = np.zeros((n_nodes, 5, cw), dtype=np.int16)
        # leaves: base state -> single-state mask; gap -> all bases, no counts
        base_masks = np.left_shift(np.int64(1), cols)
        masks[:n_leaves] = np.where(nongap, base_masks, ALL_BASES)
        onehot = (cols[:, :, None] ==
                  np.arange(4, dtype=np.int16)[None, None, :])
        freqs[:n_leaves, :4, :] = np.where(
            nongap[:, None, :], onehot.transpose(0, 2, 1).astype(np.int16), 0)
        for k, u in enumerate(internal):
            a, b = ch1[k], ch2[k]
            inter = masks[a] & masks[b]
            masks[u] = np.where(inter != 0, inter, masks[a] | masks[b])
            freqs[u] = freqs[a] + freqs[b]
        root_masks = masks[root]                                # (cw,) bitmask
        # root state = highest-frequency state WITHIN the root's Fitch set
        # (parsimony constraint); ties broken by state order A < C < G < T
        allowed = (root_masks[:, None] >> np.arange(4)) & 1     # (cw, 4)
        scores = np.where(allowed, freqs[root, :4, :].T, -1).astype(np.int64)
        root_state = np.argmax(scores, axis=1)                  # (cw,) 0-3
        for j in range(cw):
            if not active[j]:
                continue
            emitted.append(STATES[int(root_state[j])])
    return ''.join(emitted), n_columns, len(emitted), n_gated


def partition_ancestral(seq_data, coverage_threshold=DEFAULT_COVERAGE_THRESHOLD,
                       chunk=512):
    """Ancestral representative for a partition (NJ + Fitch on the MAF block);
    single-member blocks return the member's non-gap sequence; FASTA-only
    partitions fall back to the longest member sequence.

    Returns (seq_or_None, diag, newick_or_None)."""
    maf = seq_data.get('maf')
    if maf:
        width = min(len(r['seq']) for r in maf)
        if width == 0:
            return None, {'source': 'ancestral_empty'}, None
        encoded = encode_alignment(maf, width)
        n = len(maf)
        names = [rec['name'] for rec in maf]
        if n == 1:
            seq = ''.join(c for c in maf[0]['seq'][:width].upper()
                          if c in 'ACGT')
            return (seq,
                    {'source': 'ancestral_single', 'n_seqs': 1,
                     'aln_width': width},
                    f'{_sanitize_name(names[0])};')
        D, shared = compute_pairwise_distances(encoded)
        newick, children, root = neighbor_joining(D, names)
        seq, n_cols, n_emit, n_gated = fitch_ancestral(
            encoded, children, root, n, coverage_threshold, chunk=chunk)
        mean_depth = float((encoded < 4).sum(axis=0).mean()) if width else 0.0
        diag = {'source': 'ancestral_nj_fitch', 'n_seqs': n,
                'aln_width': width, 'columns': n_cols, 'emitted': n_emit,
                'gated_out_columns': n_gated,
                'mean_depth': round(mean_depth, 2)}
        return seq, diag, newick
    fasta = seq_data.get('fasta')
    if fasta:
        longest = max(fasta, key=lambda r: len(r[1]))
        return (longest[1],
                {'source': 'fasta_longest', 'n_seqs': len(fasta),
                 'longest_name': longest[0], 'longest_len': len(longest[1])},
                None)
    return None, {'source': 'none'}, None


# ─── Genome building ────────────────────────────────────────────────────────

def build_genome(path, reps, gap_size=DEFAULT_GAP_SIZE):
    """Concatenate the representatives of a sampled path in path order.

    Partitions without sequence are skipped (recorded).  Returns
    (genome_seq, genome_pids, skipped_no_seq)."""
    seq = ''
    genome_pids = []
    skipped = []
    first_emit = True
    for pid in path:
        rep = reps.get(pid)
        if not rep:
            skipped.append(pid)
            continue
        if not first_emit and gap_size > 0:
            seq += 'N' * gap_size
        seq += rep
        genome_pids.append(pid)
        first_emit = False
    return seq, genome_pids, skipped


def wrap_fasta(seq, width=80):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


# ─── Metrics ────────────────────────────────────────────────────────────────

def compute_sampling_metrics(samples, weights):
    """Per-partition sampling frequency + summary metrics over the draws."""
    n_samples = len(samples)
    n_sampled = Counter()
    lengths = []
    supports = []
    for s in samples:
        n_sampled.update(set(s['path']))
        lengths.append(s['genome_length_bp'])
        supports.append(s['support'])
    freq = {p: n_sampled.get(p, 0) / n_samples for p in weights}
    order = sorted(weights, key=lambda p: -weights[p])
    k = max(1, len(order) // 10)
    high = order[:k]
    low = order[-k:]
    high_mean = round(sum(freq[p] for p in high) / k, 4) if k else 0.0
    low_mean = round(sum(freq[p] for p in low) / k, 4) if k else 0.0
    metrics = {
        'n_samples': n_samples,
        'mean_sample_genome_length_bp': round(sum(lengths) / len(lengths), 1),
        'min_sample_genome_length_bp': min(lengths),
        'max_sample_genome_length_bp': max(lengths),
        'mean_sampling_frequency_top10pct_weight': high_mean,
        'mean_sampling_frequency_bottom10pct_weight': low_mean,
        'fraction_partitions_sampled_in_ge_half_of_runs':
            round(sum(1 for p in weights if freq[p] >= 0.5) / len(weights), 4),
    }
    return metrics, freq


# ─── Outputs ────────────────────────────────────────────────────────────────

def _community_name(prefix):
    return os.path.basename(prefix)


def write_outputs(prefix, mode, samples, best_index, weights, freq,
                  consensi, consensus_diag, ancestral, ancestral_diag,
                  trees, coverage_stats, adj, genome_seq, genome_pids,
                  skipped_no_seq, metrics, params, inputs, part_lengths):
    """Write all per-community outputs.  Returns dict of written files."""
    parent = os.path.dirname(prefix)
    if parent:
        os.makedirs(parent, exist_ok=True)
    files = {}
    partitions = sorted(weights)
    position_best = {p: i for i, p in enumerate(samples[best_index]['path'])}

    # ── (a) traversal JSON ──────────────────────────────────────────────────
    best = samples[best_index]
    best_observed = sum(1 for i in range(len(best['path']) - 1)
                        if adj.get(best['path'][i], {}).get(best['path'][i + 1], 0) > 0)
    partitions_records = {}
    for pid in partitions:
        cs = coverage_stats.get(pid, {})
        con = consensi.get(pid)
        anc = ancestral.get(pid) if mode == 'ancestral' else None
        partitions_records[str(pid)] = {
            'occurrence': cs.get('occurrence', 0),
            'fraction': cs.get('fraction', 0.0),
            'weight': round(weights[pid], 6),
            'sampling_frequency': round(freq.get(pid, 0.0), 4),
            'has_seq': con is not None or anc is not None,
            'consensus_len': len(con) if con else 0,
            'ancestral_len': len(anc) if anc else 0,
            'representative_len': part_lengths.get(pid, 0),
            'seq_source': (consensus_diag.get(pid, {}).get('source', 'none')
                           if mode == 'ml' else
                           ancestral_diag.get(pid, {}).get('source', 'none')),
            'in_best_genome': pid in position_best,
            'best_position': position_best.get(pid),
        }
    traversal_json = {
        'community': _community_name(prefix),
        'mode': mode,
        'level1': (
            'Weighted path sampling over the partition set: weight w_p = '
            '(occurrence_p / N)^alpha; first partition proportional to w; '
            'next partition b proportional to w_b * (1 + beta*adj[current][b]) '
            'among unvisited partitions fitting the remaining budget. '
            'No hard rarity threshold: every partition with occurrence >= 1 '
            'has w > 0. The stitched genome uses the sampled traversal with '
            'the highest support sum(log2 w).'),
        'level2': (
            "'ml' = majority-rule consensus of the aligned block (per-column "
            "multinomial MLE, coverage-gated); 'ancestral' = neighbor-joining "
            "tree from pairwise p-distances (gap = missing) + Fitch parsimony "
            "with frequency-based deterministic root resolution; columns "
            "below the coverage gate excluded)."),
        'parameters': params,
        'inputs': inputs,
        'n_prophages': inputs['n_prophages'],
        'n_partitions_total': len(partitions),
        'n_partitions_with_seq': sum(1 for p in partitions
                                     if part_lengths.get(p, 0)),
        'weights': {str(p): round(weights[p], 6) for p in partitions},
        'samples': [{'index': s['index'], 'path': s['path'],
                     'support': s['support'], 'adj_pairs': s['adj_pairs'],
                     'n_partitions': len(s['path']),
                     'genome_length_bp': s['genome_length_bp']}
                    for s in samples],
        'best_sample_index': best_index,
        'best_sample': {'path': best['path'], 'support': best['support'],
                        'adj_pairs': best['adj_pairs'],
                        'n_partitions': len(best['path']),
                        'genome_length_bp': best['genome_length_bp'],
                        'observed_adjacent_pairs_in_genome': best_observed,
                        'skipped_no_seq': skipped_no_seq},
        'genome': {'mode': mode, 'path': best['path'], 'pids': genome_pids,
                   'length_bp': len(genome_seq)},
        'metrics': metrics,
        'partitions': partitions_records,
    }
    json_path = prefix + '.traversal.json'
    with open(json_path, 'w') as f:
        json.dump(traversal_json, f, indent=2)
    files['traversal_json'] = json_path

    # ── (b) per-partition ML representatives ────────────────────────────────
    cons_path = prefix + '.consensus.fa'
    with open(cons_path, 'w') as f:
        for pid in partitions:
            con = consensi.get(pid)
            if not con:
                continue
            cs = coverage_stats.get(pid, {})
            src = consensus_diag.get(pid, {}).get('source', '?')
            f.write(f'>partition{pid} occurrence={cs.get("occurrence", 0)} '
                    f'len={len(con)} source={src}\n')
            f.write(wrap_fasta(con) + '\n')
    files['consensus_fasta'] = cons_path

    # ── (c) stitched ML genome ──────────────────────────────────────────────
    ml_path = prefix + '.ml.fa'
    with open(ml_path, 'w') as f:
        f.write(f'>community_{_community_name(prefix)}_ML '
                f'mode=ml n_partitions={len(genome_pids)} '
                f'len={len(genome_seq)} seed={params["seed"]}\n')
        f.write(wrap_fasta(genome_seq) + '\n')
    files['ml_genome_fasta'] = ml_path

    if mode == 'ancestral':
        # ── per-partition ancestral representatives ─────────────────────────
        anc_path = prefix + '.ancestral.fa'
        with open(anc_path, 'w') as f:
            for pid in partitions:
                anc = ancestral.get(pid)
                if not anc:
                    continue
                cs = coverage_stats.get(pid, {})
                src = ancestral_diag.get(pid, {}).get('source', '?')
                f.write(f'>partition{pid} occurrence={cs.get("occurrence", 0)} '
                        f'len={len(anc)} source={src}\n')
                f.write(wrap_fasta(anc) + '\n')
        files['ancestral_fasta'] = anc_path

        # ── ancestral stitched genome ───────────────────────────────────────
        anc_genome = ''
        first = True
        for pid in genome_pids:
            rep = ancestral.get(pid)
            if not rep:
                continue
            if not first and params.get('gap_size', 0) > 0:
                anc_genome += 'N' * params['gap_size']
            anc_genome += rep
            first = False
        ag_path = prefix + '.ancestral.genome.fa'
        with open(ag_path, 'w') as f:
            f.write(f'>community_{_community_name(prefix)}_ancestral '
                    f'mode=ancestral n_partitions={len(genome_pids)} '
                    f'len={len(anc_genome)} seed={params["seed"]}\n')
            f.write(wrap_fasta(anc_genome) + '\n')
        files['ancestral_genome_fasta'] = ag_path

        # ── per-partition NJ trees (QC) ─────────────────────────────────────
        trees_path = prefix + '.trees.nwk'
        with open(trees_path, 'w') as f:
            for pid in partitions:
                nwk = trees.get(pid)
                if nwk is None:
                    continue
                f.write(f'partition{pid}\t{nwk}\n')
        files['trees_nwk'] = trees_path

    # ── (d) coverage + sampling statistics ──────────────────────────────────
    tsv_path = prefix + '.coverage.tsv'
    with open(tsv_path, 'w') as f:
        f.write('\t'.join(['partition', 'occurrence', 'fraction', 'weight',
                           'n_sampled', 'sampling_frequency', 'has_seq',
                           'consensus_len', 'ancestral_len', 'in_best_genome',
                           'best_position']) + '\n')
        for pid in partitions:
            cs = coverage_stats.get(pid, {})
            con = consensi.get(pid)
            anc = ancestral.get(pid) if mode == 'ancestral' else None
            n_sampled = int(round(freq.get(pid, 0.0) * len(samples)))
            f.write('\t'.join([
                str(pid),
                str(cs.get('occurrence', 0)),
                f"{cs.get('fraction', 0.0):.6f}",
                f"{weights[pid]:.6f}",
                str(n_sampled),
                f"{freq.get(pid, 0.0):.4f}",
                'yes' if (con or anc) else 'no',
                str(len(con) if con else 0),
                str(len(anc) if anc else 0),
                'yes' if pid in position_best else 'no',
                str(position_best.get(pid, '')),
            ]) + '\n')
    files['coverage_tsv'] = tsv_path

    # ── stats JSON ──────────────────────────────────────────────────────────
    stats_json = {'community': _community_name(prefix), 'mode': mode,
                  'parameters': params, 'metrics': metrics,
                  'best_sample_index': best_index,
                  'genome_length_bp': len(genome_seq),
                  'n_partitions_total': len(partitions),
                  'n_partitions_with_seq': sum(1 for p in partitions
                                               if part_lengths.get(p, 0))}
    stats_path = prefix + '.stats.json'
    with open(stats_path, 'w') as f:
        json.dump(stats_json, f, indent=2)
    files['stats_json'] = stats_path

    return files


# ─── Main pipeline ──────────────────────────────────────────────────────────

def run_traversal(partitions_dir, bed_path, output_prefix, mode=DEFAULT_MODE,
                  n_samples=DEFAULT_N_SAMPLES, alpha=DEFAULT_ALPHA,
                  seed=DEFAULT_SEED, beta=DEFAULT_BETA,
                  max_length=DEFAULT_MAX_LENGTH,
                  coverage_threshold=DEFAULT_COVERAGE_THRESHOLD,
                  gap_size=DEFAULT_GAP_SIZE):
    """Run the full two-level ML pipeline.  Returns (result_dict, files_dict)."""
    t0 = time.time()

    # 1. Parse BED -> normalized per-prophage partition lists
    raw = parse_bed(bed_path)
    normalized = normalize_partition_lists(raw)
    n_prophages = len(normalized)
    if n_prophages == 0:
        raise RuntimeError(f'no prophage partition assignments found in {bed_path}')
    partitions = sorted({p for pids in normalized.values() for p in pids})
    print(f'[1] BED: {n_prophages} prophages, {len(partitions)} partitions')

    # 2. Coverage statistics + sampling weights (Level 1)
    coverage_stats = compute_coverage_stats(normalized, n_prophages)
    weights = compute_sampling_weights(coverage_stats, alpha)
    w_min = min(weights.values())
    w_max = max(weights.values())
    print(f'[2] weights: alpha={alpha}, range [{w_min:.6g}, {w_max:.6g}] '
          f'— all positive, nothing hard-excluded')

    # 3. Sequence data + Level-2 representatives
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
    print(f'[3] ML consensus: {n_with_seq}/{len(partitions)} partitions have a '
          f'representative sequence')

    ancestral = {}
    ancestral_diag = {}
    trees = {}
    if mode == 'ancestral':
        for pid in partitions:
            seq, diag, newick = partition_ancestral(seq_data[pid])
            if seq:
                ancestral[pid] = seq
            ancestral_diag[pid] = diag
            if newick:
                trees[pid] = newick
        print(f'[3b] ancestral: {len(ancestral)}/{len(partitions)} partitions '
              f'reconstructed (NJ + Fitch)')

    # 4. Adjacency graph (Level-1 support)
    adj, first_counts, last_counts, occurrence = build_adjacency_graph(normalized)
    n_edges = sum(len(t) for t in adj.values())
    print(f'[4] adjacency graph: {len(adj)} nodes, {n_edges} directed '
          f'observed edges')

    # 5. Representative lengths + weighted path sampling
    part_lengths = {pid: len(consensi[pid]) if pid in consensi else 0
                    for pid in partitions}
    samples, rng = sample_traversals(weights, adj, part_lengths, max_length,
                                     n_samples, seed, beta)
    best_index = max(range(n_samples), key=lambda i: samples[i]['support'])
    best = samples[best_index]
    print(f'[5] sampling: {n_samples} traversals (seed={seed}); best sample '
          f'#{best_index}: {len(best["path"])} partitions, '
          f'{best["genome_length_bp"]:,} bp, support={best["support"]}')

    # 6. Sampling statistics
    metrics, freq = compute_sampling_metrics(samples, weights)
    print(f'[6] sampling distribution: top-10% weight partitions sampled in '
          f'{metrics["mean_sampling_frequency_top10pct_weight"]*100:.1f}% of '
          f'runs; bottom-10% in '
          f'{metrics["mean_sampling_frequency_bottom10pct_weight"]*100:.2f}%')

    # 7. Stitched genome from the best sampled traversal
    reps = ancestral if mode == 'ancestral' else consensi
    genome_seq, genome_pids, skipped_no_seq = build_genome(best['path'], reps,
                                                           gap_size)
    print(f'[7] {mode} genome: {len(genome_seq):,} bp from '
          f'{len(genome_pids)} partitions (budget {max_length:,} bp); '
          f'{len(skipped_no_seq)} path partitions skipped (no sequence)')

    params = {
        'mode': mode, 'alpha': alpha, 'beta': beta, 'n_samples': n_samples,
        'seed': seed, 'max_length': max_length,
        'coverage_threshold': coverage_threshold, 'gap_size': gap_size,
    }
    inputs = {'bed': bed_path, 'partitions_dir': partitions_dir,
              'n_prophages': n_prophages, 'n_partitions': len(partitions)}

    # 8. Outputs
    files = write_outputs(output_prefix, mode, samples, best_index, weights,
                          freq, consensi, consensus_diag, ancestral,
                          ancestral_diag, trees, coverage_stats, adj,
                          genome_seq, genome_pids, skipped_no_seq, metrics,
                          params, inputs, part_lengths)
    print(f'[8] wrote {len(files)} output files under {output_prefix}.*')
    print(f'    completed in {time.time() - t0:.1f}s')

    result = {
        'mode': mode, 'samples': samples, 'best_sample_index': best_index,
        'metrics': metrics, 'files': files,
        'genome_length_bp': len(genome_seq),
        'genome_pids': genome_pids,
        'n_partitions_total': len(partitions),
    }
    return result, files


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Two-level ML inference for typical (ML) and ancestral '
                    'phage genomes: Level-1 weighted path sampling over the '
                    'partition adjacency graph (rare partitions less likely, '
                    'never hard-excluded) x Level-2 per-partition ML '
                    'representative (majority consensus or ancestral '
                    'reconstruction).')
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
                             '<prefix>.ml.fa, <prefix>.coverage.tsv, '
                             '<prefix>.stats.json (plus <prefix>.ancestral.fa, '
                             '<prefix>.ancestral.genome.fa, '
                             '<prefix>.trees.nwk in ancestral mode)')
    parser.add_argument('--mode', choices=['ml', 'ancestral'],
                        default=DEFAULT_MODE,
                        help='Level-2 representative per partition: '
                             '"ml" = majority-rule consensus (default); '
                             '"ancestral" = NJ tree + Fitch parsimony '
                             'ancestral reconstruction')
    parser.add_argument('--n-samples', type=int, default=DEFAULT_N_SAMPLES,
                        help='Number of weighted traversal draws '
                             '(default 10); per-partition sampling frequency '
                             'is reported over these draws')
    parser.add_argument('--alpha', type=float, default=DEFAULT_ALPHA,
                        help='Exponent of the sampling weight '
                             'w = occurrence_fraction^alpha (default 1.0); '
                             'alpha > 1 makes rare partitions even less '
                             'likely, alpha < 1 flattens the distribution')
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED,
                        help='RNG seed for the weighted sampling (default 42); '
                             'same seed -> identical samples and outputs')
    parser.add_argument('--adj-weight', type=float, default=DEFAULT_BETA,
                        help='Adjacency modulation beta in the step '
                             'distribution P(b) ~ w_b * (1 + beta*adj[a][b]) '
                             '(default 1.0)')
    parser.add_argument('--max-length', type=int, default=DEFAULT_MAX_LENGTH,
                        help='Length budget (bp) for the stitched genome '
                             '(default 150000; keeps the output at '
                             'phage-typical scale)')
    parser.add_argument('--coverage-threshold', type=float,
                        default=DEFAULT_COVERAGE_THRESHOLD,
                        help='Column coverage gate for MAF consensus and '
                             'ancestral reconstruction (fraction of block '
                             'sequences with a non-gap base; default 0.25 = '
                             'true conserved locus for staggered bundles, '
                             'matching the validated stitch_algorithm.py; '
                             '0.0 = union mosaic, byte-identical to the '
                             'committed 2363ece runs)')
    parser.add_argument('--gap-size', type=int, default=DEFAULT_GAP_SIZE,
                        help='N-runs inserted between partition blocks in the '
                             'stitched genome (default 0 = plain concatenation)')
    args = parser.parse_args(argv)

    if not os.path.isdir(args.partitions_dir):
        parser.error(f'--partitions-dir not a directory: {args.partitions_dir}')
    if not os.path.exists(args.bed):
        parser.error(f'--bed not found: {args.bed}')
    if args.mode not in ('ml', 'ancestral'):
        parser.error(f'--mode must be ml or ancestral, got {args.mode}')
    if args.n_samples < 1:
        parser.error('--n-samples must be >= 1')
    if args.alpha < 0:
        parser.error('--alpha must be >= 0')
    if args.adj_weight < 0:
        parser.error('--adj-weight must be >= 0')
    if args.coverage_threshold < 0 or args.coverage_threshold > 1:
        parser.error('--coverage-threshold must be in [0, 1]')
    if args.max_length <= 0:
        parser.error('--max-length must be positive')

    result, files = run_traversal(
        args.partitions_dir, args.bed, args.output, mode=args.mode,
        n_samples=args.n_samples, alpha=args.alpha, seed=args.seed,
        beta=args.adj_weight, max_length=args.max_length,
        coverage_threshold=args.coverage_threshold, gap_size=args.gap_size)
    return result


if __name__ == '__main__':
    main()
