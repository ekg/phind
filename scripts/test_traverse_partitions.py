#!/usr/bin/env python3
"""Tests for scripts/traverse_partitions.py — two-level ML inference
(Level 1 weighted path sampling across partitions; Level 2 ML consensus or
ancestral NJ+Fitch reconstruction per partition).

Run:  python3 scripts/test_traverse_partitions.py
      (or: pytest scripts/test_traverse_partitions.py)

Covers:
  * unit tests: BED parsing/normalization, adjacency graph, majority-rule
    consensus + coverage gate, FASTA fallback, sampling weights
    (fraction^alpha, all positive), weighted path sampling (determinism,
    common-vs-rare distribution, no hard exclusion, alpha effect),
    neighbor-joining (known matrices, determinism) and Fitch parsimony
    (conserved/majority/tie columns, tree-informative columns, gap =
    missing, column gate, 1-/2-sequence blocks, determinism).
  * integration on the archived community_3 partition set: default ML run
    emits traversal(s) + consensus + ML genome + coverage/sampling report;
    sampling-frequency table shows common partitions near 100% and rare
    partitions at a rate consistent with their weight; genome length is
    phage-typical (tens of kb .. ~150 kb); same seed -> byte-identical
    outputs.  Ancestral mode emits per-partition ancestral FASTA, an
    ancestral genome and NJ trees, also phage-typical length.
  * integration on a synthetic community (per-partition BED dir layout):
    both modes emit all outputs; FASTA-only partitions fall back to the
    longest member; determinism holds.
"""

import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import traverse_partitions as tp  # noqa: E402

COMMUNITY3_BED = os.path.join(REPO, 'research/stitching/community_3_partitions.bed')
COMMUNITY3_DIR = os.path.join(REPO, 'research/stitching/inputs/community_3_partitions')

HAS_COMMUNITY3 = all(os.path.exists(p) for p in [COMMUNITY3_BED, COMMUNITY3_DIR])


# ─── helpers / fixtures ────────────────────────────────────────────────────

def _rand_dna(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def make_synthetic_community(root):
    """Create a synthetic community mimicking `impg partition` output.

    Layout:
      root/partitions_bed/partition<N>.bed   (prophage start end, 3 columns)
      root/partitions/partition<N>.maf       (alignment block) or
      root/partitions/partition<N>.fa        (member FASTA bundle, no MAF)

    Genome model: 12 prophages; core partitions 1-8 shared by most genomes,
    medium 9-12 by a subset, rare 13-30 by a single genome each.  MAF blocks
    are staggered segment bundles (union mosaic when gated at 0.0); partitions
    29/30 have only FASTA (no MAF) to exercise the longest-seq fallback.
    """
    bed_dir = os.path.join(root, 'partitions_bed')
    part_dir = os.path.join(root, 'partitions')
    os.makedirs(bed_dir)
    os.makedirs(part_dir)

    rng = random.Random(7)
    order = list(range(1, 31))
    genomes = {}
    for g in range(1, 13):
        name = f'G{g:02d}'
        if g <= 9:
            pids = list(range(1, 9))        # 9/12 carry the core
        else:
            pids = [2, 3, 4]
        if g % 2 == 0:
            pids += list(range(9, 13))      # half carry the medium modules
        if g <= 2:
            pids += [13 + g - 1, 13 + g]    # two genomes carry rare modules
        if g == 11:
            pids += [15, 16, 17, 29]        # 29 = FASTA-only partition
        if g == 12:
            pids += [18, 19, 20, 30]        # 30 = FASTA-only partition
        genomes[name] = pids

    part_members = {}
    for name, pids in genomes.items():
        start = 0
        for pid in pids:
            part_members.setdefault(pid, []).append((name, start))
            with open(os.path.join(bed_dir, f'partition{pid}.bed'), 'a') as f:
                f.write(f'{name}\t{start}\t{start + 1000}\n')
            start += 1000

    for pid in sorted(part_members):
        members = part_members[pid]
        if pid in (29, 30):
            with open(os.path.join(part_dir, f'partition{pid}.fa'), 'w') as f:
                for name, _ in members:
                    f.write(f'>{name}\n{_rand_dna(rng, 1000)}\n')
            continue
        block_rows = []
        for name, _ in members:
            seq = _rand_dna(rng, 1000)
            block_rows.append((name, 400, seq))
        width = max(s + len(seq) for _, s, seq in block_rows)
        with open(os.path.join(part_dir, f'partition{pid}.maf'), 'w') as f:
            f.write('a score=1\n')
            for name, s, seq in block_rows:
                f.write(f's {name} {s} {len(seq)} + 10000 '
                        f'{seq}{"N" * (width - s - len(seq))}\n')
    return genomes


def _run_cli(args, cwd=REPO):
    proc = subprocess.run([sys.executable, 'scripts/traverse_partitions.py'] + args,
                          capture_output=True, text=True, cwd=cwd)
    return proc.returncode, proc.stdout + proc.stderr


def _read_fasta(path):
    seqs = {}
    name, cur = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(cur)
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line)
    if name:
        seqs[name] = ''.join(cur)
    return seqs


def _records(rows):
    return [{'name': n, 'start': 0, 'size': len(s), 'strand': '+',
             'src_size': len(s), 'seq': s} for n, s in rows]


def _freqs(samples, pids):
    n = len(samples)
    c = Counter()
    for s in samples:
        c.update(set(s['path']))
    return {p: c.get(p, 0) / n for p in pids}


# ─── BED parsing ───────────────────────────────────────────────────────────

class ParseNormalizeTests(unittest.TestCase):
    def test_combined_bed_and_normalize(self):
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, 'p.bed')
            with open(bed, 'w') as f:
                f.write('G1\t0\t1000\t1\n')
                f.write('G1\t1000\t2000\t1\n')   # contiguous run -> merge
                f.write('G1\t2000\t3000\t2\n')
                f.write('G1\t3000\t4000\t2\n')
                f.write('G1\t4000\t5000\t2\n')   # duplicate of prev -> collapse
                f.write('G2\t0\t1000\t3\n')
            raw = tp.parse_bed(bed)
            norm = tp.normalize_partition_lists(raw)
            self.assertEqual(norm['G1'], [1, 2])
            self.assertEqual(norm['G2'], [3])

    def test_bed_dir_layout(self):
        with tempfile.TemporaryDirectory() as d:
            for pid in (3, 7):
                with open(os.path.join(d, f'partition{pid}.bed'), 'w') as f:
                    f.write(f'G1\t0\t1000\nG2\t0\t1000\n')
            raw = tp.parse_bed(d)
            self.assertEqual(sorted(raw), ['G1', 'G2'])
            pids = sorted({p for v in raw.values() for p, _, _ in v})
            self.assertEqual(pids, [3, 7])


# ─── adjacency graph ───────────────────────────────────────────────────────

class AdjacencyTests(unittest.TestCase):
    def test_adjacency_graph_counts(self):
        normalized = {
            'G1': [1, 2, 3],
            'G2': [1, 2, 3],
            'G3': [1, 3],
            'G4': [2, 4],
        }
        adj, first, last, occ = tp.build_adjacency_graph(normalized)
        self.assertEqual(adj[1][2], 2)     # G1,G2
        self.assertEqual(adj[2][3], 2)
        self.assertEqual(adj[1][3], 1)     # G3
        self.assertEqual(adj[2][4], 1)     # G4
        self.assertEqual(occ[1], 3)        # G1,G2,G3
        self.assertEqual(first[1], 3)      # G1,G2,G3
        self.assertEqual(last[3], 3)       # G1,G2,G3


# ─── consensus (ML within) ─────────────────────────────────────────────────

class ConsensusTests(unittest.TestCase):
    def test_majority_rule(self):
        recs = _records([('a', 'ACGT'), ('b', 'ACGA'), ('c', 'ACGC')])
        cons, diag = tp.compute_partition_consensus(recs, 0.0)
        self.assertEqual(cons, 'ACGT')     # majority per column
        self.assertEqual(diag['n_seqs'], 3)

    def test_coverage_threshold_gate(self):
        recs = _records([('a', 'AA--'), ('b', 'AG--'), ('c', 'ATCC')])
        cons, _ = tp.compute_partition_consensus(recs, 0.5)
        # ceil(3*0.5)=2: col0 'A' (3 bases); col1 tie A/G/T broken by
        # first-encountered -> 'A'; cols 2/3 have 1 base each -> gated out
        self.assertEqual(cons, 'AA')
        cons0, _ = tp.compute_partition_consensus(recs, 0.0)
        self.assertEqual(cons0, 'AACC')    # 0.0 = any covered column

    def test_longest_seq_fallback_without_maf(self):
        seq_data = {'fasta': [('a', 'ACGT' * 2), ('b', 'ACGT' * 5)]}
        rep, diag = tp.partition_representative(seq_data, 0.0)
        self.assertEqual(rep, 'ACGT' * 5)
        self.assertEqual(diag['source'], 'fasta_longest')

    def test_maf_preferred_over_fasta(self):
        maf = _records([('a', 'ACGT')])
        seq_data = {'maf': maf, 'fasta': [('a', 'TTTTTT')]}
        rep, diag = tp.partition_representative(seq_data, 0.0)
        self.assertEqual(rep, 'ACGT')
        self.assertEqual(diag['source'], 'maf_consensus')


# ─── coverage stats + weights ──────────────────────────────────────────────

class CoverageAndWeightTests(unittest.TestCase):
    def test_fractions(self):
        normalized = {f'G{i}': [1, 2] for i in range(10)}
        normalized['G1'] = [1, 2, 3]
        stats = tp.compute_coverage_stats(normalized, 10)
        self.assertEqual(stats[1]['occurrence'], 10)
        self.assertEqual(stats[1]['fraction'], 1.0)
        self.assertEqual(stats[3]['occurrence'], 1)
        self.assertEqual(stats[3]['fraction'], 0.1)

    def test_weights_positive_and_alpha(self):
        stats = {p: {'fraction': f} for p, f in [(1, 0.8), (2, 0.5), (3, 0.001)]}
        w1 = tp.compute_sampling_weights(stats, alpha=1.0)
        self.assertAlmostEqual(w1[1], 0.8)
        self.assertAlmostEqual(w1[3], 0.001)
        # all strictly positive: nothing hard-excluded
        self.assertTrue(all(w > 0 for w in w1.values()))
        # alpha sharpens: rare partition gets even less likely
        w3 = tp.compute_sampling_weights(stats, alpha=3.0)
        self.assertLess(w3[3] / w3[1], w1[3] / w1[1])
        # alpha=0 -> uniform
        w0 = tp.compute_sampling_weights(stats, alpha=0.0)
        self.assertEqual(set(w0.values()), {1.0})


# ─── Level 1: weighted path sampling ───────────────────────────────────────

class SamplingTests(unittest.TestCase):
    WEIGHTS = {1: 0.9, 2: 0.9, 3: 0.9, 4: 0.9, 5: 0.9,
               6: 0.5, 7: 0.5, 8: 0.2, 9: 0.1, 10: 0.1,
               11: 0.05, 12: 0.05}

    def _setup(self, budget=8000, lengths=None):
        # core partitions 1-5 are mutually adjacent (observed adjacency mass),
        # medium partitions 6/7 attach to the core, so walks entering the
        # graph reliably visit the core before the budget is exhausted
        adj = {1: {2: 10}, 2: {3: 10}, 3: {4: 10}, 4: {5: 10}, 5: {1: 10},
               6: {1: 10}, 7: {2: 10}}
        part_lengths = {p: (lengths or 1000) for p in self.WEIGHTS}
        return adj, part_lengths

    def test_same_seed_identical_samples(self):
        adj, plen = self._setup()
        s1, _ = tp.sample_traversals(self.WEIGHTS, adj, plen, 8000, 50, 42)
        s2, _ = tp.sample_traversals(self.WEIGHTS, adj, plen, 8000, 50, 42)
        self.assertEqual([s['path'] for s in s1],
                         [s['path'] for s in s2])

    def test_different_seeds_differ(self):
        adj, plen = self._setup()
        s1, _ = tp.sample_traversals(self.WEIGHTS, adj, plen, 8000, 20, 1)
        s2, _ = tp.sample_traversals(self.WEIGHTS, adj, plen, 8000, 20, 2)
        self.assertNotEqual([s['path'] for s in s1],
                            [s['path'] for s in s2])

    def test_common_sampled_more_than_rare(self):
        adj, plen = self._setup()
        samples, _ = tp.sample_traversals(self.WEIGHTS, adj, plen, 8000, 500, 42)
        f = _freqs(samples, self.WEIGHTS)
        common = sum(f[p] for p in range(1, 6)) / 5
        rare = sum(f[p] for p in (9, 10, 11, 12)) / 4
        self.assertGreater(common, 0.9, 'common partitions ~ near-100%')
        self.assertLess(rare, 0.4, 'rare partitions less likely')
        self.assertGreater(common, 2 * rare)

    def test_no_hard_exclusion(self):
        # unlimited budget -> every partition (however rare) is visited in
        # every run: nothing is excluded by a threshold
        adj, plen = self._setup(budget=10_000_000)
        samples, _ = tp.sample_traversals(self.WEIGHTS, adj, plen,
                                          10_000_000, 30, 42)
        f = _freqs(samples, self.WEIGHTS)
        for p in self.WEIGHTS:
            self.assertEqual(f[p], 1.0, f'partition {p} never excluded')

    def test_alpha_sharpens_distribution(self):
        adj, plen = self._setup()
        w_flat = {p: w ** 0.5 for p, w in self.WEIGHTS.items()}
        w_peak = {p: w ** 3.0 for p, w in self.WEIGHTS.items()}
        s_flat, _ = tp.sample_traversals(w_flat, adj, plen, 5000, 500, 7)
        s_peak, _ = tp.sample_traversals(w_peak, adj, plen, 5000, 500, 7)
        f_flat = _freqs(s_flat, self.WEIGHTS)
        f_peak = _freqs(s_peak, self.WEIGHTS)
        rare_flat = sum(f_flat[p] for p in (9, 10, 11, 12)) / 4
        rare_peak = sum(f_peak[p] for p in (9, 10, 11, 12)) / 4
        self.assertLess(rare_peak, rare_flat,
                        'higher alpha -> rare partitions sampled less often')

    def test_sampler_respects_budget(self):
        adj, plen = self._setup(budget=3000)  # ~3 x 1000 bp
        samples, _ = tp.sample_traversals(self.WEIGHTS, adj, plen, 3000, 20, 42)
        for s in samples:
            self.assertLessEqual(s['genome_length_bp'], 3000)


# ─── Level 2 (ancestral): NJ + Fitch ───────────────────────────────────────

class NeighborJoiningTests(unittest.TestCase):
    def test_nj_recovers_identical_siblings(self):
        # a0==a1 and a2==a3: NJ must join the identical pair (0,1) first
        rng = random.Random(3)
        s0 = _rand_dna(rng, 400)
        s2 = _rand_dna(rng, 400)
        rows = _records([('a0', s0), ('a1', s0), ('a2', s2), ('a3', s2)])
        encoded = tp.encode_alignment(rows, len(s0))
        D, _ = tp.compute_pairwise_distances(encoded)
        newick, children, root = tp.neighbor_joining(D, [r['name'] for r in rows])
        self.assertIn('a0:0.000000,a1:0.000000', newick)

    def test_nj_deterministic(self):
        rng = random.Random(11)
        rows = _records([(f's{i}', _rand_dna(rng, 300)) for i in range(6)])
        encoded = tp.encode_alignment(rows, 300)
        D, _ = tp.compute_pairwise_distances(encoded)
        n1, c1, r1 = tp.neighbor_joining(D, [r['name'] for r in rows])
        n2, c2, r2 = tp.neighbor_joining(D, [r['name'] for r in rows])
        self.assertEqual(n1, n2)

    def test_nj_two_taxa(self):
        rows = _records([('x', 'ACGT'), ('y', 'ACGA')])
        encoded = tp.encode_alignment(rows, 4)
        D, _ = tp.compute_pairwise_distances(encoded)
        newick, children, root = tp.neighbor_joining(D, ['x', 'y'])
        self.assertIn('x:', newick)
        self.assertIn('y:', newick)
        self.assertEqual(sorted(children[root]), [0, 1])  # leaf node ids


class FitchTests(unittest.TestCase):
    def test_conserved_columns(self):
        rows = _records([('a', 'ACGT'), ('b', 'ACGT'), ('c', 'ACGT')])
        encoded = tp.encode_alignment(rows, 4)
        D, _ = tp.compute_pairwise_distances(encoded)
        newick, children, root = tp.neighbor_joining(D, ['a', 'b', 'c'])
        seq, ncol, nemit, ngated = tp.fitch_ancestral(encoded, children, root,
                                                      3, 0.25)
        self.assertEqual(seq, 'ACGT')

    def test_majority_resolution(self):
        rows = _records([('a', 'AAA'), ('b', 'AAA'), ('c', 'AAA'),
                         ('d', 'CCC'), ('e', 'CCC')])
        encoded = tp.encode_alignment(rows, 3)
        D, _ = tp.compute_pairwise_distances(encoded)
        newick, children, root = tp.neighbor_joining(D, ['a', 'b', 'c', 'd', 'e'])
        seq, _, _, _ = tp.fitch_ancestral(encoded, children, root, 5, 0.25)
        self.assertEqual(seq, 'AAA')       # 3 vs 2 -> majority base

    def test_tree_informative_column(self):
        # 6 leaves A,A,C,C,G,G on a tree where the two root subtrees have
        # sets {A,C} and {C,G} -> parsimony forces root C, while plain column
        # majority would resolve the A/G tie to A.
        children = {6: [0, 1], 7: [4, 5], 8: [6, 2], 9: [7, 3], 10: [8, 9]}
        encoded = np.array([[0], [0], [1], [1], [2], [2]], dtype=np.int8)
        seq, ncol, nemit, ngated = tp.fitch_ancestral(encoded, children, 10,
                                                      6, 0.25)
        self.assertEqual(seq, 'C')

    def test_gap_is_missing_data(self):
        # column 1 is carried by only member c (G); gap leaves a/b are
        # unconstrained (missing data) -> root state = G
        rows = _records([('a', 'A-'), ('b', 'A-'), ('c', 'AG')])
        encoded = tp.encode_alignment(rows, 2)
        D, _ = tp.compute_pairwise_distances(encoded)
        newick, children, root = tp.neighbor_joining(D, ['a', 'b', 'c'])
        seq, _, _, _ = tp.fitch_ancestral(encoded, children, root, 3, 0.25)
        self.assertEqual(seq, 'AG')        # col0 A (3x); col1 G (only carrier)

    def test_column_gate(self):
        # 5 rows: columns 0-1 covered by all, columns 2-4 covered by 1 row
        rows = _records([('a', 'AAACC'), ('b', 'AA---'), ('c', 'AA---'),
                         ('d', 'AA---'), ('e', 'AA---')])
        encoded = tp.encode_alignment(rows, 5)
        D, _ = tp.compute_pairwise_distances(encoded)
        newick, children, root = tp.neighbor_joining(D, ['a', 'b', 'c', 'd', 'e'])
        seq, ncol, nemit, ngated = tp.fitch_ancestral(encoded, children, root,
                                                      5, 0.25)
        self.assertEqual(seq, 'AA')       # cols 2-4 gated out (1/5 < 0.25)
        self.assertEqual(ngated, 3)

    def test_single_sequence_block(self):
        seq_data = {'maf': _records([('a', 'ACGT--')])}
        seq, diag, newick = tp.partition_ancestral(seq_data, 0.25)
        self.assertEqual(seq, 'ACGT')      # gaps stripped
        self.assertEqual(diag['source'], 'ancestral_single')

    def test_two_sequence_block(self):
        rows = _records([('a', 'ACGT'), ('b', 'ACGA')])
        seq, diag, newick = tp.partition_ancestral({'maf': rows}, 0.25)
        self.assertIn('a:', newick)
        self.assertIn('b:', newick)
        self.assertEqual(len(seq), 4)

    def test_fasta_fallback(self):
        seq_data = {'fasta': [('a', 'ACGT' * 3)]}
        seq, diag, newick = tp.partition_ancestral(seq_data, 0.25)
        self.assertEqual(seq, 'ACGT' * 3)
        self.assertEqual(diag['source'], 'fasta_longest')
        self.assertIsNone(newick)

    def test_ancestral_deterministic(self):
        rng = random.Random(5)
        rows = _records([(f's{i}', _rand_dna(rng, 200)) for i in range(5)])
        a1, d1, n1 = tp.partition_ancestral({'maf': rows}, 0.25)
        a2, d2, n2 = tp.partition_ancestral({'maf': rows}, 0.25)
        self.assertEqual(a1, a2)
        self.assertEqual(n1, n2)


# ─── integration: archived community_3 ─────────────────────────────────────

@unittest.skipUnless(HAS_COMMUNITY3, 'community_3 inputs not present')
class Community3MlIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='c3ml_')
        cls.prefix = os.path.join(cls.tmp, 'community_3')
        rc, out = _run_cli(['--partitions-dir', COMMUNITY3_DIR,
                            '--bed', COMMUNITY3_BED,
                            '--output', cls.prefix,
                            '--n-samples', '200', '--seed', '42'])
        cls.rc = rc
        cls.out = out
        cls.json_path = cls.prefix + '.traversal.json'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cli_runs(self):
        self.assertEqual(self.rc, 0, self.out)
        for suffix in ('.traversal.json', '.consensus.fa', '.ml.fa',
                       '.coverage.tsv', '.stats.json'):
            self.assertTrue(os.path.exists(self.prefix + suffix),
                            f'missing {suffix}')

    def test_genome_phage_typical_length(self):
        with open(self.json_path) as f:
            data = json.load(f)
        ln = data['genome']['length_bp']
        # phage-typical: tens of kb to the ~150 kb budget, NOT multi-Mbp
        self.assertGreaterEqual(ln, 30_000)
        self.assertLessEqual(ln, tp.DEFAULT_MAX_LENGTH)

    def test_sampling_frequency_common_near_100(self):
        with open(self.json_path) as f:
            data = json.load(f)
        w = data['weights']
        parts = data['partitions']
        top = sorted(w, key=lambda p: -w[p])[:10]
        mean_top = sum(parts[p]['sampling_frequency'] for p in top) / len(top)
        self.assertGreater(mean_top, 0.9,
                           'common partitions sampled in nearly every run')
        bottom = sorted(w, key=lambda p: w[p])[:10]
        mean_bottom = sum(parts[p]['sampling_frequency'] for p in bottom) / len(bottom)
        self.assertLess(mean_bottom, 0.3, 'rare partitions less likely')
        self.assertGreater(mean_top, 3 * mean_bottom)

    def test_all_weights_positive(self):
        with open(self.json_path) as f:
            data = json.load(f)
        self.assertTrue(all(w > 0 for w in data['weights'].values()))

    def test_multiple_traversals_emitted(self):
        with open(self.json_path) as f:
            data = json.load(f)
        self.assertEqual(len(data['samples']), 200)
        self.assertLessEqual(data['best_sample_index'], 199)
        self.assertEqual(len(data['best_sample']['path']),
                         len(data['genome']['path']))

    def test_coverage_report_columns(self):
        with open(self.prefix + '.coverage.tsv') as f:
            header = f.readline().strip().split('\t')
        for col in ('partition', 'occurrence', 'fraction', 'weight',
                    'sampling_frequency'):
            self.assertIn(col, header)

    def test_deterministic_cli(self):
        prefix2 = os.path.join(self.tmp, 'run2', 'community_3')
        rc, _ = _run_cli(['--partitions-dir', COMMUNITY3_DIR,
                          '--bed', COMMUNITY3_BED,
                          '--output', prefix2,
                          '--n-samples', '200', '--seed', '42'])
        self.assertEqual(rc, 0)
        for suffix in ('.traversal.json', '.ml.fa', '.coverage.tsv',
                       '.consensus.fa'):
            with open(self.prefix + suffix) as a, open(prefix2 + suffix) as b:
                self.assertEqual(a.read(), b.read(),
                                 f'{suffix} not deterministic')


@unittest.skipUnless(HAS_COMMUNITY3, 'community_3 inputs not present')
class Community3AncestralIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='c3anc_')
        cls.prefix = os.path.join(cls.tmp, 'community_3')
        rc, out = _run_cli(['--partitions-dir', COMMUNITY3_DIR,
                            '--bed', COMMUNITY3_BED,
                            '--output', cls.prefix,
                            '--mode', 'ancestral',
                            '--n-samples', '30', '--seed', '42'])
        cls.rc = rc
        cls.out = out

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cli_runs(self):
        self.assertEqual(self.rc, 0, self.out)
        for suffix in ('.traversal.json', '.consensus.fa', '.ml.fa',
                       '.ancestral.fa', '.ancestral.genome.fa', '.trees.nwk',
                       '.coverage.tsv', '.stats.json'):
            self.assertTrue(os.path.exists(self.prefix + suffix),
                            f'missing {suffix}')

    def test_ancestral_genome_phage_typical_length(self):
        with open(self.prefix + '.traversal.json') as f:
            data = json.load(f)
        ln = data['genome']['length_bp']
        self.assertGreaterEqual(ln, 10_000)
        self.assertLessEqual(ln, tp.DEFAULT_MAX_LENGTH)

    def test_trees_emitted(self):
        with open(self.prefix + '.trees.nwk') as f:
            lines = [l for l in f.read().splitlines() if l.strip()]
        self.assertGreater(len(lines), 0)
        for line in lines:
            pid, newick = line.split('\t')
            self.assertTrue(newick.endswith(';'))
            # single-member blocks have a bare-leaf newick; multi-member
            # blocks have at least one parenthesized clade
            if '(' in newick:
                self.assertGreaterEqual(newick.count('('), 1)
            else:
                self.assertEqual(newick.count(';'), 1)

    def test_ancestral_seq_per_partition(self):
        anc = _read_fasta(self.prefix + '.ancestral.fa')
        cons = _read_fasta(self.prefix + '.consensus.fa')
        self.assertGreater(len(anc), 0)
        self.assertEqual(sorted(anc), sorted(cons))  # same partition set


# ─── integration: synthetic community (per-partition BED dir) ──────────────

class SyntheticIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='synth_')
        make_synthetic_community(cls.tmp)
        cls.args = ['--partitions-dir', os.path.join(cls.tmp, 'partitions'),
                    '--bed', os.path.join(cls.tmp, 'partitions_bed')]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_ml_mode_all_outputs(self):
        prefix = os.path.join(self.tmp, 'syn_ml')
        rc, out = _run_cli(self.args + ['--output', prefix, '--n-samples', '50'])
        self.assertEqual(rc, 0, out)
        for suffix in ('.traversal.json', '.consensus.fa', '.ml.fa',
                       '.coverage.tsv', '.stats.json'):
            self.assertTrue(os.path.exists(prefix + suffix), f'missing {suffix}')
        with open(prefix + '.traversal.json') as f:
            data = json.load(f)
        self.assertGreaterEqual(data['genome']['length_bp'], 1_000)
        self.assertLessEqual(data['genome']['length_bp'], tp.DEFAULT_MAX_LENGTH)

    def test_ancestral_mode_all_outputs(self):
        prefix = os.path.join(self.tmp, 'syn_anc')
        rc, out = _run_cli(self.args + ['--output', prefix, '--mode', 'ancestral',
                                        '--n-samples', '20'])
        self.assertEqual(rc, 0, out)
        for suffix in ('.traversal.json', '.ancestral.fa',
                       '.ancestral.genome.fa', '.trees.nwk'):
            self.assertTrue(os.path.exists(prefix + suffix), f'missing {suffix}')
        with open(prefix + '.ancestral.genome.fa') as f:
            seq = ''.join(l.strip() for l in f if not l.startswith('>'))
        self.assertGreater(len(seq), 0)
        self.assertLessEqual(len(seq), tp.DEFAULT_MAX_LENGTH)

    def test_fasta_fallback_partitions_flagged(self):
        prefix = os.path.join(self.tmp, 'syn_fb')
        rc, out = _run_cli(self.args + ['--output', prefix])
        self.assertEqual(rc, 0, out)
        cons = _read_fasta(prefix + '.consensus.fa')
        self.assertIn('partition29', cons)
        self.assertIn('partition30', cons)
        with open(prefix + '.traversal.json') as f:
            data = json.load(f)
        self.assertEqual(data['partitions']['29']['seq_source'],
                         'fasta_longest')
        self.assertEqual(data['partitions']['30']['seq_source'],
                         'fasta_longest')

    def test_deterministic_cli(self):
        p1 = os.path.join(self.tmp, 'det1', 'c')
        p2 = os.path.join(self.tmp, 'det2', 'c')
        rc1, out1 = _run_cli(self.args + ['--output', p1, '--seed', '99'])
        rc2, out2 = _run_cli(self.args + ['--output', p2, '--seed', '99'])
        self.assertEqual(rc1, 0, out1)
        self.assertEqual(rc2, 0, out2)
        for suffix in ('.traversal.json', '.ml.fa', '.coverage.tsv'):
            with open(p1 + suffix) as a, open(p2 + suffix) as b:
                self.assertEqual(a.read(), b.read())

    def test_alpha_flag_accepted(self):
        prefix = os.path.join(self.tmp, 'syn_a')
        rc, out = _run_cli(self.args + ['--output', prefix, '--alpha', '2.0'])
        self.assertEqual(rc, 0, out)


if __name__ == '__main__':
    unittest.main(verbosity=2)
