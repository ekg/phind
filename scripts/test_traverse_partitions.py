#!/usr/bin/env python3
"""Tests for scripts/traverse_partitions.py (README §6 optimal typical traversal).

Run:  python3 scripts/test_traverse_partitions.py
      (or: pytest scripts/test_traverse_partitions.py)

Covers:
  * unit tests for BED parsing, adjacency graph, consensus (MAJORITY rule +
    coverage gate), longest-sequence FASTA fallback, rare-partition flags,
    k-mer-Jaccard bridging and deterministic traversal.
  * integration: the archived community_3 partition set reproduces the
    validated stitch reference (research/stitching/repro/community_3_stitched_mean.fa)
    byte-identically for the shared partition consensus and extends it into a
    phage-typical stitched genome; coverage statistics flag rare partitions.
  * the alignment-task input layout (per-partition 3-col BED dir + MAF/FASTA)
    on a synthetic community, including the no-MAF longest-seq fallback and
    genome length within the phage-typical budget.
"""

import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import traverse_partitions as tp  # noqa: E402

COMMUNITY3_BED = os.path.join(REPO, 'research/stitching/community_3_partitions.bed')
COMMUNITY3_DIR = os.path.join(REPO, 'research/stitching/inputs/community_3_partitions')
COMMUNITY3_REF = os.path.join(REPO, 'research/stitching/repro/community_3_stitched_mean.fa')

HAS_COMMUNITY3 = all(os.path.exists(p) for p in
                     [COMMUNITY3_BED, COMMUNITY3_DIR, COMMUNITY3_REF])


# ─── synthetic fixtures ────────────────────────────────────────────────────

def make_synthetic_community(root):
    """Create a synthetic community mimicking `impg partition` output.

    Layout (the alignment task's per-clade format):
      root/partitions_bed/partition<N>.bed   (prophage start end, 3 columns)
      root/partitions/partition<N>.maf       (alignment block) or
      root/partitions/partition<N>.fa        (member FASTA bundle, no MAF)

    Genome model: prophages carry an ordered partition list; core partitions
    (1-5) are shared by most genomes, medium (6-7) by a subset, rare (8-15)
    by a single genome each.  MAF blocks are staggered segment bundles so
    consensus is a union mosaic (as in real impg output); partitions 14/15
    have only FASTA (no MAF) to exercise the longest-seq fallback.
    """
    bed_dir = os.path.join(root, 'partitions_bed')
    part_dir = os.path.join(root, 'partitions')
    os.makedirs(bed_dir)
    os.makedirs(part_dir)

    rng = random.Random(7)
    order = list(range(1, 16))  # canonical typical order 1..15
    genomes = {}
    for g in range(1, 13):
        name = f'G{g:02d}'
        if g <= 9:            # 9/12 carry the core
            pids = [1, 2, 3, 4, 5]
        else:                 # 3/12 carry only part of the core (rare-ish)
            pids = [2, 3]
        if g % 2 == 0:        # half carry medium partitions
            pids += [6, 7]
        if g <= 2:            # two genomes carry rare modules
            pids += [8 + g - 1, 8 + g]  # 8,9 / 10,11
        if g == 12:
            pids += [14, 15]  # FASTA-only rare partitions
        if g == 11:
            pids += [12, 13]
        genomes[name] = pids

    # per-partition 3-col BEDs; intervals are 1000 bp windows laid out in the
    # prophage's partition order
    part_members = {}
    for name, pids in genomes.items():
        start = 0
        for pid in pids:
            part_members.setdefault(pid, []).append((name, start))
            with open(os.path.join(bed_dir, f'partition{pid}.bed'), 'a') as f:
                f.write(f'{name}\t{start}\t{start + 1000}\n')
            start += 1000

    # sequence data
    for pid in sorted(part_members):
        members = part_members[pid]
        if pid in (14, 15):
            # FASTA bundle only (no MAF) -> longest-seq fallback path
            with open(os.path.join(part_dir, f'partition{pid}.fa'), 'w') as f:
                for name, _ in members:
                    f.write(f'>{name}\n{_rand_dna(rng, 1000)}\n')
            continue
        # MAF alignment block: staggered segment bundle (each member a row),
        # with a 400 bp left-flank on rows > 0 so the block is longer than any
        # single member (union mosaic), mirroring real impg output.
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


def _rand_dna(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _run_cli(args, cwd=REPO):
    """Run the CLI in a subprocess; return (returncode, stdout)."""
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


# ─── unit tests ────────────────────────────────────────────────────────────

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


class AdjacencyTests(unittest.TestCase):
    def test_adjacency_graph_counts(self):
        normalized = {
            'G1': [1, 2, 3],
            'G2': [1, 2, 3],
            'G3': [1, 3],
            'G4': [2, 4],
        }
        adj, cooc, first, last, occ = tp.build_adjacency_graph(normalized)
        self.assertEqual(adj[1][2], 2)     # G1,G2
        self.assertEqual(adj[2][3], 2)
        self.assertEqual(adj[1][3], 1)     # G3
        self.assertEqual(adj[2][4], 1)     # G4
        self.assertEqual(occ[1], 3)        # G1,G2,G3
        self.assertEqual(first[1], 3)      # G1,G2,G3
        self.assertEqual(last[3], 3)       # G1,G2,G3
        self.assertEqual(cooc[(1, 2)], 2)  # G1,G2 (unordered pair)


class ConsensusTests(unittest.TestCase):
    def _records(self, rows):
        return [{'name': n, 'start': 0, 'size': len(s), 'strand': '+',
                 'src_size': len(s), 'seq': s} for n, s in rows]

    def test_majority_rule(self):
        recs = self._records([('a', 'ACGT'), ('b', 'ACGA'), ('c', 'ACGC')])
        cons, diag = tp.compute_partition_consensus(recs, 0.0)
        self.assertEqual(cons, 'ACGT')     # majority per column
        self.assertEqual(diag['n_seqs'], 3)

    def test_coverage_threshold_gate(self):
        # columns 2/3 are covered by only 1/3 sequences -> dropped at threshold 0.5
        recs = self._records([('a', 'AA--'), ('b', 'AG--'), ('c', 'ATCC')])
        cons, _ = tp.compute_partition_consensus(recs, 0.5)
        # ceil(3*0.5)=2: col0 'A' (3 bases); col1 is a 3-way tie A/G/T, broken
        # by first-encountered -> 'A' (same rule as the validated
        # stitch_algorithm.py, commit 2363ece); cols 2/3 have 1 base each ->
        # dropped by the coverage gate
        self.assertEqual(cons, 'AA')

    def test_longest_seq_fallback_without_maf(self):
        seq_data = {'fasta': [('a', 'ACGT' * 2), ('b', 'ACGT' * 5)]}
        rep, diag = tp.partition_representative(seq_data, 0.0)
        self.assertEqual(rep, 'ACGT' * 5)
        self.assertEqual(diag['source'], 'fasta_longest')

    def test_maf_preferred_over_fasta(self):
        maf = self._records([('a', 'ACGT')])
        seq_data = {'maf': maf, 'fasta': [('a', 'TTTTTT')]}
        rep, diag = tp.partition_representative(seq_data, 0.0)
        self.assertEqual(rep, 'ACGT')
        self.assertEqual(diag['source'], 'maf_consensus')


class CoverageStatsTests(unittest.TestCase):
    def test_rare_flags(self):
        normalized = {f'G{i}': [1, 2] for i in range(10)}
        normalized['G1'] = [1, 2, 3]          # partition 3 in 1/10
        normalized['G2'] = [1, 2, 4]          # partition 4 in 1/10
        stats = tp.compute_coverage_stats(normalized, 10, 0.2, 5)
        self.assertFalse(stats[1]['rare'])    # 10/10
        self.assertFalse(stats[2]['rare'])
        self.assertTrue(stats[3]['rare'])     # 1/10 < 0.2 and < 5
        self.assertTrue(stats[4]['rare'])


class KmerBridgeTests(unittest.TestCase):
    @staticmethod
    def _rand_seq(rng, n):
        return ''.join(rng.choice('ACGT') for _ in range(n))

    def test_identical_sequences_jaccard_one(self):
        rng = random.Random(3)
        s = self._rand_seq(rng, 400)
        seqs = {1: s, 2: s}
        sim, counts = tp.kmer_jaccard_similarity(seqs, k=15)
        self.assertEqual(sim[(1, 2)], 1.0)
        self.assertEqual(counts[1], 400 - 15 + 1)

    def test_revcomp_aware(self):
        rng = random.Random(4)
        s = self._rand_seq(rng, 400)
        seqs = {1: s, 2: tp._revcomp(s)}
        sim, _ = tp.kmer_jaccard_similarity(seqs, k=15)
        self.assertEqual(sim[(1, 2)], 1.0)  # canonical k-mers: revcomp is identical


class TraversalTests(unittest.TestCase):
    def test_deterministic_and_optimal_on_line(self):
        # prophage partition lists are pure chains; the max-weight path must
        # reproduce the canonical order 1..6 (all observed adjacencies).
        normalized = {}
        for g in range(20):
            chain = list(range(1, 7))
            normalized[f'G{g}'] = chain
        adj, cooc, first, last, occ = tp.build_adjacency_graph(normalized)
        sim = {}
        path, weight, W, joins = tp.find_optimal_traversal(
            list(range(1, 7)), adj, cooc, sim, first, last, occ, 20)
        self.assertEqual(path, list(range(1, 7)))
        self.assertEqual([j['type'] for j in joins], ['observed'] * 5)
        m = tp.compute_metrics(normalized, path, path, adj, occ)
        self.assertEqual(m['expected_path_completeness'], 1.0)
        self.assertEqual(m['fraction_prophages_fully_ordered'], 1.0)

    def test_deterministic_two_runs(self):
        normalized = {}
        for g in range(15):
            normalized[f'G{g}'] = [1, 2, 3, 4] if g % 3 else [2, 1, 4, 3]
            if g % 5 == 0:
                normalized[f'G{g}'] += [5]
        adj, cooc, first, last, occ = tp.build_adjacency_graph(normalized)
        seqs = {p: 'ACGT' * 50 for p in range(1, 6)}
        sim, _ = tp.kmer_jaccard_similarity(seqs, k=15)
        kw = dict(adj=adj, cooc=cooc, sim=sim, first_counts=first,
                  last_counts=last, occurrence=occ, n_prophages=15)
        p1, w1, _, _ = tp.find_optimal_traversal(list(range(1, 6)), **kw)
        p2, w2, _, _ = tp.find_optimal_traversal(list(range(1, 6)), **kw)
        self.assertEqual(p1, p2)
        self.assertEqual(w1, w2)


# ─── integration: archived community_3 ─────────────────────────────────────

@unittest.skipUnless(HAS_COMMUNITY3, 'community_3 inputs not present')
class Community3IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='c3_test_')
        cls.prefix = os.path.join(cls.tmp, 'community_3')
        rc, out = _run_cli(['--partitions-dir', COMMUNITY3_DIR,
                            '--bed', COMMUNITY3_BED,
                            '--output', cls.prefix])
        cls.rc = rc
        cls.out = out
        cls.json_path = cls.prefix + '.traversal.json'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cli_runs(self):
        self.assertEqual(self.rc, 0, self.out)
        self.assertTrue(os.path.exists(self.json_path))

    def test_partition51_consensus_byte_identical_to_reference(self):
        """The validated stitch result (repro/community_3_stitched_mean.fa)
        is the partition-51 majority consensus; our consensus must match it
        byte-for-byte."""
        with open(COMMUNITY3_REF) as f:
            ref = ''.join(l.strip() for l in f if not l.startswith('>'))
        cons = _read_fasta(self.prefix + '.consensus.fa')['partition51']
        self.assertEqual(len(cons), len(ref))
        self.assertEqual(cons, ref)

    def test_genome_phage_typical_length(self):
        with open(self.json_path) as f:
            data = json.load(f)
        ln = data['genome_length_bp']
        # phage-typical: tens of kb to ~150 kb budget, NOT multi-Mbp concat
        self.assertGreaterEqual(ln, 30_000)
        self.assertLessEqual(ln, tp.DEFAULT_MAX_LENGTH)
        # the validated core partitions are extended, not dropped
        genome = set(data['genome'])
        self.assertIn(51, genome)
        self.assertIn(48, genome)
        self.assertIn(241, genome)

    def test_coverage_stats_rare_flags_emitted(self):
        with open(self.prefix + '.coverage.tsv') as f:
            lines = f.read().splitlines()
        header = lines[0].split('\t')
        self.assertIn('rare', header)
        self.assertIn('occurrence', header)
        rows = [l.split('\t') for l in lines[1:]]
        rare_rows = [r for r in rows if r[header.index('rare')] == 'yes']
        self.assertGreater(len(rare_rows), 0)
        # every partition carries an occurrence count >= 1
        self.assertTrue(all(int(r[header.index('occurrence')]) >= 1 for r in rows))

    def test_traversal_keeps_rare_partitions_in_ordering(self):
        with open(self.json_path) as f:
            data = json.load(f)
        self.assertEqual(len(data['traversal']),
                         data['n_partitions_total'])
        rare_in_traversal = [p for p in data['traversal']
                             if data['partitions'][str(p)]['rare']]
        self.assertGreater(len(rare_in_traversal), 0)  # rare kept, not dropped


# ─── integration: synthetic second community (per-partition BED dir) ───────

class SyntheticLayoutIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='synth_')
        cls.genomes = make_synthetic_community(cls.tmp)
        cls.prefix = os.path.join(cls.tmp, 'community_synth')
        rc, out = _run_cli(['--partitions-dir',
                            os.path.join(cls.tmp, 'partitions'),
                            '--bed', os.path.join(cls.tmp, 'partitions_bed'),
                            '--output', cls.prefix])
        cls.rc = rc
        cls.out = out

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cli_runs_on_per_partition_bed_dir(self):
        self.assertEqual(self.rc, 0, self.out)
        for suffix in ('.traversal.json', '.consensus.fa', '.fa',
                       '.coverage.tsv', '.stats.json'):
            self.assertTrue(os.path.exists(self.prefix + suffix),
                            f'missing {suffix}')

    def test_genome_length_phage_typical(self):
        with open(self.prefix + '.traversal.json') as f:
            data = json.load(f)
        ln = data['genome_length_bp']
        self.assertGreaterEqual(ln, 10_000)
        self.assertLessEqual(ln, tp.DEFAULT_MAX_LENGTH)
        # the stitched genome contains the shared core partitions
        genome = set(data['genome'])
        self.assertTrue(genome & {1, 2, 3, 4, 5})

    def test_longest_seq_fallback_for_fasta_only_partitions(self):
        cons = _read_fasta(self.prefix + '.consensus.fa')
        self.assertIn('partition14', cons)
        self.assertIn('partition15', cons)
        with open(self.prefix + '.traversal.json') as f:
            data = json.load(f)
        self.assertEqual(data['partitions']['14']['seq_source'],
                         'fasta_longest')
        self.assertEqual(data['partitions']['15']['seq_source'],
                         'fasta_longest')
        # MAF partitions use consensus
        self.assertEqual(data['partitions']['1']['seq_source'],
                         'maf_consensus')

    def test_metrics_computed_and_bounded(self):
        with open(self.prefix + '.traversal.json') as f:
            data = json.load(f)
        m = data['metrics']
        self.assertIn('expected_path_completeness', m)
        self.assertIn('genome_coverage', m)
        self.assertIn('ordered_coverage', m)
        for k in ('expected_path_completeness', 'genome_coverage',
                  'ordered_coverage', 'fraction_prophages_fully_ordered'):
            self.assertGreaterEqual(m[k], 0.0)
            self.assertLessEqual(m[k], 1.0)

    def test_deterministic_cli(self):
        # same basename in a second directory -> byte-identical outputs
        prefix2 = os.path.join(self.tmp, 'run2', 'community_synth')
        rc, _ = _run_cli(['--partitions-dir', os.path.join(self.tmp, 'partitions'),
                          '--bed', os.path.join(self.tmp, 'partitions_bed'),
                          '--output', prefix2])
        self.assertEqual(rc, 0)
        for suffix in ('.traversal.json', '.fa', '.coverage.tsv',
                       '.consensus.fa'):
            with open(self.prefix + suffix) as a, open(prefix2 + suffix) as b:
                self.assertEqual(a.read(), b.read(),
                                 f'{suffix} not deterministic')


if __name__ == '__main__':
    unittest.main(verbosity=2)
