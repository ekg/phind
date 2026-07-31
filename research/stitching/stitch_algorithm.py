#!/usr/bin/env python3
"""
Partition Stitching Algorithm v2

Builds a properly ordered mean genome from partition alignment blocks.
Uses adjacency graph for ordering and median-position-based merging.

Algorithm:
1. For each prophage, record its ordered list of partitions (from BED file)
2. Count adjacencies between consecutive partitions
3. Build a directed graph weighted by adjacency counts
4. Find the maximum likelihood path through the graph
5. Compute consensus sequences for each partition (with coverage threshold)
6. Merge overlapping partitions in the path
7. Generate the mean genome
8. Handle accessories: partitions appearing in <50% of prophages are excluded

Usage:
    python3 stitch_algorithm.py --partition-dir <dir> --bed <bed> --output <output.fa>
    python3 stitch_algorithm.py --partition-dir <dir> --bed <bed> --output <output.fa> --ancestral <ancestral.fa>
"""

import argparse
import glob
import json
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ─── MAF parsing ────────────────────────────────────────────────────────────

def parse_partition_maf(maf_path):
    """Parse a single MAF file containing one alignment block.
    
    Returns list of dicts with name, start, size, strand, src_size, seq.
    """
    records = []
    in_block = False
    with open(maf_path) as f:
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


def get_partition_id(path):
    """Extract partition ID from path."""
    base = os.path.basename(path)
    for ext in ['.maf', '.maf.gz']:
        if base.endswith(ext):
            base = base[:-len(ext)]
    digits = ''.join(c for c in base if c.isdigit())
    return int(digits) if digits else hash(base)


# ─── Consensus computation ──────────────────────────────────────────────────

def compute_partition_consensus(records, coverage_threshold=0.25):
    """Compute majority-rule consensus from aligned sequences.
    
    Only includes positions where at least `coverage_threshold` fraction of
    sequences have a non-gap base.
    
    Returns consensus string or None if no records.
    """
    if not records:
        return None
    
    n_seqs = len(records)
    aln_len = len(records[0]['seq'])
    min_coverage = max(1, int(n_seqs * coverage_threshold))
    
    consensus = []
    for i in range(aln_len):
        counts = Counter()
        total = 0
        for rec in records:
            base = rec['seq'][i].upper()
            if base in 'ACGT':
                counts[base] += 1
                total += 1
        if total >= min_coverage and counts:
            consensus.append(counts.most_common(1)[0][0])
    
    return ''.join(consensus)


# ─── BED parsing ────────────────────────────────────────────────────────────

def parse_bed_file(bed_path):
    """Parse BED file with partition assignments.
    
    Format: prophage_name  start  end  partition_id
    
    Returns dict mapping prophage_name -> [(partition_id, start, end)]
    """
    prophage_partitions = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                name = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                pid = int(parts[3])
                prophage_partitions[name].append((pid, start, end))
    return dict(prophage_partitions)


# ─── Adjacency graph ────────────────────────────────────────────────────────

def build_adjacency_graph(prophage_partitions):
    """Build adjacency graph from prophage partition lists.
    
    Returns:
        adj: dict pid -> {target_pid: count}
        first_counts: Counter of partition appearing first
        last_counts: Counter of partition appearing last
        occurrence_counts: Counter of partition occurrences
    """
    adj = defaultdict(lambda: defaultdict(int))
    first_counts = Counter()
    last_counts = Counter()
    occurrence_counts = Counter()
    
    for name, entries in prophage_partitions.items():
        entries.sort(key=lambda x: x[1])  # Sort by start position
        pids = [e[0] for e in entries]
        if not pids:
            continue
        
        first_counts[pids[0]] += 1
        last_counts[pids[-1]] += 1
        for pid in pids:
            occurrence_counts[pid] += 1
        
        for i in range(len(pids) - 1):
            adj[pids[i]][pids[i + 1]] += 1
    
    return adj, first_counts, last_counts, occurrence_counts


def find_maximum_likelihood_path(adj, first_counts, occurrence_counts, n_prophages,
                                  accessory_threshold=0.5):
    """Find the maximum likelihood path through the partition graph.
    
    Core partitions: appear in >= accessory_threshold fraction of prophages.
    Path: start from partition that appears first most often, then greedily
    follow the highest-weight edge among core partitions.
    
    Returns:
        core_path: list of partition_ids in order
        accessory: set of accessory partition_ids
    """
    # Identify core vs accessory partitions
    core_partitions = set()
    accessory = set()
    
    for pid, count in occurrence_counts.items():
        fraction = count / max(n_prophages, 1)
        if fraction >= accessory_threshold:
            core_partitions.add(pid)
        else:
            accessory.add(pid)
    
    if not core_partitions:
        # Fallback: use all partitions
        core_partitions = set(adj.keys())
        accessory = set()
    
    # Build adjacency restricted to core partitions
    core_adj = {}
    for a in core_partitions:
        targets = {}
        for b, count in adj.get(a, {}).items():
            if b in core_partitions:
                targets[b] = count
        if targets:
            core_adj[a] = targets
    
    # Find start node: partition that appears first most often among core
    start_candidates = [(pid, count) for pid, count in first_counts.items()
                        if pid in core_partitions]
    if not start_candidates:
        start_candidates = [(pid, len(edges)) for pid, edges in core_adj.items()]
    
    if not start_candidates:
        return list(core_partitions), list(accessory)
    
    start_candidates.sort(key=lambda x: -x[1])
    start = start_candidates[0][0]
    
    # Greedy path following
    visited = {start}
    path = [start]
    current = start
    
    while True:
        targets = core_adj.get(current, {})
        unvisited = [(t, w) for t, w in targets.items() if t not in visited]
        if not unvisited:
            break
        unvisited.sort(key=lambda x: -x[1])
        next_node = unvisited[0][0]
        visited.add(next_node)
        path.append(next_node)
        current = next_node
    
    # Add remaining core partitions
    remaining = core_partitions - visited
    if remaining:
        remaining_sorted = sorted(remaining, key=lambda x: occurrence_counts.get(x, 0), reverse=True)
        path.extend(remaining_sorted)
    
    return path, list(accessory)


# ─── Overlap detection and merging ──────────────────────────────────────────

def find_overlap(seq_a, seq_b, min_overlap=50, max_overlap_check=5000):
    """Find the longest overlap between suffix of A and prefix of B.
    
    Uses a fast hash-based approach: compute hashes for all suffixes of A
    and prefixes of B, then find the longest match.
    
    Returns (overlap_len, merged_seq) or (0, None) if no overlap found.
    """
    max_possible = min(len(seq_a), len(seq_b), max_overlap_check)
    
    # Build hash set of suffix hashes for A (limited to max_possible)
    # Use Python's built-in string hashing for speed
    # Check all possible overlap lengths (from long to short)
    for overlap_len in range(max_possible, min_overlap - 1, -1):
        suffix = seq_a[-overlap_len:]
        prefix = seq_b[:overlap_len]
        if suffix == prefix:
            merged = seq_a + seq_b[overlap_len:]
            return overlap_len, merged
        
        # Allow small number of mismatches
        mismatches = 0
        fast_check = True
        for k in range(0, overlap_len, 100):
            chunk_a = suffix[k:k+100]
            chunk_b = prefix[k:k+100]
            if chunk_a != chunk_b:
                fast_check = False
                break
        if fast_check:
            # Full check
            mismatches = sum(1 for a, b in zip(suffix, prefix) if a != b)
            max_mismatches = max(1, overlap_len // 50)
            if mismatches <= max_mismatches:
                merged = seq_a + seq_b[overlap_len:]
                return overlap_len, merged
    
    return 0, None


def stitch_and_merge(path, partition_consensi, min_overlap=50):
    """Stitch partition consensus sequences in order, merging overlaps.
    
    Returns the merged genome sequence.
    """
    if not path:
        return ""
    
    result = partition_consensi.get(path[0], "")
    if not result:
        return ""
    
    for i in range(1, len(path)):
        pid = path[i]
        next_seq = partition_consensi.get(pid, "")
        if not next_seq:
            continue
        
        # Try to find overlap
        overlap_len, merged = find_overlap(result, next_seq, min_overlap)
        if merged:
            result = merged
        else:
            # No overlap: concatenate with a small gap (N's)
            result += 'N' * 10 + next_seq
    
    return result


# ─── Validation ─────────────────────────────────────────────────────────────

def compute_mash_identity(fasta_path, ancestral_path, mash_cmd='mash'):
    """Compute approximate identity between two FASTA files using MASH."""
    import subprocess
    
    # Sketch both files
    for f in [fasta_path, ancestral_path]:
        sketch = f + '.msh'
        if not os.path.exists(sketch):
            result = subprocess.run(
                [mash_cmd, 'sketch', '-i', '-o', sketch.replace('.msh', ''), f],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                return None
    
    # Compute distance
    s1 = fasta_path + '.msh'
    s2 = ancestral_path + '.msh'
    result = subprocess.run(
        [mash_cmd, 'dist', s1, s2],
        capture_output=True, text=True, timeout=300
    )
    if result.returncode == 0 and result.stdout.strip():
        for line in result.stdout.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 3:
                distance = float(parts[2])
                return 1.0 - distance
    return None


# ─── Main pipeline ──────────────────────────────────────────────────────────

def run_stitching(partition_dir, bed_path, output_path, ancestral_path=None,
                  mash_cmd='mash', accessory_threshold=0.5, 
                  coverage_threshold=0.25, min_overlap=50):
    """Run the full stitching pipeline."""
    print("=" * 60)
    print("PARTITION STITCHING ALGORITHM v2")
    print("=" * 60)
    
    t0 = time.time()
    
    # Step 1: Load partitions
    print("\n[1] Loading partitions...")
    maf_files = sorted(glob.glob(os.path.join(partition_dir, '*.maf')))
    print(f"    Found {len(maf_files)} MAF files")
    
    # Step 2: Parse BED file
    print("\n[2] Parsing BED file...")
    prophage_partitions = parse_bed_file(bed_path)
    n_prophages = len(prophage_partitions)
    print(f"    {n_prophages} prophages with partition assignments")
    
    # Step 3: Build adjacency graph
    print("\n[3] Building adjacency graph...")
    adj, first_counts, last_counts, occurrence_counts = build_adjacency_graph(prophage_partitions)
    print(f"    {len(adj)} nodes, {sum(len(t) for t in adj.values())} directed edges")
    
    # Step 4: Find maximum likelihood path
    print("\n[4] Finding maximum likelihood path...")
    core_path, accessory = find_maximum_likelihood_path(
        adj, first_counts, occurrence_counts, n_prophages,
        accessory_threshold=accessory_threshold
    )
    print(f"    Core path: {len(core_path)} partitions")
    print(f"    Accessory: {len(accessory)} partitions (<{accessory_threshold*100:.0f}% occurrence)")
    if core_path:
        print(f"    First 5: {core_path[:5]}")
        print(f"    Last 5:  {core_path[-5:]}")
    
    # Step 5: Compute consensus for each partition
    print("\n[5] Computing partition consensus sequences...")
    partition_consensi = {}
    for pid in set(core_path) | set(accessory):
        maf_path = os.path.join(partition_dir, f'partition{pid}.maf')
        if not os.path.exists(maf_path):
            continue
        records = parse_partition_maf(maf_path)
        consensus = compute_partition_consensus(records, coverage_threshold)
        if consensus:
            partition_consensi[pid] = consensus
    
    print(f"    Computed consensus for {len(partition_consensi)} partitions")
    
    # Step 6: Stitch and merge
    print("\n[6] Stitching and merging...")
    core_genome = stitch_and_merge(core_path, partition_consensi, min_overlap)
    print(f"    Core genome: {len(core_genome):,} bp")
    
    # Step 7: Write output
    print("\n[7] Writing output...")
    with open(output_path, 'w') as f:
        if core_genome:
            n_partitions = len(core_path)
            f.write(f'>community_mean_genome_core path={n_partitions}_partitions\n')
            for i in range(0, len(core_genome), 80):
                f.write(core_genome[i:i+80] + '\n')
    print(f"    Written to {output_path}")
    
    # Step 8: Compare to ancestral genome
    if ancestral_path and os.path.exists(ancestral_path):
        print("\n[8] Comparing to ancestral genome...")
        anc_records = list(SeqIO.parse(ancestral_path, 'fasta'))
        if anc_records:
            anc_len = len(anc_records[0].seq)
            print(f"    Ancestral genome: {anc_len:,} bp")
            print(f"    Core genome:      {len(core_genome):,} bp")
            print(f"    Ratio:            {len(core_genome) / max(anc_len, 1):.3f}")
            
            # Write core for MASH comparison
            import tempfile
            temp = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
            temp.write(f'>stitched_core\n{core_genome}\n')
            temp.close()
            
            identity = compute_mash_identity(temp.name, ancestral_path, mash_cmd)
            if identity is not None:
                print(f"    MASH identity vs ancestral: {identity:.4f} ({identity*100:.2f}%)")
            else:
                print(f"    MASH identity: could not compute")
            
            os.unlink(temp.name)
            for f in [temp.name + '.msh', ancestral_path + '.msh']:
                if os.path.exists(f):
                    try:
                        os.unlink(f)
                    except:
                        pass
    
    t1 = time.time()
    print(f"\n{'='*60}")
    print(f"Completed in {t1-t0:.1f}s")
    print(f"{'='*60}")
    
    return {
        'n_maf_files': len(maf_files),
        'n_prophages': n_prophages,
        'core_path_length': len(core_path),
        'n_accessory': len(accessory),
        'core_genome_length': len(core_genome),
        'ancestral_length': len(list(SeqIO.parse(ancestral_path, 'fasta'))[0].seq)
            if ancestral_path and os.path.exists(ancestral_path) else None,
        'accessory_threshold': accessory_threshold,
        'coverage_threshold': coverage_threshold,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Partition Stitching Algorithm v2: order + merge partitions into a mean genome')
    parser.add_argument('--partition-dir', required=True,
                        help='Directory containing partition MAF files')
    parser.add_argument('--bed', required=True,
                        help='BED file with partition assignments per prophage')
    parser.add_argument('--output', '-o', default='stitched_mean.fa',
                        help='Output FASTA path')
    parser.add_argument('--ancestral', '-a', default=None,
                        help='Ancestral genome FASTA for comparison')
    parser.add_argument('--mash', default='mash',
                        help='MASH command path')
    parser.add_argument('--accessory-threshold', type=float, default=0.5,
                        help='Fraction threshold for accessory vs core (default: 0.5)')
    parser.add_argument('--coverage-threshold', type=float, default=0.25,
                        help='Coverage threshold for consensus (default: 0.25)')
    parser.add_argument('--min-overlap', type=int, default=50,
                        help='Minimum overlap for merging (default: 50)')
    parser.add_argument('--json', default=None,
                        help='Path to save results as JSON')
    
    args = parser.parse_args()
    
    result = run_stitching(
        partition_dir=args.partition_dir,
        bed_path=args.bed,
        output_path=args.output,
        ancestral_path=args.ancestral,
        mash_cmd=args.mash,
        accessory_threshold=args.accessory_threshold,
        coverage_threshold=args.coverage_threshold,
        min_overlap=args.min_overlap,
    )
    
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"Results saved to {args.json}")


if __name__ == '__main__':
    main()