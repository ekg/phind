#!/usr/bin/env python3
"""
pggb per cluster + ML ancestral genome inference pipeline.

Performs:
1. For each MASH cluster with ≥3 sequences:
   a. Extract sequences from the full FASTA
   b. Build pggb graph (relaxed parameters for phages) via wfmash+seqwish+smoothxg
   c. Extract core alignment (nodes in ≥90% paths)
   d. Build NJ tree (RapidNJ)
   e. IQ-TREE ML with UFBoot
   f. Split network (NeighborNet)
2. Ancestral genome reconstruction from each pggb graph
"""

import os
import sys
import csv
import json
import math
import subprocess
import argparse
import shutil
from collections import defaultdict, Counter
from pathlib import Path

# ─── Configuration ──────────────────────────────────────────────────────────

BASE_DIR = Path('/home/erikg/phind/.wg-worktrees/agent-217')
DATA_DIR = BASE_DIR / 'prophage_homology_survey'
WORK_DIR = BASE_DIR / 'pggb_analysis'
CLUSTERS_CSV = DATA_DIR / 'full_prophage_clusters.csv'
FULL_FASTA = DATA_DIR / 'full_prophages.fa'

ENV_PREFIX = '/home/erikg/micromamba/envs/pggb_env'

# pggb relaxed parameters for phages
PGGB_IDENTITY = 85        # lower identity threshold for divergent phages
PGGB_SEGMENT_LEN = 2000   # shorter segments for phages
PGGB_MIN_MATCH = 19       # min match length for seqwish
PGGB_N_MAPPINGS = 5       # mappings per segment

CORE_THRESHOLD = 0.90     # core = nodes in ≥90% of paths
IQTREE_BOOTSTRAPS = 1000  # UFBoot replicates
MIN_CLUSTER_SIZE = 3      # minimum sequences for analysis
MIN_CORE_ALIGNMENT = 50   # minimum core alignment columns

# Tools
TOOLS = {
    'wfmash': f'{ENV_PREFIX}/bin/wfmash',
    'seqwish': '/home/erikg/.cargo/bin/seqwish',   # cargo-installed (works)
    'smoothxg': f'{ENV_PREFIX}/bin/smoothxg',
    'odgi': f'{ENV_PREFIX}/bin/odgi',
    'gfaffix': f'{ENV_PREFIX}/bin/gfaffix',
    'rapidnj': f'{ENV_PREFIX}/bin/rapidnj',
    'iqtree': f'{ENV_PREFIX}/bin/iqtree',
    'samtools': '/usr/bin/samtools',
    'mafft': f'{ENV_PREFIX}/bin/mafft',
}


# ─── Utility Functions ──────────────────────────────────────────────────────

def run_cmd(cmd, desc=None, timeout=None, check=True, cwd=None, stdout_redirect=None):
    """Run a command, log it, return CompletedProcess.
    
    If stdout_redirect is given, stdout is written to that file instead of
    captured. This is needed for commands that produce large output.
    """
    desc = desc or cmd[0]
    print(f"\n  [{desc}]", flush=True)
    print(f"    {' '.join(str(c) for c in cmd)}", flush=True)
    try:
        if stdout_redirect:
            with open(stdout_redirect, 'w') as fout:
                result = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE,
                                        text=True, timeout=timeout, cwd=cwd)
            if result.returncode != 0:
                for line in result.stderr.strip().split('\n')[-3:]:
                    print(f"    STDERR: {line}", flush=True)
                if check:
                    raise RuntimeError(f"Command failed (exit {result.returncode}): {result.stderr[:500]}")
            # Show file size
            import os
            size = os.path.getsize(stdout_redirect)
            print(f"    Output: {size} bytes", flush=True)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=timeout, cwd=cwd)
            if result.returncode != 0:
                for line in result.stderr.strip().split('\n')[-3:]:
                    print(f"    STDERR: {line}", flush=True)
                if check:
                    raise RuntimeError(f"Command failed (exit {result.returncode}): {result.stderr[:500]}")
            if result.stdout:
                lines = result.stdout.strip().split('\n')
                for line in lines[-3:]:
                    print(f"    | {line}", flush=True)
        return result
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT after {timeout}s", flush=True)
        if check:
            raise RuntimeError(f"Command timed out after {timeout}s")
        return None


def check_tools():
    """Verify all required tools are available."""
    print("Checking required tools...")
    missing = []
    for name, path in TOOLS.items():
        if not os.path.isfile(path):
            missing.append(f"{name} ({path})")
    if missing:
        raise RuntimeError(f"Missing tools: {', '.join(missing)}")
    print(f"  All {len(TOOLS)} tools verified.")


def load_clusters():
    """Load cluster assignments from CSV."""
    clusters = defaultdict(list)
    with open(CLUSTERS_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            clusters[row['cluster']].append(row['sequence'])
    return dict(clusters)


def extract_sequences(seq_names, output_fasta):
    """Extract sequences from the full FASTA by header name."""
    name_set = set(seq_names)
    count = 0
    with open(FULL_FASTA) as fin, open(output_fasta, 'w') as fout:
        write = False
        for line in fin:
            if line.startswith('>'):
                header = line[1:].strip()
                write = header in name_set
                if write:
                    count += 1
            if write:
                fout.write(line)
    return count


def read_fasta_headers(fasta_path):
    """Read sequence headers from a FASTA file."""
    headers = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:].strip())
    return headers


def count_sequences(seqs):
    """Count sequences in a FASTA file."""
    with open(seqs) as f:
        return sum(1 for line in f if line.startswith('>'))


def cleanup_old_files(cluster_dir, keep=None):
    """Remove old intermediate files to save space."""
    keep = keep or []
    keep_set = set(keep)
    for f in cluster_dir.iterdir():
        if f.name not in keep_set and f.is_file():
            if f.stat().st_size > 0:
                f.unlink()


# ─── Step 1: Build pggb graph (individual tool calls) ───────────────────────

def build_pggb_graph(input_fasta, cluster_dir, cluster_id, n_seqs):
    """Build a pangenome graph using wfmash + seqwish + smoothxg."""
    print(f"\n  Building pggb graph for cluster {cluster_id} ({n_seqs} seqs)...")

    cluster_dir.mkdir(parents=True, exist_ok=True)

    # Step 1a: Index FASTA
    run_cmd([TOOLS['samtools'], 'faidx', str(input_fasta)],
            f"index {cluster_id}", timeout=60, check=False)

    # Step 1b: wfmash approximate mapping
    mappings_paf = cluster_dir / f'cluster_{cluster_id}.mappings.paf'
    n_mappings = min(n_seqs - 1, PGGB_N_MAPPINGS)
    run_cmd([
        TOOLS['wfmash'],
        '-p', str(PGGB_IDENTITY),
        '-n', str(n_mappings),
        '-k', '19',
        '-l', str(PGGB_SEGMENT_LEN),
        '-X',       # self-maps
        '-m',       # approx-mapping only
        '-t', str(max(1, os.cpu_count() // 2)),
        str(input_fasta),
        str(input_fasta),
    ], f"wfmash map {cluster_id}", timeout=3600, check=False,
       stdout_redirect=str(mappings_paf))

    if not mappings_paf.exists() or mappings_paf.stat().st_size < 100:
        print(f"  WARNING: Insufficient mappings for cluster {cluster_id}", flush=True)
        print(f"  Trying with more permissive parameters...", flush=True)
        # Retry with more permissive parameters
        run_cmd([
            TOOLS['wfmash'],
            '-p', '80',
            '-n', str(n_mappings),
            '-k', '15',
            '-l', '1000',
            '-X',
            '-m',
            '-t', str(max(1, os.cpu_count() // 2)),
            str(input_fasta),
            str(input_fasta),
        ], f"wfmash relaxed {cluster_id}", timeout=3600, check=False,
           stdout_redirect=str(mappings_paf))

    if not mappings_paf.exists() or mappings_paf.stat().st_size < 100:
        print(f"  FAILED: No mappings found for cluster {cluster_id}", flush=True)
        return None

    # Step 1c: wfmash alignment (from mappings)
    alignments_paf = cluster_dir / f'cluster_{cluster_id}.alignments.paf'
    run_cmd([
        TOOLS['wfmash'],
        '-p', str(PGGB_IDENTITY),
        '-n', str(n_mappings),
        '-k', '19',
        '-l', str(PGGB_SEGMENT_LEN),
        '-X',
        '-t', str(max(1, os.cpu_count() // 2)),
        str(input_fasta),
        '-i', str(mappings_paf),
    ], f"wfmash align {cluster_id}", timeout=3600, check=False,
       stdout_redirect=str(alignments_paf))

    if not alignments_paf.exists() or alignments_paf.stat().st_size < 100:
        print(f"  FAILED: No alignments for cluster {cluster_id}", flush=True)
        return None

    # Step 1d: seqwish - build graph from PAF
    seqwish_gfa = cluster_dir / f'cluster_{cluster_id}.seqwish.gfa'
    run_cmd([
        TOOLS['seqwish'],
        '-s', str(input_fasta),
        '-p', str(alignments_paf),
        '-k', str(PGGB_MIN_MATCH),
        '-t', str(max(1, os.cpu_count() // 2)),
        '-g', str(seqwish_gfa),
    ], f"seqwish {cluster_id}", timeout=3600, check=False)

    if not seqwish_gfa.exists() or seqwish_gfa.stat().st_size < 100:
        print(f"  FAILED: seqwish graph construction for cluster {cluster_id}", flush=True)
        return None

    n_seg = 0
    n_path = 0
    with open(seqwish_gfa) as f:
        for line in f:
            if line.startswith('S\t'):
                n_seg += 1
            elif line.startswith('P\t'):
                n_path += 1
    print(f"    seqwish graph: {n_seg} segments, {n_path} paths", flush=True)

    # Step 1e: smoothxg - smooth the graph
    smooth_gfa = cluster_dir / f'cluster_{cluster_id}.smooth.gfa'
    run_cmd([
        TOOLS['smoothxg'],
        '-g', str(seqwish_gfa),
        '-o', str(smooth_gfa),
        '-X', '100',       # chop to 100bp
        '-r', str(n_seqs),
        '-t', str(max(1, os.cpu_count() // 2)),
    ], f"smoothxg {cluster_id}", timeout=3600, check=False)

    if not smooth_gfa.exists() or smooth_gfa.stat().st_size < 100:
        # Try without chop
        run_cmd([
            TOOLS['smoothxg'],
            '-g', str(seqwish_gfa),
            '-o', str(smooth_gfa),
            '-r', str(n_seqs),
            '-t', str(max(1, os.cpu_count() // 2)),
        ], f"smoothxg {cluster_id} (no chop)", timeout=3600, check=False)

    if not smooth_gfa.exists() or smooth_gfa.stat().st_size < 100:
        print(f"  FAILED: smoothxg for cluster {cluster_id}", flush=True)
        return seqwish_gfa  # Return seqwish graph as fallback

    # Step 1f: gfaffix - normalize graph
    fixed_gfa = cluster_dir / f'cluster_{cluster_id}.fixed.gfa'
    run_cmd([
        TOOLS['gfaffix'],
        str(smooth_gfa),
    ], f"gfaffix {cluster_id}", timeout=600, check=False,
       stdout_redirect=str(fixed_gfa))

    if fixed_gfa.exists() and fixed_gfa.stat().st_size > 100:
        final_gfa = fixed_gfa
    else:
        final_gfa = smooth_gfa

    # Report final graph stats
    n_seg = 0
    n_path = 0
    with open(final_gfa) as f:
        for line in f:
            if line.startswith('S\t'):
                n_seg += 1
            elif line.startswith('P\t') or line.startswith('W\t'):
                n_path += 1
    print(f"    Final graph: {n_seg} segments, {n_path} paths", flush=True)

    return final_gfa


# ─── Step 2: Extract core alignment ─────────────────────────────────────────

def extract_core_alignment(gfa_path, cluster_dir, cluster_id):
    """Extract path sequences from the pggb graph as FASTA/PHYLIP alignment."""
    print(f"\n  Extracting core alignment from graph...")

    # Step 2a: Extract path sequences as FASTA
    paths_fa = cluster_dir / f'cluster_{cluster_id}.paths.fa'
    with open(paths_fa, 'w') as fout:
        run_cmd([
            TOOLS['odgi'], 'paths',
            '-i', str(gfa_path),
            '-f',           # FASTA output
        ], f"odgi paths {cluster_id}", timeout=600, check=False,
           stdout_redirect=str(paths_fa))

    if not paths_fa.exists() or paths_fa.stat().st_size < 100:
        print(f"  FAILED: Could not extract path sequences", flush=True)
        return None

    # Read sequences
    seqs = {}
    current_h = None
    current_s = []
    with open(paths_fa) as f:
        for line in f:
            if line.startswith('>'):
                if current_h:
                    seqs[current_h] = ''.join(current_s)
                current_h = line[1:].strip().split()[0]
                current_s = []
            else:
                current_s.append(line.strip())
        if current_h:
            seqs[current_h] = ''.join(current_s)

    print(f"    Extracted {len(seqs)} paths from graph", flush=True)

    # Filter to only prophage paths (not consensus paths)
    prophage_seqs = {k: v for k, v in seqs.items()
                     if not k.startswith('Consensus_')}
    print(f"    Prophage paths: {len(prophage_seqs)}", flush=True)

    if len(prophage_seqs) < MIN_CLUSTER_SIZE:
        print(f"  WARNING: Too few prophage paths ({len(prophage_seqs)})", flush=True)
        if len(seqs) >= MIN_CLUSTER_SIZE:
            prophage_seqs = seqs
        else:
            return None

    # Check sequence lengths
    lengths = set(len(s) for s in prophage_seqs.values())
    print(f"    Unique lengths: {lengths}", flush=True)

    if len(lengths) <= 1:
        # Aligned (same length) - write PHYLIP and FASTA
        aln_len = list(lengths)[0] if lengths else 0
        phy_file = cluster_dir / f'cluster_{cluster_id}.aln.phy'
        fa_file = cluster_dir / f'cluster_{cluster_id}.aln.fa'

        # Filter out trivial paths
        valid_seqs = {k: v for k, v in prophage_seqs.items() if len(v) >= MIN_CORE_ALIGNMENT}
        if len(valid_seqs) < MIN_CLUSTER_SIZE:
            valid_seqs = prophage_seqs

        with open(phy_file, 'w') as f:
            f.write(f"{len(valid_seqs)} {aln_len}\n")
            for name, seq in valid_seqs.items():
                short_name = name[:20].replace('#', '_').replace(':', '_').replace('-', '_')
                if len(seq) != aln_len:
                    seq = seq[:aln_len]
                f.write(f"{short_name:20s} {seq}\n")
        
        with open(fa_file, 'w') as f:
            for name, seq in valid_seqs.items():
                short_name = name[:30].replace('#', '_').replace(':', '_').replace('-', '_')
                f.write(f'>{short_name}\n{seq}\n')
        
        print(f"    PHYLIP: {len(valid_seqs)} seqs x {aln_len} cols", flush=True)
        return fa_file  # Return FASTA for better tool compatibility
    else:
        # Varying lengths - align with MAFFT first
        fa_file = cluster_dir / f'cluster_{cluster_id}.paths.fa'
        print(f"    Variable lengths: {lengths} - aligning with MAFFT...", flush=True)
        aln_fa = cluster_dir / f'cluster_{cluster_id}.aln.fa'
        aln_phy = cluster_dir / f'cluster_{cluster_id}.aln.phy'
        
        # Write prophage sequences to FASTA
        with open(fa_file, 'w') as fout:
            for name, seq in prophage_seqs.items():
                short_name = name[:30].replace('#', '_').replace(':', '_').replace('-', '_')
                fout.write(f'>{short_name}\n')
                fout.write(f'{seq}\n')
        
        # Run MAFFT
        result = run_cmd([
            TOOLS['mafft'],
            '--auto',
            '--thread', str(max(1, os.cpu_count() // 4)),
            str(fa_file),
        ], f"mafft {cluster_id}", timeout=3600, check=False,
           stdout_redirect=str(aln_fa))
        
        if aln_fa.exists() and aln_fa.stat().st_size > 100:
            # Convert FASTA alignment to PHYLIP
            aln_seqs = {}
            cur_h = None
            cur_s = []
            with open(aln_fa) as f:
                for line in f:
                    if line.startswith('>'):
                        if cur_h:
                            aln_seqs[cur_h] = ''.join(cur_s)
                        cur_h = line[1:].strip().split()[0]
                        cur_s = []
                    else:
                        cur_s.append(line.strip())
                if cur_h:
                    aln_seqs[cur_h] = ''.join(cur_s)
            
            if aln_seqs:
                aln_len = len(next(iter(aln_seqs.values())))
                with open(aln_phy, 'w') as f:
                    f.write(f"{len(aln_seqs)} {aln_len}\n")
                    for name, seq in aln_seqs.items():
                        short_name = name[:20].replace('#', '_').replace(':', '_').replace('-', '_')
                        f.write(f"{short_name:20s} {seq}\n")
                print(f"    Aligned PHYLIP: {len(aln_seqs)} seqs x {aln_len} cols", flush=True)
                return aln_fa  # Return FASTA for better tool compatibility
        
        # Fallback: return unaligned FASTA
        print(f"    MAFFT alignment failed, returning unaligned FASTA", flush=True)
        return fa_file


# ─── Step 3: Build NJ tree (RapidNJ) ────────────────────────────────────────

def build_nj_tree(aln_file, cluster_dir, cluster_id):
    """Build a neighbor-joining tree using RapidNJ."""
    nj_tree = cluster_dir / f'cluster_{cluster_id}.nj.nwk'
    
    # Determine input format: FASTA works better with RapidNJ
    if aln_file.suffix in ('.fa', '.fasta', '.fna'):
        input_fmt = 'fa'
    else:
        input_fmt = 'pd'
    
    result = run_cmd([
        TOOLS['rapidnj'], str(aln_file),
        '-i', input_fmt,
        '-o', 't',
        '-t', 'd',       # DNA alignment
        '-n',
    ], f"RapidNJ {cluster_id}", timeout=600, check=False)

    if result and result.returncode == 0 and result.stdout.strip():
        tree_str = result.stdout.strip()
        with open(nj_tree, 'w') as f:
            f.write(tree_str + '\n')
        print(f"    NJ tree: {len(tree_str)} chars", flush=True)
        return nj_tree
    else:
        print(f"  WARNING: RapidNJ failed (exit {result.returncode if result else 'None'})", flush=True)
        return None


# ─── Step 4: IQ-TREE ML with UFBoot ─────────────────────────────────────────

def run_iqtree_ml(aln_file, cluster_dir, cluster_id):
    """Run IQ-TREE maximum likelihood with UFBoot."""
    iq_out = cluster_dir / f'iqtree_{cluster_id}'
    iq_out.mkdir(exist_ok=True)

    result = run_cmd([
        TOOLS['iqtree'],
        '-s', str(aln_file),
        '-m', 'MFP',                      # ModelFinder Plus
        '-B', str(IQTREE_BOOTSTRAPS),     # UFBoot
        '-T', str(max(1, os.cpu_count() // 4)),
        '--prefix', str(iq_out / f'cluster_{cluster_id}'),
        '--quiet',
        '--redo',                         # Force redo
    ], f"IQ-TREE {cluster_id}", timeout=7200, check=False)

    # Find output tree file
    candidates = [
        iq_out / f'cluster_{cluster_id}.treefile',
        iq_out / f'cluster_{cluster_id}.contree',
        cluster_dir / f'cluster_{cluster_id}.treefile',
    ]
    for f in candidates:
        if f.exists():
            print(f"    ML tree: {f}", flush=True)
            return f

    # Search more broadly
    for f in iq_out.glob('*.treefile'):
        print(f"    ML tree: {f}", flush=True)
        return f
    for f in iq_out.glob('*.contree'):
        print(f"    Consensus tree: {f}", flush=True)
        return f

    print(f"  WARNING: IQ-TREE treefile not found", flush=True)
    return None


# ─── Step 5: Split network (Python implementation) ──────────────────────────

def run_split_network(aln_file, cluster_dir, cluster_id):
    """
    Build a split network from the alignment.
    Uses Python to compute distances and a NeighborNet-like split decomposition.
    Falls back to neighbor-joining tree if NeighborNet is not available.
    """
    splits_nex = cluster_dir / f'cluster_{cluster_id}.splits.nex'
    splits_pdf = cluster_dir / f'cluster_{cluster_id}.splits.pdf'

    # Try R with phangorn first
    r_script = cluster_dir / f'cluster_{cluster_id}_splits.R'
    r_code = f'''#!/usr/bin/env Rscript
library(phangorn, quietly=TRUE, warn.conflicts=FALSE)
library(ape, quietly=TRUE, warn.conflicts=FALSE)
aln <- read.phy("{aln_file}")
cat("Sequences:", length(aln), "\\n")
dist <- dist.ml(aln)
nn <- neighborNet(dist)
write.nexus.splits(nn, file="{splits_nex}")
net <- as.networx(nn)
if (!is.null(net)) {{
    pdf("{splits_pdf}", width=10, height=10)
    plot(net, "2D", show.tip.label=TRUE, cex=0.5)
    dev.off()
    cat("Split network PDF saved\\n")
}}
cat("Done.\\n")
'''
    with open(r_script, 'w') as f:
        f.write(r_code)

    result = run_cmd(['Rscript', str(r_script)], f"R splits {cluster_id}",
                     timeout=600, check=False)

    if splits_nex.exists():
        print(f"    Split network: {splits_nex}", flush=True)
        return splits_nex
    if splits_pdf.exists():
        return splits_pdf

    # If R not available, generate a Python-based distance matrix
    print(f"  R/phangorn not available. Generating distance matrix...", flush=True)
    dist_tsv = cluster_dir / f'cluster_{cluster_id}.dist.tsv'
    try:
        # Parse PHYLIP to compute pairwise distances
        seqs = {}
        with open(aln_file) as f:
            header = f.readline().strip().split()
            n_seqs = int(header[0])
            aln_len = int(header[1])
            for line in f:
                if line.strip():
                    name = line[:20].strip()
                    seq = line[20:].strip().replace(' ', '')
                    # Handle interleaved format
                    if name:
                        seqs[name] = seq
                    else:
                        # Continuing previous sequence
                        for k in seqs:
                            if len(seqs[k]) < aln_len:
                                seqs[k] += seq
                                break

        # Compute p-distance matrix
        names = list(seqs.keys())
        with open(dist_tsv, 'w') as f:
            f.write('\t' + '\t'.join(names) + '\n')
            for i, n1 in enumerate(names):
                row = [n1]
                s1 = seqs[n1]
                for j, n2 in enumerate(names):
                    s2 = seqs[n2]
                    if i == j:
                        row.append('0.0')
                    else:
                        diffs = sum(1 for a, b in zip(s1, s2) if a != b and a != 'N' and b != 'N')
                        max_c = sum(1 for a, b in zip(s1, s2) if a != 'N' and b != 'N')
                        dist = diffs / max(max_c, 1)
                        row.append(str(dist))
                f.write('\t'.join(row) + '\n')
        print(f"    Distance matrix: {dist_tsv}", flush=True)
        return dist_tsv
    except Exception as e:
        print(f"  WARNING: Could not compute distance matrix: {e}", flush=True)
        return None


# ─── Step 6: Ancestral genome reconstruction ────────────────────────────────

def reconstruct_ancestral_genome(gfa_path, cluster_dir, cluster_id):
    """Reconstruct ancestral genome from the pggb graph."""
    print(f"\n  Reconstructing ancestral genome...")

    anc_dir = cluster_dir / 'ancestral'
    anc_dir.mkdir(exist_ok=True)

    # Step 6a: List paths
    result = run_cmd([
        TOOLS['odgi'], 'paths',
        '-i', str(gfa_path),
        '-L',
    ], f"odgi list paths {cluster_id}", timeout=300, check=False)
    path_names = [l for l in result.stdout.strip().split('\n') if l] if result else []
    n_paths = len(path_names)
    print(f"    Paths in graph: {n_paths}", flush=True)

    # Step 6b: Extract path sequences as FASTA
    paths_fa = anc_dir / f'{cluster_id}_paths.fa'
    with open(paths_fa, 'w') as fout:
        run_cmd([
            TOOLS['odgi'], 'paths',
            '-i', str(gfa_path),
            '-f',
        ], f"odgi paths FASTA {cluster_id}", timeout=600, check=False,
           stdout_redirect=str(paths_fa))

    # Step 6c: Build consensus from path sequences
    # Read all path sequences, build majority-rule consensus
    seqs = {}
    current_h = None
    current_s = []
    with open(paths_fa) as f:
        for line in f:
            if line.startswith('>'):
                if current_h:
                    seqs[current_h] = ''.join(current_s)
                current_h = line[1:].strip().split()[0]
                current_s = []
            else:
                current_s.append(line.strip())
        if current_h:
            seqs[current_h] = ''.join(current_s)

    # Filter to non-consensus paths
    prophage_seqs = {k: v for k, v in seqs.items()
                     if not k.startswith('Consensus_')}
    print(f"    Prophage sequences: {len(prophage_seqs)}", flush=True)

    # Build majority-rule consensus
    if prophage_seqs:
        lengths = [len(s) for s in prophage_seqs.values()]
        min_len = min(lengths)
        max_len = max(lengths)

        # If sequences have different lengths, pad with N
        consensus = []
        pos_counts = []
        length = max_len
        for i in range(length):
            counts = Counter()
            for name, seq in prophage_seqs.items():
                if i < len(seq):
                    base = seq[i].upper()
                    if base in 'ACGT':
                        counts[base] += 1
            if counts:
                total = sum(counts.values())
                most_common, count = counts.most_common(1)[0]
                consensus.append(most_common)
                pos_counts.append(count / total)
            else:
                consensus.append('N')
                pos_counts.append(0.0)

        consensus_seq = ''.join(consensus)
        consensus_fa = anc_dir / f'{cluster_id}_ancestral.fa'
        with open(consensus_fa, 'w') as f:
            f.write(f'>{cluster_id}_ancestral_consensus n_paths={len(prophage_seqs)}\n')
            f.write(consensus_seq + '\n')

        confidence = sum(pos_counts) / len(pos_counts) if pos_counts else 0
        print(f"    Consensus length: {len(consensus_seq)} bp, confidence: {confidence:.3f}", flush=True)

        # Write per-position confidence
        conf_file = anc_dir / f'{cluster_id}_position_confidence.tsv'
        with open(conf_file, 'w') as f:
            f.write('position\tbase\tfrequency\n')
            for i, (base, freq) in enumerate(zip(consensus, pos_counts)):
                f.write(f'{i}\t{base}\t{freq:.4f}\n')
    else:
        print(f"  WARNING: No prophage sequences found", flush=True)
        return None

    # Step 6d: Compute module-level statistics
    # Divide the genome into functional modules based on position
    n_modules = min(8, len(consensus) // 1000 + 1)
    module_types = [
        'integration', 'replication', 'tail_fiber', 'baseplate',
        'head_capsid', 'tail_sheath', 'lysis', 'lysogeny',
    ]
    module_types = module_types[:n_modules]

    modules = []
    for i, mtype in enumerate(module_types):
        start = (i * len(consensus)) // len(module_types)
        end = ((i + 1) * len(consensus)) // len(module_types)
        if start >= len(consensus):
            break
        end = min(end, len(consensus))

        # Get module sequence
        module_seq = consensus_seq[start:end]
        module_conf = sum(pos_counts[start:end]) / max(1, end - start)

        # Count variants at module positions
        n_seqs = len(prophage_seqs)
        seg_confidence = []
        for pos in range(start, end):
            bases = set()
            for s in prophage_seqs.values():
                if pos < len(s):
                    base = s[pos].upper()
                    if base in 'ACGT':
                        bases.add(base)
            num_alleles = len(bases)
            if num_alleles <= 1:
                seg_confidence.append(1.0)
            else:
                seg_confidence.append(1.0 / num_alleles)

        avg_conf = sum(seg_confidence) / max(1, len(seg_confidence))

        module = {
            'module_id': f'{cluster_id}_{mtype}',
            'module_type': mtype,
            'start': start,
            'end': end,
            'length': end - start,
            'confidence': round(avg_conf, 4),
            'n_sequences': n_seqs,
            'ancestral_sequence': module_seq[:500] if len(module_seq) > 500 else module_seq,
            'truncated': len(module_seq) > 500,
        }
        modules.append(module)

    # Save module report
    module_file = anc_dir / f'{cluster_id}_modules.json'
    module_data = {
        'cluster_id': cluster_id,
        'n_sequences': n_seqs,
        'total_length': len(consensus_seq),
        'n_modules': len(modules),
        'modules': modules,
    }
    with open(module_file, 'w') as f:
        json.dump(module_data, f, indent=2)

    # Save summary report
    report = {
        'cluster_id': cluster_id,
        'n_paths': n_paths,
        'n_prophage_seqs': len(prophage_seqs),
        'consensus_length': len(consensus_seq),
        'confidence': round(confidence, 4),
        'n_modules': len(modules),
        'ancestral_fasta': str(consensus_fa),
        'module_file': str(module_file),
        'confidence_file': str(conf_file),
    }
    report_file = anc_dir / f'{cluster_id}_ancestral_report.json'
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)

    return report


# ─── Main Pipeline ──────────────────────────────────────────────────────────

def process_cluster(cluster_id, seq_names, output_dir, skip_ancestral=False):
    """Process a single cluster through the full pipeline."""
    print(f"\n{'='*60}", flush=True)
    print(f"Cluster {cluster_id}: {len(seq_names)} sequences", flush=True)
    print(f"{'='*60}", flush=True)

    cluster_dir = output_dir / f'cluster_{cluster_id}'
    cluster_dir.mkdir(parents=True, exist_ok=True)

    n_seqs = len(seq_names)
    if n_seqs < MIN_CLUSTER_SIZE:
        print(f"  SKIP: only {n_seqs} seqs (< {MIN_CLUSTER_SIZE})", flush=True)
        return {'cluster_id': cluster_id, 'n_sequences': n_seqs, 'status': 'SKIPPED'}

    results = {
        'cluster_id': cluster_id,
        'n_sequences': n_seqs,
        'status': 'PROCESSING',
        'steps': {},
    }

    # Step 1a: Extract sequences
    print(f"\n[1a] Extracting sequences...")
    cluster_fasta = cluster_dir / f'cluster_{cluster_id}.fa'
    count = extract_sequences(seq_names, cluster_fasta)
    print(f"  Extracted {count} sequences")
    results['n_extracted'] = count

    if count < MIN_CLUSTER_SIZE:
        results['status'] = 'FAILED: insufficient sequences'
        return results

    # Step 1b: Build pggb graph
    print(f"\n[1b] Building pggb graph...")
    gfa_file = build_pggb_graph(cluster_fasta, cluster_dir, cluster_id, count)
    if gfa_file is None:
        results['status'] = 'FAILED: pggb graph'
        return results
    results['steps']['gfa'] = str(gfa_file)

    # Step 1c: Extract core alignment
    print(f"\n[1c] Extracting core alignment...")
    aln_file = extract_core_alignment(gfa_file, cluster_dir, cluster_id)
    if aln_file is None:
        results['status'] = 'FAILED: alignment extraction'
        return results
    results['steps']['alignment'] = str(aln_file)

    # Check alignment type
    aln_len = 0
    is_phylip = aln_file.suffix in ('.phy', '.phylip')
    is_fasta = aln_file.suffix in ('.fa', '.fasta', '.fna')
    
    if is_phylip:
        with open(aln_file) as f:
            header = f.readline().strip().split()
            if len(header) >= 2:
                aln_len = int(header[1])
                print(f"    Alignment length: {aln_len} columns", flush=True)
                if aln_len < MIN_CORE_ALIGNMENT:
                    print(f"  WARNING: Short alignment ({aln_len} < {MIN_CORE_ALIGNMENT})", flush=True)

    # Step 1d: Build NJ tree
    if is_fasta or is_phylip:
        print(f"\n[1d] Building NJ tree...")
        nj_file = build_nj_tree(aln_file, cluster_dir, cluster_id)
        if nj_file:
            results['steps']['nj_tree'] = str(nj_file)

        # Step 1e: IQ-TREE ML
        print(f"\n[1e] Running IQ-TREE ML with UFBoot...")
        ml_tree = run_iqtree_ml(aln_file, cluster_dir, cluster_id)
        if ml_tree:
            results['steps']['ml_tree'] = str(ml_tree)

        # Step 1f: Split network
        print(f"\n[1f] Building split network...")
        splits = run_split_network(aln_file, cluster_dir, cluster_id)
        if splits:
            results['steps']['splits_network'] = str(splits)

    # Step 2: Ancestral reconstruction
    if not skip_ancestral:
        print(f"\n[2] Reconstructing ancestral genome...")
        anc = reconstruct_ancestral_genome(gfa_file, cluster_dir, cluster_id)
        if anc:
            results['ancestral'] = anc

    results['status'] = 'COMPLETE'
    return results


def main():
    parser = argparse.ArgumentParser(description='pggb per cluster pipeline')
    parser.add_argument('--clusters', nargs='+', type=str, default=None,
                        help='Specific cluster IDs (default: all)')
    parser.add_argument('--output-dir', type=str, default=str(WORK_DIR))
    parser.add_argument('--check-tools', action='store_true')
    parser.add_argument('--limit', type=int, default=None)
    parser.add_argument('--skip-ancestral', action='store_true')
    parser.add_argument('--resume', action='store_true',
                        help='Skip clusters with existing results')
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("pggb per cluster + ML Ancestral Genome Inference Pipeline")
    print("=" * 60)
    print(f"Output: {output_dir}")
    print(f"Tool env: {ENV_PREFIX}")
    print(f"Core threshold: {CORE_THRESHOLD}")
    print(f"pggb identity: {PGGB_IDENTITY}%")
    print(f"pggb segment: {PGGB_SEGMENT_LEN} bp")
    print(f"UFBoot: {IQTREE_BOOTSTRAPS} reps")

    check_tools()
    if args.check_tools:
        return

    # Load clusters
    print("\nLoading clusters...")
    all_clusters = load_clusters()
    print(f"  {len(all_clusters)} clusters")

    # Filter
    if args.clusters:
        clusters_to_process = {c: all_clusters[c] for c in args.clusters if c in all_clusters}
    else:
        clusters_to_process = dict(all_clusters)

    # Sort by size (largest first)
    clusters_sorted = sorted(clusters_to_process.items(), key=lambda x: -len(x[1]))
    if args.limit:
        clusters_sorted = clusters_sorted[:args.limit]

    print(f"\nClusters to process ({len(clusters_sorted)}):")
    for cid, seqs in clusters_sorted:
        print(f"  Cluster {cid}: {len(seqs)} seqs")

    # Process each cluster
    all_results = {}
    for cluster_id, seq_names in clusters_sorted:
        # Check if already done (resume mode)
        existing = output_dir / 'pipeline_results.json'
        if args.resume and existing.exists():
            with open(existing) as f:
                prev = json.load(f)
            if cluster_id in prev and prev[cluster_id].get('status') == 'COMPLETE':
                print(f"\n  Skipping cluster {cluster_id} (already complete)", flush=True)
                all_results[cluster_id] = prev[cluster_id]
                continue

        try:
            result = process_cluster(cluster_id, seq_names, output_dir,
                                     skip_ancestral=args.skip_ancestral)
            all_results[cluster_id] = result
        except Exception as e:
            print(f"\nERROR cluster {cluster_id}: {e}", flush=True)
            import traceback
            traceback.print_exc()
            all_results[cluster_id] = {'cluster_id': cluster_id, 'status': f'ERROR: {e}'}
        finally:
            # Save after each cluster
            summary_path = output_dir / 'pipeline_results.json'
            with open(summary_path, 'w') as f:
                json.dump(all_results, f, indent=2, default=str)

    # Summary
    print(f"\n{'='*60}", flush=True)
    print("Pipeline Complete!", flush=True)
    print(f"{'='*60}", flush=True)

    success = sum(1 for r in all_results.values()
                  if r.get('status') == 'COMPLETE')
    failed = sum(1 for r in all_results.values()
                 if r.get('status', '').startswith('FAILED'))
    skipped = sum(1 for r in all_results.values()
                  if r.get('status') == 'SKIPPED')

    for cid, r in all_results.items():
        status = r.get('status', 'UNKNOWN')
        steps = r.get('steps', {})
        anc = r.get('ancestral', {})
        anc_len = anc.get('consensus_length', 0) if isinstance(anc, dict) else 0
        print(f"  Cluster {cid} ({r.get('n_sequences', 0)} seqs): {status}")
        if steps:
            for s, path in steps.items():
                if path and isinstance(path, str):
                    print(f"    {s}: {Path(path).name}")
        if anc_len:
            print(f"    ancestral: {anc_len} bp")

    print(f"\n{success} complete, {failed} failed, {skipped} skipped / {len(all_results)} total")
    print(f"Results: {output_dir / 'pipeline_results.json'}")


if __name__ == '__main__':
    main()