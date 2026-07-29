#!/usr/bin/env python3
"""
Full prophage homology survey: MASH-based analysis + PDF report.

Pipeline:
1. Read MASH edge list (or triangle matrix)
2. Build distance matrix for MDS and hierarchical clustering
3. Compute distance statistics (within/between genome)
4. Analyze sequence length distribution
5. Generate comprehensive PDF report with Typst
"""

import os
import sys
import csv
import math
import struct
import warnings
import subprocess
import tempfile
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.sparse import csr_matrix
from sklearn.manifold import MDS
from sklearn.cluster import DBSCAN
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
DATA_DIR = Path('prophage_homology_survey')
FULL_FASTA = DATA_DIR / 'full_prophages.fa'
SAMPLE_MSH = DATA_DIR / 'full_prophages_sample_10k.msh'
SAMPLE_EDGES = DATA_DIR / 'full_prophages_sample_10k_edges.tsv'
OUTPUT_DIR = DATA_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
MDS_COORDS = OUTPUT_DIR / 'full_prophage_mds_coords.csv'
CLUSTERS = OUTPUT_DIR / 'full_prophage_clusters.csv'
DIST_STATS = OUTPUT_DIR / 'full_prophage_distance_stats.csv'
REPORT_PDF = OUTPUT_DIR / 'full_prophage_homology_report.pdf'
REPORT_TYP = OUTPUT_DIR / 'full_prophage_homology_report.typ'

# Plotting
sns.set_style('whitegrid')
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
})

# ============================================================
# Helper functions
# ============================================================

def parse_sequence_lengths(fasta_path):
    """Parse sequence lengths from FASTA file."""
    seq_lens = {}
    with open(fasta_path, 'r') as f:
        current_header = None
        current_len = 0
        for line in f:
            if line.startswith('>'):
                if current_header is not None:
                    seq_lens[current_header] = current_len
                current_header = line[1:].strip()
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_header is not None:
            seq_lens[current_header] = current_len
    return seq_lens


def parse_sample_headers(fasta_path, sample_headers):
    """Parse only the sampled headers from FASTA to get lengths."""
    seq_lens = {}
    sample_set = set(sample_headers)
    with open(fasta_path, 'r') as f:
        current_header = None
        current_len = 0
        for line in f:
            if line.startswith('>'):
                if current_header is not None and current_header in sample_set:
                    seq_lens[current_header] = current_len
                current_header = line[1:].strip()
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_header is not None and current_header in sample_set:
            seq_lens[current_header] = current_len
    return seq_lens


def parse_edges(edges_path, sample_headers=None):
    """Parse MASH edge list output.
    
    Format (tab-separated): seq1  seq2  distance  p-value  shared_hashes
    Returns list of (seq1, seq2, dist) tuples.
    """
    edges = []
    seq_names = set()
    with open(edges_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                s1, s2, dist = parts[0], parts[1], float(parts[2])
                seq_names.add(s1)
                seq_names.add(s2)
                edges.append((s1, s2, dist))
    return edges, sorted(seq_names)


def extract_genome_from_header(header):
    """Extract genome accession from prophage header (GCF_XXXXX.1_prophage_N)."""
    parts = header.split('_prophage_')
    if len(parts) >= 2:
        return parts[0]
    # Fallback: try to extract GCF_... pattern
    import re
    m = re.match(r'(GCF_\d+\.\d+)', header)
    if m:
        return m.group(1)
    return header


def extract_prophage_num(header):
    """Extract prophage number from header."""
    parts = header.split('_prophage_')
    if len(parts) >= 2:
        try:
            return int(parts[1])
        except ValueError:
            pass
    return 0


# ============================================================
# Analysis Functions
# ============================================================

def analyze_sequence_lengths(seq_lens):
    """Analyze sequence length distribution."""
    lengths = np.array(list(seq_lens.values()))
    df = pd.DataFrame({
        'sequence': list(seq_lens.keys()),
        'length': list(seq_lens.values()),
    })
    
    stats = {
        'count': len(lengths),
        'mean_length': float(np.mean(lengths)),
        'median_length': float(np.median(lengths)),
        'std_length': float(np.std(lengths)),
        'min_length': float(np.min(lengths)),
        'max_length': float(np.max(lengths)),
        'q25': float(np.percentile(lengths, 25)),
        'q75': float(np.percentile(lengths, 75)),
    }
    
    return df, stats


def build_distance_matrix(edges, seq_names):
    """Build a distance matrix from edges."""
    n = len(seq_names)
    name_to_idx = {name: i for i, name in enumerate(seq_names)}
    
    # Initialize with 1.0 (max distance)
    dist_matrix = np.ones((n, n))
    np.fill_diagonal(dist_matrix, 0.0)
    
    for s1, s2, dist in edges:
        i, j = name_to_idx[s1], name_to_idx[s2]
        dist_matrix[i, j] = dist
        dist_matrix[j, i] = dist
    
    return dist_matrix, name_to_idx


def analyze_distance_distribution(edges):
    """Analyze the distribution of pairwise distances."""
    distances = np.array([e[2] for e in edges])
    
    stats = {
        'num_pairs': len(distances),
        'mean_distance': float(np.mean(distances)),
        'median_distance': float(np.median(distances)),
        'std_distance': float(np.std(distances)),
        'min_distance': float(np.min(distances)),
        'max_distance': float(np.max(distances)),
        'q25': float(np.percentile(distances, 25)),
        'q75': float(np.percentile(distances, 75)),
    }
    
    # Distance bins
    bins = np.arange(0, 0.55, 0.02)
    hist, bin_edges = np.histogram(distances, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    return distances, stats, (hist, bin_edges, bin_centers)


def within_between_genome_distances(edges, seq_names):
    """Compute within-genome and between-genome distance distributions."""
    # Map sequences to genomes
    seq_to_genome = {name: extract_genome_from_header(name) for name in seq_names}
    
    within_distances = []
    between_distances = []
    
    for s1, s2, dist in edges:
        g1, g2 = seq_to_genome[s1], seq_to_genome[s2]
        if g1 == g2:
            within_distances.append(dist)
        else:
            between_distances.append(dist)
    
    within_stats = {}
    between_stats = {}
    
    if within_distances:
        wd = np.array(within_distances)
        within_stats = {
            'count': len(wd),
            'mean': float(np.mean(wd)),
            'median': float(np.median(wd)),
            'std': float(np.std(wd)),
        }
    
    if between_distances:
        bd = np.array(between_distances)
        between_stats = {
            'count': len(bd),
            'mean': float(np.mean(bd)),
            'median': float(np.median(bd)),
            'std': float(np.std(bd)),
        }
    
    return within_distances, between_distances, within_stats, between_stats


def run_mds(dist_matrix, seq_names, n_components=2):
    """Run MDS dimensionality reduction."""
    mds = MDS(n_components=n_components, dissimilarity='precomputed',
              random_state=42, normalized_stress='auto', n_init=1, max_iter=100)
    coords = mds.fit_transform(dist_matrix)
    
    df = pd.DataFrame({
        'sequence': seq_names,
        'MDS1': coords[:, 0],
        'MDS2': coords[:, 1],
        'genome': [extract_genome_from_header(n) for n in seq_names],
    })
    
    stress = mds.stress_
    return df, stress


def run_clustering(dist_matrix, seq_names, method='ward', max_clusters=20):
    """Run hierarchical clustering and assign cluster labels."""
    # Convert distance matrix to condensed form
    condensed = squareform(dist_matrix)
    
    # Linkage
    Z = linkage(condensed, method=method)
    
    # Determine optimal number of clusters using elbow method
    # Look at the last max_clusters merges
    n = len(seq_names)
    heights = Z[-max_clusters:, 2]  # Last max_clusters merge heights
    
    # Find the largest gap in merge heights
    diffs = np.diff(heights)
    if len(diffs) > 0:
        n_clusters = np.argmax(diffs) + 2
    else:
        n_clusters = min(5, n)
    
    n_clusters = max(2, min(n_clusters, max_clusters))
    
    # Assign clusters
    labels = fcluster(Z, n_clusters, criterion='maxclust')
    
    df = pd.DataFrame({
        'sequence': seq_names,
        'cluster': labels,
        'genome': [extract_genome_from_header(n) for n in seq_names],
    })
    
    cluster_sizes = df['cluster'].value_counts().sort_index()
    
    return df, Z, n_clusters, cluster_sizes


def analyze_genome_prophage_diversity(edges, seq_names):
    """Analyze how many prophages per genome and their diversity."""
    seq_to_genome = {name: extract_genome_from_header(name) for name in seq_names}
    
    # Count prophages per genome
    genome_counts = Counter(seq_to_genome.values())
    
    # For genomes with multiple prophages, compute intra-genome diversity
    genome_diversity = {}
    genome_prophage_distances = defaultdict(list)
    
    for s1, s2, dist in edges:
        g1, g2 = seq_to_genome[s1], seq_to_genome[s2]
        if g1 == g2:
            genome_prophage_distances[g1].append(dist)
    
    for genome, dists in genome_prophage_distances.items():
        if len(dists) > 0:
            genome_diversity[genome] = {
                'mean_distance': float(np.mean(dists)),
                'max_distance': float(np.max(dists)),
                'num_prophages': genome_counts[genome],
                'num_pairs': len(dists),
            }
    
    return genome_counts, genome_diversity


# ============================================================
# Plotting Functions
# ============================================================

def plot_length_distribution(seq_df, stats, output_path):
    """Plot sequence length distribution histogram."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    ax = axes[0]
    lengths = seq_df['length']
    ax.hist(lengths / 1000, bins=100, color='steelblue', edgecolor='white', alpha=0.8)
    ax.axvline(stats['mean_length'] / 1000, color='red', linestyle='--', 
               label=f"Mean: {stats['mean_length']/1000:.1f} kb")
    ax.axvline(stats['median_length'] / 1000, color='orange', linestyle='--',
               label=f"Median: {stats['median_length']/1000:.1f} kb")
    ax.set_xlabel('Sequence Length (kb)')
    ax.set_ylabel('Count')
    ax.set_title(f'Prophage Sequence Length Distribution (n={stats["count"]:,})')
    ax.legend()
    
    # Log-scale histogram
    ax = axes[1]
    ax.hist(lengths / 1000, bins=100, color='steelblue', edgecolor='white', alpha=0.8)
    ax.set_yscale('log')
    ax.axvline(stats['mean_length'] / 1000, color='red', linestyle='--',
               label=f"Mean: {stats['mean_length']/1000:.1f} kb")
    ax.axvline(stats['median_length'] / 1000, color='orange', linestyle='--',
               label=f"Median: {stats['median_length']/1000:.1f} kb")
    ax.set_xlabel('Sequence Length (kb)')
    ax.set_ylabel('Count (log scale)')
    ax.set_title('Length Distribution (log scale)')
    ax.legend()
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_distance_distribution(distances, stats, output_path, within_dists=None, between_dists=None):
    """Plot distance distribution histogram."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Overall distance distribution
    ax = axes[0]
    ax.hist(distances, bins=50, color='steelblue', edgecolor='white', alpha=0.8)
    ax.axvline(stats['mean_distance'], color='red', linestyle='--',
               label=f"Mean: {stats['mean_distance']:.4f}")
    ax.axvline(stats['median_distance'], color='orange', linestyle='--',
               label=f"Median: {stats['median_distance']:.4f}")
    ax.set_xlabel('MASH Distance')
    ax.set_ylabel('Count')
    ax.set_title(f'Pairwise Distance Distribution\n(n={stats["num_pairs"]:,} pairs)')
    ax.legend()
    
    # Distance distribution (log scale)
    ax = axes[1]
    ax.hist(distances, bins=50, color='steelblue', edgecolor='white', alpha=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('MASH Distance')
    ax.set_ylabel('Count (log scale)')
    ax.set_title('Distance Distribution (log scale)')
    
    # Within vs between genome distances
    ax = axes[2]
    if within_dists is not None and between_dists is not None:
        ax.hist(within_dists, bins=30, alpha=0.6, color='green', 
                label=f'Within genome (n={len(within_dists):,})', density=True)
        ax.hist(between_dists, bins=30, alpha=0.6, color='gray',
                label=f'Between genome (n={len(between_dists):,})', density=True)
        ax.set_xlabel('MASH Distance')
        ax.set_ylabel('Density')
        ax.set_title('Within vs Between Genome Distances')
        ax.legend()
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_mds_scatter(mds_df, output_path, cluster_df=None):
    """Plot MDS scatter plot with optional cluster coloring."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # MDS colored by genome density
    ax = axes[0]
    scatter = ax.scatter(mds_df['MDS1'], mds_df['MDS2'], 
                         c='steelblue', alpha=0.5, s=8, edgecolors='none')
    ax.set_xlabel('MDS Dimension 1')
    ax.set_ylabel('MDS Dimension 2')
    ax.set_title(f'MDS Plot of Prophage Sequences\n(n={len(mds_df):,})')
    
    # MDS colored by cluster
    ax = axes[1]
    if cluster_df is not None:
        # Merge cluster assignments
        merged = mds_df.merge(cluster_df[['sequence', 'cluster']], on='sequence', how='left')
        colors = merged['cluster']
        n_clusters = colors.nunique()
        scatter = ax.scatter(merged['MDS1'], merged['MDS2'], 
                            c=colors, cmap='tab20', alpha=0.6, s=8, edgecolors='none')
        cbar = plt.colorbar(scatter, ax=ax, label='Cluster')
        ax.set_title(f'MDS Colored by Cluster (n={n_clusters} clusters)')
    else:
        ax.scatter(mds_df['MDS1'], mds_df['MDS2'],
                   c='steelblue', alpha=0.5, s=8, edgecolors='none')
        ax.set_title('MDS Plot')
    ax.set_xlabel('MDS Dimension 1')
    ax.set_ylabel('MDS Dimension 2')
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_cluster_analysis(seq_df, cluster_sizes, output_path):
    """Plot cluster size distribution and characteristics."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Cluster size bar plot
    ax = axes[0]
    sizes = cluster_sizes.values
    labels = [str(c) for c in cluster_sizes.index]
    colors = plt.cm.tab20(np.linspace(0, 1, len(sizes)))
    ax.bar(range(len(sizes)), sizes, color=colors, edgecolor='white')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Number of Sequences')
    ax.set_title('Cluster Size Distribution')
    ax.set_xticks(range(len(sizes)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    
    # Cluster size distribution (log)
    ax = axes[1]
    ax.hist(sizes, bins=min(50, len(sizes)), color='steelblue', edgecolor='white')
    ax.set_yscale('log')
    ax.set_xlabel('Cluster Size')
    ax.set_ylabel('Number of Clusters (log scale)')
    ax.set_title('Cluster Size Distribution (log scale)')
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_dendrogram(Z, output_path, max_d=200):
    """Plot dendrogram (truncated for large datasets)."""
    fig, ax = plt.subplots(figsize=(16, 6))
    
    dendrogram(Z, truncate_mode='lastp', p=min(30, len(Z)),
               leaf_rotation=90., leaf_font_size=8.,
               show_contracted=True, color_threshold=0.7*max(Z[:, 2]),
               ax=ax)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Distance')
    ax.set_title('Hierarchical Clustering Dendrogram (truncated)')
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_genome_prophage_count(genome_counts, output_path):
    """Plot distribution of prophage count per genome."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    counts = np.array(list(genome_counts.values()))
    
    # Histogram
    ax = axes[0]
    max_count = min(50, int(np.percentile(counts, 99)))
    ax.hist(counts, bins=range(1, max_count + 2), color='steelblue', edgecolor='white', alpha=0.8)
    ax.set_xlabel('Number of Prophages per Genome')
    ax.set_ylabel('Number of Genomes')
    ax.set_title(f'Prophage Count per Genome\n(n={len(genome_counts):,} genomes)')
    
    # Log scale
    ax = axes[1]
    ax.hist(counts, bins=range(1, max_count + 2), color='steelblue', edgecolor='white', alpha=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('Number of Prophages per Genome')
    ax.set_ylabel('Number of Genomes (log scale)')
    ax.set_title('Prophage Count per Genome (log scale)')
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_genome_diversity(genome_diversity, output_path):
    """Plot intra-genome prophage diversity."""
    if not genome_diversity:
        return output_path
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    genomes = list(genome_diversity.keys())
    mean_dists = [genome_diversity[g]['mean_distance'] for g in genomes]
    counts = [genome_diversity[g]['num_prophages'] for g in genomes]
    
    scatter = ax.scatter(counts, mean_dists, alpha=0.6, s=20, 
                         c=mean_dists, cmap='viridis', edgecolors='none')
    ax.set_xlabel('Number of Prophages in Genome')
    ax.set_ylabel('Mean Intra-Genome Distance')
    ax.set_title('Intra-Genome Prophage Diversity')
    plt.colorbar(scatter, ax=ax, label='Mean Distance')
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return output_path


# ============================================================
# Typst Report Generation
# ============================================================

def generate_typst_report(length_stats, dist_stats, within_stats, between_stats,
                          mds_df, cluster_df, cluster_sizes, n_clusters, mds_stress,
                          genome_counts, genome_diversity, seq_len, output_path):
    """Generate a Typst presentation with all figures."""
    
    total_genomes = len(genome_counts)
    genomes_with_multiple = sum(1 for c in genome_counts.values() if c > 1)
    
    def esc(s):
        return str(s).replace('<', '\\<').replace('>', '\\>')
    
    typst = []
    typst.append('#set page(width: 16in, height: 9in, margin: (x: 0.5in, y: 0.3in))')
    typst.append('#set text(size: 10pt)')
    typst.append('')
    typst.append('#align(center)[#text(size: 24pt, weight: "bold")[Full Prophage Homology Survey]]')
    typst.append('#align(center)[#text(size: 14pt)[MASH-based Analysis of 132,393 Prophage Sequences from 26,000 E. coli Genomes]]')
    typst.append('')
    typst.append('#grid(columns: (1fr, 1fr, 1fr, 1fr),')
    typst.append('  align: center,')
    typst.append(f'  [#text(size: 18pt, weight: "bold")[{seq_len:,}]\\ #text(size: 9pt)[Total Sequences]],')
    typst.append(f'  [#text(size: 18pt, weight: "bold")[{total_genomes:,}]\\ #text(size: 9pt)[Genomes]],')
    typst.append(f'  [#text(size: 18pt, weight: "bold")[{length_stats["mean_length"]/1000:.1f} kb]\\ #text(size: 9pt)[Mean Length]],')
    typst.append(f'  [#text(size: 18pt, weight: "bold")[{dist_stats["num_pairs"]:,}]\\ #text(size: 9pt)[{esc("Similar Pairs (d<0.5)")}]])')
    typst.append('')
    typst.append('#pagebreak()')
    
    # Section 1: Sequence Length Distribution
    typst.append('#text(size: 18pt, weight: "bold")[1. Sequence Length Distribution]')
    typst.append('#v(0.2cm)')
    typst.append('#grid(columns: (1.5fr, 1fr), gutter: 0.3cm,')
    typst.append('  [#image("length_distribution.png", width: 100%)],')
    typst.append('  [#table(columns: (1fr, 1fr),')
    typst.append('    [Metric], [Value],')
    typst.append(f'    [Count], [{length_stats["count"]:,}],')
    typst.append(f'    [Mean], [{length_stats["mean_length"]/1000:.1f} kb],')
    typst.append(f'    [Median], [{length_stats["median_length"]/1000:.1f} kb],')
    typst.append(f'    [Std Dev], [{length_stats["std_length"]/1000:.1f} kb],')
    typst.append(f'    [Min], [{length_stats["min_length"]/1000:.1f} kb],')
    typst.append(f'    [Max], [{length_stats["max_length"]/1000:.1f} kb],')
    typst.append(f'    [Q25], [{length_stats["q25"]/1000:.1f} kb],')
    typst.append(f'    [Q75], [{length_stats["q75"]/1000:.1f} kb])')
    typst.append('  ])')
    typst.append('')
    typst.append('#pagebreak()')
    
    # Section 2: Distance Distribution
    typst.append('#text(size: 18pt, weight: "bold")[2. MASH Distance Distribution]')
    typst.append('#v(0.2cm)')
    typst.append('#image("distance_distribution.png", width: 100%)')
    typst.append('#v(0.2cm)')
    typst.append('#grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.2cm,')
    typst.append('  [#table(columns: (1fr, 1fr), [Metric], [Overall],')
    typst.append(f'    [Pairs], [{dist_stats["num_pairs"]:,}],')
    typst.append(f'    [Mean], [{dist_stats["mean_distance"]:.4f}],')
    typst.append(f'    [Median], [{dist_stats["median_distance"]:.4f}],')
    typst.append(f'    [Std], [{dist_stats["std_distance"]:.4f}])],')
    if within_stats:
        typst.append('  [#table(columns: (1fr, 1fr), [Metric], [Within Genome],')
        typst.append(f'    [Pairs], [{within_stats["count"]:,}],')
        typst.append(f'    [Mean], [{within_stats["mean"]:.4f}],')
        typst.append(f'    [Median], [{within_stats["median"]:.4f}],')
        typst.append(f'    [Std], [{within_stats["std"]:.4f}])],')
    if between_stats:
        typst.append('  [#table(columns: (1fr, 1fr), [Metric], [Between Genome],')
        typst.append(f'    [Pairs], [{between_stats["count"]:,}],')
        typst.append(f'    [Mean], [{between_stats["mean"]:.4f}],')
        typst.append(f'    [Median], [{between_stats["median"]:.4f}],')
        typst.append(f'    [Std], [{between_stats["std"]:.4f}])],')
    typst.append(')')
    typst.append('')
    typst.append('#pagebreak()')
    
    # Section 3: MDS Analysis
    typst.append('#text(size: 18pt, weight: "bold")[3. Multidimensional Scaling (MDS)]')
    typst.append('#v(0.2cm)')
    typst.append(f'#text(size: 11pt)[MDS stress: {mds_stress:.4f}]')
    typst.append('#v(0.2cm)')
    typst.append('#image("mds_scatter.png", width: 100%)')
    typst.append('')
    typst.append('#pagebreak()')
    
    # Section 4: Cluster Analysis
    typst.append('#text(size: 18pt, weight: "bold")[4. Hierarchical Clustering]')
    typst.append('#v(0.2cm)')
    typst.append(f'#text(size: 11pt)[Optimal number of clusters: {n_clusters}]')
    typst.append('#v(0.2cm)')
    typst.append('#image("cluster_analysis.png", width: 100%)')
    typst.append('#v(0.2cm)')
    typst.append('#image("dendrogram.png", width: 100%)')
    typst.append('')
    typst.append('#pagebreak()')
    
    # Section 5: Genome-level Analysis
    typst.append('#text(size: 18pt, weight: "bold")[5. Genome-Level Prophage Diversity]')
    typst.append('#v(0.2cm)')
    typst.append('#grid(columns: (1.5fr, 1fr), gutter: 0.3cm,')
    typst.append('  [#image("genome_prophage_count.png", width: 100%)],')
    typst.append('  [#table(columns: (1fr, 1fr),')
    typst.append('    [Metric], [Value],')
    typst.append(f'    [Genomes with prophages], [{total_genomes:,}],')
    typst.append(f'    [Genomes with 2+ prophages], [{genomes_with_multiple:,}],')
    typst.append(f'    [Max prophages/genome], [{max(genome_counts.values())}],')
    typst.append(f'    [Mean prophages/genome], [{np.mean(list(genome_counts.values())):.1f}])')
    typst.append('  ])')
    typst.append('')
    typst.append('#v(0.3cm)')
    typst.append('#image("genome_diversity.png", width: 80%)')
    typst.append('')
    typst.append('#pagebreak()')
    
    # Section 6: Summary
    typst.append('#text(size: 18pt, weight: "bold")[6. Summary Statistics]')
    typst.append('#v(0.2cm)')
    typst.append('#table(columns: (3fr, 1fr, 2fr),')
    typst.append('  [Metric], [Value], [Notes],')
    typst.append(f'  [Total sequences], [{seq_len:,}], [Full prophage set],')
    typst.append(f'  [Total genomes], [{total_genomes:,}], [From 26k E. coli genomes],')
    typst.append(f'  [Mean sequence length], [{length_stats["mean_length"]/1000:.1f} kb], [Range: {length_stats["min_length"]/1000:.1f}-{length_stats["max_length"]/1000:.1f} kb],')
    typst.append(f'  [{esc("Similar pairs (d<0.5)")}], [{dist_stats["num_pairs"]:,}], [Out of ~{(seq_len*(seq_len-1)/2):.0e} total pairs],')
    typst.append(f'  [{esc("Mean distance (similar pairs)")}], [{dist_stats["mean_distance"]:.4f}], [{esc("Only pairs with d<0.5")}],')
    typst.append(f'  [MDS stress], [{mds_stress:.4f}], [Lower is better],')
    typst.append(f'  [Number of clusters], [{n_clusters}], [Hierarchical clustering],')
    typst.append(f'  [Genomes with 2+ prophages], [{genomes_with_multiple:,}], [Intra-genome diversity],')
    typst.append(')')
    typst.append('')
    typst.append('#v(1cm)')
    typst.append('#align(center)[#text(size: 10pt)[Generated by MASH v2.3 + Python analysis pipeline. \\')
    typst.append(f'Analysis date: 2026-07-29. Sample size: 10,000 sequences for pairwise analysis.]]')
    
    with open(output_path, 'w') as f:
        f.write('\n'.join(typst))
    
    return output_path


# ============================================================
# Main Pipeline
# ============================================================

def main():
    print("=" * 60)
    print("Full Prophage Homology Survey - Analysis Pipeline")
    print("=" * 60)
    
    # Step 1: Parse sequence lengths from full FASTA
    print("\n[1/6] Parsing sequence lengths from full FASTA...")
    seq_lens = parse_sequence_lengths(FULL_FASTA)
    print(f"  Parsed {len(seq_lens)} sequences")
    
    # Step 2: Analyze sequence lengths
    print("\n[2/6] Analyzing sequence length distribution...")
    seq_df, length_stats = analyze_sequence_lengths(seq_lens)
    print(f"  Mean: {length_stats['mean_length']/1000:.1f} kb, "
          f"Median: {length_stats['median_length']/1000:.1f} kb, "
          f"Range: {length_stats['min_length']/1000:.1f}-{length_stats['max_length']/1000:.1f} kb")
    
    # Step 3: Parse MASH edge list
    print("\n[3/6] Parsing MASH edge list...")
    if not SAMPLE_EDGES.exists():
        print(f"  ERROR: Edge file not found: {SAMPLE_EDGES}")
        print("  Make sure MASH triangle has completed.")
        sys.exit(1)
    
    edges, seq_names = parse_edges(SAMPLE_EDGES)
    print(f"  Parsed {len(edges):,} edges, {len(seq_names):,} sequences")
    
    # Step 3a: Distance distribution
    print("\n  Computing distance distribution...")
    distances, dist_stats, dist_hist = analyze_distance_distribution(edges)
    print(f"  Mean distance: {dist_stats['mean_distance']:.4f}, "
          f"Median: {dist_stats['median_distance']:.4f}")
    
    # Step 3b: Within/between genome distances
    print("\n  Computing within/between genome distances...")
    within_dists, between_dists, within_stats, between_stats = \
        within_between_genome_distances(edges, seq_names)
    if within_stats:
        print(f"  Within-genome pairs: {within_stats['count']:,}, "
              f"mean: {within_stats['mean']:.4f}")
    if between_stats:
        print(f"  Between-genome pairs: {between_stats['count']:,}, "
              f"mean: {between_stats['mean']:.4f}")
    
    # Step 4: Subsample for MDS and clustering (faster)
    print("\n[4/6] Subsampling for MDS and clustering...")
    
    # Use up to 5000 sequences for pairwise analysis
    max_mds = 5000
    if len(seq_names) > max_mds:
        import random
        random.seed(42)
        mds_seqs = set(random.sample(seq_names, max_mds))
        # Filter edges to only include MDS sequences
        mds_edges = [(s1, s2, d) for s1, s2, d in edges if s1 in mds_seqs and s2 in mds_seqs]
        mds_names = sorted(mds_seqs)
        print(f"  Subsampled to {len(mds_names)} sequences for MDS/clustering")
    else:
        mds_edges = edges
        mds_names = seq_names
    
    print("\n  Building distance matrix...")
    dist_matrix, name_to_idx = build_distance_matrix(mds_edges, mds_names)
    print(f"  Distance matrix: {dist_matrix.shape[0]}x{dist_matrix.shape[1]}")
    
    # MDS (faster: n_init=1, max_iter=100)
    print("  Running MDS (this may take a few minutes)...")
    mds_df, mds_stress = run_mds(dist_matrix, mds_names)
    print(f"  MDS stress: {mds_stress:.4f}")
    
    # Hierarchical clustering
    print("\n  Running hierarchical clustering...")
    cluster_df, Z, n_clusters, cluster_sizes = run_clustering(dist_matrix, mds_names)
    print(f"  Optimal clusters: {n_clusters}")
    print(f"  Cluster sizes: min={cluster_sizes.min()}, max={cluster_sizes.max()}")
    
    # Step 5: Genome-level analysis
    print("\n[5/6] Analyzing genome-level prophage diversity...")
    genome_counts, genome_diversity = analyze_genome_prophage_diversity(edges, seq_names)
    print(f"  Genomes: {len(genome_counts)}, "
          f"Multi-prophage: {sum(1 for c in genome_counts.values() if c > 1)}")
    
    # Step 6: Generate plots and report
    print("\n[6/6] Generating plots and PDF report...")
    
    # Generate all PNG figures
    print("  Generating length distribution plot...")
    plot_length_distribution(seq_df, length_stats, 
                             OUTPUT_DIR / 'length_distribution.png')
    
    print("  Generating distance distribution plot...")
    plot_distance_distribution(distances, dist_stats,
                               OUTPUT_DIR / 'distance_distribution.png',
                               within_dists, between_dists)
    
    print("  Generating MDS scatter plot...")
    plot_mds_scatter(mds_df, OUTPUT_DIR / 'mds_scatter.png', cluster_df)
    
    print("  Generating cluster analysis plot...")
    plot_cluster_analysis(seq_df, cluster_sizes, 
                          OUTPUT_DIR / 'cluster_analysis.png')
    
    print("  Generating dendrogram...")
    plot_dendrogram(Z, OUTPUT_DIR / 'dendrogram.png')
    
    print("  Generating genome prophage count plot...")
    plot_genome_prophage_count(genome_counts, 
                               OUTPUT_DIR / 'genome_prophage_count.png')
    
    print("  Generating genome diversity plot...")
    plot_genome_diversity(genome_diversity, 
                          OUTPUT_DIR / 'genome_diversity.png')
    
    # Generate Typst report
    print("  Generating Typst report...")
    generate_typst_report(
        length_stats, dist_stats, within_stats, between_stats,
        mds_df, cluster_df, cluster_sizes, n_clusters, mds_stress,
        genome_counts, genome_diversity, len(seq_lens),
        REPORT_TYP
    )
    
    # Compile Typst to PDF
    print("  Compiling Typst to PDF...")
    result = subprocess.run(
        ['typst', 'compile', str(REPORT_TYP), str(REPORT_PDF)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  WARNING: Typst compilation failed: {result.stderr}")
        print("  PDF may not be generated. Check Typst installation.")
    else:
        print(f"  PDF generated: {REPORT_PDF}")
    
    # Save output CSV files
    print("\n  Saving output CSV files...")
    mds_df.to_csv(MDS_COORDS, index=False)
    print(f"  MDS coordinates: {MDS_COORDS}")
    
    cluster_df.to_csv(CLUSTERS, index=False)
    print(f"  Cluster assignments: {CLUSTERS}")
    
    # Distance stats
    dist_stats_df = pd.DataFrame([{
        'metric': 'overall',
        **{k: v for k, v in dist_stats.items()}
    }])
    if within_stats:
        dist_stats_df = pd.concat([dist_stats_df, pd.DataFrame([{
            'metric': 'within_genome',
            **{k: v for k, v in within_stats.items()}
        }])])
    if between_stats:
        dist_stats_df = pd.concat([dist_stats_df, pd.DataFrame([{
            'metric': 'between_genome',
            **{k: v for k, v in between_stats.items()}
        }])])
    # Add length stats
    dist_stats_df = pd.concat([dist_stats_df, pd.DataFrame([{
        'metric': 'sequence_length',
        'num_pairs': length_stats['count'],
        'mean': length_stats['mean_length'],
        'median': length_stats['median_length'],
        'std': length_stats['std_length'],
        'min': length_stats['min_length'],
        'max': length_stats['max_length'],
        'q25': length_stats['q25'],
        'q75': length_stats['q75'],
    }])])
    dist_stats_df.to_csv(DIST_STATS, index=False)
    print(f"  Distance stats: {DIST_STATS}")
    
    # Summary
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print(f"  Output files:")
    print(f"    PDF Report: {REPORT_PDF}")
    print(f"    MDS Coords: {MDS_COORDS}")
    print(f"    Clusters:   {CLUSTERS}")
    print(f"    Stats:      {DIST_STATS}")
    print(f"    Typst:      {REPORT_TYP}")
    print("=" * 60)


if __name__ == '__main__':
    main()