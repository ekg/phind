#!/usr/bin/env python3
"""
Full 132k×132k clustered heatmap using community detection + binned rendering.

Pipeline:
1. Read the full edge list (36M edges, d<0.5) from main repo
2. Parse all 132k sequence names from FASTA
3. Build networkx graph from edges and run Louvain community detection
4. Assign isolated sequences (no edges at d<0.5) as singleton communities
5. Sort all 132k sequences by community, then by intra-cluster degree
6. Create 2D binned heatmap (2000×2000 bins) aggregating all 36M edges
7. Create cluster-level heatmap (mean distance between communities)
8. Generate PNG, standalone PDF
9. Update Typst report, compile, and prepare for SCP
"""

import os
import sys
import gc
import time
import warnings
import subprocess
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcolors

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
MAIN_REPO = Path('/home/erikg/phind')
DATA_DIR = Path('prophage_homology_survey')
EDGES_FILE = Path('/home/erikg/phind/prophage_homology_survey/full_prophages_edges.tsv')
FASTA_FILE = Path('/home/erikg/phind/prophage_homology_survey/full_prophages.fa')
OUTPUT_DIR = DATA_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
HEATMAP_PNG = OUTPUT_DIR / 'full_heatmap.png'
HEATMAP_STANDALONE_PDF = OUTPUT_DIR / 'full_heatmap_standalone.pdf'
CLUSTER_HEATMAP_PNG = OUTPUT_DIR / 'cluster_heatmap.png'
CLUSTER_MEMBERSHIP_CSV = OUTPUT_DIR / 'full_heatmap_clusters.csv'
REPORT_TYP = OUTPUT_DIR / 'full_prophage_homology_report.typ'

# Binning parameters
N_BINS = 2000

# Plotting
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
})


# ============================================================
# Step 0: Parse all sequence names from FASTA
# ============================================================

def parse_all_sequences(fasta_path):
    """Parse ALL sequence names from the FASTA file (132k sequences)."""
    print("\n[0/6] Parsing all sequence names from FASTA...")
    t0 = time.time()
    seq_names = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
                seq_names.append(header)
    print(f"  Parsed {len(seq_names):,} sequence names in {time.time()-t0:.1f}s")
    return seq_names


# ============================================================
# Step 1: Build graph from edge list
# ============================================================

def build_graph(edges_path, all_seq_names):
    """
    Build a networkx graph from the edge list.
    Edge list format: seq1, seq2, dist, p-val, shared-hashes (tab-separated)
    Returns: graph, set of connected nodes
    """
    print("\n[1/6] Building graph from edge list...")
    t0 = time.time()
    all_seq_set = set(all_seq_names)
    
    G = nx.Graph()
    
    # Count total lines first
    print("  Counting edges...")
    total_lines = 0
    with open(edges_path, 'r') as f:
        for _ in f:
            total_lines += 1
    print(f"  Total edges: {total_lines:,}")
    
    # Stream edges into graph
    print("  Adding edges to graph (this may take a few minutes)...")
    batch_size = 500000
    edges_batch = []
    connected_set = set()
    edge_count = 0
    max_dist = 0.0
    
    with open(edges_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                s1, s2, dist = parts[0], parts[1], float(parts[2])
                # Only include sequences that are in our FASTA
                if s1 not in all_seq_set or s2 not in all_seq_set:
                    continue
                weight = 1.0 - dist  # Louvain uses similarity weight
                edges_batch.append((s1, s2, weight))
                connected_set.add(s1)
                connected_set.add(s2)
                edge_count += 1
                if dist > max_dist:
                    max_dist = dist
            
            if len(edges_batch) >= batch_size:
                G.add_weighted_edges_from(edges_batch)
                edges_batch = []
                if edge_count % 2000000 == 0:
                    elapsed = time.time() - t0
                    print(f"    {edge_count:,} edges processed ({elapsed:.1f}s), "
                          f"{len(connected_set):,} connected seqs so far")
    
    # Add remaining batch
    if edges_batch:
        G.add_weighted_edges_from(edges_batch)
    
    print(f"  Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    print(f"  Connected sequences (with d<0.5 edges): {len(connected_set):,}")
    print(f"  Isolated sequences (no edges): {len(all_seq_names) - len(connected_set):,}")
    print(f"  Max distance in edges: {max_dist:.4f}")
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return G, connected_set


# ============================================================
# Step 2: Community detection (Louvain)
# ============================================================

def run_community_detection(G, all_seq_names, connected_set):
    """
    Run Louvain community detection on the graph.
    Assigns isolated sequences (no edges) as singleton communities.
    Returns dict mapping node -> community_id, and community size info.
    """
    print("\n[2/6] Running community detection...")
    t0 = time.time()
    
    from networkx.algorithms.community import louvain_communities
    
    # Run Louvain on the connected component
    print("  Running Louvain on connected sequences...")
    communities_iter = louvain_communities(G, weight='weight', seed=42)
    
    # Assign community IDs to connected sequences
    node_to_community = {}
    connected_comm_count = 0
    connected_comm_sizes = []
    for i, comm in enumerate(communities_iter):
        connected_comm_sizes.append(len(comm))
        for node in comm:
            node_to_community[node] = i
        connected_comm_count = i + 1
    
    print(f"  Louvain found {connected_comm_count} communities from {len(connected_set):,} connected sequences")
    print(f"  Connected community sizes: min={min(connected_comm_sizes)}, "
          f"max={max(connected_comm_sizes)}, mean={np.mean(connected_comm_sizes):.0f}")
    
    # Assign isolated sequences as singleton communities
    iso_count = 0
    for seq in all_seq_names:
        if seq not in node_to_community:
            node_to_community[seq] = connected_comm_count + iso_count
            iso_count += 1
    
    n_communities = connected_comm_count + iso_count
    print(f"  Total communities: {n_communities:,} "
          f"({connected_comm_count} connected + {iso_count:,} singletons)")
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return node_to_community, connected_comm_count, iso_count


# ============================================================
# Step 3: Order sequences by community + degree
# ============================================================

def order_sequences(node_to_community, G, all_seq_names):
    """
    Order sequences by community, then by degree within community.
    Isolated sequences (not in G) have degree 0.
    Returns dict mapping node -> position (0..n-1), and list of (node, community, degree).
    """
    print("\n[3/6] Ordering sequences by community + degree...")
    t0 = time.time()
    
    # Group nodes by community
    community_nodes = defaultdict(list)
    for node, comm in node_to_community.items():
        community_nodes[comm].append(node)
    
    # Sort communities by size (descending) for visual clarity
    sorted_communities = sorted(community_nodes.keys(), 
                                key=lambda c: len(community_nodes[c]), 
                                reverse=True)
    
    # Build ordered list: within each community, sort by degree (descending)
    ordered_nodes = []
    node_to_pos = {}
    pos = 0
    
    for comm in sorted_communities:
        nodes = community_nodes[comm]
        # Sort by degree within community (descending)
        nodes_sorted = sorted(nodes, key=lambda n: G.degree(n) if G.has_node(n) else 0, 
                             reverse=True)
        for n in nodes_sorted:
            node_to_pos[n] = pos
            degree = G.degree(n) if G.has_node(n) else 0
            ordered_nodes.append((n, comm, degree))
            pos += 1
    
    print(f"  Ordered {len(ordered_nodes):,} sequences in {len(sorted_communities):,} communities")
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return node_to_pos, ordered_nodes


# ============================================================
# Step 4: Create binned heatmap
# ============================================================

def create_binned_heatmap(edges_path, node_to_pos, n_bins=2000):
    """
    Create a 2D binned heatmap by iterating through all edges.
    Each cell shows the mean distance for all pairs in that bin region.
    Bins with no edges get distance = 1.0 (max distance).
    """
    print(f"\n[4/6] Creating binned heatmap ({n_bins}×{n_bins})...")
    t0 = time.time()
    n_nodes = len(node_to_pos)
    
    # Initialize bin accumulators
    bin_sum = np.zeros((n_bins, n_bins), dtype=np.float64)
    bin_count = np.zeros((n_bins, n_bins), dtype=np.int64)
    
    def pos_to_bin(p):
        """Convert position (0..n_nodes-1) to bin index (0..n_bins-1)."""
        return min(int(p * n_bins / n_nodes), n_bins - 1)
    
    print(f"  Processing {n_nodes:,} sequences across {n_bins}×{n_bins} bins...")
    
    # Stream edges and accumulate
    edge_count = 0
    valid_edge_count = 0
    report_interval = 5000000
    
    with open(edges_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                s1, s2, dist = parts[0], parts[1], float(parts[2])
                edge_count += 1
                
                # Get positions
                pos1 = node_to_pos.get(s1)
                pos2 = node_to_pos.get(s2)
                if pos1 is None or pos2 is None:
                    continue
                
                valid_edge_count += 1
                
                # Map to bins
                i = pos_to_bin(pos1)
                j = pos_to_bin(pos2)
                
                # Accumulate
                bin_sum[i, j] += dist
                bin_sum[j, i] += dist
                bin_count[i, j] += 1
                bin_count[j, i] += 1
                
                if edge_count % report_interval == 0:
                    elapsed = time.time() - t0
                    print(f"    {edge_count:,} edges processed ({elapsed:.1f}s)")
    
    print(f"  Total edges: {edge_count:,}, valid: {valid_edge_count:,}")
    
    # Compute mean distance for each bin
    with np.errstate(divide='ignore', invalid='ignore'):
        bin_mean = np.where(bin_count > 0, bin_sum / bin_count, 1.0)
    
    # Self-pairs in the same bin: diagonal should be 0
    for i in range(n_bins):
        if bin_count[i, i] > 0:
            # Only internal edges contributed, not self-pairs
            pass
    
    print(f"  Bins with data: {np.sum(bin_count > 0):,} / {n_bins * n_bins:,} "
          f"({100 * np.sum(bin_count > 0) / (n_bins * n_bins):.1f}%)")
    print(f"  Mean distance in filled bins: {np.mean(bin_mean[bin_count > 0]):.4f}")
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return bin_mean, bin_count


def render_heatmap(bin_mean, ordered_nodes, n_communities, 
                   output_png, output_pdf=None, title=None):
    """Render the binned heatmap as PNG and optional PDF."""
    print(f"  Rendering heatmap...")
    t0 = time.time()
    
    n_bins = bin_mean.shape[0]
    n_nodes = len(ordered_nodes)
    
    # Get community boundaries for overlay
    community_boundaries = []
    current_comm = ordered_nodes[0][1] if ordered_nodes else None
    start_pos = 0
    for i, (node, comm, degree) in enumerate(ordered_nodes):
        if comm != current_comm:
            # Convert position to bin coordinate
            bin_pos = i * n_bins / n_nodes
            community_boundaries.append((bin_pos, current_comm, i - start_pos))
            current_comm = comm
            start_pos = i
    # Last community
    bin_pos = n_nodes * n_bins / n_nodes
    community_boundaries.append((bin_pos, current_comm, n_nodes - start_pos))
    
    # Create figure for PNG
    fig, ax = plt.subplots(figsize=(10, 9))
    
    vmin, vmax = 0.0, 0.5
    im = ax.imshow(bin_mean, aspect='equal', cmap='viridis_r', 
                   vmin=vmin, vmax=vmax, interpolation='nearest',
                   extent=[0, n_bins, 0, n_bins])
    
    # Add community boundary lines (only draw for communities with >1 member)
    for bin_pos, comm, size in community_boundaries:
        if size > 1:
            pos = bin_pos
            ax.axhline(pos, color='white', linewidth=0.3, alpha=0.3)
            ax.axvline(pos, color='white', linewidth=0.3, alpha=0.3)
    
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('MASH Distance', fontsize=10)
    
    ax.set_xlabel('Sequence (ordered by community)', fontsize=10)
    ax.set_ylabel('Sequence (ordered by community)', fontsize=10)
    
    n_connected = sum(1 for _, c, _ in ordered_nodes if c < n_communities - 1)
    
    if title:
        ax.set_title(title, fontsize=11, fontweight='bold')
    else:
        ax.set_title(f'Prophage Homology Heatmap\n'
                     f'{n_nodes:,} sequences, {n_communities:,} communities, '
                     f'{n_bins}×{n_bins} bins', 
                     fontsize=11, fontweight='bold')
    
    ax.text(0.02, 0.98, f'{n_nodes:,} sequences × {n_bins} bins',
            transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches='tight', pad_inches=0.1)
    print(f"  Saved: {output_png}")
    plt.close(fig)
    
    if output_pdf:
        fig_pdf, ax_pdf = plt.subplots(figsize=(8, 8))
        im_pdf = ax_pdf.imshow(bin_mean, aspect='equal', cmap='viridis_r',
                               vmin=vmin, vmax=vmax, interpolation='nearest',
                               extent=[0, n_bins, 0, n_bins])
        
        for bin_pos, comm, size in community_boundaries:
            if size > 1:
                pos = bin_pos
                ax_pdf.axhline(pos, color='white', linewidth=0.3, alpha=0.3)
                ax_pdf.axvline(pos, color='white', linewidth=0.3, alpha=0.3)
        
        cbar_pdf = fig_pdf.colorbar(im_pdf, ax=ax_pdf, shrink=0.8, pad=0.02)
        cbar_pdf.set_label('MASH Distance', fontsize=10)
        
        ax_pdf.set_xlabel('Sequence (ordered by community)', fontsize=10)
        ax_pdf.set_ylabel('Sequence (ordered by community)', fontsize=10)
        ax_pdf.set_title(f'Prophage Homology Heatmap (n={n_nodes:,})', 
                         fontsize=11, fontweight='bold')
        
        plt.tight_layout()
        fig_pdf.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.2)
        print(f"  Saved: {output_pdf}")
        plt.close(fig_pdf)
    
    print(f"  Render time: {time.time() - t0:.1f}s")


# ============================================================
# Step 5: Cluster-level heatmap
# ============================================================

def create_cluster_heatmap(edges_path, node_to_community, n_communities, output_png,
                           connected_comm_count):
    """
    Create a cluster-level heatmap showing mean distance between communities.
    Only shows connected communities (those with >1 member), not singletons.
    """
    connected_communities = set(range(connected_comm_count))
    
    print(f"\n[5/6] Creating cluster-level heatmap...")
    t0 = time.time()
    
    # Accumulate sum and count for each pair of connected communities
    cluster_sum = np.zeros((connected_comm_count, connected_comm_count), dtype=np.float64)
    cluster_count = np.zeros((connected_comm_count, connected_comm_count), dtype=np.int64)
    
    edge_count = 0
    report_interval = 5000000
    
    with open(edges_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                s1, s2, dist = parts[0], parts[1], float(parts[2])
                
                c1 = node_to_community.get(s1)
                c2 = node_to_community.get(s2)
                if c1 is None or c2 is None:
                    continue
                
                # Only accumulate for connected communities
                if c1 in connected_communities and c2 in connected_communities:
                    cluster_sum[c1, c2] += dist
                    cluster_sum[c2, c1] += dist
                    cluster_count[c1, c2] += 1
                    cluster_count[c2, c1] += 1
                
                edge_count += 1
                if edge_count % report_interval == 0:
                    elapsed = time.time() - t0
                    print(f"    {edge_count:,} edges processed ({elapsed:.1f}s)")
    
    print(f"  Total edges: {edge_count:,}")
    
    # Compute mean distance for each community pair
    with np.errstate(divide='ignore', invalid='ignore'):
        cluster_mean = np.where(cluster_count > 0, cluster_sum / cluster_count, 1.0)
    
    # Set diagonal to 0 (within-community distance is 0 by definition)
    np.fill_diagonal(cluster_mean, 0.0)
    
    # Compute community sizes
    community_sizes = Counter(node_to_community.values())
    
    print(f"  Connected communities: {connected_comm_count}")
    print(f"  Community pairs with data: {np.sum(cluster_count > 0):,} / {connected_comm_count * connected_comm_count:,}")
    
    # Render
    print(f"  Rendering cluster heatmap...")
    
    # Sort communities by size for visual clarity
    sorted_communities = sorted(range(connected_comm_count), 
                                 key=lambda c: community_sizes.get(c, 0), 
                                 reverse=True)
    
    # Reorder matrix
    cluster_mean_ordered = cluster_mean[sorted_communities][:, sorted_communities]
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    data_max = np.max(cluster_mean_ordered[cluster_mean_ordered < 0.99])
    vmin, vmax = 0.0, min(0.5, data_max)
    
    im = ax.imshow(cluster_mean_ordered, aspect='equal', cmap='viridis_r',
                   vmin=vmin, vmax=vmax, interpolation='nearest')
    
    tick_labels = [f'C{c+1}\n({community_sizes.get(sorted_communities[c], 0)})' 
                   for c in range(connected_comm_count)]
    ax.set_xticks(range(connected_comm_count))
    ax.set_yticks(range(connected_comm_count))
    ax.set_xticklabels(tick_labels, fontsize=6, rotation=90)
    ax.set_yticklabels(tick_labels, fontsize=6)
    
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Mean MASH Distance', fontsize=10)
    
    ax.set_xlabel('Community', fontsize=10)
    ax.set_ylabel('Community', fontsize=10)
    ax.set_title(f'Cross-Community Prophage Homology\n'
                 f'{connected_comm_count} connected communities, mean pairwise distance',
                 fontsize=12, fontweight='bold')
    
    # Add text annotations
    for i in range(connected_comm_count):
        for j in range(connected_comm_count):
            if i != j and cluster_mean_ordered[i, j] < 0.99:
                val = cluster_mean_ordered[i, j]
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', 
                       fontsize=5, color='white' if val > 0.2 else 'black')
    
    plt.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches='tight', pad_inches=0.1)
    print(f"  Saved: {output_png}")
    plt.close(fig)
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return cluster_mean, sorted_communities


# ============================================================
# Typst Report Update
# ============================================================

def update_typst_report(typ_path, n_communities, n_sequences, n_edges, n_bins,
                        connected_comm_count, n_isolated):
    """Add the new heatmap slides to the Typst report."""
    print(f"\n  Updating Typst report: {typ_path}")
    
    with open(typ_path, 'r') as f:
        content = f.read()
    
    if '7. Full-Resolution Heatmap' in content:
        # Replace existing heatmap section
        print("  Heatmap section exists, replacing...")
        # Find start and end of heatmap sections
        start_marker = '#pagebreak()\n#text(size: 18pt, weight: "bold")[7. Full-Resolution Heatmap'
        end_marker = '#pagebreak()\n#text(size: 18pt, weight: "bold")[8. Cross-Community'
        
        old_start = content.find(start_marker)
        old_end = content.find(end_marker)
        
        if old_start >= 0 and old_end >= 0:
            content = content[:old_start] + content[old_end:]
    
    # Helper to escape angle brackets for Typst
    def esc(s):
        return str(s).replace('<', '\\<').replace('>', '\\>')
    
    # Build new sections
    new_sections = f"""
#pagebreak()
#text(size: 18pt, weight: "bold")[7. Full-Resolution Heatmap (Community Detection)]
#v(0.2cm)
#text(size: 11pt)[Full 132k×132k clustered heatmap using Louvain community detection on all 36M similar pairs. \
The matrix is binned to {n_bins}×{n_bins} resolution. White lines show community boundaries. \
Dark regions indicate closely related prophage clusters. Of {n_sequences:,} sequences, \
{n_sequences - n_isolated:,} have at least one similar pair {esc('(d<0.5)')} and form {connected_comm_count} communities; \
{n_isolated:,} are singletons with no similar pairs.]
#v(0.2cm)
#image("full_heatmap.png", width: 90%)

#pagebreak()
#text(size: 18pt, weight: "bold")[8. Cross-Community Homology Matrix]
#v(0.2cm)
#text(size: 11pt)[Mean pairwise MASH distance between each pair of communities. \
Values near 0 indicate closely related communities; values near 0.5 indicate distant relationships. \
Diagonal is zero (within-community). Showing {connected_comm_count} connected communities; \
{n_isolated:,} singleton communities (no similar pairs) are omitted.]
#v(0.2cm)
#image("cluster_heatmap.png", width: 80%)
"""
    
    # Insert before the summary/ending content
    insertion_point = content.rfind('#v(1cm)')
    if insertion_point == -1:
        insertion_point = len(content)
    
    # Check if we already have section 7/8 and replace
    if '#pagebreak()\n#text(size: 18pt, weight: "bold")[7. Full-Resolution Heatmap' in content:
        # Already handled above
        pass
    else:
        content = content[:insertion_point] + new_sections + content[insertion_point:]
    
    with open(typ_path, 'w') as f:
        f.write(content)
    
    print(f"  Updated Typst report")


# ============================================================
# Main Pipeline
# ============================================================

def main():
    print("=" * 60)
    print("Full 132k Prophage Heatmap Generator")
    print("=" * 60)
    print(f"Edge list: {EDGES_FILE}")
    print(f"FASTA:     {FASTA_FILE}")
    print(f"Output:    {OUTPUT_DIR}")
    
    # Step 0: Parse all sequence names
    all_seq_names = parse_all_sequences(FASTA_FILE)
    print(f"  Total sequences: {len(all_seq_names):,}")
    
    # Step 1: Build graph from edges
    G, connected_set = build_graph(EDGES_FILE, all_seq_names)
    
    # Step 2: Community detection
    node_to_community, connected_comm_count, n_isolated = \
        run_community_detection(G, all_seq_names, connected_set)
    n_communities = connected_comm_count + n_isolated
    
    # Step 3: Order sequences
    node_to_pos, ordered_nodes = order_sequences(node_to_community, G, all_seq_names)
    
    # Save cluster membership
    print("\n  Saving cluster membership...")
    cluster_records = []
    for node, comm in node_to_community.items():
        cluster_records.append({'sequence': node, 'community': comm})
    cluster_df = pd.DataFrame(cluster_records)
    cluster_df.to_csv(CLUSTER_MEMBERSHIP_CSV, index=False)
    print(f"  Saved: {CLUSTER_MEMBERSHIP_CSV} ({len(cluster_df):,} records)")
    
    # Free graph memory
    del G
    gc.collect()
    
    # Step 4: Create binned heatmap
    bin_mean, bin_count = create_binned_heatmap(EDGES_FILE, node_to_pos, N_BINS)
    
    n_nodes = len(ordered_nodes)
    n_edges = sum(bin_count.flatten()) // 2  # Each edge counted twice
    
    # Render heatmap
    render_heatmap(
        bin_mean, ordered_nodes, n_communities,
        HEATMAP_PNG, HEATMAP_STANDALONE_PDF,
        title=f'Prophage Homology Heatmap\n'
              f'{n_nodes:,} sequences, {n_communities:,} communities, {N_BINS}×{N_BINS} bins'
    )
    
    # Step 5: Cluster-level heatmap
    create_cluster_heatmap(EDGES_FILE, node_to_community, n_communities, 
                          CLUSTER_HEATMAP_PNG, connected_comm_count)
    
    # Step 6: Update Typst report
    update_typst_report(REPORT_TYP, n_communities, n_nodes, n_edges, N_BINS,
                        connected_comm_count, n_isolated)
    
    # Compile Typst to PDF
    print("\n  Compiling Typst report to PDF...")
    report_pdf = OUTPUT_DIR / 'full_prophage_homology_report.pdf'
    result = subprocess.run(
        ['typst', 'compile', str(REPORT_TYP), str(report_pdf)],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        print(f"  WARNING: Typst compilation failed:")
        print(f"  {result.stderr[:500]}")
    else:
        file_size = os.path.getsize(report_pdf)
        print(f"  PDF generated: {report_pdf} ({file_size / 1e6:.1f} MB)")
    
    # Summary
    print("\n" + "=" * 60)
    print("Pipeline Complete!")
    print("=" * 60)
    print(f"  Heatmap PNG:        {HEATMAP_PNG}")
    print(f"  Standalone PDF:     {HEATMAP_STANDALONE_PDF}")
    print(f"  Cluster heatmap:    {CLUSTER_HEATMAP_PNG}")
    print(f"  Cluster membership: {CLUSTER_MEMBERSHIP_CSV}")
    print(f"  Typst report:       {REPORT_TYP}")
    print(f"  PDF report:         {report_pdf}")
    print(f"\n  Key stats:")
    print(f"    Total sequences:            {n_nodes:,}")
    print(f"    Connected (with edges):     {n_nodes - n_isolated:,}")
    print(f"    Isolated (no edges):        {n_isolated:,}")
    print(f"    Connected communities:      {connected_comm_count}")
    print(f"    Total communities:          {n_communities:,}")
    print(f"    Total edges:                {n_edges:,}")
    print(f"    Heatmap bins:               {N_BINS}×{N_BINS}")
    print("=" * 60)


if __name__ == '__main__':
    main()