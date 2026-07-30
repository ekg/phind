#!/usr/bin/env python3
"""
Focused heatmap of the connected component only (19,638 prophages, 12 communities).

Uses the pre-computed cluster assignments from full_heatmap_clusters.csv to
avoid re-running Louvain. Filters to the 12 real communities (0-11), builds
a 2000x2000 binned heatmap ordered by community, and generates per-community
zoom heatmaps for the 4 largest communities.

Pipeline:
1. Read cluster assignments, filter to 12 connected communities
2. Get degree (connectivity) for each prophage from the edge list
3. Build 2000x2000 binned heatmap
4. Render focused heatmap with community boundaries
5. Generate per-community zoom heatmaps for top 4
6. Update Typst report
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
MAIN_REPO = Path('/home/erikg/phind')
DATA_DIR = Path('prophage_homology_survey')
EDGES_FILE = Path('/home/erikg/phind/prophage_homology_survey/full_prophages_edges.tsv')
CLUSTERS_FILE = Path('/home/erikg/phind/prophage_homology_survey/full_heatmap_clusters.csv')
OUTPUT_DIR = DATA_DIR

# Output files
HEATMAP_PNG = OUTPUT_DIR / 'connected_heatmap.png'
HEATMAP_STANDALONE_PDF = OUTPUT_DIR / 'connected_heatmap_standalone.pdf'
REPORT_TYP = OUTPUT_DIR / 'full_prophage_homology_report.typ'

# Binning parameters
N_BINS = 2000

# Connected communities (from Louvain, 0-11 are the real communities)
CONNECTED_COMMUNITIES = set(range(12))

# Community names for display
COMMUNITY_LABELS = {
    0: 'C1', 1: 'C2', 2: 'C3', 3: 'C4', 4: 'C5',
    5: 'C6', 6: 'C7', 7: 'C8', 8: 'C9', 9: 'C10',
    10: 'C11', 11: 'C12'
}

# Plotting
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
})


def load_cluster_assignments(clusters_path):
    """Load cluster assignments from CSV, filter to connected communities."""
    print(f"\n[0] Loading cluster assignments from {clusters_path}...")
    t0 = time.time()
    df = pd.read_csv(clusters_path)
    print(f"  Total sequences: {len(df):,}")
    
    # Filter to connected communities
    connected = df[df['community'].isin(CONNECTED_COMMUNITIES)]
    print(f"  Connected sequences (communities 0-11): {len(connected):,}")
    
    # Build map: sequence -> community
    node_to_community = {}
    for _, row in connected.iterrows():
        node_to_community[row['sequence']] = int(row['community'])
    
    # Community sizes
    community_sizes = Counter(node_to_community.values())
    print(f"  Community sizes:")
    for c in sorted(community_sizes.keys()):
        size = community_sizes[c]
        print(f"    Community {c} ({COMMUNITY_LABELS[c]}): {size:,} sequences")
    
    print(f"  Time: {time.time() - t0:.1f}s")
    return node_to_community, community_sizes


def compute_degrees(edges_path, connected_nodes):
    """Compute degree for each connected prophage from the edge list."""
    print(f"\n[1] Computing degrees for connected prophages...")
    t0 = time.time()
    
    degrees = defaultdict(int)
    edge_count = 0
    connected_edge_count = 0
    report_interval = 5000000
    
    with open(edges_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                s1, s2 = parts[0], parts[1]
                edge_count += 1
                if s1 in connected_nodes and s2 in connected_nodes:
                    degrees[s1] += 1
                    degrees[s2] += 1
                    connected_edge_count += 1
                    
            if edge_count % report_interval == 0:
                elapsed = time.time() - t0
                print(f"    {edge_count:,} edges scanned ({elapsed:.1f}s)")
    
    # Ensure all connected nodes have at least degree 0
    for node in connected_nodes:
        if node not in degrees:
            degrees[node] = 0
    
    print(f"  Total edges scanned: {edge_count:,}")
    print(f"  Connected edges: {connected_edge_count:,}")
    print(f"  Connected nodes: {len(connected_nodes):,}")
    print(f"  Mean degree: {np.mean(list(degrees.values())):.1f}")
    print(f"  Max degree: {max(degrees.values())}")
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return degrees, connected_edge_count


def order_sequences(node_to_community, degrees):
    """
    Order sequences by community (by size descending), then by degree within community.
    Returns ordered list of (node, community, degree) and dict node->position.
    """
    print(f"\n[2] Ordering sequences by community + degree...")
    t0 = time.time()
    
    # Group nodes by community
    community_nodes = defaultdict(list)
    for node, comm in node_to_community.items():
        community_nodes[comm].append(node)
    
    # Sort communities by size (descending) for visual clarity
    sorted_communities = sorted(community_nodes.keys(), 
                                key=lambda c: len(community_nodes[c]), 
                                reverse=True)
    print(f"  Community order (by size): {[COMMUNITY_LABELS[c] for c in sorted_communities]}")
    
    # Build ordered list
    ordered_nodes = []
    node_to_pos = {}
    pos = 0
    
    for comm in sorted_communities:
        nodes = community_nodes[comm]
        # Sort by degree within community (descending)
        nodes_sorted = sorted(nodes, key=lambda n: degrees.get(n, 0), reverse=True)
        for n in nodes_sorted:
            node_to_pos[n] = pos
            ordered_nodes.append((n, comm, degrees.get(n, 0)))
            pos += 1
    
    assert len(ordered_nodes) == len(node_to_community), \
        f"Ordered {len(ordered_nodes)} but expected {len(node_to_community)}"
    
    print(f"  Ordered {len(ordered_nodes):,} sequences in {len(sorted_communities)} communities")
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return node_to_pos, ordered_nodes, sorted_communities


def create_binned_heatmap(edges_path, connected_nodes, node_to_pos, n_bins=2000):
    """
    Create a 2D binned heatmap from edges between connected prophages.
    Each cell shows the mean distance for all pairs in that bin region.
    Bins with no edges get distance = 1.0 (max distance = no similarity).
    """
    print(f"\n[3] Creating binned heatmap ({n_bins}×{n_bins})...")
    t0 = time.time()
    n_nodes = len(node_to_pos)
    
    # Initialize bin accumulators
    bin_sum = np.zeros((n_bins, n_bins), dtype=np.float64)
    bin_count = np.zeros((n_bins, n_bins), dtype=np.int64)
    
    def pos_to_bin(p):
        return min(int(p * n_bins / n_nodes), n_bins - 1)
    
    print(f"  Processing {n_nodes:,} sequences across {n_bins}×{n_bins} bins...")
    
    edge_count = 0
    valid_edge_count = 0
    report_interval = 5000000
    
    with open(edges_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                s1, s2, dist = parts[0], parts[1], float(parts[2])
                edge_count += 1
                
                # Only process edges between connected prophages
                if s1 not in connected_nodes or s2 not in connected_nodes:
                    continue
                    
                pos1 = node_to_pos.get(s1)
                pos2 = node_to_pos.get(s2)
                if pos1 is None or pos2 is None:
                    continue
                
                valid_edge_count += 1
                
                i = pos_to_bin(pos1)
                j = pos_to_bin(pos2)
                
                bin_sum[i, j] += dist
                bin_sum[j, i] += dist
                bin_count[i, j] += 1
                bin_count[j, i] += 1
                
                if edge_count % report_interval == 0:
                    elapsed = time.time() - t0
                    print(f"    {edge_count:,} edges scanned ({elapsed:.1f}s)")
    
    print(f"  Total edges scanned: {edge_count:,}")
    print(f"  Valid edges (both in connected set): {valid_edge_count:,}")
    
    # Compute mean distance for each bin
    with np.errstate(divide='ignore', invalid='ignore'):
        bin_mean = np.where(bin_count > 0, bin_sum / bin_count, 1.0)
    
    filled_bins = np.sum(bin_count > 0)
    print(f"  Bins with data: {filled_bins:,} / {n_bins * n_bins:,} "
          f"({100 * filled_bins / (n_bins * n_bins):.1f}%)")
    print(f"  Mean distance in filled bins: {np.mean(bin_mean[bin_count > 0]):.4f}")
    print(f"  Time: {time.time() - t0:.1f}s")
    
    return bin_mean, bin_count


def render_heatmap(bin_mean, ordered_nodes, sorted_communities, community_sizes,
                   output_png, output_pdf=None, title=None):
    """Render the binned heatmap as PNG and optional PDF with community boundaries."""
    print(f"  Rendering heatmap...")
    t0 = time.time()
    
    n_bins = bin_mean.shape[0]
    n_nodes = len(ordered_nodes)
    
    # Compute community boundaries in bin coordinates
    # Build cumulative positions
    community_boundaries = []
    cum_pos = 0
    for comm in sorted_communities:
        size = community_sizes[comm]
        cum_pos += size
        # Position in bin space
        bin_pos = cum_pos * n_bins / n_nodes
        community_boundaries.append((bin_pos, comm, size))
    
    # Color palette for community boundary lines (subtle colors)
    boundary_colors = ['#e74c3c', '#e67e22', '#f1c40f', '#2ecc71', '#1abc9c',
                       '#3498db', '#9b59b6', '#e91e63', '#00bcd4', '#ff5722',
                       '#795548', '#607d8b']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 9))
    
    vmin, vmax = 0.0, 0.5
    im = ax.imshow(bin_mean, aspect='equal', cmap='inferno', 
                   vmin=vmin, vmax=vmax, interpolation='nearest',
                   extent=[0, n_bins, 0, n_bins])
    
    # Add community boundary lines
    prev_pos = 0
    for i, (bin_pos, comm, size) in enumerate(community_boundaries):
        if size > 1:
            color = boundary_colors[i % len(boundary_colors)]
            ax.axhline(bin_pos, color=color, linewidth=0.6, alpha=0.5)
            ax.axvline(bin_pos, color=color, linewidth=0.6, alpha=0.5)
            ax.axhline(prev_pos, color=color, linewidth=0.6, alpha=0.5, linestyle='--')
            ax.axvline(prev_pos, color=color, linewidth=0.6, alpha=0.5, linestyle='--')
        prev_pos = bin_pos
    
    # Add community labels
    prev_pos = 0
    for i, (bin_pos, comm, size) in enumerate(community_boundaries):
        mid_pos = (prev_pos + bin_pos) / 2
        label = f"{COMMUNITY_LABELS[comm]}\n({size:,})"
        ax.text(n_bins * 1.01, mid_pos, label, fontsize=6, va='center', ha='left',
                color=boundary_colors[i % len(boundary_colors)])
        ax.text(mid_pos, n_bins * 1.01, label, fontsize=6, ha='center', va='bottom',
                color=boundary_colors[i % len(boundary_colors)])
        prev_pos = bin_pos
    
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.08)
    cbar.set_label('MASH Distance', fontsize=10)
    
    ax.set_xlabel('Sequence (ordered by community)', fontsize=10)
    ax.set_ylabel('Sequence (ordered by community)', fontsize=10)
    
    n_communities = len(sorted_communities)
    
    if title:
        ax.set_title(title, fontsize=11, fontweight='bold')
    else:
        ax.set_title(f'Connected Prophage Homology Heatmap\n'
                     f'{n_nodes:,} sequences, {n_communities} communities, '
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
        im_pdf = ax_pdf.imshow(bin_mean, aspect='equal', cmap='inferno',
                               vmin=vmin, vmax=vmax, interpolation='nearest',
                               extent=[0, n_bins, 0, n_bins])
        
        # Add community boundary lines
        prev_pos = 0
        for i, (bin_pos, comm, size) in enumerate(community_boundaries):
            if size > 1:
                color = boundary_colors[i % len(boundary_colors)]
                ax_pdf.axhline(bin_pos, color=color, linewidth=0.6, alpha=0.5)
                ax_pdf.axvline(bin_pos, color=color, linewidth=0.6, alpha=0.5)
                ax_pdf.axhline(prev_pos, color=color, linewidth=0.6, alpha=0.5, linestyle='--')
                ax_pdf.axvline(prev_pos, color=color, linewidth=0.6, alpha=0.5, linestyle='--')
            prev_pos = bin_pos
        
        # Add community labels
        prev_pos = 0
        for i, (bin_pos, comm, size) in enumerate(community_boundaries):
            mid_pos = (prev_pos + bin_pos) / 2
            label = f"{COMMUNITY_LABELS[comm]} ({size:,})"
            ax_pdf.text(n_bins * 1.01, mid_pos, label, fontsize=6, va='center', ha='left',
                        color=boundary_colors[i % len(boundary_colors)])
            prev_pos = bin_pos
        
        cbar_pdf = fig_pdf.colorbar(im_pdf, ax=ax_pdf, shrink=0.8, pad=0.08)
        cbar_pdf.set_label('MASH Distance', fontsize=10)
        
        ax_pdf.set_xlabel('Sequence (ordered by community)', fontsize=10)
        ax_pdf.set_ylabel('Sequence (ordered by community)', fontsize=10)
        ax_pdf.set_title(f'Connected Prophage Homology (n={n_nodes:,})', 
                         fontsize=11, fontweight='bold')
        
        plt.tight_layout()
        fig_pdf.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.2)
        print(f"  Saved: {output_pdf}")
        plt.close(fig_pdf)
    
    print(f"  Render time: {time.time() - t0:.1f}s")


def create_community_zoom(edges_path, connected_nodes, node_to_community, 
                          community_id, community_size, output_png, node_to_pos_all=None):
    """
    Create a per-community zoom heatmap for a single community.
    Shows only the edges within that community.
    """
    community_label = COMMUNITY_LABELS.get(community_id, f'C{community_id+1}')
    print(f"\n[4.{community_id}] Creating zoom heatmap for {community_label} "
          f"({community_size:,} members)...")
    t0 = time.time()
    
    # Get all nodes in this community
    community_nodes = set()
    for node, comm in node_to_community.items():
        if comm == community_id:
            community_nodes.add(node)
    
    assert len(community_nodes) == community_size, \
        f"Expected {community_size} nodes, got {len(community_nodes)}"
    
    # Build positions within this community only
    n_local = len(community_nodes)
    local_positions = {}
    for i, node in enumerate(sorted(community_nodes)):
        local_positions[node] = i
    
    # Determine zoom bin size
    if n_local < 100:
        zoom_bins = n_local  # 1:1 for small communities
    elif n_local < 500:
        zoom_bins = n_local // 2
    elif n_local < 2000:
        zoom_bins = 500
    elif n_local < 4000:
        zoom_bins = 800
    else:
        zoom_bins = 1000
    
    # Initialize bin accumulators
    bin_sum = np.zeros((zoom_bins, zoom_bins), dtype=np.float64)
    bin_count = np.zeros((zoom_bins, zoom_bins), dtype=np.int64)
    
    def pos_to_bin(p, n_bins_local):
        return min(int(p * n_bins_local / n_local), n_bins_local - 1)
    
    # Stream edges and filter for this community
    edge_count = 0
    valid_edge_count = 0
    
    with open(edges_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                s1, s2, dist = parts[0], parts[1], float(parts[2])
                edge_count += 1
                
                if s1 not in community_nodes or s2 not in community_nodes:
                    continue
                    
                pos1 = local_positions.get(s1)
                pos2 = local_positions.get(s2)
                if pos1 is None or pos2 is None:
                    continue
                
                valid_edge_count += 1
                
                i = pos_to_bin(pos1, zoom_bins)
                j = pos_to_bin(pos2, zoom_bins)
                
                bin_sum[i, j] += dist
                bin_sum[j, i] += dist
                bin_count[i, j] += 1
                bin_count[j, i] += 1
    
    # Compute mean
    with np.errstate(divide='ignore', invalid='ignore'):
        bin_mean = np.where(bin_count > 0, bin_sum / bin_count, 1.0)
    
    filled_bins = np.sum(bin_count > 0)
    print(f"  Community edges: {valid_edge_count:,}")
    print(f"  Bins: {zoom_bins}×{zoom_bins}, filled: {filled_bins:,} / {zoom_bins*zoom_bins:,}")
    
    # Render
    fig, ax = plt.subplots(figsize=(8, 7))
    
    vmin, vmax = 0.0, 0.5
    im = ax.imshow(bin_mean, aspect='equal', cmap='inferno',
                   vmin=vmin, vmax=vmax, interpolation='nearest',
                   extent=[0, zoom_bins, 0, zoom_bins])
    
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('MASH Distance', fontsize=10)
    
    ax.set_xlabel('Sequence (within community)', fontsize=10)
    ax.set_ylabel('Sequence (within community)', fontsize=10)
    ax.set_title(f'Community {community_label} Zoom\n'
                 f'{n_local:,} prophages, {valid_edge_count:,} edges, '
                 f'{zoom_bins}×{zoom_bins} bins',
                 fontsize=11, fontweight='bold')
    
    # Add stats
    ax.text(0.02, 0.98, f'{n_local:,} sequences',
            transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches='tight', pad_inches=0.1)
    print(f"  Saved: {output_png}")
    plt.close(fig)
    print(f"  Time: {time.time() - t0:.1f}s")


def update_typst_report(typ_path, n_nodes, n_edges, n_bins, n_communities,
                        community_sizes, top4_communities):
    """Update the Typst report with the focused heatmap."""
    print(f"\n[5] Updating Typst report: {typ_path}")
    t0 = time.time()
    
    with open(typ_path, 'r') as f:
        content = f.read()
    
    def esc(s):
        return str(s).replace('<', '\\<').replace('>', '\\>')
    
    # Build the new section 7 content (focused heatmap)
    # We'll insert it after the existing section 7 (full heatmap) and
    # renumber old section 8 to 9, etc.
    
    # Describe community sizes
    comm_sizes_str = ', '.join(
        f'{COMMUNITY_LABELS[c]}: {community_sizes[c]:,}' 
        for c in sorted(community_sizes.keys(), key=lambda c: community_sizes[c], reverse=True)
    )
    
    top4_zoom_section = ""
    for i, comm_id in enumerate(top4_communities):
        size = community_sizes[comm_id]
        label = COMMUNITY_LABELS[comm_id]
        top4_zoom_section += f"""
#v(0.3cm)
#image("community_zoom_{comm_id}.png", width: 80%)
#v(0.1cm)
#text(size: 9pt)[*Community {label}* ({size:,} prophages) — zoom heatmap showing intra-community homology structure.]
"""
    
    # Build the focused heatmap section
    focused_section = f"""
#pagebreak()
#text(size: 18pt, weight: "bold")[8. Focused Heatmap: Connected Component Only]
#v(0.2cm)
#text(size: 11pt)[Focused heatmap of the {n_nodes:,} prophages that form the connected component \
(those with at least one similar pair at {esc('d<0.5')}). These {n_nodes:,} sequences form \
{n_communities} communities. The {n_nodes - sum(community_sizes[c] for c in top4_communities):,} \
remaining prophages form 8 smaller communities. The matrix is binned to {n_bins}×{n_bins} resolution. \
Color community boundaries are shown. Community sizes: {comm_sizes_str}.]
#v(0.2cm)
#image("connected_heatmap.png", width: 90%)

#pagebreak()
#text(size: 18pt, weight: "bold")[9. Per-Community Zoom Heatmaps]
#v(0.2cm)
#text(size: 11pt)[Zoom heatmaps for the 4 largest communities, showing detailed intra-community homology structure. \
Each zoom uses adaptive binning to reveal fine-scale relationships.]
#v(0.2cm)
{top4_zoom_section}
"""
    
    # Find the position to insert — before the final "Generated by" line
    generation_marker = content.rfind('Generated by MASH')
    if generation_marker > 0:
        gen_line_start = content.rfind('\n', 0, generation_marker) + 1
        insert_pos = gen_line_start
        content = content[:insert_pos] + focused_section + '\n' + content[insert_pos:]
    else:
        # Append at end as fallback
        content += focused_section
    
    with open(typ_path, 'w') as f:
        f.write(content)
    
    print(f"  Updated Typst report")
    print(f"  Time: {time.time() - t0:.1f}s")


def main():
    print("=" * 60)
    print("Connected Component Prophage Heatmap Generator")
    print("=" * 60)
    print(f"Edge list:  {EDGES_FILE}")
    print(f"Clusters:   {CLUSTERS_FILE}")
    print(f"Output:     {OUTPUT_DIR}")
    
    # Step 0: Load cluster assignments
    node_to_community, community_sizes = load_cluster_assignments(CLUSTERS_FILE)
    connected_nodes = set(node_to_community.keys())
    n_connected = len(connected_nodes)
    
    # Step 1: Compute degrees
    degrees, connected_edge_count = compute_degrees(EDGES_FILE, connected_nodes)
    
    # Step 2: Order sequences
    node_to_pos, ordered_nodes, sorted_communities = order_sequences(
        node_to_community, degrees
    )
    
    # Step 3: Create binned heatmap
    bin_mean, bin_count = create_binned_heatmap(
        EDGES_FILE, connected_nodes, node_to_pos, N_BINS
    )
    
    # Step 4: Render heatmap
    n_communities = len(sorted_communities)
    render_heatmap(
        bin_mean, ordered_nodes, sorted_communities, community_sizes,
        HEATMAP_PNG, HEATMAP_STANDALONE_PDF,
        title=f'Connected Prophage Homology Heatmap\n'
              f'{n_connected:,} sequences, {n_communities} communities, {N_BINS}×{N_BINS} bins'
    )
    
    # Step 5: Per-community zoom heatmaps for top 4
    # Top 4 by size: community 10 (5208), 0 (4999), 4 (3030), 2 (2826)
    top4 = sorted(community_sizes.keys(), key=lambda c: community_sizes[c], reverse=True)[:4]
    print(f"\n  Top 4 communities for zoom: {[COMMUNITY_LABELS[c] for c in top4]}")
    print(f"  Sizes: {[community_sizes[c] for c in top4]}")
    
    for comm_id in top4:
        zoom_png = OUTPUT_DIR / f'community_zoom_{comm_id}.png'
        create_community_zoom(
            EDGES_FILE, connected_nodes, node_to_community,
            comm_id, community_sizes[comm_id], zoom_png
        )
    
    # Step 6: Update Typst report
    update_typst_report(
        REPORT_TYP, n_connected, connected_edge_count, N_BINS,
        n_communities, community_sizes, top4
    )
    
    # Compile Typst to PDF
    print("\n[6] Compiling Typst report to PDF...")
    report_pdf = OUTPUT_DIR / 'full_prophage_homology_report.pdf'
    result = subprocess.run(
        ['typst', 'compile', str(REPORT_TYP), str(report_pdf)],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        print(f"  WARNING: Typst compilation failed:")
        print(f"  {result.stderr[:1000]}")
    else:
        file_size = os.path.getsize(report_pdf)
        print(f"  PDF generated: {report_pdf} ({file_size / 1e6:.1f} MB)")
    
    # Summary
    print("\n" + "=" * 60)
    print("Pipeline Complete!")
    print("=" * 60)
    print(f"  Focused heatmap PNG:     {HEATMAP_PNG}")
    print(f"  Standalone PDF:          {HEATMAP_STANDALONE_PDF}")
    for comm_id in top4:
        zoom_png = OUTPUT_DIR / f'community_zoom_{comm_id}.png'
        print(f"  Zoom C{comm_id}:              {zoom_png}")
    print(f"  Typst report:            {REPORT_TYP}")
    print(f"  PDF report:              {report_pdf}")
    print(f"\n  Key stats:")
    print(f"    Connected sequences:     {n_connected:,}")
    print(f"    Connected communities:   {n_communities}")
    print(f"    Connected edges:         {connected_edge_count:,}")
    print(f"    Heatmap bins:            {N_BINS}×{N_BINS}")
    print("=" * 60)


if __name__ == '__main__':
    main()