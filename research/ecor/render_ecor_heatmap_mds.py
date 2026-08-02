#!/usr/bin/env python3
"""
render_ecor_heatmap_mds.py — ECOR-highlighted heatmap / MDS figures.

ECOR prophage elements are all *singleton* Louvain communities in the MASH
graph (no connected-component neighbours, verified: 0/300 in communities 0-11),
so the existing connected-component heatmap cannot show them. Instead we
highlight ECOR elements in the spaces that do resolve them:

  ecor_heatmap.png + ecor_heatmap_interactive.html
      300 x 300 cophenetic distance matrix (all ECOR elements, tree-path
      distances in the all-prophage UPGMA tree), rows/cols dendrogram-ordered,
      strain labels — "ECOR in tree space".
  ecor_mds.png + ecor_mds_interactive.html
      Metric MDS (sklearn) of the same 300x300 cophenetic matrix.
  ecor_mds_subset_highlighted.png + ecor_mds_subset_interactive.html
      Reuse of the survey's existing 5,000-prophage MDS
      (prophage_homology_survey/full_prophage_mds_coords.csv) with the
      7 ECOR elements present in that subset marked in red.

Usage:
  python3 render_ecor_heatmap_mds.py
"""

import csv
import json
import os
import time
from pathlib import Path

import numpy as np
import ete3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn import manifold

HERE = Path(__file__).resolve().parent
MASH = Path('/home/erikg/phind/research/mash_tree')
SURVEY = Path('/home/erikg/phind/prophage_homology_survey')

TREE_NWK = MASH / 'full_prophages_tree.nwk'
TAGS = HERE / 'ecor_leaf_tags.tsv'
INDEX = HERE / 'ecor_inspection_index.csv'
SUBSET_MDS = SURVEY / 'full_prophage_mds_coords.csv'
CONN_CLUSTERS = SURVEY / 'full_heatmap_clusters.csv'

HEAT_PNG = HERE / 'ecor_heatmap.png'
HEAT_HTML = HERE / 'ecor_heatmap_interactive.html'
MDS_PNG = HERE / 'ecor_mds.png'
MDS_HTML = HERE / 'ecor_mds_interactive.html'
SUBSET_PNG = HERE / 'ecor_mds_subset_highlighted.png'
SUBSET_HTML = HERE / 'ecor_mds_subset_interactive.html'

ECOR_COLOR = '#d62728'


def load_tags():
    tags = {}
    with open(TAGS) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['is_ecor'] == 'TRUE':
                tags[row['prophage_id']] = row['ecor_strain']
    return tags


def load_index():
    rows = {}
    with open(INDEX) as f:
        for row in csv.DictReader(f):
            rows[row['prophage_id']] = row
    return rows


def cophenetic_matrix(leaves, pairs):
    """300x300 pairwise tree-path distances. pairs = [(leaf_a, leaf_b), ...]."""
    n = len(leaves)
    D = np.zeros((n, n))
    for i in range(n):
        D[i, i] = 0.0
        for j in range(i + 1, n):
            d = leaves[i].get_distance(leaves[j])
            D[i, j] = d
            D[j, i] = d
    return D


def render_heatmap(D, pids, strains, idx_rows, png_path, html_path):
    import plotly.graph_objects as go

    n = len(pids)
    # dendrogram order
    Z = linkage(D, method='average')
    order = leaves_list(Z)
    Ds = D[np.ix_(order, order)]
    labels = [pids[i] for i in order]
    strain_labels = [strains[p] for p in labels]

    # --- static PNG ---
    fig, ax = plt.subplots(figsize=(13, 12), dpi=140)
    im = ax.imshow(Ds, cmap='viridis_r', aspect='auto', vmin=0, vmax=1.0)
    # strain group separators
    prev = None
    for i, s in enumerate(strain_labels):
        if s != prev and i > 0:
            ax.axhline(i - 0.5, color='white', lw=0.6)
            ax.axvline(i - 0.5, color='white', lw=0.6)
        prev = s
    ax.set_xticks(range(n))
    ax.set_xticklabels(strain_labels, rotation=90, fontsize=5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(strain_labels, fontsize=5)
    ax.set_xlabel('ECOR element (dendrogram order)')
    ax.set_ylabel('ECOR element')
    ax.set_title('ECOR × ECOR cophenetic distance in the all-prophage MASH UPGMA tree '
                 f'({n} elements)', fontsize=11)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03)
    cbar.set_label('tree-path distance (UPGMA height)')
    fig.tight_layout()
    fig.savefig(png_path, bbox_inches='tight')
    plt.close(fig)
    print(f'[heatmap] wrote {png_path.name} ({os.path.getsize(png_path)/1e6:.1f} MB)')

    # --- interactive HTML (plotly) ---
    fig = go.Figure(go.Heatmap(
        z=Ds, x=strain_labels, y=strain_labels,
        colorscale='Viridis', zmin=0, zmax=1.0,
        hovertemplate='%{x}<br>%{y}<br>d = %{z:.4f}<extra></extra>',
    ))
    fig.update_layout(
        title=dict(text=f'ECOR × ECOR cophenetic distance, all-prophage MASH tree ({n} elements)',
                   font=dict(size=14)),
        height=900, xaxis=dict(title='ECOR element', tickangle=-90, tickfont=dict(size=8)),
        yaxis=dict(title='ECOR element', tickfont=dict(size=8), autorange='reversed'),
        margin=dict(l=10, r=10, t=60, b=130),
    )
    fig.write_html(html_path, include_plotlyjs='cdn', config={'displaylogo': False})
    print(f'[heatmap] wrote {html_path.name}')


def render_mds(D, pids, strains, idx_rows, png_path, html_path):
    import plotly.graph_objects as go

    n = len(pids)
    mds = manifold.MDS(n_components=2, dissimilarity='precomputed',
                       random_state=42, normalized_stress='auto', max_iter=500)
    X = mds.fit_transform(D)
    xs, ys = X[:, 0], X[:, 1]

    # --- static PNG ---
    fig, ax = plt.subplots(figsize=(11, 10), dpi=140)
    ax.scatter(xs, ys, s=22, c='#7f8c9b', alpha=0.75, edgecolors='none')
    ax.scatter(xs, ys, s=22, c=ECOR_COLOR, alpha=0.95, edgecolors='none')
    for i in range(n):
        ax.annotate(strains[pids[i]], (xs[i], ys[i]), fontsize=6,
                    ha='center', va='bottom', color='#333')
    ax.set_title(f'Metric MDS of ECOR elements in all-prophage MASH tree space ({n} elements; '
                 f'cophenetic distances)', fontsize=11)
    ax.set_xlabel('MDS1')
    ax.set_ylabel('MDS2')
    fig.tight_layout()
    fig.savefig(png_path, bbox_inches='tight')
    plt.close(fig)
    print(f'[mds] wrote {png_path.name}')

    # --- interactive HTML ---
    hover = [f'{strains[p]}<br>{p}<br>d(nearest) = {min(D[i, j] for j in range(n) if j != i):.4f}'
             for i, p in enumerate(pids)]
    fig = go.Figure(go.Scatter(
        x=xs, y=ys, mode='markers+text', text=[strains[p] for p in pids],
        textposition='top center', textfont=dict(size=9),
        marker=dict(size=9, color=ECOR_COLOR),
        customdata=hover, hovertemplate='%{customdata}<extra></extra>',
    ))
    fig.update_layout(
        title=dict(text=f'Metric MDS of {n} ECOR prophage elements (cophenetic distance, '
                        'all-prophage MASH tree)', font=dict(size=14)),
        height=800, margin=dict(l=10, r=10, t=60, b=40),
        xaxis=dict(title='MDS1'), yaxis=dict(title='MDS2'),
    )
    fig.write_html(html_path, include_plotlyjs='cdn', config={'displaylogo': False})
    print(f'[mds] wrote {html_path.name}')


def render_subset_mds():
    """Reuse the survey's 5,000-prophage MDS; mark the ECOR elements in it."""
    import plotly.graph_objects as go

    tags = load_tags()
    df_rows = []
    with open(SUBSET_MDS) as f:
        for row in csv.DictReader(f):
            df_rows.append(row)
    xs = np.array([float(r['MDS1']) for r in df_rows])
    ys = np.array([float(r['MDS2']) for r in df_rows])
    pids = [r['sequence'] for r in df_rows]

    is_ecor = [p in tags for p in pids]
    ecor_pids = [p for p, e in zip(pids, is_ecor) if e]
    print(f'[subset-mds] {len(pids)} prophages, {len(ecor_pids)} ECOR in subset: '
          f'{ecor_pids[:10]}...')

    fig, ax = plt.subplots(figsize=(12, 9), dpi=140)
    ax.scatter(xs[~np.array(is_ecor)], ys[~np.array(is_ecor)], s=3,
               c='#aab6bf', alpha=0.6, label=f'other prophages ({len(pids) - len(ecor_pids):,})')
    ax.scatter(xs[np.array(is_ecor)], ys[np.array(is_ecor)], s=42, c=ECOR_COLOR,
               edgecolors='#7a1d1d', linewidths=0.5,
               label=f'ECOR prophage ({len(ecor_pids)})')
    for p in ecor_pids:
        i = pids.index(p)
        ax.annotate(tags[p], (xs[i], ys[i]), fontsize=8, ha='center', va='bottom',
                    color=ECOR_COLOR, fontweight='bold')
    ax.set_title('Existing survey MDS (5,000-prophage subset) with ECOR elements marked',
                 fontsize=11)
    ax.set_xlabel('MDS1')
    ax.set_ylabel('MDS2')
    ax.legend(loc='best', fontsize=9)
    fig.tight_layout()
    fig.savefig(SUBSET_PNG, bbox_inches='tight')
    plt.close(fig)
    print(f'[subset-mds] wrote {SUBSET_PNG.name}')

    # interactive HTML
    hover = [f'{tags.get(p, "")}<br>{p}' if e else p
             for p, e in zip(pids, is_ecor)]
    colors = [ECOR_COLOR if e else '#aab6bf' for e in is_ecor]
    sizes = [12 if e else 4 for e in is_ecor]
    fig = go.Figure(go.Scatter(
        x=xs, y=ys, mode='markers', marker=dict(size=sizes, color=colors, opacity=0.85),
        customdata=hover, hovertemplate='%{customdata}<extra></extra>',
    ))
    fig.update_layout(
        title=dict(text='Survey MDS (5,000-prophage subset) — ECOR elements marked red',
                   font=dict(size=14)),
        height=800, margin=dict(l=10, r=10, t=60, b=40),
        xaxis=dict(title='MDS1'), yaxis=dict(title='MDS2'),
    )
    fig.write_html(SUBSET_HTML, include_plotlyjs='cdn', config={'displaylogo': False})
    print(f'[subset-mds] wrote {SUBSET_HTML.name}')


def main():
    tags = load_tags()
    idx_rows = load_index()
    pids = sorted(tags.keys())
    n = len(pids)
    print(f'[load] {n} ECOR elements')

    print('[tree] parsing ...')
    t0 = time.time()
    tree = ete3.Tree(str(TREE_NWK), format=1)
    name_idx = {}
    for node in tree.traverse('preorder'):
        name_idx[node.name] = node
    leaves = [name_idx[p] for p in pids]
    print(f'       parsed in {time.time()-t0:.1f}s')

    print('[matrix] cophenetic distances ...')
    t0 = time.time()
    D = cophenetic_matrix(leaves, None)
    print(f'       {n}x{n} matrix in {time.time()-t0:.1f}s')

    render_heatmap(D, pids, tags, idx_rows, HEAT_PNG, HEAT_HTML)
    render_mds(D, pids, tags, idx_rows, MDS_PNG, MDS_HTML)
    render_subset_mds()


if __name__ == '__main__':
    main()
