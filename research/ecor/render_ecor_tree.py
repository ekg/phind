#!/usr/bin/env python3
"""
render_ecor_tree.py — ECOR-highlighted tree visualizations.

Renders:
  ecor_tree_highlighted.png          full 132,393-leaf MASH tree, ECOR tips red
                                     (static, tall, zoomable image)
  ecor_tree_interactive.html         interactive *skeleton* tree: the minimal
                                     subtree spanning all 300 ECOR leaves, with
                                     non-ECOR-only clades collapsed to triangles
                                     (plotly, ECOR leaves labelled by strain)
  ecor_neighborhood_explorer.html    dropdown explorer: per-ECOR local clade
                                     (max <=60 leaves) rendered as an interactive
                                     tree (plotly)
  ecor_neighborhoods/<prophage_id>.png   per-ECOR static neighbourhood tree image

Inputs: research/mash_tree/full_prophages_tree.nwk, ecor_leaf_tags.tsv,
        ecor_inspection_index.csv (for neighbourhood metadata).

Usage:
  python3 render_ecor_tree.py [--full] [--skeleton] [--explorer] [--neighborhoods] [--all]
"""

import csv
import html
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import ete3
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
MASH = Path('/home/erikg/phind/research/mash_tree')
TREE_NWK = MASH / 'full_prophages_tree.nwk'
TAGS = HERE / 'ecor_leaf_tags.tsv'
INDEX = HERE / 'ecor_inspection_index.csv'
NBR_DIR = HERE / 'ecor_neighborhoods'
SKELETON_HTML = HERE / 'ecor_tree_interactive.html'
EXPLORER_HTML = HERE / 'ecor_neighborhood_explorer.html'
FULL_PNG = HERE / 'ecor_tree_highlighted.png'

NBR_MAX = 60

ECOR_COLOR = '#d62728'
NON_ECOR_COLOR = '#9aa4ad'
SKELETON_EDGE = '#5b6770'


# ------------------------------------------------------------------ layout
def layout(root):
    """Rectangular tree layout.
    Returns (leaf_order, coords, parent_map) where:
      coords[id(node)] = (x, y);  x = distance from root, y = leaf rank
      parent_map[id(node)] = id(parent) or None
    Iterative (no Python recursion; safe for 132k-leaf trees)."""
    # postorder
    stack = [(root, False)]
    post = []
    while stack:
        node, vis = stack.pop()
        if vis:
            post.append(node)
        else:
            stack.append((node, True))
            for c in reversed(node.children):
                stack.append((c, False))

    counts = {}
    for node in post:
        if node.is_leaf():
            counts[id(node)] = 1
        else:
            counts[id(node)] = sum(counts[id(c)] for c in node.children)

    x = {}
    parent = {}
    stack = [(root, None, 0.0)]
    while stack:
        node, par, xv = stack.pop()
        x[id(node)] = xv
        parent[id(node)] = par
        if node.is_leaf():
            pass
        else:
            for c in reversed(node.children):
                stack.append((c, id(node), xv + c.dist))

    y = {}
    leaf_order = []
    rank = 0
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            y[id(node)] = float(rank)
            leaf_order.append(node)
            rank += 1
        else:
            for c in reversed(node.children):
                stack.append(c)
    # internal y = midpoint of children extents (subtree contiguous in leaf order)
    for node in post:
        if not node.is_leaf():
            lo = min(y[id(c)] for c in node.children)
            hi = max(y[id(c)] for c in node.children)
            y[id(node)] = (lo + hi) / 2.0

    coords = {k: (x[k], y[k]) for k in x}
    return leaf_order, coords, parent


def edge_segments(root, coords, parent):
    segs = []
    for node in root.traverse('preorder'):
        pid = parent.get(id(node))
        if pid is None:
            continue
        x1, y1 = coords[pid]
        x2, y2 = coords[id(node)]
        segs.append(((x1, y1), (x2, y2)))
    return segs


def edge_segments_iter(nodes, coords, parent):
    """Segments for an arbitrary list of nodes sharing one coords/parent space."""
    segs = []
    for node in nodes:
        pid = parent.get(id(node))
        if pid is None:
            continue
        segs.append((coords[pid], coords[id(node)]))
    return segs


# ------------------------------------------------------------- data loading
def load_tags():
    tags = {}
    with open(TAGS) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['is_ecor'] == 'TRUE':
                tags[row['prophage_id']] = row['ecor_strain']
    return tags


def load_index():
    idx = {}
    with open(INDEX) as f:
        for row in csv.DictReader(f):
            idx[row['prophage_id']] = row
    return idx


def build_name_index(tree):
    idx = {}
    for node in tree.traverse('preorder'):
        idx[node.name] = node
    return idx


# ------------------------------------------------------------- full tree PNG
def render_full_tree_png(tree, name_idx, ecor_ids, out_path):
    print('[full] layout of full tree ...')
    t0 = time.time()
    leaf_order, coords, parent = layout(tree)
    print(f'       layout done in {time.time()-t0:.1f}s')

    ecor_set = set(ecor_ids)
    n = len(leaf_order)

    # edges: color ECOR leaf edges red
    ecor_edges = []
    plain_edges = []
    for node in tree.traverse('preorder'):
        pid = parent.get(id(node))
        if pid is None:
            continue
        seg = (coords[pid], coords[id(node)])
        (ecor_edges if node.name in ecor_set else plain_edges).append(seg)

    fig_h = max(24.0, n / 1800.0)  # ~73 inches for 132k leaves
    fig, ax = plt.subplots(figsize=(9, fig_h), dpi=180)
    if plain_edges:
        ax.add_collection(LineCollection(plain_edges, colors=NON_ECOR_COLOR,
                                         linewidths=0.4, alpha=0.6))
    if ecor_edges:
        ax.add_collection(LineCollection(ecor_edges, colors=ECOR_COLOR,
                                         linewidths=1.1, alpha=1.0))
    # ECOR tip markers
    ex, ey = [], []
    for pid in ecor_ids:
        node = name_idx.get(pid)
        if node is None:
            continue
        px, py = coords[id(node)]
        ex.append(px)
        ey.append(py)
    ax.scatter(ex, ey, s=2.0, c=ECOR_COLOR, marker='s', zorder=5)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-1, n)
    ax.set_title(f'All-prophage MASH UPGMA tree (n={n:,} leaves) — ECOR elements highlighted (red, n={len(ecor_ids)})',
                 fontsize=12, loc='left')
    ax.set_xlabel('UPGMA height (MASH distance)')
    ax.set_ylabel(f'leaf order (1..{n:,})')
    ax.legend(handles=[
        Line2D([0], [0], color=ECOR_COLOR, lw=2, label=f'ECOR prophage ({len(ecor_ids)})'),
        Line2D([0], [0], color=NON_ECOR_COLOR, lw=2, label=f'other prophages ({n - len(ecor_ids):,})'),
    ], loc='upper right', fontsize=11, framealpha=0.9)
    ax.set_yticks([])
    fig.savefig(out_path, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    print(f'       wrote {out_path.name} ({os.path.getsize(out_path)/1e6:.1f} MB)')

    # companion: ECOR tip-position strip (leaf order vs count)
    strip_path = out_path.with_name('ecor_tree_tip_strip.png')
    ranks = [coords[id(name_idx[pid])][1] for pid in ecor_ids
             if name_idx.get(pid) is not None]
    fig, ax = plt.subplots(figsize=(12, 2.6), dpi=150)
    ax.scatter(ranks, [0.5] * len(ranks), s=1.5, c=ECOR_COLOR, marker='|', zorder=3)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, n)
    ax.set_yticks([])
    ax.set_xlabel(f'leaf order (1..{n:,}); 300 red ticks = ECOR element positions')
    ax.set_title('ECOR element positions across the all-prophage MASH tree (leaf order)',
                 fontsize=11, loc='left')
    fig.tight_layout()
    fig.savefig(strip_path, bbox_inches='tight')
    plt.close(fig)
    print(f'       wrote {strip_path.name}')


# ------------------------------------------------------------- skeleton tree
class SNode:
    __slots__ = ('name', 'dist', 'children', 'is_ecor', 'strain', 'clade_size')
    def __init__(self, name, dist, children=None, is_ecor=False, strain='', clade_size=None):
        self.name = name
        self.dist = dist
        self.children = children if children is not None else []
        self.is_ecor = is_ecor
        self.strain = strain
        self.clade_size = clade_size

    def is_leaf(self):
        return len(self.children) == 0

    def traverse(self, order='preorder'):
        """Simple iterative traversal yielding self and descendants."""
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            for c in reversed(node.children):
                stack.append(c)


def build_skeleton(tree, name_idx, ecor_ids, strain_of):
    """Minimal subtree spanning the ECOR leaves; non-ECOR-only clades collapsed
    into leaves named 'clade:N'. Returns SNode root + counts."""
    t0 = time.time()
    counts = Counter()
    for pid in ecor_ids:
        node = name_idx.get(pid)
        if node is None:
            continue
        cur = node
        while cur is not None:
            counts[id(cur)] += 1
            cur = cur.up

    def convert(node):
        kids = []
        for c in node.children:
            if counts.get(id(c), 0) > 0:
                kids.append(convert(c))
            else:
                kids.append(SNode(f'clade:{len(c)}', c.dist, clade_size=len(c)))
        is_ecor = node.is_leaf() and counts.get(id(node), 0) > 0
        strain = strain_of.get(node.name, '') if is_ecor else ''
        return SNode(node.name, node.dist, kids, is_ecor=is_ecor,
                     strain=strain)

    root = convert(tree)
    n_ecor = sum(1 for n in root.traverse() if n.is_ecor)
    n_collapsed = sum(1 for n in root.traverse() if n.clade_size is not None)
    print(f'[skeleton] built in {time.time()-t0:.1f}s: ECOR leaves={n_ecor}, '
          f'collapsed clades={n_collapsed}, nodes={sum(1 for _ in root.traverse())}')
    return root


def plotly_figure_from_snode(root, title, height=2600):
    """Rectangular plotly figure for an SNode tree. Returns dict for fig.update_layout.
    Builds traces: edges (one line trace), ECOR leaves (red markers + text),
    collapsed clades (grey triangles), internal nodes (dots)."""
    import plotly.graph_objects as go

    leaf_order, coords, parent = layout(root)

    # edges
    ex, ey = [], []
    for node in root.traverse('preorder'):
        pid = parent.get(id(node))
        if pid is None:
            continue
        (x1, y1) = coords[pid]
        (x2, y2) = coords[id(node)]
        ex += [x1, x2, None]
        ey += [y1, y2, None]

    ecor_x, ecor_y, ecor_txt = [], [], []
    ecor_hover = []
    coll_x, coll_y, coll_txt = [], [], []
    for node in root.traverse('preorder'):
        x, y = coords[id(node)]
        if node.is_ecor:
            ecor_x.append(x); ecor_y.append(y)
            ecor_txt.append(node.strain)
            ecor_hover.append(f'{node.strain}<br>{node.name}')
        elif node.clade_size is not None:
            coll_x.append(x); coll_y.append(y)
            coll_txt.append(f'clade:{node.clade_size}')

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=ex, y=ey, mode='lines', line=dict(color=SKELETON_EDGE, width=1.0),
        hoverinfo='skip', name='branches'))
    if ecor_x:
        fig.add_trace(go.Scatter(
            x=ecor_x, y=ecor_y, mode='markers+text', text=ecor_txt,
            textposition='middle right', textfont=dict(size=9, color=ECOR_COLOR),
            marker=dict(color=ECOR_COLOR, size=9, symbol='circle'),
            customdata=ecor_hover, hovertemplate='%{customdata}<extra></extra>',
            name='ECOR'))
    if coll_x:
        fig.add_trace(go.Scatter(
            x=coll_x, y=coll_y, mode='markers',
            marker=dict(color='#c9d1d6', size=7, symbol='triangle-up',
                        line=dict(color='#7f8b94', width=0.6)),
            text=coll_txt, hovertemplate='%{text}<extra></extra>', name='collapsed clade'))
    ymax = len(leaf_order)
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        height=height, margin=dict(l=10, r=140, t=60, b=30),
        xaxis=dict(title='UPGMA height (MASH distance)', range=[-0.02, 1.05]),
        yaxis=dict(range=[-1, ymax + 1], showticklabels=False),
        showlegend=True, legend=dict(x=1.02, y=1.0),
        hovermode='closest',
    )
    return fig


def render_skeleton_html(tree, name_idx, ecor_ids, strain_of, out_path):
    root = build_skeleton(tree, name_idx, ecor_ids, strain_of)
    fig = plotly_figure_from_snode(
        root, 'ECOR-highlighted skeleton of the all-prophage MASH tree '
              '(minimal subtree spanning all ECOR leaves; non-ECOR clades collapsed '
              'to triangles)', height=2800)
    n_ecor = sum(1 for n in root.traverse() if n.is_ecor)
    n_collapsed = sum(1 for n in root.traverse() if n.clade_size is not None)
    fig.write_html(out_path, include_plotlyjs='cdn', full_html=True,
                   config={'displaylogo': False})
    print(f'[skeleton] wrote {out_path.name} (ECOR={n_ecor}, collapsed={n_collapsed})')


# ------------------------------------------------------- neighborhood explorer
def render_explorer_html(tree, name_idx, ecor_ids, strain_of, idx_rows, out_path):
    """One HTML with a dropdown: per-ECOR local clade (max <=60 leaves) tree."""
    import plotly.graph_objects as go

    ecor_data = {}
    for pid in ecor_ids:
        leaf = name_idx.get(pid)
        if leaf is None:
            continue
        clade = leaf
        while clade.up is not None and len(clade.up) <= NBR_MAX:
            clade = clade.up
        leaf_order, coords, parent = layout(clade)
        members = [lf.name for lf in leaf_order]
        # edge lines
        ex, ey = [], []
        for node in clade.traverse('preorder'):
            p = parent.get(id(node))
            if p is None:
                continue
            ex += [coords[p][0], coords[id(node)][0], None]
            ey += [coords[p][1], coords[id(node)][1], None]
        # tip coordinates in leaf order
        lx = [coords[id(nd)][0] for nd in leaf_order]
        ly = [coords[id(nd)][1] for nd in leaf_order]
        row = idx_rows.get(pid, {})
        ecor_data[pid] = {
            'name': pid,
            'strain': strain_of.get(pid, ''),
            'edges': [ex, ey],
            'members': members,
            'lx': lx,
            'ly': ly,
            'n': len(members),
            'taxonomy': row.get('taxonomy', ''),
            'contig': row.get('source_contig', ''),
            'coords': (row.get('contig_start', ''), row.get('contig_end', '')),
        }

    json_blob = json.dumps(ecor_data)
    print(f'[explorer] embedded {len(ecor_data)} neighborhoods '
          f'({len(json_blob)/1e6:.1f} MB JSON)')

    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>ECOR neighbourhood explorer</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<style>
 body {{ font-family: -apple-system, 'Segoe UI', Roboto, sans-serif; margin: 20px; }}
 h1 {{ font-size: 18px; }}
 select {{ font-size: 14px; padding: 5px; width: 440px; }}
 .info {{ font-size: 13px; color: #444; margin: 8px 0; }}
 #tree {{ width: 100%; height: 560px; }}
 #members {{ font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 12px;
            border: 1px solid #ddd; padding: 8px; margin-top: 8px; height: 180px;
            overflow-y: scroll; white-space: pre; }}
</style></head><body>
<h1>ECOR neighbourhood explorer — local clade (≤{NBR_MAX} leaves) of each ECOR prophage in the MASH tree</h1>
<div class="info">Red = ECOR element · grey = closest relatives in the all-prophage UPGMA tree · zoom/pan freely.</div>
<select id="sel"></select>
<div id="info" class="info"></div>
<div id="tree"></div>
<div id="members"></div>
<script>
const DATA = {json_blob};
const sel = document.getElementById('sel');
const keys = Object.keys(DATA).sort();
for (const k of keys) {{
  const o = document.createElement('option');
  o.value = k; o.textContent = `${{DATA[k].strain || '?'}} — ${{k}}`;
  sel.appendChild(o);
}}
function render() {{
  const pid = sel.value; const d = DATA[pid];
  const traceEdge = {{ x: d.edges[0], y: d.edges[1], mode: 'lines',
    line: {{ color: '#5b6770', width: 1.2 }}, hoverinfo: 'skip' }};
  const n = d.members.length;
  const fig = {{
    data: [
      traceEdge,
      {{ x: d.lx, y: d.ly, mode: 'markers+text', text: d.members,
        textposition: 'middle right', textfont: {{ size: 10,
          color: d.members.map(m => m === pid ? '#d62728' : '#555') }},
        marker: {{ size: d.members.map(m => m === pid ? 11 : 5),
          color: d.members.map(m => m === pid ? '#d62728' : '#9aa4ad') }},
        customdata: d.members,
        hovertemplate: '%{{customdata}}<extra></extra>' }}
    ],
    layout: {{
      title: {{ text: `${{d.strain}} — ${{pid}} — neighbourhood of ${{n}} leaves`,
               font: {{ size: 13 }} }},
      height: 520, margin: {{ l: 10, r: 200, t: 50, b: 30 }},
      xaxis: {{ title: 'UPGMA height (MASH distance)', range: [-0.03, 1.08] }},
      yaxis: {{ showticklabels: false, range: [-1, n + 1] }},
      showlegend: false, hovermode: 'closest'
    }}
  }};
  Plotly.react('tree', fig.data, fig.layout);
  document.getElementById('info').innerHTML =
    `${{n}} leaves · ECOR: ${{pid}} · taxonomy: ${{d.taxonomy}} · contig: ${{d.contig}} (${{d.coords[0]}}–${{d.coords[1]}})`;
  document.getElementById('members').textContent = d.members.join('\\n');
}}
sel.addEventListener('change', render);
render();
</script>
</body></html>"""
    with open(out_path, 'w') as f:
        f.write(page)
    print(f'[explorer] wrote {out_path.name}')


# ------------------------------------------------------- neighbourhood PNGs
def render_neighborhood_pngs(tree, name_idx, ecor_ids, strain_of, idx_rows, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    done = 0
    for pid in ecor_ids:
        leaf = name_idx.get(pid)
        if leaf is None:
            continue
        clade = leaf
        while clade.up is not None and len(clade.up) <= NBR_MAX:
            clade = clade.up
        leaf_order, coords, parent = layout(clade)
        n = len(leaf_order)
        segs = []
        for node in clade.traverse('preorder'):
            p = parent.get(id(node))
            if p is None:
                continue
            segs.append((coords[p], coords[id(node)]))
        fig_h = max(3.0, 0.28 * n)
        fig, ax = plt.subplots(figsize=(10, fig_h), dpi=100)
        ax.add_collection(LineCollection(segs, colors=NON_ECOR_COLOR, linewidths=1.0))
        # labels
        for nd in leaf_order:
            x, y = coords[id(nd)]
            color = ECOR_COLOR if nd.name == pid else '#444'
            fw = 'bold' if nd.name == pid else 'normal'
            ax.text(x, y, f'  {nd.name}', fontsize=8, color=color,
                    fontweight=fw, va='center', ha='left')
        # ECOR edge + tip
        e = name_idx[pid]
        ep = parent.get(id(e))
        if ep is not None:
            ax.plot([coords[ep][0], coords[id(e)][0]],
                    [coords[ep][1], coords[id(e)][1]],
                    color=ECOR_COLOR, lw=2.0)
        ax.scatter([coords[id(e)][0]], [coords[id(e)][1]], s=28, c=ECOR_COLOR, zorder=5)
        ax.set_xlim(-0.05, 1.1)
        ax.set_ylim(-1, n + 0.5)
        ax.set_yticks([])
        row = idx_rows.get(pid, {})
        ax.set_title(f'{pid}  ·  {strain_of.get(pid, "")}  ·  '
                     f'{row.get("source_contig", "")}:{row.get("contig_start", "")}-'
                     f'{row.get("contig_end", "")}  ({row.get("length_bp", "")} bp)  '
                     f'·  {n} leaves in clade', fontsize=9, loc='left')
        fig.tight_layout()
        fig.savefig(out_dir / f'{pid}.png', bbox_inches='tight', pad_inches=0.15)
        plt.close(fig)
        # palette-quantize + optimize to keep the directory git-friendly
        from PIL import Image as PILImage
        p = out_dir / f'{pid}.png'
        im = PILImage.open(p).convert('RGB')
        q = im.quantize(colors=256, method=PILImage.Quantize.FASTOCTREE)
        q.save(p, optimize=True)
        done += 1
    print(f'[neighborhoods] wrote {done} PNGs to {out_dir} in {time.time()-t0:.1f}s')


# ------------------------------------------------------------------ main
def main():
    args = sys.argv[1:]
    do_full = '--full' in args or '--all' in args or len(args) == 0
    do_skel = '--skeleton' in args or '--all' in args or len(args) == 0
    do_expl = '--explorer' in args or '--all' in args or len(args) == 0
    do_nbr = '--neighborhoods' in args or '--all' in args or len(args) == 0

    print('[load] tree ...')
    t0 = time.time()
    tree = ete3.Tree(str(TREE_NWK), format=1)
    name_idx = build_name_index(tree)
    tags = load_tags()
    idx_rows = load_index()
    ecor_ids = sorted(tags.keys())
    print(f'       tree={len(tree)} leaves, {len(ecor_ids)} ECOR ids, {time.time()-t0:.1f}s')

    if do_full:
        render_full_tree_png(tree, name_idx, ecor_ids, FULL_PNG)
    if do_skel:
        render_skeleton_html(tree, name_idx, ecor_ids, tags, SKELETON_HTML)
    if do_expl:
        render_explorer_html(tree, name_idx, ecor_ids, tags, idx_rows, EXPLORER_HTML)
    if do_nbr:
        render_neighborhood_pngs(tree, name_idx, ecor_ids, tags, idx_rows, NBR_DIR)


if __name__ == '__main__':
    main()
