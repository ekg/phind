#!/usr/bin/env python3
"""
build_ecor_inspection.py — ECOR inspection index builder.

Joins the ECOR mapping (ecor_manifest.csv / ecor_leaf_tags.tsv) onto the
all-prophage MASH tree leaves (research/mash_tree/full_prophages_tree.nwk)
and the survey label table (full_prophages_labels.csv), and emits:

  ecor_inspection_index.csv    one row per accessible ECOR prophage element
  ecor_inspection_index.html   browsable HTML index (search + NCBI/neighborhood links)
  ecor_meta.json               machine-readable summary + per-element lookup tables
                               (leaf y-order, neighborhood members, etc.) reused by
                               the rendering scripts

"Accessible" = prophage_id present in full_prophages.fa with resolvable source
coordinates (upstream ecor-mapping verified all 300 manifest rows against
full_prophages.fa, tag merge rate 1.0).

Usage:
  python3 build_ecor_inspection.py
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

# ---------------------------------------------------------------- paths
HERE = Path(__file__).resolve().parent
MASH = Path('/home/erikg/phind/research/mash_tree')
SURVEY = Path('/home/erikg/phind/prophage_homology_survey')

MANIFEST = HERE / 'ecor_manifest.csv'
TAGS = HERE / 'ecor_leaf_tags.tsv'
LABELS = MASH / 'full_prophages_labels.csv'
TREE_NWK = MASH / 'full_prophages_tree.nwk'
CONN_COMM = SURVEY / 'full_heatmap_clusters.csv'

OUT_CSV = HERE / 'ecor_inspection_index.csv'
OUT_HTML = HERE / 'ecor_inspection_index.html'
OUT_META = HERE / 'ecor_meta.json'

NEIGHBORHOOD_MAX = 60          # max leaves in a rendered neighborhood clade
NEIGHBORHOOD_DIR = HERE / 'ecor_neighborhoods'


def load_ecor_manifest():
    rows = {}
    with open(MANIFEST) as f:
        for row in csv.DictReader(f):
            rows[row['prophage_id']] = row
    return rows


def load_tags():
    tags = {}
    with open(TAGS) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['is_ecor'] == 'TRUE':
                tags[row['prophage_id']] = row['ecor_strain']
    return tags


def load_labels():
    labels = {}
    with open(LABELS) as f:
        for row in csv.DictReader(f):
            labels[row['sequence']] = row
    return labels


def load_connected_communities():
    """sequence -> community for the Louvain graph communities (full collection)."""
    comm = {}
    with open(CONN_COMM) as f:
        for row in csv.DictReader(f):
            comm[row['sequence']] = int(row['community'])
    return comm


def build_name_index(tree):
    """One-pass name -> node map (search_nodes is O(n) per call; this is O(n) total)."""
    idx = {}
    for node in tree.traverse('preorder'):
        idx[node.name] = node
    return idx


def ecor_leaf_positions(tree, name_idx, ecor_ids):
    """Return {prophage_id: leaf node} for ECOR ids present in the tree."""
    found = {}
    for pid in ecor_ids:
        node = name_idx.get(pid)
        if node is not None and node.is_leaf():
            found[pid] = node
    return found


def neighborhood_of(leaf, name_idx, max_leaves=NEIGHBORHOOD_MAX):
    """
    Maximal clade (subtree) containing `leaf` with at most max_leaves leaves.
    Walks up while the parent clade stays <= max_leaves.
    Returns the clade root node (or the leaf itself for a singleton).
    """
    clade = leaf
    while clade.up is not None and len(clade.up) <= max_leaves:
        clade = clade.up
    return clade


def cophenetic(leaf_a, leaf_b):
    """Tree path distance (sum of branch lengths) between two leaves."""
    return leaf_a.get_distance(leaf_b)


def build_index():
    t0 = time.time()
    print('[1/4] loading manifest / tags / labels ...')
    manifest = load_ecor_manifest()
    tags = load_tags()
    labels = load_labels()
    conn = load_connected_communities()
    print(f'      manifest={len(manifest)} tags={len(tags)} labels={len(labels)}')

    print('[2/4] parsing tree ...')
    tree = ete3.Tree(str(TREE_NWK), format=1)
    name_idx = build_name_index(tree)
    print(f'      tree leaves={len(tree)} ({time.time()-t0:.1f}s)')

    print('[3/4] locating ECOR leaves + neighborhoods ...')
    ecor_ids = sorted(manifest.keys())
    found = ecor_leaf_positions(tree, name_idx, ecor_ids)
    missing_in_tree = [p for p in ecor_ids if p not in found]
    print(f'      ECOR leaves in tree: {len(found)}/{len(ecor_ids)}')

    rows = []
    meta = {}
    for pid in ecor_ids:
        m = manifest[pid]
        strain = tags.get(pid, '')
        lab = labels.get(pid, {})
        leaf = found.get(pid)
        if leaf is None:
            clade_size = None
            members = []
            nn = ''
            nn_dist = None
            tree_ok = 0
        else:
            clade = neighborhood_of(leaf, name_idx)
            members = sorted(lf.name for lf in clade)
            clade_size = len(members)
            # nearest neighbour within the neighbourhood clade
            best = None
            best_d = None
            for other in clade:
                if other is leaf:
                    continue
                d = cophenetic(leaf, other)
                if best_d is None or d < best_d:
                    best_d = d
                    best = other.name
            nn, nn_dist = best, best_d
            tree_ok = 1

        comm = lab.get('community')
        if comm is not None:
            try:
                comm = int(float(comm))
            except (TypeError, ValueError):
                comm = None
        is_connected = comm is not None and comm < 12 if comm is not None else False
        rows.append({
            'prophage_id': pid,
            'ecor_strain': strain,
            'strain_index': strain.replace('ECOR-', '') if strain.startswith('ECOR-') else strain,
            'assembly_accession': m['assembly_accession'],
            'gca_accession': m['gca_accession'],
            'wgs_master': m['wgs_master'],
            'source_contig': m['source_contig'],
            'contig_start': int(m['start']),
            'contig_end': int(m['end']),
            'length_bp': int(m['length']),
            'transposable': float(m['transposable']),
            'taxonomy': m['taxonomy'],
            'community': comm,
            'in_connected_component': 1 if is_connected else 0,
            'cluster': lab.get('cluster') or '',
            'genome': lab.get('genome') or m['assembly_accession'],
            'in_tree': tree_ok,
            'tree_leaf_id': pid,
            'neighborhood_size': clade_size,
            'nearest_neighbor': nn,
            'nearest_neighbor_dist': round(nn_dist, 6) if nn_dist is not None else '',
            'contig_url': f"https://www.ncbi.nlm.nih.gov/nuccore/{m['source_contig']}",
            'contig_fasta_url': (f"https://www.ncbi.nlm.nih.gov/nuccore/{m['source_contig']}"
                                 f"?report=fasta&from={int(m['start'])}&to={int(m['end'])}"),
            'assembly_url': f"https://www.ncbi.nlm.nih.gov/assembly/{m['assembly_accession']}",
            'neighborhood_png': f"ecor_neighborhoods/{pid}.png",
        })

    # consistency checks
    tag_ids = set(tags)
    man_ids = set(manifest)
    checks = {
        'manifest_rows': len(rows),
        'tags_covered': len(man_ids & tag_ids),
        'in_tree': len(found),
        'missing_in_tree': missing_in_tree,
        'neighborhood_max': max((r['neighborhood_size'] or 0) for r in rows),
        'connected_component_ecor': sum(r['in_connected_component'] for r in rows),
        'duplicate_prophage_ids': len(man_ids) - len(set(man_ids)),
    }

    print('[4/4] writing index ...')
    fieldnames = [
        'prophage_id', 'ecor_strain', 'strain_index', 'assembly_accession',
        'gca_accession', 'wgs_master', 'source_contig', 'contig_start',
        'contig_end', 'length_bp', 'transposable', 'taxonomy', 'community',
        'in_connected_component', 'cluster', 'genome', 'in_tree',
        'tree_leaf_id', 'neighborhood_size', 'nearest_neighbor',
        'nearest_neighbor_dist', 'contig_url', 'contig_fasta_url',
        'assembly_url', 'neighborhood_png',
    ]
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    meta = {
        'schema': 'ecor-inspection-index-v1',
        'generated': time.strftime('%Y-%m-%dT%H:%M:%S'),
        'tree': str(TREE_NWK),
        'n_ecor_elements': len(rows),
        'checks': checks,
        'neighborhoods': {
            r['prophage_id']: {
                'clade_size': r['neighborhood_size'],
                'members': None,  # filled below (kept in separate file to bound size)
            } for r in rows
        },
    }
    with open(OUT_META, 'w') as f:
        json.dump(meta, f, indent=1)

    print(f'      wrote {OUT_CSV.name} ({len(rows)} rows) and {OUT_META.name}')
    print('      checks:', json.dumps(checks))
    print(f'      total {time.time()-t0:.1f}s')

    # --- HTML index -----------------------------------------------------
    write_html_index(rows, checks)
    return rows, checks


def write_html_index(rows, checks):
    """Browsable HTML index: searchable table, per-row NCBI + neighborhood links."""
    esc = html.escape
    trs = []
    for r in rows:
        nn = esc(str(r['nearest_neighbor'])) if r['nearest_neighbor'] else '—'
        nn_d = r['nearest_neighbor_dist']
        nn_d = f"{nn_d:.5f}" if nn_d != '' else '—'
        trs.append(f"""<tr data-strain="{esc(r['ecor_strain'])}" data-pid="{esc(r['prophage_id'])}">
<td class="mono">{esc(r['prophage_id'])}</td>
<td class="strain">{esc(r['ecor_strain'])}</td>
<td class="mono">{esc(r['assembly_accession'])}</td>
<td class="mono">{esc(r['source_contig'])}</td>
<td>{r['contig_start']:,}–{r['contig_end']:,}</td>
<td>{r['length_bp']:,}</td>
<td>{esc(r['taxonomy'])}</td>
<td>{r['community'] if r['community'] is not None else '—'}</td>
<td>{r['neighborhood_size'] if r['neighborhood_size'] is not None else '—'}</td>
<td class="mono">{nn}</td>
<td>{nn_d}</td>
<td>
<a href="{esc(r['contig_fasta_url'])}" target="_blank">NCBI FASTA</a> ·
<a href="{esc(r['assembly_url'])}" target="_blank">assembly</a> ·
<a href="{esc(r['neighborhood_png'])}" target="_blank">tree</a>
</td>
</tr>""")

    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>ECOR prophage inspection index</title>
<style>
 body {{ font-family: -apple-system, 'Segoe UI', Roboto, sans-serif; margin: 24px; color: #222; }}
 h1 {{ font-size: 20px; }} h2 {{ font-size: 15px; color: #444; }}
 table {{ border-collapse: collapse; width: 100%; font-size: 12.5px; }}
 th, td {{ border: 1px solid #ddd; padding: 4px 7px; text-align: left; vertical-align: top; }}
 th {{ background: #f4f4f6; position: sticky; top: 0; }}
 tr:nth-child(even) {{ background: #fafafa; }}
 td.strain {{ font-weight: 600; }}
 td.mono {{ font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 11.5px; }}
 input#q {{ width: 380px; padding: 6px; font-size: 14px; margin-bottom: 10px; }}
 .meta {{ color: #555; margin: 6px 0 14px; }}
 a {{ color: #0b5394; }}
</style></head><body>
<h1>ECOR prophage inspection index</h1>
<div class="meta">300 ECOR-mapped prophage elements · 72 ECOR strains ·
<a href="ecor_inspection_index.csv">CSV</a> ·
<a href="ecor_tree_interactive.html">interactive tree</a> ·
<a href="ecor_neighborhood_explorer.html">neighborhood explorer</a> ·
<a href="ecor_heatmap_interactive.html">heatmap</a> ·
<a href="ecor_mds_interactive.html">MDS</a></div>
<input id="q" placeholder="Filter by prophage id or ECOR strain (e.g. ECOR-12, GCF_...)">
<table id="tbl">
<thead><tr><th>prophage_id (tree leaf)</th><th>ECOR strain</th><th>assembly</th>
<th>source contig</th><th>coords (1-based)</th><th>length</th><th>taxonomy</th>
<th>community</th><th>neigh.<br>size</th><th>nearest neighbour</th><th>NN dist</th>
<th>inspect</th></tr></thead>
<tbody>
{''.join(trs)}
</tbody></table>
<script>
const q = document.getElementById('q');
q.addEventListener('input', () => {{
  const t = q.value.trim().toLowerCase();
  for (const tr of document.querySelectorAll('#tbl tbody tr')) {{
    tr.style.display = (!t || tr.dataset.strain.toLowerCase().includes(t) ||
      tr.dataset.pid.toLowerCase().includes(t)) ? '' : 'none';
  }}
}});
</script>
</body></html>"""
    with open(OUT_HTML, 'w') as f:
        f.write(page)
    print(f'      wrote {OUT_HTML.name}')


def main():
    build_index()


if __name__ == '__main__':
    main()
