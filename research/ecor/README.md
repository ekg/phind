# ECOR-highlighted inspection (all-prophage MASH tree / heatmap / MDS)

**Task:** `ecor-highlighted-inspection` — merge the ECOR tags onto the
all-prophage MASH tree/heatmap with highlighting, and produce an inspection
index so each accessible ECOR prophage element can be inspected individually.

**Inputs used**

| input | path |
|---|---|
| all-prophage MASH UPGMA tree | `research/mash_tree/full_prophages_tree.nwk` (132,393 leaves) |
| tree label table (community/cluster/genome/in_tree) | `research/mash_tree/full_prophages_labels.csv` |
| ECOR mapping manifest (300 elements, 72 strains) | `research/ecor/ecor_manifest.csv` (from `ecor-mapping`) |
| ECOR leaf tags | `research/ecor/ecor_leaf_tags.tsv` |
| survey heatmap cluster assignments (all 132,393) | `prophage_homology_survey/full_heatmap_clusters.csv` |
| survey MDS (5,000-prophage subset) | `prophage_homology_survey/full_prophage_mds_coords.csv` |
| prophage coordinates (authoritative) | `26k_prophage1.csv` |
| prophage sequences (read-only) | `prophage_homology_survey/full_prophages.fa` (via recovery dir) |

## Key finding that shapes the figures

**All 300 ECOR prophage elements are singleton Louvain communities** in the
MASH graph (verified: 0/300 in the 12 shared connected communities 0–11; each
ECOR element is its own community id 36312–36611). I.e. no ECOR element has a
close relative among the 19,638 prophages of the connected component — the
existing connected-component heatmap (`connected_heatmap_generator.py`)
contains *no* ECOR element and therefore cannot be reused to "highlight ECOR".
The ECOR elements *do* have small local clades (1–60 leaves) in the UPGMA
tree, which is where their relatives live. The figures below therefore
highlight ECOR in tree space (leaf positions, neighbourhoods, cophenetic
distance heatmap/MDS) and in the survey's existing MDS subset.

## Deliverables (all under `research/ecor/`)

### 1. Inspection index
| file | contents |
|---|---|
| `ecor_inspection_index.csv` | **one row per ECOR element (300 rows, 100% of the ecor-mapping manifest)**: `prophage_id`, `ecor_strain`, `assembly_accession`, `wgs_master`, `source_contig`, `contig_start`, `contig_end`, `length_bp`, `transposable`, `taxonomy`, `community`, `in_connected_component`, `cluster`, `genome`, `in_tree`, **`tree_leaf_id`** (cross-link to the MASH tree leaf; = prophage_id since tree leaves are named by prophage id), `neighborhood_size`, `nearest_neighbor` + `nearest_neighbor_dist` (closest relative in the tree), NCBI URLs (contig FASTA with `from/to` coords, assembly) |
| `ecor_inspection_index.html` | browsable version: search box (strain / prophage id), 3 inspect links per row (NCBI FASTA, assembly, neighbourhood tree image) |
| `ecor_meta.json` | machine-readable summary + checks |

Accessibility definition: element present in `full_prophages.fa` with
resolvable source coordinates. Upstream `ecor-mapping` verified all 300
manifest rows against `full_prophages.fa` (tag merge rate 1.0); this task
re-verified all 300 are tree leaves (`in_tree=1`) and spot-checked 3
elements independently (see below).

### 2. Highlighted tree (static + interactive)
| file | contents |
|---|---|
| `ecor_tree_highlighted.png` | full 132,393-leaf UPGMA tree, ECOR leaf edges + tips in red (~15k px tall, zoomable) |
| `ecor_tree_tip_strip.png` | leaf-order positions of the 300 ECOR elements (300 red ticks) |
| `ecor_tree_interactive.html` | **interactive skeleton tree**: the minimal subtree spanning all 300 ECOR leaves (8,045 nodes), with non-ECOR-only clades collapsed to grey triangles (hover: clade size) and ECOR leaves labelled by strain (hover: prophage id) — plotly, zoomable |
| `ecor_neighborhood_explorer.html` | dropdown explorer: for **each** ECOR element, its local clade (max ≤ 60 leaves — its closest relatives in the tree) as an interactive tree; red = the ECOR element, grey = relatives; member list below |
| `ecor_neighborhoods/<prophage_id>.png` | 300 static per-element neighbourhood tree images (linked from the index HTML) |

### 3. Highlighted heatmap / MDS
| file | contents |
|---|---|
| `ecor_heatmap.png` / `ecor_heatmap_interactive.html` | 300 × 300 **cophenetic distance** matrix of all ECOR elements (tree-path distances in the all-prophage UPGMA tree), dendrogram-ordered, strain labels, strain separators — "ECOR in tree space" |
| `ecor_mds.png` / `ecor_mds_interactive.html` | metric MDS (sklearn) of the same 300 × 300 matrix, strain labels |
| `ecor_mds_subset_highlighted.png` / `ecor_mds_subset_interactive.html` | **reuse of the existing survey MDS** (`full_prophage_mds_coords.csv`, 5,000-prophage subset) with the 7 ECOR elements present in that subset marked in red (labels + legend) |

## Reproduction

```bash
cd research/ecor
python3 build_ecor_inspection.py          # index CSV/HTML + ecor_meta.json
python3 render_ecor_tree.py --all         # full-tree PNG, skeleton HTML, explorer HTML, 300 neighbourhood PNGs
python3 render_ecor_heatmap_mds.py        # heatmap + MDS figures
```

Environment notes: `pip install --break-system-packages plotly` (system
python is PEP-668 managed); `ete3`, `pandas`, `numpy`, `scipy`, `sklearn`,
`matplotlib` were already present. ete3 is used for tree parsing only
(its Qt-based `render()` is unavailable headless); all drawing is matplotlib
(static) + plotly (interactive).

## Validation

- **Rendered highlighted tree**: `ecor_tree_highlighted.png` + interactive
  `ecor_tree_interactive.html` (browser-verified: 300 ECOR leaf labels, 3,723
  collapsed clade markers render).
- **Highlighted heatmap/MDS**: `ecor_heatmap.png` (+ interactive),
  `ecor_mds.png` (+ interactive), `ecor_mds_subset_highlighted.png`
  (+ interactive; browser-verified: 7 red ECOR markers in the reused subset).
- **Index coverage**: `ecor_inspection_index.csv` = 300 rows = 100% of the
  300 ECOR-mapped elements in `ecor_manifest.csv`; every row carries
  `tree_leaf_id` (MASH tree leaf id) + `in_tree=1`; every row links to a
  neighbourhood PNG that exists on disk.
- **Spot-check (3 elements, independent of upstream)**: for
  `GCF_003334405.1_prophage_2` (ECOR-1), `GCF_003334385.1_prophage_1`
  (ECOR-2) and `GCF_003333875.1_prophage_3` (ECOR-36):
  - manifest coordinates == `26k_prophage1.csv` begin/end, and
    `length == end − begin + 1` (19682/7511/28209 bp) — all true;
  - sequence extracted from `full_prophages.fa` length == manifest length
    (19682/7511/28209) — all true;
  - all three are leaves of `full_prophages_tree.nwk` — true;
  - neighbourhood PNG exists for each — true.
