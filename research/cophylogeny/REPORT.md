# Prophage vs host cophylogeny in *E. coli*: co-diversification or horizontal transfer?

## 1. Scope and data

| Component | File | Size |
|---|---|---|
| Prophage phylogeny | `research/mash_tree/full_prophages_tree.nwk` | 132,393 tips |
| Host (*E. coli*) phylogeny | `research/host_mash_tree/host_tree.nwk` | 26,074 tips |
| Prophage → host map | `research/cophylogeny/mapping/prophage_host_map.tsv` | **100% coverage** (132,393/132,393) |
| Tight prophage clades | `research/clades/*/tight_clades.json` | 585 clades / 12 communities, 19,638 prophages |
| Host phylogroups | `research/phylogroups/phylogroups.tsv` | A/B1/B2/C/D/E/F/G + clade/other |

Each prophage leaf `GCF_xxx_prophage_N` maps to host `GCF_xxx` by accession
prefix. 132,393 prophages → 26,074 hosts (many-to-one; median 5 prophages/host,
max 25, min 1). **The mapping covers 100% of prophage leaves and every mapped
host exists in the host tree (0 missing).**

Analyses use R `ape` + `phytools` (`cophylo`) + `paco` (PACo) +
`ape::parafit`.

## 2. Three resolutions

### 2.1 ECOR-resolution tanglegram
Prune the host tree to the 72 ECOR reference strains and the prophage tree to
the 300 prophages resident in those strains (research/ecor/). `phytools::cophylo`
aligns both phylogenies and draws the 300 many-to-one association links.
Host tips are labeled/colored by curated ECOR phylogroup (A, B1, B2, C, D, E,
F, G).
- `ecor_resolution/ecor_tanglegram.png`
- `ecor_resolution/ecor_tanglegram_interactive.html` (hover a host → highlights
  its resident prophage lineage links)

### 2.2 Clade-resolution
Collapse 132,393 prophages → 585 tight clades (within 12 communities
`research/clades/`), and hosts → 10 phylogroup groups. Produces the bipartite
**association_matrix.tsv (585 prophage clades × 10 host phylogroups; counts)**,
a host × clade presence matrix, and a compact community × phylogroup
tanglegram.
- `clade_resolution/association_matrix.tsv`
- `clade_resolution/compact_tanglegram.png`
- `clade_resolution/community_phylogroup.tsv`

### 2.3 Full-26k quantitative
Across the 3,766 hosts that carry ≥1 prophage assigned to a tight clade:
host-pair mash (patristic) distance is compared to the Jaccard dissimilarity
of their resident prophage-clade sets (Mantel), and PACo + ParaFit are run on
a matched subset of 200 hosts.
- `full_26k/cophylogeny_stats.json`
- `full_26k/matched_subset.tsv`

## 3. Results

### 3.1 Mantel test (full matched set, 3,766 hosts; 999 perms, Pearson)
**r = 0.0279, p = 0.001**

Host patristic (mash) distance vs prophage-clade-set Jaccard dissimilarity: a
very small but statistically significant positive correlation. Hosts that are
more distant in the mash tree share on average only slightly fewer of the same
prophage clades.

### 3.2 PACo (matched subset, 200 hosts × 377 clades; cailliez, r0, 499 perms)
**Procrustes sum-of-squared-residuals = 0.6851, p ≈ 0 (< 0.002)**

Two phylogenies are significantly more congruent than chance given the
observed host–prophage-clade network, but the residual is high
(ss ≈ 0.69 → only ~31% of the variance is reconciled by the Procrustes
superimposition), i.e. the congruence is real but weak.

### 3.3 ParaFit (matched subset, 200 hosts × 377 clades; D1 host, D2 clade, HP01 host × clade)
**ParaFitGlobal = 0.9101, p = 0.001 (nperm = 999)**

The global cophylogenetic fit is statistically significant, again indicating
a non-random association between prophage-clade and host phylogenies.

## 4. Interpretation: co-diversification vs horizontal transfer

All three tests reach statistical significance (Mantel p=0.001, PACo p<0.002,
ParaFit p=0.001) — there is real, non-random phylogenetic structure between the
prophage-clade and host side. Yet the *magnitude* of the effect is uniformly
weak:
- Mantel r ≈ **0.028** means host mash distance explains essentially none of the
  variance in prophage-clade-set composition (~0.08% of variance).
- PACo Procrustes ss ≈ **0.69** leaves most of the variance unexplained by
  superimposition.
- The ECOR and clade-resolution tanglegrams (Section 2) show many resident
  prophage lineages scattered across distantly related host strains, and most
  phylogroup/community bins contain a broad mix of prophage clades.

**Conclusion: horizontal transfer dominates.** The overwhelming pattern is one
of prophages moving promiscuously across the *E. coli* pangenome — closely
related prophages are found in phylogenetically distant hosts, and a given host
lineage carries a heterogeneous set of prophage clades. The small but
significant Mantel/PACo/ParaFit signals are consistent with a modest
co-phylogenetic component: some prophage lineages do track their host clade
(vertically maintained / host-adapted), contributing a detectable, but minor,
cophylogeny signal on top of a largely horizontally transferred background. This
fits the expectation for temperate phages in *E. coli*: integration into a
particular chromosomal site, plus occasional lineage-specific enrichment, gives
weak-to-moderate codivergence that is easily swamped by lateral exchange.

## 5. Caveats
- The tight-clade definition covers 19,638 of 132,393 prophages (the subset
  used for all-wave alignment); the remaining prophages lack clade assignment
  and are therefore excluded from the clade-level and quantitative analyses.
- Host phylogroup is the grouping used for the host side of the association
  matrix; finer host mash-clades could be substituted.
- Mash trees are approximate (Mash sketches), not full core-genome alignments.
