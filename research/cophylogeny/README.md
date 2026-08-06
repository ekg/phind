# Prophage vs host cophylogeny

Compares the prophage phylogeny (`research/mash_tree/full_prophages_tree.nwk`,
132,393 leaves) to the host *E. coli* phylogeny
(`research/host_mash_tree/host_tree.nwk`, 26,074 leaves). Each prophage leaf
maps to its host by accession prefix (`GCF_xxx_prophage_N` → `GCF_xxx`), a
many-to-one mapping (132,393 prophages → 26,074 hosts; median 5 prophages/host,
max 25).

Implemented with R `ape` + `phytools` (cophylo) + `paco` (PACo) +
`ape::parafit`, delivered at three resolutions.

## Mapping
`mapping/`
- `prophage_host_map.tsv` — full 132,393-row prophage → host map (covers 100%).
- `mapping_stats.json` — coverage / multiplicity statistics.

## Resolution 1 — ECOR tanglegram
`ecor_resolution/`
- `ecor_tanglegram.png` — host tree (72 ECOR strains, colored by phylogroup)
  vs prophage tree (300 ECOR prophages) with connecting links
  (`phytools::cophylo`).
- `ecor_tanglegram_interactive.html` — self-contained interactive SVG tanglegram
  (hover a host strain to highlight its prophage lineage links).
- `ecor_tanglegram_assoc.tsv` — prophage → host association used by cophylo.

## Resolution 2 — Clade resolution
`clade_resolution/`
- `association_matrix.tsv` — bipartite **prophage_clade (585) × host_phylogroup (10)
  counts**.
- `clade_memberships.tsv` — prophage → tight clade (+ community, host, phylogroup).
- `clade_meta.tsv`, `clade_host_set.tsv` — clade metadata and host sets.
- `host_phylogroup_map.tsv` — host accession → phylogroup.
- `host_clade_matrix.tsv` — host (3,766) × prophage-clade (585) 0/1 presence.
- `community_phylogroup.tsv` — 12-community × 10-phylogroup aggregation.
- `compact_tanglegram.png` — compact community × phylogroup tanglegram.

## Resolution 3 — Full-26k quantitative
`full_26k/`
- `cophylogeny_stats.json` — Mantel (host mash distance vs prophage-clade-set
  Jaccard), PACo, and ParaFit results.
- `matched_subset.tsv` — hosts used for the PACo / ParaFit subset.
- `full_26k_stats.R` — reproducible analysis script.

## Interpretation
See `REPORT.md`.
