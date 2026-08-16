# NTM v1 → v2 — what the redo changed

Generated: 2026-08-16 (task ntm-v2-release). Compares the v1 release
(`ntm/RELEASE.md`, geNomad calls, 913 genomes) against v2 (collaborator
QC-passed cohort, 2,388 genomes). Numbers below come from
`v1_v2_summary.json` / `v1_v2_clade_map.tsv` / `v1_v2_prophage_stats.tsv`
(script `ntm/v2/scripts/compare_v1_v2.py`) and the like-for-like QC re-run of
v1 (`v1_qc_baseline/qc_table.tsv`).

## Headline

| metric | v1 | v2 | change |
|---|--:|--:|---|
| source cohort | 7,352 accessions (NCBI family query, TB-complex excluded) → 7,303 genomes | 13,122 canonical genomes (collaborator QC-passed list) + 9,543 run-assembly ids pending | +79% genomes now, more to come |
| prophage caller | geNomad v1.12 end-to-end | collaborator manifest calls (shipped per accession) | caller swap |
| prophages | 10,438 (from 4,021 genomes) | 8,502 used (canonical cohort; 23,212 deduped calls total, 14,741 on pending runs) | −18% used |
| prophage clades (thr 0.25/max 100) | 913 (472 alignable + 441 singletons) | 2,388 (1,251 + 1,137) | ×2.6 |
| ML genomes | 913 | 2,388 | ×2.6 |
| ancestral genomes | 472 | 1,251 | ×2.7 |
| host clades (re-sketch, dist ≤0.05) | 342 represented (of 730) | 623 represented (of 1,044) | +82% |
| host species labelled | NTM only (TB-complex excluded by design) | TB-complex labels kept verbatim (ANI: *M. decipiens*), decision pending | scope change |

## Prophage caller agreement (collaborator calls vs geNomad)

On the **2,867 genomes** called by both (numeric accession match; geNomad
prophages there: 8,351, collaborator calls: 7,704), matched by same contig +
reciprocal interval overlap ≥50%:

| | prophages | % |
|---|--:|--:|
| matched (≥50% reciprocal overlap) | **3,448** | 41% of v1 / 45% of v2 |
| geNomad-only (dropped by collaborator calls) | 4,903 | 59% of v1 |
| collaborator-only (added vs geNomad) | 4,260 | 55% of v2 |

Stringency sensitivity: reciprocal ≥30% → 4,273 matched; **any overlap ≥1 bp →
4,697 (56% of v1 geNomad prophages)**. So ~44% of geNomad calls have **no
overlapping collaborator call at all** on shared genomes — this is genuine
caller disagreement (different tools/thresholds), not boundary jitter. Both
directions disagree at similar rates: the collaborator calls are not a
superset or subset of geNomad; they are a different view. Per-accession detail:
`v1_v2_prophage_stats.tsv`.

Cohort effects add to this: v1's cohort (7,303 genomes, GCF-heavy: 5,303 GCF /
2,000 GCA, TB-complex excluded) and v2's canonical cohort (13,122, GCA-preferred,
TB labels kept) share **1,150 exact accessions but 5,407 assemblies by numeric**
(assembly twins resolved to different copies), so most shared assemblies enter v2
with different GCA/GCF sequence copies — and 7,715 v2 assemblies were never in
v1 at all.

## Clade continuity (member overlap)

Each v1 clade was mapped through matched prophages to the v2 clade(s)
containing its members (`v1_v2_clade_map.tsv`):

| v1 clades (913) | count | % |
|---|--:|--:|
| majority-matched (≥50% of matched members in one v2 clade, ≥2 members) | **104** | 11% |
| any member matched (≥1) | 286 | 31% |
| of which ≥50% modal incl. single-member clades | 223 | 24% |
| split across multiple v2 clades | 128 | 14% |
| vanished (no member has a v2 counterpart) | **627** | 69% |

Only 272 v2 clades (11%) contain any matched v1 member — the v2 clade tree is
mostly new: a consequence of (a) the caller disagreement above, (b) the larger
canonical cohort (2,634 abscessus-dominated redundant assemblies → 736 clades
with internal mash distance 0.0), and (c) fresh clade ids (v2 clades are
recomputed, not updated in place). **Do not treat v2 clade ids as v1
continuations**; use `v1_v2_clade_map.tsv` to bridge where a bridge exists.

## Host-clade coverage

- v1: 342 host clades represented (of 875 host clades over 7,303 genomes;
  TB-complex excluded).
- v2: **623** host clades represented (of 1,044 over 13,122 genomes).
- Species emphasis shifted with the cohort: abscessus complex grew from 212 →
  **1,182** ML genomes (host clade 0002; the collaborator list is heavily
  abscessus-sampled), avium/intracellulare remain large, and TB-labelled
  genomes now contribute 23 ML genomes (host_clade_0001) that v1 excluded by
  design (kept verbatim pending user decision — see RELEASE.md cohort note).
- v1 host-clade ids are a different id space (v2 re-sketch); species-level
  correspondence via `per_ntm_clade_summary.tsv` (v2) vs v1's.

## ML genome length / GC / N-run distributions (same QC, fixed metric)

| statistic | v1 (913) | v2 (2,388) | v2 alignable only (1,251) |
|---|--:|--:|--:|
| length min / median / mean / max (kb) | 4.3 / **33.5** / 53.8 / 153 | 0.2 / **20.3** / 28.3 / 150 | 0.2 / 30.3 / 37.0 / 150 |
| <10kb | 73 | 707 | 240 |
| 10–50kb | 483 | 1,311 | 676 |
| 50–100kb | 177 | 330 | 295 |
| 100–150kb | 122 | 39 | 39 |
| ≥150kb | 58 | 1 | 1 |
| GC median | **63.98%** | **63.81%** | 63.72% |
| N-run: median / #≥100 / max | 0 / 10 / 411 | 0 / 6 / 5,507 | 0 / 0 / 0 |
| member recovery <0.5 (mosaic metric) | 461/913 (median 0.25) | 1,014/1,251 alignable (median 0.0) | same |

Interpretation:

- **GC is stable at ~64%** in both releases — the mycobacterial-host
  expectation holds after the caller and cohort swap (v1 63.98%, v2 63.81%).
- **v2 is shorter and more small-genome-heavy**: 707 genomes <10 kb (mostly
  short singleton collaborator calls pass through unchanged), and the ≥150 kb
  tail collapsed (58 → 1) because the 150-kb traversal cap now binds rarely
  (v1's large 100-member k-nearest mosaics frequently hit it).
- **N-runs improved** on alignable clades (0 genomes with N-run ≥100 in v2
  alignable vs v1's 10), but one v2 singleton carries a 5,507-N run
  (verbatim member sequence).
- **Member recovery is partial in both versions** (mosaic property of
  `traverse_partitions.py` output over redundant k-nearest partition sets).
  v1 was never measured before because the QC script's completeness/identity
  math was inverted (fixed in this task); v1's baseline in `v1_qc_baseline/`
  is the first honest look, and it shows the same phenomenon (median 0.25).
  v2's alignable median is lower (0.0) mainly because 100-member clades are
  dominated by near-identical redundant members (median internal mash 0.006),
  so each member is a smaller fraction of the stitched mosaic.

## What to carry forward

1. v2 supersedes v1 for the canonical cohort; keep v1 for the ~1,050
   accessions absent from the collaborator list.
2. When run-assembly FASTAs arrive (5,907 prophage-bearing runs), re-run
   clades → ML → release; expect the biggest movement in the abscessus
   complex and the singleton tail.
3. The user decision on TB-labelled genomes (13,227 accessions, ANI
   *M. decipiens*) determines whether host_clade_0001's 23 ML genomes stay.
4. The prophage-caller disagreement (44% of geNomad calls unmatched by any
   collaborator call) is large enough that a third-caller arbitration
   (PhiSpy/Phigaro, as already recommended in v1's caveats for *M. ulcerans*)
   would be worthwhile before any biological conclusion rests on either
   caller alone.
