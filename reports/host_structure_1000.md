# Frozen host population structure — exact 1,000-assembly pilot

**Release:** `host-structure-1000-v1-3e16e725f70d0fdd`  
**Verdict:** **PASS**  
**External immutable release:** `/home/erikg/phind-data/ecoli26k/v1/releases/run-host-structure-1000/host-structure-1000-v1-3e16e725f70d0fdd`  
**Scope:** host-only, phage-blind; frozen before any downstream phage association.

## Automatic input and safety gates

The run consumed the immutable selection `pilot-cohorts-v1-8afc0ea03d9e50dc`
(`release.json` `d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99`), canonical exact-N=1,000
release `canonical-cohort-1000-v1-4bc3e029e6e0be44` (`release.json`
`14a39b424f2a23de6fa52c173b00e03b167e897baf3a9dbcd9876e31e999740c`), and host consumer certification
`consumer-compatibility-v1-78d7e93f19fa3d87` (`release.json`
`021719ddadd7bb7fa2932d2ef9cb25da9c666ebe0389988691283011ee12f4c7`). Cohort SHA-256 is `265a1e7784a4d5db3ea3577892feba8173290518b6c621f7e5091dbad66bfe77`.
All 1,000 exact assembly revisions, rows, object/index checksums, PanSN names,
contig counts and bases were accounted in immutable cohort order. No genome was
acquired. The global union is the exact frozen 1,000-revision set and never
exceeded 1,000.

Both root inputs matched at start and finish:

- `26k_ecoli_accession.txt`: `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`
- `26k_prophage1.csv`: `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`

The second file was **opaque-hashed only**. Host computation did not parse or
read any prophage presence, count, taxonomy, cluster, coordinates, sequence, or
derived trait. Prophage source/coordinate and extraction semantics are
`NOT_APPLICABLE_HOST_ONLY_ANALYSIS_EXTRACTION_BLOCKED`, never fabricated as PASS. The declared host-only allow-list and
negative unit test enforce this boundary.

## All-host Mash overview

Pinned Mash 2.3 used whole assembly files (never per-contig `-i`) with baseline
`k=21,s=10,000,seed=42`. Exact validation found
**499,500** unordered
off-diagonal pairs, **1,000,000**
directed records, 1,000 zero diagonals, exact triangle/direct agreement and
symmetry. Pinned RapidNJ 2.3.2 retained exactly 1,000 tips. This is an unrooted
whole-genome similarity dendrogram, not a substitution-model phylogeny.

| configuration | sampled Spearman vs baseline | nearest-neighbor agreement | split Jaccard |
|---|---:|---:|---:|
| `k21_s1000_seed42` | 0.990474 | 80.400% | 54.454% |
| `k21_s10000_seed42` | 1.000000 | 100.000% | 100.000% |
| `k21_s50000_seed42` | 0.998803 | 93.000% | 74.606% |
| `k31_s10000_seed42` | 0.996447 | 89.900% | 65.066% |

Six independently seeded `k=21,s=10,000` sketch trees measured sketch-sampling
stability. Of 997 baseline splits,
611 were present under all four parameter settings
and 641 had at least 95% seed support.
Branches below 95% combined support were collapsed, not forced. The recorded
Mash sketch/distance plus RapidNJ command wall time was **171.9 s**;
per-command argv, wall time, `/usr/bin/time -v`, stderr and output digests are
in the external provenance logs.

## Representatives, duplicates, placement and clades

Every eligible host has an exact tip/assembly/biological-unit, sequence class,
near-duplicate class, host-genetic representative, nearest-neighbor tie set,
sampling medoid, placement and membership row. Because no BioSample/isolate
identity table was supplied as a pinned host input, the biological analysis
unit is explicitly the exact assembly revision; downstream code must not
silently merge it.

- exact sequence classes: **1,000**
- exact duplicate non-primary assemblies: **0**
- near-duplicate classes (`Mash D <= 0.0001`): **989**
- near-duplicate non-representatives: **11**
- frozen supported unrooted clades: **7**
- supported fixed memberships: **517**
- deliberately ambiguous memberships: **483**

Only mutually disjoint, 20–400-host unrooted splits present in every Mash
parameter tree and at ≥95% seed support became fixed clades. Other tips say
`AMBIGUOUS_UNROOTED_BACKBONE_OR_UNSUPPORTED`; no association-friendly topology
was chosen. Complete host-only medoid alternatives at k=12, 16 and 20 remain in
the release (3,000 rows).

## High-fidelity host core ensemble

Sixteen host-genetic sampling partitions selected medoids, maximin-diverse and
boundary cases, and fragmentation-QC extremes without phage data. Digest-pinned
minimap2 2.31 built assembly-to-medoid reference-coordinate core alignments. The initial 90%-of-entire-reference
breadth screen correctly failed where lineage-specific accessory sequence was
not shared; those values remain a diagnostic and were not relabeled PASS. The
host-only pilot therefore uses an auditable core-genome denominator: the core
must span at least 50% of the reference, every selected sample must call at
least 95% of that core, and mean core missingness must be at most 5%. Dense
local panels have primary and diverse alternative-reference trees, missingness
checks, a deterministic SNP-density recombination-candidate mask, at least 100
non-recombinant parsimony-informative sites, 100 SNP-site bootstrap replicates,
and 95% support collapse. **6/16** lineages passed every
primary/alternative scientific and reference gate; **10/16** sparse or failing lineages are explicitly blocked from clade inference and all their memberships remain ambiguous. No failed lineage was used to rescue a clade. Among passing lineages, the minimum primary-vs-alternative split concordance was **0.667**. The minimap2 mapping commands took **158.3 s** total wall time; Python alignment/mask/bootstrap work is included in the outer run timing.

The density mask is explicitly a conservative host-core diagnostic, **not** a
claim that Gubbins ran. Unmasked/masked outputs and alternative references make
recombination and reference bias auditable. No independently verified outgroup
exists inside the frozen cohort and out-of-cohort acquisition is prohibited,
so every primary biological tree remains unrooted; midpoint rooting was not
used.

## Resource, restart, determinism and publication

The declared allocation was 68,719,476,736
B RAM, 4,000,000,000,000 B scratch,
30,000,000,000 B durable and
500,000 inodes. Live `findmnt`,
ownership/write probes, bytes/inodes, quotas, and unfinished-write reserve were
recorded before every batch/stage. Peak RSS was
297,144,320 B (0.432%); scratch
upper bound was 10,350,901,329 B and
26,249 files. OOM, cgroup swap growth and system swap
growth were zero; every ≤70% RAM/disk, ≤50% inode and 2x unfinished-write gate
PASS.

A real SIGKILL after five committed materialized views exposed no release.
The same run ID resumed only checksum-valid units and atomically published the
release. `SHA256SUMS` covers the exact input manifest, outputs, ledgers,
provenance and resource evidence; `COMPLETE` was fsynced last before a
same-filesystem rename. Independent semantic validation reports SHA-256
`814a2fcdfb1157e1c7b3c42a7b81a513679eb64ad5dd9d2244abe0c137620d04`. A post-publication rerun performs no
recompute or mutation.

## Interpretation boundary

This release freezes host topology, support-collapsed clades, all-host mappings
and alternative partitions **before** phage association. Mash branches describe
whole-assembly k-mer similarity; local core-SNP branches describe the sampled
host clonal frame conditional on mapping, mask and reference. Unsupported or
reference-sensitive structure is ambiguous, not evidence of absence and not a
license to choose whichever partition strengthens a phage result.
