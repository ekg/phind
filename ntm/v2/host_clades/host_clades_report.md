# NTM v2 — host clades on the QC-passed cohort

Method and results for `ntm/v2/host_clades/` (task ntm-v2-host).

## Method

Per the chat-agent update (2026-08-16), subsetting v1 host distances is wrong
for v2 (only ~1,057 of 13,122 NCBI assemblies overlap the v1 cohort; v1 has no
run-assemblies). v2 uses a **full re-sketch** over the QC-passed cohort:

1. filelist = every `*.pansn.fa.gz` in `genomes/canonical_objects/`
   (+ `genomes/run_assemblies/` when present — none yet, so this ships the
   NCBI-subset version; the pipeline extends unchanged).
2. `mash sketch -k 21 -s 10000 -l filelist.txt -o host -p 32`
   (identical parameters to v1 `ntm/scripts/host_clades_mash.py`).
3. `mash dist -p 64 host.msh host.msh > host.dist`
   (full all-vs-all; N² lines).
4. Cluster with the same logic as v1: union-find connected components at
   `mash dist <= 0.05` (~95% ANI, species level) → components sorted by size
   → fresh clade ids `host_clade_0001..` for v2.
5. Annotations from the QC-passed accession list
   (`NTM_QC_passed_accession_list.tsv`): `species` verbatim, `genus` = first
   token of species, `organism` = species (the v2 accession list has no
   separate organism column). v1's manifest-based annotations would not cover
   the new genomes.

Script: `ntm/v2/scripts/host_clades_mash_v2.py` (reuses v1 clustering code).

## Cohort

- **13,122** canonical NCBI assemblies (GCA-preferred; 332 GCF-only), all
  downloaded/linked, all present in the QC-passed accession list (0 missing).
- run-assemblies: **0** present in `genomes/run_assemblies/` — table schema is
  designed so a re-run over an extended filelist (canonical + runs) just works;
  run-id rows (9,543) are already annotated in the accession list.

### Composition (accession-list species calls, canonical NCBI cohort)

- TB-complex (species containing 'tuberculosis', excl. paratuberculosis): **7244 (55.2%)** — kept verbatim per collaborator QC list (user decision pending; no filtering applied)
- abscessus complex (species containing 'abscessus'): **2336 (17.8%)**
- distinct species strings: 264

Top 15 species:

| n | species |
|---:|---|
| 7237 | Mycobacterium tuberculosis |
| 1046 | Mycobacterium ulcerans |
| 1046 | Mycobacteroides abscessus subsp. abscessus |
| 703 | Mycobacteroides abscessus |
| 413 | Mycobacteroides abscessus subsp. massiliense |
| 282 | Mycobacterium sp. |
| 234 | Mycobacterium avium subsp. hominissuis |
| 215 | Mycobacterium avium subsp. paratuberculosis |
| 185 | Mycobacterium kansasii |
| 174 | Mycobacteroides abscessus subsp. bolletii |
| 138 | Mycolicibacterium sp. |
| 83 | Mycobacterium marinum |
| 77 | Mycobacteroides chelonae |
| 67 | Mycobacterium intracellulare |
| 61 | Mycolicibacterium smegmatis |

## Results

Clustering from `host.dist` (39,575,994 pairs; union ops: 21,023,097) at
dist ≤ 0.05 (~95% ANI, single-linkage connected components — v1 logic):

- **genomes clustered: 13,122** (all filelist accessions assigned)
- **clades: 1,044** (fresh ids `host_clade_0001`–`host_clade_1044`, size-desc)
- singletons: 1,006 (mostly singletons/small clades — no join partner ≤ 0.05)

Clade size distribution (top):

| n | # clades |
|---:|---:|
| 7465 | 1 |
| 2336 | 1 |
| 1139 | 1 |
| 507 | 1 |
| 187 | 1 |
| 152 | 1 |
| 61 | 2 |
| 25 | 1 |
| 18 | 1 |
| 17 | 1 |
| 15 | 1 |
| 10 | 1 |
| 9 | 1 |
| 8 | 2 |
| 7 | 4 |
| ... | 1,024 more clades (mostly singletons) |

Top clades (dominant species):

| clade | n | dominant |
|---|---|---|
| host_clade_0001 | 7465 (56.9%) | Mycobacterium tuberculosis (TB complex — kept verbatim per collaborator QC list) |
| host_clade_0002 | 2336 (17.8%) | Mycobacteroides abscessus subsp. abscessus (abscessus complex) |
| host_clade_0003 | 1139 | Mycobacterium ulcerans (+ M. marinum) |
| host_clade_0004 | 507 | Mycobacterium avium subsp. hominissuis (+ M. avium complex) |
| host_clade_0005 | 187 | Mycobacterium kansasii |
| host_clade_0006 | 152 | Mycobacterium intracellulare (+ MAC) |

## Validation

- [x] **Row completeness:** 13,122 rows == 13,122 filelist accessions; 0
  duplicate accessions; inner-join over the NCBI cohort is exact (100%).
  The full accession list (33,082 rows) is NOT fully represented by design:
  9,543 run-assembly ids have no genomes ingested yet (0 in
  `genomes/run_assemblies/`), and GCA/GCF twins (10,417 numerics) collapse to
  one canonical genome each — both documented above. Missing=0, so nothing
  was dropped silently.
- [x] **Annotation join:** all 13,122 accessions resolve to a species row in
  `NTM_QC_passed_accession_list.tsv`; 0 missing.
- [x] **Clade size sanity:** top clade host_clade_0001 = 56.9% (>40%) but it is
  TB-complex (M. tuberculosis-labelled), which the collaborator QC list keeps
  verbatim per chat-agent update 2026-08-16 — flagged, explained, not filtered.
  No other clade exceeds 40% (next is abscessus complex at 17.8%).
- [x] **Method parity:** identical mash parameters (k=21, s=10000) and
  clustering threshold (dist ≤ 0.05) as v1 `ntm/scripts/host_clades_mash.py`;
  fresh clade ids for v2.
- [x] **Script validation:** `host_clades_mash_v2.py` reports `validation: PASS`
  (see `cluster.log`).

Summary table (from task log):

```
clades: 1044; assigned: 13122 accessions
TB-complex: 7244 (55.2%) — kept verbatim
abscessus complex: 2336 (17.8%)
top 15 species: M. tuberculosis 7237 | M. ulcerans 1046 |
  M. abscessus subsp. abscessus 1046 | M. abscessus 703 |
  M. abscessus subsp. massiliense 413 | Mycobacterium sp. 282 |
  M. avium subsp. hominissuis 234 | M. avium subsp. paratuberculosis 215 |
  M. kansasii 185 | M. abscessus subsp. bolletii 174 |
  Mycolicibacterium sp. 138 | M. marinum 83 | M. chelonae 77 |
  M. intracellulare 67 | M. smegmatis 61
```

## Files

| file | contents |
|---|---|
| `host_clades.tsv` | accession \| host_clade_id \| species \| genus \| organism (one row per cohort genome) |
| `host_clade_summary.tsv` | host_clade_id \| count \| dominant_species \| species_distribution |

## Extending to run-assemblies

When `genomes/run_assemblies/{acc}/{acc}.pansn.fa.gz` arrive (ingest via
`ntm/v2/scripts/ingest_run_assemblies.py`), regenerate the filelist
(`find canonical_objects run_assemblies -name '*.pansn.fa.gz' | sort`),
re-sketch + dist, and re-run `host_clades_mash_v2.py` unchanged. Clade ids are
recomputed fresh; annotations resolve via the accession list (run ids already
present, data_source=ASSEMBLY).
