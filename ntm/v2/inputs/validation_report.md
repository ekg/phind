# NTM v2 collaborator input validation report

Generated: (see coverage.tsv mtime) — script `ntm/v2/scripts/validate_inputs.py`

## 1. File integrity (sha256)

| file | size (bytes) | sha256 |
|---|---|---|
| `NTM_QC_passed_accession_list.tsv` | 2599713 | `da26584badaa4fca7cd8434abaf0bdf70116152012a1e9217991bf20092241ce` |
| `NTM_QC_passed_prophage_master_manifest.tsv` | 17507506 | `c88686661928c66ab9638da77dd0d615e87bc70975ad859a81e7644142c83386` |

## 2. Schema — NTM_QC_passed_accession_list.tsv

- rows (excl. header): **33082**
- columns (4): `genome_id`, `accession`, `species`, `data_source`

| column | type | notes |
|---|---|---|
| genome_id | str | identifier, identical to accession for run-assemblies |
| accession | str | NCBI accession: GCA_/GCF_ assembly or ERR/SRR/DRR run id |
| species | str | collaborator species call (may include `[tuberculosis]`) |
| data_source | str | `NCBI` (assembly) or `ASSEMBLY` (run-assembly) |

| accession prefix | count |
|---|---|
| GCA | 12790 |
| GCF | 10749 |
| ERR | 5620 |
| SRR | 3915 |
| DRR | 8 |

| data_source | count |
|---|---|
| NCBI | 23539 |
| ASSEMBLY | 9543 |

- unique NCBI assemblies (by numeric part): **13122**
  - with both GCA+GCF twins: 10417
  - GCA-only: 2373
  - GCF-only: 332
- unique run ids (ERR/SRR/DRR): **9543**
- total unique genome ids (numeric + runs): **22665**

### Species composition (accession list, all rows)

| species | count |
|---|---|
| Mycobacterium tuberculosis | 13216 |
| Mycobacteroides abscessus | 6291 |
| Mycobacteroides abscessus subsp. abscessus | 2077 |
| Mycobacterium avium | 1932 |
| Mycobacterium ulcerans | 1854 |
| Mycobacterium intracellulare | 1310 |
| Mycobacteroides abscessus subsp. massiliense | 803 |
| Mycobacterium sp. | 610 |
| Mycobacteroides abscessus subsp. bolletii | 609 |
| Mycobacterium avium subsp. hominissuis | 445 |
| Mycobacterium kansasii | 375 |
| Mycolicibacterium smegmatis | 259 |
| Mycobacterium avium subsp. paratuberculosis | 254 |
| Mycobacteroides chelonae | 245 |
| Mycolicibacterium sp. | 234 |
| Mycobacterium marinum | 170 |
| Mycolicibacterium fortuitum | 102 |
| Mycobacterium paraintracellulare | 70 |
| Mycobacterium colombiense | 64 |
| Mycobacterium canetti | 61 |
| Mycobacterium intracellulare subsp. chimaera | 54 |
| Mycolicibacterium nivoides | 54 |
| Mycolicibacterium septicum | 53 |
| Mycolicibacterium senegalense | 43 |
| Mycobacterium riyadhense | 42 |
| Mycobacterium genavense | 42 |
| Mycobacterium intracellulare subsp. intracellulare | 39 |
| Mycolicibacterium neoaurum | 33 |
| Mycobacteroides immunogenum | 32 |
| Mycobacterium avium subsp. avium | 31 |

Total distinct species strings: 274

### TB-complex subset (species containing 'tuberculosis', excl. paratuberculosis)

- rows: 13227 (of 33082)
- unique NCBI assemblies: 7377
- unique run ids: 0

> v1 excluded MTBC by design; v2 default = keep collaborator list verbatim (pending user confirmation). Documented here for the scope decision.

## 3. Schema — NTM_QC_passed_prophage_master_manifest.tsv

- rows (excl. header): **51004**
- columns (23): `genome_id`, `species`, `species_confidence`, `data_source`, `accession`, `fasta_path`, `busco_complete_pct`, `busco_fragmented_pct`, `checkm2_completeness_pct`, `checkm2_contamination_pct`, `assembly_contigs`, `assembly_length_bp`, `assembly_n50_bp`, `ani_pct_to_closest_type_strain`, `ani_closest_type_strain_species`, `has_prophage`, `prophage_id`, `prophage_contig`, `prophage_start`, `prophage_end`, `prophage_length_bp`, `prophage_transposable_element`, `prophage_family`

| column | type | notes |
|---|---|---|
| genome_id | str | `{accession}_{assembly_name}_genomic` for NCBI; run id for ASSEMBLY |
| species | str | collaborator species call |
| species_confidence | str | empty / MODERATE / HIGH |
| data_source | str | NCBI or ASSEMBLY |
| accession | str | GCA_/GCF_ assembly or ERR/SRR/DRR run id |
| fasta_path | str | collaborator Iridis cluster path (NOT public) |
| busco_complete_pct | float | QC metric |
| busco_fragmented_pct | float | QC metric |
| checkm2_completeness_pct | float | QC metric |
| checkm2_contamination_pct | float | QC metric |
| assembly_contigs | int | contig count |
| assembly_length_bp | int | assembly size |
| assembly_n50_bp | int | N50 |
| ani_pct_to_closest_type_strain | float | ANI to nearest type strain |
| ani_closest_type_strain_species | str | type strain species |
| has_prophage | bool (True/False) | flag; False rows have empty prophage_* cols |
| prophage_id | str | `{genome_id}_prophageN` (N 1-based) |
| prophage_contig | str | WGS contig accession (GCA rows bare, GCF rows NZ_-prefixed); NODE_* for runs |
| prophage_start | float (int) | 1-based INCLUSIVE start |
| prophage_end | float (int) | 1-based INCLUSIVE end |
| prophage_length_bp | float (int) | end-start+1 |
| prophage_transposable_element | bool | transposable element flag |
| prophage_family | str | geNomad family call |

- rows with wrong column count: **0** (of 51004)

- has_prophage value counts: {'False': 20839, 'True': 30165}

- False rows: all prophage_* cols empty: True
- True rows with any empty prophage_* col: **0** (row accounting)

## 4. Unique genome & prophage counts (dedup rule: prefer GCA rows)

| source | prophage rows (True) |
|---|---|
| GCA | 8183 |
| GCF | 7241 |
| ERR | 9938 |
| SRR | 4789 |
| DRR | 14 |

| kept prophage rows after dedup | source | count |
|---|---|---|
| | GCA | 8183 |
| | GCF | 319 |
| | ERR | 9938 |
| | SRR | 4789 |
| | DRR | 14 |
- total kept prophage rows: **23243**
- GCF-only numeric assemblies contributing prophages: **143**
- unique prophage-bearing genomes after dedup: **9434** (3527 NCBI numerics + 5907 run-assemblies)

- prophage_id uniqueness after dedup: **23243 unique ids for 23243 rows**
  - no duplicate prophage_ids ✔

## 5. Prophage length distribution & coordinate sanity

- rows with begin>end or unparseable coords: **0**
- rows where length != end-start+1 (1-based inclusive): **0**

| length stat | bp |
|---|---|
| min | 237 |
| median | 20311 |
| mean | 21933.0 |
| max | 96653 |

| length bin | count |
|---|---|
| 0-5000 | 3318 |
| 5000-10000 | 6428 |
| 10000-20000 | 5267 |
| 20000-40000 | 11357 |
| 40000-60000 | 3619 |
| 60000-100000 | 176 |
| 100000-1000000000 | 0 |

| prophage_family | count |
|---|---|
| Siphoviridae | 26633 |
| Unknown | 2313 |
| Myoviridae | 881 |
| Siphoviridae / Podoviridae | 167 |
| Myoviridae / Siphoviridae | 94 |
| Podoviridae | 63 |
| Fuselloviridae | 10 |
| Siphoviridae / Fuselloviridae | 2 |
| Myoviridae / Siphoviridae / Podoviridae | 2 |

## 6. Coverage vs v1 (7,352 local genomes)

- v1 canonical_objects dirs: 7352
  - with pansn.fa.gz: 7303
  - empty dirs: 49
- v1 unique numeric accessions: 7325

- v2 NCBI assemblies already local, **GCA-preferred exact-match** (dir name == preferred accession): **1150**
- v2 NCBI assemblies MISSING locally (need download): **11972** (of which GCF-only numerics: 219)
- prefix-agnostic overlap (numeric id present locally under either GCA/GCF): **5423**
- v1 assemblies not in v2 list: **1902**

> The download task (ntm-v2-download-12-030) uses the GCA-preferred exact-match definition: link local dirs where the preferred accession exists; download the rest.

> Reconciliation vs earlier chat-agent estimate (13,123 unique numerics; 1,057 local; 12,030 to download): precise recount gives **13,122** unique numerics, **1,150** GCA-preferred local, **11,972** to download. Small deltas from the estimate; coverage.tsv is the ground truth.

- v1 host_clades.tsv unique accessions: 7303 (numeric: 7295)
- host_clades ∩ v2 NCBI assemblies: **5407**
- host_clades ∩ v1-local-and-v2: 5407

- v1 prophage CSV: 10438 rows, 4064 unique genomes
- v2 prophage-bearing genomes (after dedup): 3527 NCBI + 5907 runs
- overlap (v1 ∩ v2 prophage-bearing numerics): **2866**

## 7. Per-genome prophage count: v1 geNomad vs v2 manifest

- v1 genomes with prophages also in v2 (comparable): 2867
- per-genome prophage count identical: **1021**
- count differs: 1846

| v1 genome | v1 count | v2 count |
|---|---|---|
| GCF_014263335.1 | 3 | 13 |
| GCF_010722915.1 | 1 | 9 |
| GCF_000184435.1 | 1 | 8 |
| GCF_015689175.1 | 2 | 9 |
| GCF_025263565.1 | 3 | 10 |
| GCF_900136205.1 | 9 | 16 |
| GCF_001296255.1 | 4 | 10 |
| GCF_001296275.1 | 4 | 10 |
| GCF_001296415.1 | 5 | 11 |
| GCF_010723135.1 | 3 | 9 |
| GCF_900130755.1 | 9 | 3 |
| GCF_900141685.1 | 3 | 9 |

### Interpretation of count differences

v1 ran geNomad v1.12.0 (find-proviruses default, end-to-end mode) on the v1 cohort. The collaborator ran their own geNomad pipeline (version/params unknown) on their assemblies. Both call prophages per genome, but counts differ for 1846 of 2867 shared prophage-bearing genomes. Contributing factors:

- collaborator manifest retains very short predictions (min 237 bp; v1 min was 3,220 bp) — e.g. GCF_001296255.1: 4 (v1) vs 10 (v2, incl. 1,428 bp and 4,800 bp calls);
- geNomad version/parameter drift between runs (e.g. GCF_900130755.1: 9 (v1) vs 3 (v2));
- assembly version differences (GCA vs GCF mirror or updated versions).

v2 keeps the collaborator's calls verbatim (their QC-passed list is the v2 denominator); v1 counts are recorded here only as a cross-check.

## 8. Scaffold-name mapping decision (collaborator contig → PanSN)

v1 PanSN header: `{assembly_acc}#1#{contig}` (e.g. `GCA_000523695.1#1#JAOB01000032.1`).
v1 `ntm_prophages.csv` stored `source_seq` = the full PanSN header.

Collaborator `prophage_contig` values:

- run: NODE_*: 14741
- GCA-row: bare WGS: 8183
- GCF-row: NZ_-prefixed WGS: 7213
- GCF-row: NC_ chromosome (complete genomes): 28

**Decision:** for NCBI rows, map collaborator contig → PanSN using the LOCAL genome dir's prefix: `{local_acc}#1#{contig}` where the contig suffix matches the dir prefix convention (GCA dirs → bare WGS accession; GCF dirs → `NZ_`-prefixed). v1 GCA dirs (canonical) match GCA-row contigs directly; GCF rows drop `NZ_` when mapped onto a GCA dir and keep it on a GCF dir.

Worked example on ≥3 genomes (coordinate check vs local .fai):

| genome | prophage_id | contig (collab) | PanSN (mapped) | fai len | start | end | within? |
|---|---|---|---|---|---|---|---|
| GCF_041063225.1 | GCA_041063225.1_ASM4106322v1_genomic_prophage1 | JBFUXV010000022.1 | `GCF_041063225.1#1#NZ_JBFUXV010000022.1` | 73764 | 11273 | 44377 | OK |
| GCF_041063225.1 | GCF_041063225.1_ASM4106322v1_genomic_prophage1 | NZ_JBFUXV010000022.1 | `GCF_041063225.1#1#NZ_JBFUXV010000022.1` | 73764 | 11273 | 44377 | OK |
| GCF_039023525.1 | GCA_039023525.1_ASM3902352v1_genomic_prophage1 | JBCHKV010000010.1 | `GCF_039023525.1#1#NZ_JBCHKV010000010.1` | 211844 | 153591 | 174537 | OK |
| GCF_039023525.1 | GCF_039023525.1_ASM3902352v1_genomic_prophage1 | NZ_JBCHKV010000010.1 | `GCF_039023525.1#1#NZ_JBCHKV010000010.1` | 211844 | 153591 | 174537 | OK |
| GCF_025823205.1 | GCA_025823205.1_ASM2582320v1_genomic_prophage1 | JACKVH010000022.1 | `GCF_025823205.1#1#NZ_JACKVH010000022.1` | 101242 | 2017 | 6988 | OK |
| GCF_025823205.1 | GCA_025823205.1_ASM2582320v1_genomic_prophage2 | JACKVH010000013.1 | `GCF_025823205.1#1#NZ_JACKVH010000013.1` | 213945 | 115524 | 124593 | OK |
| GCF_025823205.1 | GCA_025823205.1_ASM2582320v1_genomic_prophage3 | JACKVH010000012.1 | `GCF_025823205.1#1#NZ_JACKVH010000012.1` | 1907946 | 1701227 | 1703986 | OK |
| GCF_025823205.1 | GCF_025823205.1_ASM2582320v1_genomic_prophage1 | NZ_JACKVH010000022.1 | `GCF_025823205.1#1#NZ_JACKVH010000022.1` | 101242 | 2017 | 6988 | OK |

- spot-checks performed on overlapping v1-local genomes: 8 (fai-resolved: 8)
- contig missing from .fai: 0 | outside bounds: 0

## 9. Row accounting (no silent drops)

- manifest rows total: 51004
- True rows: 30165 → kept after dedup: 23243
  - GCA rows kept: 8183 (GCA preferred)
  - GCF rows kept: 319 (GCF-only numerics: 143)
  - run rows kept: 14741
- False rows (genomes without prophage): 20839 — retained in manifest, not prophage set
- True rows dropped by dedup (GCF twins of GCA): 6922
- rows with wrong column count: 0
- True rows with empty prophage cols: 0

Every manifest row is accounted for: kept (prophage set), retained as no-prophage genome, or dropped by explicit dedup rule (GCF twin where GCA preferred).

## 10. Species composition (manifest)

| species | genomes (rows) |
|---|---|
| Mycobacteroides abscessus | 15887 |
| Mycobacterium tuberculosis | 13238 |
| Mycobacteroides abscessus subsp. abscessus | 5599 |
| Mycobacterium avium | 2063 |
| Mycobacteroides abscessus subsp. massiliense | 1939 |
| Mycobacterium ulcerans | 1854 |
| Mycobacterium intracellulare | 1740 |
| Mycobacteroides abscessus subsp. bolletii | 940 |
| Mycobacteroides chelonae | 718 |
| Mycobacterium sp. | 678 |
| Mycobacterium avium subsp. hominissuis | 570 |
| Mycobacterium kansasii | 517 |
| Mycolicibacterium sp. | 355 |
| Mycolicibacterium smegmatis | 274 |
| Mycobacterium avium subsp. paratuberculosis | 254 |
| Mycolicibacterium fortuitum | 246 |
| Mycobacterium marinum | 237 |
| Mycobacteroides immunogenum | 195 |
| Mycolicibacterium nivoides | 148 |
| Mycolicibacterium senegalense | 135 |
| Mycobacterium intracellulare subsp. chimaera | 96 |
| Mycobacterium paraintracellulare | 91 |
| Mycobacterium kubicae | 86 |
| Mycobacterium paragordonae | 81 |
| Mycobacterium colombiense | 74 |
| Mycobacterium riyadhense | 69 |
| Mycobacterium canetti | 65 |
| Mycolicibacterium neoaurum | 63 |
| Mycolicibacterium goodii | 62 |
| Mycolicibacterium gilvum | 61 |

