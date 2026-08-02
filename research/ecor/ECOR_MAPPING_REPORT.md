# ECOR mapping — E. coli Reference collection in the all-prophage collection

**Task:** `ecor-mapping` — map/highlight the ECOR elements within the
all-prophage collection.

**Status:** COMPLETE — 72/72 ECOR strains resolved to RefSeq assemblies present
in the 26k cohort; 300 prophage elements mapped; tag merge rate 1.0 (300/300).

## 1. What ECOR is

The **E. coli Reference (ECOR) collection** is the historical reference set of
72 *Escherichia coli* isolates (ECOR-1 … ECOR-72) assembled by Ochman &
Selander (1984) to represent the known genetic diversity of the species
(J Bacteriol 157:690–693). The canonical reference genome set for all 72
strains is:

> Patel IR, Gangiredla J, Lacher DW, Mammel MK, et al. (2018) **Draft Genome
> Sequences of the *Escherichia coli* Reference (ECOR) Collection.**
> Microbiol Resour Announc 7:e01133-18. doi:10.1128/mra.01133-18
> — BioProject **PRJNA230969**, WGS accessions **QOWM00000000.1 … QOZE00000000.1**
> (ECOR-59 was deposited out-of-order as QOZF00000000.1).

## 2. Crosswalk method (strain → GCF accession)

The prophage collection (`26k_prophage1.csv`, `full_prophages.fa`) is keyed on
RefSeq **GCF_** assembly accessions. The ECOR reference genomes were deposited
as GenBank WGS assemblies, each of which NCBI has paired with an **identical**
RefSeq assembly. The crosswalk:

1. **Strain number source:** NCBI `assembly_summary` (taxid 562) rows for
   BioProject PRJNA230969, whose `infraspecific_name` is the authoritative
   `strain=MOD1-ECOR<n>` label (all 72 strains, exactly once).
2. **Paper cross-check:** every NCBI `wgs_master` was matched against the
   paper's Table 1 (embedded in `build_ecor_manifest.py` as
   `ECOR_WGS_PAPER`). 0 discrepancies.
3. **GCA → GCF:** each GenBank assembly (`GCA_0033xxxxx.1`) has an identical
   paired RefSeq assembly (`GCF_0033xxxxx.1`, `paired_asm_comp=identical`),
   confirmed present in the RefSeq assembly summary.
4. **Cohort membership:** all 72 GCF accessions are present in
   `26k_ecoli_accession.txt` (26,077 accessions), `all_valid_accessions.txt`,
   and the frozen pilot frame `manifests/pilot-cohorts-v1/frame.tsv.gz`
   (72/72).
5. **Prophage join:** `26k_prophage1.csv` rows where `genome == GCF_...`
   provide `prophage_id`, `source_contig` (scaffold), and 1-based inclusive
   `begin`/`end` coordinates (per prophage-semantics v2 policy
   `C1_RAW_1_BASED_CLOSED`; `length = end − begin + 1`).
6. **FASTA verification:** every mapped `prophage_id` was looked up in
   `prophage_homology_survey/full_prophages.fa` (132,393 sequences); sequence
   length equals the manifest `length` column for **all 300** elements
   (0 mismatches).

## 3. Results

### Reconciliation (72 canonical strains)

| status | count |
|---|---|
| matched (GCF resolved, in cohort, prophage elements mapped) | **72** |
| gcf-resolved-not-in-cohort | 0 |
| not-found-by-accession | 0 |
| **total** | **72** |

Coverage vs the canonical ~72-strain ECOR list: **72/72 = 100%**.

### Elements

| metric | value |
|---|---|
| ECOR prophage elements (manifest rows) | **300** |
| elements present in `full_prophages.fa` | **300** |
| **tag merge rate** | **1.0** (300/300) |
| sum of element base pairs | 6,879,652 |
| elements per strain | 1 … 9 (mean 4.17) |

Bidirectional consistency: the set of `full_prophages.fa` headers whose GCF is
an ECOR GCF (300) equals the set of `26k_prophage1.csv` prophage ids for ECOR
GCFs (300) — no orphan headers, no missing elements.

### Ambiguity: alternate ECOR-labeled assemblies in the cohort

64 of 72 strain numbers have **additional** ECOR-labeled GCF assemblies in the
cohort (66 total alternates). These are almost all BioProject **PRJNA224116**
(Messerer M, Fischer W, Schubert S, PLoS One 2017; the 2018 reference paper
**explicitly cautions against this dataset**, flagging possible contamination —
multiple O/H serotyping loci — in 32 of its ECOR assemblies), plus
`GCF_001865905.1` (Karolinska) and `GCF_002197975.1` (UMass) for ECOR-31/13.

**Decision:** the canonical mapping uses the 2018 PRJNA230969 reference
assemblies for all 72 strains. Alternates are recorded (never merged into the
manifest) in `ecor_alternate_assemblies.tsv` with the paper's contamination
flag, so downstream consumers can decide whether to consider them.

## 4. Deliverables (under `research/ecor/`)

| file | content |
|---|---|
| `ecor_manifest.csv` | one row per ECOR prophage element: `ecor_strain, assembly_accession (GCF), gca_accession, wgs_master, prophage_id, source_contig, start, end, length, transposable, taxonomy` (300 rows) |
| `ecor_strain_reconciliation.tsv` | per-strain status, GCF, element counts, tag merge rate (72 rows) |
| `ecor_reconciliation_summary.json` | exact counts + tag merge rate (machine-readable) |
| `ecor_leaf_tags.tsv` | **ECOR-tag per MASH-tree leaf**: all 132,393 `full_prophages.fa` headers with `is_ecor` (TRUE/FALSE) and `ecor_strain` — direct join key for highlighting the MASH tree (`ecor-highlighted-inspection`) |
| `ecor_alternate_assemblies.tsv` | alternate ECOR-labeled GCF assemblies in the cohort, with contamination flags |
| `build_ecor_manifest.py` | reproducible builder (reads NCBI assembly summaries + cohort inputs) |

## 5. Caveats

- The 2018 ECOR reference assemblies are **draft (Contig-level)** and their
  RefSeq copies are "identical" pairs (`refseq_category=na`), not annotated
  Reference genomes. Coordinates therefore refer to the draft scaffolds
  (`NZ_QOWM01…`).
- Coordinates are 1-based inclusive per the project's prophage-semantics v2
  policy (`C1_RAW_1_BASED_CLOSED`).
- If another ECOR strain's genome was sequenced at higher quality elsewhere
  (e.g., complete chromosomes), it is **not** the 2018 reference set and is out
  of scope for this canonical mapping; alternates are documented above.
- The 11 extraction FAILs recorded in `prophage_homology_survey/extraction_log.txt`
  (GCF_001291365.1, GCF_003886435.1, GCF_012029685.1) involve **no** ECOR
  accessions.
