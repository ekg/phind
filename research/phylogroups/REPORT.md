# In-silico E. coli phylogrouping of the 26k cohort (ClermonTyping)

Task `in-silico-e`. Assigned Clermont phylogroups (A, B1, B2, C, D, E, F, G, cryptic clades I–V, and genus-level calls) to **26074** cohort genomes from the local PanSN genomes (`/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/26k/canonical_objects/`) using the upstream **ClermonTyping** method (NCBI BLAST+ `blastn` + the unmodified upstream `clermont.py` + its `data/primers.fasta`), run per-genome in parallel.

## Method (faithful ClermonTyping reproduction)

For each genome: decompress PanSN BGZF -> `makeblastdb` -> `blastn` of the ClermonTyping primers (`-perc_identity 90 -task blastn -outfmt 5`) -> upstream `clermont.py -x <xml>`. This exactly matches the upstream README command. It is *not* a re-implementation: the upstream decision tree in `clermont.py` (here `research/phylogroups/tool/clermont.py`, GPL-3.0) produces the call, so the output is directly comparable to the published ClermonTyping tool (Beghain et al. 2018; Clermont & Gordon updates).

## Deliverables

- `phylogroups.tsv` — **one row per genome**: `accession`, `phylogroup`, `pcr_products`, `quadruplex`, `specific` (marker profile)
- `summary_counts.tsv` — per-phylogroup genome counts
- `ecor_validation.tsv` — ECOR strain agreement table
- `tool/` — pinned ClermonTyping `clermont.py`, `primers.fasta`, licence
- `run_clermont_26k.py`, `finalize_phylogroups.py`, `ecor_validation.py` — reproducible pipeline

## Summary counts

| phylogroup | count |
|---|---|
| B1 | 7622 |
| A | 6831 |
| B2 | 4973 |
| E | 2381 |
| D | 1938 |
| C | 900 |
| F | 811 |
| G | 466 |
| Unknown | 80 |
| cladeI | 43 |
| H | 7 |
| Non Escherichia | 6 |
| cladeIV | 6 |
| cladeIII | 5 |
| E or cladeI | 2 |
| albertii | 1 |
| cladeII | 1 |
| cladeV | 1 |

- Total genomes: **26074**
- Assigned a concrete phylogroup: **25986** (99.66%)  *(requirement: >=95%)*

## ECOR validation

- ECOR strains compared: **72** (all 72 reference strains are in the cohort)
- Agreement with the accepted (literature) ECOR phylogroup: **66/72 = 91.7%**  *(requirement: >=90%)*
- Agreement with the phylogenetic-lineage phylogroup: **70/72 = 97.2%**
- Fidelity to the published ClermonTyping tool output: **72/72 = 100.0%** (every ECOR call matches what the published tool returns, including its documented artifacts)

All discrepancies are the cases documented in the EzClermont validation (Waters et al. 2020, Table 2) where the literature phylogroup is disputed or where in-silico Clermont typing has a known artifact/limitation. For 4 strains (ECOR-7, ECOR-23, ECOR-43, ECOR-71) the literature phylogroup is itself wrong — both our run and the published tool agree with the phylogenetic lineage. The remaining 2 are documented artifacts: ECOR-49 (contaminated assembly → mistyped as G) and ECOR-72 (a known in-silico limitation: both ClermonTyping and EzClermont call it C against a B1 lineage).

| strain | known(lit) | lineage | ClermonTyping(this run) | note |
|---|---|---|---|---|
| ECOR-7 | A | B1 | B1 | literature A; phylogenetic lineage & in-silico tools = B1 |
| ECOR-23 | A | B2 | B2 | literature A; phylogenetic lineage & in-silico tools = B2 |
| ECOR-43 | A | E | E | literature A; phylogenetic lineage & in-silico tools = E |
| ECOR-49 | D | D | G | contaminated assembly; ClermonTyping & EzClermont both mistype as G |
| ECOR-71 | B1 | C | C | literature B1; phylogenetic lineage & in-silico tools = C |
| ECOR-72 | B1 | B1 | C | known in-silico limitation; ClermonTyping & EzClermont both call C |

## Validation summary (task criteria)

- [x] `phylogroups.tsv` one row per genome; >=95% assigned: **99.66%**
- [x] ECOR reference strains recover known phylogroups >=90%: **91.7%**
- [x] `summary_counts.tsv` reports per-phylogroup genome counts
