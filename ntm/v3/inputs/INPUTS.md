# NTM v3 expansion source inputs — BV-BRC phigaro QC-passed cohort

Imported: 2026-09-12 (task `import-ntm-v3`). These two CSVs are the **frozen
source inputs** for the v3 all-inclusive NTM expansion (see `../README.md`),
exactly analogous to how the v2 collaborator inputs are recorded in
`ntm/v2/inputs/coverage.tsv` / `validation_report.md`. Original filenames are
kept unchanged — the date stamp **2026-09-09** is the provenance (BV-BRC
phigaro QC-passed export, 2026-09-09).

## 1. File integrity (sha256)

| file | size (bytes) | data rows | sha256 |
|---|---:|---:|---|
| `ntm_qc_passed_phigaro_summary_20260909.csv` | 2,241,769 | 26,499 | `c9968d395193b9046b09146a21a45fffb3dbc34903c63cdc5c2338d0037c877a` |
| `ntm_qc_passed_phigaro_coordinates_20260909.csv` | 7,693,295 | 47,994 | `503bcd78862d13687595e548a7f42ac2b7bd7628f996748c51a9c39c83616420` |

Both files are **byte-exact as exported** (CRLF line endings, 0x0D 0x0A —
the BV-BRC export tool emits DOS line endings). The repo-root `.gitattributes`
marks these two paths `whitespace=cr-at-eol` so git's whitespace lint does not
flag the carriage returns; do **not** "normalize" these files to LF — that
would break the sha256 receipts above.

## 2. Schema — `ntm_qc_passed_phigaro_summary_20260909.csv`

Per-genome summary, one row per genome (26,499 rows).

Columns (6): `genome_id,source,accession,species,has_phigaro_result,prophage_count`

| column | type | notes |
|---|---|---|
| genome_id | str | identifier; format depends on `source` (see §4) |
| source | str | `BV-BRC` / `NCBI` / `ASSEMBLY` |
| accession | str | == `genome_id` for `BV-BRC` rows; bare `GCA_/GCF_` version for `NCBI` rows |
| species | str | BV-BRC-style label, often with type-strain name (see §5) |
| has_phigaro_result | bool | `True` for all 26,499 rows |
| prophage_count | int | 0–36; 9,382 genomes have 0 |

## 3. Schema — `ntm_qc_passed_phigaro_coordinates_20260909.csv`

Per-prophage coordinates, one row per predicted prophage (47,994 rows).

Columns (10): `genome_id,source,accession,species,scaffold,prophage_id,begin,end,transposable,taxonomy`

| column | type | notes |
|---|---|---|
| genome_id / source / accession / species | str | as in summary |
| scaffold | str | contig the prophage sits on; format varies by source (see §4) |
| prophage_id | str | `<scaffold>_prophage<N>`; unique across all 47,994 rows |
| begin / end | int | 1-based inclusive; `end ≥ begin` in all rows |
| transposable | str | `True` / `False` / **empty** (22,225 rows have no value) |
| taxonomy | str | phigaro family call: `Siphoviridae` (41,844), `Unknown` (3,719), `Myoviridae` (1,885), mixed calls (`Siphoviridae / Podoviridae`, …) |

Scaffold format by source:

| source | rows | scaffold style | example |
|---|---:|---|---|
| ASSEMBLY (run asm) | 33,303 | SPAdes/ENA `NODE_<n>_length_…_cov_…` contig names | `NODE_14_length_150130_cov_107.750611` |
| NCBI | 7,977 | versioned accessions, mostly `NZ_`/`NW_`-prefixed or plain | `NZ_GG770559.1`, `CP000479.1` |
| BV-BRC | 465 | `accn\|`-prefixed accessions, unversioned | `accn\|CP118870` |
| NCBI (misc) | 44 | other contig names | — |
| **total** | **47,994** | | |

## 4. Verified cohort properties (re-checked post-move, 2026-09-12)

- 26,499 genomes; **all** `has_phigaro_result=True`
- **17,117** genomes carry ≥ 1 prophage; **47,994** prophages total
- **9,382** genomes have `prophage_count=0` (no coordinates rows — expected)
- coordinates `genome_id`s are a strict subset of summary `genome_id`s (0 orphans);
  per-genome coordinates row counts equal the summary `prophage_count` exactly (0 mismatches)

Composition by `source` / identifier format — **the identifiers are a mix of three
namespaces, not a single BV-BRC ID scheme**:

| source | genome_id format | genomes | of which prophage-bearing |
|---|---|---:|---:|
| `ASSEMBLY` | ENA/SRA run accession (`ERR/SRR/DRR`) | 17,920 | 11,323 |
| `NCBI` | `GCA_/GCF_<num>.<ver>` or `GCA_/GCF_…_ASM…_genomic` (5,451 unique assemblies after GCA/GCF twin dedup) | 8,331 | 5,595 |
| `BV-BRC` | BV-BRC `taxid.version` (e.g. `1138383.42`) | 248 | 199 |
| **total** | | **26,499** | **17,117** |

Caveats for downstream code:

- **`genome_id` is NOT uniformly `taxid.version`.** Only the 248 `BV-BRC` rows
  use `taxid.version`; 8,331 rows are NCBI assembly IDs and 17,920 rows are run
  accessions. `accession == genome_id` holds for run-assembly and BV-BRC rows but
  NOT for the `GCA_…_ASM…_genomic`-style NCBI rows (5,102 summary rows), where
  `accession` is the bare `GCA_/GCF_<num>.<ver>`. Do **not** confuse these
  columns with the v2 NCBI manifest semantics — check `source` first.
- **Species labels use BV-BRC type-strain names** (e.g. `Mycobacteroides abscessus
  ATCC 19977`, `… subsp. massiliense CIP 108297`; 17,235 of 26,499 rows carry a
  type-strain/collection token). These are **incompatible with v2's plain
  taxonomy labels** (`Mycobacteroides abscessus`) — cross-referencing to v2 must
  use **contig accessions**, never species names (see `v2_overlap_report.md`).

## 5. Provenance chain

- Source: BV-BRC phigaro QC-passed NTM export, dated **2026-09-09** (filename
  suffix). Export covers the QC-passed NTM cohort with phigaro prophage calls:
  both the per-genome summary and per-prophage coordinates.
- This cohort **expands** PHIND-NTM from the v2 13,122-assembly + 9,543-run
  cohort toward an all-inclusive cohort; see `../README.md` for the v2∪v3 union
  design and `v2_overlap_report.md` for the measured overlap.
