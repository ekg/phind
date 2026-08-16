# NTM v2 — PanSN genome download report

Generated: 2026-08-16T13:27:46.575498Z

## Counts

| metric | count |
|---|---|
| cohort total (unique assemblies) | 13122 |
| downloaded | 11972 |
| linked from v1 | 1150 |
| failed | 0 |
| resolved (downloaded+linked+failed) | 13122 |
| total bp | 63,808,478,470 |
| GCF-only accessions (gcf_only.tsv) | 332 |

## Validation

- [x] resolved == cohort total
- [x] samtools faidx ok on all downloaded
- [x] gzip -t passes on sampled files
- [x] spot-check contig names match manifest

## Failures

None.

## Spot-check (contig names vs manifest prophage_contig)

- GCA_000157895.2: manifest prophage_contig CP006836.1 -> normalized CP006836.1 MATCH (fai has 2 contigs)
- GCA_000157895.2: manifest prophage_contig CP006835.1 -> normalized CP006835.1 MATCH (fai has 2 contigs)
- GCA_000157895.2: manifest prophage_contig CP006835.1 -> normalized CP006835.1 MATCH (fai has 2 contigs)
- GCA_000164135.1: manifest prophage_contig GG770559.1 -> normalized GG770559.1 MATCH (fai has 124 contigs)
- GCA_000164135.1: manifest prophage_contig GG770554.1 -> normalized GG770554.1 MATCH (fai has 124 contigs)
- GCA_000164135.1: manifest prophage_contig GG770561.1 -> normalized GG770561.1 MATCH (fai has 124 contigs)
- GCA_000174035.1: manifest prophage_contig ACFI01000164.1 -> normalized ACFI01000164.1 MATCH (fai has 258 contigs)

## GCF-only accessions (recorded in gcf_only.tsv)

332 accessions — see `gcf_only.tsv`.