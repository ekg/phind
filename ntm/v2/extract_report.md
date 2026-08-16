# NTM v2 — prophage FASTA extraction report

Generated: 2026-08-16T13:42:05Z

## Counts

| metric | count |
|---|---|
| manifest rows (excl. header) | 51004 |
| has_prophage=True rows | 30165 |
| kept after GCA-dedup | 23243 |
| NCBI prophage rows (extract) | 8502 |
| run-assembly rows (skipped, pending v2.1) | 14741 |
| NCBI extracted into FASTA | 8502 |
| NCBI missing (no local fasta / contig) | 0 |
| NCBI error (faidx / mismatch) | 0 |

## Length distribution (extracted NCBI prophages)

| stat | bp |
|---|---|
| min | 237 |
| median | 18725 |
| mean | 21952 |
| max | 93713 |
| total | 186637140 |

Coordinates are 1-based inclusive (end-start+1 == length_bp, verified in ntm-v2-acquire). Extraction used samtools faidx with the raw start/end (off-by-one delta = 0) except for the start=0 rows below, where start is 0-based (length_bp == end - start + 1 still holds), so the region was extracted as [1, end] and the extracted length is length_bp - 1 (deliberate off-by-one). Every non-zero-start extracted sequence length equals its prophage_length_bp. Columns checked: 8502 unique (acc,contig,start,end) keys; duplicate same-coord rows: 0.

prophage_id duplicates among NCBI rows: 0.

Rows with 0-based start (start==0, extracted as [1, end], len == length_bp - 1): 13

| prophage_id | contig | end | length_bp | extracted |
|---|---|---|---|---|
| GCA_900960235.1_PRJEB31972-3_genomic_prophage2 | CAAHFM010000055.1 | 13200 | 13201 | 13200 |
| GCA_000270825.1_ASM27082v1_genomic_prophage5 | AKUX01000007.1 | 28307 | 28308 | 28307 |
| GCA_030330585.1_ASM3033058v1_genomic_prophage1 | JAQPLK010000002.1 | 15835 | 15836 | 15835 |
| GCA_030330595.1_ASM3033059v1_genomic_prophage3 | JAQPLL010000005.1 | 15835 | 15836 | 15835 |
| GCA_049066535.1_ASM4906653v1_genomic_prophage2 | JBMEWF010000012.1 | 31976 | 31977 | 31976 |
| GCA_047739115.1_ASM4773911v1_genomic_prophage1 | JBLLUQ010000023.1 | 28111 | 28112 | 28111 |
| GCA_900132295.1_12163_2_73_genomic_prophage2 | FSBM01000004.1 | 28105 | 28106 | 28105 |
| GCA_900133575.1_12082_5_56_genomic_prophage1 | FSDN01000012.1 | 29572 | 29573 | 29572 |
| GCA_900134085.1_10625_4_52_genomic_prophage4 | FVIT01000018.1 | 20171 | 20172 | 20171 |
| GCA_900134115.1_10625_4_53_genomic_prophage1 | FVIU01000017.1 | 21022 | 21023 | 21022 |
| GCA_900139745.1_10660_1_19_genomic_prophage1 | FVUW01000021.1 | 24419 | 24420 | 24419 |
| GCA_000987455.1_ASM98745v1_genomic_prophage1 | LBEW01000009.1 | 28674 | 28675 | 28674 |
| GCA_001545925.1_ASM154592v1_genomic_prophage1 | LMVQ01000168.1 | 32509 | 32510 | 32509 |

## Missing / error rows

None — every NCBI prophage extracted.

## Run-assembly skip

14741 run-assembly prophage rows (ERR/SRR/DRR) were skipped because their FASTAs are not yet local (pending the collaborator shipping them to ntm-v2-run; they arrive in a v2.1 extension).
