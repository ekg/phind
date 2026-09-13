# NTM v3 — unified prophage manifest + full_prophages.fa extraction report

Generated: 2026-09-13T01:39:50Z

## Design

v3 is all-inclusive: BV-BRC phigaro coordinates are primary for genomes
present in the 2026-09-09 BV-BRC phigaro QC-passed export (any row
namespace: ASSEMBLY runs, NCBI assemblies, BV-BRC taxid.version); v2
collaborator calls are kept only for genomes absent from the export
(v2-only NCBI assemblies + the 135 v2-only runs). No genome carries
calls from both callers — on overlap the v2 set is superseded (counted
below, not emitted). Coordinates are used as given; the caller is never
re-run (v2 convention).

## Manifest reconciliation

| metric | count |
|---|---:|
| coordinates CSV rows | 47,994 |
| stub rows (empty prophage_id) in CSV | 22,225 |
| stub rows dropped as exact duplicates of named rows | 5,927 |
| stub rows kept with synthesized `<scaffold>_prophage<N>` id | 16,298 |
| unique prophages after stub dedup | 42,067 |
| export genomes whose summary prophage_count double-counts | 2,549 |
| canonical assemblies with >1 export genome/object | 2,218 (GCA/GCF twin objects: 1,870) |
| export genomes superseded export-internally | 2,252 |
| export rows superseded export-internally (one genome per assembly kept, GCA-preferred) | 5,746 |
| export genomes kept (phigaro calls) | 14,865 |
| v2 manifest rows / has_prophage / kept after GCA-dedup | 51,004 / 30,165 / 23,243 |
| v2 numeric rows superseded (genome in export) | 7,966 (3,060 genomes) |
| v2 run rows superseded (run in export) | 14,658 |
| v2 rows kept (genomes absent from export) | 619 (536 NCBI + 83 runs) |
| **unified manifest rows** | **36,940** |
| — source=BV-BRC (caller=phigaro) | 36,321 |
| — source=V2 (caller=v2-collab) | 619 |

Genomes in both cohorts (reconciliation basis): every v2 kept
numeric assembly whose 9-digit numeric is in the export (3,060 prophage-bearing
genomes) and every v2 run accession present in the export
(14,658) had its v2 prophage set superseded by the phigaro set;
536 v2 NCBI rows and 83 v2 run rows are the v2-only retention.

## Extraction counts

| metric | count |
|---|---:|
| unified manifest rows | 36,940 |
| extractable rows (genome local) | 9,446 (BV-BRC 8,910 + V2 536) |
| extracted into FASTA | 9,446 |
| errors (scaffold/end/length) | 0 |
| blocked rows (run assemblies, collaborator pending) | 27,494 (BV-BRC 27,411 + V2 83) |
| distinct canonical objects read | 3,994 |
| FASTA records written | 9,446 == extractable − errors |

Run-assembly accounting (all-inclusive, none dropped): every v2 run
prophage row is present in the union — 14,658 via the export phigaro
rows (blocked_run: FASTAs pending collaborator delivery, ntm/v2/run_assemblies/
REQUEST.md PENDING) and 83 as source=V2 rows (the 36
prophage-bearing v2-only runs; the other 99 v2-only runs carry no
prophages). ENA substitute assemblies are content-mismatched and are
not used (v2 ena_backfill study).

## Length distribution (extracted prophages)

| stat | bp |
|---|---:|
| min | 237 |
| median | 18,276 |
| mean | 21,838 |
| max | 93,713 |
| total | 206,284,318 |

## Coordinate conventions

Coordinates are 1-based inclusive: `length == end − begin + 1` holds for
all 36,940 manifest rows (asserted at build time), and every extracted
sequence length equals its manifest length — except the begin==0 rows
(16 extractable of 52 manifest-wide, same assemblies and
convention as the v2 extract report): there begin is 0-based while end
stays 1-based inclusive, so the region is extracted as [1, end] and the
extracted length is length − 1. The transposable flag is metadata only
and does not shift coordinates (verified across every row: extracted
length == end − max(1,begin) + 1).

Scaffold→contig resolution methods (normalized both sides:
accn\| stripped, unversioned, NZ_/NW_ prefix added or stripped — both
directions occur: RefSeq objects rename contigs the export lists
unprefixed, and GCA twin objects store NZ_ contigs unprefixed):

| method | rows |
|---|---:|
| exact | 9,240 |
| exact+accn_strip | 107 |
| nz_add | 72 |
| version_suffix+accn_strip | 17 |
| nz_strip | 10 |

## Flagged warnings (kept, not dropped)

Prophages outside [1,000, 100,000] bp: **91** — all extracted and retained:

| prophage | source | length |
|---|---|---:|
| 36809.1962#1#accn|JBOZLN010000004_prophage1 | BV-BRC | 680 |
| 36809.1967#1#accn|JBOZLT010000005_prophage1 | BV-BRC | 680 |
| 36809.1970#1#accn|JBOZLV010000004_prophage1 | BV-BRC | 680 |
| 36809.1974#1#accn|JBOZLX010000007_prophage1 | BV-BRC | 680 |
| 36809.1979#1#accn|JBOZMD010000005_prophage1 | BV-BRC | 680 |
| 36809.1983#1#accn|JBOZMH010000005_prophage1 | BV-BRC | 680 |
| 36809.2012#1#accn|JBOZNK010000004_prophage1 | BV-BRC | 680 |
| 36809.2030#1#accn|JBGKBR010000004_prophage1 | BV-BRC | 680 |
| GCA_001946105.1#1#GCA_001946105.1_ASM194610v1_genomic_prophage1 | V2 | 813 |
| GCA_001948695.1#1#GCA_001948695.1_ASM194869v1_genomic_prophage1 | V2 | 813 |
| GCA_002101775.1#1#LQPI01000027.1_prophage1 | BV-BRC | 516 |
| GCA_002800615.1#1#NQRO01000009.1_prophage1 | BV-BRC | 731 |
| GCA_002800805.1#1#NQSE01000008.1_prophage1 | BV-BRC | 680 |
| GCA_002800815.1#1#NQSF01000008.1_prophage1 | BV-BRC | 680 |
| GCA_002800845.1#1#NQSG01000009.1_prophage1 | BV-BRC | 680 |
| GCA_002800945.1#1#NQSS01000010.1_prophage1 | BV-BRC | 978 |
| GCA_002800965.1#1#NQSU01000009.1_prophage1 | BV-BRC | 680 |
| GCA_002802225.1#1#NQSD01000008.1_prophage1 | BV-BRC | 680 |
| GCA_002802285.1#1#NQSH01000009.1_prophage1 | BV-BRC | 680 |
| GCA_003072295.1#1#QDEW01000008.1_prophage1 | BV-BRC | 731 |
| GCA_003582525.1#1#QXBL01000001.1_prophage1 | BV-BRC | 731 |
| GCA_003582645.1#1#QXBR01000001.1_prophage1 | BV-BRC | 731 |
| GCA_003582765.1#1#QXBX01000004.1_prophage1 | BV-BRC | 680 |
| GCA_003582845.1#1#QXCB01000004.1_prophage1 | BV-BRC | 680 |
| GCA_004105485.1#1#GCA_004105485.1_ASM410548v1_genomic_prophage1 | V2 | 813 |
| GCA_004106125.1#1#GCA_004106125.1_ASM410612v1_genomic_prophage1 | V2 | 813 |
| GCA_015023895.1#1#JACVDN010000001.1_prophage1 | BV-BRC | 680 |
| GCA_015355655.1#1#JACLAQ010000007.1_prophage2 | BV-BRC | 744 |
| GCA_015355675.1#1#JACLAP010000005.1_prophage2 | BV-BRC | 744 |
| GCA_015499795.1#1#JACDRJ010000027.1_prophage1 | BV-BRC | 731 |
| GCA_016756035.1#1#AP024240.1_prophage7 | BV-BRC | 650 |
| GCA_017176265.1#1#JADWXJ010000010.1_prophage1 | BV-BRC | 680 |
| GCA_017183555.1#1#CP063318.1_prophage1 | BV-BRC | 680 |
| GCA_020055525.1#1#WEHM01000002.1_prophage1 | BV-BRC | 641 |
| GCA_028211145.1#1#JAQLTS010000015.1_prophage2 | BV-BRC | 638 |
| GCA_030463595.1#1#CP119724.1_prophage1 | BV-BRC | 680 |
| GCA_030513445.1#1#JAROLG010000008.1_prophage1 | BV-BRC | 680 |
| GCA_030513785.1#1#JAROLS010000006.1_prophage1 | BV-BRC | 731 |
| GCA_049066265.1#1#JBMEVX010000015.1_prophage1 | BV-BRC | 731 |
| GCA_049066725.1#1#JBMEWK010000008.1_prophage1 | BV-BRC | 731 |
| GCA_900133475.1#1#FVHJ01000008.1_prophage1 | BV-BRC | 853 |
| GCA_900133695.1#1#FVHZ01000002.1_prophage1 | BV-BRC | 853 |
| GCA_900133725.1#1#FVHQ01000008.1_prophage1 | BV-BRC | 853 |
| GCA_900133745.1#1#FVHX01000001.1_prophage1 | BV-BRC | 853 |
| GCA_900133755.1#1#FVIG01000008.1_prophage1 | BV-BRC | 853 |
| GCA_900135145.1#1#FSGF01000001.1_prophage1 | BV-BRC | 731 |
| GCA_900135205.1#1#FVLA01000022.1_prophage1 | BV-BRC | 731 |
| GCA_900136385.1#1#FSIR01000005.1_prophage1 | BV-BRC | 680 |
| GCA_900136595.1#1#FVNL01000004.1_prophage1 | BV-BRC | 731 |
| GCA_900136615.1#1#FVNJ01000005.1_prophage1 | BV-BRC | 731 |
| … 41 more (see manifest; filter length column) | | |

end > contig length rows: **0**

None — every prophage fits inside its resolved contig.

## begin==0 rows (extracted as [1, end], extracted length = length − 1)

| header | genome | scaffold | end | length |
|---|---|---|---:|---:|
| GCA_000270825.1#1#AKUX01000007.1_prophage1 | GCA_000270825.1 | AKUX01000007.1 | 28,307 | 28,308 |
| GCA_001545925.1#1#GCA_001545925.1_ASM154592v1_genomic_prophage1 | GCA_001545925.1_ASM154592v1_genomic | LMVQ01000168.1 | 32,509 | 32,510 |
| GCA_030330585.1#1#JAQPLK010000002.1_prophage1 | GCA_030330585.1_ASM3033058v1_genomic | JAQPLK010000002.1 | 15,835 | 15,836 |
| GCA_030330595.1#1#JAQPLL010000005.1_prophage1 | GCA_030330595.1_ASM3033059v1_genomic | JAQPLL010000005.1 | 15,835 | 15,836 |
| GCA_049066535.1#1#GCA_049066535.1_ASM4906653v1_genomic_prophage2 | GCA_049066535.1_ASM4906653v1_genomic | JBMEWF010000012.1 | 31,976 | 31,977 |
| GCA_900134085.1#1#FVIT01000018.1_prophage1 | GCA_900134085.1 | FVIT01000018.1 | 20,171 | 20,172 |
| GCA_900134115.1#1#FVIU01000017.1_prophage1 | GCA_900134115.1 | FVIU01000017.1 | 21,022 | 21,023 |
| GCA_900139745.1#1#FVUW01000021.1_prophage1 | GCA_900139745.1 | FVUW01000021.1 | 24,419 | 24,420 |
| GCA_900960235.1#1#CAAHFM010000055.1_prophage1 | GCA_900960235.1 | CAAHFM010000055.1 | 13,200 | 13,201 |
| GCF_000271105.1#1#NZ_AKUP01000011.1_prophage1 | GCF_000271105.1 | NZ_AKUP01000011.1 | 28,400 | 28,401 |
| GCF_000987455.1#1#NZ_LBEW01000009.1_prophage1 | GCF_000987455.1_ASM98745v1_genomic | NZ_LBEW01000009.1 | 28,674 | 28,675 |
| GCF_047739115.1#1#NZ_JBLLUQ010000023.1_prophage1 | GCF_047739115.1_ASM4773911v1_genomic | NZ_JBLLUQ010000023.1 | 28,111 | 28,112 |
| GCF_054074565.1#1#NZ_JBOZMN010000019.1_prophage1 | GCF_054074565.1_ASM5407456v1_genomic | NZ_JBOZMN010000019.1 | 35,441 | 35,442 |
| GCF_057397965.1#1#NZ_JBXXZF010000018.1_prophage1 | GCF_057397965.1_ASM5739796v1_genomic | NZ_JBXXZF010000018.1 | 27,824 | 27,825 |
| GCF_900132295.1#1#NZ_FSBM01000004.1_prophage1 | GCF_900132295.1_12163_2_73_genomic | NZ_FSBM01000004.1 | 28,105 | 28,106 |
| GCF_900133575.1#1#NZ_FSDN01000012.1_prophage1 | GCF_900133575.1_12082_5_56_genomic | NZ_FSDN01000012.1 | 29,572 | 29,573 |

## Independent spot check (10 genomes, seed 42)

Each row re-queried via a fresh `samtools faidx` call and byte-compared
to the record parsed back out of the written FASTA:

| header | scaffold | contig | begin | end | len | result |
|---|---|---|---:|---:|---:|---|
| GCA_905219475.1#1#CAJMWI010000001.1_prophage1 | CAJMWI010000001.1 | GCA_905219475.1#1#CAJMWI010000001.1 | 4,098,292 | 4,110,816 | 12,525 | PASS |
| GCA_002304075.1#1#MBFN01000161.1_prophage1 | MBFN01000161.1 | GCA_002304075.1#1#MBFN01000161.1 | 2,880 | 9,186 | 6,307 | PASS |
| GCA_000270865.1#1#AKUF01000004.1_prophage1 | AKUF01000004.1 | GCA_000270865.1#1#AKUF01000004.1 | 14,902 | 20,346 | 5,445 | PASS |
| GCF_002013725.1#1#NZ_MAER01000011.1_prophage1 | NZ_MAER01000011.1 | GCF_002013725.1#1#NZ_MAER01000011.1 | 3,482 | 19,500 | 16,019 | PASS |
| GCA_021559935.1#1#CP060055.1_prophage2 | CP060055.1 | GCA_021559935.1#1#CP060055.1 | 1,321,851 | 1,340,041 | 18,191 | PASS |
| GCA_017176255.1#1#JADWXG010000002.1_prophage1 | JADWXG010000002.1 | GCA_017176255.1#1#JADWXG010000002.1 | 489,260 | 499,228 | 9,969 | PASS |
| GCA_015499625.1#1#JACDRB010000004.1_prophage2 | JACDRB010000004.1 | GCA_015499625.1#1#JACDRB010000004.1 | 404,531 | 416,506 | 11,976 | PASS |
| GCA_002802425.1#1#NQSX01000002.1_prophage1 | NQSX01000002.1 | GCA_002802425.1#1#NQSX01000002.1 | 525,177 | 528,528 | 3,352 | PASS |
| GCF_001954125.1#1#NZ_MBGE01000274.1_prophage1 | NZ_MBGE01000274.1 | GCF_001954125.1#1#NZ_MBGE01000274.1 | 46,382 | 53,491 | 7,110 | PASS |
| GCA_002102225.1#1#LQPB01000023.1_prophage1 | LQPB01000023.1 | GCA_002102225.1#1#LQPB01000023.1 | 3,883 | 33,112 | 29,230 | PASS |

## Output

- FASTA: `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/full_prophages.fa` — 9,446 records, 206,284,318 bp, sha256 `5adfd2aee6e1e03d109a49a3b1fbba837a21649851fc39ee9c9ff07217f0480c`
- per-record status: `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/full_prophages.fa.manifest.tsv`
- manifest (repo): `ntm/v3/inputs/v3_prophage_manifest.tsv.gz` (36,940 rows; plain copy on NVMe at `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/inputs/v3_prophage_manifest.tsv`)

FASTA stays on NVMe (repo holds only code, manifests and small
reports — v1/v2 rule).

