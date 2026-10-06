# NTM v3 — unified prophage manifest + full_prophages.fa extraction report

Generated: 2026-10-06T16:51:28Z

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
| extractable rows (genome local) | 36,857 (BV-BRC 36,321 + V2 536) |
| extracted into FASTA | 36,857 |
| errors (scaffold/end/length) | 0 |
| blocked rows (run assemblies, collaborator pending) | 83 (BV-BRC 0 + V2 83) |
| distinct canonical objects read | 15,317 |
| FASTA records written | 36,857 == extractable − errors |

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
| median | 18,226 |
| mean | 21,315 |
| max | 100,585 |
| total | 785,633,429 |

## Coordinate conventions

Coordinates are 1-based inclusive: `length == end − begin + 1` holds for
all 36,940 manifest rows (asserted at build time), and every extracted
sequence length equals its manifest length — except the begin==0 rows
(46 extractable of 52 manifest-wide, same assemblies and
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
| exact | 36,651 |
| exact+accn_strip | 107 |
| nz_add | 72 |
| version_suffix+accn_strip | 17 |
| nz_strip | 10 |

## Flagged warnings (kept, not dropped)

Prophages outside [1,000, 100,000] bp: **263** — all extracted and retained:

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
| DRR317567#1#NODE_6_length_321285_cov_65.078585_prophage1 | BV-BRC | 100,585 |
| DRR444226#1#NODE_4_length_230282_cov_15.436384_prophage1 | BV-BRC | 537 |
| ERR115000#1#NODE_15_length_117014_cov_9.211228_prophage1 | BV-BRC | 731 |
| ERR115004#1#NODE_24_length_79134_cov_19.873835_prophage1 | BV-BRC | 731 |
| ERR115005#1#NODE_14_length_164903_cov_20.224534_prophage1 | BV-BRC | 731 |
| ERR115031#1#NODE_4_length_372540_cov_25.784955_prophage2 | BV-BRC | 237 |
| ERR115032#1#NODE_4_length_372540_cov_24.629134_prophage1 | BV-BRC | 237 |
| ERR115033#1#NODE_4_length_372540_cov_27.174960_prophage2 | BV-BRC | 237 |
| ERR115040#1#NODE_24_length_79133_cov_25.091315_prophage1 | BV-BRC | 731 |
| ERR115046#1#NODE_22_length_79118_cov_22.203673_prophage1 | BV-BRC | 731 |
| ERR115047#1#NODE_4_length_338841_cov_24.056573_prophage2 | BV-BRC | 237 |
| ERR115079#1#NODE_24_length_79124_cov_29.069180_prophage1 | BV-BRC | 731 |
| ERR119103#1#NODE_19_length_77853_cov_38.103525_prophage1 | BV-BRC | 680 |
| ERR119107#1#NODE_14_length_164826_cov_37.801773_prophage1 | BV-BRC | 731 |
| ERR13148610#1#NODE_51_length_30047_cov_59.265432_prophage1 | BV-BRC | 554 |
| ERR16089468#1#NODE_9_length_234012_cov_173.900524_prophage1 | BV-BRC | 680 |
| ERR16089539#1#NODE_9_length_234012_cov_214.094885_prophage1 | BV-BRC | 680 |
| ERR16089633#1#NODE_9_length_234012_cov_238.695963_prophage1 | BV-BRC | 680 |
| ERR16914737#1#NODE_18_length_120521_cov_63.284448_prophage1 | BV-BRC | 883 |
| ERR2524309#1#NODE_16_length_140867_cov_64.265659_prophage1 | BV-BRC | 731 |
| ERR2524342#1#NODE_12_length_163579_cov_66.657292_prophage1 | BV-BRC | 731 |
| ERR2524354#1#NODE_3_length_687324_cov_85.222345_prophage2 | BV-BRC | 731 |
| ERR2759453#1#NODE_1_length_787494_cov_26.280716_prophage1 | BV-BRC | 731 |
| ERR2759457#1#NODE_1_length_1934245_cov_24.731971_prophage1 | BV-BRC | 731 |
| ERR2759471#1#NODE_4_length_636007_cov_24.246119_prophage1 | BV-BRC | 680 |
| ERR2759478#1#NODE_3_length_525046_cov_28.504138_prophage1 | BV-BRC | 680 |
| ERR2759479#1#NODE_4_length_525046_cov_28.251785_prophage1 | BV-BRC | 680 |
| ERR2759480#1#NODE_6_length_402313_cov_23.476363_prophage1 | BV-BRC | 680 |
| ERR2759483#1#NODE_4_length_449205_cov_25.609736_prophage1 | BV-BRC | 680 |
| ERR2759484#1#NODE_3_length_449196_cov_26.050283_prophage1 | BV-BRC | 680 |
| ERR2759485#1#NODE_4_length_449196_cov_28.750766_prophage1 | BV-BRC | 680 |
| ERR2759487#1#NODE_5_length_449205_cov_22.316992_prophage1 | BV-BRC | 680 |
| ERR2759488#1#NODE_3_length_525046_cov_26.168073_prophage1 | BV-BRC | 680 |
| ERR2759489#1#NODE_6_length_355569_cov_26.791237_prophage1 | BV-BRC | 680 |
| ERR2759492#1#NODE_3_length_449196_cov_27.255298_prophage1 | BV-BRC | 680 |
| ERR2759494#1#NODE_6_length_402363_cov_25.181540_prophage1 | BV-BRC | 680 |
| ERR2759495#1#NODE_4_length_449205_cov_25.920053_prophage1 | BV-BRC | 680 |
| ERR2759496#1#NODE_4_length_449205_cov_24.488883_prophage1 | BV-BRC | 680 |
| ERR2759497#1#NODE_5_length_449205_cov_26.911192_prophage1 | BV-BRC | 680 |
| ERR2759500#1#NODE_4_length_422143_cov_26.497761_prophage1 | BV-BRC | 680 |
| ERR2759522#1#NODE_3_length_525046_cov_25.635079_prophage1 | BV-BRC | 680 |
| ERR3198401#1#NODE_3_length_515028_cov_12.823854_prophage1 | BV-BRC | 731 |
| … 213 more (see manifest; filter length column) | | |

end > contig length rows: **0**

None — every prophage fits inside its resolved contig.

## begin==0 rows (extracted as [1, end], extracted length = length − 1)

| header | genome | scaffold | end | length |
|---|---|---|---:|---:|
| ERR1413183#1#NODE_29_length_53278_cov_57.370012_prophage1 | ERR1413183 | NODE_29_length_53278_cov_57.370012 | 26,489 | 26,490 |
| ERR15535446#1#NODE_31_length_32374_cov_23.438740_prophage1 | ERR15535446 | NODE_31_length_32374_cov_23.438740 | 22,031 | 22,032 |
| ERR3012659#1#NODE_3_length_182958_cov_24.190121_prophage1 | ERR3012659 | NODE_3_length_182958_cov_24.190121 | 23,101 | 23,102 |
| ERR330893#1#NODE_67_length_21375_cov_4.526220_prophage1 | ERR330893 | NODE_67_length_21375_cov_4.526220 | 14,937 | 14,938 |
| ERR337788#1#NODE_23_length_24386_cov_35.715384_prophage1 | ERR337788 | NODE_23_length_24386_cov_35.715384 | 23,808 | 23,809 |
| ERR484975#1#NODE_27_length_26202_cov_713.379776_prophage1 | ERR484975 | NODE_27_length_26202_cov_713.379776 | 16,828 | 16,829 |
| ERR5412607#1#NODE_37_length_50832_cov_23.344892_prophage1 | ERR5412607 | NODE_37_length_50832_cov_23.344892 | 21,878 | 21,879 |
| ERR5412799#1#NODE_29_length_42306_cov_7.886523_prophage1 | ERR5412799 | NODE_29_length_42306_cov_7.886523 | 40,257 | 40,258 |
| ERR5413109#1#NODE_71_length_21642_cov_17.558345_prophage1 | ERR5413109 | NODE_71_length_21642_cov_17.558345 | 19,701 | 19,702 |
| ERR5413703#1#NODE_36_length_45043_cov_33.375834_prophage1 | ERR5413703 | NODE_36_length_45043_cov_33.375834 | 38,240 | 38,241 |
| ERR5413755#1#NODE_45_length_37154_cov_28.307178_prophage1 | ERR5413755 | NODE_45_length_37154_cov_28.307178 | 21,538 | 21,539 |
| ERR5414153#1#NODE_38_length_50696_cov_11.219802_prophage1 | ERR5414153 | NODE_38_length_50696_cov_11.219802 | 24,750 | 24,751 |
| ERR7253671#1#NODE_23_length_95955_cov_39.359071_prophage1 | ERR7253671 | NODE_23_length_95955_cov_39.359071 | 20,507 | 20,508 |
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
| SRR14719107#1#NODE_1_length_802290_cov_100.357398_prophage1 | SRR14719107 | NODE_1_length_802290_cov_100.357398 | 27,571 | 27,572 |
| SRR14719143#1#NODE_25_length_21996_cov_227.234865_prophage1 | SRR14719143 | NODE_25_length_21996_cov_227.234865 | 13,115 | 13,116 |
| SRR14719277#1#NODE_20_length_43016_cov_1427.344652_prophage1 | SRR14719277 | NODE_20_length_43016_cov_1427.344652 | 40,005 | 40,006 |
| SRR17774879#1#NODE_27_length_68999_cov_14.921478_prophage1 | SRR17774879 | NODE_27_length_68999_cov_14.921478 | 37,660 | 37,661 |
| SRR19543348#1#NODE_71_length_23702_cov_418.039704_prophage1 | SRR19543348 | NODE_71_length_23702_cov_418.039704 | 19,175 | 19,176 |
| SRR21939558#1#NODE_4_length_392523_cov_25.262316_prophage1 | SRR21939558 | NODE_4_length_392523_cov_25.262316 | 27,070 | 27,071 |
| SRR21939562#1#NODE_27_length_50832_cov_25.472485_prophage1 | SRR21939562 | NODE_27_length_50832_cov_25.472485 | 21,878 | 21,879 |
| SRR21939572#1#NODE_14_length_102707_cov_205.917802_prophage1 | SRR21939572 | NODE_14_length_102707_cov_205.917802 | 27,070 | 27,071 |
| SRR21939590#1#NODE_2_length_392528_cov_26.212032_prophage1 | SRR21939590 | NODE_2_length_392528_cov_26.212032 | 27,070 | 27,071 |
| SRR32740559#1#NODE_5_length_292346_cov_257.915191_prophage1 | SRR32740559 | NODE_5_length_292346_cov_257.915191 | 24,401 | 24,402 |
| SRR3394782#1#NODE_33_length_54486_cov_106.615486_prophage1 | SRR3394782 | NODE_33_length_54486_cov_106.615486 | 51,257 | 51,258 |
| SRR6045021#1#NODE_20_length_89275_cov_12.265757_prophage1 | SRR6045021 | NODE_20_length_89275_cov_12.265757 | 14,001 | 14,002 |
| SRR6045402#1#NODE_58_length_21880_cov_55.905013_prophage1 | SRR6045402 | NODE_58_length_21880_cov_55.905013 | 13,979 | 13,980 |
| SRR6046714#1#NODE_43_length_21279_cov_30.978445_prophage1 | SRR6046714 | NODE_43_length_21279_cov_30.978445 | 20,150 | 20,151 |
| SRR6046890#1#NODE_37_length_42337_cov_7.427094_prophage1 | SRR6046890 | NODE_37_length_42337_cov_7.427094 | 35,777 | 35,778 |
| SRR6871408#1#NODE_21_length_91858_cov_29.861039_prophage1 | SRR6871408 | NODE_21_length_91858_cov_29.861039 | 27,070 | 27,071 |
| SRR8291522#1#NODE_53_length_27688_cov_15.311072_prophage1 | SRR8291522 | NODE_53_length_27688_cov_15.311072 | 25,851 | 25,852 |

## Independent spot check (10 genomes, seed 42)

Each row re-queried via a fresh `samtools faidx` call and byte-compared
to the record parsed back out of the written FASTA:

| header | scaffold | contig | begin | end | len | result |
|---|---|---|---:|---:|---:|---|
| SRR11090511#1#NODE_62_length_36083_cov_32.100817_prophage1 | NODE_62_length_36083_cov_32.100817 | SRR11090511#1#NODE_62_length_36083_cov_32.100817 | 803 | 7,907 | 7,105 | PASS |
| ERR16914762#1#NODE_22_length_83412_cov_33.337625_prophage1 | NODE_22_length_83412_cov_33.337625 | ERR16914762#1#NODE_22_length_83412_cov_33.337625 | 41,384 | 47,429 | 6,046 | PASS |
| DRR317494#1#NODE_1_length_532610_cov_73.703132_prophage1 | NODE_1_length_532610_cov_73.703132 | DRR317494#1#NODE_1_length_532610_cov_73.703132 | 139,198 | 149,807 | 10,610 | PASS |
| SRR24581785#1#NODE_21_length_109192_cov_123.982862_prophage1 | NODE_21_length_109192_cov_123.982862 | SRR24581785#1#NODE_21_length_109192_cov_123.982862 | 58,719 | 64,580 | 5,862 | PASS |
| ERR484990#1#NODE_14_length_148115_cov_26.614021_prophage1 | NODE_14_length_148115_cov_26.614021 | ERR484990#1#NODE_14_length_148115_cov_26.614021 | 39,479 | 67,130 | 27,652 | PASS |
| ERR4022328#1#NODE_19_length_79838_cov_18.375760_prophage1 | NODE_19_length_79838_cov_18.375760 | ERR4022328#1#NODE_19_length_79838_cov_18.375760 | 64,999 | 73,680 | 8,682 | PASS |
| ERR369324#1#NODE_2_length_948644_cov_27.668791_prophage1 | NODE_2_length_948644_cov_27.668791 | ERR369324#1#NODE_2_length_948644_cov_27.668791 | 440,404 | 483,184 | 42,781 | PASS |
| ERR3142121#1#NODE_1_length_1750669_cov_21.989051_prophage1 | NODE_1_length_1750669_cov_21.989051 | ERR3142121#1#NODE_1_length_1750669_cov_21.989051 | 498,818 | 508,786 | 9,969 | PASS |
| SRR22333106#1#NODE_1_length_222142_cov_111.984892_prophage1 | NODE_1_length_222142_cov_111.984892 | SRR22333106#1#NODE_1_length_222142_cov_111.984892 | 159,262 | 202,924 | 43,663 | PASS |
| ERR15501063#1#NODE_34_length_67779_cov_40.821822_prophage1 | NODE_34_length_67779_cov_40.821822 | ERR15501063#1#NODE_34_length_67779_cov_40.821822 | 35,985 | 42,030 | 6,046 | PASS |

## Output

- FASTA: `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/full_prophages.fa` — 36,857 records, 785,633,429 bp, sha256 `a949c653b104521f9d52aa95af4a233a3e35418c10cf3c782e26a6e3317441bf`
- per-record status: `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/full_prophages.fa.manifest.tsv`
- manifest (repo): `ntm/v3/inputs/v3_prophage_manifest.tsv.gz` (36,940 rows; plain copy on NVMe at `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/inputs/v3_prophage_manifest.tsv`)

FASTA stays on NVMe (repo holds only code, manifests and small
reports — v1/v2 rule).
