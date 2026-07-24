# Prophage table distribution and interval audit

## Scope, provenance, and reproducibility

This is a local, offline audit of the two root inputs only. It did **not** inspect mounted storage, download or resolve sequence, rename identifiers, assess IMPG, build/query a pangenome, or compute a tree. Both inputs are opened read-only. `reproduce.sh` checks their expected SHA-256 values before analysis and verifies that the pre/post hashes match.

| Input | Bytes | Physical lines | SHA-256 |
|---|---:|---:|---|
| `26k_prophage1.csv` | 12,393,209 | 132,405 | `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996` |
| `26k_ecoli_accession.txt` | 417,239 | 26,078 | `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5` |

Regenerate from any working directory with `artifacts/prophage_summary/reproduce.sh`. The implementation uses Python's standard library only, makes no network calls, writes LF-terminated TSVs, fixes ordering explicitly, and places no timestamp in outputs. Quantiles use linear interpolation at `(n-1)p` (R type 7 / NumPy `linear`).

## Headline result and the meaning of “tagged”

The table has no column named `tag`, `tagged`, `status`, `quality`, or `completeness`. Therefore there is **no uniquely evidenced generic ‘tagged prophage’ filter**. Rather than silently choosing one, this audit carries two competing, exactly reproducible interpretations alongside all records:

1. **Flag-positive:** `Decimal(transposable)` is the integer `1` (the observed raw value is `1.0`). This is the only explicit binary, status-like field. The column name supports ‘transposable’, but no supplied documentation establishes a broader meaning such as quality or completeness.
2. **Taxonomy-assigned:** trimmed, case-folded `taxonomy` is neither empty nor `Unknown`. This is defensible only if ‘tagged’ meant taxonomically labeled; it is not equivalent to the binary flag.

All rows have a nonempty `prophage_id`, so ‘has an ID’ selects all records and is not a distinct subset. No completeness/quality stratification is possible from these inputs.

All genome fractions below use **26,077 normalized unique denominator genomes**; locus/record fractions state their own denominator.

| Scope (exact rule) | Source records | Unique loci | Genomes zero | Genomes one | Genomes multiple | Total bp (`end-begin+1`) | Median locus bp |
|---|---:|---:|---:|---:|---:|---:|---:|
| All records | 132,404 / 132,404 | 132,404 / 132,404 | 0 / 26,077 (0.00%) | 1,200 / 26,077 (4.60%) | 24,877 / 26,077 (95.40%) | 3,262,327,228 / 132,404 valid loci | 21,817.50 / 132,404 valid loci |
| Flag-positive (transposable = 1.0) | 7,695 / 132,404 | 7,695 / 7,695 | 20,435 / 26,077 (78.36%) | 4,446 / 26,077 (17.05%) | 1,196 / 26,077 (4.59%) | 230,813,644 / 7,695 valid loci | 25,364 / 7,695 valid loci |
| Taxonomy-assigned (taxonomy != Unknown) | 115,442 / 132,404 | 115,442 / 115,442 | 171 / 26,077 (0.66%) | 1,790 / 26,077 (6.86%) | 24,116 / 26,077 (92.48%) | 3,092,341,893 / 115,442 valid loci | 25,575 / 115,442 valid loci |

The accession denominator is coextensive with the CSV genome keys: every normalized denominator genome has at least one all-record locus. Thus the all-record zero count is genuinely 0 / denominator, but the list cannot estimate prevalence in an independent collection. The subset zero counts remain informative for the two explicit filters.

## Row, parsing, normalization, and denominator reconciliation

| Stage | Result | Denominator / rule |
|---|---:|---|
| CSV physical lines | 132,405 | file bytes split into lines |
| CSV records | 132,405 | 1 header + 132,404 data records |
| Seven-field parse failures | 0 | 132,404 source data records |
| Normalization failures | 0 | 132,404 source data records |
| Normalized records | 132,404 | 132,404 source data records |
| Valid coordinates | 132,404 | 132,404 normalized records |
| Unique locus keys | 132,404 | 132,404 normalized records |
| Duplicate-locus groups / extra rows | 0 / 0 | key `(genome, scaffold, begin, end)` |
| Exact duplicate-record groups | 0 | all seven source fields |
| Duplicate `prophage_id` groups | 0 | 132,404 unique IDs |

The accession file has:

- **26,078 physical / 26,078 nonblank / 26,078 exact-unique tokens**; duplicate extra lines: 0.
- **26,077 normalized unique assembly accessions**, defined strictly by `^GC[AF]_[0-9]+\.[0-9]+$`. `per_genome.tsv` has exactly this many data rows, in accession-file order.
- **1 invalid token:** `genome` at line 26078. It is not a genome and is excluded from the normalized genome denominator, but is retained in reconciliation metrics.
- Exact-key CSV join failures: **0 / 26,077 unique CSV genome keys; 0 / 132,404 records**. Versionless join failures: **0 / 26,077 keys; 0 / 132,404 records**.
- Denominator genomes without an exact CSV hit: **0 / 26,077**.

A versionless key is used for diagnostics only, never to merge rows into the per-genome table. Stripping assembly versions creates **0 collision groups / 26,077 denominator base keys** and **0 / 26,077 CSV base keys**. Stripping scaffold versions creates **0 / 105,793 base-key collision groups**. These are current-data observations, not permission to drop versions.

## Schema and semantic evidence

The header is exactly `end,genome,scaffold,begin,transposable,taxonomy,prophage_id` (`26k_prophage1.csv`, line 1). The first data row demonstrates float-formatted integer coordinates/flag and versioned keys (`26k_prophage1.csv`, line 2).

| Column | Observed values / constraints | Supported interpretation | Unsupported or unknown |
|---|---|---|---|
| `end` | nonmissing; integer-like decimals; range 1,366–7,512,837 | interval end coordinate | indexing convention not explicitly declared |
| `genome` | 26,077 exact values; all match `GCF_...version` | RefSeq assembly-style source key and exact denominator join key | organism quality/completeness |
| `scaffold` | 105,793 exact values; prefixes `NC_` 447, `NZ_` 131,957 records | versioned sequence/contig source key | chromosome vs plasmid vs unplaced status cannot be inferred from prefix |
| `begin` | nonmissing; integer-like decimals; range 1–7,479,841; minimum 1 | interval begin coordinate | indexing convention not explicitly declared |
| `transposable` | `0.0` 124,709/132,404 (94.19%), `1.0` 7,695/132,404 (5.81%) | binary flag named ‘transposable’ | method, evidence threshold, quality, completeness, and meaning of 0 beyond not flag-positive |
| `taxonomy` | 21 exact strings; `Unknown` 16,962/132,404 | class-like source label; exact strings preserved | classifier/method/confidence; mixed labels are not resolved |
| `prophage_id` | 132,404 unique; every value matches `<genome>_prophage_<integer>` and prefix agrees with `genome` | source record identifier | biological stable ID across releases |

The source contains no explicit tag/status/quality/completeness/chromosome-role column. `transposable` is treated as flag-like and `taxonomy` as class-like only because of their names and observed values; unsupported semantics remain unknown.

## Coordinates, locus definition, overlaps, and coordinate QC

All begin/end values parse as finite integers despite `.0` formatting; all have `begin >= 1` and `end >= begin`. The minimum begin of 1 is consistent with **1-based closed** coordinates, so `end-begin+1` is the primary bp sensitivity. This is not conclusive without producer metadata or sequence; consequently every per-genome row and summary also reports `end-begin`. The two totals differ by exactly one bp per valid unique locus. No contig lengths are supplied, so right-edge clipping and out-of-bounds coordinates cannot be assessed.

A locus is the exact normalized tuple **(`genome`, `scaffold`, integer `begin`, integer `end`)**. Deduplication keeps the first source row as representative; subset selection is applied to rows first and then the same locus key is deduplicated. Taxonomy and flag are deliberately absent from the key. Overlap is assessed only within exact `(genome, scaffold)`, on valid unique loci, as closed-interval intersection (`next.begin <= prior.end`). Nesting includes full containment; exact duplicates are handled before overlap.

| Scope | Overlap pairs | Participating loci | Pairwise overlap bp | Nested pairs | Participating nested loci |
|---|---:|---:|---:|---:|---:|
| All records | 0 / 132,404 valid loci | 0 / 132,404 | 0 / 0 pairs | 0 / 132,404 valid loci | 0 / 132,404 |
| Flag-positive (transposable = 1.0) | 0 / 7,695 valid loci | 0 / 7,695 | 0 / 0 pairs | 0 / 7,695 valid loci | 0 / 7,695 |
| Taxonomy-assigned (taxonomy != Unknown) | 0 / 115,442 valid loci | 0 / 115,442 | 0 / 0 pairs | 0 / 115,442 valid loci | 0 / 115,442 |

Diagnostic, not biological-quality, flags:

| Scope | begin <= 3 | length <1,000 bp | length >200,000 bp | Denominator |
|---|---:|---:|---:|---:|
| All records | 12,279 (9.27%) | 817 (0.62%) | 1 (0.00%) | 132,404 valid unique loci |
| Flag-positive (transposable = 1.0) | 365 (4.74%) | 0 (0.00%) | 0 (0.00%) | 7,695 valid unique loci |
| Taxonomy-assigned (taxonomy != Unknown) | 11,897 (10.31%) | 120 (0.10%) | 1 (0.00%) | 115,442 valid unique loci |

`begin <= 3` may indicate a left-contig-edge/truncated call; short/long cutoffs are transparent review thresholds, not evidence-based quality labels. `interval_qc.tsv` identifies every flagged source row.

## Count, bp, and length distributions

Per-genome distributions include all normalized denominator genomes, including zeros. Bp uses valid unique loci and the primary `end-begin+1` sensitivity.

| Scope | Loci/genome median / p95 / max | bp/genome median / p95 / max | Locus bp p05 / median / p95 |
|---|---:|---:|---:|
| All records | 5 / 10 / 25 (N=26,077 genomes) | 106,467 / 275,634.00 / 903,677 (N=26,077) | 4,434 / 21,817.50 / 47,813 (N=132,404 loci) |
| Flag-positive (transposable = 1.0) | 0 / 1 / 11 (N=26,077 genomes) | 0 / 48,392.20 / 442,805 (N=26,077) | 6,787 / 25,364 / 67,666.40 (N=7,695 loci) |
| Taxonomy-assigned (taxonomy != Unknown) | 4 / 9 / 23 (N=26,077 genomes) | 99,992 / 268,494.20 / 892,706 (N=26,077) | 6,192.20 / 25,575 / 48,860.95 (N=115,442 loci) |

`summary_metrics.tsv` adds min, p01, p05, p25, median, p75, p95, p99, max, mean, and population SD for per-genome locus counts, hit-scaffold counts, both bp conventions, and locus lengths for all three scopes.

### Extreme genomes

Top genomes are descriptive extremes, not automatically errors. Every rank is among the same N=26,077 denominator genomes.

| Scope | Rank | Genome exact key | Unique loci | Total inclusive bp |
|---|---:|---|---:|---:|
| All records | 1 | `GCF_000194335.1` | 25 | 402,376 |
| All records | 2 | `GCF_001039215.2` | 24 | 571,418 |
| All records | 3 | `GCF_014217035.1` | 23 | 614,708 |
| All records | 4 | `GCF_000703225.1` | 23 | 560,341 |
| All records | 5 | `GCF_002863685.1` | 22 | 726,011 |
| All records | 6 | `GCF_001039155.2` | 22 | 584,066 |
| All records | 7 | `GCF_001039075.2` | 22 | 504,366 |
| All records | 8 | `GCF_000172015.1` | 21 | 472,315 |
| All records | 9 | `GCF_014058445.2` | 20 | 816,525 |
| All records | 10 | `GCF_002156845.1` | 20 | 726,863 |
| Flag-positive (transposable = 1.0) | 1 | `GCF_002863665.1` | 11 | 442,805 |
| Flag-positive (transposable = 1.0) | 2 | `GCF_003018395.1` | 10 | 329,748 |
| Flag-positive (transposable = 1.0) | 3 | `GCF_002863685.1` | 10 | 325,023 |
| Flag-positive (transposable = 1.0) | 4 | `GCF_000026345.1` | 9 | 297,512 |
| Flag-positive (transposable = 1.0) | 5 | `GCF_014623405.1` | 9 | 272,060 |
| Flag-positive (transposable = 1.0) | 6 | `GCF_004358365.1` | 8 | 389,650 |
| Flag-positive (transposable = 1.0) | 7 | `GCF_003018075.1` | 8 | 349,832 |
| Flag-positive (transposable = 1.0) | 8 | `GCF_001420935.1` | 8 | 336,859 |
| Flag-positive (transposable = 1.0) | 9 | `GCF_003019015.1` | 8 | 323,355 |
| Flag-positive (transposable = 1.0) | 10 | `GCF_003017805.1` | 8 | 312,457 |
| Taxonomy-assigned (taxonomy != Unknown) | 1 | `GCF_000703225.1` | 23 | 560,341 |
| Taxonomy-assigned (taxonomy != Unknown) | 2 | `GCF_001039215.2` | 22 | 555,483 |
| Taxonomy-assigned (taxonomy != Unknown) | 3 | `GCF_001039075.2` | 22 | 504,366 |
| Taxonomy-assigned (taxonomy != Unknown) | 4 | `GCF_000194335.1` | 22 | 376,904 |
| Taxonomy-assigned (taxonomy != Unknown) | 5 | `GCF_910593925.1` | 20 | 684,213 |
| Taxonomy-assigned (taxonomy != Unknown) | 6 | `GCF_000220005.1` | 20 | 607,462 |
| Taxonomy-assigned (taxonomy != Unknown) | 7 | `GCF_001039155.2` | 20 | 565,133 |
| Taxonomy-assigned (taxonomy != Unknown) | 8 | `GCF_000172015.1` | 20 | 466,436 |
| Taxonomy-assigned (taxonomy != Unknown) | 9 | `GCF_003018515.1` | 19 | 883,763 |
| Taxonomy-assigned (taxonomy != Unknown) | 10 | `GCF_014058445.2` | 19 | 815,535 |

Top 10 by inclusive bp for every scope (each rank is among N=26,077 denominator genomes):

| Scope | Rank | Genome exact key | Unique loci | Total inclusive bp |
|---|---:|---|---:|---:|
| All records | 1 | `GCF_003112225.1` | 19 | 903,677 |
| All records | 2 | `GCF_003018515.1` | 19 | 883,763 |
| All records | 3 | `GCF_003018455.1` | 19 | 878,326 |
| All records | 4 | `GCF_003018055.1` | 18 | 827,949 |
| All records | 5 | `GCF_003018035.1` | 19 | 824,146 |
| All records | 6 | `GCF_014058445.2` | 20 | 816,525 |
| All records | 7 | `GCF_003966465.1` | 18 | 787,138 |
| All records | 8 | `GCF_001695515.1` | 19 | 772,916 |
| All records | 9 | `GCF_017165395.1` | 19 | 756,818 |
| All records | 10 | `GCF_002156845.1` | 20 | 726,863 |
| Flag-positive (transposable = 1.0) | 1 | `GCF_002863665.1` | 11 | 442,805 |
| Flag-positive (transposable = 1.0) | 2 | `GCF_004358365.1` | 8 | 389,650 |
| Flag-positive (transposable = 1.0) | 3 | `GCF_005037715.1` | 7 | 387,715 |
| Flag-positive (transposable = 1.0) | 4 | `GCF_003018515.1` | 6 | 379,954 |
| Flag-positive (transposable = 1.0) | 5 | `GCF_003018035.1` | 5 | 359,543 |
| Flag-positive (transposable = 1.0) | 6 | `GCF_003112225.1` | 5 | 351,238 |
| Flag-positive (transposable = 1.0) | 7 | `GCF_003018075.1` | 8 | 349,832 |
| Flag-positive (transposable = 1.0) | 8 | `GCF_003018835.2` | 6 | 339,502 |
| Flag-positive (transposable = 1.0) | 9 | `GCF_001420935.1` | 8 | 336,859 |
| Flag-positive (transposable = 1.0) | 10 | `GCF_003018775.1` | 5 | 336,185 |
| Taxonomy-assigned (taxonomy != Unknown) | 1 | `GCF_003112225.1` | 17 | 892,706 |
| Taxonomy-assigned (taxonomy != Unknown) | 2 | `GCF_003018515.1` | 19 | 883,763 |
| Taxonomy-assigned (taxonomy != Unknown) | 3 | `GCF_003018455.1` | 18 | 875,432 |
| Taxonomy-assigned (taxonomy != Unknown) | 4 | `GCF_003018055.1` | 16 | 818,397 |
| Taxonomy-assigned (taxonomy != Unknown) | 5 | `GCF_003018035.1` | 18 | 817,129 |
| Taxonomy-assigned (taxonomy != Unknown) | 6 | `GCF_014058445.2` | 19 | 815,535 |
| Taxonomy-assigned (taxonomy != Unknown) | 7 | `GCF_003966465.1` | 18 | 787,138 |
| Taxonomy-assigned (taxonomy != Unknown) | 8 | `GCF_001695515.1` | 18 | 766,078 |
| Taxonomy-assigned (taxonomy != Unknown) | 9 | `GCF_017165395.1` | 16 | 734,979 |
| Taxonomy-assigned (taxonomy != Unknown) | 10 | `GCF_015687505.1` | 15 | 723,314 |

Shortest and longest unique loci (primary inclusive sensitivity):

| Tail | Rank | Genome | Scaffold | Begin–end | bp | Prophage ID |
|---|---:|---|---|---:|---:|---|
| shortest | 1 | `GCF_000619605.1` | `NZ_JHHM01000131.1` | 3,327–3,587 | 261 | `GCF_000619605.1_prophage_9` |
| shortest | 2 | `GCF_000194215.1` | `NZ_AEZJ02000027.1` | 30,691–31,017 | 327 | `GCF_000194215.1_prophage_6` |
| shortest | 3 | `GCF_002134695.1` | `NZ_NEOP01000023.1` | 82,233–82,579 | 347 | `GCF_002134695.1_prophage_3` |
| shortest | 4 | `GCF_002510595.1` | `NZ_NLWG01000021.1` | 48,104–48,450 | 347 | `GCF_002510595.1_prophage_1` |
| shortest | 5 | `GCF_020554025.1` | `NZ_JAJCGM010000053.1` | 2,091–2,440 | 350 | `GCF_020554025.1_prophage_5` |
| shortest | 6 | `GCF_020554065.1` | `NZ_JAJCGN010000052.1` | 2,091–2,440 | 350 | `GCF_020554065.1_prophage_5` |
| shortest | 7 | `GCF_000511425.1` | `NZ_AYJX01000017.1` | 61,043–61,419 | 377 | `GCF_000511425.1_prophage_1` |
| shortest | 8 | `GCF_001614875.1` | `NZ_LVMG01000240.1` | 3,891–4,267 | 377 | `GCF_001614875.1_prophage_1` |
| shortest | 9 | `GCF_005389745.1` | `NZ_BFVM01000156.1` | 3,484–3,860 | 377 | `GCF_005389745.1_prophage_3` |
| shortest | 10 | `GCF_006350215.1` | `NZ_VESX01000063.1` | 23,661–24,037 | 377 | `GCF_006350215.1_prophage_5` |
| longest | 1 | `GCF_014216695.1` | `NZ_CP050202.1` | 1,745,719–1,968,002 | 222,284 | `GCF_014216695.1_prophage_2` |
| longest | 2 | `GCF_014216715.1` | `NZ_CP050203.1` | 5,145,073–5,333,645 | 188,573 | `GCF_014216715.1_prophage_6` |
| longest | 3 | `GCF_900449745.1` | `NZ_UGEE01000003.1` | 2,373–188,770 | 186,398 | `GCF_900449745.1_prophage_1` |
| longest | 4 | `GCF_014216755.1` | `NZ_CP050205.1` | 5,045,580–5,218,270 | 172,691 | `GCF_014216755.1_prophage_4` |
| longest | 5 | `GCF_001677475.2` | `NZ_CP015229.1` | 1,794,569–1,960,962 | 166,394 | `GCF_001677475.2_prophage_4` |
| longest | 6 | `GCF_013182975.1` | `NZ_SNOJ01000006.1` | 1,783–166,566 | 164,784 | `GCF_013182975.1_prophage_1` |
| longest | 7 | `GCF_013183005.1` | `NZ_SNOL01000012.1` | 1,783–166,566 | 164,784 | `GCF_013183005.1_prophage_2` |
| longest | 8 | `GCF_013183115.1` | `NZ_SNOQ01000009.1` | 3,185–166,542 | 163,358 | `GCF_013183115.1_prophage_2` |
| longest | 9 | `GCF_002196645.1` | `NZ_MVOU01000011.1` | 3,479–164,706 | 161,228 | `GCF_002196645.1_prophage_1` |
| longest | 10 | `GCF_000599665.1` | `NZ_CP007392.1` | 2,086,491–2,246,127 | 159,637 | `GCF_000599665.1_prophage_2` |

## Contig/scaffold and category distributions

There are **105,793 exact original scaffold keys / 132,404 records** and **0 / 105,793 keys assigned to multiple exact genome keys**. Prefix counts are identifier-form distributions only:

| Prefix | Records | Unique scaffold keys | Record denominator |
|---|---:|---:|---:|
| `NZ_` | 131,957 (99.66%) | 105,729 / 105,793 | 132,404 all records |
| `NC_` | 447 (0.34%) | 64 / 105,793 | 132,404 all records |

No chromosome/contig role field exists. In particular, `NC_`/`NZ_` prefixes alone are not used to label a sequence chromosome, plasmid, complete, or draft. The most reused exact scaffold keys (reuse here means multiple prophage records on the same sequence) are:

| Scaffold | Records | Record denominator |
|---|---:|---:|
| `NZ_CP050218.1` | 23 | 132,404 |
| `NZ_CP015831.1` | 20 | 132,404 |
| `NZ_CP024618.1` | 20 | 132,404 |
| `NZ_AGTD01000001.1` | 19 | 132,404 |
| `NZ_CP076230.1` | 19 | 132,404 |
| `NZ_CP016625.1` | 19 | 132,404 |
| `NZ_CP021339.1` | 19 | 132,404 |
| `NZ_CP015853.1` | 19 | 132,404 |
| `NZ_CP040305.1` | 19 | 132,404 |
| `NZ_CP038300.1` | 19 | 132,404 |

Exact taxonomy strings (mixed labels are kept intact):

| Taxonomy exact string | Records | Denominator |
|---|---:|---:|
| Myoviridae | 49,342 (37.27%) | 132,404 all records |
| Siphoviridae | 36,897 (27.87%) | 132,404 all records |
| Podoviridae | 20,420 (15.42%) | 132,404 all records |
| Unknown | 16,962 (12.81%) | 132,404 all records |
| Siphoviridae / Podoviridae | 4,314 (3.26%) | 132,404 all records |
| Myoviridae / Podoviridae | 2,395 (1.81%) | 132,404 all records |
| Myoviridae / Siphoviridae | 1,025 (0.77%) | 132,404 all records |
| Myoviridae / Siphoviridae / Podoviridae | 844 (0.64%) | 132,404 all records |
| Siphoviridae / Podoviridae / Bicaudaviridae | 150 (0.11%) | 132,404 all records |
| Myoviridae / Siphoviridae / Bicaudaviridae | 15 (0.01%) | 132,404 all records |
| Myoviridae / Siphoviridae / Podoviridae / Fuselloviridae | 13 (0.01%) | 132,404 all records |
| Siphoviridae / Bicaudaviridae | 7 (0.01%) | 132,404 all records |
| Microviridae | 5 (0.00%) | 132,404 all records |
| Siphoviridae / Microviridae | 4 (0.00%) | 132,404 all records |
| Bicaudaviridae | 3 (0.00%) | 132,404 all records |
| Myoviridae / Bicaudaviridae | 2 (0.00%) | 132,404 all records |
| Inoviridae | 2 (0.00%) | 132,404 all records |
| Myoviridae / Siphoviridae / Microviridae | 1 (0.00%) | 132,404 all records |
| Myoviridae / Siphoviridae / Podoviridae / Bicaudaviridae | 1 (0.00%) | 132,404 all records |
| Myoviridae / Microviridae | 1 (0.00%) | 132,404 all records |
| Podoviridae / Bicaudaviridae | 1 (0.00%) | 132,404 all records |

`category_counts.tsv` provides the same dimensions for all three scopes, plus field missingness, coordinate validity, diagnostic length class, per-genome zero/one/multiple, exact hit-scaffold-count frequencies, and records-per-scaffold-key bins. Every count has an explicit denominator and fraction. All seven fields have 0 missing values in all-record rows; subset missingness is also tabulated.

## Lossless current-key crosswalk requirements

This task does not rename identifiers. For the future PanSN/BGZF crosswalk, coordinate survival requires retaining these exact source components:

- assembly key: `genome`, including `GCF_` prefix and version suffix;
- sequence key: `scaffold`, including `NC_`/`NZ_` prefix and version suffix;
- source `begin` and `end` strings plus their normalized integer values and an explicit coordinate-convention field;
- `prophage_id` and `source_row` for row-level traceability;
- the source input SHA-256 and the composite locus key `(genome, scaffold, begin, end)`.

Do not join on version-stripped assembly or scaffold keys merely because this snapshot has zero collision groups; future releases may contain multiple versions. Do not rely on `prophage_id` alone as a coordinate key. `interval_qc.tsv` is the lossless tabular handoff for every current normalized row: it stores both exact original keys, version-separated diagnostics, original numeric strings, normalized coordinates, row ID, join flags, and QC flags. The dedicated `pansn-bgzip-genome-layout` task owns the canonical old-to-new crosswalk and future names.

## Likely artifacts and limits

- The terminal accession token `genome` and the exact coextensiveness of 26,077 normalized accessions with 26,077 CSV genome keys suggest the denominator list may have been derived from the table (or at least selected to match it). Therefore 0 / 26,077 all-record zero-hit genomes is not an independent prevalence estimate.
- 12,279 / 132,404 unique loci begin at <=3 and 817 / 132,404 are <1,000 bp. These may reflect contig-edge/truncated or very short calls, but sequence/contig lengths and caller output are absent, so the cause is unverified.
- 1 / 132,404 loci exceed 200,000 bp; inspect the interval row before downstream use. Long calls are not automatically artifacts.
- Exact duplicate, overlap, and nesting checks found 0 extra duplicate-locus rows, 0 overlapping pairs, and 0 nested pairs / 132,404 all unique loci. This does not detect homology or duplicate biology across different assemblies/scaffolds.
- No sequence lengths, caller name/version/settings, taxonomy method/confidence, status definition, sample metadata, assembly completeness, or quality field is supplied. Out-of-bounds checks, right-edge clipping, biological validation, and completeness/quality stratification remain unsupported.

## Machine-readable outputs and plots

- `per_genome.tsv`: one row per normalized denominator genome, in input order; record/locus/valid-locus/hit-scaffold counts and both bp sensitivities for all three scopes.
- `summary_metrics.tsv`: reconciliation, joins, distributions, interval QC, semantic-field presence, and explicit denominators.
- `category_counts.tsv`: categorical distributions, missingness, fractions, and denominators.
- `interval_qc.tsv`: one row per normalized source record with exact keys, coordinates, dedup ranks, joins, and QC flags.
- `plots/per_genome_locus_counts.svg`: zero-inclusive genome-count histogram (N printed).
- `plots/per_genome_bp_ecdf.svg`: zero-inclusive bp ECDF (N printed; log10(1+bp) display).
- `plots/locus_length_distribution.svg`: scope-specific length-bin percentages with locus denominators.
- `plots/taxonomy_distribution.svg`: top exact taxonomy strings plus Other, with record denominator.

SVGs are deterministic, dependency-free, contain accessible `<title>`/`<desc>` text, and expose exact counts in bar/line tooltips where applicable.
