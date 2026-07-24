# Prophage source semantics and extraction gate

## Release verdict

**`EXTRACTION_BLOCKED`** — consumer action **`REJECT`**.

The immutable CSV has no attributable producer/caller version or transformation
provenance. Phigaro v2.3.0 is a strong **format hypothesis**, but Phigaro v2.4.0
changed the same tabular boundaries and added an ID column. The available
annotation-boundary diagnostic favors the raw-coordinate/1-based hypothesis but
is not a producer declaration or known-base oracle. Selecting a coordinate
convention would therefore be a guess.

Current release:

- ID: `prophage-semantics-v1-f5619e221ff272ae`
- durable path: `/home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source/prophage-semantics-v1-f5619e221ff272ae/`
- tracked reference: `artifacts/prophage_semantics/release_reference.json`
- policy: `artifacts/prophage_semantics/semantics_policy_v1.json`
- schema: `workflow/prophage_semantics/semantics-policy-v1.schema.json`

An extraction-dependent consumer must require exactly `EXTRACTION_GO`; missing,
`CONDITIONAL`, unknown, or `EXTRACTION_BLOCKED` is rejected. The demonstrated
strict consumer check exits 2 for this release.

## Automatic immutable-input gate

All gates were evaluated automatically; no approval, override, substituted
object, new assembly download, or reduced cohort was used.

| Input/evidence | Required identity | Result |
|---|---|---|
| `26k_ecoli_accession.txt` | SHA-256 `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5` | PASS at start and finish |
| `26k_prophage1.csv` | SHA-256 `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996` | PASS at start and finish |
| integrated plan | SHA-256 `fb58d25a6f4971137ab0dcb82dae09eac5d177e37ba34f857845ce2d7e0a6da8` | PASS |
| source audit | SHA-256 `feb9b687fb0722a4f073f7105088f19a2dabebcbc0378dbe5c4799b0b7f29fdc` | PASS |
| canonical N=10 manifest | SHA-256 `4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23` | PASS |
| predecessor release | `canonical-cohort-010-v1-e71484de9994fc28`, exact order of 10 revisions | PASS |
| predecessor external release | `COMPLETE` plus every row of its `SHA256SUMS` | PASS |
| global exact-assembly union | exactly 10 revisions, declared cap 1,000 | PASS |

The predecessor's applicable identity, checksum, row, BGZF/index/name,
coordinate-policy, compatibility, restart, and publication gates were checked
before its objects were read. In the semantics release, every applicable
engineering gate is exactly `PASS`, including the source-policy gate (the
fail-closed dual policy is valid) and pinned-consumer compatibility gate (the
strict consumer correctly rejects non-GO). The separate scientific
`extraction_eligibility` verdict is `EXTRACTION_BLOCKED`; this is the hard stop,
not an engineering-gate bypass. The N=10 order is immutable and appears in
`sentinel_summary.json`. One predecessor GFF alias view was quarantined for an
out-of-range source feature; the immutable source GFF remains checksum-pinned and
was used only as an explicitly non-decisive boundary diagnostic.

## Evidence, in priority order

### 1. Source provenance

No producer documentation or sidecar is present. The CSV header is exactly
`end,genome,scaffold,begin,transposable,taxonomy,prophage_id`; it contains no
caller, version, command/config, output format, strand, topology, completeness,
edge, duplicate-policy, or `tagged` field. Repository reports preserve that
negative result (`reports/prophage_distribution.md`, “Schema and semantic
evidence”; integrated plan D4/D5 and §4.2). This is the highest-priority result.

A final bounded metadata-only retry inspected the file's complete Git history,
its creation-commit tree, all current tracked text references, and three exact
filename/header/field web queries. Commit
`7381847222cd7cd25aa3042ec7eeeacfb528092b` added only the accession and CSV
files with the message “Add ecoli accession and prophage data files”; it added no
pipeline or sidecar. The web queries produced no relevant exact attribution.
This non-decisive attempt is recorded in
`artifacts/prophage_semantics/targeted_provenance_retry.json`; it read no bases,
made no assembly download, and does not alter the immutable release or verdict.

### 2. Version-matched candidate caller code

The CSV's five semantic fields fingerprint Phigaro ≤2.3 TSV after plausible
aggregation:

- Phigaro v2.3.0 writes `scaffold,begin,end,transposable,taxonomy` (plus optional
  `vog`) and assigns TSV begin/end directly from boundary gene coordinates
  ([v2.3.0 output header and boundaries](https://github.com/bobeobibo/phigaro/blob/aea9469d09cdbfbb528998ebc43232ee9f44decd/phigaro/batch/task/run_phigaro.py#L88-L105),
  [lines 139–142](https://github.com/bobeobibo/phigaro/blob/aea9469d09cdbfbb528998ebc43232ee9f44decd/phigaro/batch/task/run_phigaro.py#L139-L142),
  [written rows](https://github.com/bobeobibo/phigaro/blob/aea9469d09cdbfbb528998ebc43232ee9f44decd/phigaro/batch/task/run_phigaro.py#L227-L254)).
- In that version, `transposable=True` is returned when the retained ordered pVOG
  records contain a second `Integration` group hit without an intervening group
  other than `Other`; it is not a completeness or quality flag
  ([candidate algorithm](https://github.com/bobeobibo/phigaro/blob/aea9469d09cdbfbb528998ebc43232ee9f44decd/phigaro/to_html/preprocess.py#L223-L237)).
- Candidate taxonomy is the joined set of tied maximum-frequency taxonomy codes
  among mapped pVOGs; no mapped code yields `Unknown`. Slash-separated labels are
  ties, not validated modern taxonomy
  ([candidate classifier](https://github.com/bobeobibo/phigaro/blob/aea9469d09cdbfbb528998ebc43232ee9f44decd/phigaro/data.py#L256-L271)).

This fingerprint does **not** attribute the CSV. The aggregate's `genome` and
`prophage_id` could have been added or rewritten by unversioned post-processing.
Critically, v2.4.0 added IDs and computes TSV boundaries as boundary-gene
`begin-1,end-1`
([v2.4 TSV header](https://github.com/bobeobibo/phigaro/blob/1ff5f85cee31e418bce24e4cd51c7528c43bc968/phigaro/batch/task/run_phigaro.py#L87-L102),
[boundary computation](https://github.com/bobeobibo/phigaro/blob/1ff5f85cee31e418bce24e4cd51c7528c43bc968/phigaro/batch/task/run_phigaro.py#L135-L139),
[written rows](https://github.com/bobeobibo/phigaro/blob/1ff5f85cee31e418bce24e4cd51c7528c43bc968/phigaro/batch/task/run_phigaro.py#L247-L256)).
Its release notes also document a coordinate fix affecting ≤2.3 GFF/BED output
([v2.4.0 tracker](https://github.com/bobeobibo/phigaro/blob/1ff5f85cee31e418bce24e4cd51c7528c43bc968/version_tracker.md#L1-L8)).
The exact candidate commits and full-file SHA-256 values are pinned in
`artifacts/prophage_semantics/evidence_inventory.json`.

### 3. Bounded sentinel diagnostic

Exactly 56 CSV rows belong to the frozen ten assemblies (per-assembly counts
5,3,8,5,3,4,7,10,6,5 in cohort order). Source GFF feature boundaries were read
inside the already pinned package ZIPs; **zero FASTA bases** and zero new
assemblies were read or downloaded.

| Add delta to both raw boundaries | Begin matches | End matches | Both match | Denominator |
|---:|---:|---:|---:|---:|
| -1 | 2 | 0 | 0 | 56 rows |
| 0 | 49 | 52 | 45 | 56 rows |
| +1 | 0 | 0 | 0 | 56 rows |

Both coordinate candidates are in contig range for all 56 rows. Thus raw
boundaries align strongly with NCBI annotation boundaries, while raw+1 does not.
The result is **`NON_DECISIVE`**: NCBI GFF is not the missing caller's Prodigal
output or producer provenance, 11 calls do not match both boundaries, all 1,223
predecessor contigs have topology `unknown`, and there is no independently
specified expected first/last base. Merely hashing the two candidate slices
would prove they differ, not which is correct. Therefore no “known-base” result
is claimed and no source sequence was copied.

## Versioned semantic decisions

Every dimension has a status, confidence, evidence list, and extraction-critical
flag in the machine policy.

| Dimension | Status / confidence | Supported policy |
|---|---|---|
| producer/caller/version | `UNRESOLVED` / `NONE` | Phigaro ≤2.3 remains an unattributed hypothesis; blocks extraction |
| user's term `tagged` | `UNRESOLVED` / `NONE` | no generic tagged subset; it is not mapped silently |
| `transposable` | `HYPOTHESIS_NOT_ATTRIBUTABLE` / `LOW` | preserve raw 0.0/1.0; candidate Phigaro meaning may be reported only as hypothesis; never quality/completeness |
| taxonomy labels | `HYPOTHESIS_NOT_ATTRIBUTABLE` / `LOW` | preserve exact strings; mixed labels and `Unknown` are not rewritten or modernized |
| coordinate base/end | `UNRESOLVED` / `NONE` | no canonical interval is emitted; blocks extraction |
| strand/orientation | `ABSENT_UNKNOWN` / `HIGH` | never infer strand or reverse-complement; blocks oriented extraction |
| topology/circularity | `UNRESOLVED` / `NONE` | no wrap/rotation inference from ordered begin/end; blocks extraction |
| contig-edge behavior | `UNRESOLVED` / `NONE` | begin≤3 is a diagnostic, not truncation/completeness; blocks extraction |
| completeness | `ABSENT_UNKNOWN` / `HIGH` | no completeness subset or label exists |
| duplicate-locus rules | `UNRESOLVED` / `NONE` | producer behavior unknown; preserve every row; blocks extraction even though this snapshot has no duplicate loci |

### The three scopes remain distinct

| Scope ID | Reversible rule | Rows |
|---|---|---:|
| `all_records` | every immutable source row | 132,404 |
| `transposable_flag_positive` | `Decimal(raw transposable) == 1`; observed member text `1.0` | 7,695 |
| `taxonomy_assigned` | trimmed/case-folded taxonomy is neither empty nor `Unknown` | 115,442 |

All 132,404 rows, 132,405 physical lines, 26,077 genome keys, and every raw field
were accounted. There are zero duplicate exact locus groups, extra duplicate
locus rows, exact-record duplicate groups, and duplicate `prophage_id` groups.
That observation is not promoted into a producer deduplication rule. No normalized
row table was materialized; the immutable CSV remains the lossless row store.

## Dual-convention release plan

Two candidates remain active and neither is selected:

1. **C1, raw 1-based closed**: `[b,e] -> [b-1,e)`.
2. **C2, raw 0-based inclusive**: `[b,e] -> [b,e+1)`.

Release may change to `EXTRACTION_GO` only when digest-pinned producer/caller and
post-processing provenance identifies the relevant version **and** a
version-matched native TSV/BED/GFF/saved-FASTA fixture supplies declared expected
boundary bases. For each frozen sentinel, the selected policy must pass contig
bounds, expected first/last base, expected length, whole-slice digest,
source-to-canonical round-trip, topology/wrap, edge, and strand-state checks.
Unknown strand remains unknown; no test may infer it from sequence plausibility.

## Release engineering and resource evidence

The external release contains its exact input manifest, policy, evidence,
profile/diagnostic, tools and argv provenance, append-only state/failure/resource
ledgers, resource and restart summaries, `release.json`, and complete
`SHA256SUMS`. An injected exit 75 occurred after static units but before
publication; `COMPLETE` was absent. Restart independently checked each completed
unit's SHA-256, appended resume events, and published no mixed files. `COMPLETE`
was created and fsynced last, followed by same-filesystem atomic rename.

Resource allocations were nonblank: 8 GiB RAM, 1 GB durable, 4 TB scratch, and
100,000 inodes; predicted peaks were 50 MB/50 MB and 100 files with 10 MB
unfinished-write reservation. Live durable/scratch floors, write probes,
ownership, `findmnt`, bytes/inodes, ≤70% RAM/disk, ≤50% inode, ≥2× unfinished
writes, no swap growth, and promotion preflight all passed. This N=10 metadata
boundary diagnostic is non-scale-bearing, so exponent/slope gates are not
applicable.

Eight task tests cover root/cardinality accounting, policy completeness,
manifest/checksum mismatch, blank/overcommitted resource refusal, interrupted
unit validation and mixed-resume refusal, absent-`COMPLETE` rejection, exact
predecessor inventory/cap/order, and bounded sentinel counts. A deterministic
rerun left every release byte unchanged. No extraction, clustering, pangenome
processing, or new *E. coli* assembly download occurred.
