# Canonical cohort 250

**Verdict: PASS.** The exact frozen N=250 rung was prepared as immutable release
`canonical-cohort-250-v1-a6184d7d6ee08bda` at:

```text
/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-250/canonical-cohort-250-v1-a6184d7d6ee08bda/
```

The release has 250 terminal `VALIDATED` genome rows, 41,050 lossless contig
crosswalk rows, and 1,276,442,466 bases. It stores 150 new source/canonical
object pairs and resolves the validated first 100 read-only through immutable
object references. It did not copy, download, or recompress any predecessor
object.

## Gate applicability and immutable inputs

The graph's explicit dependency clarification makes canonical acquisition and
canonicalization independent of unresolved prophage coordinates. Therefore the
integrated N=100 `GO_250` verdict is recorded, not fabricated, as
`NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY`; integrated extraction/query remains
separately blocked. No human wait, waiver, conditional verdict, or manual
override was used. Every gate applicable to this canonical preparation is
unqualified `PASS`.

Execution verified tracked and external inventories, exact release IDs,
immutable/PASS verdicts, and exact bytes for:

- selection release `pilot-cohorts-v1-8afc0ea03d9e50dc`, tracked `release.json`
  SHA-256 `d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99`;
- exact 91,475-byte, 250-row `cohort-0250.tsv`, SHA-256
  `ba2cf2909ccf62a0c1944a76b522edc5600953511ec355479117b4a419acbc9f`;
- predecessor `canonical-cohort-100-v1-6be4c0dde65f31d0`, tracked and external
  `release.json` SHA-256
  `3b91b24e23323ef971a13f22825e512a233bb592ed641ea9b270a2f1fd683795`;
- consumer release `consumer-compatibility-v1-78d7e93f19fa3d87`, tracked and
  external `release.json` SHA-256
  `021719ddadd7bb7fa2932d2ef9cb25da9c666ebe0389988691283011ee12f4c7`,
  including all 19 consumer rows `PASS`.

The first 100 N=250 rows equal the ordered N=100 inventory byte-for-byte by
exact accession revision; all 250 rows are unique and retain frozen order. The
URL constructor was narrowed in-process to those 250 revisions, and only rows
101-250 reached acquisition. The graph-wide sequence-bearing union was exactly
100 at start and exactly the frozen 250 at finish, a subset of the single
frozen collection and below the hard 1,000-revision cap.

Both root inputs matched at start and finish and were not edited:

- `26k_ecoli_accession.txt`:
  `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`
- `26k_prophage1.csv`:
  `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`

## Reuse, acquisition, and canonical QC

Rows 1-10 retain their physical N=10 storage references and rows 11-100 retain
N=100 storage references; the N=250 manifest resolves both forms through the
N=100 predecessor. Every reused source and canonical `SHA256SUMS`, manifest,
and eight role-specific artifact digests matched the N=100 assembly table.
There were zero predecessor acquisition or canonicalization events.

The 150 new source packages total 303,448,257 bytes; all 250 referenced source
packages total 507,185,175 bytes. Every new package passed exact accession
identity in path/catalog/report, ZIP structure and CRC, upstream MD5 coverage,
and local package/member SHA-256.

FASTA members were streamed through the unchanged pinned workflow and pinned
`bgzip (htslib) 1.19`; no routine plain FASTA was retained. The 150 new
canonical BGZF files total 218,160,437 bytes, and all 250 referenced BGZF files
total 364,650,284 bytes. Every assembly passed:

- rename-only sequence digest, length, and record-order equality;
- exact `assembly_accession.version#1#reversible_contig` PanSN grammar,
  reversibility, and cohort-wide uniqueness;
- BGZF integrity and compressed/content SHA-256;
- `.fai` name/length/order, structural `.gzi`, and every-contig
  `samtools faidx` prefix-region round trip;
- source/canonical manifest, crosswalk, and role-specific checksum accounting;
- atomic object inventory and last-marker publication.

All 237 publishable GFF alias views independently reproduced their exact
seqid aliases and 1-based-closed coordinate summaries. Thirteen optional alias
views (four inherited, nine new) reproducibly failed source-coordinate bounds
and remain explicitly quarantined; genome objects remain fully validated, and
no clipped or partial annotation view was published. Exact rows and failure
reasons are in `artifacts/canonical_cohort_250/exceptions.tsv`.

## Batches, transfer, retry, and forced restart

The compact `batch_metrics.tsv` has 25 deterministic batches of ten. Each row
records per-batch and cumulative validated transfer/canonical bytes, wall and
driver CPU, peak RSS, stage/partial/final bytes, files, live free bytes/inodes,
requests, retries, failures, and restart events. The release records 399 live
preflights: before every batch, acquisition, canonicalization, batch end, and
promotion.

Across all invocations there were 151 bounded GET requests for 150 completed
new packages, one unsafe-identity partial restart, one retry event, and zero
transport failure-ledger rows. The restart sequence was:

1. acquisition self-`SIGKILL` after fsyncing 131,072 partial bytes;
2. final release and `COMPLETE` remained absent;
3. resume rejected unsafe remote range identity, deleted the partial, and
   acquired the exact package from byte zero;
4. conversion self-`SIGKILL` after 262,217 streamed bases;
5. final release and `COMPLETE` again remained absent;
6. resume discarded the 20,124-byte interrupted conversion stage, restreamed,
   validated, and atomically committed it.

The ten predecessor-only batches transparently recorded that reserved partial
until its owning batch 11 resumed it; batch 11 and all later batches ended with
zero partial bytes. The final state ends in `READY_TO_PROMOTE`. `COMPLETE` was
created and fsynced last, followed by same-filesystem stage-directory rename.
No staging directory survived. The final release is 581,575,254 bytes and
2,269 files. External `SHA256SUMS` SHA-256 is
`45fd42b76bf1c1ace3a2e882fe6a9a8f6af2457c0f5d4bc28011cb99b521c5b7`,
which is the digest recorded by `COMPLETE`.

The successful resume used 414.48 s wall, 311.59 s user, and 53.52 s system
CPU by `/usr/bin/time`; the driver recorded 411.74 s wall and 180.22 s process
CPU. Including both deliberately killed invocations, cumulative execution was
421.35 s wall and 370.40 s user+system CPU, with 243,048,448-byte peak RSS and
zero process swap events.

## Resource and scale evidence

Every live preflight records `findmnt`, filesystem identity, owner/mode,
successful write probe, free bytes/inodes, swap, explicit allocations, and the
unfinished-write reservation.

| Resource | Allocation | Predicted upper-95% / projected | Measured peak | Gate |
|---|---:|---:|---:|---|
| RAM | 17,179,869,184 B | — | 238,854,144 B (1.39%) | ≤70% PASS |
| Durable disk | 4,000,000,000 B | 727,775,238 B (18.19%) | 581,297,080 B (14.53%) | ≤70% PASS |
| Scratch disk | 4,000,000,000,000 B | 1,200,000,000 B | no biological scratch payload | ≤70% PASS |
| Inodes | 200,000 | 5,000 files (2.5%) | 2,264 files (1.132%) | ≤50% PASS |
| Unfinished write | 1,000,000,000 B | 2× retained in allocations/live space | 20,124 B interrupted stage | PASS |

Durable free space remained 2,630,974,799,872 bytes at successful-run start and
2,630,272,970,752 bytes at promotion, with more than 447 million free inodes.
Scratch remained above 5.506 TB and 1.492 billion
free inodes. Thus durable 2 TB/1 million and scratch 4 TB preflight plus 2 TB/5
million stop floors all passed. Sampled task-lifetime swap free never fell
below its start, process swap events were zero, and no OOM occurred.

The modeled upper-95 disk peak is the N=100 measured peak scaled by new-object
count plus a conservative 25% allowance. Observed N=250 peak was 99.84% of the
linear estimate and below the modeled bound. Descriptive last-two preparation
slopes show no unexplained growth jump:

| Slope | N=100 | N=250 | Change |
|---|---:|---:|---:|
| successful wall / new object | 2.964 s | 2.763 s | -6.78% |
| successful wall / new source byte | 1.455 µs | 1.366 µs | -6.15% |
| stage peak / new object | 3,881,468 B | 3,875,314 B | -0.16% |
| stage files / new object | 15.156 | 15.093 | -0.41% |

The new-object wall-time point exponent is 0.863. These are preparation
telemetry, not a substitute for the integrated biological time-exponent upper
bound or per-base slope gate; those remain
`NOT_APPLICABLE_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS` here.

## Validation, determinism, and compact handoff

The independent validator re-hashes the complete external inventory and every
referenced predecessor object, revalidates all 250 archives, all 41,050 BGZF
index/name/sequence region checks, all crosswalks and GFF alias decisions,
resource/restart ledgers, exact manifest rows, root digests, and finish union.
The project regression command passes 63/63 tests; the task plus pinned
acquisition suite passes 21/21, including manifest/row and tracked-inventory
checksum mismatch, blank/overallocated and live-floor resource refusal,
bounded URL/retry behavior, invalid completed object discard, and interrupted
promotion.

A deterministic existing-release rerun made zero additional network requests
and zero canonical recomputations. External and tracked trees, state ledger,
file/byte counts, and semantic validation JSON remained byte-identical.

Git handoff is limited to:

- `manifests/canonical-cohort-250-v1/`: exact input, compact assembly/state/
  checksum/reference/batch tables, deterministic compressed contigs, external
  inventory, release JSON, and tracked `SHA256SUMS`;
- `artifacts/canonical_cohort_250/`: compact runner/validator/tests, metrics,
  exceptions, validation/restart/resource/time evidence;
- this report.

No package, source FASTA/GFF, canonical BGZF, `.fai`, `.gzi`, raw cache, or
per-hit biological payload is in Git, and no task-owned Git artifact exceeds
10 MiB.
