# Canonical cohort 1,000

**Verdict: PASS.** The exact frozen N=1,000 ceiling is immutable release
`canonical-cohort-1000-v1-4bc3e029e6e0be44` at:

```text
/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-1000/canonical-cohort-1000-v1-4bc3e029e6e0be44/
```

It contains 1,000 terminal `VALIDATED` genome rows,
171,659 lossless contig rows, and 5,104,867,881 bases.
The first 500 objects are read-only digest references; only
rows 501-1,000 were acquired and canonicalized. No 1,001st revision was
requested, analyzed, or published, and this release does not authorize any
beyond-1,000 projection.

## Automatic immutable gates

Execution consumed selection `pilot-cohorts-v1-8afc0ea03d9e50dc`, tracked `release.json`
SHA-256 `d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99`, and the exact 365,970-byte,
1,000-row cohort SHA-256 `265a1e7784a4d5db3ea3577892feba8173290518b6c621f7e5091dbad66bfe77`. It consumed predecessor
`canonical-cohort-500-v1-c208c111708b6435`, release SHA-256
`1df34fa84f2d4223242cad20529bea331c67289653f082a8f6da4d5036764dd5`, plus checksum-pinned canonical scale
trend SHA-256 `11edf88a6a9b5b126917cb4d734ed683aeb01a56647d0631067e4162c84ae686`. The predecessor's
`applicable_gates.canonical_scale_trend`, `scale_trend.json` verdict, every
scale check, and every N=1,000 projection check were unqualified `PASS`.
All 19 pinned-consumer gates were also `PASS`. Missing, mismatched,
`CONDITIONAL`, or failed inputs are rejected by code and tests; no prompt,
wait, waiver, shrink, reconstruction, or substitution path exists.

Canonical preparation records prophage source/coordinate policy and integrated
extraction as `NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED`.
It does not fabricate a coordinate or integrated-analysis PASS.

Both immutable roots matched at start and finish and were not edited:

- `26k_ecoli_accession.txt`: `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`
- `26k_prophage1.csv`: `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`

The graph-wide sequence-bearing union was exactly 500 at start and the exact
frozen 1,000 at finish, wholly within the selected collection and hard cap.

## Reuse, acquisition, canonicalization, and quarantine

All 500 predecessor references and all eight role-specific digests per row
matched the immutable N=500 chain. The request ledger contains no predecessor
acquisition or canonicalization event. The 500 incremental packages total
1,015,169,053 bytes and
incremental canonical BGZF files total
729,360,986 bytes.

Every referenced assembly passed exact accession/version package identity,
ZIP/CRC/upstream-MD5/local-SHA checks; rename-only sequence identity; exact,
unique and reversible PanSN names; BGZF integrity; `.fai`/`.gzi`; every-contig
`samtools faidx` prefix-region round trip; crosswalk; source-GFF alias/bounds;
object inventory; and atomic marker validation. There are zero genome-object
quarantines. 52 optional source-GFF
alias views reproducibly failed their own bounds and remain explicit in
`artifacts/canonical_cohort_1000/exceptions.tsv`; no partial alias view was
published. The other 948 annotation
views reproduced exactly. This annotation QC makes no prophage-coordinate claim.

## Batches, retry, and injected restart

The release has 100 deterministic ten-row batches and
5,076 append-only live resource records. The
injected acquisition process was SIGKILLed after fsyncing 131,072 partial bytes
(9.89 s); neither final release nor `COMPLETE`
appeared. Resume rejected unsafe range identity and reacquired that exact row.
Conversion was then SIGKILLed after crossing the configured 262,144-base threshold
(11.00 s), again with no final publication. Resume
discarded the interrupted conversion stage and independently validated the
completed unit. Batch evidence records 4
restart events and 262,144 partial
bytes; final partial bytes are zero.

The payload-completion attempt used 1420.59 s but
correctly stopped before promotion when sampled system swap-free decreased;
every other resource and scale check passed. That failure is append-only.
Two checksum-only promotion retries also refused on continued system swap
growth after 235.25
s total. Once live swap-free stabilized, the next checksum-only retry revalidated
all units and observed zero growth, using 130.52 s wall and
129.98 s user+system. Across six
measured invocations, wall was 1807.25
s, CPU was 1682.69 s, outer peak
RSS was 583,913,472 bytes, and process
swap events were zero. The append ledger additionally exposes
1
checksum-only code-correction interruption(s), with no payload work or
publication. There were 501
bounded GET requests for 500 completed
incremental packages and 3 explicit
failure-ledger event(s), compactly mirrored in
`artifacts/canonical_cohort_1000/failures.tsv`.

`COMPLETE` was written/fsynced last and the full staging directory was promoted
by same-filesystem rename. No stage survived. External `SHA256SUMS` SHA-256 is
`4acc432f408300fa22df39532489b3d9f210ab0e71883bbfafb513cddfe5ce9f` and `COMPLETE` SHA-256 is `c94134fcd8fe04f9c4545f3ccb1bdb8d0679d0e16cfe9bf26f18876026530443`.

## Live resources and final scaling model

| Resource | Allocation | Measured / upper model | Gate |
|---|---:|---:|---|
| RAM | 34,359,738,368 B | 583,913,472 B outer peak (1.699%) | <=70% PASS |
| Durable disk | 15,000,000,000 B | 1,991,202,841 B measured; 2,462,105,290 B modeled upper-95 | <=70% PASS |
| Scratch disk | 4,000,000,000,000 B | 3,000,000,000 B reserved | <=70% PASS |
| Inodes | 400,000 | 7,515 files; 15,000 configured | <=50% PASS |
| Unfinished write | 2,000,000,000 B | at least 2x retained live/allocation | PASS |

Every stage/batch record includes `findmnt`, mount/source/fstype, ownership and
write probe, live free bytes/inodes, assigned quotas, swap, and unfinished-write
reservation. Durable free space remained at least
2,586,436,419,584
bytes with at least
447,658,270
free inodes; scratch remained at least
5,505,572,864,000
bytes and
1,492,409,893
inodes. All durable 2 TB/1M-inode and scratch 4 TB preflight plus 2 TB/5M stop
floors passed. The successful promotion retry had no OOM, process swap, or
system swap growth; every earlier system-wide swap-growth refusal remains
explicit in the immutable failure/resource ledgers.

The final adjacent N=500 to N=1,000 time exponent is
**1.035320**; the empirical upper
bound is **1.035320** (limit 1.3).

| Per-new-base slope | N=500 | N=1,000 | Change |
|---|---:|---:|---:|
| Wall seconds | 5.42775e-07 | 5.551285e-07 | +2.276% |
| Source bytes | 0.3971922 | 0.3974875 | +0.074% |
| Stage bytes | 0.7727522 | 0.7766756 | +0.508% |
| Stage files | 2.953407e-06 | 2.942092e-06 | -0.383% |
| Peak RSS bytes | 0.2772929 | 0.2084257 | -24.836% |

Every absolute change is <=25%. All measured N=1,000 values stayed within the
checksum-pinned predecessor upper-95 projection. `scale_trend.json` publishes
descriptive final fits over N=10/100/250/500/1,000 and explicitly records
`projection_beyond_ceiling=NOT_AUTHORIZED_NOT_COMPUTED`.

## Independent validation, determinism, and compact handoff

Two full independent validators produced byte-identical PASS JSON SHA-256
`fc8dc803369f82ee39fe974e7d2a12d9c699dce1aba58077e8e528edfada7397`. Each rehashed inventories;
validated all 1,000 archives/object references and
171,659 contig/index/name/region checks; recomputed every
annotation decision and final scale check; audited restart/resources/roots/
exact cap; and rejected partial, unlisted, symlink, or plain-FASTA release data.
Task tests cover manifest/inventory mismatch, non-PASS gates, live-resource
refusal, a 1,001-row union, scale regression, bounded URL/retry, corrupt
completed units, and interrupted promotion. The task suite passes 14/14 and the
full project regression suite passes 95/95.

A post-validation existing-release rerun made zero network requests and zero
canonical recomputations. External tree SHA-256
`222625d03436ecc74e602f90f63fe26bbddcf5646dd7de83f7c45ba76b92c6d8` and tracked-manifest tree
SHA-256 `f61eff6fbbed499dc0d7c73c0e96a0b8d77d1aa16ccbfb201b9feece13282e9d` were unchanged,
as were state bytes and request/commit counts.

Git contains only compact files under `manifests/canonical-cohort-1000-v1/`,
`artifacts/canonical_cohort_1000/`, and this report. No package, genome/prophage
sequence, source FASTA/GFF, canonical BGZF, `.fai`, `.gzi`, raw cache, whole
index, or per-hit biological output is tracked; each task-owned file is under
10 MiB. Bulky objects remain solely in this task's external namespace, and no
dependent analysis is authorized beyond N=1,000.
