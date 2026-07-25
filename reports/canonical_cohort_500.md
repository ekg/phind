# Canonical cohort 500

**Verdict: PASS.** The exact frozen N=500 rung was prepared as immutable release
`canonical-cohort-500-v1-c208c111708b6435` at:

```text
/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-500/canonical-cohort-500-v1-c208c111708b6435/
```

The release has 500 terminal `VALIDATED` genome rows, 83,096 lossless contig
crosswalk rows, and 2,550,902,887 bases. It stores 250 new source/canonical
object pairs and resolves the validated first 250 read-only through immutable
digest references. No predecessor object was copied, downloaded, or
recompressed.

## Automatic gate and immutable inputs

Execution automatically consumed the canonical-only audit verdict
`CANONICAL_GO_500`, SHA-256
`699d14b32c010771280b193b2373968dcae0c0c130a87a91270413eadd9c03e5`.
That verdict pins the canonical N=250 audit result SHA-256
`8c5ce43261dd99d6c4ba6b52bf0e5e4e6a9c86168f5c95b452f2ab192cd0d8d1`
and the exact predecessor release tree. Every applicable audit and release gate
was unqualified `PASS` before acquisition began.

This is canonical preparation only. The release explicitly records both
prophage source/coordinate interpretation and integrated extraction as
`NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED`; it does not
invent an integrated GO or a coordinate-policy PASS. No human wait, waiver,
conditional verdict, shrink, substitution, or manual override was used.

Execution verified tracked and external inventories, exact release IDs,
immutable/PASS verdicts, row counts, order, and exact bytes for:

- selection release `pilot-cohorts-v1-8afc0ea03d9e50dc`, tracked `release.json`
  SHA-256 `d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99`;
- exact 182,605-byte, 500-row `cohort-0500.tsv`, SHA-256
  `bb6497bff230ecb6987dc5cda865307524a7fd63653e12e0d54d0808afe15ecb`;
- predecessor `canonical-cohort-250-v1-a6184d7d6ee08bda`, tracked/external
  `release.json` SHA-256
  `dcf2b887afa51e4e0e739ae2fef9b5a9d72fb8bc9a4d698a161a99673aaf504a`;
- consumer release `consumer-compatibility-v1-78d7e93f19fa3d87`,
  `release.json` SHA-256
  `021719ddadd7bb7fa2932d2ef9cb25da9c666ebe0389988691283011ee12f4c7`,
  including all 19 consumer rows `PASS`.

The first 250 N=500 rows equal the ordered N=250 accession inventory exactly.
All 500 exact revisions are unique and retain frozen order. The URL constructor
was narrowed in-process to only those revisions, and only rows 251-500 reached
acquisition. The global sequence-bearing union was exactly 250 at start and the
exact frozen 500 at finish, remains a subset of the frozen N=1,000 collection,
and is below the hard 1,000-revision cap.

Both root inputs matched before and after work and were not edited:

- `26k_ecoli_accession.txt`:
  `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`
- `26k_prophage1.csv`:
  `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`

## Reuse, acquisition, and canonical QC

Rows 1-10 retain N=10 physical references, rows 11-100 retain N=100
references, rows 101-250 retain N=250 references, and rows 251-500 are
self-contained in this task's release. Every reused object inventory and all
eight role-specific digests matched the audited predecessor chain. There were
zero predecessor acquisition or canonicalization events.

The 250 new source packages total 506,205,744 bytes; all 500 referenced source
packages total 1,013,390,919 bytes. Every new package passed exact accession
identity in path/catalog/report, ZIP structure and CRC, upstream MD5 coverage,
and local package/member SHA-256.

FASTA members were streamed through the unchanged pinned workflow and pinned
`bgzip (htslib) 1.19`; no routine plain FASTA was retained. The 250 new
canonical BGZF files total 363,855,997 bytes, and all 500 referenced BGZF files
total 728,506,281 bytes. Every assembly passed:

- rename-only sequence digest, length, and record-order equality;
- exact `assembly_accession.version#1#reversible_contig` PanSN grammar,
  reversibility, and cohort-wide uniqueness;
- BGZF integrity and compressed/content SHA-256;
- `.fai` name/length/order, structural `.gzi`, and every-contig
  `samtools faidx` prefix-region round trip;
- source/canonical manifest, crosswalk, and role-specific checksum accounting;
- atomic object inventory and final-marker publication.

All 473 publishable source-GFF alias views independently reproduced their exact
seqid aliases and 1-based-closed bounds. Twenty-seven optional alias views
reproducibly failed source-GFF bounds and remain quarantined; no clipped or
partial alias view was published. This annotation validation does not assert
prophage coordinate semantics. Exact rows and failure reasons are in
`artifacts/canonical_cohort_500/exceptions.tsv`.

## Batches, bounded retry, and forced restart

The compact `batch_metrics.tsv` has 50 deterministic batches of ten. It records
per-batch/cumulative validated transfer and canonical bytes, wall/CPU, RSS,
stage/partial/final bytes, files, live free bytes/inodes, requests, retries,
failures, and restart events. The append-only external ledger contains 709 live
preflights: before every batch and acquisition/canonicalization stage, at batch
ends, and before promotion.

Across all invocations there were 251 bounded GET requests for 250 completed
new packages, one unsafe-identity partial restart, one retry, and zero failure
ledger rows. The injected restart sequence was:

1. acquisition self-`SIGKILL` after fsyncing 131,072 partial bytes;
2. final release and `COMPLETE` remained absent;
3. resume rejected unsafe remote range identity, removed the partial, and
   reacquired the exact package from byte zero;
4. conversion self-`SIGKILL` after 262,195 streamed bases;
5. final release and `COMPLETE` again remained absent;
6. resume discarded the interrupted conversion stage, restreamed it, and
   independently validated the completed unit.

The 25 predecessor-only batches transparently recorded the reserved interrupted
conversion stage; its owning batch 26 discarded and resumed it, and batch 26
through 50 ended with zero partial bytes. The final state event is
`READY_TO_PROMOTE`; `COMPLETE` was written and fsynced last, followed by a
same-filesystem stage-directory rename. No staging directory survived. The
final release is 985,310,126 bytes and 3,770 files. External `SHA256SUMS`
SHA-256 is
`214cf5b1885004f1faa5ddda5e9423c73f78c2cc85e6c36aaca23c0151bc3b18`,
and its digest is pinned by `COMPLETE`.

The successful resume used 698.55 s wall, 543.19 s user, and 99.82 s system by
`/usr/bin/time`; the driver recorded 691.75 s wall and 321.46 s process CPU.
Including both intentionally killed invocations, execution used 714.43 s wall
and 657.29 s user+system CPU. Outer peak RSS was 396,283,904 bytes and all
three invocations recorded zero process swap events.

## Live resources and canonical scale gates

Every preflight records `findmnt`, filesystem identity, owner/mode, successful
write probe, free bytes/inodes, swap, positive assigned quotas, and unfinished
write reservation.

| Resource | Allocation | N=500 modeled/measured | Gate |
|---|---:|---:|---|
| RAM | 34,359,738,368 B | 396,283,904 B outer peak (1.153%) | ≤70% PASS |
| Durable disk | 10,000,000,000 B | 1,211,035,583 B modeled upper-95; 984,842,116 B measured | ≤70% PASS |
| Scratch disk | 4,000,000,000,000 B | 3,000,000,000 B reserved prediction; no biological scratch payload | ≤70% PASS |
| Inodes | 400,000 | 12,000 configured; 3,764 measured | ≤50% PASS |
| Unfinished write | 2,000,000,000 B | at least 2× retained live and in allocations | PASS |

Durable free space remained 2,628,507,979,776 bytes at successful-run start and
2,627,618,496,512 bytes at promotion, with over 447 million free inodes.
Scratch remained at 5.505 TB and 1.492 billion free inodes. Durable 2 TB/1
million and scratch 4 TB preflight plus 2 TB/5 million stop floors all passed.
Swap free increased from 13,516,800 to 14,188,544 bytes over the measured run;
process swap events were zero and no OOM occurred.

The canonical scale gate compares the adjacent N=100→250 and N=250→500
incremental workloads. Its empirical time-exponent upper bound is **1.015672**
(limit 1.3). Every last-two-rung per-new-base slope change is within 25%:

| Per-new-base slope | N=250 | N=500 | Change |
|---|---:|---:|---:|
| wall seconds | 5.3891e-7 | 5.4277e-7 | +0.717% |
| source bytes | 0.397173 | 0.397192 | +0.005% |
| stage bytes | 0.760839 | 0.772752 | +1.566% |
| stage files | 2.9633e-6 | 2.9534e-6 | -0.333% |
| peak RSS bytes | 0.312628 | 0.277293 | -11.303% |

The N=1,000 upper-95 projection fits all available N=10, 100, 250, and 500
rung observations, then takes the maximum of the fit and linear N=500 scaling
plus 25%. The legacy N=10 Stage-B point is used for conservative projection but
is not mislabeled as a canonical scale-bearing adjacent transition.

| N=1,000 projected measure | Upper-95 | Comparison allocation/gate |
|---|---:|---|
| Wall | 2,865 s | 7,200 s: PASS |
| Source transfer | 1,269,009,070 B | bounded frozen 500-object increment |
| Stage disk | 2,489,053,920 B | 3.0 GB configured projection; 10 GB allocation: PASS |
| Peak RSS | 883,496,960 B | 34.36 GB assigned, ≤70%: PASS |
| Files | 9,410 | 12,000 configured; 400,000 inodes, ≤50%: PASS |

All projection reservation, disk, RAM, file, inode, and time checks are `PASS`.
This is canonical preparation telemetry, not an integrated biological result;
the integrated branch remains blocked.

## Validation, determinism, and compact handoff

The independent validator re-hashes the exact external and tracked inventories;
revalidates all 500 archives, all 83,096 BGZF/index/name/sequence region checks,
all recursive references, crosswalks and source-GFF alias decisions; recomputes
the scale projection; audits resource/restart ledgers; verifies root digests and
the exact global N=500 finish union; and rejects unlisted, partial, symlink, or
plain-FASTA release files. Task failure tests cover manifest/row and audit
checksum mismatch, blank/overallocated/live-floor resource refusal, canonical
scale refusal, bounded URL/retry behavior, invalid completed object cleanup,
and interrupted promotion. The full project-plus-audit suite passes 81/81.

An existing-release rerun made zero new network requests and zero canonical
recomputations. External tree SHA-256
`711e2720d2890dd391a74f8c757b389bf9277f0202bff3ed8682a2d328a96936`
and tracked-manifest tree SHA-256
`c1b0af06f26f1b7808716fd4c09143d85b8610928305095f8c165f13ae4e7780`
were unchanged, as were the append-only state ledger and request/commit counts.
Two full independent semantic validations produce byte-identical compact JSON.

The Git handoff is limited to:

- `manifests/canonical-cohort-500-v1/`: exact input, compact assembly/state/
  checksum/reference/batch tables, three deterministic compressed contig
  manifest parts, external inventory, release JSON, and tracked `SHA256SUMS`;
- `artifacts/canonical_cohort_500/`: compact runner/validator/tests, metrics,
  exceptions, validation/restart/resource/scale/time evidence;
- this report.

No package, genome/prophage sequence, source FASTA/GFF, canonical BGZF, `.fai`,
`.gzi`, raw cache, per-hit biological output, or whole index is in Git. Every
new task-owned Git artifact is below 10 MiB.
