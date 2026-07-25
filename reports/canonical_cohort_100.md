# Canonical cohort 100

**Verdict: PASS.** The frozen N=100 rung was prepared as immutable release
`canonical-cohort-100-v1-6be4c0dde65f31d0` at:

```text
/home/erikg/phind-data/ecoli26k/v1/releases/prepare-canonical-cohort-100/canonical-cohort-100-v1-6be4c0dde65f31d0/
```

The release contains 100 terminal `VALIDATED` genome rows, 18,098 lossless
contig crosswalk rows, and 512,421,261 bases. It contains 90 newly acquired
source/canonical object pairs and resolves the validated first 10 read-only by
predecessor release ID plus source/canonical object-inventory digest. It did
not download, copy, or recompress those first 10 objects.

## Immutable automatic input gates

Execution proceeded without a human wait only after verifying all of the
following tracked and external inventories, exact release IDs, immutable/PASS
verdicts, and manifest bytes:

- frozen selection release `pilot-cohorts-v1-8afc0ea03d9e50dc`, tracked
  `release.json` SHA-256
  `d134f5a31deff39ac1614df0ecf20ce91a1388f1e9673c0f41efd231d2b5eb99`;
- exact 36,573-byte, 100-row `cohort-0100.tsv`, SHA-256
  `13e203961a9fcec18a8a09e690582652d8085b2a386811e6c6a03184b9489182`;
- predecessor `canonical-cohort-010-v1-e71484de9994fc28`, tracked and external
  `release.json` SHA-256
  `4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23`;
- compatibility release `consumer-compatibility-v1-78d7e93f19fa3d87`,
  tracked and external `release.json` SHA-256
  `021719ddadd7bb7fa2932d2ef9cb25da9c666ebe0389988691283011ee12f4c7`,
  with all 19 required consumer rows `PASS`;
- the N=100 first ten accession revisions exactly equal the ordered N=10
  release inventory; all 100 exact revisions are unique and in frozen order.

Both immutable root inputs matched before and after work:

- `26k_ecoli_accession.txt`:
  `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`
- `26k_prophage1.csv`:
  `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`

The committed graph-wide sequence-bearing union was 10 at start and exactly
the frozen 100 at finish. The URL constructor was narrowed in-process to only
rows 1-100 and acquisition was invoked only for rows 11-100. The global union
is a subset of the frozen collection and below the hard cap of 1,000.

## Source and canonical objects

The task reused the unchanged production primitives in
`workflow/acquisition_canonicalization/pilot.py`; it did not modify that pinned
workflow. Each of the 90 new source objects passed exact accession directory
and catalog identity, ZIP structure/CRC, upstream MD5 coverage, and local
package/member SHA-256. The 100 source packages total 203,736,918 bytes;
183,308,530 bytes belong to the newly acquired 90.

FASTA members were streamed directly from their packages through pinned
`bgzip (htslib) 1.19`, with no retained plain FASTA. Every output uses exact
`assembly_accession.version#1#reversible_contig` PanSN naming and passed:

- source/output sequence digest, length, and record-order equality;
- cohort-wide PanSN uniqueness and reversible source-header crosswalks;
- BGZF integrity and compressed/content SHA-256;
- `.fai` exact name/length/order, structural `.gzi`, and every-contig
  `samtools 1.19.2 faidx` prefix-region round trip;
- atomic object `SHA256SUMS` plus last `COMPLETE` publication.

The 100 canonical BGZF files total 146,489,847 bytes. All source GFF files
remain unchanged in their checksum-validated source packages. Ninety-six
assemblies have validated alias-only annotation views after every GFF seqid
and 1-based-closed coordinate passed. Four optional alias views were explicitly
withheld, never clipped or partially published:

| Order | Exact revision | Terminal annotation action | Evidence |
|---:|---|---|---|
| 8 | `GCF_000167895.3` | `QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW` | row 622, `90076-90333` out of range |
| 23 | `GCF_005885955.1` | `QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW` | row 11359, `93426-102935` out of range |
| 68 | `GCF_001617565.1` | `QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW` | row 10490, `82553-82984` out of range |
| 75 | `GCF_001577325.1` | `QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW` | row 10312, `113369-114226` out of range |

The genome objects for these four assemblies remain rename-only validated;
only invalid coordinate-bearing optional annotation views are quarantined.
Exact reasons are in `artifacts/canonical_cohort_100/exceptions.tsv` and the
per-object manifests.

## Batches, retries, and forced restart

Work ran in ten deterministic batches of ten, with live preflight before each
batch and every acquisition/canonicalization stage. The successful resume took
265.06 seconds wall and 116.33 seconds process CPU (the enclosing
`/usr/bin/time` measured 266.78 seconds wall, 198.27 seconds user, and 37.85
seconds system including child tools). Batch 2-10 wall times were 25.97-31.95
seconds; the first batch only validated the ten predecessor references.

There were 92 bounded payload GET attempts for 90 completed new packages.
One initial chunked-transfer `IncompleteRead(847 bytes read)` stopped without
publication and was retained in the append-only failure ledger. The bounded
transport retry wrapper then re-entered the pinned primitive's unsafe-partial
identity check rather than appending untrusted bytes.

The required real restart sequence then demonstrated:

1. acquisition self-`SIGKILL` after fsyncing 131,072 partial bytes;
2. no final release or `COMPLETE`, followed by an identity-unsafe delete and
   exact byte-zero retry because NCBI supplied no strong resumable identity;
3. source-object atomic commit followed by conversion self-`SIGKILL` after
   262,150 bases;
4. no final publication, followed by
   `INTERRUPTED_CONVERSION_STAGE_DISCARDED`, full restream, and validated
   canonical-object commit.

The final release state ends with `READY_TO_PROMOTE`; `COMPLETE` was created
and fsynced last, and the whole stage was promoted by same-filesystem rename.
The release has 1,369 files and 349,501,982 bytes. Its `SHA256SUMS` SHA-256 is
`06cd0c25a23bb9c48e73f5175def14e97f5a008d41ff65d8c4d8ccecc15e0176`,
which is the digest recorded by `COMPLETE`. No staging directory survived.

## Resource gates

Every one of 205 live preflight records is `PASS` and records `findmnt`, owner,
mode, successful write probe, free bytes/inodes, swap, explicit allocations,
and unfinished-write reservation. Allocations and observed peaks were:

| Resource | Allocation / bound | Predicted upper-95% | Measured peak | Gate |
|---|---:|---:|---:|---|
| RAM | 8,589,934,592 B | — | 142,684,160 B (1.66%) | ≤70% PASS |
| Durable disk | 2,000,000,000 B | 500,000,000 B (25.0%) | 349,332,114 B (17.47%) | ≤70% PASS |
| Scratch disk | 4,000,000,000,000 B | 500,000,000 B | no biological scratch payload | ≤70% PASS |
| Inodes | 100,000 | 2,000 files (2.0%) | 1,364 files (1.364%) | ≤50% PASS |
| Unfinished write | 500,000,000 B | 2× retained in allocations and live space | PASS | PASS |

Durable free space remained above 2 TB and 447 million free inodes. Scratch
preflight remained above 4 TB and 1.492 billion free inodes, with the 2 TB/5
million stop floors preserved. Process swap events and observed system swap
growth were zero. This task prepares canonical inputs and performs no
scale-bearing integrated biological analysis, so the time-exponent and
last-two-rung per-base-slope gate is explicitly
`NOT_APPLICABLE_PREPARATION_ONLY_NO_INTEGRATED_ANALYSIS`; downstream integrated
N=100 analysis must evaluate it.

## Validation and compact handoff

The independent validator re-hashes the complete release and all referenced
N=10 objects, revalidates every source package, BGZF, index, crosswalk, alias
policy, role-specific checksum, resource record, restart event, root digest,
and finish union. Task and predecessor regression tests cover manifest and
checksum mismatch, blank/over-allocated resources, payload-cap refusal,
transport retry, archive corruption, safe/unsafe partial handling, conversion
restart, and interrupted promotion.

Tracked handoff files are limited to:

- `manifests/canonical-cohort-100-v1/`: exact input, compact assembly/state/
  checksum/reference/batch tables, deterministic compressed contigs,
  external inventory, release JSON, and tracked `SHA256SUMS`;
- `artifacts/canonical_cohort_100/`: compact runner/validator/tests, metrics,
  exceptions, validation/restart/resource/time evidence;
- this report.

No package, source FASTA/GFF, BGZF, `.fai`, `.gzi`, or other biological payload
is in Git, and no tracked task artifact exceeds 10 MiB.
