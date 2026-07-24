# Stage-B 10-assembly acquisition and canonicalization

**Verdict: PASS for the bounded pilot.** Immutable release `canonical-cohort-010-v1-e71484de9994fc28` contains exactly the ten frozen Stage-B exact assembly revisions, ten checksum-complete source objects, and ten checksum-complete PanSN/BGZF canonical objects. No scale-bearing rung or eleventh assembly was requested.

External release:

```text
/home/erikg/phind-data/ecoli26k/v1/releases/run-10-assembly-acquisition/canonical-cohort-010-v1-e71484de9994fc28/
```

Its 169 files occupy 37,615,933 bytes. `SHA256SUMS` hashes to `96a40035c15684d4c3c12c88f8134c32c4df421eb9d138119581ab7473badc44`; `COMPLETE` contains that digest and was fsynced last before the same-filesystem directory rename. The independent validator rechecked the entire inventory and all source/canonical semantics.

## Immutable input and bound gates

The workflow verified tracked and external predecessor `collection-v1-f7494b4b89d1382b`, tracked `release.json` SHA-256 `59c6907e2c053e9d8ac3df8d5eb820bab0097030a9259ca2c9354c47cb6642bf`, and the exact 2,246-byte, 10-row Stage-B manifest SHA-256 `0d179cbafce2ba1fa14d1929a4acd6621810a335f25bcd7ec67dd2083eb101f6`. Its order and row hashes were not reconstructed or changed. The documented predecessor status `EXACT_VERSION_VALID_METADATA_UNAVAILABLE` for `GCF_000167895.3` remains explicit.

Both root inputs matched at start and finish:

- `26k_ecoli_accession.txt`: `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`
- `26k_prophage1.csv`: `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`

The committed graph-wide sequence-bearing union was zero before this pilot and exactly these ten revisions afterward. It is a subset of all 26,077 frozen candidates and below the hard cap of 1,000. The URL constructor rejects any versionless, substituted, or non-Stage-B accession before network I/O. Retries transferred the same exact first revision twice after the injected interruption, but the distinct sequence-bearing identity union never exceeded ten.

## Source acquisition

Acquisition was sequential, rate-limited, and bounded to single-accession NCBI Datasets v2 packages requesting genomic FASTA, GFF3, and sequence report. The server returned neither ETag, Last-Modified, Content-Length, nor byte-range support. Consequently, the interrupted 131,072-byte first partial was **not** appended: the restart recorded `ACQUISITION_PARTIAL_RESTART_IDENTITY_UNSAFE`, deleted it, and fetched that exact revision from byte zero. This is the fail-closed remote-identity behavior; range resume is separately covered by a strong-ETag unit test.

Every source object passed ZIP structure/CRC, traversal/symlink, exact package accession directory/catalog identity, all upstream MD5 entries, and local package/member SHA-256 gates. NCBI omits only generic top-level `README.md` from `md5sum.txt`; ZIP CRC and local archive SHA-256 cover it, while all accession data, reports, FASTA, GFF, and catalog members are upstream-checksummed. Each source directory has its manifest, receipt/remote identity, `SHA256SUMS`, and last `COMPLETE`, and was atomically renamed before canonicalization consumed it. The compact assembly table's `download_attempts` is the final successful receipt invocation; the append-only state ledger is authoritative for the first assembly's two total GET attempts across `SIGKILL`.

## Canonical objects and annotation policy

Conversion streamed the ZIP FASTA member directly into pinned `bgzip -@ 2 -l 6 --binary -c`. No source or cohort-wide plain FASTA was extracted or retained. Each source record became exact `resolved_accession.version#1#reversibly_encoded_source_token`; safe bytes passed unchanged and all unsafe UTF-8 bytes use uppercase `%HH` escaping (with the policy's digest alias only above the length bound). Each canonical file passed:

- per-contig source/output sequence SHA-256, length, and record-order equality;
- cohort-wide unique PanSN names and complete 1,223-row crosswalk accounting;
- `bgzip -t`, canonical-content SHA-256, and compressed-artifact SHA-256;
- `.fai` exact names/lengths/order, structural `.gzi` validation, and `samtools faidx` PanSN region checks for every contig;
- atomic per-assembly manifest/inventory/`COMPLETE` publication.

Source GFF remains unchanged inside the source package. A transformed GFF was intentionally not published because no pinned transformed-GFF consumer gate exists. Nine assemblies instead received validated explicit alias tables after every GFF seqid resolved to a FASTA/sequence-report alias and every 1-based-closed interval was in range.

The documented metadata-omission assembly provided an exact payload but exposed a real source inconsistency: `GCF_000167895.3` GFF row 622 uses `NZ_CP021212.1:90076-90333` although that source sequence is 90,229 bases. The workflow stopped that attempt, then applied the already specified “annotation view only where validated” policy: its entire derived annotation alias view is quarantined/empty with the exact reason retained; its unchanged source GFF remains checksummed; its genome FASTA/canonical BGZF remains rename-only validated. Thus no invalid coordinate-bearing annotation was promoted or silently clipped. The source-coordinate/annotation-policy gate is PASS because the invalid optional view was withheld, not bypassed.

| Order | Exact assembly | Package bytes | Bases | Contigs | Canonical BGZF bytes | Annotation-view state |
|---:|---|---:|---:|---:|---:|---|
| 1 | `GCF_000005845.2` | 1,790,946 | 4,641,652 | 1 | 1,314,927 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 2 | `GCF_000812325.1` | 2,066,665 | 5,205,271 | 237 | 1,501,230 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 3 | `GCF_002302315.1` | 2,085,137 | 5,319,464 | 3 | 1,504,775 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 4 | `GCF_004664255.1` | 2,034,978 | 5,096,076 | 144 | 1,454,522 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 5 | `GCF_015644385.1` | 2,046,069 | 5,121,256 | 139 | 1,465,307 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 6 | `GCF_020829045.1` | 2,027,488 | 5,182,045 | 8 | 1,463,145 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 7 | `GCF_921380995.1` | 2,063,121 | 5,232,387 | 159 | 1,492,512 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 8 | `GCF_000167895.3` | 2,250,915 | 5,654,428 | 265 | 1,618,499 | `QUARANTINED_SOURCE_GFF_NO_ALIAS_VIEW` |
| 9 | `GCF_001881595.4` | 1,970,974 | 4,957,459 | 111 | 1,413,336 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |
| 10 | `GCF_000498835.2` | 2,092,095 | 5,321,624 | 156 | 1,517,172 | `ALIASES_VALIDATED_NO_TRANSFORMED_GFF` |

Totals are 20,428,388 packaged bytes, 52,500,819 source FASTA bytes, 51,731,662 bases, 1,223 contigs, and 14,745,425 canonical BGZF bytes. BGZF is 0.285 bytes/base; package transfer is 1.385 times canonical BGZF. These are measured pilot ratios, not a production extrapolation confidence interval.

## Forced interruption, resume, and failure evidence

The acquisition invocation self-sent `SIGKILL` after fsyncing 131,072 partial bytes. No source object, release directory, or `COMPLETE` was visible. The second injected invocation checksum-committed the source object, then self-sent `SIGKILL` after 262,160 FASTA bases; it left only a truncated BGZF in a hidden canonical stage. Restart logged `INTERRUPTED_CONVERSION_STAGE_DISCARDED`, never appended the BGZF, streamed it again, and committed a validated object. At both interruption points the overall final release and its `COMPLETE` were absent.

The append-only failure ledger also retains three fail-closed development/restart discoveries before the successful immutable publication: NCBI's generic README MD5 omission, overly broad sequence-report assembly-accession aliasing, and the real GFF bound failure above. Each invocation stopped without publication. The first two implementation issues were covered by tests before resume; the third now produces an explicit annotation-view quarantine rather than partial aliases. No independently checksum-complete source/canonical object was redownloaded or recompressed.

A post-publication rerun made zero network requests and did not acquire or recompress anything. External tree digest remained `fad738a5ae9bef6370fdc36e9366a8d552e4712896ecb9946146374dfe16b041`; tracked release-tree digest remained `4475115a500fb342678661195d74b1e9848c049297fd6ac28ae86d8e0bf86afa`. This establishes deterministic checksum-complete resume and byte-stable compact republication.

## Resource evidence

Before initial work, every acquisition/conversion unit, and promotion, the workflow recorded live mount identity, ownership/mode, write probe, exact free bytes/inodes, allocations, and reservations. There are 52 resource records across forced/restarted invocations. Mounts remained root ext4 `/dev/nvme0n1p2` for durable and XFS `/dev/nvme3n1` at `/mnt/nvme3n1` for scratch.

Explicit allocations were 8,589,934,592 bytes RAM, 10,000,000,000 bytes durable, 4,000,000,000,000 bytes scratch, and 100,000 inodes; predicted peaks were 500,000,000 bytes on each tier, 500 files, with a 500,000,000-byte unfinished-write reservation. At the successful publication invocation's start, durable/scratch free bytes were 2,635,671,732,224 / 5,506,473,230,336, with 447,716,010 / 1,492,409,952 free inodes. Finish values were 2,635,660,890,112 / 5,506,473,230,336 bytes and 447,715,965 / 1,492,409,952 inodes. All hard floors and 2× unfinished-write checks passed.

Peak RSS was 37,793,792 bytes (0.44% of assigned RAM); process swaps, `/usr/bin/time` swaps, and observed system swap growth were zero. Measured staged release peak was 37,588,400 bytes and 164 files, far below 70% disk and 50% inode allocations. The non-scale-bearing Stage-B pilot has no time exponent or last-two-rung per-base slope, so that gate is explicitly `NOT_APPLICABLE_STAGE_B_NON_SCALE_BEARING`, not PASS.

## Validation and tracked handoff

`python -m unittest -v workflow.acquisition_canonicalization.test_pilot` passes 11 tests covering manifest/checksum mismatch, resource refusal, the exact-accession cap, archive/upstream corruption, safe/unsafe resume identity, streamed naming/index/GFF aliases, annotation quarantine, and interrupted conversion/promotion. `validate_release.py` independently reports PASS for 12 gate groups, all 10 assembly rows, 1,223 contig rows, 80 role-specific checksum rows, 51,731,662 bases, exact root hashes, and the finish global union of ten.

The tracked handoff is `manifests/canonical-cohort-010-v1/`: exact Stage-B input, assembly/state/checksum tables, deterministic compressed contig inventory, external inventory, release JSON, and tracked `SHA256SUMS`. Compact validation/resource/restart/time evidence is under `artifacts/acquisition_canonicalization_10/`. Git contains no package, FASTA, GFF, BGZF, `.fai`, `.gzi`, or other biological payload.
