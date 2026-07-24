# Storage and genome-acquisition inventory

**Report/evidence date:** 2026-07-24 UTC.  Filesystem topology was captured at
2026-07-24T19:49:14.323652111Z; byte/inode and write-probe evidence at
2026-07-24T19:49:36.559153454Z.  Remote metadata checks ran between
2026-07-24T19:51:54Z and 19:56:20Z.  Values are point-in-time observations,
not reservations.

## Scope and decision

This is an inventory and capacity estimate only.  It did **not** download any
sequence payload, build indexes for the production collection, convert a
production FASTA, run a tree/pangenome job, or inspect prophage-table
statistics.  The exact PanSN names, FASTA-header transformation, and BGZF
conversion policy remain owned by `pansn-bgzip-genome-layout`.

A later, separately authorized acquisition can initially use the root
filesystem for durable inputs:

* **durable input:** `/home/erikg/phind-data/ecoli26k/v1`
* **separate work/scratch:**
  `/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1`
* **recommended checksum replica:**
  `/mnt/nvme2n1/erikg/phind-genome-backup/ecoli26k/v1`

These paths are recommendations; they were not created.  The durable path is a
project-associated sibling on the current `/` filesystem, not inside the Git
working tree; its existing parent `/home/erikg` is erikg:erikg 0750.  This
avoids making Git scan or accidentally stage genome payload.  At the
observation time, `/` meets the conservative 2.4 TB preflight guardrail.  Its 2.647 TB
user-available space leaves 647 GB above the hard 2.0 TB stop threshold.  That
fits the high acquisition envelope of 392 GB, including a modeled replica and
25% margin.  An all-collection plain-FASTA transient would add as much as 153
GB and still only leave about 102 GB above the stop threshold; therefore it is
explicitly not recommended.  Stream each source file into BGZF instead.
Potentially much larger IMPG/tree/pangenome scratch is a separate estimate and
must not be inferred from the small acquisition footprint.

## Input inventory

### Integrity and reconciliation

`26k_ecoli_accession.txt` was read as bytes.  It is ASCII, contains LF line
terminators only, is LF-terminated, and has no NULs or CRLFs.

| Measure | Count |
|---|---:|
| bytes | 417,239 |
| physical lines / inventory data rows | 26,078 |
| nonblank lines | 26,078 |
| blank lines | 0 |
| unique normalized nonblank values | 26,078 |
| distinct duplicate groups | 0 |
| rows in duplicate groups | 0 |
| duplicate excess rows | 0 |
| valid assembly accessions | 26,077 |
| malformed/header-like lines | 1 |

The sole malformed value is line 26,078, `genome`.  No whitespace changed on
normalization.  The version distribution among valid accessions is `.1`:
25,149; `.2`: 911; `.3`: 14; `.4`: 3.

Rules used by `artifacts/genome_input_inventory.tsv`:

* `normalized_identifier` is the ASCII raw line with surrounding whitespace
  removed.  `raw_identifier` excludes only the physical LF terminator.
* Blank means the normalized identifier is empty.
* Malformed means blank or failure of the fully anchored regex
  `^(GCF|GCA)_[0-9]{9}\.[1-9][0-9]*$`.
* Duplicate means exact equality of nonblank normalized identifiers.  Every
  member of a repeated value would receive the same stable `dup_NNNNNN`
  group.  Because there are none, `duplicate_group` is `.` in every row.
* `is_blank_or_malformed=true` is therefore present only on the last row.

Checksums before the audit and after writing the row inventory were identical:

```text
SHA-256  1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5
BLAKE2b  660bc051b28c52bba3519dce8fdf0312941d3a736ef994ba12b54595f69a6cc0452752660abb89fb7b9451d5ce2732dea084156efffff5295d6a95309473e611
```

Thus the input remained byte-for-byte unchanged.  A final validation repeats
SHA-256 rather than relying only on this observation.

### Identifier kind and metadata resolution

All 26,077 syntactically valid values have `GCF_` prefixes and are **versioned
NCBI RefSeq assembly accessions**, not nucleotide-sequence accessions.  There
are no `GCA_` GenBank assembly identifiers and no nucleotide-record identifiers
in the file.  NCBI documents `GCF_` + nine digits + version as the RefSeq
assembly format and explains that versions change when underlying sequences
change [NCBI assembly versioning and status](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/policies-annotation/genome-processing/version-status/)
(accessed 2026-07-24T19:50Z).  NCBI also distinguishes latest, replaced, and
suppressed assemblies; a valid accession is not proof that its status will
remain current at later acquisition time.

A metadata-only POST to
`https://api.ncbi.nlm.nih.gov/datasets/v2/genome/check` checked all 26,077
accessions in 27 ordered batches of at most 1,000.  From
2026-07-24T19:56:04.479873Z through 19:56:20.729923Z it returned all 26,077 in
`valid_assemblies`, zero invalid/no-result accessions, 496,030 response bytes,
and required zero retries.  Each TSV row records its POST batch and access
window.  `api_check_valid_assembly` means that this endpoint recognized the
exact accession.version on that date; it does not promise later payload
availability or establish a per-record organism/status snapshot.

### Ten-accession size/metadata preview (no sequence bytes)

Ten deterministic representatives were chosen by input position/version:

| Selection | Input line | Accession |
|---|---:|---|
| first | 1 | GCF_000005845.2 |
| 10th percentile | 2,609 | GCF_000812325.1 |
| 25th percentile | 6,520 | GCF_002302315.1 |
| median | 13,039 | GCF_004664255.1 |
| 75th percentile | 19,558 | GCF_015644385.1 |
| 90th percentile | 23,469 | GCF_020829045.1 |
| last valid | 26,077 | GCF_921380995.1 |
| first version 3 | 80 | GCF_000167895.3 |
| first version 4 | 4,531 | GCF_001881595.4 |
| median version 2 | 1,506 | GCF_000498835.2 |

Exact metadata URL/command shape (the CSV value was the ten accessions above,
in that order):

```bash
sample_csv='GCF_000005845.2,GCF_000812325.1,GCF_002302315.1,GCF_004664255.1,GCF_015644385.1,GCF_020829045.1,GCF_921380995.1,GCF_000167895.3,GCF_001881595.4,GCF_000498835.2'
curl -fsSL --retry 3 --retry-delay 1 --connect-timeout 10 --max-time 60 \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$sample_csv/dataset_report"
curl -fsSL --retry 3 --retry-delay 1 --connect-timeout 10 --max-time 60 \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$sample_csv/download_summary"
```

At 2026-07-24T19:53:09Z, one combined NCBI Datasets v2 `download_summary`
reported `record_count=10`, `all_genomic_fasta` 10 files / 14.721868 MB, and
`genome_gff` 10 files / 4.694811 MB.  The corresponding `dataset_report`
response was 34,617 bytes and returned nine detailed records; all nine said
current RefSeq assemblies and *Escherichia coli*, with mean total sequence
length 5,119,692.667 bases and mean 106.56 contigs.  It omitted detailed record
GCF_000167895.3, although the separate `/check` endpoint recognized that exact
version as valid and `download_summary` counted all ten.  This discrepancy is
preserved as uncertainty, not silently filled.

The size-preview response was 1,582 bytes and uses the API's `size_mb` fields.
Those fields are treated as decimal-MB packaged-transfer estimates here; they
are not direct measurements of canonical BGZF output.  NCBI describes genome
packages, genomic FASTA, GFF3, the catalog, and decompressed-file `md5sum.txt`
in its [genome package reference](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/)
(accessed 2026-07-24T19:59Z).  The official
[large-download guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download/)
(accessed 2026-07-24T19:59Z) recommends the Datasets CLI and dehydrated package
workflow for at least 1,000 genomes or packages over 15 GB.

**Payload-sampling accounting:** 0 accessions, 0 sequence bytes, and 0
annotation payload bytes.  Only JSON metadata/size responses were transferred.
Temporary JSON response files used for the ten-record preview were removed.
No bounded payload directory exists to clean up.

## Filesystem evidence

### Devices, capacity, bytes, and inodes

`available` is `statvfs.f_bavail`, i.e. bytes available to this unprivileged
user; root's `free` includes filesystem-reserved blocks.  Decimal TB are shown
only for readability; exact bytes are authoritative.

| Path | source; filesystem | device classification | capacity bytes | free bytes | user-available bytes | free inodes |
|---|---|---|---:|---:|---:|---:|
| `/` | `/dev/nvme0n1p2`; ext4 | local NVMe partition (root/OS) | 15,239,634,198,528 | 3,415,259,312,128 | 2,647,149,449,216 | 447,717,990 |
| `/mnt/nvme1n1` | `/dev/nvme1n1`; XFS | local NVMe | 15,360,854,417,408 | 3,497,784,344,576 | 3,497,784,344,576 | 1,499,670,672 |
| `/mnt/nvme2n1` | `/dev/nvme2n1`; XFS | local NVMe | 15,360,854,417,408 | 4,171,527,573,504 | 4,171,527,573,504 | 1,499,001,374 |
| `/mnt/nvme3n1` | `/dev/nvme3n1`; XFS | local NVMe | 15,360,854,417,408 | 5,506,658,492,416 | 5,506,658,492,416 | 1,492,409,972 |

`lsblk -b` identified all four physical devices as non-rotating NVMe Samsung
MZQL215THBLA-00A07 devices.  Relevant `findmnt` sources are `/dev/nvme*`; none
of the four is NFS, MooseFS, or an unknown/network mount.  XFS data mounts use
`rw,nodev,noatime,nodiratime,...,inode64,noquota`; root is `rw,relatime` ext4.
No benchmark was run, so performance conclusions are qualitative: local NVMe
should offer much lower latency and higher throughput than NFS/MooseFS, while
root may contend with the OS and other work.  Each observed filesystem is on a
single device with no RAID/redundancy evidenced by `lsblk`; a same-host replica
is protection against accidental corruption/deletion, not site or chassis
failure.  The mounted data devices and root are shared, so free space can
change concurrently.

Evidence commands (read-only) were:

```bash
date --iso-8601=ns --utc
findmnt -T PATH -o TARGET,SOURCE,FSTYPE,OPTIONS,SIZE,AVAIL,USE% --bytes
df -B1 / /mnt/nvme1n1 /mnt/nvme2n1 /mnt/nvme3n1
df -i / /mnt/nvme1n1 /mnt/nvme2n1 /mnt/nvme3n1
stat -f -c '%T %S %b %f %a %c %d' PATH
lsblk -b -o NAME,KNAME,TYPE,TRAN,FSTYPE,SIZE,ROTA,MODEL,MOUNTPOINTS
```

### Ownership and writeability

| Mount/path | mountpoint owner/mode | tested Erik-owned directory owner/mode | write result |
|---|---|---|---|
| `/` | root:root 0755 | `/home/erikg/phind`, erikg:erikg 0775 | pass |
| `/mnt/nvme1n1` | root:root 0755 | `/mnt/nvme1n1/erikg`, erikg:erikg 0755 | pass |
| `/mnt/nvme2n1` | root:root 0755 | `/mnt/nvme2n1/erikg`, erikg:erikg 0755 | pass |
| `/mnt/nvme3n1` | root:root 0777 | `/mnt/nvme3n1/erikg`, erikg:erikg 0755 | pass |

Each probe used `set -o noclobber`, a unique name containing user, PID, random
number, and nanosecond timestamp, wrote 59 bytes, checked nonzero size, removed
it, then verified absence.  The exact basenames were:

```text
.wg-storage-probe.erikg.2219243.18273.1784922576578252923
.wg-storage-probe.erikg.2219243.21800.1784922576588401498
.wg-storage-probe.erikg.2219243.7056.1784922576597542102
.wg-storage-probe.erikg.2219243.1636.1784922576607270147
```

A trap removed every accumulated probe on normal exit, HUP, INT, or TERM.  A
post-probe `find` restricted to depth one in each exact candidate directory
reported zero `.wg-storage-probe.erikg.*` files.

## Bounded existing-data search

No whole-filesystem crawl was performed and no existing file was changed.
The search had three bounded stages:

1. `find` of each exact `/mnt/nvme*n1/erikg` root at `-xdev -mindepth 1
   -maxdepth 1`, capped to 250 displayed entries.  There were only 19, 9, and 9
   immediate entries on NVMe 1, 2, and 3 respectively.
2. Case-insensitive candidate-name matching (`genom`, `ecoli`, `coli`,
   `assembl`, `refseq`, `fasta`, `pangenom`, `phage`) to depth two under the
   project and those three Erik-owned roots, capped to 300 displayed matches.
3. A read-only Python filename scanner rooted only at
   `/home/erikg/phind` and the three `/mnt/.../erikg` paths.  It did not follow
   symlinks, required the same device, descended through directory depth three
   (thus inspected files inside those directories at relative path depth four),
   and imposed both 100,000 entries and 20 seconds per root.  It recognized
   only `.fa`, `.fna`, or `.fasta` with optional `.gz`, `.bgz`, or `.bgzf`,
   extracted versioned GCA/GCF accessions from basenames, and intersected them
   with this input list.

The shell commands for the listing/name stages were, for each explicitly
listed root (with output caps shown):

```bash
find "$root" -xdev -mindepth 1 -maxdepth 1 \
  -printf '%y\t%u\t%g\t%m\t%s\t%p\n' | sort | head -250
find "$root" -xdev -mindepth 1 -maxdepth 2 \
  \( -iname '*genom*' -o -iname '*ecoli*' -o -iname '*coli*' \
     -o -iname '*assembl*' -o -iname '*refseq*' -o -iname '*fasta*' \
     -o -iname '*pangenom*' -o -iname '*phage*' \) \
  -printf '%y\t%u\t%s\t%p\n' | sort | head -300
```

The here-document Python command used constants
`roots=[/home/erikg/phind,/mnt/nvme1n1/erikg,/mnt/nvme2n1/erikg,/mnt/nvme3n1/erikg]`,
`max_entries=100000`, `max_seconds=20`, `max_directory_depth=3`,
`follow_symlinks=False`, and regex
`(?i)\\.(fa|fna|fasta)(\\.(gz|bgz|bgzf))?$`; it counted every inspected entry
before applying the suffix/accession intersection.  This records the actual
bounds without implying a persistent helper script exists.

All four scans finished without hitting a cap:

| Root | entries examined | sequence-like candidates | exact input-accession filename hits | truncated |
|---|---:|---:|---:|---|
| `/home/erikg/phind` | 292 | 0 | 0 | no |
| `/mnt/nvme1n1/erikg` | 2,701 | 0 | 0 | no |
| `/mnt/nvme2n1/erikg` | 42,305 | 0 | 0 | no |
| `/mnt/nvme3n1/erikg` | 31,580 | 54 | 0 | no |

Relevant project files at that time were the 417,239-byte accession input and
12,393,209-byte `26k_prophage1.csv`; worktree copies of these small inputs are
not genome payloads.  NVMe 1's immediate entries were training/checkpoint
projects.  NVMe 2 held large language-model/text data and a `phrs` directory.
NVMe 3's 54 sequence-like hits were visibly human/primate, HPRC, VGP, and
Drosophila data by path/name, for example
`HPRCv2/chm13v2.0.psn.fa.gz`, `11way/GCA_028858775.2.fa.gz`, and
`drosophila/assemblies/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz`.
None contained an assembly accession from this 26,077-accession set in its
basename.  **Reusable existing input genomes found within the bounded scope:
zero.**  This does not prove none exist deeper than the depth bound, under
other users, elsewhere in `/home/erikg`, behind symlinks, in generic filenames,
or on unexamined mounts.

## Capacity and transfer model

All GB/TB in the estimate are decimal (10^9/10^12 bytes).  The complete
machine-readable model is `artifacts/genome_acquisition_estimate.tsv`.
Low/central/high are scenarios, not confidence intervals or promises.

| Component | low | central | high | Included in 392 GB high acquisition envelope? |
|---|---:|---:|---:|---|
| packaged genomic-FASTA transfer | 31.292 GB | 38.390 GB | 57.369 GB | transfer, not additive disk residence |
| compressed GFF3 transfer/retention | 7.823 GB | 12.243 GB | 26.077 GB | yes, retained copy |
| canonical BGZF FASTA | 31.292 GB | 41.723 GB | 62.585 GB | yes |
| `.fai` | 0.050 GB | 0.200 GB | 2.000 GB | yes |
| `.gzi` | 0.020 GB | 0.040 GB | 0.080 GB | yes |
| metadata/manifests/checksums | 0.050 GB | 0.250 GB | 1.000 GB | yes |
| **primary durable subtotal** | **39.236 GB** | **54.456 GB** | **91.742 GB** | yes |
| checksum replica | 39.236 GB | 54.456 GB | 91.742 GB | yes |
| indexes/intermediates scratch | 10 GB | 30 GB | 100 GB | yes |
| streaming conversion scratch | 2 GB | 5 GB | 15 GB | yes |
| bounded source-chunk staging | 2 GB | 5 GB | 15 GB | yes |
| pre-margin acquisition capacity | 92.471 GB | 148.912 GB | 313.484 GB | yes |
| 25% safety margin | 23.118 GB | 37.228 GB | 78.371 GB | yes |
| **acquisition capacity** | **115.589 GB** | **186.140 GB** | **391.855 GB** | yes |
| optional all-set plain FASTA | 122.353 GB | 137.511 GB | 153.333 GB | **no; stream instead** |
| later compute scratch | 0.5 TB | 2 TB | 8 TB | **no; separate scope/path** |

Arithmetic is explicit: primary durable is BGZF + retained GFF3 + `.fai` +
`.gzi` + metadata.  Pre-margin capacity is two primary copies +
indexes/intermediates + conversion scratch + source staging.  Acquisition
capacity is 1.25 times that subtotal.  Source transfer is not added as a full
resident copy because acquisition is chunked; its bounded staging peak is.
The optional plain-FASTA row is also excluded because source gzip/plain FASTA
can be streamed one file at a time into BGZF.

The file/inode model is 78,331 low (three files per assembly plus 100), 104,808
central (four per assembly plus 500), and 157,462 high (six per assembly plus
1,000).  This includes canonical FASTA, `.fai`, `.gzi`, optional GFF3, and
allowances for metadata/source records; it is intentionally much larger than
current manifest counts.

`bgzip`/htslib documents BGZF as concatenated gzip-compatible blocks smaller
than 64 KiB and `.gzi` as a count followed by pairs of 64-bit compressed and
uncompressed offsets [bgzip manual](https://www.htslib.org/doc/bgzip.html)
(accessed 2026-07-24T19:59Z).  `samtools faidx` supports BGZF FASTA and `.fai`
and `.gzi` index paths [samtools-faidx manual](https://www.htslib.org/doc/samtools-faidx.html)
(accessed 2026-07-24T19:59Z).  Locally installed evidence was bgzip/htslib 1.19
and samtools 1.19.2.  This supports the sidecar categories, but exact options,
compression level, headers, and names await the dedicated policy task.

### Transfer-duration scenarios

Transfer bytes are genomic FASTA plus GFF3: 39.116 / 50.633 / 83.446 GB.  For
an explicit sustained end-to-end rate `R` Mbps:

`hours = decimal_GB * 10^9 * 8 / (R * 10^6 * 3600)`.

| Sustained payload rate | low | central | high |
|---:|---:|---:|---:|
| 10 Mbps | 8.69 h | 11.25 h | 18.54 h |
| 50 Mbps | 1.74 h | 2.25 h | 3.71 h |
| 100 Mbps | 0.87 h | 1.13 h | 1.85 h |
| 500 Mbps | 0.17 h | 0.23 h | 0.37 h |

These are arithmetic transfer-time estimates, not wall-clock promises.  They
exclude API/package preparation, per-file latency, rate limiting, backoff,
validation, decompression/recompression, server-side changes, and retries.
With 26,077 small records, latency and packaging can dominate even when raw
bandwidth is high.

## Guardrails and layout

### Durable root

At every run and before every chunk, use `statvfs`/`df` against the actual
path, not a cached mount table.  Create the proposed directory as Erik with
0750/02750 permissions only after authorization; do not place payload in the
Git worktree.

* Start only when `/` has at least **2,400,000,000,000 available bytes** and
  **1,000,000 available inodes**.
* Stop before creating/downloading the next chunk if `/` has less than
  **2,000,000,000,000 available bytes** or 1,000,000 available inodes.
* Also stop if `available - next_chunk_peak < 2 TB`; checking only the current
  value can overrun the threshold during an atomic write.
* Current available bytes exceed preflight by 247,149,449,216 bytes and current
  available inodes exceed the guardrail by 446,717,990.  The durable path
  therefore passes now, but concurrent usage can invalidate that result.

Proposed durable logical layout (exact object filenames deferred):

```text
/home/erikg/phind-data/ecoli26k/v1/
  manifests/             # immutable input, resolution, object, and run snapshots
  source-metadata/       # Datasets catalog/reports and HTTP provenance
  objects/               # one atomically committed assembly directory per key
  annotations/           # retained annotation if authorized
  quarantine/            # checksum/status failures, never treated as complete
```

### Scratch and moveability

Use `/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1` through a runtime
`WORK_ROOT` variable.  Durable manifests must identify assemblies and logical
artifact roles, not embed scratch absolute paths.  Per-run logs may record the
resolved `WORK_ROOT`; changing it later then does not change immutable object
manifests.

* Scratch preflight: at least **4 TB available** and **5,000,000 available
  inodes**.
* Stop dispatching new jobs below **2 TB available** or 5,000,000 available
  inodes, and reserve each job's measured peak before launch.
* NVMe 3 currently passes the central guardrail by 1.507 TB and has 1.492
  billion free inodes.  The placeholder 8 TB high compute scenario would need
  at least 10 TB including the 2 TB reserve; current free space is short by
  about **4.493 TB**.  Pilot evidence or another scratch tier is mandatory for
  that case.

Keep a verified replica on NVMe 2 only after a primary object is committed.
Because NVMe 2 and root are separate devices but the same host/chassis, later
arrange an off-host or managed durable backup if recovery from host loss is a
requirement.

## Later acquisition and resume procedure

This procedure is a design for a separately authorized worker.

1. **Freeze inputs.** Copy the exact accession bytes and their SHA-256 into a
   new immutable manifest version.  Reject/quarantine line 26,078 rather than
   interpreting it as an accession.  Never silently replace a requested
   accession.version with a newer version.
2. **Refresh metadata.** Re-run exact-version validity/status metadata and
   archive the raw JSON, response timestamp, API version/headers, request
   batches, and resolution status.  Record latest/replaced/suppressed and
   missing annotation explicitly.
3. **Plan deterministic chunks.** Preserve input order and use fixed chunk IDs.
   Start with a small pilot, then keep packages below both a chosen accession
   count (for example 100–500) and a measured 5 GB staging cap.  NCBI's large
   package guidance supports dehydrated metadata followed by rehydration; pin
   and record the later-installed Datasets CLI version.
4. **Preflight capacity/inodes.** Apply the next-chunk-aware thresholds above
   to durable, scratch, and replica filesystems.  Confirm expected devices with
   `findmnt`, so an unmounted path cannot accidentally write to root.
5. **Rate limit.** The API response at 2026-07-24T19:51:54Z advertised
   `x-ratelimit-limit: 5`; without an API key, target at most **3 metadata
   requests/s** and **2 concurrent payload transfers**, reducing concurrency
   on 429/5xx or elevated errors.  Honor `Retry-After`.  Use exponential
   backoff with jitter (for example 1, 2, 4, 8, 16, 32, 60 seconds), at most
   eight attempts per object, then quarantine for operator review.
6. **Resume downloads safely.** Write a uniquely named `.partial` on the same
   filesystem as its final staging object.  Persist URL, exact accession,
   expected size, ETag/Last-Modified when provided, byte count, and upstream
   catalog/checksum.  Resume a range only when remote identity metadata still
   agrees; otherwise discard/quarantine the partial and restart that one
   object.  Never append to a committed file.
7. **Verify source.** Validate archive/container structure, NCBI's decompressed
   MD5 entries, expected member accession, exact byte counts, and a locally
   computed SHA-256.  MD5 is retained for upstream compatibility; SHA-256 is
   the durable local identity.
8. **Stream conversion.** Feed an NCBI gzip or plain FASTA stream directly to
   the policy-approved `bgzip` invocation.  Write BGZF, `.gzi`, and `.fai`
   inside a uniquely named staging directory on the destination filesystem.
   Do not materialize the whole collection as plain FASTA.  Record source and
   canonical compressed SHA-256, canonical uncompressed-content digest,
   bytes, sequence/contig counts, headers/crosswalk, tool versions, command,
   and times.
9. **Validate and atomically commit.** Run BGZF integrity validation, build and
   test both indexes, confirm FASTA/manifest sequence counts and uniqueness,
   `fsync` files and staging directory, then rename the whole per-assembly
   staging directory into its final path.  An immutable manifest commit record
   written last is the completeness marker; multi-file existence alone is not.
10. **Idempotent resume.** On restart, skip only an object whose commit record,
    exact requested accession.version, sizes, and all checksums validate.
    Re-validate but do not overwrite committed objects.  Retry partial states
    independently; log every state transition (`planned`, `downloaded`,
    `source_verified`, `converted`, `indexed`, `committed`, `replicated`, or
    `quarantined`) in append-only JSONL and periodically seal a sorted TSV
    snapshot with SHA-256.
11. **Replicate and audit.** Copy content-addressed committed objects to the
    replica device, verify SHA-256 after copy, and periodically scrub both
    copies.  A new upstream release creates a new manifest version; it never
    mutates the old snapshot.

Illustrative single-object conversion shape only (not executed and not the
final naming/BGZF policy):

```bash
# Guard deliberately prevents accidental list-wide use.
test -n "${ONE_ACCESSION_ONLY:-}" && test "$ONE_ACCESSION_ONLY" = "$accession"
stage="$DURABLE_PARENT/.stage.${accession}.$(uuidgen)"
mkdir -m 0750 "$stage"
if test "$source_encoding" = gzip; then
  gzip -dc -- "$source_file"
else
  cat -- "$source_file"
fi | bgzip -@ "$threads" -c >"$stage/canonical.bgz.partial"
bgzip -t "$stage/canonical.bgz.partial"
bgzip -r "$stage/canonical.bgz.partial"       # illustrative .gzi creation
samtools faidx "$stage/canonical.bgz.partial" # illustrative .fai creation
# Verify/checksum/fsync, then rename the directory and append commit record.
```

Do not run this template against the accession file.  The dedicated PanSN/BGZF
task must replace the illustrative names/options and approve compatibility
before a pilot.

## Uncertainty and risks

* The ten-record size preview is small and chosen by position/version, not by
  genome size or assembly fragmentation; it does not bound rare outliers.
  Low/high per-assembly assumptions are deliberately wider than the preview.
* `download_summary` is server metadata, not a measured transfer or BGZF
  conversion.  Canonical BGZF could vary with compression level, block/text
  behavior, headers, and tool version.  A separately authorized 10-genome
  pilot should replace ratios before bulk authorization.
* All 26,077 accessions passed the validity endpoint on the access date, but a
  later worker must capture current/replaced/suppressed status and annotation
  availability.  The detailed-report omission for GCF_000167895.3 shows why
  check-validity alone is insufficient.
* Annotation is optional and format choice matters: this estimate retains
  compressed GFF3 only.  GBFF, proteins, CDS, RNA, or multiple annotation
  formats would add transfer, files, and inodes not included here.
* Free space is changing on shared local devices.  Root is acceptable for the
  hundreds-of-GB acquisition only with the hard stop; it is not authorization
  for the 0.5–8 TB compute-scratch scenario.
* No observed filesystem provides evidenced RAID, remote durability, or
  snapshots.  Same-chassis copies do not cover host-wide loss.
* The bounded filename search can miss deeper, symlinked, generically named,
  or externally indexed reusable genomes.  Its zero exact-hit result applies
  only to the documented scope.
