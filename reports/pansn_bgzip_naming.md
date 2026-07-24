# PanSN naming, identifier crosswalk, and canonical per-assembly BGZF policy

**Status:** authoritative policy for this collection's naming/crosswalk/BGZF boundary; not a collection-wide acquisition, storage-capacity, IMPG design, Mash phylogeny design, prophage-distribution analysis, or integrated execution plan.

**Evidence access date:** 2026-07-24 UTC.

**Policy version:** `pansn-bgzip-policy-v1`.

Only the first five nonblank accession lines and the prophage CSV header plus three rows were inspected to identify join fields. No real genome was downloaded, renamed, compressed, or indexed. The compatibility probe generated two fictional assemblies (four retained FASTA representations, 144,115 retained bytes) under `artifacts/pansn_bgzip_probe/synthetic/`.

## 1. Decision summary

1. A canonical FASTA sequence ID is exactly `sample#haplotype#contig`, with two literal `#` delimiters.
2. For an NCBI assembly, `sample` is the **resolved, versioned assembly accession**, including namespace and dot-version (for example, the symbolic form `GCF_<nine digits>.<version>`). It is not a strain label and not a BioSample accession.
3. Nominally haploid *E. coli* uses haplotype `1`. Chromosomes, plasmids, and unplaced/scaffold contigs of one assembly all remain under that same haplotype.
4. `contig` retains the source FASTA identifier token/accession.version when it is in the safe alphabet. Replicon role, plasmid name, topology, strain, isolate, and BioSample remain crosswalk metadata; they are not packed into the name.
5. The canonical per-assembly file is BGZF, named `<pansn_sample>.pansn.fa.gz`, plus `.fai` and `.gzi`. The storage task may choose its parent/sharding path; this policy fixes the basename and contents.
6. Renaming changes headers only. It must never reorder bases within a record, reverse-complement, trim, rotate a circular replicon, split/join contigs, or silently coerce coordinates.
7. Consumers may use an assembly as a phylogeny/IMPG/pangenome sample (`pansn_sample`) and a contig path as `pansn_sequence_name`. Human strain/isolate text is a display label only.

## 2. Published PanSN specification versus collection policy

### 2.1 Version-specific primary source

The current upstream release is **PanSN-spec v0.1.0**, tag and `main` both at commit [`166a7a4b1fe9ea691402f0fe421886ad8a8aeabc`](https://github.com/pangenome/PanSN-spec/tree/166a7a4b1fe9ea691402f0fe421886ad8a8aeabc), verified against `refs/heads/main` and accessed 2026-07-24. The repository has no later tag at access time.

The specification states the three-field pattern and field types—sample string, delimiter character, numeric haplotype, and contig/scaffold string—in [README lines 13–26](https://github.com/pangenome/PanSN-spec/blob/166a7a4b1fe9ea691402f0fe421886ad8a8aeabc/README.md#L13-L26). It recommends `#`, says a delimiter must not occur in sample/contig identifiers, recommends configurable delimiter support, and requires the sample/haplotype hierarchy to be unique across the analyzed pangenome in [lines 28–33](https://github.com/pangenome/PanSN-spec/blob/166a7a4b1fe9ea691402f0fe421886ad8a8aeabc/README.md#L28-L33). It also says PanSN is not a container for generic metadata and directs such data to a table/database in [lines 72–80](https://github.com/pangenome/PanSN-spec/blob/166a7a4b1fe9ea691402f0fe421886ad8a8aeabc/README.md#L72-L80).

The upstream document uses “suggest” rather than RFC 2119 keywords. For interoperability, this report treats the ordered three-field grammar, numeric haplotype, delimiter exclusion, and globally unique sample/haplotype hierarchy as the **published normative core**; `#` is an upstream recommendation that this project adopts.

### 2.2 What upstream does *not* specify

PanSN v0.1.0 does not choose a strain, BioSample, or assembly accession; define “number” more tightly; prescribe a safe alphabet/escaping/maximum length; define a haploid convention; distinguish chromosomes/plasmids/unplaced contigs; handle assembly revisions; define GFF lexical escaping; or require BGZF/checksums. Every rule below is therefore explicitly **collection policy**, not an upstream PanSN claim.

## 3. Collection identifier policy

### 3.1 Sample field: versioned assembly accession

NCBI describes assembly accession.version as a stable, unique identifier for a set of assembly sequence records and recommends it over an assembly name; NCBI GFF files likewise use accession.version for unambiguous seqids ([NCBI GFF3 documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/), accessed 2026-07-24). NCBI also states that sequence updates increment a version while metadata-only updates do not ([NCBI assembly data model](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/data-model/), accessed 2026-07-24).

| Candidate | Assessment | Policy |
|---|---|---|
| Strain/isolate label | Human-readable but submitter-supplied, non-unique, missing, unsafe, and mutable. Duplicate labels are expected. | **Reject as identifier**; retain raw bytes as metadata/display text. |
| BioSample accession | Stable identifier for biological source material, but does not identify a particular assembly revision; one material/sample may lead to multiple assemblies. It can also be absent from imported/local records. | **Metadata/grouping only**. Never a fallback assembly ID. |
| Unversioned assembly accession | Namespace-aware but resolves to a moving/latest version. | **Reject**. Resolution must pin the dot-version. |
| Resolved `GCA_….<v>` or `GCF_….<v>` | Globally names the exact assembly record/revision and is already safe ASCII. | **Primary sample field**. Preserve `GCA` versus `GCF`, underscore, zero padding, dot, and version exactly. |
| Locally minted identity | Needed only for data with no authoritative versioned assembly accession. | Use the full-digest fallback below; never use strain/BioSample as an implicit fallback. |

Production NCBI sample tokens must match `^GC[AF]_[0-9]{9}\.[1-9][0-9]*$`. A raw unversioned accession is resolved to an exact version and the request, response, timestamp, and resolver version are recorded before naming. Suppressed/replaced versions remain traceable with `assembly_status` and `supersedes_assembly_accession_version`; they are not silently relabeled as a newer version. If both paired GCA and GCF accessions exist, the requested/resolved namespace is the sample and the paired accession is metadata. Stripping `GCA/GCF` or the version is forbidden.

**Local deterministic fallback.** When no authoritative accession exists, compute a domain-separated SHA-256 over length-prefixed UTF-8 authority, raw source accession (or explicit `NO_ACCESSION`), declared source revision, source-record order/IDs, lengths, and uppercase sequence bytes; use `LOCAL_SHA256_<64 lowercase hex>`. Persist the complete preimage fields and sequence digest. Reuse is allowed only when those fields and content are byte-identical. A digest reuse with a different preimage is a hard collision error, not a suffixing opportunity. This gives deterministic, collision-checked identity; no claim of mathematical collision impossibility is made for SHA-256.

The sample choice means collection denominators are **assembly revisions**, not unique biological materials or strain strings. Any later deduplication/grouping by BioSample is explicit metadata logic and must not rewrite PanSN IDs.

### 3.2 Haplotype field

* Collection value for each nominally haploid assembly: `1`.
* Project grammar: `^[1-9][0-9]*$`; canonical decimal, no sign, whitespace, leading zero, or `0`.
* All chromosome, plasmid, and unplaced records in an assembly use the same `sample#1` prefix. A plasmid is not haplotype `2`.
* If a source explicitly represents multiple phased genome copies, allocate `1..N` only from authoritative phase metadata and record `haplotype_basis`. Mixed/ambiguous cultures or alternate loci are quarantined rather than guessed.

This is consistent with the pggb haploid tutorial, which assigns `1` to every haploid assembly ([pggb commit `4225c6ce…`, tutorial lines 51–64](https://github.com/pangenome/pggb/blob/4225c6ce010553ef353c5ea13805e38c2016b503/docs/rst/tutorials/divergence_estimation.rst#L51-L64), accessed 2026-07-24), but it remains this collection's policy rather than PanSN v0.1.0 syntax.

### 3.3 Contig field and replicons

1. Parse and preserve the complete source FASTA defline bytes separately. The source sequence identifier is the token after `>` through the first ASCII whitespace. Empty IDs or duplicate tokens within an assembly fail import.
2. Prefer the exact identifier present in the source FASTA being canonicalized. If it is a versioned sequence accession, retain its dot-version. Parse base/version into separate crosswalk fields without removing them from `pansn_contig`.
3. Store assembly-report sequence name, GenBank/RefSeq aliases, sequence role, assigned molecule, plasmid name, topology, and length separately. Never substitute `chromosome`, `plasmid`, or `unplaced` for a source accession.
4. A complete chromosome, each plasmid, every unplaced contig, and every scaffold in a multi-contig assembly receives one PanSN name. Record order may remain source order, but role is metadata and does not affect haplotype.
5. A new assembly dot-version produces a new `pansn_sample` and canonical file even when some contig accessions are unchanged. Old files/crosswalk rows are immutable.

### 3.4 Safe alphabet, reversible escaping, and length

The collection delimiter is literal `#`. Production sample IDs are already restricted as above. The canonical contig encoder is byte-injective:

* bytes in `[A-Za-z0-9._-]` pass unchanged;
* every other byte, including `%`, `#`, whitespace, slash, pipe, colon, shell metacharacters, and non-ASCII UTF-8 bytes, becomes uppercase `%HH`;
* decoding accepts only the canonical uppercase form and exactly reverses one layer; a raw `%23` becomes `%2523`, so it cannot collide with raw `#` encoded as `%23`;
* no trimming, case-folding, whitespace replacement, or lossy “sanitize” operation is allowed.

The full ID must contain exactly two literal `#`, nonempty fields, and at most **240 ASCII bytes**. FASTA itself has no dependable cross-tool identifier limit; 240 is a conservative project gate, not an upstream standard or guarantee. Do not truncate. If reversible percent expansion would exceed 240, use `CTGSHA256_<64 lowercase hex>` and set `contig_id_encoding=SHA256_ALIAS_V1`; the raw ID bytes remain in `source_fasta_id_token_b64`, making reversal a mandatory crosswalk lookup. A digest/preimage mismatch is a hard failure.

Before promotion enforce unique `pansn_sequence_name` across the collection, unique `(pansn_sample,pansn_haplotype,pansn_contig)`, unique FASTA IDs within a file, and one canonical path per versioned assembly. A duplicate input line is retained as a separate `INPUT` occurrence but resolves to the same assembly object; it does not create a second file. Duplicate strain names never collide because they do not participate in IDs.

### 3.5 Examples (all fictional)

`EXAMPLE_…` is a reserved probe/documentation namespace and is rejected by the production NCBI resolver. These examples assert syntax/policy behavior only and do not assert real accessions, BioSamples, strains, or biology.

| Case | Synthetic input/metadata | Valid canonical behavior |
|---|---|---|
| Chromosome | sample `EXAMPLE_ASM_A_v1`, source `EXAMPLE_CHROMOSOME_v1` | `EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1`; role=`chromosome` |
| Plasmid | same assembly, source `EXAMPLE_PLASMID_v1` | `EXAMPLE_ASM_A_v1#1#EXAMPLE_PLASMID_v1`; role=`plasmid`; still haplotype `1` |
| Unplaced/multicontig | `EXAMPLE_UNPLACED_01`, `EXAMPLE_UNPLACED_02` | two names under the same `EXAMPLE_ASM_B_v2#1#` prefix; role metadata differs as needed |
| Duplicate strain | assemblies A and B both display `EXAMPLE_DUPLICATE_STRAIN` | different assembly sample IDs; duplicate display value retained verbatim |
| Missing strain | strain field absent | sample remains assembly identity; `strain_raw_b64=.` |
| Unsafe source ID | raw bytes for `ctg#A/plasmid|β` | contig `ctg%23A%2Fplasmid%7C%CE%B2`; exactly two literal delimiters remain |
| Assembly revision | symbolic `GCF_<base>.1` replaced by `GCF_<base>.2` | distinct sample/file/tip IDs; `supersedes…` links revisions; no overwrite |

Structurally invalid/project-rejected examples:

| Name | Reason |
|---|---|
| `EXAMPLE_ASM_A_v1#EXAMPLE_CONTIG` | missing numeric haplotype/one delimiter |
| `EXAMPLE_ASM_A_v1#hapA#EXAMPLE_CONTIG` | haplotype is not numeric |
| `EXAMPLE_ASM_A_v1#0#EXAMPLE_CONTIG` or `#01#` | project requires positive canonical decimal |
| `EXAMPLE_ASM_A_v1#1#ctg#raw` | unescaped delimiter creates four fields |
| `EXAMPLE_ASM_A_v1#1#` | empty contig |
| `GCF_<base>#1#contig` | unversioned assembly is not a production sample (symbolic illustration) |
| a strain label such as `EXAMPLE_DUPLICATE_STRAIN#1#contig` | PanSN-shaped but rejected as a production collection identity |
| any ID over 240 bytes or with literal whitespace | project gate; alias/escape rather than truncate |

## 4. Lossless identifier crosswalk

The template is `artifacts/pansn_identifier_crosswalk_template.tsv`. Every data row currently has `record_status=EXAMPLE`; the values correspond only to the synthetic probe. Production generation must refuse `EXAMPLE` rows as inputs.

### 4.1 Record model and encoding

* `INPUT`: one row per physical input line/occurrence, including duplicates.
* `ASSEMBLY`: one row per unique resolved assembly revision.
* `CONTIG`: one row per source FASTA record; parent is its `ASSEMBLY`.
* `FEATURE_JOIN`: one row per annotation/prophage feature requiring a coordinate-bearing join; parent is its `CONTIG`.

`.` is the null token. All potentially unsafe or byte-significant raw values are standard Base64 of original bytes in `*_b64`; display columns are not the lossless authority. JSON is compact UTF-8 with sorted keys in production. `row_sha256` is SHA-256 of columns 1–95 joined by tabs plus a final LF. IDs and hashes are lowercase hex unless an external accession defines case.

### 4.2 Header dictionary

| Columns | Meaning |
|---|---|
| `schema_version`, `record_status`, `record_type`, `record_id`, `parent_record_id` | Schema/control fields, explicit example/production status, normalized row kind, deterministic row key, and parent join. |
| `input_source_file`, `input_line_number`, `input_raw_line_b64`, `input_raw_accession`, `input_occurrence_sha256` | Exact origin and occurrence of an accession line. The Base64 raw line is authoritative; occurrence hash includes its original LF policy as declared in provenance. |
| `resolution_status`, `requested_assembly_accession_version`, `resolved_assembly_accession_base`, `assembly_version`, `resolved_assembly_accession_version`, `paired_assembly_accession_version`, `assembly_status`, `supersedes_assembly_accession_version` | Request-to-resolution audit, parsed base/version, exact GCA/GCF identity, paired namespace alias, current/suppressed/replaced state, and revision edge. |
| `biosample_accession`, `strain_raw_b64`, `isolate_raw_b64` | Biological/display metadata. Never identity fallbacks. Raw Base64 preserves missing/unsafe/duplicate values without normalization. |
| `sample_id_basis`, `pansn_sample`, `haplotype_basis`, `pansn_haplotype` | Why the assembly/local token and haplotype were chosen, plus exact PanSN fields. |
| `source_fasta_uri`, `source_fasta_header_b64`, `source_fasta_id_token_b64`, `source_contig_id_display` | Source object, complete defline, exact source ID-token bytes, and non-authoritative display form. |
| `source_contig_accession`, `source_contig_version`, `source_contig_accession_version`, `genbank_contig_accession_version`, `refseq_contig_accession_version`, `assembly_report_sequence_name` | Parsed source sequence accession/version and assembly-report aliases needed to join FASTA, GFF, and assembly report without guessing. |
| `contig_id_encoding`, `pansn_contig`, `pansn_sequence_name` | Identity/percent/SHA-alias rule, encoded contig field, and exact `sample#haplotype#contig`. |
| `replicon_role`, `plasmid_name_raw_b64`, `topology`, `contig_length`, `contig_sequence_sha256` | Chromosome/plasmid/unplaced/scaffold role, lossless plasmid label, `linear/circular/unknown`, sequence length, and wrapping/header-independent uppercase sequence digest. |
| `canonical_bgzf_relpath`, `fasta_seqid` | Canonical BGZF location relative to the chosen storage root and literal FASTA ID (normally `pansn_sequence_name`). |
| `source_gff_file`, `source_gff_row`, `source_gff_seqid_lexical_b64`, `source_gff_seqid_decoded`, `canonical_gff_seqid_lexical`, `canonical_gff_seqid_decoded` | Original GFF origin and raw column 1, its semantic decoded ID, and the canonical GFF lexical/semantic PanSN forms. |
| `prophage_source_file`, `prophage_source_row`, `prophage_source_sha256`, `prophage_id_raw_b64`, `prophage_genome_key_raw_b64`, `prophage_contig_key_raw_b64`, `prophage_composite_locus_key_sha256` | Lossless provenance, input-file digest, and original keys for each prophage row. The composite digest is SHA-256 of ASCII `prophage-locus-v1`, one NUL byte, then each exact genome/scaffold/begin/end UTF-8 component prefixed by its unsigned 64-bit big-endian byte length; it accelerates joins but never replaces the components. |
| `source_coordinate_convention`, `source_begin_raw`, `source_begin_integer`, `source_end_raw`, `source_end_integer`, `source_strand_raw` | Declared source convention, untouched textual coordinate/strand values, and exact integer parses. Parsed integers do not imply a 0/1-based or closed/open convention. |
| `canonical_coordinate_convention`, `canonical_start_0based`, `canonical_end_0based_exclusive`, `canonical_intervals_0based_halfopen_json`, `canonical_strand`, `wraps_origin`, `touches_left_boundary`, `touches_right_boundary` | Normalized 0-based half-open representation, including disjoint intervals for circular origin crossing, strand, and physical sequence-boundary flags. |
| `phylogeny_tip_id`, `impg_sample_id`, `impg_path_id`, `pangenome_sample_id`, `pangenome_path_id` | Assembly-level tip/sample IDs (`pansn_sample`) and sequence/path IDs (`pansn_sequence_name`). These are explicit rather than assumed equal. |
| `source_blob_sha256`, `source_decompressed_sha256`, `canonical_fasta_content_sha256`, `canonical_bgzf_sha256`, `fai_sha256`, `gzi_sha256`, `source_gff_sha256` | Raw transfer, decompressed source, exact canonical FASTA bytes, BGZF bytes, both index files, and source annotation checksums. |
| `source_url`, `acquired_at_utc`, `source_etag`, `source_last_modified` | Retrieval/provenance fields; absence is explicit, not fabricated. |
| `rename_policy_version`, `bgzip_version`, `bgzip_threads`, `bgzip_level`, `samtools_version`, `transformation_command_sha256`, `provenance_json_b64`, `row_sha256` | Exact transformation policy/tool parameters, command/script digest, extensible lossless provenance, and row integrity. |

### 4.3 Coordinate and annotation invariants

The canonical feature convention is **0-based, half-open**. A source explicitly documented as 1-based closed maps `[begin,end]` to `[begin-1,end)`; no conversion occurs until the source convention is known. GFF3 v1.26 defines seqid escaping, 1-based closed coordinates, strand values, and circular-origin representation in [lines 24–58](https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md#L24-L58).

The supplied prophage CSV header exposes `genome`, `scaffold`, `begin`, `end`, and `prophage_id` but no strand column; the few inspected values do not establish whether coordinates are 0/1-based or closed/open. Therefore production prophage rows remain `source_coordinate_convention=UNRESOLVED` and `source_strand_raw=.` until producer documentation or a sentinel with known bases resolves the convention. Decimal-looking integers must be preserved as raw text and must parse as exact integers before normalization. Strand must never be inferred.

For circular records, retain source topology and rotation. Represent a wrap as ordered intervals such as `[[start,L],[0,end]]`, set `wraps_origin=true`, and retain the original raw values. `touches_*_boundary` refers to serialized sequence boundaries; it does not claim a biological edge on a circular molecule. Base-sequence digest, contig length, strand, interval length, and fetched boundary bases must agree before accepting a join.

GFF3's allowed unescaped seqid set does **not** include `#`; therefore a semantic PanSN seqid has lexical delimiters `%23` in a standards-compliant transformed GFF. Example: FASTA `EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1`, GFF lexical `EXAMPLE_ASM_A_v1%231%23EXAMPLE_CHROMOSOME_v1`, GFF decoded value equal to the FASTA ID. A `%` already present in a percent-encoded PanSN contig is escaped again as `%25` at the GFF layer. Preserve original GFF unchanged and record both layers. A downstream tool must pass a toy GFF/FASTA decode-and-coordinate pilot; raw string equality between `%23` and `#` must never be assumed.

## 5. Canonical BGZF workflow

HTSlib 1.19 documents BGZF as gzip-compatible blocks smaller than 64 KiB, stdin/stdout operation, `-@` threads, `-l` level, indexing, and `-t` integrity testing ([bgzip 1.19 manual](https://www.htslib.org/doc/1.19/bgzip.html), accessed 2026-07-24). Samtools 1.19 documents BGZF FASTA input and `faidx` index/query behavior ([samtools-faidx 1.19](https://www.htslib.org/doc/1.19/samtools-faidx.html), accessed 2026-07-24).

### 5.1 Canonical bytes and streaming conversion

* Basename: `<pansn_sample>.pansn.fa.gz`; sidecars append `.fai` and `.gzi`.
* FASTA output: LF line endings, no description after the canonical ID, uppercase sequence, deterministic source record order, 60 sequence bases per line, final LF. Ambiguous IUPAC bases are preserved under an explicitly versioned validation policy; no base replacement.
* Default compressor: installed/pinned `bgzip -@ 2 -l 6 --binary -c`. Threads may change only through an explicit job parameter and must be recorded; this report makes no fleet-capacity claim.
* Plain source: `pansn_rename ... < source.fa | bgzip -@ 2 -l 6 --binary -c > "$part"`.
* Ordinary gzip source: first `gzip -t -- source.fa.gz`, then `gzip -cd -- source.fa.gz | pansn_rename ... | bgzip ... > "$part"`.
* Use `set -o pipefail` and check every producer/renamer/compressor status. Never retain a routine uncompressed canonical copy. Scratch consists of the compressed `.part`, sidecars, small checksum/mapping files, and bounded logs.

`pansn_rename` here denotes a deterministic streaming implementation of Sections 3–4, not an unchecked text substitution. It emits the crosswalk and per-contig sequence digests while copying sequence bases. Input and output per-contig length/digest equality proves rename-only behavior.

### 5.2 Checksums and reproducibility semantics

Record distinct SHA-256 values:

1. `source_blob_sha256`: exact acquired bytes (plain, ordinary gzip, or BGZF).
2. `source_decompressed_sha256`: exact decompressed source FASTA bytes.
3. per-contig sequence digest and length: wrapping/header-independent base identity used to prove coordinate preservation.
4. `canonical_fasta_content_sha256`: exact `bgzip -cd canonical.fa.gz` bytes. This is the logical canonical-content checksum and changes when naming/wrapping policy changes.
5. `canonical_bgzf_sha256`: exact compressed bytes, for transfer/cache integrity only.
6. `.fai`, `.gzi`, mapping, annotation, command/script, and provenance checksums.

A BGZF byte checksum is **not** a stable biological/content identity. It may change with HTSlib, zlib/libdeflate, compression level, thread/block behavior, `--binary`, or platform even when decompressed content is identical. Pinning and recording these settings helps reproduce an attempt but does not justify requiring byte-identical recompression across environments. Compare canonical-content and sequence digests for semantic equality; compare the BGZF digest only to a specific stored artifact.

### 5.3 Attempt, validation, promotion, retry, and resume

1. Acquire/validate the source separately. On the destination filesystem create a collision-resistant hidden attempt basename that retains the compression suffix, e.g. `.<sample>.attempt-<uuid>.part.fa.gz`, with an attempt JSON containing source digest, policy, intended ID, tool hashes, and start time. Refuse symlinks and an unexpected existing final.
2. Stream conversion once into the part. A failed pipeline deletes its part/sidecars (or moves them to bounded diagnostic quarantine) and records the failure.
3. Validate before promotion:
   * `bgzip -t "$part"` succeeds (CRC/truncation/EOF check);
   * recompute `bgzip -cd "$part" | sha256sum`;
   * strict FASTA parser verifies exactly two `#`, field grammar, uniqueness, allowed length, legal bases, expected record count/order/lengths, and input/output per-contig sequence digests;
   * `samtools faidx "$part"` succeeds and creates nonempty `"$part.fai"` and `"$part.gzi"`;
   * `.fai` names/lengths exactly match crosswalk; quoted `samtools faidx "$part" 'sample#1#contig:1-N'` spot queries agree with streamed base digests;
   * all checksums and the crosswalk/manifest row validate.
4. Fsync part and sidecars. **Preferred when the storage layout assigns one directory per assembly:** stage the BGZF, sidecars, crosswalk fragment, and checksum-complete `COMPLETE` marker in a hidden sibling directory, fsync it, then atomically rename that directory to the final assembly directory and fsync the parent. If the selected layout cannot rename a whole directory, rename the three files on the same filesystem and publish the `COMPLETE` marker/manifest transaction **last**; consumers ignore files without it. There is no atomic three-file rename, so the final marker is the fallback set-level commit point.
5. A retry never appends to or resumes compression from an incomplete BGZF stream. Up to three policy-controlled attempts may restart from the validated source with backoff. A fully written part can resume at validation only if `bgzip -t`, source digest, intended name, policy/tool parameters, and attempt metadata all match. Otherwise remove/quarantine it and restart.
6. If the final BGZF is valid but an index is missing after a crash, regenerate only `.fai/.gzi` and republish the marker; do not recompress. If final content disagrees with the manifest, quarantine and require operator reconciliation—never overwrite silently. Cleanup removes stale attempt files only when no live attempt owns them and retains bounded failure metadata.

### 5.4 Required manifest fields

At minimum include: schema/policy version; state (`STAGING|COMPLETE|QUARANTINED`); attempt ID/timestamps/retry count; input occurrence and raw/resolved/paired assembly accessions; BioSample/strain/isolate metadata; sample/haplotype basis; source URL/ETag/Last-Modified and all source digests; every contig raw header/token, alias/accession.version, PanSN fields/full name, role/topology/length and sequence digest; canonical relative path; canonical content/BGZF/FAI/GZI digests and byte sizes; record count/total bases; bgzip executable hash/version/threads/level/`--binary`; samtools hash/version; rename script/command/config hashes; validation statuses/timestamps; acquisition host/run ID; and supersession/collision/quarantine reason. The TSV crosswalk supplies the joinable subset; the attempt manifest may add operational detail.

### 5.5 Tools that require plain FASTA

Do not replace the canonical BGZF. First seek a documented stream/stdin mode. If a tool needs a seekable plain path, materialize `bgzip -cd` into a quota-bounded, attempt-specific scratch file, verify its SHA-256 equals `canonical_fasta_content_sha256`, run the tool, and delete it with an EXIT trap. Record tool/version, decompression command, scratch path/peak bytes, content digest, and deletion status. A tool that merely accepts ordinary `.gz` has not thereby proven BGZF random-access or literal-`#` safety.

## 6. Bounded compatibility evidence

### 6.1 Probe inventory

Reproduction entry point: `artifacts/pansn_bgzip_probe/commands.sh`. Exact host/environment and binary versions/hashes are in `environment.txt` and `tool_versions.txt`. Every synthetic command, stdout, stderr, and exit status is under `synthetic/logs/`; assertions/results are under `synthetic/results/`.

The script creates only:

* fictional assembly A: chromosome + plasmid, initially plain FASTA;
* fictional assembly B: chromosome + unplaced contig, initially ordinary gzip;
* two canonical BGZF files.

It enforces at most five retained FASTA representations and 100 MiB; the recorded run has 2 assemblies, 4 FASTA files, and exactly 144,115 retained bytes. Content checksums before/after conversion match. `status_matrix.tsv` records all cases as exit 0/PASS.

### 6.2 Exact installed-tool matrix

| Boundary tested | Exact installed version | BGZF | Literal `#` | Result and scope |
|---|---|---:|---:|---|
| `bgzip -t`; stream plain and ordinary gzip to BGZF | bgzip/HTSlib 1.19, binary SHA in `tool_versions.txt` | yes | headers retained by byte/content hash | PASS for both fictional files. |
| `samtools faidx`, `.fai/.gzi`, quoted region fetch | samtools 1.19.2 using HTSlib 1.19 | yes | exact names in `.fai`; `sample#1#contig:11-40` fetch | PASS. This establishes this binary's FASTA/index/query boundary, not every samtools release. |
| wfmash target/query to PAF | `v0.24.1-24-gf64becc4`, binary SHA recorded | yes | exact query and target names in PAF columns 1/6 | PASS. This is the installed MashMap-family mapper. It does not prove classic Mash. |
| IMPG syng index and map to PAF | `impg 0.4.1`, binary SHA recorded | yes | exact target names in `.syng.names`; exact query/target names in PAF | PASS. This is an index/map input-output boundary only; no graph/pangenome was built. |
| classic `mash`, `mashmap`, `skani` | not installed | unknown | unknown | NOT TESTED; no transparent-compression claim. |
| `pggb`, `seqwish`, `odgi`, `vg`, `minigraph(-cactus)` | not installed | local behavior unknown | local behavior unknown | NOT TESTED. Current pggb docs do explicitly recommend bgzip + faidx + PanSN ([pggb commit `4225c6ce…` README lines 32–33](https://github.com/pangenome/pggb/blob/4225c6ce010553ef353c5ea13805e38c2016b503/README.md#L32-L33)), but that is not evidence for an absent local toolchain. |
| gene-oriented Panaroo/Roary/Pangraph/Panta candidates | not installed | unknown/not necessarily their direct input | GFF lexical behavior unknown | NOT TESTED; selection is outside this naming probe. |

wfmash's current upstream README separately documents PanSN `-Y '#'` grouping and BGZIP-indexed FASTA ([commit `e040aa10…` lines 61–64 and 123–142](https://github.com/waveygang/wfmash/blob/e040aa10e87cab44ed5a4db005e784be62b0bd21/README.md#L61-L64), [indexing lines](https://github.com/waveygang/wfmash/blob/e040aa10e87cab44ed5a4db005e784be62b0bd21/README.md#L123-L142)). The installed binary is a different exact build, so the local probe—not the current README—is the compatibility evidence.

### 6.3 Mandatory pilot gates for untested boundaries

Use the existing synthetic inputs only; a combined panel would be the fifth and last FASTA and total retained output must remain under 100 MiB.

**Classic Mash after pinning a version.** The phylogeny task owns sketch strategy; its delegated production boundary is one whole assembly file per sketch (no `-i`) and an external crosswalk tip equal to `pansn_sample`. First test that exact mode, then run a separate diagnostic `-i` sketch solely to observe whether literal contig names survive:

```bash
mash --version
mash sketch -o synthetic/results/mash_whole_a \
  synthetic/canonical/EXAMPLE_ASM_A_v1.pansn.fa.gz
mash info -t synthetic/results/mash_whole_a.msh \
  > synthetic/results/mash_whole_a.info.tsv
mash dist synthetic/results/mash_whole_a.msh \
  synthetic/results/mash_whole_a.msh > synthetic/results/mash_whole_a.self.tsv

# Compatibility diagnostic only; not the production sketch mode.
mash sketch -i -o synthetic/results/mash_ids_a \
  synthetic/canonical/EXAMPLE_ASM_A_v1.pansn.fa.gz
mash info -t synthetic/results/mash_ids_a.msh \
  > synthetic/results/mash_ids_a.info.tsv
grep -F 'EXAMPLE_ASM_A_v1#1#EXAMPLE_CHROMOSOME_v1' \
  synthetic/results/mash_ids_a.info.tsv
```

The production boundary passes only if BGZF is accepted, exit statuses are zero, exactly one whole-assembly sketch is emitted, and self-distance is zero; the phylogeny tip is assigned from the crosswalk rather than inferred from a contig header or file display name. The diagnostic passes only if every expected sequence name is exact. A diagnostic name failure prohibits Mash `-i` mode but does not by itself prove the whole-assembly mode changed bases; a whole-mode BGZF failure requires a checksum-verified ephemeral plain FASTA. Never silently change IDs.

**Selected graph toolchain after installation/authorization:** stream A+B into one BGZF panel, index it, run the smallest one-thread pggb/seqwish fixture, and list paths through every chosen consumer:

```bash
{ bgzip -cd synthetic/canonical/EXAMPLE_ASM_A_v1.pansn.fa.gz;
  bgzip -cd synthetic/canonical/EXAMPLE_ASM_B_v2.pansn.fa.gz; } |
  bgzip -@ 1 -l 6 --binary -c > synthetic/input/EXAMPLE_PANEL.pansn.fa.gz
samtools faidx synthetic/input/EXAMPLE_PANEL.pansn.fa.gz
pggb -i synthetic/input/EXAMPLE_PANEL.pansn.fa.gz \
  -o synthetic/results/pggb -t 1 -p 80 -s 1000 -n 2
odgi build -g synthetic/results/pggb/SELECTED_FINAL.gfa -o synthetic/results/pggb.og
odgi paths -i synthetic/results/pggb.og -L > synthetic/results/odgi.paths.txt
vg paths -x synthetic/results/pggb/SELECTED_FINAL.gfa -L > synthetic/results/vg.paths.txt
```

Replace `SELECTED_FINAL.gfa` only with the actual output reported by the pinned pggb version; do not glob an arbitrary intermediate. Pass only if input BGZF is accepted; every expected full name appears exactly, once per input path, in final GFA and both path listings; `#` is neither split nor stripped; sequences spelled by paths match source digests; and outputs remain within the cap. Failure blocks that tool/version and triggers an ephemeral plain-input or explicit translation design. This is a future compatibility gate, not authorization to build the production pangenome.

**GFF consumer gate:** create a two-feature toy GFF whose seqid lexical form uses `%23`, run the exact selected parser with assembly A BGZF, and require that it resolves the decoded PanSN ID, preserves 1-based GFF coordinates/strand, and round-trips to the same canonical 0-based intervals. Also test a negative literal/raw mismatch. A parser that compares `%23` and `#` as raw unequal strings requires a documented tool-specific alias file; it must not prompt lossy global renaming.

## 7. Release-blocking invariants

A canonical assembly is publishable only when: the versioned sample identity resolved; every source header/alias has a reversible crosswalk; PanSN and path uniqueness constraints pass; strain/BioSample never substitute for assembly identity; all sequence digests/lengths match before and after header conversion; coordinate convention is declared before any feature normalization; GFF/prophage original keys and raw coordinate/strand values remain present; topology/origin/edge flags are explicit; `bgzip -t`, content hash, `faidx`, `.fai`, `.gzi`, and quoted `#` region checks pass; and the complete marker/manifest is committed last. Untested tool/version boundaries remain blocked behind Section 6.3 rather than being assumed transparent.
