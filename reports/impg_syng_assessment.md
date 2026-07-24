# Installed IMPG/SYNG capability assessment

**Evidence date:** 2026-07-24 UTC  
**Scope:** installed `/home/erikg/.cargo/bin/impg`, ordinary tiny synthetic
FASTA/SYNG behavior, and a design assessment for a 24,000–26,000 *E. coli*
workflow. This report does **not** inventory storage, acquire genomes, compute
current prophage statistics, design the host tree, or authorize production.
BGZF and literal PanSN `#` compatibility are owned by
`pansn-bgzip-genome-layout` and are deliberately not duplicated here.

## Decision summary

1. The installed program is **`impg` 0.4.1**, “implicit pangenome graph,” with
   an embedded **`syng` backend/subcommand**. “IMPG-SYNG” is understandable
   shorthand for that pairing, but is not the executable name. There is no
   installed or authoritative evidence for “INPG-SYNG” or “INFG-SYNG”; treat
   those as speech/transcription errors, not separate tools.
2. There are two different whole-cohort representations:
   * an alignment-backed IMPG index (`.impg`, optionally one per alignment
     file), which requires PAF/1ALN/TPA; and
   * one logical SYNG index prefix, physically six files, built directly from
     FASTA or AGC.
   They are alternative query backends, not two mandatory layers. For the
   proposed first pilot, one cohort-wide SYNG prefix is the simpler hypothesis;
   it is an **implicit syncmer graph**, not a stored base-complete GFA.
3. Known prophage **coordinates already on indexed paths** should use
   `impg query -b`. Extracted/novel prophage **sequences** should use
   `impg map`; its PAF is a syncmer-anchor coordinate projection, explicitly
   not a base-level alignment/CIGAR. Exact FASTA/GFA spelling still needs the
   original sequence collection.
4. IMPG does not directly create prophage biological clusters or a
   genome-by-cluster core/accessory matrix. Those are separate stages. A
   cohort-wide index can find candidate homologous intervals; cluster
   assignment, validated presence rules, and matrix aggregation remain
   downstream work. Unrelated prophage families should not be forced into one
   graph merely to obtain a matrix.
5. There is **no evidence that this exact binary/workflow has handled 26k
   bacterial assemblies**. A SYNG route avoids the explicit all-pairs count,
   but current `--parallel-dictionary` still replays GBWT paths serially and has
   no shard-merge or build-resume operation. A staged benchmark is mandatory.
6. “Hundreds of GB” is currently defensible only as a rough **canonical
   sequence-base** scale: `N × mean assembly bp` is about 120–130 GB at a
   stated 5 Mb planning assumption. It does **not** include SYNG sidecars,
   temporary/uncompressed staging, query products, regional GFAs, logs, or
   safety headroom. Total peak storage is unknown until the pilot measures
   expansion factors.

## 1. Exact installed identity and terminology

### Installed evidence

| Item | Observed value |
|---|---|
| Resolved command | `/home/erikg/.cargo/bin/impg` |
| File type | regular executable (not a symlink); ELF 64-bit x86-64 PIE, dynamically linked, not stripped |
| Size / mtime | 17,172,496 bytes; 2026-07-24 19:17:35 UTC |
| SHA-256 | `509296fb5c052be291a1841ea41f9bd4eb98e49b58b5f22cd69603729a94285f` |
| CLI identity | `impg 0.4.1` |
| Cargo install record | `impg 0.4.1 (path+file:///home/erikg/impg)`, release profile, no optional features, x86_64 Linux |
| Embedded toolchain comment | GCC 14.3.0; rustc 1.93.0 (`254b59607`, 2026-01-19) |
| Identical local build artifact | `/home/erikg/impg/target/guix/release/impg`, same SHA-256 |
| Current source checkout at inspection | clean `1bc9cf2aa9ddb46b9a1dd0cc8fbd2f094ab9de86`, described as `v0.4.1-13-g1bc9cf2` |

The installed binary does **not** embed a source Git SHA. It predates the current
checkout commit timestamp, and the Cargo registry records only a path, package
version, profile, and Rust compiler. Therefore the precise source commit used to
produce these bytes is **not recoverable** from installed metadata. The binary's
exact auditable build identity is the SHA-256 plus ELF/toolchain/package facts
above; claiming that it was built from current checkout SHA `1bc9cf2` would be
unjustified. The official source snapshot inspected for behavior was
`9138e906bad1fe03fe7e3435c5c75639c9d5c63b`, still package version 0.4.1
([Cargo.toml lines 1–4](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/Cargo.toml#L1-L4)).
Installed help and the probe are the authority wherever source applicability is
uncertain.

Full command output, dynamic libraries, statuses, and environment are in:

* [`artifacts/impg_probe/version.txt`](../artifacts/impg_probe/version.txt)
* [`artifacts/impg_probe/environment.txt`](../artifacts/impg_probe/environment.txt)
* [`artifacts/impg_probe/help.txt`](../artifacts/impg_probe/help.txt)

### Terminology resolution

* The upstream title is **“impg: implicit pangenome graph”** and describes an
  implicit graph made from pairwise alignments
  ([README lines 1–15](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/README.md#L1-L15)).
* Official documentation says IMPG **embeds `syng` as a second backend** that
  builds a syncmer GBWT from FASTA/AGC
  ([README lines 17–33](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/README.md#L17-L33)).
* Upstream SYNG calls its product “syng” and represents paths over a conceptual
  syncmer graph
  ([SYNG README lines 1–10](https://github.com/richarddurbin/syng/blob/28050a1270fb93458cca5091b5924ee73562776f/README.md#L1-L10)).
* Installed help contains commands `impg` and `impg syng`; `impg version` is an
  invalid subcommand (exit 2), while `impg --version` succeeds. Exact searches
  found no product named INPG or INFG. Thus use **IMPG**, **SYNG**, or **IMPG's
  SYNG backend**, not INPG-SYNG/INFG-SYNG.

### Installed command vocabulary

Installed top-level help reports:

`index`, `lace`, `partition`, `query`, `refine`, `similarity`, `genotype`,
`project`, `infer`, `stats`, `graph`, `normalize-self-loops`, `crush`,
`gfa2vcf`, `graph-report`, `render`, `align`, `map`, `read-index`, `syng`,
`syng2gfa`, and `syng-repair`.

`genotype` currently exposes `cos` (alias `cosigt`); top-level alias `gt` is
source-documented. Every installed command's exhaustive help and exit status is
captured in `help.txt` (24 help invocations exited 0; intentional `impg version`
exited 2).

## 2. Representations, formats, and semantics

### The two index families are not interchangeable

| Representation | Build input | Physical product | What it means | Main consumers |
|---|---|---|---|---|
| Alignment-backed IMPG | PAF, 1ALN, TPA, or list of those; PAF requires `=`/`X` CIGAR operations | one `.impg`, or per-input `.impg` files in `per-file` mode | cache-oblivious interval-tree projection over pairwise alignments; “implicit graph” | `query`, `partition`, `refine`, `similarity`, `stats` |
| SYNG-backed IMPG | exactly one FASTA (`-f`) or one AGC archive (`--agc`) per invocation | one **logical prefix**: `.1khash`, `.1gbwt`, `.names`, `.pstep`, `.spos`, `.meta` | canonical syncmer dictionary plus signed GBWT-like paths and coordinate sidecars | `query`, `partition`, `map`, `render`, `genotype`, `infer`, `syng2gfa` |
| Explicit sequence graph | FASTA/AGC plus generated or supplied PAF, or a queried SYNG region | GFA (and optionally VCF/report/render products) | materialized sequence graph; engine may be pggb/seqwish/poa/syng | graph-specific downstream tools |

Alignment formats and the `=`/`X` requirement are authoritative in
[README lines 72–85](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/README.md#L72-L85).
The SYNG prefix contents and their functions are specified in
[the version-matched design lines 135–160](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L135-L160).
All six files are one index and must be versioned/promoted together.

A SYNG node is a canonical syncmer; each input sequence is represented as a
forward and reverse-complement signed-node path. Inter-syncmer sequence is **not
stored in full** and must be fetched again from FASTA/AGC to emit exact FASTA or
GFA gap sequence
([design lines 31–54](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L31-L54)).
Therefore a SYNG index is not a substitute for the canonical assembly files.

### Accepted inputs and outputs relevant to this plan

| Command/stage | Accepted input documented by installed 0.4.1 | Output / caveat |
|---|---|---|
| `syng` | FASTA **or** AGC; one argument, not a FASTA-list option | six-file prefix; default total syncmer `k=63`, inner `s=8`, seed 7, position sample rate 256 |
| `index` / alignment queries | PAF, 1ALN, TPA, directly or text list | `.impg`; `auto` uses per-file at ≥100 alignment files |
| `query` | an alignment backend or SYNG prefix; `seq:start-end` or BED-like file | alignment backend: auto/BED/BEDPE/PAF/GFA/VCF/MAF/FASTA variants; SYNG backend: BED, BEDPE, FASTA, regional GBWT, GFA, VCF (not PAF) |
| `map` | SYNG prefix plus query FASTA **or FASTQ** | GAF, PAF-like, binary/text pack, or projection bundle |
| `graph` | FASTA/AGC file(s) or a sequence list; optional PAF | GFA through pggb, seqwish, poa, syng/syng-local engines |
| `syng2gfa` | SYNG prefix; sequence files optional | GFA 1.0 P-lines or 1.1 W-lines; raw overlap or blunt mode; absent source DNA is replaced by `N` |
| `partition` | either index family | BED; separate GFA/VCF/MAF/FASTA; graph outputs can create many files |
| `lace` | GFA or VCF files/list | combined GFA/VCF; this is a regional graph/VCF combiner, **not** a SYNG index merge |
| `project` | explicit GFA plus GAF | graph pack/projection evidence |
| `genotype cos` / `infer` | compatible SYNG or GFA feature-space evidence (`pack`/`proj`, ranges/partitions) | cosine-ranked candidate genotypes / interval calls; not a probabilistic imputation engine |
| GFA utilities | GFA | normalized/crushed GFA, VCF, topology report, render bundle |

The installed help is more exhaustive than this workflow table. Version-matched
source confirms that SYNG query outputs differ from map outputs
([design lines 257–277](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L257-L277))
and that map PAF is only a syncmer-anchor projection, not a base alignment
([design lines 502–510](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L502-L510)).

**Compression boundary:** version-matched documentation explicitly promises
native `.paf.gz` for alignment input
([README lines 263–282](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/README.md#L263-L282)).
Installed help only says “FASTA,” not BGZF. The delegated PanSN/BGZF task's tiny
matrix subsequently demonstrated two specific installed-binary boundaries:
`impg syng -f` read one BGZF FASTA with literal `#` path names, and `impg map -q`
read a BGZF FASTA and retained literal `#` in PAF query/target names. See the
cross-task reconciliation in section 8. This does not prove every
multi-file/sequence-retrieval boundary.

### Query versus map

* **Coordinate query:** `query -r PATH:start-end` resolves `PATH` exactly through
  `.names`, walks the source range through `.pstep`, locates shared syncmers via
  `.spos`, groups/chains target hits, and sorts output. Default SYNG query then
  performs boundary refinement; `--syng-raw` is debug-only and emits padded
  syncmer-resolution ranges. Algorithm detail is documented in
  [design lines 279–301](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L279-L301)
  and [328–359](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L328-L359).
* **Sequence map:** `map -q probes.fa/fq` does not require query names to exist in
  the panel. GAF records a syncmer-node walk; PAF-like output loads coordinate
  sidecars and reports candidate panel intervals. Pack/proj are read-support
  feature products, not presence matrices. Supported map modes and feature
  semantics are documented in
  [design lines 463–489](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L463-L489)
  and [512–550](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L512-L550).
* **Regional graph:** `query -o gfa` first finds intervals, then materializes a
  regional GFA; `query -o gbwt` creates a new regional six-file index. These are
  additional per-region products, not needed merely to get coordinates
  ([design lines 361–375](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L361-L375)).

## 3. Upstream contracts and identifier rules

### Required project contract

1. **Assembly sequence:** one canonical record per contig in FASTA/AGC. A
   per-assembly FASTA list is accepted by graph/query sequence retrieval, but
   `impg syng -f` itself accepts one FASTA path. A pilot-validated aggregation
   or AGC contract is required before a cohort-wide build.
2. **Identifier:** IMPG stores only the first whitespace-delimited FASTA header
   token. The same exact, case-sensitive token must join FASTA, `.fai`, GFF
   `seqid`, prophage contig key, interval BED, query output, and the canonical
   crosswalk. This behavior is explicitly documented
   ([design lines 73–88](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L73-L88)).
   Reject blank or duplicate primary tokens before build; duplicate behavior is
   not a supported identity scheme.
3. **Coordinates:** maintain a declared canonical internal convention of
   **0-based, half-open `[start,end)`**, which is the SYNG path-walk convention
   ([design lines 211–225](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L211-L225)).
   Convert GFF3/GenBank feature coordinates once, recording original start/end,
   convention, strand, contig circularity, and conversion method. Validate
   `0 ≤ start < end ≤ contig_length`.
4. **Colon rule:** range parsing splits on the **last** colon; paths that already
   contain coordinate-like colons still need a trailing `:start-end`
   ([README lines 405–410](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/README.md#L405-L410)).
5. **Annotation:** GFF/GenBank are **not direct SYNG inputs**. They supply
   interval metadata. Give `query` a converted BED-like interval file or give
   `map` extracted FASTA sequences.
6. **Ambiguity risk:** current version-matched source maps non-ACGT input bases
   to `A` during SYNG construction
   ([design lines 100–110](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L100-L110)).
   Bacterial assembly N-runs therefore require explicit pilot tests and QC;
   do not interpret anchors through long ambiguous runs as biological evidence.
7. **Manifest:** retain raw accession, assembly accession **with version**,
   source contig ID/version, canonical path ID, sample ID, canonical FASTA path,
   content checksum, GFF checksum, contig length/topology, and coordinate
   conversion provenance. SYNG files contain no sufficient accession crosswalk.

PanSN field policy, escaping, BGZF checksums, and the lossless crosswalk schema
are delegated to `pansn-bgzip-genome-layout`.

## 4. Separation of the four requested analytical stages

| Stage | Required inputs | IMPG operation | Index/graph count | Product and limitation |
|---|---|---|---|---|
| **A. Whole-genome cohort representation** | canonical assembly FASTA/AGC and unique path IDs | preferably pilot `impg syng`; alternative pairwise alignment + `impg index` or `impg graph` | one logical cohort SYNG prefix (six files), **or** one logical alignment backend; explicit GFA is optional | queryable implicit index; SYNG does not store all bases |
| **B1. Tagged interval query** | converted prophage BED keyed to indexed contig paths | `impg query -b` | reuse the one whole-cohort index | homologous BED/FASTA/optional regional GFA; must preserve source interval ID externally |
| **B2. Extracted sequence mapping** | prophage FASTA with stable feature IDs | `impg map -q ... -o paf` | reuse the same whole-cohort SYNG index | candidate genome/path coordinates; PAF-like syncmer projection needs independent coverage/identity validation |
| **C. Distinct prophage pangenomes/clusters** | QC'd extracted sequences and/or B-stage homologs | no direct IMPG clustering command; optional `impg graph` per validated cluster or regional `query -o gfa` | normally **multiple cluster-specific graphs** or no graph; not another whole-genome index | cluster membership and graph paths; external clustering/threshold logic required |
| **D. Core/accessory presence matrix** | canonical sample crosswalk, validated cluster membership, B-stage mappings, callable/ambiguous status | external aggregation; `similarity` can describe a region but does not emit the matrix | no extra IMPG index required | genome/sample × cluster table with `present/absent/uncallable`; core/accessory thresholds and denominator must be declared |

A single whole-genome SYNG index is therefore intended for coordinate lookup and
probe mapping. “Multiple indexes” only arise if choosing alignment per-file mode,
materializing a regional `gbwt` per locus, deliberately constructing separate
cluster indexes, or versioning successive cohort releases. Regional/clustering
products must not be confused with the canonical whole-genome index.

For matrix construction, collapse contig/path hits through the canonical
assembly/sample crosswalk. Define presence using pilot-validated minimum query
coverage, anchor count, coordinate consistency, and ambiguity handling. Count at
most one presence per biological sample/cluster; retain copy count separately.
“Core” must use a stated callable-sample denominator, not raw contig count.

## 5. Parallelism, determinism, restart, and failure behavior

### Parallelism and controls

* Global `-t/--threads` defaults to 4
  ([source lines 2026–2037](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/src/main.rs#L2026-L2037)).
* `syng --parallel-dictionary` parallelizes extraction and deterministic
  sort/dedup, but path insertion remains ordered/serial. The current limitation
  and absent GBWT merge API are explicit in
  [the construction design lines 54–88](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/syng-parallel-construction.md#L54-L88).
* SYNG position-sidecar repair and map workers use Rayon threads. Position
  sample rate trades disk/checkpoint count against at most more path walking.
  There is no SYNG memory cap, temp-directory flag, shard merge, or resume flag.
* Alignment/graph routes expose `--batch-bytes`, `--batch-size`, `--max-disk`,
  `--zstd`, `--zstd-level`, `--temp-dir`, seqwish `--disk-backed`, and
  `--transclose-batch`. `--max-disk` constrains alignment temporary planning,
  not total graph/index/output disk.
* Pair work can be selected/sparsified and emitted as a joblist; joblist
  execution has `--jobs`. Do not oversubscribe: effective runnable threads are
  approximately `jobs × threads_per_job`, plus I/O helpers.
* Graph pair controls include `--pairs`, `--pairs-done`, `--pairs-remaining`,
  `--max-pairs`, `--pair-start`, `--shuffle-pairs`, and `--shuffle-seed`.
  These can checkpoint pair generation/alignment, not SYNG construction or
  final graph induction.

### Determinism

* Parallel dictionary construction sorts/deduplicates all packed syncmers and is
  intended to avoid thread-order-dependent node IDs
  ([source lines 131–151](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/src/syng_parallel.rs#L131-L151)).
* The probe's two identical 2-thread parallel-dictionary builds were byte-for-byte
  identical across all six files. A legacy serial build was semantically usable
  but assigned different node IDs/bytes. Numeric node IDs must not be treated as
  stable across construction modes, parameter changes, or cohort versions.
* Source documents deterministic query sorting. This does not prove byte
  determinism for large graph engines, random sparsification, compressed
  outputs, or interrupted/resumed pair schedules. Use explicit shuffle seeds
  and test both byte and semantic equivalence in the pilot.

### Atomicity and restart

SYNG is only **partially atomic**. `.1gbwt`, `.1khash`, and `.names` are written
directly; `.pstep` and `.spos` use temporary-file rename; `.meta` is written last
([save implementation lines 2673–2762](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/src/syng.rs#L2673-L2762)).
A killed build can leave a mixed/partial prefix. There is no append/resume.
`syng-repair` can rebuild positional sidecars from intact core files; it cannot
resume a failed core build. Normal loads reject missing/stale position sidecars
([design lines 162–166](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/docs/designs/syng-integration.md#L162-L166)).

Production wrapper policy must build inside a same-filesystem
`cohort.release.partial.<run-id>/`, validate all six files and a sentinel query,
write checksums/manifest, then atomically rename the **directory** to a versioned
release. Never publish based only on `.meta` existing. Query/graph output files
should likewise use staging paths because several commands create/truncate final
paths directly.

Alignment-backed `per-file` indexing is incrementally rebuildable at file
boundaries and `auto` selects it at ≥100 files
([README lines 263–278](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/README.md#L263-L278)).
That feature does not merge/resume a SYNG prefix.

## 6. Scale assessment for 24k–26k assemblies

### Input arithmetic and all-pairs warning

Let `N` be assemblies, `L_i` their total contig bases, and
`B = Σ L_i`. At a clearly labeled planning assumption `mean(L)=5.0 Mb`:

* `B_24k = 24,000 × 5,000,000 = 120,000,000,000 bp` (120 GB of bases);
* `B_26k = 26,000 × 5,000,000 = 130,000,000,000 bp` (130 GB of bases).

This is not a measured collection size and excludes FASTA headers/wrapping.
Canonical BGZF byte size is a separate measured manifest value.

Explicit unordered all-pairs are:

* `C(24,000,2) = 287,988,000` genome pairs;
* `C(26,000,2) = 337,987,000` genome pairs.

Thus a naive pairwise whole-genome graph route is `Θ(N²)` alignments and can
create hundreds of millions of jobs/files. Installed help says default
`sparsify=none` means all pairs. Sparsification/joblists can reduce selected
edges and distribute work, but no installed evidence establishes graph
connectivity, biological recall, or final induction feasibility at 26k. Do not
launch the default whole-genome `impg graph` route.

### SYNG growth model

Let `M` be syncmer occurrences, `U` unique canonical syncmers, `P` contig paths,
and `r` the position sample rate.

* sequence scan and ordered path replay: approximately `O(B + M)` work;
* `--parallel-dictionary`: holds input sequence vectors and occurrence words,
  then sort/dedup is `O(M log M)` comparison work and `O(M)` occurrence memory
  before reducing to `U`; source code collects all extracted occurrences before
  sorting
  ([source lines 137–151](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/src/syng_parallel.rs#L137-L151));
* `.1khash`: data-dependent `O(U)`;
* `.1gbwt`: compressed, data-dependent `O(M)`;
* position checkpoints: approximately `O(M/r + P)` records, with both `.pstep`
  and `.spos` representations;
* `.names`: `O(P + identifier bytes)`.

This model exposes the likely bottlenecks: RAM for the parallel dictionary,
serial GBWT replay, compressed-index write bandwidth, path/contig count, and
checkpoint sorting/writes. Many small contigs increase `P`, fixed metadata, and
file-index work even when `B` matches a lower-path benchmark.

Upstream SYNG reports 92 human genomes / 277 Gbp producing 4 GB `.1khash` plus
5.8 GB `.1gbwt` in ~1.5 hours, but this is a different upstream executable,
hardware, path count, sidecar set, and dataset
([SYNG README lines 10–12](https://github.com/richarddurbin/syng/blob/28050a1270fb93458cca5091b5924ee73562776f/README.md#L10-L12)).
It is useful evidence that core syncmer/GBWT machinery can reach hundreds of
Gbp, **not** evidence that installed IMPG 0.4.1 can build/query 26k assemblies.
IMPG's own README only claims query behavior across “hundreds” of genomes
([README lines 7–15](https://github.com/pangenome/impg/blob/9138e906bad1fe03fe7e3435c5c75639c9d5c63b/README.md#L7-L15)).

### Disk, intermediates, and file count

Do not collapse storage into one “hundreds of GB” number. Measure:

`D_peak = D_canonical_BGZF + D_plain_or_AGC_stage + D_index_final +`
`D_index_partial + D_alignment_temp + D_query_products + D_regional_graphs +`
`D_logs + D_safety`.

For the recommended SYNG-first pilot, alignment terms can be zero, but canonical
BGZF may require decompression/aggregation depending on the delegated
compatibility result. During a staged rebuild both old and new six-file prefixes
may coexist. If every prophage BED row emits FASTA/GFA/VCF separately, output
file count can dominate metadata/inodes; batch tables and archive/version
regional products intentionally.

Define measured expansion factors for each pilot point:

* `E_index = D_index_final / B_input`;
* `E_partial = max_t(D_partial(t)) / B_input`;
* `E_query = D_query_products / B_queried_or_extracted`;
* `F_query = output_file_count / queried_interval_count`.

For a proposed production base total `B_prod`, fit `D = a + bB` across pilot
sizes and use the upper 95% prediction, not the 3.6 kb smoke ratio. Migration is
required when

`D_peak,95 > 0.70 × D_allocatable`

or projected output files exceed 50% of the allocated inode budget. During a
run, stop before the next batch if used bytes exceed 80% of the explicit project
allocation or if remaining allocation is less than twice the predicted
unfinished write volume. `D_allocatable` must come from the storage owner; the
machine-wide `df` value is not a project allocation.

## 7. Ordinary synthetic probe

The reproducible script is
[`artifacts/impg_probe/commands.sh`](../artifacts/impg_probe/commands.sh).
It generated exactly four ordinary-name sequences: three 1,200 bp panel records
and one 360 bp query, all deterministic synthetic DNA. It did not read production
lists, test BGZF, or use `#`.

| Probe item | Result |
|---|---|
| Serial SYNG build, 1 thread | exit 0; 3 paths, 3,600 bp, 315 syncmer occurrences; 12,288 KiB peak RSS |
| Parallel-dictionary build, 2 threads | exit 0; 139 unique nodes; 12,288 KiB peak RSS |
| Exact repeat of parallel build | exit 0; all six files byte-identical |
| Raw coordinate query `genomeA:180-540` | exit 0; 3 BED rows (one per homologous path) |
| Query-sequence map | exit 0; 1 GAF walk and 10 capped PAF-like candidate rows |
| Whole-index raw `syng2gfa` | exit 0; 141 segments, 153 links, 3 paths |
| FASTA interval output | exit 0; 3 records |
| Total retained synthetic files | <62 KiB (61,259 bytes after the recorded run); entire probe directory ~308 KiB including help/logs |

The parallel six-file prefix is 6,926 bytes / 3,600 panel bases = 1.92 bytes per
input base in this tiny, fixed-overhead, highly similar example. The raw GFA is
9,085 bytes / 3,600 bases = 2.52 bytes/base. **Neither ratio is a production
estimate.** The probe is far too small to model unique-syncmer growth,
checkpoint compression, contig counts, allocator behavior, or temporary files.
It only proves ordinary FASTA build/load/query/map/GFA plumbing in this binary.

`--syng-raw` intentionally produced broad padded intervals (for example
`genomeA 63 648` for source `180–540`); it did not test default BiWFA boundary
precision. It also did not test GFF/GenBank conversion, ambiguity/N runs,
interrupt/restart, clustering, matrices, graph smoothing, production-scale
parallelism, BGZF, or literal `#`. Those capabilities remain unknown until the
specified gates.

## 8. BGZF and PanSN delegation/reconciliation

`pansn-bgzip-genome-layout` owns the canonical compatibility matrix for:

* plain FASTA vs ordinary gzip vs BGZF;
* `.gzi`/`.fai` interactions and streaming/decompression policy;
* literal `sample#haplotype#contig` preservation through IMPG build, names,
  query, map, GFA, and other tools.

This assessment did not duplicate that matrix. Reconciliation with the other
owner's completed tiny run is:

* Cross-task artifact
  `artifacts/pansn_bgzip_probe/synthetic/results/status_matrix.tsv` records
  PASS/exit 0 for `impg_syng_bgzf_hash`, `impg_map_bgzf_hash`, and its
  verification case.
* The SYNG build log shows installed `impg` reading the canonical `.fa.gz`
  (BGZF), indexing two literal PanSN-like paths, and writing them unchanged to
  `.syng.names`; the map result retains `#` in both PAF query and target names.
  Evidence paths are
  `artifacts/pansn_bgzip_probe/synthetic/logs/impg_syng_bgzf_hash.stderr.txt`,
  `.../impg_index/example_a.syng.names`, and
  `.../results/impg_map_bgzf_hash.paf`.
* Applicability is deliberately narrow: one BGZF panel file at `syng -f` and
  one BGZF query file at `map -q`. The matrix does not establish a direct
  per-assembly file-list input to `syng` (none exists), concatenated BGZF,
  default `query --sequence-files/--sequence-list` boundary refinement/FASTA
  retrieval, or all graph engines. Those remain Gate 0 pilot checks.

Accordingly, use canonical BGZF directly only at those passed boundaries. At an
untested boundary, retain canonical BGZF and stream/decompress into bounded
staging or use a tested AGC aggregation; never rename identifiers during
transformation. A `#` range must still be shell-quoted. The naming policy and
full matrix authority remain `reports/pansn_bgzip_naming.md` and its probe
artifacts.

## 9. Non-executed production templates — **DO NOT RUN**

These are templates, not authorization. Variables must resolve to a versioned
manifest and pass section 10. They were not executed.

```bash
# DO NOT RUN — whole-cohort SYNG build after BGZF/input and RAM gates.
# Choose exactly one pilot-validated input form.
impg syng \
  -f "${PILOT_VALIDATED_SINGLE_PANEL_FASTA}" \
  -o "${STAGING_DIR}/ecoli.${RELEASE}.syng" \
  --syncmer-length 63 --smer-length 8 --syncmer-seed 7 \
  --position-sample-rate 256 --parallel-dictionary \
  -t "${THREADS}"

# Alternative only if the AGC contract/tooling is independently validated:
impg syng --agc "${PANEL_AGC}" \
  -o "${STAGING_DIR}/ecoli.${RELEASE}.syng" \
  --position-sample-rate 256 --parallel-dictionary -t "${THREADS}"
```

```bash
# DO NOT RUN — query tagged, converted 0-based half-open intervals.
impg query -a "${INDEX_DIR}/ecoli.${RELEASE}.syng" \
  -b "${PROPHAGE_BED}" -d "${VALIDATED_MERGE_DISTANCE}" \
  -o bed -t "${THREADS}" > "${STAGING_DIR}/prophage_homologs.bed.partial"

# DO NOT RUN — extract query-selected sequence only where canonical files passed
# the BGZF/retrieval gate.
impg query -a "${INDEX_DIR}/ecoli.${RELEASE}.syng" \
  -b "${PROPHAGE_BED}" -d "${VALIDATED_MERGE_DISTANCE}" \
  -o fasta --sequence-list "${CANONICAL_FASTA_LIST}" \
  -t "${THREADS}" > "${STAGING_DIR}/prophage_homologs.fa.partial"
```

```bash
# DO NOT RUN — map extracted/novel prophage sequences to the same whole index.
impg map -a "${INDEX_DIR}/ecoli.${RELEASE}.syng" \
  -q "${EXTRACTED_PROPHAGE_FASTA}" -o paf \
  --min-anchors "${PILOT_MIN_ANCHORS}" --max-hits "${PILOT_MAX_HITS}" \
  -t "${THREADS}" -O "${STAGING_DIR}/prophage_to_genomes.paf"
```

```bash
# DO NOT RUN — optional graph per already-validated prophage cluster, not the
# clustering operation itself.
impg graph --sequence-files "${ONE_CLUSTER_FASTA}" \
  --gfa-engine "${PILOT_VALIDATED_ENGINE}" \
  --batch-bytes "${PILOT_BATCH_BYTES}" --max-disk "${PILOT_TEMP_BUDGET}" \
  --temp-dir "${SCRATCH_DIR}" -t "${THREADS}" \
  -g "${STAGING_DIR}/${CLUSTER_ID}.gfa.partial"
```

```bash
# DO NOT RUN — alignment-backed explicit whole graph is blocked unless a
# sparsification/connectivity pilot succeeds. Never omit --sparsify here.
impg align --sequence-list "${CANONICAL_FASTA_LIST}" \
  --format joblist --sparsify "${PILOT_VALIDATED_SPARSE_STRATEGY}" \
  --output-dir "${PAIR_PAF_DIR}" -t "${THREADS_PER_JOB}" \
  > "${STAGING_DIR}/pair_jobs.sh.partial"
# Pair execution/validation/merge is a separate resumable scheduler stage.
impg graph --sequence-list "${CANONICAL_FASTA_LIST}" \
  --paf-file "${VALIDATED_MERGED_PAF}" \
  --gfa-engine "${PILOT_VALIDATED_ENGINE}" \
  --disk-backed --temp-dir "${SCRATCH_DIR}" -t "${THREADS}" \
  -g "${STAGING_DIR}/whole.gfa.partial"
```

The final presence matrix is not an IMPG CLI template: aggregate validated PAF
or BED evidence through the canonical crosswalk with explicit thresholds,
callable/uncallable status, and sample-level denominators.

## 10. Bounded pilot and hard pass/fail gates — **DO NOT LAUNCH HERE**

### Design

Use a frozen, non-production pilot manifest stratified by assembly size,
contig count, N content, phylogenetic diversity, and prophage length. Suggested
stages are 100, 250, 500, then at most 1,000 assemblies (roughly ≤5 Gbp at the
planning mean); stop at the first failed gate. Preserve one fixed 50-assembly
validation subset at every stage. Include known source intervals, exact extracted
source sequences, mutated positives, unrelated negatives, multi-copy/repeat
cases, contig-edge/circular cases, and ambiguous-base cases. Gate 0 imports the
PanSN/BGZF compatibility verdict.

At every stage capture `/usr/bin/time -v`, elapsed phase logs, CPU utilization,
peak RSS/swap, per-file and peak partial disk, input/output checksums, file count,
inodes, syncmer `M/U`, path count, query latency distribution, and result counts.
Fit `T(B)=cB^p` and `D(B)=a+bB` with upper prediction intervals; do not extrapolate
from the tiny smoke probe.

### Gates

| Gate | PASS | FAIL / action |
|---|---|---|
| **0. Input compatibility** | PanSN/BGZF task shows exact ID round-trip and readable chosen input at every needed IMPG boundary | block; use bounded plain/AGC staging and rerun |
| **1. Manifest/IDs** | 100% unique primary tokens; 100% FASTA↔GFF↔prophage↔crosswalk joins; all coordinates in range; checksums frozen | any collision, orphan, lossy rename, or coordinate error: stop |
| **2. Index integrity** | six expected files; nonzero and checksummed; fresh process loads prefix and sentinel queries every pilot path class; no stale sidecar warning | any missing/stale/mixed file or load error: stop |
| **3. Correctness** | exact-source controls recover originating path and overlap ≥95% of source interval; fetched source spelling matches canonical sequence for reported coordinates; strand correct in 100% of controls | any systematic coordinate/strand/sequence mismatch; investigate N/edge cases |
| **4. Mapping coverage** | ≥95% of non-ambiguous positive probes have a hit covering ≥80% of query span to the originating assembly; all negative controls lack hits meeting that same rule | threshold not met or negative false positives exceed 1%; tune only on training controls, then re-test held-out controls |
| **5. Query performance** | p95 single-interval and batch latency meet a predeclared operational SLA; `T(B)` upper 95% prediction at production is within the allocated wall-time window; exponent upper bound `p ≤ 1.3` | superlinear growth beyond gate or production prediction exceeds window: stop/shard query batches or change backend |
| **6. RAM** | peak RSS ≤70% of the explicitly assigned RAM and zero OOM/swap growth; 1- vs many-thread result semantics agree | exceed cap or rising RSS per repeated batch: stop; prefer serial build/AGC or migrate |
| **7. Disk and file count** | `D_peak,95 ≤0.70×D_allocatable`; projected files ≤50% inode allocation; remaining allocation ≥2× unfinished predicted writes | stop before next stage and migrate/reduce products |
| **8. Determinism** | two identical `--parallel-dictionary` builds, including 1-thread and intended-thread runs, are byte-identical for six files; query/map outputs are byte-identical after declared sorting | differences: freeze mode/thread count and use semantic checks only after explaining provenance |
| **9. Interrupt/restart** | forced termination at ~25%, 50%, and sidecar-write phase never publishes a prefix; wrapper detects/cleans partial directory and a clean restart reproduces output; pair-job route skips validated completed pairs without duplicates | any mixed release, silent reuse, or corrupt output: wrapper is not production-ready |
| **10. Scale trend** | per-base index/partial factors stabilize (last two slopes differ ≤25%), no unexplained `M/U`, latency, or output-count jump | do not extrapolate; add an intermediate stage or change strategy |

Production authorization requires all gates, storage-owner confirmation of
`D_allocatable`, a stated runtime SLA, and a versioned immutable manifest. The
pilot may demonstrate that the current filesystem is adequate, but this report
cannot conclude that merely because canonical data are “hundreds of GB.”

## Citation register and applicability

All URLs were accessed **2026-07-24**.

* IMPG source/docs use official repository commit
  [`9138e906bad1fe03fe7e3435c5c75639c9d5c63b`](https://github.com/pangenome/impg/tree/9138e906bad1fe03fe7e3435c5c75639c9d5c63b),
  whose `Cargo.toml` is version 0.4.1. This snapshot matches installed command
  behavior observed in help/probe, but the binary lacks an embedded SHA; exact
  commit identity remains unknown as stated in section 1.
* SYNG context uses upstream commit
  [`28050a1270fb93458cca5091b5924ee73562776f`](https://github.com/richarddurbin/syng/tree/28050a1270fb93458cca5091b5924ee73562776f).
  It is architectural/benchmark context only; installed IMPG vendors a separate
  SYNG revision and adds IMPG sidecars and orchestration.
* Installed evidence and probe results are primary local evidence in the four
  required artifacts and `synthetic/probe.log`; every command includes merged
  stdout/stderr and `[exit_status=…]`.

### Known unknowns

* exact installed source commit and reproducible build recipe beyond captured
  Cargo/toolchain/local-artifact facts;
* BGZF/`#` passed delegated `syng -f` and `map -q` boundaries, but multi-file,
  query-retrieval, and graph-engine boundaries remain untested;
* full ambiguity/N-run impact, default boundary-refinement accuracy, and
  circular-contig semantics;
* 26k build time/RAM/disk, query hit explosion, and graph connectivity;
* byte determinism for large graph engines and scheduler resume;
* biologically valid prophage clustering and presence/core thresholds.

No bulk transfer, production indexing/query, BGZF conversion, host-tree build,
or full-data computation was performed.
