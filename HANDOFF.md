# Immediate handoff: *E. coli* tagged-prophage pangenomes

**Purpose.** This is an operational handoff for the phage-pangenome project. It is intentionally conservative: it records the project objective, decisions already made, a restartable workflow, and the active WorksGood (WG) research graph. It does **not** replace the evidence-backed reports that the active tasks will produce.

**Handoff state:** 2026-07-24 UTC, while the five research tasks in Section 8 were still active.

## Status language used here

| Label | Meaning |
|---|---|
| **OBSERVED** | A fact already seen in the repository, environment, or WG graph. |
| **DECIDED** | A project policy to carry forward unless explicitly changed with rationale. |
| **PROVISIONAL** | A working choice or estimate that must be tested or reconciled with an upstream report. |
| **PENDING / TBD** | Do not infer a value or behavior; the owning research task or pilot must determine it. |

## Safety boundary

- Do **not** alter `26k_ecoli_accession.txt` or `26k_prophage1.csv`.
- Do **not** download the cohort, run a 24k–26k all-vs-all distance calculation, construct a production host tree, recompress/index the production collection, or build/query a production IMPG pangenome yet.
- Production work requires completed research reports, the integrated plan, successful resource/correctness pilots, and explicit approval.
- Treat all exact cohort counts, input schemas, tagged-prophage rules, coordinate conventions, tool behavior, storage budgets, and biological thresholds as **PENDING** until their owning tasks report them. Do not promote progress-log snippets or partial worker output to findings.

# 1. Objective and scientific rationale

## 1.1 Objective

Find and characterize phage pangenomes represented by **tagged prophages** in the supplied *E. coli* cohort. The intended program is to:

1. obtain and normalize every resolvable cohort genome and its required metadata/annotation;
2. construct a host-genome phylogeny independently of prophage content; host clades must be supported by host-derived sequence rather than prophage presence or composition;
3. build and query the installed IMPG-SYNG whole-genome pangenome, subject to the installed-tool capability report and pilots;
4. map the supplied tagged-prophage records back to exact assembly contigs, validate their coordinates, and extract them losslessly;
5. cluster related prophages using sequence, gene-content, and/or synteny evidence appropriate to mosaic phages;
6. construct prophage pangenomes for biologically defensible phage clusters and, where supported, compare or stratify them by independently defined host clades;
7. estimate core, soft-core, and accessory modules with explicit denominators and detection uncertainty; and
8. identify **candidate** ancestral components using an explicit evolutionary model and uncertainty analysis rather than equating present-day frequency with ancestry.

The whole-genome graph, host tree, prophage clusters, and prophage-specific pangenomes are distinct products. One does not automatically substitute for another.

## 1.2 Terms that must remain distinct

| Concept | Operational meaning | Must not be confused with |
|---|---|---|
| **Host clade** | A group of bacterial isolates defined from host-genome evidence, with clade rules fixed without using prophage presence, phage cluster, or tagged status. | A group of similar phages found in those hosts. |
| **Phage cluster/clade** | A group of prophages supported by phage sequence, shared gene families, and/or conserved synteny. Because phages are mosaic and recombine, a network or multiple locus/module histories may be more defensible than one bifurcating “phage tree.” | A host clade, even if the phages occur mostly within one host clade. |
| **Host phylogeny** | A tree or relatedness structure based on host-derived genomic signal rather than prophage content. Mash-family whole-assembly distances may be a scalable screen, but clade claims require higher-fidelity host-core validation plus sensitivity to excluding prophage/mobile sequence and to recombination. | A prophage gene-content or synteny network. |
| **Prophage similarity analysis** | Pairwise or graph-based comparison of extracted prophages/modules by sequence, homology, gene content, and synteny. | Evidence for the bacterial species tree. |
| **Prevalence** | Observed fraction of eligible units containing a called module, with a stated denominator. | Evolutionary age or ancestral presence. |
| **Core / soft core** | Modules above predeclared prevalence/detection criteria within a precisely defined analysis set. Thresholds are **TBD** and require sensitivity analysis. | Proof that a module existed in an ancestor. |
| **Ancestral state** | A model-based estimate on an explicit tree/network with gain, loss, horizontal transfer, detection error, and topology uncertainty considered. | “Common today,” “core,” or “found in the deepest sampled branch.” |

# 2. Known inputs and environment

## 2.1 Repository inputs

| Path | Role | Current knowledge boundary |
|---|---|---|
| `26k_ecoli_accession.txt` | Supplied cohort/accession source and intended genome denominator/join reference. It seeds the immutable manifest and later metadata/assembly resolution. | **OBSERVED:** file exists at repository root. **PENDING:** physical/nonblank/unique counts, duplicates, malformed values, accession kinds, versions, and resolution outcomes. These belong to `storage-genome-acquisition-inventory`; do not guess them. |
| `26k_prophage1.csv` | Supplied prophage records from which defensible “tagged prophage” subset(s), genome/contig joins, and source intervals must be derived. | **OBSERVED:** file exists at repository root. **PENDING:** schema, row/locus counts, tag interpretation, coordinate convention, deduplication, overlaps, missingness, and accession/contig join behavior. These belong to `prophage-table-distribution`; do not guess them. |

Both inputs are immutable source evidence. Record content checksums before and after any local analysis, retain original row/line identity, and never silently rewrite or normalize them in place.

## 2.2 Installed pangenome executable

- **OBSERVED:** an executable is installed at `/home/erikg/.cargo/bin/impg`.
- **PENDING:** its exact version/build, whether and how the project term “IMPG-SYNG” maps to installed terminology, subcommands, accepted formats, graph/index/query semantics, deterministic behavior, resume behavior, resource controls, and production suitability. The owner is `impg-syng-capability`.
- Until that report is complete, do not write production commands from memory and do not assume one graph, multiple graphs, direct interval mapping, BGZF support, or literal PanSN `#` support.

## 2.3 Initial filesystem snapshot

The planning context already recorded the following **approximate initial snapshot**:

| Filesystem/location | Approximate free space observed | Other observation |
|---|---:|---|
| `/` | about 2.5 TB | Root/current filesystem. |
| `/mnt/nvme1n1` | about 3.2 TB | Erik-owned directories exist. |
| `/mnt/nvme2n1` | about 3.8 TB | Erik-owned directories exist. |
| `/mnt/nvme3n1` | about 5.1 TB | Erik-owned directories exist. |

These values are **OBSERVED but approximate**, not a timestamped capacity commitment. Exact free bytes, inodes, mount/device types, ownership/writeability, bounded existing-data inventory, recommended durable/scratch paths, acquisition size/time estimates, and safety reserves are **PENDING** from `storage-genome-acquisition-inventory`. Free space can change while other agents or users run jobs.

# 3. Decisions already made

1. **DECIDED — canonical compression and indexes.** Canonical per-assembly FASTAs should be BGZF/bgzip compressed. Produce `.fai` and `.gzi` where the relevant toolchain supports them. Consumers that do not accept BGZF must receive an explicit, bounded staged stream/copy rather than causing the canonical store to drift.
2. **DECIDED with provisional details — PanSN naming.** After authoritative verification, use PanSN sequence names, provisionally `sample#haplotype#contig`. For haploid bacteria the haplotype field is likely `1`, but the stable sample token, contig token, escaping rules, revision behavior, and exact haplotype policy are **PENDING** `pansn-bgzip-genome-layout`. Free-text strain/isolate names are metadata, not canonical IDs, until duplicate, missing, unsafe-character, and stability behavior is resolved.
3. **DECIDED — lossless identity and coordinates.** Preserve a crosswalk from each physical accession-file line and each original prophage-table row/locus through resolved assembly/version, source contig, renamed FASTA/GFF seqid, phylogeny tip, IMPG sample/path, extracted prophage, and phage cluster. Preserve raw coordinate fields and the source coordinate convention; every conversion must be explicit and reversible.
4. **DECIDED with a resource condition.** The current filesystems are acceptable for acquisition and a pilot if total canonical inputs remain in the **hundreds of GB**. This is not authorization to fill a filesystem: set a hard free-byte and free-inode floor before writing. Use an Erik-owned `/mnt/nvme*` location for larger graph/tree scratch if pilot-measured expansion warrants it. Exact paths and floors are **TBD** from the storage report and pilot.
5. **DECIDED — production gate.** Do not start production all-vs-all distances or production IMPG builds until pilot correctness and resource gates pass. At `n` genomes, Mash-family all-pairs work contains `n(n-1)/2` unordered pairs and is therefore O(n²):
   - `n = 24,000`: **287,988,000 pairs**;
   - `n = 26,000`: **337,987,000 pairs**.
6. **DECIDED — provenance and atomicity.** Immutable inputs plus versioned manifests are the source of truth. Write partial outputs under unique temporary names, validate them, and atomically promote them. Record tool/build, parameters, input/output hashes, timestamps, and status so reruns can skip validated products and safely replace only failed products.

# 4. Proposed end-to-end workflow

Each phase must have a versioned run ID and a machine-readable ledger. A phase consumes only validated outputs of its dependencies. Failed or interrupted work leaves no promoted output; resume checks hashes and completion markers rather than relying on filenames alone.

| # | Phase: inputs and dependencies | Operations | Outputs | Required checks and restart behavior |
|---:|---|---|---|---|
| 1 | **Immutable intake manifest.** Inputs: the two supplied files. | Record physical line/row identities, source paths, byte checksums, parsing status, and provenance without modifying sources. Define append-only correction/version policy. | Versioned intake manifest and raw-to-record keys. | All physical source records are accounted for, including blanks/parse failures. Rehash sources. On restart, reuse only if source checksums and manifest schema/version match. Depends on completed input audits. |
| 2 | **Metadata resolution and assembly-version selection.** Input: accession manifest. | Resolve identifiers to authoritative assembly accession/version, BioSample, strain/isolate metadata, component contigs, annotation availability, and replacement/suppression status. Apply a documented deterministic version policy. | Resolution table, frozen acquisition manifest, unresolved/ambiguous queue. | No silent one-to-many or many-to-one collapse; duplicate and obsolete/replaced accessions remain traceable. Cache source response/provenance. Restart unresolved/transient failures only; never silently upgrade a frozen version. Policy is **PENDING** storage/integration reports. |
| 3 | **Resumable acquisition.** Input: approved frozen manifest and resource budget. | Rate-limited retrieval of selected assembly FASTA and required annotation/metadata into per-object temporary paths; verify source size/checksum where available; atomically promote. | Immutable source objects, download receipts, status/error ledger. | Preflight free bytes/inodes and continuously enforce the hard floor. Retry transient failures with backoff; quarantine checksum/semantic failures. Resume partial transfers only when the source supports safe validation. Production launch remains unapproved. |
| 4 | **BGZF + PanSN canonicalization and synchronized GFF seqids.** Inputs: verified source FASTA/GFF and finalized naming policy. | Stream to canonical per-assembly BGZF without retaining unnecessary plain FASTA; rename each sequence through the crosswalk; apply exactly the same seqid mapping to GFF; create compatible `.fai`/`.gzi`; record source-content, canonical-content, and compressed-byte hashes separately. | Canonical `.fa.gz`, indexes, synchronized GFF, seqid crosswalk, conversion provenance. | Validate BGZF integrity; canonical decompressed bases must equal source bases except headers; index lengths/random access must agree; every GFF seqid maps exactly once. Use unique partials and atomic promotion. PanSN/tool details remain **PENDING** the naming task. |
| 5 | **Assembly and annotation QC.** Inputs: canonical assemblies and resolved metadata. | Measure completeness/contamination proxies, assembly length/contiguity, ambiguous bases, contig topology, taxonomic consistency, duplicates/near duplicates, annotation consistency, and known assembly-status flags. Freeze exclusion/flag rules before analysis. | QC table, inclusion set, exclusion/review reasons, strata for the pilot. | Every assembly has an explicit status; exclusions never delete source records. Recompute only when input, QC tool, or threshold version changes. Exact thresholds are **TBD**. |
| 6 | **Scalable host distance/tree plus higher-fidelity validation.** Inputs: QC-passing host assemblies; no prophage-derived labels. | Use the Mash-family design selected by `mash-phylogeny-design` for scalable screening/staging. Define host clades from non-prophage host genomic evidence, not prophage counts, sequence content, or clusters. Validate key topology/clades with a higher-fidelity host-core route; assess recombination, near duplicates, sketch parameters, representatives versus all genomes, and exclusion/masking of annotated mobile sequence; and recheck promoted clades after tagged intervals are validated in Phase 8. | Sketches, staged distances, host tree/relatedness result, support/stability summaries, host-clade assignments and method version. | Pilot O(n²) materialization, runtime, RAM, I/O, tree stability, and all-genome/representative mapping first. Whole-assembly Mash topology alone is not sufficient proof of host phylogeny or independence from prophage content. Restart at sketch/sample partitions with hash-validated shards. |
| 7 | **Whole-genome IMPG pilot, then gated build.** Inputs: canonical QC-passing genomes and installed capability/naming reports. | Test accepted inputs, IDs, graph/index build, query semantics, determinism, sharding/merge if any, and interruption/restart on the pilot. Only after approval, build the whole-genome pangenome design chosen by integration. | Pilot evidence; later, versioned whole-genome graph/index, sample/path map, build metrics and provenance. | Compare known pilot sequences/paths and queries against direct sequence truth; measure runtime, peak RAM, disk/files, expansion, failure atomicity, and resume. Number of graphs/indexes and production commands are **PENDING**. |
| 8 | **Prophage key and coordinate validation.** Inputs: audited prophage rows, canonical crosswalk, source and canonical assemblies. | Resolve each defensible tagged-subset row to one assembly version and source contig; determine source coordinate convention; validate bounds, length, strand, overlap/nesting, contig-edge, circular wrap, and duplicate-locus rules. Retain raw and normalized intervals. | Locus manifest, resolved/unresolved/anomalous queues, interval segments, coordinate conversion provenance. | Round-trip normalized intervals back to the exact raw representation where defined; compare extracted bases from source and canonical names; never clip or repair silently. Rules are **PENDING** the table audit. |
| 9 | **Prophage extraction, mapping, and whole-genome graph query.** Inputs: validated loci and versioned whole-genome index. | Extract exact interval sequence(s), preserving strand/topology and multi-segment circular cases; assign stable locus IDs and hashes; map/query sequences or coordinates by the installed IMPG-supported route. | Extracted prophage BGZF/FASTA as specified, locus-to-host/path mappings, query results, unmapped/ambiguous diagnostics. | Direct extraction must match independently reconstructed expected bases. Every input locus has a terminal status. Restart by locus/query shard; reuse only against the same graph build and sequence hash. |
| 10 | **Prophage clustering.** Inputs: validated extracted prophages and gene calls/annotations. | Compare nucleotide/protein homology, gene-family content, order/orientation, and synteny; account for fragments and mosaics. Use networks or multiple resolutions where a single tree is misleading. Freeze method/threshold version. | Prophage similarity matrices/networks, optional supported trees, stable cluster assignments and uncertainty/outlier labels. | Sensitivity across defensible thresholds and clustering routes; inspect incomplete/edge calls and reference phages. Host clade labels must not determine clusters. Restart from versioned feature/matrix shards. |
| 11 | **Cluster-specific and host-stratified prophage pangenomes.** Inputs: cluster assignments, extracted prophages, independently defined host clades. | Build a pangenome within each sufficiently supported phage cluster. Separately compare cluster prevalence/composition among host clades or analyze phage-cluster × host-clade cells where sample size allows. A bag of unrelated prophages from one host clade is not automatically one “phage pangenome.” | Per-cluster gene/module catalogs, synteny graphs, presence/absence/copy-number matrices, host-clade stratifications. | Declare units, denominators, minimum sample requirements, gene calling and homology parameters. Track fragments, split/fused genes, paralogs, missing calls, and uncertain cluster membership. |
| 12 | **Core, soft-core, and accessory estimates.** Inputs: versioned presence/copy matrices and callable-status masks. | Compute prevalence under predeclared strict-core/soft-core/accessory definitions for each analysis set; distinguish biological absence from unavailable/uncallable sequence. | Module prevalence tables with numerators, denominators, confidence/sensitivity summaries, and classifications. | Thresholds and denominators are explicit and varied in sensitivity analysis. Repeat after excluding low-quality/incomplete loci and with alternative homology/synteny rules. Exact cutoffs are **TBD**, not implied here. |
| 13 | **Candidate ancestral-component inference.** Inputs: module matrices, host tree uncertainty, phage relationships/networks, quality/detection model. | Fit or compare explicit gain/loss/transfer-aware ancestral-state models where assumptions are defensible; evaluate alternate host topologies, phage clusterings, and detection-error treatments. | Posterior/likelihood or otherwise calibrated ancestral-state estimates, uncertainty intervals/scores, alternative models, and a clearly labeled candidate list. | Never infer ancestry from prevalence alone. Report non-identifiability and sensitivity to HGT, recombination, mosaicism, sampling, and topology. Validate against external/reference phages where appropriate. |
| 14 | **Final delivery.** Inputs: all versioned products and ledgers. | Reconcile identifiers and produce analysis-ready matrices, trees, networks/graphs, reports, methods, software lockfiles, and provenance. | Final manifest/crosswalk, host tree, phage networks/clusters, pangenomes, matrices, QC/resource reports, and reproducible workflow description. | Referential-integrity audit from every final cell/tip/path back to immutable input; rerun spot checks; no orphan IDs; archive parameters/logs/checksums. |

**Restart contract across all phases:** tasks are idempotent at the manifest-row or deterministic shard level; statuses distinguish planned/running/validated/failed/quarantined; validation creates a completion record tied to input hashes and tool/parameter versions; reruns do not overwrite valid products from another version.

# 5. Data model and naming

## 5.1 Proposed manifest/crosswalk fields

Use one normalized relational model or a losslessly joinable set of TSV/Parquet tables rather than forcing all one-to-many relationships into one row. Stable keys should join these field groups:

| Group | Proposed fields |
|---|---|
| Record/provenance | `record_status`, `manifest_schema_version`, `manifest_run_id`, `input_file`, `input_file_sha256`, `input_line_number`, `raw_input_accession`, `normalized_input_accession`, `created_at`, `source_url_or_database`, `source_retrieved_at` |
| Assembly resolution | `assembly_accession`, `assembly_version`, `assembly_accession_version`, `assembly_resolution_status`, `assembly_replacement_status`, `biosample_accession`, `strain_label_raw`, `isolate_label_raw`, `metadata_source`, `metadata_record_hash` |
| Source sequence | `source_fasta_path`, `source_object_checksum`, `source_contig_id`, `source_contig_version`, `source_contig_length`, `source_contig_topology`, `source_contig_sequence_checksum` |
| **Provisional PanSN** | `pansn_policy_version`, `pansn_sample`, `pansn_haplotype`, `pansn_contig`, `pansn_full_name`, `name_escape_map`, `name_collision_status` |
| Canonical FASTA | `canonical_bgzf_path`, `canonical_content_checksum`, `canonical_compressed_checksum`, `fai_path`, `gzi_path`, `canonicalization_tool_version`, `canonicalization_status` |
| Annotation | `source_gff_path`, `source_gff_checksum`, `source_gff_seqid`, `canonical_gff_path`, `canonical_gff_seqid`, `gff_seqid_mapping_status` |
| Prophage source/locus | `prophage_table_file_checksum`, `prophage_table_row_id`, `prophage_id_raw`, `tagged_subset_rule_version`, `prophage_locus_id`, `prophage_genome_key_raw`, `prophage_contig_key_raw`, `source_start_raw`, `source_end_raw`, `source_strand_raw`, `source_coordinate_system`, `normalized_start`, `normalized_end`, `normalized_coordinate_system`, `wraps_origin`, `interval_segment_ids`, `coordinate_conversion_version`, `coordinate_qc_status` |
| Analysis identities | `phylogeny_tip_id`, `host_clade_id`, `host_clade_method_version`, `impg_build_id`, `impg_sample_id`, `impg_path_id`, `extracted_prophage_path`, `extracted_prophage_checksum`, `prophage_cluster_id`, `prophage_cluster_method_version` |
| QC/workflow | `assembly_qc_status`, `prophage_qc_status`, `inclusion_status`, `status_reason`, `current_phase`, `attempt_count`, `last_error`, `tool_and_parameter_record`, `validated_at` |

Required integrity constraints:

- `(input_file_sha256, input_line_number)` and `(prophage_table_file_checksum, prophage_table_row_id)` remain immutable source keys.
- Assembly and contig **versions** are retained; normalization never erases the original token.
- Each renamed FASTA header has an explicit old-header → PanSN mapping, and each GFF seqid uses the same mapping.
- Raw prophage coordinates are never overwritten. Normalized coordinates include the named coordinate system and conversion version. Wrapped/circular loci use explicit ordered segment records rather than an ambiguous start/end repair.
- `phylogeny_tip_id`, IMPG sample/path IDs, extracted-locus IDs, and cluster IDs are foreign keys, not independently rederived display labels.
- Human-readable strain names remain metadata unless the naming report proves a deterministic, collision-proof canonical policy.

## 5.2 PanSN status

**PROVISIONAL pending `pansn-bgzip-genome-layout`:** the expected shape is `sample#haplotype#contig`, likely with haplotype `1` for nominally haploid bacterial assemblies. The authoritative PanSN version, normative field semantics, stable sample basis, component/contig representation, escaping, assembly-revision behavior, length limits, and downstream handling of literal `#` are not settled in this handoff. Do not generate production names until that report and its compatibility probe are integrated.

# 6. Scientific caveats and required sensitivity analyses

| Caveat | Risk to result | Required response |
|---|---|---|
| Incomplete prophages and assembly breaks | Missing modules can look accessory; fragments can form false clusters. | Track completeness/callability, repeat prevalence after excluding or separately modeling incomplete/edge loci, and retain uncertainty. |
| Contig-edge and circular cases | A locus can be truncated or wrap the origin; naive intervals can extract the wrong sequence. | Record topology, source convention, strand and ordered segments; validate exact coordinate round trips and sequence hashes. |
| Duplicate, overlapping, and nested calls | Counts and denominators can be inflated; one biological locus can be represented multiple ways. | Predeclare locus/deduplication rules, preserve all source rows, and report sensitivity to competing defensible rules. |
| Gene-family thresholds | Identity/coverage parameters change cluster, core, and accessory calls. | Evaluate multiple documented thresholds and orthology methods; version all family assignments. |
| Gene splitting/fusion and annotation error | One module may appear as two genes, or vice versa. | Compare protein/domain and synteny evidence; represent split/fusion relationships rather than forcing one-to-one orthology. |
| Paralogs and copy number | Binary presence can hide duplication and ambiguous orthology. | Retain copy number and locus context; distinguish orthogroups from individual copies. |
| Synteny and rearrangement | Gene presence alone loses module structure; strict order can miss rearranged homologs. | Provide both gene-content and synteny-aware analyses and test order/orientation tolerance. |
| HGT, mosaicism, and recombination | Different phage modules have different histories; one bifurcating tree may be misleading. | Use networks/module histories where appropriate, compare sequence/content/synteny views, and avoid forcing unsupported monophyly. |
| Host recombination and mobile DNA | Whole-genome distances can distort the host tree. | Validate host clades with host-core methods, recombination-aware checks, and sensitivity to prophage/mobile regions and sketch parameters. |
| Assembly and annotation quality | Contamination, fragmentation, and inconsistent callers bias both host and phage results. | Apply versioned QC, callable masks, common reannotation where justified, and quality-stratified analyses. |
| Sampling bias and near duplicates | Overrepresented outbreaks/lineages can dominate prevalence and ancestral estimates. | Report raw and weighted/dereplicated sensitivity analyses while preserving mappings to all genomes. |
| Host-clade circularity | Defining host clades with prophage content guarantees apparent clade association. | Freeze host clades from host-derived evidence before testing prophage distribution; never use phage features to define them. |
| “Common equals ancestral” fallacy | High prevalence may reflect recent sweep/HGT; ancestral modules may be frequently lost. | Separate descriptive prevalence/core calls from model-based ancestral inference and report uncertainty/non-identifiability. |

Recommended validation is not yet an established result: compare alternative QC, homology, clustering, prevalence, topology, and gain/loss assumptions; map uncertainty through downstream analyses; and compare recovered clusters/modules against curated or otherwise appropriate external/reference phages. External references are validation evidence, not a license to assign taxonomy or ancestry without support.

# 7. Pilot and production gates

## 7.1 Small stratified pilot

Select and freeze a **small, bounded** pilot manifest only after the input audits. Exact `N_pilot` is **TBD**; do not fabricate it before knowing the strata and resource evidence. The pilot should intentionally cover, where the audited data support them:

- resolvable, ambiguous, versioned, duplicated, and failed accession joins;
- genomes with zero, one, and multiple candidate/tagged loci;
- typical and extreme assembly sizes/contiguity/QC states;
- chromosome, plasmid, unplaced/multicontig, and unusual source seqids;
- ordinary, overlapping/nested, strand-reversed, contig-edge, and circular/wrapped intervals;
- duplicate/missing/unsafe strain labels relevant to naming;
- near-duplicate and divergent host genomes; and
- complete-looking and fragmentary prophages spanning preliminary similarity groups.

The pilot manifest, rationale, expected outcomes, and resource cap must be approved before sequence acquisition. Pilot membership must not be chosen to make joins or tools look artificially successful.

## 7.2 Measurable gates

| Gate | Pilot measurement / pass condition | Still TBD before production |
|---|---|---|
| Accession/metadata joins | Every pilot input line is accounted for as exactly resolved, explicitly ambiguous, or explicitly unresolved; no silent collapse or version drift. | Allowed unresolved fraction and manual-resolution policy for production. |
| Identifier integrity | Every canonical FASTA/GFF seqid, phylogeny tip, IMPG sample/path, locus, extraction, and cluster foreign key resolves through the crosswalk; collision count is zero for promoted IDs. | Final PanSN tokens and per-tool transformations. |
| Coordinate round trip | For every promoted pilot locus, raw → normalized → raw conversion is exact where defined; independent source/canonical extraction yields identical expected bases and length; anomalies are quarantined, not clipped. | Source CSV convention and policy for ambiguous/invalid calls. |
| BGZF/index compatibility | BGZF integrity passes; decompressed content hash and random-access slices agree with source; `.fai`/`.gzi` are usable; synchronized GFF seqids all resolve; each actual consumer is tested with BGZF and literal PanSN `#`. | Consumer-specific pass/fallback matrix and tool versions. |
| IMPG graph/query correctness | Known pilot sequence/path membership and positive/negative queries match direct truth; sample/path IDs round-trip; deterministic behavior is measured. | Exact subcommands, one-versus-multiple graph design, acceptable query metrics. |
| Runtime and peak RAM | Capture wall/CPU time and peak RSS per phase and per input/base/graph unit under fixed concurrency; extrapolation formula and uncertainty are recorded. | Maximum allowed runtime/RAM and production concurrency based on reports and available hardware. |
| Disk/file expansion | Measure source, BGZF, indexes, annotation, temp, graph/tree, logs, backups, and inode growth separately; include failed-attempt scratch. | Allowed expansion factors and finalized durable/scratch placement. |
| Restart/resume | Deliberately interrupt safe pilot stages; rerun leaves no corrupt promoted file, does not duplicate records, and either resumes safely or restarts the smallest documented shard. Validated hashes are stable when determinism is promised. | Tool-specific checkpoint/sharding strategy where unsupported. |
| Host-tree stability | Compare selected sketch/parameters, all-genome versus representative placement, and higher-fidelity host-core validation; report clade/topology support and discordance. | Numeric stability/support threshold, rooting/outgroup, final clade rule. |
| Prophage clustering/pangenome | Reference/known pilot relationships, alternative homology/synteny thresholds, fragment handling, and presence matrices behave as expected; every call has a denominator/status. | Final thresholds, cluster minimum size, core/soft-core cutoffs. |
| Free-space stop | Before and during every write, projected remaining free bytes **and inodes** must stay above a hard reserve; abort cleanly before crossing it. Predicted peak plus reserve must fit the selected filesystem. | Exact byte/inode floors from timestamped storage inventory and approved concurrent workload. |

## 7.3 Production authorization

Production remains **NO-GO** until all of the following hold:

1. all five research reports and artifacts in Section 8 are complete and reconciled by `integrate-phage-pangenome-plan`;
2. stable acquisition, version, PanSN, coordinate, host-clade, graph/index, clustering, and denominator policies are recorded;
3. the stratified pilot passes correctness, compatibility, restart, and resource gates;
4. measured peak storage/RAM/runtime plus a documented reserve fit the chosen filesystems/hardware, including concurrent jobs and failure scratch;
5. the O(n²) host-distance plan has an approved materialization/staging strategy rather than an implicit full run; and
6. an authorized reviewer/user approves the frozen production manifest and thresholds.

Any failed correctness, identifier, coordinate, free-space, or restart gate is a stop condition. Exact evidence-dependent thresholds remain **TBD** rather than being invented here.

# 8. Current WG graph and receiving-agent instructions

## 8.1 Graph at handoff

The reviewed planning structure is:

```text
.quality-pass-phage-pangenome-planning (done)
    ├── storage-genome-acquisition-inventory
    ├── prophage-table-distribution
    ├── impg-syng-capability
    ├── mash-phylogeny-design
    └── pansn-bgzip-genome-layout
          \  all five research tasks  /
           integrate-phage-pangenome-plan
```

The five research tasks run in parallel after the quality pass and own disjoint artifacts. The integration task waits on **all five**. Internal `.assign-*` and `.flip-*` lifecycle tasks may also appear in WG; inspect them but do not mistake them for scientific deliverables.

| Task ID | State observed at handoff | Intended artifact(s) | Dependency/role |
|---|---|---|---|
| `.quality-pass-phage-pangenome-planning` | Done | WG task metadata/descriptions and validation log; no repository report was promised. | Predecessor quality/safety review for all five research tasks. |
| `storage-genome-acquisition-inventory` | Active | `reports/storage_and_acquisition.md`; `artifacts/genome_input_inventory.tsv`; `artifacts/genome_acquisition_estimate.tsv` | After quality pass; owns accession inventory, timestamped storage evidence, acquisition estimates, paths, and resource guardrails. Feeds integration. |
| `prophage-table-distribution` | Active | `reports/prophage_distribution.md`; `artifacts/prophage_summary/per_genome.tsv`; `summary_metrics.tsv`; `category_counts.tsv`; `interval_qc.tsv`; plots; `reproduce.sh` and its source under `artifacts/prophage_summary/` | After quality pass; owns offline schema/tag/distribution/join/interval audit. Feeds integration. |
| `impg-syng-capability` | Active | `reports/impg_syng_assessment.md`; `artifacts/impg_probe/commands.sh`; `environment.txt`; `version.txt`; `help.txt`; optional bounded synthetic evidence under `artifacts/impg_probe/synthetic/` | After quality pass; owns installed IMPG terminology, formats, semantics, probe, scalability, and pilot gates. Feeds integration. |
| `mash-phylogeny-design` | Active | `reports/phylogeny_design.md`; `artifacts/phylogeny_probe/commands.sh`; `environment.txt`; `tool_versions.txt`; optional bounded benchmark under `artifacts/phylogeny_probe/benchmark/` | After quality pass; owns Mash-family host-tree design, pair/resource arithmetic, higher-fidelity validation, and gates. Feeds integration. |
| `pansn-bgzip-genome-layout` | Active | `reports/pansn_bgzip_naming.md`; `artifacts/pansn_identifier_crosswalk_template.tsv`; `artifacts/pansn_bgzip_probe/commands.sh`; `environment.txt`; `tool_versions.txt`; bounded synthetic compatibility evidence under `artifacts/pansn_bgzip_probe/synthetic/` | After quality pass; owns authoritative PanSN policy, lossless crosswalk, BGZF/index workflow, and BGZF × literal-`#` compatibility. Feeds integration. |
| `integrate-phage-pangenome-plan` | Open/waiting | `reports/phage_pangenome_project_plan.md`; `artifacts/project_manifest_template.tsv` with example rows clearly marked | Depends on all five research tasks. It must consume and cite their completed evidence, reconcile conflicts/unknowns, and define the approved downstream task graph without launching production. |

Task states can change after this document is written. Trust `wg show`, not the state label above. Do not consume an unfinished worker’s worktree or progress log as a published result. Once tasks complete and their artifacts are available on the working branch, use the reports in `reports/` and machine-readable evidence in `artifacts/`.

## 8.2 Safe orientation commands

From the WG-managed repository/worktree, run:

```bash
wg quickstart
wg service status
wg viz --no-tui
wg show <task-id>
```

Use `wg show` for every task listed above before acting. If the dispatcher service is running, do not manually claim or spawn active work. Do not redo an active task; send a WG message if coordination or clarification is needed.

## 8.3 Next-agent checklist

- [ ] Confirm the two root inputs are unchanged; do not edit them.
- [ ] Run the four orientation commands and inspect all seven task IDs above.
- [ ] Check which research tasks are done; consume only completed `reports/` and `artifacts/` outputs.
- [ ] Do not reproduce unfinished audits or turn provisional notes into facts.
- [ ] If all five reports are complete, work through `integrate-phage-pangenome-plan` rather than inventing an independent synthesis.
- [ ] Reconcile each report’s versions, assumptions, denominators, coordinate rules, resource estimates, and contradictions.
- [ ] Freeze a small stratified pilot manifest and explicit pass/fail/resource gates before any acquisition or compute.
- [ ] Keep host clades, phage clusters/networks, core prevalence, and ancestral inference separate in schemas and reports.
- [ ] Require lossless accession/contig/coordinate/PanSN/GFF/tree/IMPG/phage-cluster joins at every phase.
- [ ] Obtain approval before bulk download, production compression/indexing, O(n²) distance work, or IMPG/pangenome builds.

## 8.4 Open decisions

| Decision | Owner/evidence needed |
|---|---|
| Exact accession-file counts, identifier kinds, assembly resolution/version policy, acquisition sizes/times, reusable local data, durable/scratch paths, and hard byte/inode floor | `storage-genome-acquisition-inventory`, then integration/pilot approval |
| Prophage CSV schema, tagged-subset rule(s), denominator, deduplication/locus key, coordinate convention, joins, overlap/nesting and anomaly policy | `prophage-table-distribution`, then coordinate pilot |
| Installed IMPG version/terminology, accepted inputs, graph/path/query semantics, one versus multiple graphs/indexes, sharding/resume, resource controls, and scale feasibility | `impg-syng-capability`, naming compatibility probe, then IMPG pilot |
| Primary Mash-family route, sketch parameters, all-genome versus representative stages, matrix materialization, higher-fidelity host-core validation, rooting, support, and host-clade definition | `mash-phylogeny-design`, then tree-stability pilot |
| Normative PanSN version and canonical sample/haplotype/contig policy; assembly revisions; escaping/collisions; BGZF settings/checksums/indexing; direct consumer compatibility and fallbacks | `pansn-bgzip-genome-layout`, then integration |
| Pilot sample size, strata quotas, throughput/peak-RAM/disk limits, free-space stop values, and production concurrency | Integrated evidence and measured pilot; exact values are TBD |
| Prophage gene caller, family/homology and synteny methods, clustering resolutions, minimum supported cluster/cell size, core and soft-core thresholds, and missing-call model | Later approved analysis design plus sensitivity pilot |
| Ancestral-state model, treatment of HGT/mosaicism/detection error, topology uncertainty, external/reference-phage set, and reporting calibration | Later evolutionary-analysis design; prevalence alone cannot decide this |

The immediate safe next step is therefore **evidence integration and a bounded pilot design**, not production acquisition or computation.
