# Pilot consumer compatibility certification

**Verdict: PASS / GO for the pinned N≤1,000 pilot consumers below.**  This is
an engineering compatibility result on the exact validated ten-assembly
release, not a biological-scale analysis and not authorization for consumers or
parameters absent from this matrix.

## Immutable inputs and release

The workflow consumed predecessor
`canonical-cohort-010-v1-e71484de9994fc28`; tracked `release.json` SHA-256 is
`4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23`.
Before invoking a consumer it verified the tracked manifest inventory, the
external `COMPLETE`, every external `SHA256SUMS` row, all 80
`checksums.tsv` rows, all canonical BGZF/FAI/GZI digests, decompressed content,
exact cohort order, 10 assembly revisions, 1,223 unique PanSN paths and
51,731,662 bases. The graph-wide sequence-bearing union remained exactly those
10 revisions, below the immutable cap of 1,000.

Final release:

```text
consumer-compatibility-v1-78d7e93f19fa3d87
/home/erikg/phind-data/ecoli26k/v1/releases/certify-pilot-consumer/consumer-compatibility-v1-78d7e93f19fa3d87
```

It has 19/19 unqualified consumer gates `PASS`, 36 checksum-inventory rows,
an append-only state/failure/command/resource ledger, `COMPLETE` written last,
and same-filesystem atomic directory promotion. Independent validation passed
13 groups. `manifests/consumer-compatibility-v1/` is the only downstream
selection authority; compact earlier attempts are isolated under external
`quarantine/` and cannot validate at their new paths.

Both immutable root inputs matched before and after:

| Input | SHA-256 |
|---|---|
| `26k_ecoli_accession.txt` | `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5` |
| `26k_prophage1.csv` | `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996` |

## Exact tool identities

The disposable host/gene environment is frozen by explicit URL/build/MD5 lock
`0b8db59d01eed5762db2bcb52e581e66b61750988c6db98cc36cc4721e53ccd4`
and an independent per-package SHA-256 inventory
`b473c14076ebf01c35e3a496a4fde403ce7bf90a874bb20ec4c17afb01c8c34a`.
The existing graph environment is frozen by explicit lock
`51c39b35707d7c722a8cc814ade76b7ab968335b021946ef5a6433e9b214a899`
and package-SHA inventory
`ac10fafe863c509c3c8ae0a27928455e70e5009a370b762315477c66b93a839f`.
The former was installed only below task scratch and deleted; no system package
changed. Exact version **and help** output, absolute executable, bytes and full
SHA-256 are in `artifacts/consumer_compatibility/tool_versions.json`.

| Tool | Exact identity | Executable SHA-256 |
|---|---|---|
| bgzip | HTSlib 1.19 | `e1ca105c4785d70fa8ee21d3a3278605f44d1670585e3c271a835dab29ba5dc9` |
| samtools | 1.19.2 / HTSlib 1.19 | `e795122c091f1795179a57904087dbacc3106c90e8b2f65ec86c6714586b6861` |
| IMPG | 0.4.1 installed build | `509296fb5c052be291a1841ea41f9bd4eb98e49b58b5f22cd69603729a94285f` |
| Mash | 2.3 | `83c85c063118c8c12659baa5d990aa8821f65cde40d1c82c8eaedd144dfec205` |
| RapidNJ | 2.3.2 package | `73bbf9615f3592d540084241634b4f09f1205e75107bc529efc4359453dd0208` |
| skani | 0.3.1 | `b1d20cb7170fe40a964526eadeb6fcdf61eefa748b88ed701ec3ff8dfbc07f5f` |
| QUAST | 5.3.0 | `3da08970cfc859473880823b1c6865396202f28ff9410509b4b014e03f26bc7c` |
| gffread | 0.12.7 | `b82dc8654802b2973d1b7df3538f61c301eefce67f4ffeb8e11b8293ad5e4c6d` |
| pggb | 0.6.0 | `42c710a28788fa910e2c9b4fe0e367464f9c4bd143ca856f7b99aa8daf204706` |
| wfmash (pggb) | 0.13.0-0-gd7b6960 | `6594211403359e47b262535b30ab9507c35fc83a1bfcca836bcf02a843e531f7` |
| seqwish | 0.7.10-0-g75e807c | `04cea039dc4aa861313a6360e1872f6b5f1ef56efd4ffbdf1ea25393cd56f22e` |
| smoothxg | 0.7.4-0-ge91d6b3 | `ddb4271b5d388a5f09050f92303993f5ddb4bafb551e64d21ec8399665142d6b` |
| odgi | 0.8.6-0-ge647844f | `1079640e20de96b38a568310e131f121c998cf71a79ce629c161210d07639980` |
| vg | 1.40.0 “Suardi” | `1e6b8cd6473395a37d4c78d2de7c1b7a6a125400f91b71f5736280c0beb15a20` |
| Prodigal | 2.6.3 | `b595f2855e7407dc4a2e3b6aa60704ab6478170630327d41802f9a64491a9e16` |
| MMseqs2 | 18.8cc5c | `d6d1bf1b50c095eb0cb66f685a7fb1d9b8b69497e1e88416fd48f28e313cc773` |
| HMMER hmmbuild / hmmsearch | 3.4 | `cb102977399e55adaf28d13d8cdf5caf0ce0980872f520c5f42c9ff19afdc411` / `afbcad8930a4b1252023fcb2f9d35e3a5a5c92ce768370ea6576a3644dff51ab` |
| MCL | 22-282 | `991813a1183ef3d4836e92cb439124c93dde06225dbe2a079120ccc756b7bad7` |

The installed IMPG binary has no recoverable embedded source commit; its exact
artifact identity is the binary digest, not a guessed checkout SHA.

## Required consumer contracts

The exact invocation, input form, view ID, output-name behavior and status for
all 19 logical gates are machine-readable in
`manifests/consumer-compatibility-v1/consumers.tsv` and each gate JSON. Summary:

| Boundary | Input/view contract | Observed result |
|---|---|---|
| bgzip + samtools faidx | canonical BGZF plus predecessor `.fai/.gzi`, direct | All 10 integrity checks; quoted literal-`#` lookup; name, bases and 1-based region-to-sequence round trip; corrupt stream/unknown path rejected. **PASS** |
| IMPG `syng` build | `combined-bgzf-all10-v1`: streamed 10-file ordered BGZF, identity names, 200 MB quota | Six inseparable files, all 1,223 names exact; missing input rejected; byte-identical rebuild. **PASS** |
| IMPG `query` | SYNG prefix + quoted exact 0-based half-open path range, direct | Source range recovered; unknown path rejected; raw BED is correctly treated as syncmer-padded rather than exact alignment. **PASS** |
| IMPG `map` | SYNG prefix + bounded BGZF query, direct | Query/target names exact, origin candidate present, missing query rejected, byte-identical rerun. PAF is a **syncmer-anchor projection, not a base alignment**. **PASS** |
| Mash | ten canonical BGZF files, whole-file mode only (never `-i`) | Exactly one sketch/assembly, 45 off-diagonal pairs, missing input rejected, byte-identical triangle. **PASS** |
| RapidNJ | `mash-to-rapidnj-full-v1`: strict lower-to-full PHYLIP adapter | RapidNJ 2.3.2 rejects Mash lower PHYLIP directly. The adapter validates row `i` has `i` distances, expands symmetry, inserts zero diagonal, changes no labels, then produces exactly 10 quoted Newick tips. Malformed matrix rejected. **PASS only with this staged adapter.** |
| skani | canonical BGZF, direct | Self ANI approximately 100, filename retained, missing input rejected. **PASS** |
| QUAST | canonical BGZF, direct | Metrics report produced; missing input rejected. Assembly identity must still join through manifest, not display label. **PASS** |
| gffread | `gff-semantic-alias-v1`: bounded plain FASTA + column-1 semantic alias view | gffread compares raw seqids: lexical `%23` does **not** match FASTA `#`. Strict one-layer decode plus recorded lexical↔semantic map succeeds; 1-based closed `101..300` equals canonical `[100,300)` and exact bases. Orphans rejected. **PASS only with adapter.** |
| pggb → seqwish → smoothxg | `graph-fixture-v1`: two-path BGZF derived from one authorized canonical path, 100 MB quota | Literal `#` retained into final GFA; missing inputs rejected. **PASS** |
| odgi + vg | pggb final GFA, direct | Listed path names and path-spelled FASTA bases exactly equal both inputs; missing inputs rejected. **PASS** |
| Prodigal | bounded BGZF graph fixture, direct | BGZF works and protein/gene semantics exactly equal plain input; literal `#` survives protein IDs and GFF seqids; missing input returns exact exit 5. **PASS** |
| MMseqs2 | bounded Prodigal proteins | Cluster TSV retains full `#` identifiers; missing input rejected. **PASS** |
| HMMER | bounded protein alignment/HMM + proteins | Profile build/search recovers exact `#` target; missing alignment returns exact exit 6. **PASS** |
| MCL | ABC network over exact protein identifiers | Cluster labels round-trip exactly; missing input rejected. **PASS** |

### Identifier and coordinate findings

- Literal PanSN `sample#haplotype#contig` survived every applicable tested
  build/query/map/graph/gene/clustering output.
- Unsafe raw `ctg#A/plasmid|β and space` became the byte-injective PanSN field
  `ctg%23A%2Fplasmid%7C%CE%B2%20and%20space`; at the GFF lexical layer,
  semantic `#` becomes `%23` and the existing `%23` bytes become `%2523`.
  Lowercase/noncanonical escapes and unknown semantic IDs fail closed.
- Samtools establishes exact base/coordinate retrieval. IMPG `query` uses
  indexed 0-based half-open coordinates but `--syng-raw` boundaries are padded.
  IMPG `map` is deliberately not used to claim base-exact coordinates.
- GFF3 coordinates remain 1-based closed until the declared one-time conversion;
  the test recovered the exact canonical `[100,300)` sequence.

## Determinism, restart, resources and cleanup

The ordered combined BGZF, six IMPG sidecars, IMPG query, IMPG map and Mash
triangle were byte-identical on exact rerun. A real SIGKILL left only an
unpublished staging directory; the wrapper detected/removed it, and a clean
restart promoted only after checksums and `COMPLETE`. Manifest/checksum
mutation, traversal, blank/overallocated resources, malformed Mash matrix,
unsafe IDs and interrupted promotion are permanent tests.

Live preflights recorded `findmnt`, ownership/write probe, free bytes/inodes,
explicit allocations and unfinished-write reserve before install/consumers and
at finish. Assigned RAM was 8 GiB; peak RSS was 906,067,968 B (10.55%), with no
OOM and zero swap growth. Scratch peaked at 2,719,568,941 B and 23,309 files:
0.00075 of the 4 TB allocation and 0.23309 of 100,000 inodes. Projected files
were conservatively reset to 25,000 (0.25). Durable upper-95 planning fraction
was 0.20; unfinished-write reserve factors were 8× durable and 7,994× scratch.
Scale-exponent/slope gates are not applicable to this non-scale-bearing
compatibility task.

Every sequence/index/tool-environment view was under the task scratch namespace
and is absent after success. The final external release is 328 KiB; Git-owned
outputs are compact logs/manifests/code only, with no FASTA, graph, index, hit
or cache payload.

## Explicit exclusions

- **CheckM2 is not selected** for this isolate-assembly pilot and has no frozen
  model-database digest; downstream may not assume it compatible.
- **FastANI is not selected**; skani 0.3.1 is the pinned ANI route.
- **Panaroo/Roary are not selected**; the concrete phage pilot route is
  Prodigal → MMseqs2/HMMER/MCL plus cluster-scoped pggb graph tooling.
- This report does not turn Mash/RapidNJ output into a phylogeny. It is only an
  **unrooted host genomic-similarity dendrogram** pending higher-fidelity
  validation.
- It does not conflate IMPG build, coordinate query and sequence map, and it
  does not claim that IMPG clusters genes or emits the final presence matrix.
