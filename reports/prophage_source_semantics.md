# Prophage source semantics and extraction gate

## Release verdict

**`EXTRACTION_GO`** — consumer action **`ALLOW`** — policy **v2**.

The historical `26k_prophage1.csv` is **decisively attributable** to **Phigaro
2.3.0** native TSV output, with a reversible `genome`/`prophage_id` and scaffold
aggregation layered on top. All **56** frozen-cohort rows are reproduced
**exactly** (begin/end/transposable/taxonomy/scaffold) by a version-pinned,
independently re-verified Phigaro v2.3.0 rerun; the coordinate convention is
**1-based inclusive boundary-gene coordinates**. Extraction is authorized under
the selected candidate `C1_RAW_1_BASED_CLOSED`.

Current published semantics release:

- ID: `prophage-semantics-v2-7dc695b85e5fd229`
- durable path: `/home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source/prophage-semantics-v2-7dc695b85e5fd229/`
- tracked reference: `artifacts/prophage_semantics/release_reference.json`
- policy: `artifacts/prophage_semantics/semantics_policy_v2.json`
- schema: `workflow/prophage_semantics/semantics-policy-v1.schema.json`

Historical v1 record (kept immutable): `prophage-semantics-v1-f5619e221ff272ae`
with `EXTRACTION_BLOCKED` / policy v1 — the source-alone investigation before the
pinned-caller dependency resolved.

An extraction-dependent consumer must require exactly `EXTRACTION_GO`;
missing, `CONDITIONAL`, unknown, or `EXTRACTION_BLOCKED` is rejected. The
demonstrated strict consumer check (`validate --require-extraction-go`) exits 0.

### Mandatory pinned-caller dependency: resolved

The independent `rerun-pinned-phigaro` task published the immutable external
release `phigaro-version-comparison-v1-e7cfa43b9231aee5` under
`rerun-phigaro-version-comparison/`, with `COMPLETE` and an 84-file `SHA256SUMS`
inventory that round-trips. This task independently verified its release ID,
`COMPLETE`, full SHA-256 inventory, exact N=10 input identity/order, the
version-pinned tools/database/config, the official fixture gate, the
engineering gates, and the **two separate** machine verdicts. It then
**re-derived the attribution evidence itself** rather than trusting the
predecessor's comparison code or its `DECISIVE` string.

The machine gate in `artifacts/prophage_semantics/pinned_caller_input_gate.json`
is therefore **`PASS`**: `historical_csv_attribution=DECISIVE`,
`modern_v2_4_pilot=GO` (strictly separate), `decisive_evidence_independently_sound=true`,
and `historical_csv_extraction=EXTRACTION_GO` / consumer action `ALLOW`.

`workflow/prophage_semantics/pinned_caller_gate.py` implements the independent
verifier and decision rule; its tests prove that a missing predecessor, a
checksum-mismatched inventory, a non-decisive attribution, or an unsound
independent re-derivation fails closed, and — independently — that
`modern_v2_4_pilot=GO` plus `historical_csv_attribution=NON_DECISIVE` still
yields `EXTRACTION_BLOCKED`. No modern result can stand in for historical
attribution.

## Automatic immutable-input gate

All gates were evaluated automatically; no approval, override, substituted
object, new assembly download, or reduced cohort was used.

| Input/evidence | Required identity | Result |
|---|---|---|
| `26k_ecoli_accession.txt` | SHA-256 `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5` | PASS at start and finish |
| `26k_prophage1.csv` | SHA-256 `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996` | PASS at start and finish |
| integrated plan | SHA-256 `fb58d25a6f4971137ab0dcb82dae09eac5d177e37ba34f857845ce2d7e0a6da8` | PASS |
| source audit | SHA-256 `feb9b687fb0722a4f073f7105088f19a2dabebcbc0378dbe5c4799b0b7f29fdc` | PASS |
| canonical N=10 manifest | SHA-256 `4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23` | PASS |
| predecessor release | `canonical-cohort-010-v1-e71484de9994fc28`, exact order of 10 revisions | PASS |
| pinned-caller release | `phigaro-version-comparison-v1-e7cfa43b9231aee5`, `COMPLETE` + 84-row inventory | PASS |
| global exact-assembly union | exactly 10 revisions, declared cap 1,000 | PASS |

## Evidence, in priority order

### 1. Source provenance (the CSV itself)

The CSV header is exactly `end,genome,scaffold,begin,transposable,taxonomy,prophage_id`;
the file carries no caller/version/argv/config/output-format/strand/topology
declaration. **However**, the missing header provenance is superseded by
reproduction-based attribution below, which is stronger than a header claim. A
bounded metadata retry (file Git history, creation-commit tree, exact
filename/header/field web queries) found no sidecar or pipeline; commit
`7381847` added only the two data files. Recorded in
`artifacts/prophage_semantics/targeted_provenance_retry.json` (read no bases,
no downloads).

### 2. Pinned-caller reproduction — DECISIVE

Pinned Phigaro **v2.3.0** (commit `aea9469`) and **v2.4.0** (commit `1ff5f85`)
were rerun on the exact frozen N=10 cohort under pinned Prodigal 2.6.3 /
HMMER 3.3.2 / pVOG database / config (byte-identical tool binaries across both
envs); the official fixture gate passes for both versions.

Independently re-derived (not trusting the predecessor's comparison code):

| Check | Result |
|---|---|
| v2.3.0 vs CSV exact rows (all fields) | **56 / 56** (begin_delta=0, end_delta=0) |
| v2.4.0 vs CSV | 56 / 56 with begin_delta=-1, end_delta=-1 (off-by-one) |
| v2.3 − v2.4 boundary signature | **+1 on both ends for all 56** → uniquely identifies v2.3.0 |
| saved-FASTA length = end−begin+1 | **56 / 56** (end inclusive; 0 half-open) |
| CSV lacks the v2.4.0 `id` column | yes → matches v2.3.0, not v2.4.0 |

Post-processing is reversible: CSV `genome` = PanSN scaffold first component;
CSV `scaffold` = last component; `prophage_id` = `<genome>_prophage_<N>`;
`transposable` `1.0`/`0.0` = Phigaro `True`/`False`.

### 3. Upstream producer source (code-backed coordinate convention)

Fetched directly from the Phigaro repository at the pinned commits:

- **v2.3.0** `run_phigaro.py`: `begin = genes[phage.begin].begin; end = genes[phage.end].end; writer.writerow((scaffold.name, begin, end, transposable, taxonomy))` — raw boundary-gene coordinates, **1-based inclusive**, no offset.
- **v2.4.0** `run_phigaro.py`: `begin = genes[phage.begin].begin - 1; end = genes[phage.end].end - 1`, plus an added `id` column — **0-based**.

This code reference confirms both the convention (1-based inclusive for the
historical CSV) and the version signature (the off-by-one). Permalinks and
full-file SHA-256 values are pinned in
`artifacts/prophage_semantics/evidence_inventory.json`.

### 4. Bounded source-GFF diagnostic (NON_DECISIVE, retained for transparency)

The earlier N=10 source-GFF boundary diagnostic (read inside pinned package
ZIPs, zero host FASTA bases) showed raw boundaries align with NCBI annotation
for 45/56 rows while raw+1 aligns for 0/56. This is deliberately
**NON_DECISIVE** provenance — NCBI GFF is not the caller's Prodigal output and
is not a known-base oracle. It is retained as corroboration only; the decisive
evidence is the pinned-caller reproduction above.

## Versioned semantic decisions (v2)

Every dimension has a status, confidence, evidence list, and extraction-critical
flag in the machine policy.

| Dimension | Status / confidence | Decision |
|---|---|---|
| producer/caller/version | `RESOLVED` / `HIGH` | **Phigaro 2.3.0**; extraction-critical, now resolved |
| user's term `tagged` | `UNRESOLVED` / `NONE` | explicit unresolved user/source term; **not** silently relabeled; not extraction-critical |
| `transposable` | `RESOLVED` / `HIGH` | Phigaro `if_transposable()` Integration-group flag; raw `1.0`/`0.0` preserved |
| taxonomy labels | `RESOLVED` / `HIGH` | Phigaro `define_taxonomy()` modal pVOG codes; exact strings preserved |
| coordinate base/end | `RESOLVED` / `HIGH` | **1-based inclusive** boundary-gene coordinates; extraction-critical, now resolved |
| strand/orientation | `RESOLVED_ABSENT` / `HIGH` | absent from source by design (no TSV strand column); **not required** for lossless interval extraction; never inferred |
| topology/circularity | `RESOLVED_NONWRAPPING` / `HIGH` | no topology marker; all 132,404 rows are begin≤end (non-wrapping); extraction is a forward slice, topology-independent |
| contig-edge behavior | `RESOLVED_DECLARED_INTERVAL` / `HIGH` | extract declared `[begin,end]` verbatim; no clipping/extension; edge callability is interpretation, not extraction |
| completeness | `ABSENT_UNKNOWN` / `HIGH` | no completeness field; not required to extract the declared interval |
| duplicate-locus rules | `RESOLVED_NO_DUPLICATES` / `HIGH` | zero duplicate `(genome,scaffold,begin,end)` loci in the snapshot; no dedup; every row preserved |

### The three scopes remain distinct

| Scope ID | Reversible rule | Rows |
|---|---|---:|
| `all_records` | every immutable source row | 132,404 |
| `transposable_flag_positive` | `Decimal(raw transposable) == 1`; observed member text `1.0` | 7,695 |
| `taxonomy_assigned` | trimmed/case-folded taxonomy is neither empty nor `Unknown` | 115,442 |

All 132,404 rows, 132,405 physical lines, 26,077 genome keys, and every raw field
were accounted. There are zero duplicate exact locus groups, extra duplicate
locus rows, exact-record duplicate groups, and duplicate `prophage_id` groups.
No normalized row table was materialized; the immutable CSV remains the lossless
row store. Extraction is permitted for all three scopes alike; none is silently
relabeled.

## Coordinate policy (v2)

- **SELECTED — C1, raw 1-based inclusive**: `[b,e] -> extract contig[b-1:e]` (0-based half-open); length = e−b+1.
- **REJECTED — C2, raw 0-based inclusive**: this is the Phigaro **v2.4.0** convention (`begin-1`/`end-1`); the historical CSV matches v2.3.0 exactly, not v2.4.0.

The historical `EXTRACTION_BLOCKED` was changed to `EXTRACTION_GO` only because
`historical_csv_attribution=DECISIVE` **and** the evidence uniquely establishes
original caller/version/post-processing plus exact coordinate semantics — all
independently re-verified. `NON_DECISIVE`, mere similarity, or a modern-only `GO`
would have preserved `EXTRACTION_BLOCKED`.

## Release engineering and resource evidence

The external release contains its exact input manifest, policy, evidence,
profile/diagnostic, tools and argv provenance, append-only state/failure/resource
ledgers, resource and restart summaries, `release.json`, and complete
`SHA256SUMS`. An injected exit 75 occurred after static units but before
publication; `COMPLETE` was absent. Restart independently checked each completed
unit's SHA-256, appended resume events, and published no mixed files. `COMPLETE`
was created and fsynced last, followed by same-filesystem atomic rename.

Resource allocations were nonblank: 8 GiB RAM, 1 GB durable, 4 TB scratch, and
100,000 inodes; predicted peaks were 50 MB/50 MB and 100 files with 10 MB
unfinished-write reservation. Live durable/scratch floors, write probes,
ownership, `findmnt`, bytes/inodes, ≤70% RAM/disk, ≤50% inode, ≥2× unfinished
writes, no swap growth, and promotion preflight all passed. This N=10 metadata
diagnostic is non-scale-bearing, so exponent/slope gates are not applicable.

Nineteen task tests cover the pure decision rule (DECISIVE/NON_DECISIVE/sound
× modern GO/NO_GO), missing/checksum-mismatched predecessor fail-closed,
root/cardinality accounting, v1-BLOCKED and v2-GO policy validation,
manifest/checksum mismatch, blank/overcommitted resource refusal, interrupted
unit validation and mixed-resume refusal, absent-`COMPLETE` rejection, gate-
derived fail-closed verdict, exact predecessor inventory/cap/order, the live
pinned-caller PASS/GO integration check, and v2-release `--require-extraction-go`
validation. A deterministic rerun left every release byte unchanged. No
production extraction, clustering, pangenome processing, or new *E. coli*
assembly download occurred; the exact global distinct-assembly union remains 10
of at most 1,000.
