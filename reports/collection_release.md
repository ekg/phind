# Collection v1 accession/assembly freeze

**Verdict: PASS (metadata-only).** Release `collection-v1-f7494b4b89d1382b` is immutable and complete. No genome, annotation, index, or other sequence-bearing payload was requested or downloaded.

## Immutable contract

| Item | SHA-256 / value |
|---|---|
| `26k_ecoli_accession.txt` (start and finish) | `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5` |
| `26k_prophage1.csv` (start and finish) | `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996` |
| `artifacts/project_manifest_template.tsv` | `f99d90ed2a5aedc9595081911db2a4f9331923b90018aeb018a70c6007557afd` |
| external `SHA256SUMS` | `bd27dcca6c1e33eda8ccaafd30224aae5a2b5770ee6a87246474d6d55844f0e9` |
| external completeness marker | same digest, written last as `COMPLETE` |

The sole predecessor, `quality-pass-for`, was complete before execution and its graph/cardinality/resource-contract audits were logged as PASS; it declared no collection release or manifest artifact to consume. The required integrated template above was present and verified. There was no human wait, override, substituted input, or conditional verdict.

The durable release is:

`/home/erikg/phind-data/ecoli26k/v1/releases/freeze-collection-release/collection-v1-f7494b4b89d1382b/`

Consumers must reject it unless `COMPLETE` exists, hashes to the `SHA256SUMS` digest above, every inventory entry verifies, and tracked `release.json` has `verdict=PASS`. The external inventory check passed for all 174 files. The 146,349,493-byte release includes the lossless plain manifests and 109,227,103 bytes across 162 raw request/response/receipt files. Git contains only deterministic compressed views and compact evidence.

## Reconciliation and exact-version policy

All 26,078 physical LF-terminated input occurrences are represented. Lines 1–26,077 are unique exact versioned `GCF_#########.#` values. Line 26,078 is preserved losslessly and explicitly rejected as `REJECTED_MALFORMED_TERMINAL_TOKEN` for token `genome`. Occurrence IDs hash a domain tag, input digest, physical line number, and Base64 raw bytes; assembly IDs hash the exact accession.version. No versionless or `latest` query is used.

Fresh NCBI Datasets v2 metadata was retrieved in deterministic batches from only:

- `POST /datasets/v2/genome/check`
- `POST /datasets/v2/genome/dataset_report`

The exact-version check recognized all 26,077 accessions. Detailed reports were returned for 26,005; 72 exact versions were valid at `/check` but omitted by `dataset_report`, so they remain `EXACT_VERSION_VALID_METADATA_UNAVAILABLE` rather than being replaced. This includes Stage-B `GCF_000167895.3`, reproducing the earlier documented omission.

Among detailed reports, 25,291 have `assembly_status=current` and 714 have the authoritative terminal status `suppressed`. The manifest retains every supplied identity and its terminal evidence; downstream selection/acquisition must not treat a suppressed row as eligible. Thus the input-candidate identity accounting remains 26,077, while the non-suppressed candidate pool is 25,363 (25,291 current plus 72 exact-check-valid rows lacking detailed terminal status). Any downstream policy requiring affirmative `current` status has only 25,291 eligible rows. There are 26,005 paired accessions, BioSamples, and annotations; 25,011 strain values and 1,303 isolate values. Raw strings are Base64, not identifiers. No reported exact accession differed from its request, and no silent supersession occurred.

## Published manifests

External immutable plain manifests:

| Manifest | Rows | Bytes | SHA-256 |
|---|---:|---:|---|
| `manifests/occurrences.tsv` | 26,078 | 12,455,807 | `468bb37b2a1f2eb13d84ec6068b5e4ffad32697c6583591c283defdc238c3bcd` |
| `manifests/assemblies.tsv` | 26,077 | 24,570,044 | `4c2ce6558894996e68dd2f9daa6bcd12370fc300205d9fc68511ee4a47db8d24` |
| `manifests/stage_b_10.tsv` | 10 | 2,246 | `0d179cbafce2ba1fa14d1929a4acd6621810a335f25bcd7ec67dd2083eb101f6` |

Tracked deterministic views in `manifests/collection-v1/` are checksummed by its `SHA256SUMS`: `occurrences.tsv.gz` `7ef0563c…bbcc32`, `assemblies.tsv.gz` `5d72d583…a5041d`, identical Stage-B TSV `0d179cba…101f6`, and identical release JSON `59c6907e…6642bf`. Decompression exactly reproduces the external plain manifest bytes.

The Stage-B order is immutable: `GCF_000005845.2`, `GCF_000812325.1`, `GCF_002302315.1`, `GCF_004664255.1`, `GCF_015644385.1`, `GCF_020829045.1`, `GCF_921380995.1`, `GCF_000167895.3`, `GCF_001881595.4`, `GCF_000498835.2`. Exact IDs, input positions, selection rules, statuses, and per-row SHA-256 values are in `stage_b_10.tsv`.

## Provenance, resource, and restart evidence

The external release contains exact input manifest, raw API cache, response-level receipts and digests, `provenance.json` (tool, Python/platform, argv, bounded environment), append-only `state.jsonl`, empty append-only `failures.jsonl`, start/end resource records, output inventory, and final marker. Requests were limited to batches of 1,000, paced at no more than approximately 3/s, retried with bounded exponential backoff, and resumed only when request/response/receipt digests matched.

A deliberate kill after 55 seconds left no final directory and no `COMPLETE`. Restart with the same run ID recorded `RESUME_PREFLIGHT_PASS`, reused independently checksum-valid batches, completed the remaining batches, rechecked root inputs, and atomically renamed the staging directory. The state ledger has 49 starts, 48 validations, and one visibly interrupted start, followed by `ALL_METADATA_VALIDATED`, finish immutability PASS, and `READY_TO_PROMOTE`; no mixed output was published.

Resource evidence was captured before every batch. Durable root `/dev/nvme0n1p2` (ext4) started with 2,638,392,553,472 available bytes and 447,716,494 free inodes; scratch `/dev/nvme3n1` (XFS) had 5,506,473,230,336 bytes and 1,492,409,957 inodes. Write/ownership/mount gates passed. Explicit allocations were 4 GiB RAM, 10 GiB durable, 4 TB scratch, and 1,000,000 inodes; predicted peak disk/files were 250 MB/200 with 100 MB unfinished-write reservation. Peak RSS was 493,629,440 bytes (11.5% of assigned RAM, below 70%), swap growth was zero, and the external timed invocation reported zero swaps. Disk, inode, 2× unfinished-write, durable, and scratch floors all passed.

BGZF/index/name round-trip, coordinate policy, consumer compatibility, and scale-slope gates are not applicable to this metadata-only phase; they are not labeled PASS. Source/checksum, exact identity, row accounting, provenance, resource, deterministic resume, atomic promotion, and global-cap gates are PASS. Distinct sequence-bearing assemblies acquired or analyzed here: **0**, a subset of the frozen cohort and below the global cap of 1,000.

## Reproduction and fail-closed behavior

Run `python -m unittest -v workflow.collection_release.test_resolver`. Tests cover immutable-input mismatch, blank/insufficient resource refusal, unrequested/version-drift API response rejection, checksum-complete immutable resume, and interrupted staging cleanup/restart. `workflow/collection_release/README.md` documents invocation. The resolver has no sequence URL or payload-download implementation. Any source, cardinality, exact-version, response provenance, raw-cache checksum, RAM/swap, filesystem, row-count, inventory, or atomic-completeness failure exits `NO_GO` and cannot produce a consumer-visible release.
