# Canonical N=250 branch contract audit

## Verdict and scope

**`CANONICAL_GO_500` — canonical acquisition/canonicalization only.** This is
not an extraction GO, an integrated-analysis GO, or a claim that prophage
coordinates are interpretable. The audited source release is the already
published immutable `canonical-cohort-250-v1-a6184d7d6ee08bda`; this audit did
not alter that release or any predecessor.

Machine handoff:

- `artifacts/canonical_cohort_250_contract_audit/verdict.json`
- verdict SHA-256: `699d14b32c010771280b193b2373968dcae0c0c130a87a91270413eadd9c03e5`
- `artifacts/canonical_cohort_250_contract_audit/audit.json`
- audit SHA-256: `8c5ce43261dd99d6c4ba6b52bf0e5e4e6a9c86168f5c95b452f2ab192cd0d8d1`
- both are pinned by
  `artifacts/canonical_cohort_250_contract_audit/SHA256SUMS` and can be consumed
  fail-closed with `python3 artifacts/canonical_cohort_250_contract_audit/consume_verdict.py --repo-root .`.

## Contract inconsistency and ruling

The historical `prepare-canonical-cohort-250` task still says canonical
preparation requires an integrated N=100 automatic `GO_250`. That is the stale
clause identified by the evaluator. The subsequently authorized graph contract
separates canonical acquisition/canonicalization scaling from the failed
prophage-source/extraction branch. The published N=250 report already records
the integrated gate as `NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY`, rather than
fabricating a GO (`reports/canonical_cohort_250.md:16-24`).

This audit applies that corrected canonical-only contract and preserves history:

1. the stale WG prose is identified verbatim and remains unchanged;
2. gates needed to acquire, checksum, canonicalize, index, annotate, resume,
   resource, and hand off canonical objects are evaluated;
3. prophage extraction/source-coordinate and integrated biological gates are
   explicitly `NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_*`;
4. the failed source-semantics branch and its strict rejection remain in force.

The graph check is implemented independently in
`artifacts/canonical_cohort_250_contract_audit/audit.py:539-588`. The machine
verdict is constructed with a canonical-only scope and an explicit
`integrated_go_claimed=false` at
`artifacts/canonical_cohort_250_contract_audit/audit.py:1058-1088`.

## Independent method and immutable-byte result

The audit code does not import the N=250 runner, its validator, or the pinned
acquisition implementation. It uses Python standard-library parsers plus the
pinned local `bgzip`/`samtools` executables. Inventory coverage and COMPLETE
pins are independently checked by
`artifacts/canonical_cohort_250_contract_audit/audit.py:102-140`; exact nested
selection identity is re-derived at lines 179-190; archive/MD5/accession checks
start at lines 229-302; BGZF/FAI/GZI/faidx checks start at lines 357-433; and
GFF alias/coordinate recomputation starts at lines 472-529.

Start and finish snapshots were identical:

| Immutable view | Files | Bytes | SHA-256 |
|---|---:|---:|---|
| external N=250 release | 2,269 | 581,575,254 | `5a11927a00486d327d0a33f42ec6e4361b1716118a4e623e1da99401038cad79` |
| tracked N=250 manifest view | 10 | 6,892,308 | `37e595ef8d745f0b0c9074f8f9a88aab73b117b671e0fde35e97dff37ff6853d` |

Additional immutable pins independently matched:

- tracked and external `release.json`:
  `dcf2b887afa51e4e0e739ae2fef9b5a9d72fb8bc9a4d698a161a99673aaf504a`;
- external `SHA256SUMS`:
  `45fd42b76bf1c1ace3a2e882fe6a9a8f6af2457c0f5d4bc28011cb99b521c5b7`;
- external `COMPLETE`:
  `026ba58d2865915284df8e32abb05ecb1ace862650f5f784144e751dc4b9bce0`;
- exact 250-row cohort:
  `ba2cf2909ccf62a0c1944a76b522edc5600953511ec355479117b4a419acbc9f`;
- N=100 release JSON:
  `3b91b24e23323ef971a13f22825e512a233bb592ed641ea9b270a2f1fd683795`;
- recursive N=10 release JSON:
  `4cf1e5f7abb11d13dbae886543a343b0a57a389b46aa3df4ebc4fb14d280ff23`.

The two root inputs matched before and after:

- `26k_ecoli_accession.txt`:
  `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`;
- `26k_prophage1.csv`:
  `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`.

## Applicable canonical gates

| Gate | Independent result | Evidence summary |
|---|---|---|
| Exact selection and nesting | PASS | 250 unique exact accession revisions in frozen order; ordered N=100 prefix and N=10 prefix match on every frozen identity/stratum/occurrence field (rung-specific design weights are separately row-hashed). |
| Recursive object identity | PASS | Rows 1-10 resolve to immutable N=10, rows 11-100 to immutable N=100, rows 101-250 to SELF; 250 refs and all 500 source/canonical object inventories rehashed. |
| SHA inventories | PASS | 2,267 external release rows, 9 tracked release rows, 24 upstream artifact rows, 2,000 role-specific checksum rows, recursive predecessor inventories, exact coverage, no symlinks. |
| Archive/upstream identity | PASS | All 250 ZIPs passed structure/CRC, bounded expansion, exact accession catalog/path/report identity, upstream MD5 coverage, and local member/package SHA-256 replay. |
| BGZF/PanSN/index | PASS | All 250 BGZF streams passed `bgzip -t`; source-to-canonical rename-only length/order/sequence hashes, 41,050 globally unique reversible PanSN names, every FAI/GZI, and all 41,050 batched `samtools faidx` prefix regions passed. |
| Crosswalk/annotation | PASS | All object/cohort crosswalk rows reproduced. All 237 published alias views and their 1-based-closed bounds reproduced; all 13 empty quarantined views reproduced the exact source-GFF failure and exposed no partial alias view. |
| Restart/atomicity | PASS | Acquisition SIGKILL → unsafe partial restart and conversion SIGKILL → stage discard were ordered in the append-only ledger; requests were confined to rows 101-250; final event was `READY_TO_PROMOTE`; no partial/stage/symlink/plain FASTA survived. |
| Deterministic rerun | PASS | Published zero-network/zero-download/zero-recompression rerun pins matched; two independent audit executions produced byte-identical semantic audit JSON while immutable release/root snapshots stayed identical. |
| Resource and projection | PASS | All 399 preflights and 25 batch rows passed allocation, write, free-byte/inode, swap, reservation, measured, and upper-95 checks. See below. |
| Global cap | PASS | Live release scan union was exactly the frozen N=250 set, a subset of frozen N=1,000 and below cap 1,000. |
| Compact Git | PASS | No canonical task-owned Git file exceeds 10 MiB; no ZIP/BGZF/FAI/GZI/GFF/plain FASTA payload is tracked; this task adds only compact audit code/evidence and this report. |

The per-object deep loop is at
`artifacts/canonical_cohort_250_contract_audit/audit.py:735-856`; resource
checks begin at lines 869-954; and immutable finish equality is enforced at
lines 984-988.

### Resource floors, trend, and N=500 projection

All recorded durable samples remained at least 2 TB free and 1,000,000 free
inodes. All scratch samples remained at least 4 TB free and 5,000,000 free
inodes. Allocations were positive; 2× unfinished-write reservation, ≤70% RAM
and disk, ≤50% inode, zero process swap events, and zero sampled swap growth
all passed.

The corrected canonical branch makes canonical preparation trend evidence
applicable without pretending it is integrated biological evidence. Independent
N=100→N=250 recomputation gave:

- canonical preparation time exponent: **0.862526** (limit 1.3);
- wall/new object change: **−6.78%**;
- wall/new canonical base change: **−6.32%**;
- wall/new source byte change: **−6.15%**;
- stage bytes/new object change: **−0.16%**;
- stage files/new object change: **−0.41%**.

Every last-two-rung slope change is within 25%. Conservatively scaling the
measured N=250 incremental workload to 250 new N=500 objects and adding 25%
gives 1,211,035,584 stage bytes, 4,717 files, 497,612,800-byte RSS, and 863.5 s
wall. Those projections remain below 70% disk/RAM and 50% inode comparison
allocations. N=500 must still repeat its own live preflight; this audit grants
no waiver. The projection calculation is at
`artifacts/canonical_cohort_250_contract_audit/audit.py:930-947`.

## Extraction remains blocked; no gate relabeling

`resolve-prophage-source` is still terminal `failed` with
`EXTRACTION_BLOCKED`; immutable semantics release
`prophage-semantics-v1-f5619e221ff272ae` still says consumer action `REJECT`,
and its strict consumer still exits 2. WG reports each of
`run-integrated-100-genome`, `run-integrated-250-genome`,
`run-integrated-500-genome`, and `run-integrated-syng` transitively blocked by
that failed source branch.

The N=250 release field `source_coordinate_annotation_policy=PASS` means only
that retained source **GFF annotation alias views** obey their documented GFF3
coordinate bounds. It is not a prophage extraction/source-coordinate verdict.
The latter is recorded here as
`NOT_APPLICABLE_CANONICAL_PREPARATION_ONLY_EXTRACTION_BLOCKED`, never PASS.
Likewise integrated extraction GO and integrated biological scale trend are
NOT_APPLICABLE to this canonical preparation. The complete separation is
machine-readable in `audit.json` and `verdict.json`, and is rechecked by the
consumer (`artifacts/canonical_cohort_250_contract_audit/consume_verdict.py:35-80`).

## Validation and side effects

- Independent audit: PASS twice, deterministic semantic bytes.
- Audit failure-mode tests: 6/6 PASS (digest mismatch, inventory coverage,
  row hashing, nested identity, reversible PanSN encoding, and scope-safe
  verdict construction).
- Existing project regressions: 63/63 PASS.
- Audit network requests: 0; sequence downloads: 0; canonical objects written
  or recomputed: 0; release writes: 0.
- The audit read local sequence/archive bytes to recompute validation hashes and
  index/annotation invariants; it did not reacquire, recanonicalize, copy, or
  publish sequence.
