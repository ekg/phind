# Independent Phigaro v2.4.0 N=10 pilot certification

## Machine verdict

**`REJECTED`** — not qualified for scaling. No later scaling decision may use this
attempt. The machine-readable decision is
`artifacts/phigaro_v2_4_certification/certification.json`; the one-line handoff
is `artifacts/phigaro_v2_4_certification/verdict.txt`.

This is a fail-closed evidence decision, not a claim that Phigaro itself failed.
There is no published pilot to certify: the predecessor supplied no tracked
release reference, producer validator, commit, or registered artifact, and the
specified external release parent does not exist. Consequently there are no
release bytes from which to confirm a Phigaro process, native outputs, or ten
outcomes.

## Independent method

The audit did not trust the predecessor `done` state. It:

1. attempted the producer's required strict-validator path exactly and recorded
   exit 2 (`strict_validator_run.txt`); the path is absent;
2. searched tracked Git paths for the promised reference/code/report and found
   none;
3. checked the exact external release parent and a bounded, depth-two set of
   Phigaro-named release directories; both searches found none;
4. ran the independent read-only fail-closed certifier twice; both executions
   exited 2 and produced byte-identical `REJECTED` JSON;
5. independently parsed the frozen N=10 cohort and hashed the two root inputs
   before and after the audit, comparing both to `HEAD` bytes.

The frozen-order constant and parser begin at
`workflow/phigaro_v2_4_certification/certify.py:14-25` and
`:148`; reference/producer/release discovery is enforced at lines 156-224;
release marker and exact-inventory checks are at lines 31-45 and 120-141; deep
semantic gates remain non-PASS without a digest-valid release at lines 266-278;
and `CERTIFIED_GO` is impossible unless every deep gate is `PASS` at lines
292-297. Failure-mode tests are at
`workflow/phigaro_v2_4_certification/test_certify.py:40-79`.

## Gate results

| Required gate | Result | Independent evidence |
|---|---|---|
| Tracked release reference | **FAIL** | `artifacts/phigaro_v2_4_pilot/release_reference.json` is absent and untracked. |
| Absolute external path | **FAIL** | Required parent `/home/erikg/phind-data/ecoli26k/v1/releases/complete-phigaro-v2-4-n10-pilot` does not exist; bounded search found no Phigaro-named candidate. |
| Producer strict validator | **FAIL** | `workflow/phigaro_v2_4_pilot/validate_release.py` is absent; exact attempted command exited 2. |
| `COMPLETE`, `SHA256SUMS`, complete inventories/digests | **NOT EVALUATED / blocking** | No referenced release directory or bytes exist. Missing evidence is not PASS. |
| Exact ordered release N=10 manifest | **NOT EVALUATED / blocking** | No release manifest exists. The independent *basis* passed: 10 rows, order 1..10, SHA-256 `83f8925eb261ea915192e3923fe8d858f09abfa6791a27dbda20530da44ec146`. |
| Phigaro/dependency/database/config pins | **NOT EVALUATED / blocking** | No tool, database, config, or environment ledger exists. |
| Real argv/process logs for all ten genomes | **NOT EVALUATED / blocking** | No argv ledger, stdout/stderr, exit status, PID/process evidence, or ten run records exist. A real Phigaro v2.4.0 process cannot be confirmed for any genome. |
| Ten outcomes/call counts versus native outputs | **NOT EVALUATED / blocking** | No outcomes ledger or TSV/BED/GFF/saved-FASTA output exists, so zero of ten outcomes or call counts can be reconciled. |
| Interval/base round trips | **NOT EVALUATED / blocking** | No native intervals or saved sequences exist. Historical CSV semantics were not inspected or reinterpreted. |
| Resource evidence | **NOT EVALUATED / blocking** | No producer elapsed/RSS/disk/inode evidence exists. The audit's own bounded `/usr/bin/time -v` record is not substituted. |
| Deterministic/restart evidence | **NOT EVALUATED / blocking** | No producer evidence exists. The independent rejection rerun was byte-identical, but that does not certify a missing pilot. |
| Root-input immutability | **PARTIAL PASS, overall blocking** | Current start/finish/`HEAD` hashes match (below), but the missing release has no producer start/finish root ledger. |

The exact ordered basis was:

1. `GCF_000005845.2`
2. `GCF_000812325.1`
3. `GCF_002302315.1`
4. `GCF_004664255.1`
5. `GCF_015644385.1`
6. `GCF_020829045.1`
7. `GCF_921380995.1`
8. `GCF_000167895.3`
9. `GCF_001881595.4`
10. `GCF_000498835.2`

Current root-input bytes were unchanged throughout the audit and matched
`HEAD`:

- `26k_ecoli_accession.txt`:
  `1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5`
- `26k_prophage1.csv`:
  `6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996`

## Predecessor and side-effect evidence

`predecessor_inventory.txt` records zero commits ahead on the predecessor
branch, only the WG worktree-management symlink/cleanup marker as untracked,
source evaluation 0.12, FLIP 0.64, and `No artifacts available from
dependencies yet.` This corroborates but does not replace the direct release
checks.

The audit wrote nothing below the external release tree, made no network
request, downloaded no host assembly, ran no Phigaro process, and dispatched or
executed no larger cohort. It treated `26k_prophage1.csv` only as opaque bytes
for SHA-256 immutability; no historical CSV semantics were interpreted.

## Validation

- Independent certifier: expected fail-closed exit 2, `REJECTED`.
- Deterministic certifier rerun: PASS, byte-identical JSON.
- Failure-mode unit tests: 3/3 PASS (missing reference, missing release, unsafe
  inventory entries).
- Full project Python test discovery: 58/58 PASS, recorded in `project_tests.txt`.
- Audit inventory: `artifacts/phigaro_v2_4_certification/SHA256SUMS`.

The only machine decision is **`REJECTED`**. It is not conditional GO and does
not support scaling.
