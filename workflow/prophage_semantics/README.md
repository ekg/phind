# Prophage source-semantics gate

This directory implements the versioned, fail-closed policy for the immutable
`26k_prophage1.csv`. It performs no production extraction, clustering, sequence
download, or sequence-bearing analysis. The only bounded diagnostics read
source-GFF members inside the already-validated ten-assembly predecessor
release and read the predecessor's already-extracted native prophage
TSV/saved-FASTA outputs; it reads zero host-genome FASTA bases and copies no
packages.

## Contract

- `semantics-policy-v1.schema.json` defines the policy shape.
- `artifacts/prophage_semantics/semantics_policy_v1.json` is the **historical
  v1 policy**: the source-alone investigation, kept as an immutable BLOCKED
  record (`prophage-semantics-v1-f5619e221ff272ae`).
- `artifacts/prophage_semantics/semantics_policy_v2.json` is the **current v2
  policy**: it authorizes `EXTRACTION_GO` on the basis of decisive, independently
  re-verified pinned-caller attribution.
- Every raw source row remains in the immutable CSV. The policy stores raw-field
  and row-identity requirements; it does not materialize a lossy normalized table.
- `all_records`, `transposable_flag_positive`, and `taxonomy_assigned` are three
  distinct scopes. None is aliased to the unresolved term `tagged`; `tagged`
  remains an explicit, non-extraction-critical unresolved user/source term.
- v2 selects coordinate candidate `C1_RAW_1_BASED_CLOSED` (1-based inclusive;
  canonical transform `[begin-1, end)` 0-based half-open) and rejects
  `C2_RAW_0_BASED_INCLUSIVE` (the Phigaro v2.4.0 convention).
- The extraction verdict is **derived from the pinned-caller consumption gate**,
  never hardcoded. A missing/non-PASS gate, a non-DECISIVE historical
  attribution, or an independent re-verification that is not sound leaves the
  release BLOCKED (the v1 BLOCKED record then remains active).

## Verdict derivation (fail-closed)

1. `pinned_caller_gate.py` consumes the immutable predecessor external release
   `phigaro-version-comparison-v1-e7cfa43b9231aee5`, verifying COMPLETE, the
   84-file SHA-256 inventory, exact N=10 identity/order, version-pinned
   tools/database/config, the official fixture gate, and the two **separate**
   machine verdicts (`historical_csv_attribution`, `modern_v2_4_pilot`).
2. It then **re-derives the attribution evidence itself**
   (`independent_rerun_verification.py`) — 56/56 v2.3-vs-CSV exact reproduction,
   the v2.3−v2.4 +1 boundary signature, and inclusive saved-FASTA lengths —
   without trusting the predecessor's comparison code or its `DECISIVE` string.
3. Historical extraction becomes `EXTRACTION_GO` only when
   `historical_csv_attribution == DECISIVE` **and** the independent
   re-derivation is sound. A modern v2.4 `GO` alone never authorizes historical
   extraction.

## Reproduce

```bash
python -m pytest -v workflow/prophage_semantics/

# 1. independently verify + derive the pinned-caller verdict
python workflow/prophage_semantics/pinned_caller_gate.py \
  --repo . \
  --output artifacts/prophage_semantics/pinned_caller_input_gate.json

# 2. build/validate the versioned semantics release (verdict derived from gate)
python workflow/prophage_semantics/release.py run \
  --repo . \
  --durable-root /home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source \
  --scratch-root /mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/resolve-prophage-source \
  --run-id semantics-v2-final-20260725 \
  --assigned-ram-bytes 8589934592 \
  --durable-allocation-bytes 1000000000 \
  --scratch-allocation-bytes 4000000000000 \
  --inode-allocation 100000 \
  --predicted-durable-peak-bytes 50000000 \
  --predicted-scratch-peak-bytes 50000000 \
  --predicted-files 100 \
  --unfinished-write-bytes 10000000
```

A missing/absent predecessor release, missing `COMPLETE`, inventory mismatch,
changed N=10 identity/order, incomplete tool/database/config pins, a failed
fixture/engineering gate, or an independent re-derivation that is not sound
exits 2 and yields `NO_GO` / `EXTRACTION_BLOCKED`. The current published release
is `prophage-semantics-v2-7dc695b85e5fd229` (`EXTRACTION_GO`).

For the injected-interruption test, add `--inject-stop-before-complete` against
a throwaway durable root; exit 75 is required, no `COMPLETE` may exist, and the
same command without that flag must validate every existing static unit before
resuming. Publication writes `COMPLETE` last, fsyncs it, then atomically renames
the whole staging directory.

Validate the release and demonstrate the consumer gate:

```bash
python workflow/prophage_semantics/release.py validate \
  /home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source/prophage-semantics-v2-7dc695b85e5fd229 \
  --require-extraction-go   # exits 0 with consumer_action ALLOW
```

An existing complete release is never rewritten. An identical rerun validates
its full checksum inventory and returns
`EXISTING_IMMUTABLE_RELEASE_VALIDATED`.
