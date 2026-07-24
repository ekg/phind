# Prophage source-semantics gate

This directory implements the versioned, fail-closed policy for the immutable
`26k_prophage1.csv`. It performs no extraction, clustering, sequence download,
or sequence-bearing analysis. The only bounded diagnostic reads source GFF
members inside the already validated ten-assembly predecessor release; it reads
zero FASTA bases and does not copy the packages.

## Contract

- `semantics-policy-v1.schema.json` defines the policy shape.
- `artifacts/prophage_semantics/semantics_policy_v1.json` is the current policy.
- Every raw source row remains in the immutable CSV. The policy stores raw-field
  and row-identity requirements; it does not materialize a lossy normalized table.
- `all_records`, `transposable_flag_positive`, and `taxonomy_assigned` are three
  distinct scopes. None is aliased to the unresolved term `tagged`.
- Both plausible caller-family coordinate conventions remain explicit candidates:
  1-based closed `[b,e] -> [b-1,e)` and 0-based inclusive `[b,e] -> [b,e+1)`.
- `EXTRACTION_BLOCKED` is an expected scientific result, but a hard downstream
  gate. Consumers must invoke `validate --require-extraction-go` and reject it.

## Reproduce

Run tests first:

```bash
python -m unittest -v workflow.prophage_semantics.test_release
```

The production command requires nonblank resource allocations. The values used
for the v1 release were:

```bash
python workflow/prophage_semantics/release.py run \
  --repo . \
  --durable-root /home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source \
  --scratch-root /mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/resolve-prophage-source \
  --run-id semantics-v1-final-20260724 \
  --assigned-ram-bytes 8589934592 \
  --durable-allocation-bytes 1000000000 \
  --scratch-allocation-bytes 4000000000000 \
  --inode-allocation 100000 \
  --predicted-durable-peak-bytes 50000000 \
  --predicted-scratch-peak-bytes 50000000 \
  --predicted-files 100 \
  --unfinished-write-bytes 10000000
```

For the injected interruption test, add `--inject-stop-before-complete`; exit 75
is required, no `COMPLETE` may exist, and the exact command without that flag
must validate every existing static unit before resuming. Publication writes
`COMPLETE` last, fsyncs it, then atomically renames the whole staging directory.

Validate the release and demonstrate consumer refusal:

```bash
python workflow/prophage_semantics/release.py validate \
  /home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source/prophage-semantics-v1-f5619e221ff272ae

# Intentionally exits 2 while the verdict is not EXTRACTION_GO.
python workflow/prophage_semantics/release.py validate \
  /home/erikg/phind-data/ecoli26k/v1/releases/resolve-prophage-source/prophage-semantics-v1-f5619e221ff272ae \
  --require-extraction-go
```

An existing complete release is never rewritten. An identical rerun validates
its full checksum inventory and returns
`EXISTING_IMMUTABLE_RELEASE_VALIDATED`.
