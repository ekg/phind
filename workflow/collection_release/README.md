# Collection release resolver

`resolver.py` freezes the supplied accession occurrences and exact assembly revisions. It calls only two NCBI Datasets v2 JSON metadata endpoints (`genome/check` and `genome/dataset_report`); it contains no sequence download path.

Key properties:

- verifies both immutable root-input SHA-256 values and the manifest template before work;
- rejects anything except exact versioned `GCF_#########.#` values and preserves raw lines as Base64;
- uses deterministic occurrence, assembly, release, batch, and Stage-B identities;
- archives request bodies, responses, headers/receipts, and response SHA-256 values externally;
- rate limits requests, retries transient failures, and reuses only checksum-valid cache units;
- fails closed on cardinality, exact-version identity, provenance, resource, RAM/swap, or publication errors;
- writes `COMPLETE` last and atomically renames a same-filesystem staging directory.

## Tests

```bash
python -m unittest -v workflow.collection_release.test_resolver
```

## Production invocation

All resource allocations are required integer arguments; omitted or zero allocations are `NO_GO`. See the argv captured in the external release's `provenance.json` for the executed command. Consumers must verify `COMPLETE`, `SHA256SUMS`, `release.json`, and tracked manifest digests before use. `resolution_status=EXACT_VERSION_VALID_METADATA_UNAVAILABLE` means `/genome/check` validated the requested exact revision but the detailed-report endpoint omitted it; it is not permission to substitute `latest`.
