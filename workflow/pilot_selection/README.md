# Pilot cohort selection

`selection.py` freezes a metadata-only, nested cohort of exact assembly revisions. It never opens `26k_prophage1.csv` during selection, has no sequence-download endpoint, and performs no biological analysis.

## Frozen design

1. Verify the immutable collection release, its complete 26,077-row exact-version frame, the byte-identical Stage-B ten, and the validated Stage-B acquisition release.
2. Exclude only the predecessor's authoritative terminal `assembly_status=suppressed` rows. The 25,363 eligible rows are split into the exact ten Stage-B certainty units and 25,353 main units.
3. Keep Stage-B order 1–10. Rank every remaining main unit by `SHA256(seed NUL assembly_id NUL exact_accession.version)`, with `assembly_id` as a collision tie-breaker.
4. Freeze the first 1,000 identities; rungs 10/100/250/500/1,000 are exact prefixes.
5. Attach engineering-control labels only after order is frozen. Synthetic fixtures are excluded from host/prevalence inference; observed controls appear once with the manifest's design weight.

The literal seed, allowed and forbidden selection fields, exact stratum sizes, and inclusion probabilities are published in `selection_policy.json`. The full `frame.tsv` maps every one of the 26,077 predecessor candidates to eligibility, stratum, random key, cohort order (or `.`), first rung, and per-rung probability.

The current source-semantics release is consulted only to label controls. Its engineering gates are PASS, but its extraction status is `EXTRACTION_BLOCKED`; all future interval controls remain blocked and consumers must reject extraction.

## Tests

```bash
python -m unittest -v workflow.pilot_selection.test_selection
```

Failure tests cover pinned-manifest mismatch, resource allocation/refusal, exact-version and Stage-B drift, forbidden-trait independence, duplicate/nesting logic, and interrupted atomic promotion.

## Execution

The mandatory restart check is two automatic invocations using one run ID. The first intentionally dies after staging `SHA256SUMS` but before `COMPLETE`; no final release exists. The restart discards the incomplete release stage, reuses only the independently checksum-complete selection unit, then publishes atomically.

```bash
python workflow/pilot_selection/selection.py --repo-root . \
  --run-id pilot-cohorts-v1-execution --inject-kill-before-complete
# expected: SIGKILL / exit 137 and no final COMPLETE
python workflow/pilot_selection/selection.py --repo-root . \
  --run-id pilot-cohorts-v1-execution
```

Default explicit allocations are 4 GiB RAM, 1 GB durable, 4 TB scratch, and 100,000 inodes; predicted peaks are 100 MB per filesystem, 100 files, and 50 MB unfinished-write reservation. Live `findmnt`, ownership/write probes, bytes/inodes, allocations, reservations, RSS, and swap evidence are recorded before selection/resume and publication.

Validate the immutable external release after publication:

```bash
python workflow/pilot_selection/validate_release.py --repo-root . \
  --external-release /home/erikg/phind-data/ecoli26k/v1/releases/select-freeze-1k/<release_id>
```

Consumers reject a missing `COMPLETE`, inventory mismatch, non-PASS applicable gate, changed manifest byte/hash/count/order, or source-semantics use beyond post-selection controls.
