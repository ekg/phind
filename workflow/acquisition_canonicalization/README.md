# Stage-B acquisition and canonicalization

`pilot.py` is the bounded production workflow for the immutable ten-row Stage-B manifest. It can construct payload URLs only for the ten exact assembly revisions in `manifests/collection-v1/stage_b_10.tsv`; a versionless, substituted, or eleventh accession is rejected before I/O.

## Contracts

- Verifies the tracked and external predecessor release, `release.json` SHA-256, exact Stage-B bytes/SHA-256/order/row hashes, and both immutable root inputs before acquisition.
- Runs live `findmnt`, ownership/write, byte, inode, explicit allocation, unfinished-write, RAM, and swap gates before every acquisition/conversion stage.
- Downloads one NCBI Datasets v2 package at a time with bounded retry/rate limiting. A partial is range-resumed only when strong remote identity and byte-range support still match; otherwise it is discarded and safely restarted.
- Validates ZIP structure/CRC, complete upstream `md5sum.txt` coverage, exact accession directory/catalog identity, and local SHA-256 before atomically committing a source object.
- Streams the package FASTA directly through the versioned rename policy into `bgzip`; no extracted or cohort-wide plain FASTA is retained. Compression restarts rather than appends after an interruption.
- Publishes exact `assembly#1#contig` PanSN IDs using byte-reversible uppercase percent encoding (or the policy's digest alias), `.fai`, `.gzi`, per-contig rename-only digests, and annotation alias tables only after GFF3 seqid and 1-based-closed bounds validation. The source GFF remains unchanged inside its checksum-validated package; no transformed GFF is claimed.
- Commits each source/canonical object with `COMPLETE` last, then commits the entire external release by same-filesystem rename. Consumers must verify `COMPLETE` and `SHA256SUMS`.

The production release is deterministic in identity:

```text
canonical-cohort-010-v1-e71484de9994fc28
```

under:

```text
/home/erikg/phind-data/ecoli26k/v1/releases/run-10-assembly-acquisition/
```

The approved run scratch namespace is:

```text
/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-10-assembly-acquisition/stage-b-010-v1/
```

## Invocation

The run used explicit allocations of 8 GiB RAM, 10 GB durable, 4 TB scratch, and 100,000 inodes; predicted durable/scratch peaks were 500 MB each, 500 files, and a 500 MB unfinished-write reservation.

```bash
python -m workflow.acquisition_canonicalization.pilot run \
  --durable-task-root /home/erikg/phind-data/ecoli26k/v1/releases/run-10-assembly-acquisition \
  --scratch-root /mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-10-assembly-acquisition/stage-b-010-v1 \
  --run-id stage-b-010-v1 \
  --assigned-ram-bytes 8589934592 \
  --durable-allocation-bytes 10000000000 \
  --scratch-allocation-bytes 4000000000000 \
  --inode-allocation 100000 \
  --predicted-durable-peak-bytes 500000000 \
  --predicted-scratch-peak-bytes 500000000 \
  --predicted-files 500 \
  --unfinished-write-bytes 500000000
```

The first two invocations add, respectively, `--inject-kill acquisition` and `--inject-kill conversion`; each self-sends `SIGKILL` only after flushing its partial and append-only event. The final invocation omits injection and must observe both safe-restart paths before publication.

## Validation

```bash
python -m unittest -v workflow.acquisition_canonicalization.test_pilot
python -m workflow.acquisition_canonicalization.validate_release \
  --external-release /home/erikg/phind-data/ecoli26k/v1/releases/run-10-assembly-acquisition/canonical-cohort-010-v1-e71484de9994fc28 \
  --output artifacts/acquisition_canonicalization_10/validation.json
```

Tests cover manifest/checksum mismatch, blank/over-allocated resources, accession-cap refusal, archive/upstream checksum corruption, strong-identity range resume, unsafe partial restart, streamed BGZF/PanSN/GFF alias round-trip, interrupted conversion, and interrupted atomic promotion. Re-running `pilot.py` against the immutable complete release performs no network request and republishes only checksum-verified compact tracked views.
