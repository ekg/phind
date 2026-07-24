# Canonical consumer compatibility certification

This workflow consumes only the immutable ten-assembly release
`canonical-cohort-010-v1-e71484de9994fc28`. It verifies the tracked release,
every external `SHA256SUMS` entry, all 80 object/index rows, decompressed
canonical content, exact cohort order, 1,223 PanSN paths, 51,731,662 bases, the
global 1,000-assembly cap, and both immutable root-input hashes before invoking
a consumer.

## Pinned environments

- `environment-linux-64.explicit.lock` freezes exact package URLs/builds with
  repository MD5s (the format consumed by micromamba 1.5.9), while
  `environment-package-sha256.tsv` independently pins SHA-256 for every package.
  It provides Mash 2.3, RapidNJ 2.3.2, skani 0.3.1, QUAST 5.3.0, gffread 0.12.7,
  Prodigal 2.6.3, MMseqs2 18.8cc5c, HMMER 3.4, and MCL 22.282. The environment
  is installed below task scratch and deleted after validation; no system
  package is changed.
- `graph-linux-64.explicit.lock` and `graph-package-sha256.tsv` freeze the
  already installed pggb 0.6.0 graph environment
  (pggb/wfmash/seqwish/smoothxg/odgi/vg). Every executable is also
  pinned by SHA-256 in the release.
- System bgzip/samtools and installed IMPG 0.4.1 are consumed only at their
  captured executable SHA-256. IMPG build provenance cannot be reconstructed
  beyond those bytes, so the release never invents a source commit.

`environment.yml` is a readable intent file. The explicit lock, not a fresh
solve, is the execution authority.

## Input/view rules

Canonical per-assembly BGZF is direct for bgzip, samtools, Mash, skani and
QUAST. IMPG SYNG uses a streamed, all-ten combined BGZF because `impg syng -f`
accepts one FASTA path, not a list. pggb uses a bounded two-path BGZF fixture
derived from one authorized canonical path. Prodigal 2.6.3 unexpectedly but
reproducibly accepts that BGZF directly and matches the plain-input gene/protein
semantics; MMseqs2/HMMER/MCL consume its bounded protein derivative. gffread compares
seqids lexically, so a strict adapter converts standards-compliant GFF `%23`
seqids to semantic FASTA `#` values in column 1 only while recording a
reversible map.

Every staged view records sources, order/range, checksum, quota, reversible
name mapping and cleanup. All sequence/index-bearing views and the installed
host-tool environment remain below the task scratch namespace and are removed
on success or failure.

## Execution

The following allocations were used for the certification (decimal disk units):

```bash
DURABLE=/home/erikg/phind-data/ecoli26k/v1/releases/certify-pilot-consumer
SCRATCH=/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/certify-pilot-consumer/compat-v1
COMMON="--durable-task-root $DURABLE --scratch-root $SCRATCH \
 --assigned-ram-bytes 8589934592 --durable-allocation-bytes 5000000000 \
 --scratch-allocation-bytes 4000000000000 --inode-allocation 100000 \
 --predicted-durable-peak-bytes 1000000000 \
 --predicted-scratch-peak-bytes 3000000000 --predicted-files 25000 \
 --unfinished-write-bytes 500000000"

python -m workflow.compatibility.pilot interrupt-test $COMMON
python -m workflow.compatibility.pilot run $COMMON --run-id compat-v1
```

Each stage records live `findmnt`, ownership/write probe, bytes/inodes, explicit
allocations and swap. Blank allocations or threshold/floor failures are
`NO_GO`. The interrupt test uses a real SIGKILL: the partial staging directory
cannot publish, is cleaned, and only a fresh checksum-complete directory is
atomically promoted.

## Validation

```bash
python -m unittest -v workflow.compatibility.test_compatibility
python -m workflow.compatibility.validate_release \
  --external-release /home/erikg/phind-data/ecoli26k/v1/releases/certify-pilot-consumer/consumer-compatibility-v1-RELEASE
```

The validator requires one unqualified `PASS` machine gate per selected
consumer, exact inventory coverage, immutable predecessor/root digests, compact
Git ownership, resource thresholds, cleanup, and kill/restart evidence. Missing,
`BLOCKED`, `CONDITIONAL`, or substituted inputs fail closed.

Mash+RapidNJ is certified only as an **unrooted genomic-similarity
dendrogram**, never a phylogeny. IMPG `query` and `map` remain distinct: query
projects an indexed coordinate range; map emits a syncmer-anchor projection,
not a base alignment.
