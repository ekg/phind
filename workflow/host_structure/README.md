# Frozen N=1,000 host-only structure workflow

This workflow consumes only the digest-pinned frozen pilot selection, exact
canonical N=1,000 assembly release, certified host consumer environment, and
opaque root-file hashes. It never parses the prophage CSV and rejects any
phage-derived biological key or path from its run/release manifests.

## Methods

1. Verify exact release IDs, SHA-256 pins, immutable cohort order, 1,000/1,000
   rows, object/index checksums, global cap evidence, root hashes, and every
   applicable predecessor gate.
2. Record a live resource preflight before each stage/batch. Materialize only
   validated canonical BGZF objects into task scratch, checking decompressed
   content, FAI/GZI, PanSN names, contig count, and bases.
3. Run pinned Mash 2.3 whole-assembly sketches (`k=21,s=1k/10k/50k` and
   `k=31,s=10k`), exact 499,500 lower-triangle pairs, all 1,000,000 directed
   pairs for symmetry/matrix validation, pinned RapidNJ 2.3.2, and six
   independent sketch-seed trees. Mash outputs are genomic-similarity
   dendrograms, not substitution-model phylogenies.
4. Choose medoids, diverse/boundary cases, and QC extremes only from host
   distances and assembly QC. Build reference-coordinate core alignments with
   digest-pinned minimap2, deterministic SNP-density recombination-candidate
   masks, core-SNP NJ trees, 100 site bootstraps, support collapse, and an
   alternative-reference topology check in each of 16 limited-diversity
   sampling partitions. The heuristic mask is not misrepresented as Gubbins.
5. Keep trees unrooted because no independently verified outgroup occurs in
   the frozen cohort and acquisition outside it is blocked. Freeze only
   parameter- and seed-stable unrooted Mash splits; leave every other all-host
   membership ambiguous and retain k=12/16/20 host-genetic alternatives.
6. Write `COMPLETE` last and atomically rename the external staging directory.
   Alignment/sketch/matrix/PAF/plain FASTA payloads never enter git.

## Commands

```bash
python -m unittest -v workflow.host_structure.test_host_structure

# Required real SIGKILL attempt; expected shell status is 137/9 and no release.
python -m workflow.host_structure.runner --inject-kill materialize

# Same run ID resumes only checksum-valid views and publishes atomically.
python -m workflow.host_structure.runner --inject-kill none

python -m workflow.host_structure.validate_release \
  --external /home/erikg/phind-data/ecoli26k/v1/releases/run-host-structure-1000/<release_id>
```

The runner defaults to a declared 64 GiB RAM allocation, 30 GB durable
allocation, 4 TB scratch allocation, 500,000 inodes, 10/20 GB durable/scratch
upper-95% predictions, 10,000 projected files, and a 5 GB unfinished-write
reservation. Blank or threshold-violating allocations are hard failures.
