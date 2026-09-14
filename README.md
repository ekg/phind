# PHIND — ancestral and mean-field phage genome mining from *E. coli* prophages

**Objective.** Mine ancestral and mean-field phage genomes from the prophage
elements of *E. coli*. We do **not** build trees from whole bacterial genomes:
all work is on the ~132k prophage elements extracted from ~26k *E. coli*
assemblies (`prophage_homology_survey/full_prophages.fa`, headers
`GCF_<assembly>_prophage_N`).

**Pipeline at a glance.**

```
MASH sketch + triangle + tree          (done → research/mash_tree/)
  └─ cluster prophages                 (done → 5001 clusters / 12 communities)
       └─ per-cluster all-to-all alignment with fastga/sweepga     (TO RUN)
            └─ impg partition → small partitions (500–1000 bp max) (TO RUN)
                 └─ render partitions                              (TO RUN)
                      └─ optimal traversal / ordering of partitions
                         (rare partitions tolerated)               (CODE TO WRITE)
                           ├─ ancestral genome estimation          (code exists, needs repair/port)
                           └─ major-allele mean genome
                              (most-similar seq per partition via impg similarity) (TO WRITE)
```

This README documents the full plan — from initial MASH to ancestral genome
estimation — including what already exists, where it lives, what must still be
written, and the exact commands/tools to use.

---

## 1. Tools (verified installed)

| Tool | Path / binary | Role in pipeline |
|---|---|---|
| `mash` | system | sketching, distance triangle, identity checks |
| **`allwave`** | `~/.cargo/bin/allwave` (v0.1.0; src at `/home/erikg/allwave`, repo `pangenome/allwave`) | **THE all-vs-all wave aligner — biWFA pairwise alignment with sparsification, PAF out** |
| `impg` | `~/.cargo/bin/impg` | `partition`, `similarity`, `render`, `syng`, `map`, `refine`, `query`, `index` (NOT the aligner — see §4) |
| `FastGA` | `~/.cargo/bin/FastGA` | alternative pairwise aligner (not the historical wave path) |
| `sweepga` | via `impg graph` / `pggb-env` | graph induction from allwave PAF (the historical `allwave/sweepga` pairing) |
| `wfmash` | `~/.cargo/bin/wfmash` | mapper (used by pggb path; not the primary path) |
| `seqwish` | `~/.cargo/bin/seqwish` | graph induction (used by `impg graph`) |

Key `impg` subcommands for this project:

- `impg graph` — "Build a pangenome graph from FASTA sequences using sweepga+seqwish" (per-cluster alignment; accepts `--sequence-list` and emits GFA; `-a/--paf-file` skips realignment when a PAF exists).
- `impg partition --window-size <N>` — partition the alignment into homologous regions. **`--window-size` is the knob that controls partition size** (see §4 for the 500–1000 bp target).
- `impg similarity -a <paf>` — pairwise similarity between sequences in a region (basis for the major-allele mean genome, §7).
- `impg render --index <syng> --target-range <seq:start-end>` — render a syng-backed region into a GBZ-style bundle (partition rendering, §5).
- `impg syng` — build a GBWT syncmer index (fallback when all-to-all alignment is not feasible, §8).

---

## 2. Step 1 — MASH distance triangle + tree (DONE)

**Result.** `research/mash_tree/`:

- `full_prophages.msh` — MASH sketch of all 132,393 elements (k=21, s=10000; exact params in `research/mash_tree/COMMANDS.log`).
- `full_prophages_mash.dist` — full pairwise distance triangle, 8.76e9 pairs / 35.05 GB (git-ignored derived artifact). Verified: exact expected size, 50-pair spot-check vs direct `mash dist` → 0 mismatches (`triangle_verify.json`).
- `full_prophages_tree.nwk` / `full_prophages_tree_labeled.nwk` — UPGMA tree, 132,393 leaves, 0 dropped, ultrametric (`tree_verify.json`).
- `full_prophages_labels.csv`, `label_merge_stats.json` — community/cluster labels merged onto every leaf (100% coverage).
- `README.md`, `COMMANDS.log`, `scripts/` — reproducible methodology (sketch → chunked pairwise → sequential merge → UPGMA → labels).

**Why this step exists.** The MASH distance triangle is the input to clustering:
each prophage element gets a leaf, and similarity of prophage *content* drives
the clusters (not whole-genome relatedness of the host).

---

## 3. Step 2 — Clustering (DONE, two views)

- **5001 clusters**: `prophage_homology_survey/full_prophage_clusters.csv` (from `c6ce454`, `full_prophage_analysis.py`), with MDS coords in `full_prophage_mds_coords.csv`.
- **12 communities** (denser connected-component view): `prophage_homology_survey/full_heatmap_clusters.csv` (19,638 leaves in shared communities; 112,755 isolates), heatmaps in `prophage_homology_survey/*.png` + interactive/PDF outputs.

**Decision point.** Per-cluster work should use the community partition
(12 communities) as the primary grouping — the historical per-community runs
(`mean_genomes/`, §6) already used it. The 5001-cluster view is a finer grain
for small, highly-related groups; a cluster that is too large to align (see
§4) can be sub-clustered with MASH first.

---

## 4. Step 3 — Per-cluster all-to-all alignment + small partitions (TO RUN)

**Alignment — the wave step is `allwave` (`pangenome/allwave`), NOT `impg
align`.** AllWave is a high-performance pairwise aligner using **bidirectional
wavefront alignment (biWFA, WFA2-lib)** with sparsification and mash-based
orientation detection. Installed: `~/.cargo/bin/allwave` v0.1.0 (source at
`/home/erikg/allwave`).

```bash
allwave -i <community.fa> -o alignments.paf -p auto -t <threads>
# historical project usage (recovery records): allwave with k-nearest=5, k-farthest=2
allwave -i <community.fa> -o alignments.paf -p tree:5:2:0.1 -t <threads>
```

`allwave` options (installed binary): `-i/--input` (FASTA, gz ok), `-o/--output`
(PAF), `-s/--scores` (match,mismatch,gap_open,gap_ext[,gap_open2,gap_ext2];
default `0,5,8,2,24,1`, proven for 85–95% ANI), `-x/--preset` (ANI presets:
98%/90%/80%/70%/60%), `-t/--threads`, `-p/--sparsification`,
`--no-progress`, `--mash-matrix` (dump mash distances), `--wfa-orientation`,
`-k/--keep-prefixes`, `-e/--exclude-prefixes`.

Sparsification strategies (`-p`, default `giant:0.99`):

- `none` — all n² pairs (small sets only, < ~100 seqs)
- `auto` — giant-component connectivity at 95% probability (== `giant:0.95`)
- `random:<frac>` — deterministic-hash random fraction of pairs
- `giant:<prob>` — Erdős–Rényi graph kept connected w/ probability <prob> (edge prob ≈ (log n − log(−log x))/n)
- `tree:<near>:<far>:<random>[:<kmer>]` — k-mer similarity tree: <near> k-nearest,
  <far> k-farthest (stranger-joining), <random> fraction, kmer size (default 15).
  **Historical project config: k-nearest=5, k-farthest=2** (`tree:5:2:…`).

**"All-wave with auto sparsification" is the intended per-clade scaling
mechanism**: `-p auto` (or `tree:…` for the historical config) prunes the
pair space while keeping a connected alignment graph, so full all-vs-all
biWFA remains tractable across all 12 clades. The failed `run-all-wave` task
(agent-225, TUI: `failed · →589k ←42k ◎3.6M §2.4k, 13h`) is the prior attempt
at exactly this; any re-run must use `allwave` and record strategy+params per
community in the manifest.

**The allwave PAF then feeds graph induction (sweepga/seqwish) → `impg
partition`.** `impg align` exists but is a *different, newer* sparsified
pairwise aligner — it is NOT the historical allwave path and should not be
confused with it.

**Partitioning.** Apply `impg partition` to the alignment:

```bash
impg partition -a cluster.alignments.paf --window-size 1000 -o cluster.partitions
```

**Partition size target: 500–1000 bp maximum, smaller if needed.**
Historical partition sizes are **too large** for the current goal. Measured on
`community_3_partitions.bed` (archived `mean_genomes/`, §6):

```
n=7329  min=4  median=1000  max=5955  mean=1513
```

i.e. mean 1.5 kb, max ~6 kb — while the requirement is 500–1000 bp max.
`impg partition --window-size` must be tuned (e.g. 500) and the resulting
BED/MAF size distribution re-checked so that **no partition exceeds ~1 kb**.
Small partitions are what make the traversal step (§6) tractable: each
partition is a compact homologous block that can be ordered and consensused
independently.

**Historical artifacts that show the shape of this step** (archived, not on
main):

- `research_outputs/` in archived worktree agent-221 (commit `0c8373c`): `cluster_3.alignments.paf`, `cluster_3.syng.*` (syng index), `partitions/partition*.bed` + `partition*.maf` (55 partitions), `partition_run.log`, `partition_maf_run.log`, `full_pipeline.py` (heaviest-bundle consensus per partition), `analyze_partitions.py`, `compare_ancestral.py` / `compare_ancestral_v2.py`.
- The `mean_genomes/` directory (archived agent-223, commit `64ddebf`): per-community BEDs + MAF partition dirs for 11 communities (0,1,2,3,4,5,7,8,9,10,11).

**Import decision (from artifact audit).** These outputs are REUSABLE evidence
but are not canonical on `main`; the audit recommends re-deriving partitions
with the corrected (smaller) window size rather than importing the old,
too-large partitions.

**Caution — the historical "mean genomes" are naive concatenations, not
phage genomes.** The archived `mean_genomes_report.md` shows per-community
"mean genome" lengths of 2.8–4.5 Mbp (e.g. community 0: 28,301 partitions →
3,977,131 bp; community 3: 7,329 partitions → 1,398,681 bp). Real phages are
~50–150 kb. These outputs were produced by concatenating *every* partition,
including rare/noise blocks — exactly the failure mode the traversal step
(§6) exists to fix. Treat them as evidence of what not to do: the final
community genome must come from an *optimized typical traversal* of small
partitions, not a full concatenation.

---

## 5. Step 4 — Render the partitions (TO RUN)

For each partition, produce a rendered alignment bundle so the block can be
visually inspected:

```bash
impg render --index cluster.syng --target-range <seq_name>:<start>-<end> -O partition_render/
```

(For a non-syng path, render from the PAF/MAF directly via the scripts in
`research_outputs/` / `prophage_homology_survey/`.) Rendered partitions are
the unit of inspection: each shows which prophages cover the block and how
they align. This is also where rare partitions (present in only a few
genomes) become visible and are flagged for the traversal step.

---

## 6. Step 5 — Optimal traversal / ordering of partitions (CODE TO WRITE — MISSING)

**This is the central piece of code that does not yet exist.**

Problem: given the set of partitions for a community (many of which are rare —
present in only a few prophages, or missing entirely from most), find the
**best typical traversal**: sample the most common paths through the partition
set with probability proportional to each partition's occurrence/adjacency
support, so the resulting genome is the typical/majority-like (ML) genome.

**Rare partitions are sampled, not excluded and not forced in** (user
correction 2026-08-03): every partition gets a nonzero sampling chance
weighted by its support — common partitions dominate the traversal, rare
modules can still appear occasionally (biologically relevant), but they do
not bloat the genome. Two earlier designs were explicitly rejected: (v1)
deterministically retaining rare partitions by bridging, and (v2) excluding
partitions below a hard occurrence threshold.

What exists today (partial, needs extension):

- `research/stitching/stitch_algorithm.py` (commit `2363ece`, **validated
  byte-identical**): builds an adjacency graph from consecutive partition
  pairs across all prophages (`build_adjacency_graph`), finds a maximum
  likelihood path (`find_maximum_likelihood_path`), computes a majority-rule
  consensus per partition (`compute_partition_consensus`), and stitches
  (`stitch_and_merge` with overlap detection). Validation: `research/stitching/validation_report.md` — 50% threshold → 1 core partition (53,886 bp); 45% threshold → 124,935 bp at 78.05% MASH identity to the ancestral genome.

Gaps to fill in the new traversal code:

1. **Weighted sampling of rare/absent partitions** — current code drops partitions
   below a coverage threshold; the new code should *sample* the traversal by
   probability proportional to partition weight (occurrence^alpha), so rare
   partitions appear with low (nonzero) probability rather than being
   excluded or forced in.
2. **Optimization objective** — sample the most common paths with probability
   proportional to their support (a weighted path-sampling problem over the
   partition adjacency graph; deterministic greedy ML path exists, but the
   weighted sampling with an explicit objective + validation metric is missing).
3. **Per-partition representative choice** — for each partition in the
   traversal, pick the representative sequence (majority-rule consensus today;
   see §7 for the most-similar-to-all alternative).
4. **Reproducibility** — deterministic ordering, documented params, and a
   stats/report output per community.

Deliverable: `scripts/traverse_partitions.py` consuming
`partitions/*.bed` + `partitions/*.maf` and emitting (a) the sampled
traversal order(s) with per-partition weights, (b) per-partition
consensus/representative sequence, (c) the stitched community genome per
sample, (d) coverage + sampling statistics per partition (occurrence,
fraction, weight, sampling frequency). CLI:
`traverse_partitions.py --partitions-dir <dir> --bed <bed> --output <prefix> [--n-samples N] [--alpha A] [--seed S]`.

---

## 7. Step 6 — Ancestral genome estimation + major-allele mean genome

Two-level ML inference (this is the core of the traversal deliverable,
`scripts/traverse_partitions.py`):

**(a) ML across partitions (Level 1).** Weighted path sampling over the
partition adjacency graph: each partition is sampled with probability
proportional to its occurrence/adjacency support (weight = fraction^alpha),
so the traversal concentrates on the most common paths while rare partitions
are **sampled, just less likely** (never excluded, never forced in).

**(b) ML within each partition (Level 2).** For each sampled partition, pick
the maximum-likelihood representative of the aligned block (majority-rule
consensus over the MAF, or the sequence closest to the consensus).

**(c) Ancestral twist (alternate mode).** Same Level-1 traversal, but Level 2
is ancestral state reconstruction within each partition (per-partition
NJ/ML tree + root reconstruction, parsimony or likelihood), yielding an
**ancestral genome** per clade as a distinct output from the ML genome.

Deliverables per clade: ML genome + ancestral genome (mode flag
`--mode ml|ancestral`), sampled traversal orders with per-partition weights,
coverage + sampling statistics (occurrence, fraction, weight, sampling
frequency).

Two complementary outputs per community (and finally, an overall synthesis):

**(a) Ancestral genome.** Majority-rule consensus per partition, ordered by the
traversal (§6), stitched into one sequence. Partial implementation exists in
`research/stitching/stitch_algorithm.py` (MSA majority vote per partition) and
in archived `research_outputs/full_pipeline.py` (heaviest-bundle consensus).
Historical validation: community 3 stitched at 45% threshold = 124,935 bp,
78.05% identity vs the pggb ancestral genome (which itself came from
`pggb_analysis/cluster_6/ancestral/`, 69,153 bp, 8 modules, conf 0.5523 —
pggb path now deprecated, use fastga/sweepga + impg instead).

**(b) Major-allele mean genome — most similar to all.** The genome formed by
taking, for each partition, the **most similar sequence** (not the consensus).
This can be done with `impg similarity`:

```bash
impg similarity -a cluster.alignments.paf --region <partition_bed> ...
```

per partition, pick the element with the highest mean similarity to all other
elements of the partition; concatenate the chosen elements in traversal order.
This yields the observed "mean-field" phage genome — the real sequence most
representative of the community — vs the constructed consensus in (a).

Deliverable: `scripts/mean_field_genome.py` producing
`<community>_meanfield.fa` (+ a comparison table between consensus and
mean-field genomes: identity, length, module content).

---

## 8. Fallback — impg syng if all-to-all alignment is not feasible

For communities too large/divergent for a single all-to-all alignment,
`impg syng` (GBWT syncmer index) is the fallback:

```bash
impg syng -i cluster.fa -o cluster.syng          # build index (done historically for cluster_3)
impg map -i cluster.syng ...                     # map sequences to index
impg partition -i cluster.syng ...               # partition via index
impg refine ...                                  # refine loci to maximize sample span
```

For prophages and fastga/sweepga, direct all-to-all alignment should be
feasible; syng is the escape hatch if a cluster blows up. A historical syng
index already exists: `research_outputs/cluster_3.syng.*` (archived agent-221).

---

## 9. Status summary

| Step | Status | Location |
|---|---|---|
| MASH sketch + triangle + UPGMA tree | ✅ done | `research/mash_tree/` |
| Clustering (5001 / 12 communities) | ✅ done | `prophage_homology_survey/` |
| ECOR mapping + highlight + inspection index | ✅ done | `research/ecor/` |
| Artifact audit (reuse vs failed) | ✅ done | `research/artifact_audit.md` |
| Stitching validation (`2363ece`) | ✅ done (REUSABLE) | `research/stitching/` |
| **`allwave` all-wave + sparsification per clade (`-p auto` / historical `tree:5:2:…`)** | ⏳ **paused — regrouped; must use allwave, not impg align/pggb** | `run-all-wave` failed historically (agent-225) |
| `impg partition` at ≤1 kb windows | ⏳ to run (old partitions too large) | archived `mean_genomes/`, `research_outputs/` |
| Partition rendering | ⏳ to run | `impg render` |
| Traversal / ordering optimizer (rare partitions) | ❌ **code to write** | `scripts/traverse_partitions.py` |
| Ancestral genome estimation | 🟡 partial | `stitch_algorithm.py` + archived `full_pipeline.py` |
| Major-allele mean-field genome (impg similarity) | ❌ **code to write** | `scripts/mean_field_genome.py` |

---

## 10. Working rules (from recovery)

- Do **not** modify `/home/erikg/phind.recovery-20260801T164014Z` or
  `/mnt/nvme3n1/erikg/phind-recovery-20260801T164014Z` (read-only evidence).
- Do not push/rewrite git history; oversized artifacts stay git-ignored
  (3.1 GB `full_prophages.fa`, 35 GB `full_prophages_mash.dist`).
- Do not trust historical task terminal labels — verify artifacts (see
  `research/artifact_audit.md`).
- Small active WG graph only; `.wg` stays out of git.
