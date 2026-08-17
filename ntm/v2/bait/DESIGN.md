# NTM prophage bait & junction panel — preregistered design (v1)

**Task.** `design-ntm-prophage` · **Date.** 2026-08-17 · **Status.** preregistered
before the pilot run; thresholds below are fixed before looking at bait outcomes.

Purpose: a **reproducible bait panel** for bounded Logan/SRA validation of the
NTM v2 reconstructed ML phage genomes (task `execute-bounded-ntm`). Logan is
queried with **31-mers**, so every bait is designed, masked and reported at
k = 31 granularity.

## Inputs (all checksummed into `provenance.json`)

| Input | Path (external root `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/`) | Use |
|---|---|---|
| Functional QC per genome | `annotation/report/per_genome_functional_qc.tsv` | pilot eligibility (tier/cohort/flag) |
| Release manifest | `ml/release_manifest.tsv` | clade → status, length, flags |
| ML genome FASTA | `ml/all_ntm2_ml_phage_genomes.fa` | bait sequence source |
| Prophage distance matrix | `clades/0/distances.npz` (+ `members.json`) | **diversity-aware selection (phage-side only)** |
| Tight clades | `clades/0/tight_clades.json` | clade membership |
| Traversal (junctions) | `clades/<cid>/ml.traversal.json`, `clades/<cid>/partitions.bed` | partition joins, breakpoints, observed adjacency |
| Member prophages | `full_prophages.fa` | observed member-prophage representatives |
| Collaborator prophage master manifest | `inputs/NTM_QC_passed_prophage_master_manifest.tsv` | prophage → accession/contig/coords; host intervals |
| Pharokka CDS table (ML genomes) | `annotation/pharokka_out/pharokka_cds_final_merged_output.tsv` | module-gene targeting; IS/plasmid/AMR/rRNA composition |
| Host assemblies (PanSN bgzip) | `genomes/canonical_objects/<acc>/<acc>.pansn.fa.gz` | host-like 31-mer masking; host control baits |

Bulky intermediates, ledgers and logs stay external under
`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/bait/`. Only the small panel
FASTA + manifests + this design are committed in-repo (`ntm/v2/bait/panel/`).

## Stage 1 — Pilot selection (deterministic, diversity-aware, host-blind)

**Eligibility** (all must hold):

1. `cohort == ntm2_ml_reconstructed` (release `status=ml` — excludes all
   singleton pass-through fragments by construction);
2. functional-QC tier ∈ {`A_strong_candidate`, `B_moderate_candidate`};
3. functional-QC `flag == ok`;
4. release length ≥ 5,000 bp (enough for interior + junction baits);
5. release flags exclude `too_short` and `long_n_run` (mosaic flags
   `low_completeness`/`low_identity` are **allowed** — they measure the
   expected mosaic property, not bait usability);
6. clade has ≥ 2 members (implied by `status=ml`).

**Host-blindness rule.** Host clade ids and species labels are **never**
selection features. They are recorded in the selection manifest purely as
provenance. Diversity is measured in **prophage mash space** only.

**Selection algorithm** (fixed, no randomness):

- For each eligible clade, the representative is its **medoid member**
  (member minimizing mean mash distance to co-members; ties → lexicographic
  prophage id). Missing (NaN) distances — mash triangles without shared
  hashes — count as 1.0.
- **Tier-A pool:** greedy farthest-point sampling over medoid-to-medoid
  distances, seeded with the tier-A candidate having the lexicographically
  smallest clade id; each next pick maximizes the minimum distance to all
  already-selected representatives (ties → smaller clade id). Stop at 15.
- **Tier-B pool:** same farthest-point procedure over tier-B candidates but
  measured against **all** already-selected representatives (A picks first),
  stop at 5. → **20 ML genomes**.
- If a pool is short, fill from the other tier; every deviation is recorded.

**Member-prophage representatives.** For each selected clade, up to 2 observed
members: the medoid member and the member farthest from the medoid
(diversity), requiring member length ≥ 600 bp. Deterministic tie-breaking by
prophage id.

**Output:** `panel/pilot_selection.tsv` (committed, 20 rows + member
representatives) and the full eligibility ledger (external,
`work/eligibility_ledger.tsv`) giving every candidate an explicit
include/exclude reason.

## Stage 2 — Bait classes and coordinates

All baits are **0.5–2.5 kb**, 1-based inclusive coordinates on the named
source molecule (forward orientation; masking handles both strands).

| Class | Source | Window | Count cap |
|---|---|---|---|
| `interior_module` | ML genome | 1,200 bp centered on a key-module CDS midpoint (min of start/stop to max), priority: terminase large → portal → major head/capsid → tail tape measure → endolysin → integrase; module genes located in Pharokka CDS table (`annot` regex + PHROG category, mirroring the functional-QC module definitions); overlapping windows deduplicated keeping higher priority | 3 / genome |
| `interior_generic` | ML genome | 1,200 bp at evenly spaced offsets (25%, 55%, 85% of genome length), only used to fill a genome to 3 interior baits | fill to 3 |
| `junction` | ML genome | 1,000 bp centered on the exact breakpoint between consecutive stitched partitions (500 bp each side; clamped at genome ends and shifted inward to preserve width) | 3 / genome |
| `member_interior` | observed member prophage (`full_prophages.fa`) | 1,200 bp centered at prophage midpoint (width = min(1200, member length); skipped if < 500) | 1 / member representative |
| `control_shuffled_negative` | derived | fixed-seed (42) first-order Markov sequence fitted to the first accepted `interior_module` bait of 3 distinct genomes: matches length and mono-/dinucleotide composition, destroys 31-mer content; gated on ≤ 5 canonical 31-mers shared with the whole panel (deterministic redraws) with composition deviation reported | 3 total |
| `control_host_negative` | source host assembly | 1,000 bp window centered in the largest prophage-free stretch of the prophage-bearing contig of a source host (3 distinct clades; deterministic tie-breaks); expected to be flagged host-like by our own filter — doubles as a **filter validation control** | 3 total |
| `control_public_reference` | reserved | curated public mycobacteriophage references (sibling task `curate-public-mycobacteriophage`); **0 rows in v1**, class reserved so the panel schema is forward-compatible | pending |

**Junction provenance.** Every junction bait records both contributing
partition ids (`partition_a` = upstream, `partition_b` = downstream), the exact
breakpoint coordinate (start of `partition_b` in the ML genome = cumulative
`representative_len` of preceding pids), each partition's occurrence count and
ML-genome coordinates, `co_occurrence_count` (members carrying both
partitions, from `partitions.bed`) and `observed_adjacency_count` (members in
which the two partitions are exactly contiguous). `adjacency_observed = false`
marks joins that exist only in the reconstruction — the primary scientific
target of the panel.

Junction ranking: highest `min(occurrence_a, occurrence_b)` first (best-supported
joins), then lower partition id; both partitions must have sequence.

## Stage 3 — 31-mer masking and bait QC gates

k = 31, canonical form = min(code, reverse-complement code) with exact 62-bit
2-bits-per-base encoding (no hashing collisions). Per-bait statistics over all
positional 31-mers:

| Mask | Rule |
|---|---|
| `ambiguous` | 31-mer contains any non-ACGT base |
| `low_complexity` | 31-mer overlaps a ≥ 7-bp homopolymer run **or** base-composition Shannon entropy < 1.2 bits |
| `host_like` | canonical 31-mer present in the **host panel** (see below) |
| `internal_dup` | canonical 31-mer occurs > 1× within the bait; all but the first occurrence masked |

**Host panel.** Source host genomes of the selected clades' members and
representatives (the realistic host-contamination risk), ≤ 4 accessions per
clade, ≤ 60 total, deterministic order (member-representative sources first,
then medoid sources, then accession). Host sequences are scanned **with every
collaborator-called prophage interval on that accession masked out** (from the
master manifest), so prophage-derived 31-mers are not spuriously host-like.
Panel composition + per-file SHA-256 recorded in `provenance.json`.

**Region-level rejection** (a bait, not a k-mer, property): a window is
*dominated* by non-phage markers when ≥ 50% of its CDS-covered bases (with ≥
30% CDS coverage overall) are CDS matching the reject-marker regex
`transposase|insertion|plasmid|resistance|ribosomal|rRNA|tRNA` or carrying a
CARD/VFDB hit, or the PHROG category `DNA, RNA and nucleotide metabolism`
(dominated only — normal module baits are unaffected). Reject with reason
`dominated_nonphage_markers`.

**Acceptance gates** (checked in fixed order; first failure is the recorded
reason): length ∈ [500, 2500]; `dominated_nonphage_markers` false;
ambiguous_frac ≤ 0.10; host_like_frac ≤ 0.10; low_complexity_frac ≤ 0.30;
usable (unmasked, unique) 31-mers ≥ 200; exact-sequence duplicate of an
already-accepted bait → `duplicate_sequence`; > 50% of the bait's canonical
31-mers shared with an already-accepted bait → `cross_bait_kmer_dup` (higher
priority baits kept). `control_host_negative` baits are exempt from the
host_like gate and from the host-mask subtraction in the usable-kmer
computation by design (`expected_behavior=host_like`) — their actual
host_like fraction must be reported and is itself a filter validation
signal.

**Rejection ledger:** every rejected candidate bait (and every masked-out
module window) is recorded with bait id, class, source coordinates, failure
reason and 31-mer statistics — no silent drops.

## Stage 4 — Outputs & provenance

Committed in-repo (`ntm/v2/bait/panel/`, byte-deterministic reruns):

- `baits.fa` — `>NTMBAIT_<clade>_<CLASS>_<nn>` headers, one line per record;
- `pilot_selection.tsv` — 20 selected genomes + member representatives,
  with host-clade/species provenance columns explicitly marked not-used;
- `bait_manifest.tsv` — full per-bait provenance: id, class, genome/clade,
  member prophage (member baits), contig, start/end/length, module gene
  (phrog/annot/coords) or junction fields (partitions, breakpoint,
  occurrences, co-occurrence/observed adjacency), sequence SHA-256, and all
  31-mer statistics (total/ambiguous/low_complexity/host_like/internal_dup/
  cross-bait-shared/usable), GC, acceptance status;
- `rejection_ledger.tsv`, `filter_stats.tsv` (mask composition summary),
  `provenance.json` (input SHA-256s, host-panel checksums, parameters,
  deterministic `run_id` = SHA-256 over inputs+parameters — no wall-clock in
  committed outputs), `README.md`.

External (`…/ntm/v2/bait/`): eligibility ledger, host-panel staging list and
k-mer scan log, run log with timestamps.

## Validation contract (task `## Validation` mapping)

- Deterministic, diversity-aware pilot of reconstructed-only A/B candidates;
  controls separately labeled (`control_*`) — tested by rerun equality and by
  unit test on the greedy selector.
- Every bait 0.5–2.5 kb with unique ID + SHA-256, source coordinates, class
  and 31-mer stats; ambiguous bases and duplicates reported — enforced by a
  panel self-audit step that fails the run on any violation.
- Junction baits record both partitions + exact breakpoint; negative controls
  included (shuffled + host-derived).
- Host/plasmid/repeat filtering reproducible: input + host-panel checksums in
  provenance; rejection ledger complete.
- Unit tests (`scripts/test_bait_design.py`): reverse complements, boundary
  coordinates, duplicated 31-mers, masking (ambiguous/low-complexity/host-like),
  deterministic selection, window clamping, Markov-control composition gates.

## Explicit caveats

1. Host masking covers **source-host assemblies outside called prophages** —
   not all 13,122 hosts (bounded, checksummed panel). Residual host-like
   31-mers from unsampled hosts remain possible; the pilot queries Logan at
   accession level anyway and downstream confirmation is by alignment.
2. IS/plasmid/AMR/rRNA filtering is annotation-driven (Pharokka/PHROG + CARD/
   VFDB columns); unannotated mobile elements are not caught — the
   low-complexity/uniqueness gates are the second line of defense.
3. Module-gene regexes mirror the functional-QC definitions; a module bait is
   *centered on* a module gene, not an amplicon of the gene alone.
4. `control_public_reference` baits are deferred to
   `curate-public-mycobacteriophage`; the panel schema accepts them without
   redesign.
5. Species/host-clade labels on baits are collaborator calls recorded as
   provenance only; they are never used as evidence of phage relatedness
   (host-blind selection rule above).
