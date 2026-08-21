# E. coli prophage bait & junction panel — preregistered design (v1)

**Task.** `design-e-coli` · **Date.** 2026-08-21 · **Status.** preregistered
before any screening; thresholds below are fixed before looking at bait
outcomes. This design mirrors `ntm/v2/bait/DESIGN.md` (NTM v1) and reuses the
generalized bait pipeline (`ntm/v2/bait/scripts/baitlib.py`) delivered by
`implement-bait-design`.

Purpose: a **reproducible, frozen bait panel** over the 585-genome *E. coli*
ML phage release (`research/ml_phage_genomes/`) for the stratified Logan
SRA screen (task `local-e-coli`). Logan is queried with **31-mers**, so every
bait is designed, masked and reported at k = 31 granularity.

## Inputs (all checksummed into `provenance.json`)

| Input | Path | Use |
|---|---|---|
| ML genome FASTA (585) | `research/ml_phage_genomes/all_ml_phage_genomes.fa.gz` | bait sequence source + traversal verification target |
| Release manifest | `research/ml_phage_genomes/release_manifest.tsv` | clade → status, length, members |
| Release QC table | `research/ml_phage_genomes/qc_table.tsv` (+ `qc_flagged.tsv`) | release-flag exclusion (`too_short`, `long_n_run`) |
| Clade → source genomes | `research/ml_phage_genomes/clade_genomes.json` | provenance |
| Functional QC (per genome) | `research/phage_annotation/per_genome_annotation_qc.tsv` (`source=ecoli_ml`; identical sha256 to the NVMe annotation-report copy) | pilot eligibility / tier derivation |
| Tight clades | `research/clades/<community>/tight_clades.json` | clade membership (prophage ids) |
| Prophage MASH triangle | `research/mash_tree/full_prophages_mash.dist` (float32 condensed, 132,393 ids) + `data/ids.txt` | **diversity-aware selection (phage-side only)**, medoids |
| Partition data (per clade) | main tree `research/clades/<community>/<clade>/{partitions.bed,partitions/,manifest.json}` (staged to NVMe) | traversal regeneration → junction provenance |
| Traversal tool | `scripts/traverse_partitions.py` (in-repo, seeded) | deterministic `ml.traversal.json` regeneration |
| Member prophages | `prophage_homology_survey/full_prophages.fa` (main tree; selected ids staged to NVMe) | member-prophage representative baits |
| Prophage calls | `26k_prophage1.csv` (132,404 calls: accession/scaffold/begin/end) | host prophage-interval masking, member coordinates |
| Pharokka CDS table | `/mnt/nvme3n1/…/annotation/pharokka_out/pharokka_cds_final_merged_output.tsv` (219,370 CDS; contigs `ecoli_ml\|clade_<cid>_ML`) | module-gene targeting; non-phage-marker domination |
| Host genomes (PanSN bgzip) | `/mnt/nvme3n1/…/ecoli26k/v1/26k/canonical_objects/<acc>/<acc>.pansn.fa.gz` (26,074 E. coli) | host-like 31-mer masking; host negative controls |
| Public references | `panel/public_refs/{NC_001416.1.fa (λ, 48,502 bp), NC_002167.1.fa (HK97, 39,732 bp)}` fetched from NCBI E-utilities | positive controls; known-phage calibration |

Bulky intermediates, staged inputs, ledgers and logs stay external under
`/mnt/nvme3n1/erikg/phind-genome-work/ecoli_bait/`. Only the small panel
FASTA + manifests + public-reference slices' sources + this design are
committed in-repo (`research/ecoli_bait/panel/`).

## Junction provenance rule (STOP-and-report clause)

NTM's junction baits carry `observed_adjacency_count` from per-clade
`ml.traversal.json` files. **No traversal JSON files were retained for the
E. coli ML release.** Rather than approximating, this design **regenerates
them deterministically** with the exact released parameters and **verifies
each regenerated genome byte-identical** to the frozen release:

- parameters (from `research/ml_phage_genomes/RELEASE.md`):
  `--mode ml --n-samples 5 --seed 42 --alpha 1 --adj-weight 1
  --coverage-threshold 0.25 --max-length max(30000, int(1.6 × max member
  length))` (max member length from the clade `manifest.json`);
- verification: regenerated `<cid>.ml.fa` sequence **must equal** the
  `clade_<cid>_ML` record of `all_ml_phage_genomes.fa.gz` byte-for-byte;
  verified pre-registered on 4 clades spanning the size range (845 bp →
  106,696 bp, incl. a 100-member clade) — all matched;
- a clade whose regeneration does **not** reproduce the released genome
  gets **zero junction baits** (ledger reason `traversal_regeneration_mismatch`)
  — no approximated junctions are ever emitted; every other bait class for
  that clade is unaffected;
- if traversal inputs are absent for a selected clade (no `partitions.bed`),
  the same zero-junction rule applies (`traversal_inputs_missing`), loudly
  recorded in `provenance.json`.

## Stage 1 — Pilot selection (deterministic, diversity-aware, host-blind)

**Eligibility** (all must hold):

1. functional-QC row exists (`source=ecoli_ml`) and derived tier ∈ {A, B}
   (tier derivation below);
2. release `status == ml` (excludes all 106 singleton pass-through fragments
   by construction);
3. release length ≥ 5,000 bp; `n_members` ≥ 2;
4. release QC flags exclude `too_short` and `long_n_run` (mosaic flags are
   allowed — they measure the expected mosaic property).

**Tier derivation (fixed before outcomes).** The E. coli functional-QC
`flag` is a `;`-joined composite; tiers mirror the NTM A/B intent:

- **Tier A (strong candidate):** `flag == ok`;
- **Tier B (moderate candidate):** flags ⊆ {`host_genes`, `low_completeness`},
  non-empty, ≠ `ok` — i.e. every clade whose only defects are host-gene
  carriage and/or completeness below the reporting bar, with **no**
  `missing_core_structural`, `many_truncated` or `contamination`;
- anything else is ineligible.

Measured at preregistration: pool A = 70 clades, pool B = 279 clades.

**Host-blindness rule.** Host/phylogroup labels are never selection
features; they are recorded in the selection manifest purely as provenance.
Diversity is measured in **prophage mash space** only.

**Selection algorithm** (fixed, no randomness; identical to NTM):

- representative per clade = **medoid member** (minimising mean mash
  distance to co-members over the condensed float32 triangle; NaN — no
  shared hashes — counts as 1.0; ties → lexicographic prophage id);
- **Tier-A pool:** greedy farthest-point sampling over medoid-to-medoid
  distances seeded with the tier-A candidate whose medoid has the
  lexicographically smallest clade id; stop at **18**;
- **Tier-B pool:** same procedure over tier-B candidates measured against
  **all** already-selected representatives; stop at **6** → **24 genomes**
  (mirrors NTM stage-1's 20; A:B = 3:1 as in NTM 15:5);
- pool shortfall → fill from the other tier; every deviation recorded in
  the eligibility ledger.

**Member representatives.** Per selected clade: the medoid member (bait
donor) and the member farthest from the medoid (recorded for provenance),
both requiring ≥ 600 bp. Deterministic tie-breaks by prophage id.

**Output:** `panel/pilot_selection.tsv` (24 rows, committed) + full
eligibility ledger (external `work/eligibility_ledger.tsv`) giving every
clade an explicit include/exclude reason.

## Stage 2 — Bait classes and coordinates

All baits are **500–2,500 bp**, 1-based inclusive coordinates on the named
source molecule, forward orientation (masking is strand-symmetric).

| Class | Source | Window | Count |
|---|---|---|---|
| `interior_module` | ML genome | 1,200 bp centred on a key-module CDS midpoint; priority terminase large → portal → major head/capsid → tail tape measure → endolysin → integrase (Pharokka annot regex, same module definitions as NTM) | **1 / clade** |
| `interior_generic` | ML genome | 1,200 bp at genome midpoint; only when no module gene is annotated | fill to 1 / clade |
| `junction` | ML genome | 1,000 bp centred on the exact breakpoint between consecutive stitched partitions (500 bp each side; clamped/shifted inward at genome ends) | ≤ 3 / clade |
| `member_interior` | observed member prophage (medoid representative) | 1,200 bp centred at prophage midpoint (width = min(1200, len); skipped if < 500) | 1 / clade |
| `control_public_reference` | λ `NC_001416.1`, HK97 `NC_002167.1` | 1,200 bp mid-genome slice (centered window at `L//2+1`); exact accessions, slice coordinates and reference sha256 recorded | **2 total** |
| `control_shuffled_negative` | derived | fixed-seed (42) first-order Markov sequence fitted to the first accepted `interior_module` baits of 2 distinct genomes; length + mono-/dinucleotide composition matched; ≤ 5 canonical 31-mers shared with the whole panel (deterministic redraws ≤ 10); composition deviation reported | **2 total** |
| `control_host_negative` | source host assembly | 1,000 bp centred in the largest prophage-free stretch of any contig (deterministic tie-breaks); one from **K-12 MG1655 (`GCF_000005845.2`)** and one from a selected clade's source host; expected host-like (filter validation) | **2 total** |

**Junction selection (E. coli-specific).** Candidate joins are consecutive
partition pairs in the regenerated traversal genome; both partitions must
have sequence. **Only purely-reconstructed joins — `observed_adjacency_count
== 0` (no member in which the two partitions are exactly contiguous) — are
eligible**; this is the primary scientific target (the reconstruction
hypothesis the screen tests). Ranking: `co_occurrence_count` descending,
then `min(occurrence_a, occurrence_b)` descending, then lower `partition_a`.
Ties beyond that are impossible (partition ids unique per clade). Every
junction bait records both partition ids, exact breakpoint, each partition's
occurrence count and ML-genome coordinates, `co_occurrence_count`,
`observed_adjacency_count`, and `adjacency_observed=no` (by construction).
Clades with fewer than 3 purely-reconstructed joins simply contribute fewer
junction baits — no fallback to observed joins, no approximation.

## Stage 3 — 31-mer masking and bait QC gates

Identical to NTM v1 (`baitlib.py`, k = 31, exact 62-bit 2-bit encoding,
canonical = min(code, revcomp code)):

| Mask | Rule |
|---|---|
| `ambiguous` | 31-mer contains any non-ACGT base |
| `low_complexity` | ≥ 7-bp homopolymer overlap or base entropy < 1.2 bits |
| `host_like` | canonical 31-mer present in the host panel |
| `internal_dup` | canonical 31-mer occurs > 1× within the bait |

**Host panel.** **K-12 MG1655 (`GCF_000005845.2`) always first**, then ≤ 4
source accessions per selected clade (member-representative sources first,
then sorted member sources), ≤ 60 total, deterministic order (selection
rank, then accession). Host sequences are scanned **with every called
prophage interval on that accession masked out** (from `26k_prophage1.csv`),
so prophage-derived 31-mers are not spuriously host-like. Panel composition
+ per-file SHA-256 recorded in `provenance.json`.

**Region-level rejection** (NTM rule): a window is *dominated* by
non-phage markers when ≥ 50% of its CDS-covered bases (≥ 30% CDS coverage
overall) match the reject-marker regex
`transposase|insertion|plasmid|resistance|ribosomal|rRNA|tRNA`, carry a
CARD/VFDB hit, or are PHROG category `DNA, RNA and nucleotide metabolism`
→ reject with `dominated_nonphage_markers`.

**Acceptance gates** (fixed order; first failure is the recorded reason):
length ∈ [500, 2500]; not marker-dominated; ambiguous_frac ≤ 0.10;
host_like_frac ≤ 0.10; low_complexity_frac ≤ 0.30; usable (unmasked,
unique) 31-mers ≥ 200; exact-sequence duplicate → `duplicate_sequence`;
> 50% canonical 31-mers shared with an already-accepted bait →
`cross_bait_kmer_dup` (higher-priority baits kept). Host-derived negative
controls are exempt from the host gate by design
(`expected_behavior=host_like`) with their actual host-like fraction
reported. Public-reference positive controls pass the same gates (they are
phage DNA; host-like fraction expected ≈ 0).

**Rejection ledger:** every rejected candidate (and each masked-out module
window) recorded with id, class, coordinates, failure reason and 31-mer
statistics — no silent drops.

## Stage 4 — Outputs & provenance

Committed in-repo (`research/ecoli_bait/panel/`, byte-deterministic):

- `baits.fa` — `>ECBAIT_<class>_<nn>` (+ clade-tagged for per-clade baits);
- `pilot_selection.tsv` — 24 selected genomes (tier, medoid, reps);
- `bait_manifest.tsv` — full per-bait provenance (id, class, clade/genome,
  contig, start/end/length, module gene or junction fields incl.
  `observed_adjacency_count`, sequence sha256, all 31-mer statistics, GC,
  acceptance status, expected behavior);
- `rejection_ledger.tsv`, `filter_stats.tsv`, `provenance.json` (input
  SHA-256s incl. staged traversal inputs + regenerated traversal checksums
  + per-clade verification results + public-ref checksums; parameters;
  deterministic `run_id` = SHA-256 over inputs+parameters — no wall-clock
  in committed outputs), `README.md`, `public_refs/` (the two reference
  FASTAs, < 100 KB total).

External (`…/ecoli_bait/`): staged clade partitions, regenerated traversal
JSONs + ML FASTAs + verification log, eligibility ledger, member-prophage
staging, host-panel staging + scan log, run log.

## Validation contract (task `## Validation` mapping)

- Every bait 500–2,500 bp, unique ID + sha256, source coordinates, class
  and 31-mer stats — enforced by a panel self-audit that **fails the run**
  on any violation (including: junction rows span their breakpoint,
  `adjacency_observed == no`, positive-control slices byte-match their
  reference).
- Junction baits carry `observed_adjacency_count == 0` provenance +
  co-occurrence counts; regeneration verified byte-identical per clade;
  no approximated junctions (missing/mismatch → zero junction baits,
  loudly recorded).
- Positive controls are exact mid-genome 1,200 bp slices of λ
  `NC_001416.1` and HK97 `NC_002167.1` with reference sha256s; shuffled
  negatives composition-matched (mono/dinucleotide L1 reported) and gated
  on ≤ 5 shared canonical 31-mers with the panel.
- Host-masked: no bait exceeds 10% host-like 31-mer fraction (gate), host
  panel checksummed; residual risk documented (bounded 60-genome panel,
  not all 26k hosts).
- Deterministic rerun byte-identical (`--verify-determinism`); pytest
  passes (`scripts/test_ecoli_bait_design.py` + the shared
  `ntm/v2/bait/scripts/test_bait_design.py`).
- Design preregistered (this file committed) **before** the panel is
  generated and before any screening.

## Explicit caveats

1. Host masking covers a bounded, checksummed ≤ 61-genome panel (K-12 +
   ≤ 60 diverse source hosts) outside called prophages — not all 26,074 E.
   coli hosts. Residual host-like 31-mers from unsampled hosts remain
   possible; Logan hits are confirmed by alignment downstream.
2. Junction provenance is **regenerated, not archived**: byte-identity of
   each regenerated ML genome to the frozen release is the verification;
   a mismatch removes junctions for that clade rather than weakening the
   rule. Regeneration depends on the untracked main-tree partition data,
   which is staged + checksummed before use.
3. IS/plasmid/AMR/rRNA filtering is annotation-driven (Pharokka/PHROG +
   CARD/VFDB); unannotated mobile elements rely on the low-complexity and
   uniqueness gates as second line of defense.
4. Module baits are *centred on* a module gene, not amplicons of the gene
   alone; module regexes mirror the functional-QC definitions.
5. Tier B admits `host_genes`-flagged clades deliberately (274/585 carry
   the flag; excluding them would bias the panel toward laboratory-typical
   phages); their baits still pass the same marker-domination and
   host-like gates.
6. The λ/HK97 positive controls validate the screen's recovery of
   *known* temperate phages; they do not calibrate detection of the
   reconstructed *novel* clades — that is the junction/module bait
   classes' job.
