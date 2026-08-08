# PHIND-NTM — ML & ancestral phage genomes from non-tuberculous *Mycobacterium*

**Deliverable.** Maximum-likelihood (and ancestral) phage genomes mined from the
prophages of the NTM assembly universe, **one ML genome per prophage tight clade,
labeled by NTM host clade**.

**Date.** 2026-08-08

## Cohort & pipeline (end-to-end on NTM)

1. **Accessions** — family *Mycobacteriaceae* minus TB-complex/leprae → **7,352**
   NTM assemblies (covers reclassified genera: *Mycobacteroides* abscessus,
   *Mycolicibacterium* fortuitum/smegmatis, …). `ntm/scripts/acquire_ntm_accessions.py`.
2. **Download** — 7,303 valid PanSN bgzip genomes (49 NCBI-suppressed).
3. **Prophage calling** — geNomad v1.12 `end-to-end` (find-proviruses, NN on,
   `PYTHONNOUSERSITE=1`) → **10,438 prophages** (`ntm_prophages.csv`).
4. **Extract** → `full_prophages.fa` (10,438 seqs, median 34 kb).
5. **MASH sketch + triangle + UPGMA** over all prophages.
6. **Tight clades** (`scripts/build_tight_clades.py`, thr 0.25 / max 100) → **913
   clades** (472 alignable ≥2 members + 441 singletons; internal median mash 0.074).
7. **Per-clade allwave + impg partition** (`scripts/per_clade_alignment_pipeline.py`,
   k-nearest `tree:5:0:0.0`, window 500) → ≤1 kb partitions (median 500 bp).
8. **Two-level ML + ancestral traversal** (`scripts/traverse_partitions.py`).

## Files (`ntm/v1/release/`)

| File | Contents |
|---|---|
| `all_ntm_ml_phage_genomes.fa` | **913 ML phage genomes** (472 ML-reconstructed + 441 singletons) |
| `all_ntm_ancestral_phage_genomes.fa` | 472 ancestral-state genomes |
| `release_manifest.tsv` | per-clade: members, status, length, source genomes, host clades, species |
| `host_range.tsv` | prophage clade → host clade ids + species it spans |
| `per_ntm_clade_summary.tsv` | host clade → # ML genomes |

Headers: `>ntm_<cid>_ML status=ml|singleton n_members=N length=L host_clades=… species=…`

## QC (913 ML genomes, ~49 MB)

- Lengths: min 4.3 kb, **median 33.5 kb**, mean 53.8 kb, max 153 kb. Buckets:
  <10kb 73 · 10–50kb 483 · 50–100kb 177 · 100–150kb 122 · ≥150kb 58.
- **GC median 64.0%** (matches mycobacterial hosts ~64%); median max-N-run 0
  (10 genomes with an N-run ≥100, max 411).
- **342 NTM host clades** represented.

## ML genomes per NTM clade (top host clades)

| host clade | species | ML genomes |
|---|---|--:|
| 0001 | *Mycobacteroides abscessus* | 212 |
| 0006 | *Mycolicibacterium fortuitum* | 54 |
| 0003 | *Mycobacterium avium* | 45 |
| 0005 | *Mycobacterium intracellulare* | 38 |
| 0013 | *Mycolicibacterium senegalense* | 29 |
| 0015 | *M. arcueilense* | 22 |
| 0019 | *Mycobacteroides immunogenum* | 21 |
| 0002 | *Mycobacterium marinum* | 16 |

## Caveats

1. **Bimodal prophage yield.** Rapidly-growing species yield E. coli-like counts
   (abscessus 2.8/genome, chelonae 2.3, fortuitum 2.2); slow-growers lower
   (avium/intracellulare/marinum ~0.6–1.0). ***M. ulcerans* returned ~0/1,040
   genomes** — consistent with its known genome reduction, but a geNomad
   marker-coverage blind spot cannot be ruled out; a PhiSpy/Phigaro cross-check
   is recommended before concluding ulcerans is phage-free.
2. **k-nearest sparsification** in large clades → sparse pairwise evidence; ML
   genomes from capped 100-member clades are statistically grounded but should
   be treated as candidates for the largest constructs.
3. Rare partitions are *sampled* (not excluded); each ML genome is one draw
   (seed 42). `--mode ancestral` gives the ancestral-state alternative per clade.
4. Mash trees are approximate (sketch-based), not core-genome alignments.
