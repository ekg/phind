# Non-tuberculous *Mycobacterium* (NTM) prophage-derived ML phage genomes

**913 maximum-likelihood (ML) phage genomes** (+ 472 ancestral-state genomes)
mined from prophages found in **7,303 NTM genome assemblies** (10,438 prophage
elements), spanning **342 NTM host clades** (species/lineages).

Companion files:
- `phage_ml_genomes.fa.gz` — one ML phage genome per prophage tight clade (913 entries, **bgzipped + faidx-indexed**)
- `phage_ancestral_genomes.fa.gz` — ancestral-state alternative per partitioned clade (472, bgzipped)
- `phage_ml_genomes.tsv` — per-genome metadata
- `host_range.tsv` — prophage clade → NTM host clade ids + species it spans

> Files are block-gzipped (`.gz`). Decompress with `bgzip -d` / `gunzip`, or use
> `samtools faidx phage_ml_genomes.fa.gz <name>` for random access (`.fai`/`.gzi`
> indexes are included for the ML set).

## What these genomes are

Each entry is a **maximum-likelihood (typical) genome** reconstructed from a
coherent clade of related prophages:

1. **NTM cohort** = family *Mycobacteriaceae* minus the TB complex & *M. leprae*
   (covers reclassified genera: *Mycobacteroides* abscessus, *Mycolicibacterium*
   fortuitum/smegmatis, …). 7,303 assemblies downloaded as PanSN bgzip FASTA.
2. **Prophage calling** with **geNomad v1.12** (`end-to-end`, find-proviruses,
   neural-net classifier) → 10,438 prophages.
3. **MASH distance triangle** over all prophages; **tight clades** (913) so every
   member is within ~75%+ ANI of its representative (leader clustering,
   threshold 0.25, clade size ≤ 100). 441 singletons = divergent prophages with
   no close relatives.
4. Within each clade, **all-wave alignment** (`allwave`, biWFA, k-nearest only),
   then **`impg partition`** into small homologous blocks (≤ 1 kb; median 500 bp).
5. A **weighted traversal** samples the most common path through each clade's
   partitions (probability ∝ partition occurrence); the **ML genome** is the
   majority-rule consensus of that typical path. `--mode ancestral` gives an
   NJ + Fitch-parsimony ancestral-state alternative.

Headers encode provenance:
```
>ntm_<clade>_ML status=ml|singleton n_members=N length=L host_clades=HC1,HC2,... species=Sp1,Sp2,...
```

## TSV columns (`phage_ml_genomes.tsv`)

| column | meaning |
|---|---|
| `clade_id` | prophage tight-clade id (`0_NNNN`) |
| `n_members` | prophages in the clade (1 for singletons) |
| `status` | `ml` (reconstructed) or `singleton` (observed sequence) |
| `ancestral` | ancestral-state genome produced? (`yes`/`n/a`) |
| `length` | genome length (bp) |
| `n_source_genomes` | distinct NTM genomes carrying this clade |
| `host_clades` | NTM host-clade ids spanned (see `host_range.tsv`) |
| `species` | NCBI species of source genomes |

## QC summary (913 ML genomes)

| length | count |
|---|--:|
| < 10 kb | 73 |
| 10–50 kb | 483 |
| 50–100 kb | 177 |
| 100–150 kb | 122 |
| ≥ 150 kb | 58 |

Median 33.5 kb, mean 53.8 kb. **GC median 64.0%** (matches mycobacterial hosts
~64%). Median max-N-run 0 (10 genomes with an N-run ≥ 100, max 411). All
non-singleton clades: internal median pairwise MASH ≤ 0.25.

## ML genomes per NTM clade (top host clades)

| host clade | species | ML genomes |
|---|---|--:|
| 0001 | *Mycobacteroides abscessus* | 212 |
| 0006 | *Mycolicibacterium fortuitum* | 54 |
| 0003 | *Mycobacterium avium* | 45 |
| 0005 | *Mycobacterium intracellulare* | 38 |
| 0013 | *Mycolicibacterium senegalense* | 29 |
| 0019 | *Mycobacteroides immunogenum* | 21 |
| 0002 | *Mycobacterium marinum* | 16 |

## Caveats

1. **Bimodal prophage yield.** Rapidly-growing species yield E. coli-like counts
   (abscessus 2.8/genome, chelonae 2.3, fortuitum 2.2); slow-growers lower.
   ***M. ulcerans* returned ~0 prophages across 1,040 genomes** — consistent with
   its known genome reduction, but a geNomad marker-coverage blind spot cannot
   be ruled out; a PhiSpy/Phigaro cross-check is recommended before concluding
   ulcerans is phage-free.
2. **Sparse pairwise evidence in large clades.** k-nearest alignment aligns a
   minority of possible pairs in capped 100-member clades; ML genomes from these
   are statistically grounded (weighted sampling) but should be treated as
   candidates for the largest constructs.
3. **Rare modules are sampled, not excluded.** Each ML genome is one draw
   (seed 42); minor variants arise across seeds.
4. **Mash trees are approximate** (sketch-based), not core-genome alignments.

## Reproducibility

End-to-end scripts live in `ntm/scripts/` (`acquire_ntm_accessions.py`,
`call_ntm_prophages.py`, `extract_full_ntm_prophages.py`, `full_ntm_mash_clades.py`,
`build_ntm_release.py`) reusing the shared `scripts/` pipeline
(`build_tight_clades.py`, `per_clade_alignment_pipeline.py`, `segment_paf.py`,
`traverse_partitions.py`). Design & status: `ntm/README.md`, `ntm/RELEASE.md`.
