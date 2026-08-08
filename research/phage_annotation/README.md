# Phage genome annotation + QC report (all genomes)

Per-genome gene annotation and QC for every phage genome generated —
**E. coli ML (585), NTM ML (913), NTM ancestral (472) = 1,970 genomes**.

## How made (≤64 threads)

1. **Pharokka** (`pharokka run -m --mmseqs2_only --skip_extra_annotations
   --skip_mash -g prodigal-gv`, PHROG DB) → gene calls (prodigal-gv) + PHROG
   function/category per CDS + VFDB/CARD. (Fast mmseqs2 path; the structural
   `phold`/ProstT5 path was tried but is far too slow on CPU for this scale.)
2. **CheckV** `end_to_end` (CheckV DB v1.5) → completeness (AAI/HMM-based),
   contamination, viral/host gene counts, quality category.
3. **build_annotation_report.py** merges into a per-genome table.

Driver: `ntm/scripts/build_annotation_report.py`. Full Pharokka/CheckV outputs
live on nvme (`annotation/pharokka_out`, `annotation/checkv_out`); only the
report TSVs are committed.

## Files

- `per_genome_annotation_qc.tsv` — one row per genome (1,970)
- `summary_by_source.tsv` — aggregates per source (ecoli_ml / ntm_ml / ntm_anc)

### `per_genome_annotation_qc.tsv` columns
genome_id, source, clade_id, length, gc_perc, cds_density, gene_count,
truncated_genes (prodigal partial ≠ 00), **terminase/portal/capsid/tail/
integrase/lysis** (0/1, detected from PHROG `annot` + `category`),
n_key_proteins, has_head_packaging/has_tail_cat/has_lysis_cat/has_integration_cat,
hypothetical, hypothetical_frac, checkv_miuvig (quality), completeness_pct,
contamination, viral_genes, host_genes, **flag**.

`flag` ∈ ok | missing_core_structural (no terminase+portal+capsid) |
low_completeness (CheckV Low-quality/Not-determined) | contamination |
host_genes | many_truncated (>20% genes partial).

## Summary (by source)

| source | n | median genes | %terminase | %capsid | %missing core | %low-compl | %contam | med hypoth |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| ecoli_ml | 585 | 83 | 62 | 57 | 31 | 28 | 4.1 | 0.57 |
| ntm_ml | 913 | 55 | 42 | 35 | 42 | 52 | 5.6 | 0.82 |
| ntm_anc | 472 | 129 | 42 | 36 | 42 | 33 | 1.7 | 0.86 |

CheckV quality (all): Complete 1, High-quality 753, Medium-quality 379,
Low-quality 792, Not-determined 45.

## Reading the numbers

- **Core machinery is detected in the large clades.** The 100-member clades
  (e.g. `ntm_0_0000_ML`, `clade_0_0000_ML`) carry terminase + portal + capsid +
  tail + lysis (5–6/6 key proteins). `missing_core_structural` concentrates in
  short/singleton fragments and the most divergent lineages.
- **NTM genomes are more fragmented** than E. coli (52% vs 28% Low-quality) and
  more hypothetical (~82% vs ~57%) — expected for (a) reconstructed ML genomes
  from diverse, ORFan-rich *Mycobacteroides/Mycolicibacterium* phages and (b) the
  shorter singleton fragments included as their own genomes.
- **host_genes (1,110)** reflects these being **prophage-derived** — flanking
  host genes at integration boundaries are expected, not a pipeline defect.

## Caveats

1. Annotation is **mmseqs2-only** (sequence-based); remote/novel ORFans stay
   "hypothetical" — the high hypothetical rate for NTM is partly real
   (mycobacteriophage ORFans) and partly annotation depth. A structural pass
   (phold, on GPU) would reduce it.
2. CheckV completeness is **AAI/reference-based** — it under-calls novelty
   (45 Not-determined) and is conservative for reconstructed genomes.
3. `missing_core` is a **detection** flag, not proof of absence — a genome may
   have an intact capsid that didn't hit a PHROG profile.
