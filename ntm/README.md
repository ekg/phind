# PHIND-NTM — ancestral and ML phage genomes from non-tuberculous *Mycobacterium* (NTM) prophages

**Objective.** Extend the PHIND pipeline from *E. coli* to **NTM** (non-tuberculous
mycobacteria): download the NTM assembly universe, call prophages, cluster them,
build per-clade prophage pangenomes with `allwave` + `impg partition`, and run the
two-level ML / ancestral traversal to emit **maximum-likelihood phage genomes per
NTM clade**. This is the same proven pipeline as the *E. coli* release
(`research/ml_phage_genomes/`), re-pointed at a multi-species host cohort.

This README is the design + plan of record. The *E. coli* machinery
(`scripts/build_tight_clades.py`, `per_clade_alignment_pipeline.py`,
`traverse_partitions.py`, `qc_ml_phage_genomes.py`) is reused ~90% as-is; only
the NTM-specific front-end (accessions, download path, prophage caller, host
clades) is new.

---

## 1. What counts as NTM (correct taxonomic scope)

NTM = family **Mycobacteriaceae** minus the *M. tuberculosis* complex (MTC) and
*M. leprae* / *M. lepromatosis*. The family scope is required because the 2018
Gupta reclassification split the old *Mycobacterium* into several genera still
used by NCBI:

| New genus (NCBI) | Notable NTM species |
|---|---|
| *Mycobacterium* (s.s.) | *M. avium, M. intracellulare, M. kansasii, M. marinum, M. ulcerans, M. gordonae, M. malmoense, M. xenopi, M. simiae, M. haemophilum* |
| *Mycobacteroides* | *M. abscessus* (3 subspp: abscessus/massiliense/bolletii), *M. chelonae, M. immunogenum* |
| *Mycolicibacterium* | *M. fortuitum, M. smegmatis, M. goodii, M. poriferae, M. senegalense, M. neoaurum* |
| *Mycolicibacter*, *Mycolicibacillus* | *M. kubicae*, etc. |

Querying only `"Mycobacterium"` **misses ~half** the NTM (abscessus, chelonae,
fortuitum, smegmatis, …). Scope: `"Mycobacteriaceae"[Organism] NOT MTC NOT
leprae NOT lepromatosis`.

## 2. Scoping result (verified 2026-08-06)

`acquire_ntm_accessions.py` (E-utilities) → `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/accessions/`:

- **7,352 NTM assemblies** (after organism-name exclusion + RefSeq/GenBank dedup), **5,331 RefSeq + 2,021 GenBank**.
- **38.97 Gb total sequence**, median genome **5.39 Mb**.
- Assembly levels: 550 Complete, 69 Chromosome, 2,785 Scaffold, 3,948 Contig.
- **735 species**; dominant (clade-defining) ones:

| Species (NCBI genus) | assemblies | ~genome |
|---|---:|---:|
| *Mycobacteroides abscessus* | 2,634 | 5.1 Mb |
| *Mycobacterium ulcerans* | 1,055 | 5.3 Mb |
| *Mycobacterium sp.* (unclassified) | 677 | 3.9 Mb |
| *Mycobacterium avium* | 363 | 5.3 Mb |
| *Mycobacterium kansasii* | 209 | 6.5 Mb |
| *Mycobacterium intracellulare* | 160 | 5.6 Mb |
| *Mycobacterium marinum* | 103 | 6.2 Mb |
| *Mycobacteroides chelonae* | 91 | 5.1 Mb |
| *Mycolicibacterium fortuitum* | 89 | 6.6 Mb |
| *Mycolicibacterium smegmatis* | 72 | 6.8 Mb |
| … (~25 more with ≥10 assemblies) | | |

Abscessus alone splits into its 3 clinical subspecies (abscessus 1,122 /
massiliense ~450 / bolletii ~119), which are themselves natural clades.

## 3. Pipeline at a glance

```
NTM accessions (E-utilities, done)                 ntm/v1/accessions/
  └─ download PanSN bgzip (acquire_remaining.py)   ntm/v1/genomes/   [TO RUN]
       └─ call prophages (geNomad/Phigaro)         ntm/v1/prophages/ [TO RUN]
            └─ extract prophage FASTA              ntm/v1/full_prophages.fa
                 ├─ MASH sketch+triangle+tree       ntm/v1/mash_tree/  (prophages)
                 │     └─ tight prophage clades     ntm/v1/clades/     (build_tight_clades)
                 │          └─ per-clade allwave + segment + impg partition
                 │               └─ ML + ancestral traversal (traverse_partitions.py)
                 │                    └─ prophage-clade ML phage genomes
                 └─ host-genome MASH dist          ntm/v1/host_clades/ (7.3k genomes)
       JOIN prophage-clade genome × host-clade distribution
            └─ ML phage genomes PER NTM CLADE  +  QC + release   ntm/v1/release/
```

## 4. Two-level clade scheme (the key design choice)

The *E. coli* run was single-species, so "clade" meant a prophage-content
community. NTM spans hundreds of species, so we use **two nested levels**, both
driven by MASH distance, and both reusing existing code:

1. **Host clades** (the "NTM clade" reporting unit). Whole-genome `mash dist`
   over the 7.3k assemblies → species/lineage clades (e.g. *M. abscessus*
   subsp. *abscessus*, *M. avium*, *M. ulcerans*, …). A phage genome is reported
   "per NTM clade" via its host-clade membership distribution.
2. **Prophage tight clades** (the alignment unit). `mash dist` over all NTM
   prophages → tight prophage clades via `build_tight_clades.py` (every member
   within `--threshold` MASH of its representative, `--max-size` capped) so
   `allwave` only aligns genuinely related phages. Each tight clade → one ML
   phage genome.

**Join.** Each prophage-clade ML genome is annotated with the host clade(s) its
members come from → delivers "phage ML genomes **per NTM clade**", while
preserving the host-range signal (a prophage clade spanning several NTM species
is itself a finding). This mirrors how the *E. coli* genomes were labeled by
source assembly.

## 5. Prophage calling (the one new dependency)

*E. coli* came with a pre-called prophage CSV; for NTM we call them ourselves.
`micromamba` is available, so a caller can be installed via bioconda.
Candidates (decision in `ntm-prophage-caller-capability`):

- **geNomad** (recommended): marker-gene + BATH, robust on novel/high-GC genomes,
  actively maintained, gives coordinates + taxonomy + plasmid separation.
- **Phigaro v2.4** (project-certified for *E. coli*, `workflow/phigaro_v2_4_certification/`):
  Prodigal CDS + PHROG HMMs; consistent with the *E. coli* calls; needs
  hmmer + PHROG db. NTM genomes are ~64% GC and phage DNA often differs in GC,
  which is actually a useful signal for these callers.

Mycobacteria are high-GC (~64%); both callers handle this. The chosen caller is
smoke-tested on a handful of NTM genomes before the full run.

## 6. Tools (unchanged from the *E. coli* pipeline)

| Tool | Status | Role |
|---|---|---|
| `mash` 2.3 | installed | sketch, dist triangle, host + prophage trees |
| `allwave` 0.1.0 | `~/.cargo/bin/allwave` | per-clade all-vs-all biWFA alignment (PAF) |
| `impg` 0.4.1 | `~/.cargo/bin/impg` | `partition` (small ≤1 kb partitions) |
| `samtools`/`bgzip` 1.19 | installed | PanSN FASTA, faidx extraction |
| `micromamba` | installed | prophage-caller environment |

Reused scripts (parametrized for NTM paths): `scripts/build_tight_clades.py`,
`scripts/extract_clade_fasta.py`, `scripts/per_clade_alignment_pipeline.py`,
`scripts/segment_paf.py`, `scripts/traverse_partitions.py`,
`scripts/qc_ml_phage_genomes.py`. New: `ntm/scripts/acquire_ntm_accessions.py`
(done), `ntm/scripts/download_ntm_genomes.py` (thin wrapper on acquire logic),
`ntm/scripts/call_ntm_prophages.py`, `ntm/scripts/host_clades_mash.py`.

## 7. Status

| Step | Status | Output |
|---|---|---|
| NTM accession scoping (7,352) | ✅ done | `ntm/v1/accessions/` |
| Prophage-caller capability pick | ⏳ next | decision doc |
| Download 7,352 genomes → PanSN bgzip | ⏳ to run | `ntm/v1/genomes/` |
| Call prophages on all genomes | ⏳ to run | `ntm/v1/prophages/` |
| Extract full prophage FASTA | ⏳ to run | `ntm/v1/full_prophages.fa` |
| Prophage MASH sketch+triangle+tree | ⏳ to run | `ntm/v1/mash_tree/` |
| Host-genome MASH dist → NTM host clades | ✅ done | `ntm/v1/host_clades/` |
| Tight prophage clades (`build_tight_clades`) | ✅ done | `ntm/v1/clades/` (913 clades, `ntm/v1/mash_clades/clades/`) |
| Per-clade allwave + partition | ✅ done | `ntm/v1/clades/<c>/` (913 manifests, `pipeline_results.json`) |
| ML + ancestral traversal | ⏳ to run | prophage-clade genomes |
| Join by host clade → genomes per NTM clade | ⏳ to run | `ntm/v1/release/` |
| QC + release manifest | ⏳ to run | `ntm/v1/release/RELEASE.md` |

## 8. Working rules

- NTM data lives under `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/`; the repo
  holds only code + small manifests (oversized FASTAs/triangles stay git-ignored).
- Do not touch the *E. coli* inputs (`26k_*`, `prophage_homology_survey/full_prophages.fa`).
- `ntm/v1/accessions/ntm_accessions.txt` is the immutable accession manifest
  (checksum before analysis); the manifest TSV carries full metadata.
- Same allwave rule as *E. coli*: sparsified k-nearest (`-p tree:5:0:0.0` /
  `tree:10:0:0.0`), k-farthest = 0 (no stranger-joining), scores `0,5,8,2,24,1`.
