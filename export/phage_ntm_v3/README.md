# Mycobacteriaceae prophage-derived ML phage genomes — NTM v3 (unified)

**1,304 maximum-likelihood (ML) phage genomes** (+ **893** ancestral-state
genomes) mined from **36,857 prophages** across the **union cohort** of
**34,163 Mycobacteriaceae genomes** — spanning **443 host clades**, of which
**292** are represented in the ML set.

This is the v3 all-inclusive expansion: the 2026-09-22 collaborator delivery
(26,499 genomes: 17,920 SRA-assembled runs + 8,331 NCBI + 248 BV-BRC) unioned
with v2 local holdings. Unlike the v1 NTM export, the cohort is **not restricted
to NTM**: MTC (TB complex) genomes are included and **labelled**, not filtered.

Companion files:
- `phage_ml_genomes.fa.gz` — one ML phage genome per prophage tight clade (1,304 entries, **bgzipped + faidx-indexed**)
- `phage_ancestral_genomes.fa.gz` — ancestral-state alternative per partitioned clade (893, bgzipped)
- `phage_ml_genomes.tsv` — per-clade metadata (superset of the v1 columns; adds `host_scope`, `host_clade_ids`, `n_host_clades`, `source_breakdown`)
- `per_scope_summary.tsv` — NTM / MTC / MIXED roll-up
- `cross_boundary_report.tsv` — prophage clades whose hosts span NTM↔MTC
- `host_clade_scope.tsv` — audit of the scope derivation per host clade
- `release_scope_statement.md` — the scoping statement

> Files are block-gzipped (`.gz`). Decompress with `bgzip -d` / `gunzip`, or use
> `samtools faidx phage_ml_genomes.fa.gz <name>` for random access (`.fai`/`.gzi`
> indexes are included for the ML set).

## What these genomes are

Each entry is a **maximum-likelihood (typical) genome** reconstructed from a
coherent clade of related prophages:

1. **Cohort** — union of the 26,499-genome delivery and v2 holdings → 34,163
   canonical PanSN bgzip objects (26,890 NTM + 7,273 MTC). Prophage coordinates
   are the collaborator's BV-BRC phigaro QC-passed calls (47,994 rows).
2. **Prophage extraction** — unified manifest, **36,857 prophages** extracted, 0 errors.
3. **Host clades** — whole-genome MASH over all 34,163 genomes → **443 host clades**.
4. **Prophage tight clades** — MASH triangle over all prophages; leader
   clustering (threshold 0.25, max clade size 100, community 0) → **1,304 clades**
   (893 alignable + 411 singletons). Member partition is exact; median internal
   MASH ≤ 0.25 with 0 violations.
5. **Per-clade alignment** — `allwave` (biWFA, sparsified k-nearest
   `tree:k:0:0`, no stranger-joining), then `impg partition` (window 500) →
   **726,121 partitions / 1,847,004 intervals**.
6. **ML + ancestral traversal** — weighted path sampling over each clade's
   partition graph; ML = majority-rule consensus of the typical path;
   `--mode ancestral` = NJ + Fitch parsimony. `--n-samples 5 --seed 42`.
7. **Annotation** — Pharokka v1.10.1 (`-m --mmseqs2_only`) + CheckV v1.1.1:
   2,197 per-genome GFF3 (272,974 CDS). Index: `ntm/v3/annotation_gff_index.tsv`.

Headers encode provenance:
```
>ntm3_<clade>_ML status=ml|singleton n_members=N length=L host_clade_ids=HC1,HC2,... host_scope=NTM|MTC|MIXED species=Sp1,Sp2,...
```

## Scope — label, not filter

The cohort deliberately spans the whole family. `host_scope` is a **label**:

| `host_scope` (per prophage clade) | clades | ML genomes | members |
|---|---:|---:|---:|
| NTM | 1,285 | 1,285 | 36,228 |
| MTC | 12 | 12 | 202 |
| MIXED | 7 | 7 | 427 |
| **total** | **1,304** | **1,304** | **36,857** |

- Per **genome**, `host_scope` ∈ {NTM, MTC} (binary). A host clade is MTC when
  the majority of its genomes carry an MTC species prefix; members inherit it.
  Only `host_clade_0002` (7,273/7,467 = 97.4% MTC) qualifies.
- **7 prophage clades span the boundary**; two strongly mixed:
  `0_0306` (93 MTC / 7 NTM, median MASH 0.035) and `0_0338` (90/10, 0.045).
  See `cross_boundary_report.tsv` for per-clade counts + confidence.
- For a **strictly-NTM** analysis, take `host_scope == NTM` plus the NTM members
  of `MIXED` clades. Nothing needs re-running.

## QC summary (1,304 ML genomes)

| length | count |
|---|--:|
| < 10 kb | 381 |
| 10–50 kb | 452 |
| 50–100 kb | 96 |
| 100–150 kb | 292 |
| ≥ 150 kb | 83 (at the 150 kb budget cap) |

Median 24.0 kb, mean 57.2 kb. 411 singletons pass through as observed sequences.
All alignable clades: internal median pairwise MASH ≤ 0.25.

## ML genomes per host clade (top)

| host clade | ML genomes |
|---|---:|
| `host_clade_0001` (*Mycobacteroides abscessus*) | 435 |
| `host_clade_0008` (*Mycolicibacterium smegmatis*) | 163 |
| `host_clade_0004` (*Mycobacterium paraintracellulare*) | 134 |
| `host_clade_0003` (*M. avium* subsp. *paratuberculosis*) | 86 |
| `host_clade_0006` (*Mycobacterium kansasii*) | 66 |
| `host_clade_0009` (*M. ulcerans*) | 47 |

## Caveats

1. **"Tight" means the median**, not the max: 311 alignable clades contain at
   least one pair above 0.25 (some saturated at 1.0). Documented criterion.
2. **Cap-limited clustering**: 275 of 1,304 clades sit exactly at the
   `--max-size 100` cap, so the largest prophage families are split at 100.
3. **Rare modules are sampled, not excluded** — each ML genome is one draw
   (seed 42); minor variants arise across seeds.
4. **135 v2-only run accessions** have no delivered FASTA (36 objects, 83
   prophage rows) and are genuinely unavailable.
5. Host species labels for v2-only genomes come from the v2 QC list; the MTC
   block is v2-only and was externally spot-verified.
6. Ancestral FASTA headers carry `status=ml` (inherited); cosmetic.

## Reproducibility

Stage drivers: `ntm/v3/scripts/` (`ingest_delivery_raw.sh`,
`build_v3_pansn_objects.py`, `build_v3_prophage_manifest.py`,
`extract_v3_prophages.py`, `full_ntm_mash_clades_v3.py`,
`host_clades_mash_v3.py`, `build_v3_ml_genomes.py`, `run_annotation_v3.py`,
`build_v3_release.py`, `validate_v3_partitions.py`, `validate_v3_ml.py`),
reusing the shared `scripts/` pipeline. Design, status and receipts:
`ntm/v3/RELEASE.md`, `ntm/v3/PLAYBOOK.md`, `ntm/v3/README.md`.

Verified end-to-end by independent recomputation (`ntm/v3/s8_release_validation.md`):
counts, **100% sequence byte-identity vs source**, scope reconciliation
(391/12/7/1285/19), cross-boundary report and annotation index all reproduce.
