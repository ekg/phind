# NTM v3 — unified Mycobacteriaceae prophage phage-genome catalog

**Deliverable.** Maximum-likelihood (ML) and ancestral "typical" phage genomes,
one per prophage tight clade, derived from prophages called across the whole
**Mycobacteriaceae** family — non-tuberculous mycobacteria (NTM) *and* the
*M. tuberculosis* complex (MTC) together, as a single catalog with scope labels.

**Date:** 2026-10-06 · **Cohort:** union of the 2026-09-22 collaborator delivery
(26,499 genomes) and v2 local holdings.

---

## ⚠️ CORRECTION (2026-10-07) — the ML genomes are mosaics, not genomes

The v3 ML genomes below are **not single phage genomes**. The reconstruction
walk (`scripts/traverse_partitions.py`, default free mode) treats partitions as
independent tiles and only *biases* toward observed adjacency
(`w * (1 + beta*adj)`); it never *requires* it. It therefore accumulates
high-support blocks until the `--max-length 150000` budget stops it.

Measured consequences:

| | |
|---|---|
| clades whose ML genome is >1.5× the median member length | **747 / 1,304 (57%)**, up to **6.2×** |
| prophages covered by those mosaic clades | **35,747 / 36,857 (97%)** |
| CheckV "High-quality" genomes that are mosaics | **387 / 395 (98%)** |

**CheckV completeness rewards the mosaic**: concatenating every common block
includes every marker gene, so a 150 kb chimera scores 99–100% complete while a
real genome-length reconstruction scores ~81%. The quality label in this
release is therefore an artifact of the defect, not evidence against it.
Per-clade audit: `ntm/v3.1/v3_mosaic_audit.tsv` (length ratio, CheckV, flag).

**Do not use the v3 ML/ancestral FASTAs as genome sequences.** Use the v3.1
genome-path reconstructions (`--path-mode observed`, see `ntm/v3.1/SPEC.md` and
the A/B report), which walk real observed adjacencies to a natural terminus and
land at ~1.0× the member length. v3.1 replaces the **ML step only** — clades,
partitions, host clades and annotation are unchanged.

---

## Headline numbers

| stage | count |
|---|---:|
| cohort genomes (union) | 34,163 canonical objects |
| — NTM | 26,890 |
| — MTC | 7,273 |
| prophages extracted | 36,857 |
| prophage tight clades | 1,304 (893 alignable + 411 singletons) |
| partitions | 726,121 (1,847,004 intervals) |
| **ML phage genomes** | **1,304** |
| **ancestral phage genomes** | **893** |
| annotated (Pharokka GFF3) | 2,197 (1,304 ML + 893 ancestral), 272,974 CDS |
| cross-boundary prophage clades (NTM↔MTC hosts) | 7 (2 high-confidence) |

## End-to-end pipeline

1. **Ingest** — 26,499 delivered FASTAs (17,920 SRA-assembled runs + 8,331 NCBI
   + 248 BV-BRC) → 25,951 canonical objects; union with 8,212 v2-only objects
   = 34,163.
2. **Host clades** — MASH over 34,163 genomes → 443 host clades.
3. **Prophage extraction** — unified manifest, 36,857 prophages, 0 errors.
4. **Prophage MASH clades** — thr 0.25 / max-size 100 / community 0 → 1,304
   clades (member partition exact; median internal MASH ≤ 0.25, 0 violations).
5. **Per-clade alignment + partition** — allwave `tree:k:0:0` + impg window 500;
   893 alignable with full outputs, 0 failures.
6. **ML + ancestral traversal** — `traverse_partitions.py`, both modes,
   `--n-samples 5 --seed 42`; determinism re-run 10 clades × 2 modes, 0 mismatches.
7. **Annotation** — Pharokka v1.10.1 (`-m --mmseqs2_only`) + CheckV v1.1.1;
   per-genome GFF3 + committed index.
8. **Release** — host-clade join, scope labels, cross-boundary report.

## Files

`$NVME/ntm/v3/release/` (`$NVME = /mnt/nvme3n1/erikg/phind-genome-work`):

| file | sha256 |
|---|---|
| `all_ntm_v3_ml_phage_genomes.fa` | `7856859bf93524f19e04e5d0cb073f11f561ba278fcec18bc8c98bde6af6da16` |
| `all_ntm_v3_ancestral_phage_genomes.fa` | `97ec5e8750fc616ca0c2c6c5fc76e10e6bf1b79afdb00a7c156ee332872f8c48` |
| `release_manifest.tsv` | `e7a6fd028d107c0229a9d27d7e7e827c01cbe9efecc74e422d08919b611b82aa` |
| `per_scope_summary.tsv` | `03bb9b0b74fe86b6ec5b877a456b827b9b2a5e2b809201c116054e191010efa9` |
| `cross_boundary_report.tsv` | `6ff725f76a5a13b9683274d45d9d4b20fbbeb74c775407dd04e556bf168612f6` |
| `host_clade_scope.tsv` | `21d366b416465e29fd9a3aa8d650ed9cb76ec550a757de0eb6087a7adfc79e38` |
| `release_scope_statement.md` | `c5169656ce87e0dec637ceb1e781123f7bbd464079cc52499df9078a9018c234` |

Upstream: `full_prophages.fa` sha256
`a949c653b104521f9d52aa95af4a233a3e35418c10cf3c782e26a6e3317441bf` (787 MB, 36,857 records).

Annotation: `$NVME/ntm/v3/annotation/pharokka_out/per_genome_gff/` (2,197 GFF3);
committed index `ntm/v3/annotation_gff_index.tsv`.

## Scope — label, not filter

Per the binding scope decision (`ntm/v3/PLAYBOOK.md`):

- Per genome `host_scope` ∈ {NTM, MTC} — a **label**, nothing is dropped.
- Per prophage clade label ∈ {NTM, MTC, MIXED}.
- A host clade is MTC when the **majority** of its genomes carry an MTC species
  prefix; members inherit the clade's scope. Only `host_clade_0002` (7,273/7,467
  = 97.4% MTC) qualifies. The rule is auditable in `host_clade_scope.tsv`, and
  `release_scope_statement.md` records the 9 `Mycobacterium sp.` genomes that
  are clade-derived MTC.

### Subsetting for a strictly-NTM analysis

```
host_scope == NTM   → 1,285 clades, 36,228 prophage members
host_scope == MTC   →    12 clades,    202 prophage members
host_scope == MIXED →     7 clades,    427 prophage members (span both)
```

A strictly-NTM deliverable is `host_scope == NTM` **plus the NTM members of MIXED
clades** (`cross_boundary_report.tsv` gives per-clade NTM/MTC member counts).
Nothing needs re-running.

## Findings

- **MTC prophages form a distinct island** — 12 wholly-MTC clades.
- **3 prophage lineages cross the NTM↔MTC boundary at tight distance**, two
  strongly mixed: `0_0306` (93 MTC / 7 NTM, median MASH 0.035) and `0_0338`
  (90/10, 0.045). Only visible because the cohort is unified.
- Prophage clustering is **cap-limited**: 275 of 1,304 clades sit exactly at the
  `--max-size 100` cap, so the largest families are split at 100.
- "Tight" means **median** pairwise MASH ≤ 0.25; 311 alignable clades contain at
  least one pair above 0.25 (documented criterion, not a defect).

## Caveats

- 135 v2-only run accessions have no delivered FASTA (36 objects, 83 prophage
  rows) and remain genuinely unavailable.
- Host species labels for v2-only genomes come from the v2 QC list; the MTC block
  is v2-only and was externally spot-verified (e.g. `GCA_000154585.2` = *M.
  tuberculosis* KZN 4207).
- Ancestral FASTA headers carry `status=ml` (inherited); cosmetic.

## Verification

Every stage was independently recomputed by a separate read-only lane; the final
release passed a full verification (`ntm/v3/s8_release_validation.md`): counts,
**100% sequence byte-identity vs source**, scope reconciliation (391/12/7/1285/19),
cross-boundary report, and annotation index all reproduce exactly.
