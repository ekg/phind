# Integrated 250-Genome Pilot Report

**Release ID:** `integrated-pilot-250-v1-31cd0f3658837296`  
**Task ID:** `run-integrated-250-genome`  
**Date:** 2026-07-26  
**Verdict:** PASS  
**Scale-up authorization:** `GO_500`

---

## Executive Summary

This report documents the **scale-bearing** integrated end-to-end pilot on the frozen
250-assembly rung. It consumed only validated canonical objects
(`canonical-cohort-250-v1-a6184d7d6ee08bda`) and a PASS coordinate/source-semantics
policy (`prophage-semantics-v2-7dc695b85e5fd229`, EXTRACTION_GO/ALLOW), and reused the
immutable N=100 integrated pilot (`integrated-pilot-100-v1-0a11eda244a9def8`) read-only
as the prior rung of record.

Unlike the N=100 rung (non-scale-bearing), this rung applies the **scale-trend gate**:
the per-assembly build runtime scales with time exponent **1.016** (limit ≤1.3) and the
runtime per-base slope changes **+1.8%** (limit ≤25% unexplained). All applicable gates
passed, and the pairwise scaling authorizes the next rung.

**Automatic downstream scale-up to N=500: AUTHORIZED (`GO_500`).**

---

## 1. Inputs and Predecessor Verification

| Input | Release ID | SHA-256 (release.json) | Verdict |
|-------|------------|------------------------|---------|
| Canonical Cohort 250 | canonical-cohort-250-v1-a6184d7d6ee08bda | dcf2b887afa51e4e0e739ae2fef9b5a9d72fb8bc9a4d698a161a99673aaf504a | PASS |
| Prior Integrated Pilot 100 (reused) | integrated-pilot-100-v1-0a11eda244a9def8 | 6816c4e24f6511e45196d91112da96ab7f56082732c4c363b74ae7010a80e273 | PASS/ALLOW |
| Prophage Semantics v2 | prophage-semantics-v2-7dc695b85e5fd229 | 5d8403eb070d8a62140adfe7260b7fde6897598f72ac1c536879e78e8ea2b992 | EXTRACTION_GO/ALLOW |
| Root: 26k_ecoli_accession.txt | — | 1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5 | PASS (start & finish) |
| Root: 26k_prophage1.csv | — | 6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996 | PASS (start & finish) |
| Integrated Plan | — | fb58d25a6f4971137ab0dcb82dae09eac5d177e37ba34f857845ce2d7e0a6da8 | PASS |
| Prophage Distribution Audit | — | feb9b687fb0722a4f073f7105088f19a2dabebcbc0378dbe5c4799b0b7f29fdc | PASS |

**Global distinct assembly cap:** 1,000 (enforced)  
**Actual distinct assemblies in this rung:** 250  
**New assembly downloads:** 0 (all reused from the canonical cohort)

---

## 2. Reuse vs. Rebuild (release-scoped, checksum-proven)

The immutable N=100 integrated pilot is reused **read-only** as the prior rung of
record (release ID, methodology, coordinate policy, and contract structure). Only
products that legitimately depend on N were rebuilt at N=250:

| Decision | Object | Evidence |
|----------|--------|----------|
| REUSE | N=100 integrated release | COMPLETE + SHA256SUMS (18 rows) verified by digest |
| REUSE | N=100 cohort/bases | frozen canonical-cohort-100 (512,421,261 bases) |
| REUSE | coordinate policy / methodology | integrated plan + C1_RAW_1_BASED_CLOSED (unchanged) |
| REBUILD | assembly QC reconciliation | 250 assemblies via cross-predecessor resolution |
| REBUILD | host sketches / SYNG prefix / joins / extraction / query / clustering / matrix | scoped to the 250 cohort |
| REBUILD | scale-trend | N=100 prior rung vs N=250 current rung |

Gate: `prior_integrated_reuse` = PASS

---

## 3. Assembly QC Reconciliation

All 250 assemblies in the frozen cohort were reconciled via **cross-predecessor
canonical-object resolution** (`object_refs.tsv`): 10 objects reused from the N=10
release, 90 from the N=100 release, and 150 self-contained in the N=250 release. Every
assembly reached the `VALIDATED` terminal state with COMPLETE marker, manifest,
contigs, and the BGZF + `.fai` + `.gzi` index triple. No fuzzy/versionless joins or
guessed coordinates/strands.

| Metric | Value |
|--------|-------|
| Total assemblies | 250 |
| Distinct sequence-bearing assemblies | 250 |
| Cohort contigs | 41,050 |
| Cohort total bases | 1,276,442,466 |
| Terminal state: VALIDATED | 250 |
| Terminal state: QUARANTINED | 0 |
| BGZF / .fai / .gzi triples present | 250 / 250 / 250 |

Gate: `assembly_qc_reconciliation` = PASS

---

## 4. Host-Only Sketches (Engineering Validation)

Host-only Mash sketches/distances were computed for engineering validation of the
pipeline over 250 BGZF inputs. **No biological clade was defined or selected from
phage traits**; phage-positive engineering controls were held out from any clade
definition. This preserves the host/phage methodological separation.

Gate: `host_sketches_engineering` = PASS; `phage_blind_host` = PASS

---

## 5. Whole-Cohort SYNG Prefix (Six-File)

One staged whole-cohort SYNG prefix was built over the 250-cohort BGZF Pansn fastas.
The six inseparable files (`.meta .names .1khash .pstep .spos .1gbwt`) are staged
together; **partial and final prefixes never coexist** in a published location, and a
killed/partial build is never published.

| Component | Files |
|-----------|-------|
| Core | cohort.1khash, cohort.1gbwt |
| Sidecars | cohort.meta, cohort.names, cohort.pstep, cohort.spos |

Gate: `syng_prefix_integrity` = PASS

---

## 6. Lossless Prophage-Row Joins

All in-scope source rows were joined losslessly to the 250-assembly cohort by exact
accession, with the explicit coordinate policy **C1_RAW_1_BASED_CLOSED** selected and
C2_RAW_0_BASED_INCLUSIVE rejected. Every scoped row carries non-empty identity and
coordinate fields (genome, scaffold, begin, end). The three source scopes are preserved
separately.

| Scope | Rows in Cohort |
|-------|----------------|
| all_records | 1,249 |
| transposable_flag_positive | 68 |
| taxonomy_assigned | 1,096 |

Gate: `prophage_joins_lossless` = PASS

---

## 7. Bounded Extraction with Controls

Extraction used the resolved 1-based-closed convention (C1_RAW_1_BASED_CLOSED →
0-based half-open [begin-1, end)) with explicit contig-edge, circular-wrap, and
unknown-strand controls; exact-source digest verification and `samtools faidx`
round-trip enforced.

Gate: `extraction_controls` = PASS

---

## 8. IMPG Interval Query and Sequence Map

IMPG interval query (`impg query -b`) and sequence map (`impg map`) were exercised with
independent origin-recovery (≥95% interval overlap, 100% strand/spelling), coverage
positive controls (≥95% cover, ≥80% on origin), and negative controls (false positives
≤1%). Query partitions and hits are bounded by the scoped rows.

Gate: `impg_query_correctness` = PASS

---

## 9. Preliminary Clustering and Matrix

Clustering was performed at whole-element, protein/domain-family, and syntenic-module
levels, preserving unit type, copy number, callability, evidence IDs, and separate
source scopes. The long-form present/absent/uncallable/ambiguous matrix carries copy
counts and evidence IDs over 250 analysis units.

Gate: `clustering_preliminary` = PASS; `matrix_states` = PASS

---

## 10. Scale-Trend (the new scale-bearing gate)

This is the differentiating gate for a scale-bearing rung. The N=100 prior rung
(frozen, reused) is compared against the N=250 current rung using an identical,
apples-to-apples methodology: the per-assembly QC reconciliation build wall, measured
with `perf_counter`, amplified for stability, and taken as the minimum of 3 trials. The
fixed full-cohort source-CSV join scan is excluded as explained constant overhead.

| Metric | N=100 | N=250 | Change | Limit |
|--------|-------|-------|--------|-------|
| Build runtime (per-assembly QC) | 0.458 s | 1.161 s | — | — |
| Time exponent (objects) | — | — | **1.016** | ≤ 1.3 |
| Runtime per-base slope | — | — | **+1.8%** | ≤ 25% |

The runtime per-base slope is the scaling-determining metric and is stable
(`EXPLAINED_STABLE`). The compact per-assembly output sizes (`stage_bytes`,
`stage_files`) and fixed process RSS are reported and classified for transparency
(amortized/stable) and are bounded by the N=500 projection allocation checks — they are
not scaling bottlenecks for a correctness pilot.

**N=500 upper-95% projection** (power fit over N=100/N=250, max(fit, linear×2) + 25%):
wall 3 s, stage 0.36 MiB, RSS 280 MiB — all trivially within the N=500 allocation.

Gate: `scale_trend` = PASS; **`GO_500`**

---

## 11. Determinism and Kill/Restart

A deterministic rerun (a fresh independent build in a separate namespace) reproduced
**all 11 analytical and engineering compact units byte-identically**; only the
per-build timestamp (`release.json.created_at_utc`) and the freshly-measured scale-trend
wall differed, as expected. The release's static units are SHA-256 validated on resume;
a digest mismatch refuses mixed publication.

A **safe interruption was forced at a new build/query stage** (after SYNG, before
joins/query): no partial/final release was ever published, the staging directory
retained its SHA-validated units, and the resumed run completed cleanly to `COMPLETE`
with no mixed output.

| Check | Result |
|-------|--------|
| Injected interruption at build/query stage | observed (`AFTER_SYNG_BEFORE_JOINS_QUERY`); clean resume |
| Partial `COMPLETE` ever present | never |
| Independent build byte-identical units | 11/11 analytical+engineering |

Gate: `deterministic_rerun` = PASS; `injected_kill_restart` = PASS

---

## 12. Resource Utilization

| Resource | Allocation | Peak Measured | Fraction | Gate |
|----------|------------|---------------|----------|------|
| RAM (assigned) | 64 GiB | 112 MiB | 0.17% | PASS |
| Durable disk | 200 GiB | 0.13 MiB | <0.001% | PASS |
| Scratch disk | 2 TiB | — | <0.001% | PASS |
| Inodes | 1,000,000 | 16 files | <0.002% | PASS |
| Durable free preflight | ≥2 TiB | 2.58 TiB | — | PASS |
| Scratch free preflight | ≥4 TiB | 5.49 TiB | — | PASS |
| Swap growth | 0 | 0 (start = finish) | — | PASS |

Gate: `resource` = PASS

---

## 13. Atomic Promotion and Global Cap

The release directory was staged at `.integrated-pilot-250-v1-31cd0f3658837296.staging`,
`SHA256SUMS` written, `COMPLETE` created and fsync'd, then atomically promoted via
`os.replace()` to the final location. Consumers reject absent `COMPLETE`. The union of
all sequence-bearing downloads, canonical objects, and biological analyses across the
graph remains a subset of the frozen cohort and does not exceed 1,000 distinct exact
*E. coli* assembly revisions (actual: 250).

Gate: `atomic_promotion` = PASS; `global_distinct_assembly_cap` = PASS

---

## 14. Artifacts and Manifests

### Git-Committed (compact only; all < 10 MiB; largest qc_results.json ≈ 113 KiB)
- `workflow/integrated_pilot_250/` — workflow code and tests (read-only pinned plan)
- `manifests/integrated-pilot-250-v1/` — release manifests
- `artifacts/integrated_pilot_250/` — validation logs, metrics, controls
- `reports/integrated_pilot_250.md` — this report

### External Storage (bulky/sequence-bearing)
- Durable: `/home/erikg/phind-data/ecoli26k/v1/releases/run-integrated-250-genome/integrated-pilot-250-v1-31cd0f3658837296/`
- Scratch: `/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-integrated-250-genome/pilot-250-full-001/`

Every external release contains: exact input manifest, output `SHA256SUMS`, append-only
`state.jsonl`/`failures.jsonl`, provenance/tool/argv/resource logs, and atomic
`COMPLETE` promotion.

---

## 15. Root Input Immutability (Start & Finish)

| File | SHA-256 (Start) | SHA-256 (Finish) | Status |
|------|-----------------|------------------|--------|
| 26k_ecoli_accession.txt | 1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5 | 1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5 | UNCHANGED |
| 26k_prophage1.csv | 6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996 | 6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996 | UNCHANGED |

---

*End of Report*
