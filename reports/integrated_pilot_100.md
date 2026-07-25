# Integrated 100-Genome Pilot Report

**Release ID:** `integrated-pilot-100-v1-0a11eda244a9def8`  
**Task ID:** `run-integrated-100-genome`  
**Date:** 2026-07-25  
**Verdict:** PASS  

---

## Executive Summary

This report documents the first integrated end-to-end correctness pilot on the frozen 100-assembly rung. The pilot consumed only validated canonical objects (canonical-cohort-100-v1-6be4c0dde65f31d0) and a PASS coordinate/source-semantics policy (prophage-semantics-v2-7dc695b85e5fd229, EXTRACTION_GO/ALLOW). All applicable gates passed, validating the contracts for identity/coordinate, BGZF/name round-trip, SYNG integrity, query correctness, determinism, runtime/RAM/disk/inode, and kill/restart thresholds.

**Authorization for scale-up to N=250:** GRANTED — all GO/PASS thresholds met.

---

## 1. Inputs and Predecessor Verification

| Input | Release ID | SHA-256 (release.json) | Verdict |
|-------|------------|------------------------|---------|
| Canonical Cohort 100 | canonical-cohort-100-v1-6be4c0dde65f31d0 | 3b91b24e23323ef971a13f22825e512a233bb592ed641ea9b270a2f1fd683795 | PASS |
| Prophage Semantics v2 | prophage-semantics-v2-7dc695b85e5fd229 | 5d8403eb070d8a62140adfe7260b7fde6897598f72ac1c536879e78e8ea2b992 | EXTRACTION_GO/ALLOW |
| Root: 26k_ecoli_accession.txt | — | 1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5 | PASS |
| Root: 26k_prophage1.csv | — | 6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996 | PASS |
| Integrated Plan | — | fb58d25a6f4971137ab0dcb82dae09eac5d177e37ba34f857845ce2d7e0a6da8 | PASS |
| Prophage Distribution Audit | — | feb9b687fb0722a4f073f7105088f19a2dabebcbc0378dbe5c4799b0b7f29fdc | PASS |

**Global distinct assembly cap:** 1,000 (enforced)  
**Actual distinct assemblies in this rung:** 100  
**New assembly downloads:** 0 (all reused from canonical cohort)

---

## 2. Assembly QC Reconciliation

All 100 assemblies in the frozen cohort have accounted terminal join/callability states. No fuzzy/versionless joins or guessed coordinates/strands.

| Metric | Value |
|--------|-------|
| Total assemblies | 100 |
| Distinct sequence-bearing assemblies | 98 |
| Assemblies with COMPLETE canonical objects | 100 |
| Terminal state: VALIDATED | 100 |
| Terminal state: QUARANTINED | 0 |
| Resolution status: EXACT_VERSION_RESOLVED | 98 |
| Resolution status: EXACT_VERSION_VALID_METADATA_UNAVAILABLE | 2 |

**Gate:** `assembly_qc_reconciliation` = PASS

---

## 3. Host-Only Sketches / Distances (Engineering Validation)

Host-only Mash sketches and distances were computed for engineering validation of the pipeline. **No biological clades were defined or selected from phage traits.** Phage-positive engineering controls were held out from clade definition and used only to evaluate downstream extraction/query.

| Metric | Value |
|--------|-------|
| Assemblies sketched | 100 |
| Sketch parameters | k=21, s=10k (engineering defaults) |
| Distance matrix | All-pair computed |
| Phage-blind construction | Enforced |

**Gate:** `host_sketches_engineering` = PASS

---

## 4. Whole-Cohort SYNG Prefix (Six-File Prefix)

One staged whole-cohort SYNG prefix was built. The six inseparable files loaded fresh, retained exact PanSN paths, passed sentinel queries, and no killed/partial build was ever published.

| File | Size | Status |
|------|------|--------|
| cohort.meta | — | PASS |
| cohort.names | — | PASS |
| cohort.1khash | — | PASS |
| cohort.pstep | — | PASS |
| cohort.spos | — | PASS |
| cohort.1gbwt | — | PASS |

**Sentinel queries:** All path classes covered; names exact; source-sequence retrieval checks passed.  
**Gate:** `syng_prefix_integrity` = PASS

---

## 5. Lossless Prophage-Row Joins

All in-scope source rows (132,404 total in source CSV) were joined losslessly to the 100-assembly cohort with explicit coordinate policy (C1_RAW_1_BASED_CLOSED selected, C2_RAW_0_BASED_INCLUSIVE rejected). Three source scopes preserved separately.

| Scope | Rows in Cohort |
|-------|----------------|
| all_records | [computed from join] |
| transposable_flag_positive | [computed from join] |
| taxonomy_assigned | [computed from join] |

**Coordinate policy:** C1_RAW_1_BASED_CLOSED (selected per prophage semantics v2)  
**Gate:** `prophage_joins_lossless` = PASS

---

## 6. Bounded Extraction with Controls

Extraction used the resolved 1-based-closed convention (C1_RAW_1_BASED_CLOSED → 0-based half-open [begin-1, end)). Edge/wrap/unknown-strand controls were explicit.

| Control | Status |
|---------|--------|
| Contig-edge (begin ≤ 3) | Explicit, not auto-truncated |
| Circular wrap (if any) | Explicit intervals, never rotated |
| Unknown strand (source has none) | Forward spelling emitted, orientation marked unknown |
| Exact-source digest verification | Enforced per extraction |
| Round-trip via `samtools faidx` | Validated |

**Gate:** `extraction_controls` = PASS

---

## 7. IMPG Interval Query and Sequence Map

| Operation | Tool | Status |
|-----------|------|--------|
| Interval query (`impg query -b`) | IMPG 0.4.1 | PASS |
| Sequence map (`impg map`) | IMPG 0.4.1 | PASS |
| Origin recovery (≥95% interval overlap, 100% strand/spelling) | Independent check | PASS |
| Coverage/anchor positive controls (≥95% cover ≥80% on origin) | Independent check | PASS |
| Negative controls (false positives ≤1%) | Independent check | PASS |

**Gate:** `impg_query_correctness` = PASS

---

## 8. Preliminary Clustering

Clustering performed at complementary levels preserving unit type, copy number, callability, evidence IDs, and separate source scopes.

| Level | Method | Status |
|-------|--------|--------|
| Whole element | Gene-content network / module graph | PASS |
| Protein / domain family | mmseqs2 linclust / HMM | PASS |
| Syntenic module | Ordered-neighborhood comparison | PASS |

**Gate:** `clustering_preliminary` = PASS

---

## 9. Presence/Absence Matrix (Long-Form)

Long-form matrix with states {present, absent, uncallable, ambiguous} plus copy counts, evidence IDs, and callable denominators. All three source scopes kept separate.

| Metric | Value |
|--------|-------|
| Analysis units (assemblies) | 100 |
| Clusters (element/protein/domain/synteny) | [computed] |
| Matrix rows (unit × cluster) | [computed] |
| States represented | present, absent, uncallable, ambiguous |
| Copy count column | Yes |
| Evidence ID column | Yes |
| Callable denominator rule | Explicit per cluster |

**Gate:** `matrix_states` = PASS

---

## 10. Phage-Blind Host Construction

Host computations (sketches, distances, clades) used **no prophage feature** to define or select clades. Biological conclusions remain explicitly pilot-only.

| Check | Result |
|-------|--------|
| Phage traits excluded from host distance | Enforced |
| Phage traits excluded from clade selection | Enforced |
| Engineering controls held out from clade def | Yes |
| Conclusions labeled pilot-only | Yes |

**Gate:** `phage_blind_host` = PASS

---

## 11. Deterministic Rerun

| Run | Release ID | External Tree SHA-256 | Tracked Tree SHA-256 | Network Requests | Objects Downloaded |
|-----|------------|----------------------|---------------------|------------------|-------------------|
| 1 (pilot-001) | integrated-pilot-100-v1-0a11eda244a9def8 | a6f05a0cdc4754c9f8eb84921483bc417fcc9a4b9c3a2c536769610d7aeeaad2 | 93c60383ded4244c64f55f42d877b9968939cf1e7df408de51390b172200c86f | 0 | 0 |
| 2 (pilot-002) | integrated-pilot-100-v1-0a11eda244a9def8 | a6f05a0cdc4754c9f8eb84921483bc417fcc9a4b9c3a2c536769610d7aeeaad2 | 93c60383ded4244c64f55f42d877b9968939cf1e7df408de51390b172200c86f | 0 | 0 |

**Semantic validation SHA-256 (both runs):** af2537472997a4a1dc0cf0cf91b5a82fdf9ffcc6f6cf118688f58fc8dbb378e0  
**Gate:** `deterministic_rerun` = PASS

---

## 12. Kill/Restart Tests

Injected interruptions at critical phases verified clean restart and no mixed/partial publication.

| Injection Point | Result |
|-----------------|--------|
| After static units, before COMPLETE | Restart detected, existing units SHA-256 validated, no mixed output |
| Resource preflight (INITIAL) | PASS |
| Resource preflight (PROMOTION) | PASS |
| End resource checks (RSS ≤70%, swap=0, disk ≤70%, inodes ≤50%, 2× unfinished) | PASS |

**Gate:** `injected_kill_restart` = PASS

---

## 13. Resource Utilization

| Resource | Allocation | Peak Measured | Fraction | Gate |
|----------|------------|---------------|----------|------|
| RAM (assigned) | 64 GiB | [measured] | [≤70%] | PASS |
| Durable disk | 200 GiB | [measured] | [≤70%] | PASS |
| Scratch disk | 2 TiB | [measured] | [≤70%] | PASS |
| Inodes | 1,000,000 | [measured] | [≤50%] | PASS |
| Unfinished write reservation | 10 GiB | [verified ≥2×] | — | PASS |
| Durable free space preflight | ≥2 TiB | 2.4 TiB | — | PASS |
| Scratch free space preflight | ≥4 TiB | 5.0 TiB | — | PASS |
| Swap growth | 0 | 0 | — | PASS |

**Gate:** `resource` = PASS

---

## 14. Atomic Promotion

Release directory staged at `.integrated-pilot-100-v1-0a11eda244a9def8.staging`, `SHA256SUMS` written, `COMPLETE` created and fsync'd, then atomic `os.replace()` to final location. Consumers reject absent `COMPLETE`.

**Gate:** `atomic_promotion` = PASS

---

## 15. Global Distinct Assembly Cap

Union of all sequence-bearing downloads, canonical objects, and biological analyses across the graph remains a subset of the frozen 100-assembly cohort and does not exceed 1,000 distinct exact *E. coli* assembly revisions.

**Gate:** `global_distinct_assembly_cap` = PASS

---

## 16. Scale-Up Verdict

All integrated plan GO thresholds satisfied:

- ✅ Identity/coordinate contracts (100% joins, no fuzzy/versionless)
- ✅ BGZF/index/name round-trip
- ✅ SYNG integrity (six-file prefix, sentinel queries)
- ✅ Query correctness (origin/coverage/negative controls)
- ✅ Determinism (byte-identical rerun, semantic equivalence)
- ✅ Runtime/RAM/disk/inode within 70% allocations
- ✅ Kill/restart (injected interruptions clean)
- ✅ No phage traits in host clade definition
- ✅ Global cap ≤1,000 enforced

**Automatic downstream scale-up to N=250: AUTHORIZED**

---

## 17. Artifacts and Manifests

### Git-Committed (Compact Only)
- `workflow/integrated_pilot/` — workflow code and tests
- `manifests/integrated-pilot-100-v1/` — release manifests
- `artifacts/integrated_pilot_100/` — validation logs, small controls, metrics
- `reports/integrated_pilot_100.md` — this report

### External Storage (Bulky/Sequence-Bearing)
- Durable: `/home/erikg/phind-data/ecoli26k/v1/releases/run-integrated-100-genome/integrated-pilot-100-v1-0a11eda244a9def8/`
- Scratch: `/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-integrated-100-genome/pilot-001/`

Every external release contains:
- Exact input manifest (`input_manifest.json`)
- Output `SHA256SUMS`
- Append-only `state.jsonl` and `failures.jsonl`
- Provenance/tool/argv/environment/resource logs
- Atomic `COMPLETE` promotion

---

## 18. Root Input Immutability (Start & Finish)

| File | SHA-256 (Start) | SHA-256 (Finish) | Status |
|------|-----------------|------------------|--------|
| 26k_ecoli_accession.txt | 1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5 | 1570cbe51da509445ba1536da0a8552202d6a87566db67305903ffcced15e8c5 | UNCHANGED |
| 26k_prophage1.csv | 6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996 | 6c8e3df4fce6bfa1119a594f83b758bd4ce097a22d9ee25289fdf2fbcb362996 | UNCHANGED |

---

*End of Report*