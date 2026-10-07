# NTM v3.1 — genome-path reconstruction: results and honest verdict

**Date:** 2026-10-07 · **Scope:** replaces the **ML step only**. Clades, host
clades, partitions and annotation are unchanged from v3.

## Why

v3's ML genomes are mosaics, not genomes (see the correction in
`ntm/v3/RELEASE.md` and `v3_mosaic_audit.tsv`):

- the reconstruction walk only *biased* toward observed adjacency
  (`w * (1 + beta*adj)`), never *required* it, so it accumulated high-support
  blocks until the 150 kb budget stopped it;
- 747/1,304 clades (57%) exceeded 1.5× the median member length, up to 6.2×,
  covering 97% of prophages;
- **387 of the 395 CheckV "high-quality" genomes were mosaics (98%)** — CheckV
  completeness *rewards* the chimera, because concatenating every common block
  includes every marker gene.

## Changes (all backward-compatible; `--path-mode free` is byte-identical to v3)

| flag | effect |
|---|---|
| `--path-mode observed` | walk restricted to **observed adjacency** edges; terminates at a natural **terminus**; no budget fallback |
| `--target-length N` | per-clade genome length target (CheckV `aai_expected_length` median over members), [0.8, 1.2]× guards, closest-to-target selection |
| `--close-circle` | optional circular closure, guarded (off by default) |

`--split-mode medoids` (structure-driven clustering, replacing the arbitrary
`--max-size 100`) is implemented and tested but **not applied to production** —
the mosaic defect lives in the ML step, not the clustering.

Evidence: free mode byte-identical to HEAD; 61 tests pass.

## Length: how we know it

We don't know phage length a priori. Measured against CheckV's marker-based
`aai_expected_length`, the member-length proxies are biased:

| estimator | vs CheckV expected |
|---|---|
| member **median** | **0.78× — biased LOW** (fragment calls) |
| member **max** | overshoots (rule v2 → 1.45× expected) |
| **CheckV `aai_expected_length`** | the target adopted |

## A/B (25 worst clades, same partitions, only the walk varies)

| rule | length ÷ expected | median completeness |
|---|---:|---:|
| v3 free (mosaic) | ~2.7× | 99.6% *(artifact)* |
| median-length match | 0.78× | 80.6% |
| max-member bound | 1.45× | 94.1% |
| **v3.1 (CheckV target)** | **0.99×** (19/25 in ±20%) | 82.5% |

## Full set (1,304 clades)

**High completeness and true genome length are mutually exclusive.** Reaching
94% requires 1.45× the expected length — i.e. padding. At genome scale the
honest number is lower:

| group | n | v3 completeness | v3 "HQ" | **v3.1 completeness** | v3.1 HQ |
|---|---:|---:|---:|---:|---:|
| singleton (pass-through) | 411 | 13.3% | 7 | 13.3% | 7 |
| small (2–10) | 432 | 24.9% | 35 | 22.6% | 26 |
| mid (11–50) | 125 | 80.5% | 52 | **30.2%** | 8 |
| large (51–100) | 336 | **100.0%** | 301 | **73.3%** | 83 |

ML length: median **14,324 bp** (v3: ~150 kb); median 1.00× target; 820/1,304
within ±20%; only **2/1,304** trip CheckV's ">1.5× expected" warning.

### Verdict

**v3.1 does not recover more complete genomes — it reveals that few were ever
complete.** v3's 100% in large clades was the mosaic artifact; the real number is
73.3%. v3.1 yields 124 high-quality genomes vs v3's 395, but 98% of those 395
were chimeras. v3.1 is the honest catalog; it is not a larger one.

## Remaining lever

**697/1,304 clades used the `member_seq_median` fallback target**, measured at
0.78× the true expected length — so those reconstructions are ~22% too short,
which is most of the mid-clade drop (80.5% → 30.2%). The fallback fires where
members are too fragmentary for CheckV AAI. A two-pass refinement (reconstruct →
CheckV AAI on the product → retarget) is the obvious next improvement.

Separately: many source prophage calls are themselves fragments, which bounds
what any stitching method can recover.

## Chimerism — the metric that matters for synthesis

`path_adj_fraction = 1.0` proves each step is an observed adjacency *somewhere*,
not that the whole path came from *one* genome. `ntm/v3.1/chimerism.py` measures
that directly against `partitions.bed` (per-member partition order):

| | |
|---|---|
| **single-member** (path == one genome's order) | **408 / 893 = 45.7%** |
| longest single-member run, median | **85.7%** of the path |
| members needed to cover the path | median **2**, p90 12, max 69 |

So 45.7% of reconstructions are *observed* sequences (a single prophage's partition
order, safe to synthesise as-is), and the remainder join 2+ genome contexts —
median 86% of their length still one contiguous run, but with at least one junction
no single genome contains. Per-clade verdicts: `$NVME/ntm/v3.1/chimerism.tsv`
(`verdict` ∈ {single_member, chimeric}, `longest_run_frac`, `n_members_to_cover`).

**Release implication:** carry `verdict` and `longest_run_frac` per genome, and
flag chimeric junctions explicitly — for a synthesis customer a novel junction is a
different product from an observed fragment.

## Deliverables

`$NVME/ntm/v3.1/catalog/` (`$NVME = /mnt/nvme3n1/erikg/phind-genome-work`):

| file | |
|---|---|
| `all_ntm_v31_ml_phage_genomes.fa.gz` (+ `.fai`/`.gzi`) | 1,304 ML genomes, bgzip, sha256 `7bd3e150…` vs uncompressed |
| `all_ntm_v31_ancestral_phage_genomes.fa.gz` | 893 ancestral genomes |
| `release_manifest.tsv` | 1,304 rows: `target_len`, `target_source`, `length_ratio_vs_target`, `checkv_expected_length`, `checkv_completeness`, `checkv_quality`, `n_partitions`, `terminated_by` |

Repo: this report, `SPEC.md`, `ab_summary_v3_vs_v31.tsv`, `v3_mosaic_audit.tsv`,
and the drivers (`ab_subset.py`, `build_v31_targets.py`, `run_v31_ml.py`).
