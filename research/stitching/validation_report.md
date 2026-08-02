# Validation Report: Partition Stitching Algorithm (v2 — validated re-audit)

**Generated:** 2026-08-02 (v2; v1 dated 2026-07-31 at commit `2363ece`)
**Task:** `stitching-validation` (agent-6)
**Community:** 3 (1,234 sequences)
**Commit under audit:** `2363ece` (build-partition-stitching)
**Inputs (now committed in-repo):** `research/stitching/inputs/community_3_partitions/`
(502 MAF alignment blocks, sha256-verified identical to the archived originals),
`research/stitching/inputs/3_ancestral.fa` (113,502 bp),
`research/stitching/community_3_partitions.bed` (973 prophages)

---

## 0. Executive summary

| Question | Verdict |
|---|---|
| Reproducible run from BED+MAF? | ✅ **YES** — end-to-end re-run is byte-identical (sha256) to the committed outputs at both thresholds (50% and 45%). Entry point: `research/stitching/run_community3_validation.sh`. |
| Threshold semantics (50% vs 45%) | ✅ **Documented** — the threshold is a biological parameter (`accessory_threshold` = min fraction of prophages carrying a partition); the occurrence distribution has a natural knee at ~45–49%; **both** thresholds produce genomes far longer than the true community mean (9.6 kb) because partition-block consensi are union mosaics (see §4, §8). The 45% "fallback" was justified by a **wrong length target** (see §7). |
| Overlap handling (no suffix/prefix merge) | ✅ **CORRECT for these partitions** — all adjacent core-partition intervals are exactly contiguous (0 bp overlap, 0 bp gap) in every co-occurring prophage (70/70, 52/52, 86/86 pairwise adjacencies). Partitions are sequential alignment blocks; concatenation is the right model. The 10-N join markers are cosmetic (see §5). |
| Identity claim (78.05 % MASH) | ❌ **NOT SUPPORTED by direct alignment.** MASH 78.05 % is a uniform-divergence model extrapolation from a 0.5 % k-mer Jaccard. Direct alignment (minimap2/wfmash) finds only ~6–9.5 kb (5–8 % of the ancestral) at 30–44 % identity; partition 48 (33.5 kb of the stitched genome) shares **zero** detectable homology with the ancestral (see §6). |
| Ancestral genome provenance | ✅ **Provenance confirmed; validity rejected.** `3_ancestral.fa` = pggb cluster-3 ancestral from `scripts/pggb_per_cluster_pipeline.py` (`5a22618`), built from the **same 1,234 sequences** as community 3 — but it is a naive index-position consensus padded to the maximum path length (502→113,456 bp paths), i.e. dominated by the single longest prophage (conf 0.484). It is **not** a valid community mean genome (see §7). |
| Verdict | **REPAIRABLE** — algorithm code is deterministic, reproducible, and mechanically correct (REUSABLE as a pipeline); the **scientific claims** in v1 (mean genome ~111 kb, 78.05 % identity, 1.10× ancestral) are **not supported** and require repair: honest threshold/coverage semantics, alignment-based identity reporting, and a valid reference (see §10). |

---

## 1. Reproducibility (verified end-to-end)

Re-ran `stitch_algorithm.py` (after additive fixes, §9) from the committed in-repo
inputs with the exact parameters of the original committed runs
(`coverage_threshold=0.0`; accessory threshold 0.5 and 0.45):

| Output | Committed sha256 (v1) | Re-run sha256 | Match |
|---|---|---|---|
| `community_3_stitched_mean.fa` (50 %) | `f4e1c0f379fe8b62c8479e0cac4cae2168aa3358cd2c6e206a59d28c5de616fd` | `f4e1c0f379fe8b62c8479e0cac4cae2168aa3358cd2c6e206a59d28c5de616fd` | ✅ byte-identical |
| `community_3_stitched_mean_45pct.fa` (45 %) | `ad09fe7140701668c49438c04fa4e7c5184b9337603a37c0290c25e7d3ed9601` | `ad09fe7140701668c49438c04fa4e7c5184b9337603a37c0290c25e7d3ed9601` | ✅ byte-identical |
| `stitching_results.json` / `_45pct.json` | (v1 schema) | extended schema (§9), same core values | ✅ values identical, fields added |

Runtime ~14 s per threshold. This closes the reproducibility gap flagged by
`audit-recovered-artifacts` ("the stitch algorithm … cannot currently be re-run
end-to-end from main alone") — the 502 MAF inputs and the ancestral FASTA are
now committed under `research/stitching/inputs/`.

## 2. What the algorithm actually does

1. Parse the BED (`prophage → ordered partition intervals`).
2. Count adjacencies between consecutive partitions per prophage → directed graph.
3. Greedy max-likelihood path over partitions with occurrence ≥ threshold.
4. Per-partition consensus from the partition MAF alignment block
   (majority rule per column, `coverage_threshold` gate on non-gap depth).
5. Concatenate core-partition consensi in path order (`N`×10 join markers).

The output is best described as an **ordered pangenome-core path**: the
concatenation of the consensus of the most-common partition *bundles*, not a
"mean genome" in the sense of a single typical member's genome (§8).

## 3. Threshold semantics — decision (concern 1)

**Semantics (unchanged, now documented in the code header):**
`accessory_threshold` = minimum fraction of prophages with BED assignments
(973) in which a partition must appear to be included in the core path.
0.5 = strict-majority "core genome"; lower values relax toward a
pan-genome-style union of regions.

**Occurrence distribution (973 prophages with assignments):**

| Partition | Occurrence | Fraction |
|---|---|---|
| 51 | 746 | 76.7 % |
| 239 | 476 | 48.9 % |
| 241 | 468 | 48.1 % |
| 48 | 451 | 46.4 % |
| 55 | 428 | 44.0 % |
| …261 partitions | <20 % | long tail |

**Decision.** There is a natural knee between 48.9 % and 76.7 % (fractions are
BED-occurrence / 973 prophages — the metric the algorithm uses):
- threshold ∈ (48.9 %, 76.7 %] → **only partition 51** (53,886 bp; 0.47× ancestral) — the strict majority core;
- threshold ∈ (48.1 %, 48.9 %] → partitions **51, 239**;
- threshold ∈ (46.4 %, 48.1 %] → partitions **51, 239, 241**;
- threshold ∈ (44.0 %, 46.4 %] → partitions **48, 51, 239, 241** (124,935 bp; the v1 "45 %" result);
- threshold ≤ 44.0 % adds partition 55, then a long tail of low-occurrence partitions.

**Correct use:** the threshold is a **parameter to be set from the biology or a
stated objective**, not a constant. The v1 report's justification ("45 % chosen
to hit ~111 kb / >70 % identity") is **invalid** because the target was derived
from a broken reference (§7) — and, independently, **no threshold on this
partition set yields a mean-genome-sized result**: the true community mean is
9,550 bp (§8) while the outputs range 53.9–1,164 kb. If the objective is a
*strict core genome*, 50 % is correct (1 partition). If the objective is an
*extended core / pan-genome path*, pick the knee (~45 %) and say so. The 45 %
output remains the best-documented "extended core path" for community 3.

## 4. Overlap handling — verdict (concern 2)

**Verdict: the no-overlap (concatenation) model is CORRECT for these
partitions.**

- For every prophage carrying both members of an adjacent core pair, the two
  BED intervals are **exactly contiguous**: 48→51 (70/70 prophages), 51→239
  (52/52), 239→241 (86/86) — 0 bp overlap, 0 bp gap, in both orientations.
- The partitions are therefore **sequential alignment blocks**, as the v1
  report claimed; a suffix/prefix overlap merger is unnecessary.
- The `find_overlap()` code never triggers (0 overlap joins at 45 %; confirmed
  by the join counters added in this re-audit). The `N`×10 join markers are
  cosmetic bookkeeping between independent block consensi; the junction
  sequence itself is not represented in either consensus. `--gap-size 0`
  gives a plain concatenation if desired (default stays 10 → byte-identical).

## 5. Identity re-check — MASH vs direct alignment (concern 3)

**Claim under test:** "MASH identity vs ancestral 78.05 %" (45 % threshold).

**What MASH actually measured:** `mash dist` returned distance 0.219531
(5 / 1000 shared hashes). The script reports `1 − dist = 78.05 %`. With
k=21, this is the MASH model extrapolation from a **k-mer Jaccard of only
0.5 %** and is only interpretable as a divergence estimate under a uniform
substitution model — which does not hold here (see below).

**Independent exact-k-mer check (k=21, no sketching):** only 1,103 of
~106,000 stitched-genome 21-mers (1.04 %) occur in the ancestral. These are
not spread uniformly — they are concentrated in one window:
stitched 110–120 kb ↔ ancestral 5–15 kb (i.e. almost entirely partition 241),
as short collinear runs of 3–47 bp — the signature of a moderately-diverged
patch, not genome-wide 78 % identity.

**Direct alignment (minimap2 `map-ont`; wfmash p=60/50 as cross-check):**

| Query (stitched) | Aligned to ancestral | Identity | Notes |
|---|---|---|---|
| whole stitched genome (124,935 bp) | ~6.2 kb (minimap2) – 9.5 kb (wfmash) = 5–8 % of ancestral | ~34–59 % on aligned portion | mostly partition 241 ↔ anc 3.8–8 kb |
| partition 48 (33,533 bp) | **0 bp** | — | no detectable homology |
| partition 51 (53,886 bp) | ~1.0–1.4 kb | ~11–97 % | one fragment, anc ~12.3 kb |
| partition 239 (21,817 bp) | ~1.1 kb (3 fragments) | ~18 % | anc ~15–21 kb |
| partition 241 (15,669 bp) | ~4.2 kb | ~44 % | anc ~3.8–8.0 kb |

**Conclusion: the 78.05 % MASH claim is NOT supported.** Direct alignment
shows the stitched genome and the ancestral share only a small homologous
patch (partition 241 ↔ ancestral ~5–15 kb) plus tiny fragments; 90 % of the
stitched genome (partitions 48+51+239) has no credible homology to the
ancestral. MASH's number is an artifact of interpreting sparse k-mer overlap
under a model that does not apply to these mosaic sequences. **Any future
identity claim must be alignment-based** (PAF matches/aligned-bases), with
coverage (fraction of the reference aligned) reported alongside identity.

## 6. Ancestral genome provenance (concern 4)

**Provenance (confirmed):** `3_ancestral.fa` (header `3_ancestral_consensus`,
113,502 bp, single contig, sha256 `da210d97fb341e69f65507ecf3023422f8750d7ff1b85e9112020b9c5c17ee00`)
is the pggb **cluster-3 ancestral** from the `5a22618` per-cluster pipeline
(`scripts/pggb_per_cluster_pipeline.py` → `reconstruct_ancestral_genome()`),
produced by `odgi paths -f` followed by a **position-wise majority consensus
over paths padded to max length**. The 1,234 sequences in pggb cluster 3 are
identical (by ID) to community 3's `community_3.fa`.

**Validity: REJECTED as a mean-genome reference.** The reconstruction aligns
paths **by index**, which is invalid for variable-length paths: cluster-3 path
lengths range 502 → 113,456 bp, so the consensus is padded to the longest path
and dominated by the single 113 kb prophage (that prophage maps to the
"ancestral" at 80.8 % over 42 kb; a typical 7.5 kb prophage maps at 22 % over
35 % of its length). Its own reported confidence is 0.484. Consequences:

- The task expectation "mean genome ~111 kb / >70 % identity to ancestral"
  was **derived from this broken reference** — the real community-3 mean
  prophage is ~9.6 kb (§8).
- The v1 report's length ratio (1.10× ancestral) and identity target were
  chasing an artifact.

**Recommendation:** do not use `3_ancestral.fa` (or any of the
`ancestral_genomes/*` consensus files built the same way) as a community mean
genome reference. The pggb *graphs* themselves (`5a22618` output) may be
reusable; the naive positional consensi are not.

## 7. Community-3 reality check — why the lengths are what they are

- 1,234 prophages: **mean 9,550 bp, median 7,527 bp** (755 seqs < 10 kb,
  413 in 10–20 kb, ~65 > 20 kb, max 113,456 bp). The "~111 kb" expectation
  matched only the single longest sequence.
- 973/1,234 (78.8 %) have partition intervals in the BED; 261 prophages
  (21 %) are absent from the partition graph.
- Partition MAF blocks are **staggered-segment bundles**, not standard MSAs:
  e.g. partition 51 = 701 segments of median 2.8 kb scattered across a
  54 kb block (mean column depth 33/701; only 3.1 % of columns ≥ 25 % depth).
  The consensus at `coverage_threshold=0.0` therefore emits the **union
  mosaic** of all segment positions (every column with ≥ 1 base), which is why
  a 7.5 kb prophage's 3 partitions sum to 92 kb of "consensus".
- Net effect: the stitched 124,935 bp (45 %) ≈ 13× the true mean genome, and
  even the strict-core 53,886 bp (50 %) ≈ 5.6× the mean. **The stitched
  output is an ordered pangenome-core path, not a mean genome.**

## 8. Code changes applied (`stitch_algorithm.py`, all additive/backward-compatible)

1. **Threshold & consensus semantics documented** in the module docstring
   (what `accessory_threshold` and `coverage_threshold` mean, and the
   staggered-bundle caveat).
2. **Column-depth diagnostics** in `compute_partition_consensus()`: per-partition
   n_seqs, width, mean/max depth, fraction of columns covered / passing the
   threshold; `run_stitching()` flags blocks whose mean depth ≪ n_seqs as
   staggered bundles. (All 4 core partitions are flagged.)
3. **Join accounting**: `stitch_and_merge()` now returns (seq, n_overlap_joins,
   n_gap_joins); the stderr/stdout log and JSON report overlap vs gap joins.
4. **`--gap-size`** (default 10) controls the N-run join marker; 0 = plain
   concatenation.
5. **Richer JSON**: adds `core_path`, `n_overlap_joins`, `n_gap_joins`,
   `mash_identity_vs_ancestral`, and per-partition `consensus_stats`. Existing
   fields unchanged (values identical to v1 JSONs).
6. **`run_community3_validation.sh`**: reproducible entry point regenerating
   both FASTA outputs and both JSONs from the committed inputs.

FASTA outputs are **byte-identical** to v1 at the committed parameter sets
(verified, §1); only the JSON schema gained fields.

## 9. Verdict

**REPAIRABLE.**

- **REUSABLE (as-is, mechanically):** `stitch_algorithm.py` is deterministic,
  correct for its stated mechanics (BED → adjacency → path → block consensus →
  concatenation), and reproduces byte-identically; the overlap-handling model
  is correct for sequential impg partition blocks; the 45 % / 50 % outputs are
  reproducible artifacts.
- **NOT REUSABLE (as v1 interpreted them):** the v1 report's scientific claims
  — "mean genome ~111 kb" (wrong target; true mean 9.6 kb), "MASH identity
  78.05 %" (contradicted by direct alignment), "1.10× ancestral" (meaningless
  against a broken reference), and the 45 %-to-hit-111 kb threshold rationale —
  are **not supported** and must be replaced with the corrected framing above.
- **Repairs shipped in this re-audit:** threshold/coverage semantics documented;
  staggered-bundle diagnostics added; join bookkeeping added; identity claim
  re-checked by alignment; ancestral reference flagged as invalid; inputs
  committed so the pipeline re-runs end-to-end.

## 10. Generalization to other communities

**Mechanically: yes.** The pipeline needs, per community, (a) a partition BED
(`prophage → ordered partition intervals`) and (b) per-partition MAF alignment
blocks — exactly the outputs of the `impg partition` step. `audit-recovered-
artifacts` confirmed agent-223's `run-impg-partition` produced per-partition
consensus FASTA + BEDs for 12 communities (committed on branch `64ddebf`),
so the same stitching can be run community-wide. A FASTA-consensus input mode
(drop-in for the MAF reader) would be a small addition if only consensi, not
MAFs, are available.

**Scientifically: conditional.** The same caveats apply everywhere:
- The stitched output is a **pangenome-core path**, not a mean genome; its
  length must be compared against the community's *actual* prophage length
  distribution, not an assumed target.
- `coverage_threshold` must be chosen with the bundle structure in mind (the
  new diagnostics show whether a block is a true MSA or a mosaic).
- Identity must be reported from **alignments** (with coverage), never from
  MASH distance alone on such mosaic sequences.
- The `pggb_per_cluster_pipeline` positional ancestral consensi must **not**
  be used as references (naive index alignment over variable-length paths).

## 11. Coordination with `audit-recovered-artifacts`

| Artifact | audit-recovered-artifacts | This re-audit |
|---|---|---|
| `stitch_algorithm.py`, FASTA/JSON outputs | REUSABLE (byte-identical repro) | ✅ confirmed; JSON schema extended (values identical) |
| `validation_report.md` v1 | REPAIRABLE (2 doc bugs) | replaced by this report (4 substantive claims corrected) |
| 502 MAF inputs + `3_ancestral.fa` | not on main (recovery dir only) | **now committed** under `research/stitching/inputs/` |
| pggb cluster-3 ancestral (`3_ancestral.fa`) | REUSABLE (uncommitted, import) | imported, provenance confirmed, **validity rejected as mean-genome reference** |
| agent-223 `run-impg-partition` mean genomes (e.g. community 3 = 1,398,681 bp) | REPAIRABLE/partial, crude concatenation | consistent: both approaches over-inflate vs the 9.6 kb true mean; stitching (124.9 kb) is closer to the ancestral-span but still not a mean genome |
| agent-220 `run-pggb-ancestral`, agent-224, agent-225 | FAILED | unchanged |

## 12. Files

- `research/stitching/stitch_algorithm.py` — fixed (additive changes, §9)
- `research/stitching/run_community3_validation.sh` — reproducible entry point
- `research/stitching/inputs/community_3_partitions/` (502 MAF) + `inputs/3_ancestral.fa` — committed inputs
- `research/stitching/community_3_stitched_mean.fa` / `_45pct.fa` — regenerated (byte-identical)
- `research/stitching/stitching_results.json` / `_45pct.json` — regenerated (extended schema)
- `research/stitching/validation_report.md` — this report
