# Artifact Audit: Recovered PHIND Work (2026-08-02)

**Audit task:** `audit-recovered-artifacts`
**Auditor:** agent-1 (Programmer)
**Date:** 2026-08-02

This audit verifies which recovered PHIND artifacts/commits are genuinely
reusable vs failed attempts. Historical task terminal labels were **not**
trusted — every artifact was inspected and, where possible, reproduced.

## Sources inspected

| Source | Path | Role |
|---|---|---|
| Fresh checkout (main) | `/home/erikg/phind` @ `a5070e7` (20 commits ahead of origin) | read/write (report lives here) |
| Archived old checkout + worktrees | `/home/erikg/phind.recovery-20260801T164014Z` | read-only |
| Recovery evidence | `/mnt/nvme3n1/erikg/phind-recovery-20260801T164014Z` (manifests, task-agent-session-map.tsv, git bundle `phind-all-refs.bundle`, sha256 `234c74cff8169034bfa5f977a4c211445672ee7de10e474b48664bb7066729fd` — **verified, matches**) | read-only |

Cross-reference: `PHIND_RECOVERY_HANDOFF.md` labels `run-pggb-ancestral`,
`run-impg-partition`, `benchmark-find-fastest`, `run-all-wave` as **failed**
and `research-impg-partition` as **done** (but only in the archived worktree).
This audit confirms those labels *except* that `run-impg-partition`
(agent-223) left real committed outputs (11 community mean genomes + report)
that are absent from main.

---

## Verification Matrix

Status legend: **REUSABLE** = verified working/consistent; **REPAIRABLE** =
real but incomplete or needs minor fix/import; **FAILED** = no usable output;
**JUNK** = do not import.

| # | Artifact | Commit / ref | Status | Checksum evidence | Reuse recommendation |
|---|---|---|---|---|---|
| 1 | `research/stitching/stitch_algorithm.py` | `2363ece` (main) | **REUSABLE** | sha256 `0057f3f16ff21724b04bdb3344ba4518d429a1acd75cad28ce60d48fc3913ec8` (main == archived agent-226) | Use as-is; deterministic, 14 s for community 3 |
| 1 | `research/stitching/community_3_stitched_mean.fa` (50% result, 53,886 bp, path=1) | `2363ece` (main) | **REUSABLE** | sha256 `f4e1c0f379fe8b62c8479e0cac4cae2168aa3358cd2c6e206a59d28c5de616fd` — byte-identical repro | Reference output for 50% threshold |
| 1 | `research/stitching/community_3_stitched_mean_45pct.fa` (124,935 bp, path=4) | `2363ece` (main) | **REUSABLE** | sha256 `ad09fe7140701668c49438c04fa4e7c5184b9337603a37c0290c25e7d3ed9601` — byte-identical repro | Best current "mean genome" for community 3 |
| 1 | `research/stitching/community_3_partitions.bed` | `2363ece` (main) | **REUSABLE** | sha256 `d3259a15adf6391adb42e17ead6c5511d94067fdfe6623798d614c7382821b60` | identical to `prophage_homology_survey/mean_genomes/community_3_partitions.bed` |
| 1 | `research/stitching/stitching_results.json`, `stitching_results_45pct.json` | `2363ece` (main) | **REUSABLE** | sha256 `9353121b5b3b1471c80c700777bc9c30403093d57e7115cb80f94cf344f6d71d` / `241abe2f8a6c0f2e80f14f2b551be703c904733d88fa568dfe3a7297ea991017` — repro identical | attach |
| 1 | `research/stitching/validation_report.md` | `2363ece` (main) | **REPAIRABLE** | — | Two doc bugs (see §1.2); scientific claims all verified |
| 2 | `pggb_analysis/pipeline_results.json` | `5a22618` (main) | **REUSABLE** | sha256 `5b1b49476fea73d5255c6d066d99db429dc7b2b5f4760e12371ea196e6af90d5` | cluster 6 COMPLETE; no cluster 1 entry (incomplete) |
| 2 | `pggb_analysis/cluster_6/**` (gfa, aln.phy, iqtree_6, dist.tsv, ancestral/) | `5a22618` (main) | **REUSABLE** | `6_ancestral.fa` sha256 `f2745f42996aa7c84f6a3ace1aed0ec19cdd31321dee2ef1f4b217cb1b74f101` (main == archived) | full pipeline output; 69,153 bp / 8 modules / conf 0.5523 verified |
| 2 | `pggb_analysis/cluster_1/cluster_1.aln.fa` | `5a22618` (main) | **FAILED** | 0 bytes | MAFFT never finished; do not treat as alignment |
| 2 | `pggb_analysis/cluster_1/*.gfa, *.paf, *.paths.fa` | `5a22618` (main) | **REPAIRABLE** | — | graph built (189,663 seg / 8,235 paths) but 679 MB + 75 MB PAFs are git-tracked (see §5) |
| 2 | `scripts/pggb_per_cluster_pipeline.py` | `5a22618` (main) | **REUSABLE** | — | rerun with phage-appropriate params (handoff: `p=75 s=250 l=500 k=11 ani-diff=80`) |
| 2 | `pggb_analysis/ancestral_genomes/{3,5,7,8,9,11}_ancestral.fa` (+ modules.json) | **uncommitted** (archived agent-220 worktree + recovery root, untracked) | **REUSABLE (uncommitted — import ASAP)** | `3_ancestral.fa` sha256 `da210d97fb341e69f65507ecf3023422f8750d7ff1b85e9112020b9c5c17ee00` (both copies identical) | 6 real ancestral genomes that will be LOST when recovery dir is cleaned; cluster 3 (113,502 bp) is the reference used by stitching validation |
| 3 | `research_outputs/` scripts + 55 BED + 52 MAF + consensus + comparison (agent-221) | `0c8373c` (archived only) | **REUSABLE** | committed in git object `0c8373c` (present in fresh repo) | **IMPORT** (cherry-pick); see §3 |
| 3 | `research_outputs/cluster_3.syng.*` (index), `cluster_3.fa`, `cluster_3.alignments.paf` | uncommitted (agent-221 worktree) | **REPAIRABLE** | `cluster_3.fa` sha256 `d6a4e48901783fa7596ce63a8dbb97ea253d4d43340962dd9732965890895a20` (== agent-224 `community_3.fa`) | re-derive: syng index 458 ms; `cluster_3.fa` from `full_prophages.fa`; PAF regenerable |
| 3 | `research_outputs/partition_run.log` (109 MB), `partition_maf_run.log` (102 MB) | uncommitted (agent-221) | **JUNK** | — | do not import |
| 4 | agent-220 `run-pggb-ancestral` (branch head == main commit) | no commits | **FAILED** | its `pipeline_results.json`: cluster 1 `FAILED: pggb graph` | nothing reusable from the re-run itself; 11 `seqwish-*` temp dirs + `target/` are junk |
| 4 | agent-223 `run-impg-partition` — `prophage_homology_survey/mean_genomes/` (11 mean.fa, 2,721 partition .fasta, 11 .bed, report.md) | `5125aee` + `64ddebf` (archived only) | **REPAIRABLE** | committed on branch only; **not on main** | mean genomes are crude concatenations (community 3 = 1,398,681 bp vs 113,502 bp ancestral, 2.3 % N) — inferior to stitching; per-partition FASTA consensus files reusable as inputs; PDF/md reports fine as records |
| 4 | agent-224 `benchmark-find-fastest` | no commits | **FAILED** | `community_3.fa` sha256 `d6a4e489…` (regenerable) | only regenerable inputs; nothing to salvage |
| 4 | agent-225 `run-all-wave` | no commits | **FAILED** | — | empty worktree; nothing to salvage |
| 5 | `prophage_homology_survey/full_prophages.fa` | `ba841f4`/`6f96e89` (main) | **REUSABLE** | 3,265,823,570 B (3.04 GiB); git blob `af0fc973ec51d5b062e26339d5f4de8c74db74f1` matches working tree; sha256 `ed85b2fb549be18bc638d8485f5b5add7c2d394f3822efe66a90ca6d979758d3` | integrity fully verified; 132,393 seqs consistent with extraction log |
| 5 | `prophage_homology_survey/extraction_log.txt` | `ba841f4` (main) | **REUSABLE** | — | 132,404 expected − 11 FAIL = 132,393 ✓ (exact ID match) |
| 5 | `extract_prophages.py` | `ba841f4` (main) | **REUSABLE** | sha256 `7af1da4e6c13a34da85a08513d08a9720ea9744909e911917d2b01f2291f888c` (main == ba841f4) | reuse for re-derivation |

---

## 1. `2363ece` build-partition-stitching — REUSABLE (reproduced byte-identical)

### 1.1 Reproduction

Ran `research/stitching/stitch_algorithm.py` from fresh main with the exact
parameters of the committed runs (coverage_threshold=0.0; accessory 0.5 and
0.45) against the archived inputs:

- `--partition-dir /home/erikg/phind.recovery-20260801T164014Z/prophage_homology_survey/mean_genomes/community_3_partitions/` (502 MAF files)
- `--bed research/stitching/community_3_partitions.bed` (973 prophages; BED identical to the homology-survey copy)
- `--ancestral /home/erikg/phind.recovery-20260801T164014Z/pggb_analysis/ancestral_genomes/3_ancestral.fa` (113,502 bp)

**Result: SUCCESS in ~14 s per run.** All four outputs (`community_3_stitched_mean.fa`,
`community_3_stitched_mean_45pct.fa`, both `stitching_results*.json`) are
**byte-identical** (sha256 above) to the committed files. Reproduction outputs
attached in `research/stitching/repro/`.

> Note: the 502 MAF partition files and the ancestral FASTA are **not** on
> fresh main — they exist only in the archived recovery dir. The stitch
> algorithm itself is on main but cannot currently be re-run end-to-end from
> main alone.

### 1.2 Claims in `validation_report.md` — checked

| Claim | Verdict |
|---|---|
| 50 % threshold → 1 core partition (51, 76.7 % occurrence), 53,886 bp, 0.47× ancestral | **VERIFIED** (repro identical: path=[51], 53,886 bp, ratio 0.475) |
| 50 % row MASH identity "N/A" | **PARTIAL** — the script actually computes and prints `0.0000` (mash dist 1.0, 0/1000 shared hashes). "Too short for meaningful comparison" is the correct interpretation; "N/A" is imprecise. |
| 45 % threshold → 4 core partitions (48, 51, 239, 241), 124,935 bp, 1.10× ancestral | **VERIFIED** (repro identical: path=[48, 51, 239, 241], 124,935 bp, ratio 1.101) |
| MASH identity 78.05 % | **VERIFIED** (independent mash: dist 0.219531 → identity 78.05 %; 5/1000 shared hashes, p=3.7e-27) |
| Adjacency path 48→51→239→241 | **VERIFIED** (repro prints exactly this path; edges 48→51: 52, 51→239: 28, 239→241: 44 prophages) |
| 50 % vs 45 % threshold semantics | **VERIFIED** — `accessory_threshold` = fraction of 973 prophages; 50 % → 487 prophages → only partition 51 passes; 45 % → 438 → partitions 48/51/239/241 pass |

**Doc bugs found (REPAIRABLE):**
1. "Output Files" section says `community_3_stitched_mean.fa` is the "124,935 bp at 45 % threshold" — **wrong**: that file is the 50 % result (header `path=1_partitions`, 53,886 bp). The 45 % result is `community_3_stitched_mean_45pct.fa`.
2. The report's "50 % | N/A" identity cell is imprecise (script returns 0.00 %, see above).

---

## 2. `5a22618` pggb per-cluster — cluster_6 REUSABLE, cluster_1 INCOMPLETE

`pggb_analysis/pipeline_results.json` on main contains **only** cluster 6
(`"status": "COMPLETE"`). Verified against the actual files:

- **cluster_6 COMPLETE** — all steps present: `cluster_6.fixed.gfa`, `cluster_6.aln.phy` (8.3 MB), `iqtree_6/cluster_6.treefile`, `cluster_6.dist.tsv`, plus `ancestral/`:
  - `6_ancestral_report.json`: `consensus_length 69153`, `confidence 0.5523`, `n_modules 8`, `n_prophage_seqs 53` — **matches the claimed 69,153 bp / 8 modules / conf 0.5523**.
  - `6_modules.json`: 8 modules (integration, replication, tail_fiber, baseplate, head_capsid, tail_sheath, lysis, lysogeny), total_length 69,153.
  - `6_ancestral.fa`: 69,153 bp + newline; `6_position_confidence.tsv`: 69,153 positions.
- **cluster_1 INCOMPLETE** — graph pipeline ran to completion (`cluster_1.fixed.gfa` 37 MB, "Final graph: 189663 segments, 8235 paths"; `cluster_1.paths.fa` extracted, 577 prophage paths) but **`cluster_1.aln.fa` is 0 bytes** and `pipeline.log` ends at `[mafft 1] / mafft --auto --thread 64 cluster_1.paths.fa`. No `.aln.phy`, no `iqtree_1`, no ancestral. Root cause: MAFFT `--auto` on 577 sequences of ~45–110 kb never completed (intractable). `pipeline_results.json` correctly omits cluster 1.

**Additional finding:** six more ancestral genomes (clusters **3, 5, 7, 8, 9, 11**)
were produced by the pipeline but **never committed** — they exist only as
untracked files in the archived agent-220 worktree and the recovery root
`pggb_analysis/ancestral_genomes/` (identical content, e.g. `3_ancestral.fa`
sha256 `da210d97…` in both). These are small (a few hundred KB total), real
outputs, and will be lost when the recovery dir is deleted.

| cluster | seqs | ancestral len | conf |
|---|---|---|---|
| 3 | 1,234 | 113,502 | 0.484 |
| 5 | 34 | 8,888 | 0.4985 |
| 7 | 31 | 19,183 | 0.6577 |
| 8 | 5 | 51,243 | 0.5684 |
| 9 | 27 | 24,689 | 0.5775 |
| 11 | 65 | 39,154 | 0.546 |

---

## 3. `0c8373c` research-impg-partition — outputs only in archived worktree; DECISION: IMPORT

**What exists (agent-221 worktree `research_outputs/`):**

- **Committed in `0c8373c`** (objects present in fresh main repo; branch `wg/agent-221/research-impg-partition` only): `analyze_partitions.py`, `compare_ancestral.py`, `compare_ancestral_v2.py`, `full_pipeline.py`; `partitions/` = **55 BED + 52 MAF** (3 largest MAFs — partitions 4, 5, 16, 259/174/40 MB — were never generated, they timed out, so they are *not* in the commit and *not* anywhere); `research_output/` = `RESEARCH_REPORT.md`, `ancestral_comparison.json`, `comparison_report.md`, `heaviest_bundle_consensus.fa`, `heaviest_bundle_consensus_maf.fa`, `partition_stats.json`.
- **Untracked (worktree filesystem only)**: `cluster_3.fa` (1,234 seqs), `cluster_3.syng.*` (6-file syncmer index), `cluster_3.alignments.paf` (268 MB) + `.impg`, `ancestral_to_community.paf`, two 100+ MB run logs.
- **Not on main**: nothing from this task is on main.

**Verified content** (from `RESEARCH_REPORT.md` / JSONs): 55 partitions (window 10000, merge 5000), 52 MAF consensus, heaviest-bundle consensus 608,861 bp for community 3, community→ancestral identity 79.10 %, partition-consensus→ancestral 42.76 %.

**Decision: IMPORT, not re-derive.** The durable scientific artifacts are
already in git commit `0c8373c`; re-derivation would repeat ~2 min of
partitioning + 3 MAF timeouts for zero gain. Concretely:

1. Cherry-pick `0c8373c` onto main (or the working feature branch):
   `git cherry-pick 0c8373c` (objects are present locally; do not push).
2. Do **not** import: `cluster_3.alignments.paf` (268 MB, regenerable), the two
   ~100 MB run logs (junk), or the syng index binaries (rebuild in **458 ms**:
   `impg syng`).
3. If the 3 missing MAFs (partitions 4, 5, 16 — the core homologous regions,
   973/817/451 seqs) are needed, they must be **re-derived** with a larger
   timeout / more memory (report: POA for 973 seqs × 2.25 MB used 12+ GB RAM).

---

## 4. Failed attempts — what exists, what is junk

| Agent / task | Branch head | Commits | Working-tree evidence | Verdict |
|---|---|---|---|---|
| agent-220 `run-pggb-ancestral` | `fcff092` (== main) | **none** | `pggb_analysis/pipeline_results.json` → cluster 1 `"FAILED: pggb graph"`; cluster_1 files deleted/modified uncommitted; 11 `seqwish-*` temp dirs; `target/` dir | **FAILED** — no committed output. The only valuable files in its tree are the *uncommitted* ancestral genomes from the earlier `5a22618` run (§2), which should be imported separately. |
| agent-223 `run-impg-partition` | `64ddebf` | `5125aee`, `64ddebf` (branch only, not on main) | committed: 11 community mean genomes + 2,721 per-partition consensus FASTA + 11 BEDs + `mean_genomes_report.md`; `64ddebf` adds `mean_genomes_report.pdf` (54 KB) | **REPAIRABLE / PARTIAL** — real committed outputs, but quality is dubious: e.g. `community_3_mean.fa` = 1,398,681 bp (12.3× the 113,502 bp ancestral, 2.3 % Ns + ambiguous bases) — a crude concatenation, inferior to `2363ece`'s 124,935 bp / 78.05 % stitching result. Reuse the per-partition FASTA consensus files as inputs; treat the assembled "mean genomes" as draft-quality. Community 6 produced no mean genome (⚠ in report). |
| agent-224 `benchmark-find-fastest` | `fcff092` (== main) | **none** | untracked `community_3.fa` (== agent-221 `cluster_3.fa`, sha256 `d6a4e489…`), `.fai`, `community_3_names.txt` | **FAILED** — no committed output; the only files are regenerable inputs. |
| agent-225 `run-all-wave` | `fcff092` (== main) | **none** | empty worktree | **FAILED / JUNK** — nothing to salvage. |

---

## 5. `full_prophages.fa` integrity — REUSABLE (fully verified)

| Check | Result |
|---|---|
| Sequence count | **132,393** ✓ (matches task claim) |
| Header uniqueness | 132,393 unique — **no duplicates** ✓ |
| vs source CSV (`26k_prophage1.csv`, 132,404 prophage_id rows) | all 132,393 fasta headers present in CSV (0 orphans); exactly **11** CSV ids missing |
| vs extraction log (`extraction_log.txt`, committed in `ba841f4`) | log: "132,404 prophage records"; **11 FAIL** lines; the 11 missing ids **exactly match** the 11 FAILs (GCF_001291365.1 ×2, GCF_003886435.1 ×7, GCF_012029685.1 ×2 — "FASTA not found") → 132,404 − 11 = 132,393 ✓ |
| File size / git | 3,265,823,570 B (3.04 GiB ≈ "3.1 GB"); working tree matches committed blob `af0fc973ec51d5b062e26339d5f4de8c74db74f1` (added in `6f96e89` rescue checkpoint of `extract-full-132k`) |
| sha256 | `ed85b2fb549be18bc638d8485f5b5add7c2d394f3822efe66a90ca6d979758d3` |
| Extraction script | `extract_prophages.py` on main == `ba841f4` version (sha256 `7af1da4e…`) ✓ |

**Verdict:** the file is intact and exactly consistent with the documented
extraction (11 genuinely-missing source genomes, all accounted for).

⚠ **Storage note:** `full_prophages.fa` (3.11 GiB blob) plus
`cluster_1.alignments.paf` (647 MiB) and `cluster_1.mappings.paf` (75.8 MiB)
are **tracked in git** on main (see `manifests/git-blobs-over-50MiB.tsv`).
This is why GitHub refuses the push (per `PHIND_RECOVERY_HANDOFF.md`). An
LFS/migration decision is out of scope here, but downstream tasks should know
the repo cannot be pushed until this is resolved.

---

## 6. Commits on fresh main vs only in archived worktrees

**On fresh main (`a5070e7`, all 20 commits ahead of origin):**
`6d2c822`, `d6f666e`, `ba841f4` (extract-full-132k), `6f96e89`, `61ac2b4`,
`c6ce454` (full-prophage-homology), `5a22618` (pggb-per-cluster), `063a4dd`,
`51a7ea8`, `b012524` (full-132k-prophage), `ec9a754`, `4d6328e`,
`c52553c` (focused-heatmap-connected), `062bae5`, `fcff092`, `f7abf5b`,
`9f9a440`, `2363ece` (build-partition-stitching), `2f991af`, `a5070e7`.

**Only in archived worktrees (NOT on main):**
- `wg/agent-221/research-impg-partition` → **`0c8373c`** (impg partition research, community 3)
- `wg/agent-223/run-impg-partition` → **`5125aee`** (run impg partition + mean genome on 12 communities), **`64ddebf`** (mean genomes report PDF)

**No commits at all (branch head == main commit):**
`wg/agent-220/run-pggb-ancestral`, `wg/agent-224/benchmark-find-fastest`,
`wg/agent-225/run-all-wave`. (`wg/agent-226/build-partition-stitching`'s commit
`2363ece` IS on main.)

All commit objects for `0c8373c`/`5125aee`/`64ddebf` are present in the fresh
repo's object store, so cherry-picking from the archived refs works without
the bundle.

---

## 7. Downstream notes

- **`stitching-validation`**: stitch_algorithm.py + inputs on main, but the
  502 partition MAFs + ancestral FASTA are only in the recovery dir — copy them
  to main (or re-derive from `full_prophages.fa` + `impg`/`wfmash`) before the
  recovery dir is removed. The community-3 ancestral (113,502 bp, sha256
  `da210d97…`) and `heaviest_bundle_consensus.fa` are available there.
- **`ecor-mapping`**: full_prophages.fa is verified intact (132,393 seqs);
  header naming is `>GCF_<acc>_prophage_<n>`. Community assignment for every
  sequence is on main in `prophage_homology_survey/full_heatmap_clusters.csv`
  (sequence,community; community 3 = 1,234 seqs; communities 0/1/10 = 4,999/
  2,128/5,208; long tail of singletons — filter when mapping).
- **`all-prophage-mash`**: full_prophages.fa is the right input (verified);
  expect 132,393 leaves, 11 source genomes absent from the original 26k cohort.
- Import candidates that will be lost with the recovery dir: uncommitted
  ancestral genomes (clusters 3,5,7,8,9,11), agent-221 research outputs,
  agent-223 mean genomes (draft quality).
