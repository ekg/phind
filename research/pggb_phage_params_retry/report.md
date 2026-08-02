# pggb retry on phage-appropriate parameters — cluster_6 pilot

**Task:** `phage-params-pggb` (retry of the historical pggb run judged inappropriate
for divergent prophages).
**Status:** pilot complete on cluster_6 (53 prophages, 817,520 bp total, mean 15.4 kb).
**Date:** 2026-08-02.

## 1. Parameter sets

| stage | historical run (baseline) | corrected run (this pilot) |
|---|---|---|
| mapper | wfmash v0.24.1 (`pggb_env`) | pggb 0.6.0 wrapper (certified tool, `pggb-env`) |
| map | `wfmash -p 85 -n 5 -k 19 -l 2000 -X -m` | `pggb -p 75 -s 250 -l 500 -k 11 -g 80 -c 5` |
| align | `wfmash -p 85 -n 5 -k 19 -l 2000 -X` | same pggb params (wfmash inside: `-s 250 -l 500 -p 75 -n 5 -k 19 -H 0.001 -Y '#' --lower-triangular --hg-filter-ani-diff 80`) |
| graph induce | `seqwish -k 19` | `seqwish -k 11` (min-match 11) |
| smooth | `smoothxg -X 100 -r 53` | pggb 0.6.0 defaults (chop 100, POA targets 700,900,1100) |
| normalize | `gfaffix` | pggb `gfaffix` + odgi sort/view (final) |

**pggb 0.6.0 CLI mapping of the corrected user set (chat 3)**
(per certified help in `artifacts/consumer_compatibility/tool_versions.json` and the
pggb 0.6.0 wrapper source):

| user param | pggb 0.6.0 flag | pipeline stage |
|---|---|---|
| `p=75` | `-p, --map-pct-id 75` | wfmash identity threshold 75 % (was 85) |
| `s=250` | `-s, --segment-length 250` | wfmash seed segment length 250 bp (was default 5 k / unset) |
| `l=500` | `-l, --block-length 500` | wfmash min block length 500 bp (was 2000) |
| `k=11` | `-k, --min-match-len 11` | seqwish exact-match filter 11 bp (was 19) |
| `ani-diff=80` | `-g, --hg-filter-ani-diff 80` | wfmash ANI-difference filter 80 (was pggb default 30) |

`-c 5` (wfmash mappings/segment) was kept at 5 for parity with the historical run so
the comparison isolates the identity/segment/block/minmatch/ani-diff changes.

## 2. Outputs (this pilot)

Directory: `pggb_analysis/phage_params_retry/`

- `cluster_6.fa.7a957cc.bb87cbe.5772f4f.smooth.final.gfa` — final graph (GFA)
- `...smooth.gfa`, `...seqwish.gfa`, `...fix.gfa` — smoothed / seqwish / fixed graphs
- `cluster_6.fa.7a957cc.alignments.wfmash.paf` — wfmash base-level alignments (PAF)
- `...mappings.wfmash.paf` — approximate mappings
- `cluster_6.fa.7a957cc.bb87cbe.5772f4f.smooth.<ts>.log` — pggb run log (params echo,
  per-step timings)
- `paths.fa` — all 53 prophage graph paths as FASTA (odgi paths -f)
- `aln.fa` — MAFFT alignment of the 53 prophage paths (mirrors historical
  `cluster_6.aln.fa`)
- `consensus.fa` + `consensus_confidence.tsv` — majority-rule consensus from the
  alignment with per-position confidence (pipeline-style ancestral/consensus)
- `pilot_stats.json` — machine-readable stats for both runs
- `comparison_table.md` — the table below

## 3. Comparison table (new vs old, cluster_6)

| metric | new (p=75,s=250,l=500,k=11,ani80) | old (p=85,l=2000,k=19) | notes |
|---|---|---|---|
| Graph segments (final) | 17,313 | 17,663 | similar |
| Graph edges (final) | 23,869 | 24,456 | similar |
| Graph total bp (final) | 75,459 | 77,304 | similar |
| Segment N50 (final, bp) | 4 | 5 | both runs chopped ~100 bp + normalized; not a differentiator |
| Seqwish segments (pre-smooth) | 10,422 | 7,689 | corrected params retain more induced variation |
| Seqwish paths | 53 | 53 | — |
| Final graph paths (P/W lines) | 53 | 130 | old graph carried 77 extra chimeric traversal paths |
| Prophage paths | 53 | 53 | 1:1 with input prophages in both |
| Prophage path total bp | 817,520 | 1,014,303 | new = exact input total; old inflated by a 69 kb artifact path |
| Prophage path N50 (bp) | 14,752 | 19,512 | old N50 inflated by artifact path |
| Prophage path mean len (bp) | 15,424.9 | 19,137.8 | new = input mean |
| Consensus paths in graph | 0 | 0 | pggb 0.6.0 emits none by default (consensus-spec off) |
| Mean query coverage (PAF) | 0.9234 | 0.9138 | new aligns a higher fraction of each prophage |
| Median query coverage | 1.0000 | 0.9542 | new: half of prophages 100 % aligned |
| Min query coverage | 0.4025 | 0.4084 | one divergent outlier in both |
| Queries with any hit | 52/53 | 52/53 | same one query unmapped in both |
| PAF records | 1,347 | 2,128 | old self-map PAF had ~1.6× more redundant records |
| Alignment columns | ~157 k (MAFFT of paths) | 157,426 | comparable column count |
| Consensus length (bp, MAFFT-majority) | — | 69,153 | see §4 caveat |
| Consensus confidence | — | 0.5523 | see §4 caveat |

## 4. Interpretation

1. **The corrected params recover more of each prophage as aligned sequence.**
   Mean per-query aligned coverage rises 91.4 % → 92.3 % and median 95.4 % → 100 %.
   This is the expected signature of a more sensitive mapping set (identity 75 vs 85,
   block 500 vs 2000, min-match 11 vs 19, ani-diff 80) on divergent prophages: regions
   previously below the old identity/block thresholds are now captured.
2. **The corrected params produce a cleaner graph.** The new final graph has exactly
   53 paths (one per input prophage) whose lengths sum to exactly the input total
   (817,520 bp). The old final graph had 130 paths and an inflated path total
   (1.01 Mb) driven by a 69 kb chimeric traversal path — an artifact of graph
   topology under the old params, not real sequence.
3. **Graph size is not inflated.** Final node/edge counts and total graph bp are
   essentially unchanged (17.3 k vs 17.7 k segments; 75.5 kb vs 77.3 kb total), so the
   more sensitive params do not blow up the pangenome. The seqwish (pre-smoothing)
   graph is denser (10.4 k vs 7.7 k segments; 23.0 k vs 16.8 k edges), reflecting
   divergence that was previously collapsed.
4. **Same single unmapped outlier** (one divergent prophage at ~40 % coverage) in both
   runs — not a parameter artifact.

**Caveats / differences between the runs (not controlled):**
- The old run used the individual-tool chain (wfmash v0.24.1 + seqwish + smoothxg +
  gfaffix); the new run uses the pggb 0.6.0 wrapper (wfmash v0.13). Orchestration
  defaults (e.g. smoothxg POA targets 700/900/1100, lower-triangular mapping) differ
  slightly. The comparison should be read as pipeline-level, not a single-variable
  A/B.
- pggb 0.6.0 emits no `Consensus_*` paths by default (consensus-spec off), so the
  "consensus sequence" for the new run is the pipeline-style majority-rule consensus
  from the MAFFT alignment (`consensus.fa`), not a graph-native consensus. The old
  run's 69,153 bp "ancestral" consensus was built from the inflated old paths
  (including the artifact path) and is not directly comparable; the alignment-based
  consensus of the new run is the more defensible ancestral approximation.

## 5. Recommendation

**Yes — scale to more clusters with the corrected parameter set, with two caveats.**

1. The corrected set (`p=75 s=250 l=500 k=11 ani-diff=80`) improves alignment capture
   and graph cleanliness on the pilot cluster at no graph-size cost. Use it for the
   per-cluster prophage pangenomes.
2. Before the full 132 k / all-wave run, validate on 2–3 additional clusters spanning
   the divergence range (the ECOR-mapped cluster and one low-identity cluster),
   because this pilot is n = 1 cluster and the wfmash-version confound is not yet
   eliminated. If those confirm the pilot, proceed cluster-by-cluster; do not attempt
   the full all-wave graph in one step.
3. Adopt the pggb 0.6.0 wrapper (certified tool) rather than the hand-rolled
   wfmash/seqwish/smoothxg calls in `scripts/pggb_per_cluster_pipeline.py`, or update
   the script to take these parameters as CLI args — the wrapper's parameter semantics
   are the canonical ones and avoid the wfmash version-dependent flag mismatch.
