# NTM v2 functional-QC pipeline (Pharokka + CheckV)

Per-genome gene annotation, module detection and CheckV quality control for the
**current NTM v2 generated genomes — 2,388 ML + 1,251 ancestral = 3,639
records** from `ntm/v2/release/` (task `build-ntm-v2`; executed by
`run-ntm-v2`). Adapts the established v1 E. coli/NTM workflow
(`research/phage_annotation/`) with pinned tools, a separate v2 analysis root,
restart safety, and a graded candidate-functionality classification.

## Scientific framing

- **ML genomes are the primary target**; ancestral genomes are an explicit
  **comparator**. All summaries are broken out by cohort:
  - `ntm2_ml_reconstructed` — release `status=ml` (1,251 clade ML genomes)
  - `ntm2_ml_singleton` — release `status=singleton` (1,137 short pass-through
    fragments; expected to dominate low tiers and are reported separately so
    they don't drown the reconstructed-clade signal)
  - `ntm2_anc` — 1,251 ancestral-state genomes
- **Graded tiers, never binary function claims** (below).
- **Missing PHROG/mmseqs2 hits are "not detected", not proof of absence.**

## Tools (pinned, read-only shared DBs from the v1 evidence tree)

| Tool | Version | Binary | Database |
|---|---|---|---|
| Pharokka | v1.10.1 | `/home/erikg/micromamba/envs/pharokka/bin/pharokka` | `annotation/pharokka_db` (PHROG, 9 Sep 2025) |
| CheckV | v1.1.1 | `/home/erikg/micromamba/envs/phage-annot/bin/checkv` | `annotation/checkv_db/checkv-db-v1.5` |

(paths under `/mnt/nvme3n1/erikg/phind-genome-work/`). Threads capped at
**64** (`--threads`, hard cap). Commands — identical fast configuration to v1:

```bash
pharokka run -m --mmseqs2_only --skip_extra_annotations --skip_mash \
  -g prodigal-gv -i ROOT/input/all_v2_phage_genomes.fa -o ROOT/pharokka_out \
  -d annotation/pharokka_db -t 64 --locustag NTMV2
checkv end_to_end ROOT/input/all_v2_phage_genomes.fa ROOT/checkv_out \
  -d annotation/checkv_db/checkv-db-v1.5 -t 64
```

## Driver: `ntm/v2/scripts/run_annotation_v2.py`

Analysis root (default, external): `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/annotation`.
The v1 tree (`.../phind-genome-work/annotation`) is **only read** (databases);
the driver refuses any root inside it.

- **`--dry-run`** prepares and validates inputs, prints the pinned commands and
  writes `run_state/tool_provenance.json` (versions, DB identity, release
  FASTA SHA-256s) without running the tools.
- **Input preparation** streams both release FASTAs into
  `input/all_v2_phage_genomes.fa` with **bare, exact genome IDs** (first header
  token only, e.g. `ntm2_0_0000_ML`, `ntm2_0_0000_ANCESTRAL`) plus
  `input/genome_index.tsv` (`genome_id, source ∈ {ntm2_ml, ntm2_anc}, cohort,
  clade_id, length, status, host_clades, origin_fa`). The stage **aborts** on
  any duplicate ID, ID-set mismatch vs the release, length mismatch, or a
  count ≠ 2,388 ML + 1,251 ancestral (round-trip is exact and duplicate-free).
- **Restart-safe stages** `prepare → pharokka → checkv → report`: each writes
  `run_state/<stage>.done.json` (command, versions, DBs, threads, timestamp,
  exit code); a stage is skipped when its marker exists and its outputs are
  present (`--force` redoes, `--stages a,b` selects).
- Provenance: `run_state/tool_provenance.json`; release FASTA SHA-256s, tool
  versions, DB file inventories (name/size/mtime).

Validated (dry run, this machine): 3,639 unique records = 2,388 ML + 1,251
ancestral, round-trip exact, duplicate-free; v1 evidence tree untouched.

## Report: `ntm/v2/scripts/build_annotation_report_v2.py`

Merges Pharokka per-CDS (`pharokka_cds_final_merged_output.tsv`,
`pharokka_length_gc_cds_density.tsv`) + CheckV (`quality_summary.tsv`) into
`ROOT/report/`:

- `per_genome_functional_qc.tsv` — one row per genome (sorted, deterministic;
  reruns are byte-identical): gene/truncation/hypothetical counts, the six key
  modules, PHROG category coverage, temperate signal, CheckV
  grade/completeness/contamination/viral/host genes, hygiene flags, and the
  `candidate_functionality` tier + `tier_rationale`.
- `summary_by_cohort.tsv` / `summary_by_source.tsv` — cohort- and source-level
  aggregates (medians, %-module-detection, tier percentages).
- `candidate_functionality_by_cohort.tsv` — tier × cohort matrix.

Key-module detection matches v1 (PHROG `annot` regexes + category mapping:
`tail`, `integration and excision`→integrase, `lysis`), with one targeted fix:
the v1 lysis regex's bare `R` alternative matched any capital R (e.g. "RNA
polymerase"); it now requires a word-bounded `R`/`Rz` (spanin subunits).

### Graded candidate-functionality tiers (precedence top-down)

| Tier | Criteria |
|---|---|
| `A_strong_candidate` | CheckV Complete/High/Medium **and** coherent head-packaging (≥2 of terminase/portal/capsid) **and** tail **and** lysis **and** truncated_frac ≤ 0.20 **and** contamination = 0 |
| `B_moderate_candidate` | CheckV Complete/High/Medium **and** ≥1 head-packaging module **and** (tail or lysis) **and** truncated_frac ≤ 0.35 **and** contamination = 0 |
| `C_partial_module_evidence` | ≥1 key module detected but tier-A/B criteria unmet (e.g. CheckV Low/Not-determined, head-packaging gap, contamination >0, heavy truncation). `tier_rationale` names the binding constraint. |
| `D_modules_not_detected` | No key-module annotation detected — a **detection** result under mmseqs2-only PHROG matching (novel/divergent ORFans, short fragments, annotation depth), **not proof of biological absence** |
| `no_annotation` | Pharokka called 0 CDS — identifier/tool failure to investigate, never a biological claim |

**Integrase is supportive, not required**: recorded (`integrase`,
`temperate_signal`) since these genomes are prophage-derived, but deliberately
excluded from the tier ladder — temperate signal is not a universal
requirement of a functional phage genome.

Each row carries a human/machine-readable `tier_rationale`
(e.g. `checkv=High-quality; head=3/3[terminase+portal+capsid]; tail=1;
lysis=1; integrase=1(support-only); trunc_frac=0.050; contam=0; all tier-A
criteria met`), so no classification is a black box and annotation failure is
never read as biological absence.

## Tests

`ntm/v2/scripts/test_annotation_v2.py` (pytest, 25 tests): key-module
detection (text + category + the RNA-polymerase-is-not-lysis regression),
every tier boundary (contamination, truncation, CheckV grade, integrase
optionality), cohort derivation, summary/matrix correctness, byte-determinism
of reruns, CheckV-coverage abort semantics, driver round-trip on a synthetic
release (duplicates/counts abort), thread-cap enforcement, v1-tree write
guard, and the real-release round-trip (3,639 records) when the release is
mounted.

```bash
python3 -m pytest ntm/v2/scripts/test_annotation_v2.py -v
```

## Caveats (inherited from v1)

1. Annotation is **mmseqs2-only** (sequence similarity); remote/novel ORFans
   stay "hypothetical". The high hypothetical fraction expected for NTM is
   partly real (mycobacteriophage ORFans) and partly annotation depth — a
   structural pass (phold on GPU) would reduce it.
2. CheckV completeness is AAI/reference-based — conservative for reconstructed
   genomes and under-calls novelty (Not-determined).
3. `D_modules_not_detected` and `missing_core_structural` are detection
   statements, not absence claims.
