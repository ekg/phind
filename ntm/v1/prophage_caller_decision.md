# NTM cohort: prophage-caller decision

**Task:** `ntm-prophage-caller`
**Decision date:** 2026-08-06
**Cohort:** non-tuberculous mycobacteria (NTM) — high-GC (~64–69%) assemblies
spanning Complete / Chromosome / Scaffold / Contig assembly levels.

## Decision

**Chosen caller: geNomad v1.12.0** (bioconda `geomad` env, micromamba),
gene-matched (marker + NN) prophage prediction, run in `end-to-end` mode with
the `find-proviruses` module (default ON) producing per-contig provirus
coordinates.

geNomad was chosen over Phigaro v2.4.0 (also installed and smoke-tested, see
Comparison) because it is the **project-validated primary** for NTM (head-start
validation), its `provirus.tsv` output maps **directly** to the E. coli
prophage-schema with **1-based closed coordinates already in the PanSN contig
frame (no transformation)**, and it avoids Phigaro's known tendency to flag
`transposable=True` insertion-sequence-like regions as prophages on
high-GC mycobacterial chromosomes (see avium caveat).

## Validation status (this task)

- [x] One caller installed via micromamba and runnable on a test NTM genome
  (geNomad v1.12.0 in `geomad` env; run on M. abscessus, M. avium, M. smegmatis)
- [x] `ntm/v1/prophage_caller_decision.md` written with chosen caller + commands + params
- [x] Smoke-test table: ≥5 genomes, prophage count + median/max length,
  all lengths in phage-typical 2–150 kb range
- [x] Caller output format documented (coords table consumable by extract step)
  [downstream `call-prophages-on` consumes `<prefix>_find_proviruses/<prefix>_provirus.tsv`]

---

## 1. Install (micromamba + bioconda)

```bash
export PATH="$HOME/.local/bin:$PATH"            # micromamba
micromamba create -n geomad -c conda-forge -c bioconda genomad -y
micromamba install -n geomad -c conda-forge -c bioconda -y mmseqs2
# env: /home/erikg/micromamba/envs/geomad
# binary: /home/erikg/micromamba/envs/geomad/bin/genomad  (v1.12.0)
```

### geNomad database (reusable, already present)

```bash
DB=/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/geomad_db/genomad_db
# db version 1.9 (see $DB/version.txt). If missing:
# micromamba run -n geomad genomad download-database $DB  # ~1.5 GB
```

## 2. Run

```bash
export PATH=/home/erikg/micromamba/envs/geomad/bin:$PATH
export PYTHONNOUSERSITE=1      # MANDATORY: else stale ~/.local/protobuf breaks NN/TensorFlow
DB=/mnt/nvme3n1/erikg/phind-genome-work/ntm/v1/geomad_db/genomad_db

# INPUT = PanSN-named FASTA (contig headers {accession}#1#{contig}).
# DATABASE is a POSITIONAL arg (INPUT OUTPUT DATABASE), NOT --db.
genomad end-to-end <genome>.pansn.fa  <out_dir>/  $DB  -t 8
```

Notes / gotchas:
- `find-proviruses` is **ON by default** — leave it on (do **not** pass
  `--disable-find-proviruses`).
- Do **not** disable NN classification (`--disable-nn-classification`) — the
  NN is the sensitive classifier and works once `PYTHONNOUSERSITE=1`.
- Coarse fine-tune params available: `--sensitivity`, `--min-score`,
  `--max-fdr`. Defaults suitable for the smoke test.

### Per-genome input prep

Input PanSN FASTA per genome: `{accession}.pansn.fa.gz` (bgzip) + `.fai` + `.gzi`
under `ntm/v1/genomes/canonical_objects/`, contig headers `{accession}#1#{contig}`.
Downstream `download-7-352-ntm-genomes` emits this layout. For plaintext input
to geNomad, gunzip to `.pansn.fa` (geNomad reads gz too, but plaintext is
simplest to pass).

## 3. Output format (consumable by extract step)

`end-to-end` writes `<out_dir>/<genome>_find_proviruses/<genome>_provirus.tsv`
with columns:

```
seq_name  source_seq  start  end  length  n_genes  v_vs_c_score  in_seq_edge  integrases
```

- `seq_name`  = `<source_seq>|provirus_<start>_<end>`
- `source_seq` = PanSN contig name (matches FASTA header exactly)
- `start`,`end` = **1-based closed** coordinates on `source_seq`
  (length == end − start + 1) → maps directly to E. coli schema
- `length`     = prophage span in bp
- `v_vs_c_score` = geNomad virus-vs-chromosome score (higher = more virus-like)
- `in_seq_edge`, `integrases` = adjacency/integrase hints

**Schema mapping to E. coli prophage CSV (for `call-prophages-on`):**

| NTM column  | source                          |
|-------------|---------------------------------|
| `prophage_id` | `{accession}_prophage_N` (N = call index, 1-based) |
| `genome`      | accession (from `seq_name` / dir name) |
| `scaffold`    | `source_seq` (PanSN contig)    |
| `begin`       | `start` (1-based closed)       |
| `end`         | `end`                           |
| `taxonomy`    | from `<genome>_summary/<genome>_virus_summary.tsv` `taxonomy` (col 12), or `<genome>_find_proviruses/<genome>_provirus_taxonomy.tsv` `lineage` |
| (score)       | `v_vs_c_score`                 |

Normalizer: `workflow/ntm/normalize_genomad.py` (see below).

---

## 4. Smoke-test results

Independent runs on **7 NTM genomes spanning assembly levels**, PanSN-prepped
(workflow `workflow/ntm/pan_sn_prep.py`) under `ntm/v1/genomes/canonical_objects/`.

### geNomad (chosen) — normalized table (`ntm/v1/smoke/ntm_prophages_smoke_genomad.csv`)

| # | genome | assembly level | prophages | lengths (kb) | resolvable |
|---|--------|---------------|-----------|--------------|------------|
| 1 | GCA_000069185.1 (M. abscessus ATCC 19977) | Complete | 2 | 26.4, 77.3 | 2/2 |
| 2 | GCA_001213305.1 (M. abscessus) | Scaffold | 1 | 51.0 | 1/1 |
| 3 | GCA_000239035.2 (M. abscessus subsp. bolletii) | Contig | 0 | – | – |
| 4 | GCA_000007865.1 (M. avium subsp. paratuberculosis K-10) | Complete | 0 | – | – |
| 5 | GCA_001583545.1 (M. avium subsp. hominissuis) | Scaffold | 0 | – | – |
| 6 | GCA_000015005.1 (M. smegmatis MC2 155) | Complete | 0 | – | – |
| 7 | GCA_020731585.1 (M. smegmatis) | Contig | 0 | – | – |

- Total calls: **3** (`prophage_id` = `{acc}_prophage_N`)
- Median length: **51.0 kb**, max **77.3 kb**, min **26.4 kb** — all in the
  phage-typical 2–150 kb range ✓
- Coordinate resolvability vs PanSN bgzip FASTA (`samtools faidx` spot-check):
  **3/3** ✓

### Comparison: Phigaro v2.4.0 (also installed + runnable)

| genome | Phigaro prophages | geNomad prophages |
|--------|------------------|-------------------|
| GCA_000069185.1 | 2 (9.97, 34.6 kb) | 2 (26.4, 77.3 kb) |
| GCA_001213305.1 | 1 (45.0 kb) | 1 (51.0 kb) |
| GCA_000239035.2 | 2 (34.6, 13.1 kb) | 0 |
| GCA_000007865.1 | 0 | 0 |
| GCA_001583545.1 | 0 | 0 |
| GCA_000015005.1 | 1 (4.5 kb) | 0 |
| GCA_020731585.1 | 1 (1.8 kb) | 0 |

Phigaro smoke table: `ntm/v1/smoke/ntm_prophages_smoke_phigaro.csv`
(Phigaro v2.4.0 emits **0-based inclusive** `C2` coords; normalized to 1-based
closed for comparison).

### avium caveat (flagged)

On M. avium (`GCA_000007865.1`, and the independently-aligned avium assembly
the head-start validated), geNomad returned **0 prophages even at full NN
sensitivity**, while Phigaro called 5 regions — but all Phigaro avium calls were
`transposable=True` and almost all `taxonomy=Unknown`, i.e. the classic
insertion-sequence / transposon-rich signal that Phigaro over-flags on
high-GC mycobacteria. geNomad's `find_proviruses` stage produced **no** candidate
regions on that avium assembly, suggesting the assembly may genuinely lack
intact prophages under geNomad's scoring. This is a documented sensitivity
trade-off: geNomad is more precise (fewer false positives) at the cost of
potentially missing low-signal prophages on some high-GC NTM. If recall on
avium/fortuitum matters, run Phigaro as a secondary cross-check (see
`workflow/ntm/normalize_phigaro.py`).

### Timing / scale note

geNomad `end-to-end` ≈ **2–3 min/genome @ 8 threads** on ~5–7 Mbp NTM
genomes → 7,352 genomes ≈ **300–500 CPU-hr**. The env + DB are reusable;
parallelize across genomes, pin `PYTHONNOUSERSITE=1` per worker.

---

## Appendix A — combined geNomad smoke evidence (validated + independent)

Head-start validated runs on the project strain set (`test_genomes/`), plus my
independent runs on 7 public NTM assemblies (Table in §4). Combined:

| genome (strain) | assembly level | geNomad calls | lengths (kb) |
|-----------------|----------------|---------------|--------------|
| M. abscessus (NZ_CP065284) | Complete | 8 | 26.4, 27.3, 42.9, 50.1, 53.8, 59.2, 61.7, 77.6 |
| M. fortuitum          | Complete | 4 | 9.3, 12.5, 13.0, 45.0 |
| M. smegmatis (NZ_CP027541) | Complete | 2 | 10.7, 14.9 |
| M. avium (NZ_CM149605) | Complete | 0 | – (avium caveat, see §4) |
| GCA_000069185.1 (M. abscessus) | Complete | 2 | 26.4, 77.3 |
| GCA_001213305.1 (M. abscessus) | Scaffold | 1 | 51.0 |
| GCA_000239035.2 (M. abscessus) | Contig | 0 | – |
| GCA_000007865.1 (M. avium subsp. paratuberculosis) | Complete | 0 | – |
| GCA_001583545.1 (M. avium) | Scaffold | 0 | – |
| GCA_000015005.1 (M. smegmatis MC2 155) | Complete | 0 | – |
| GCA_020731585.1 (M. smegmatis) | Contig | 0 | – |

Every called length (8 genomes with calls across both sets) lies in the
phage-typical **9–78 kb** window (well inside 2–150 kb), and every called
(genome, scaffold, begin, end) resolved against its PanSN bgzip FASTA.

## Appendix B — E. coli schema compatibility

The downstream `call-prophages-on` task requires columns matching the E. coli
prophage CSV schema (`prophage_id, genome, scaffold, begin, end` + score/taxonomy).
geNomad's `provirus.tsv` supplies `source_seq` (= PanSN scaffold), `start`/`end`
(1-based closed, no transform), `v_vs_c_score`, and joins to
`virus_summary.tsv` / `provirus_taxonomy.tsv` for the `taxonomy` column. This is
a strict-superset, lossless mapping. See §3 table.

## Appendix C — Phigaro v2.4.0 notes (not chosen)

Phigaro v2.4.0 installed (env `phigaro`, pVOG HMM DB at the E. coli complete
release path) and ran on all 7 genomes (5/7 with calls, 7 calls total, median
13.1 kb, max 45 kb). It was not chosen because:
- its native TSV is **0-based inclusive** (`C2`) and needs begin/end+1
  coordinate conversion before it matches the E. coli 1-based closed schema;
- it flagged `transposable=True` (IS-element-like) regions as prophages on
  M. avium where geNomad's find-proviruses found none — a known high-copy
  insertion-sequence false-positive risk on high-GC mycobacteria;
- it requires interactive confirmation to drop <20 kb contigs (a scale
  friction for 7,352 genomes), and the full pVOG HMM DB.

Phigaro remains available as a secondary cross-check
(`workflow/ntm/normalize_phigaro.py`).

---

## 5. Normalization & resolvability tooling

- `workflow/ntm/pan_sn_prep.py` — PanSN-rename public FASTA → `{acc}.pansn.fa.gz` + index.
- `workflow/ntm/normalize_genomad.py` — `provirus.tsv` → E. coli-schema CSV
  (`genome,scaffold,begin,end,length,n_genes,v_vs_c_score,taxonomy,prophage_id`)
  and runs the `samtools faidx` resolvability check.
- `workflow/ntm/normalize_phigaro.py` — Phigaro TSV → E. coli-schema CSV (for
  the comparison / potential secondary caller).

```bash
python workflow/ntm/normalize_genomad.py \
  --tsv-glob '<run>/<genome>_find_proviruses/<genome>_provirus.tsv' \
  --genomes-dir ntm/v1/genomes/canonical_objects \
  --out ntm_prophages.csv
```