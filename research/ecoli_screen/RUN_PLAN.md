# RUN_PLAN — Local E. coli Logan-subset screen (stratified Tier E1+E2) via S3 assemblies

**Task:** `local-e-coli` · **Preregistered:** 2026-08-22T00:25:00Z (UTC)
**Toolchain:** reuse of the validated `local-ntm-logan` pipeline
(`ntm/v2/external_validation/local_screen/*.py`, minimap2 2.31 +
back_to_sequences 0.8.4, proven at 17.4k accessions / 97/97 confirmations)
with **path/generalization diffs only** (documented §9); the screen, gate,
confirmation and ladder logic are unchanged. Zero Logan Search spend.

Everything below is fixed **before any assembly download**. The metadata
sweep (§4) is design-time input gathering required by the task ("exact
BioProject accessions chosen from RunInfo metadata at design time") and
downloads **no assemblies**. Amendments require a logged entry in
`AMENDMENTS.md` (none at preregistration time).

## 1. Objective

Screen the frozen 113-bait E. coli prophage panel
(`research/ecoli_bait/panel/baits.fa`, sha256 in §2) against Logan's freely
downloadable per-accession assemblies (S3 `logan-pub`, contigs) for a
**preregistered stratified subset** of E. coli SRA WGS runs, then calibrate
the method against the NTM run (`ntm/v2/external_validation/local_screen/REPORT.md`):

- **Tier E1 — depth tier**: all WGS runs from 2–4 high-value
  studies/populations (ST131 + global-diversity + longitudinal archetypes;
  exact BioProjects chosen mechanically from RunInfo metadata per §5, frozen
  with justification + sha256). Target ~8–12k runs.
- **Tier E2 — breadth tier**: stratified random sample across all remaining
  E. coli WGS runs, stratified by BioProject with a per-project cap (NTM
  showed ERP001039-type clone-complex collapse). Target ~10–15k runs.

## 2. Inputs (frozen — no redesign)

| Input | Path | Integrity |
|---|---|---|
| Bait panel (113 rows: 24 interior_module, 59 junction, 24 member_interior, 2 poscon, 2 host-ctrl, 2 shuffled-neg) | `research/ecoli_bait/panel/baits.fa` + `bait_manifest.tsv` | per-bait `seq_sha256` in manifest; panel provenance `research/ecoli_bait/panel/provenance.json` (run id `2c0055272b71cc48`, deterministic, byte-identical reruns) |
| Positive-control references | `research/ecoli_bait/panel/public_refs/NC_001416.1.fa` (λ, 48,502 bp), `NC_002167.1.fa` (HK97, 39,732 bp) | sha256 recorded in panel provenance; poscon baits are exact mid-genome slices of these |
| Toolchain | `ntm/v2/external_validation/local_screen/*.py` | generalized copies under `scripts/` (§9 diff list) |
| Thresholds | 0.7 primary screen threshold; minimap2 ladder `sr → asm5 → asm20`; no threshold fishing | mirrors `logan/PILOT_PLAN.md` §0,§2,§4 as executed in the NTM run |

## 3. Universe query (calibrated 2026-08-22T00:18–00:24Z, frozen)

```
"Escherichia coli"[Organism] AND wgs[strategy]      → 538,919 experiments
```

Calibration notes (receipts in `logs/eutils_calibration.log`):
`strategy wgs` as bare tokens misparses (1,050); `AND strategy:wgs[sb]`
(538,919) and `AND wgs[strategy]` (538,919) agree and match the task's
"≈ 539k runs". The frozen form is `wgs[strategy]` (flat query; nested-paren
defect from the NTM run avoided). Uids retrieved by `usehistory` +
`retstart` pagination at `retmax=100000` (deep paging verified at
retstart=500000 before preregistration).

## 4. Metadata sweep (design-time; no assemblies downloaded)

1. `esearch` (query §3, usehistory, retmax 100k, retstart paging) →
   `metadata/ecoli_universe_uids.txt` (+ sha256 of the frozen uid list).
2. `esummary` (db=sra, retmode=json, chunks of 400 uids, ≤ 1 req/s,
   3 retries / 5-20-60 s backoff — identical protocol to NTM Amendment A2)
   → `metadata/ecoli_universe_runinfo.tsv`: uid, runs, scientific_name,
   library_strategy, library_source, platform, bioproject, n_runs.
   Every response sha256-logged (`logs/esummary.log`).
3. **Frozen filters** applied to the sweep before any tier design:
   - keep `library_strategy == "WGS"` exactly (server query is
     strategy-side; belt-and-braces);
   - keep `scientific_name` startswith `Escherichia coli` (E. coli incl.
     O157:H7 etc.; drops e.g. mixed/metagenome/synthetic naming);
   - one row per RUN accession (runs deduped; multi-run experiments expand;
     per-run tier attribution kept).
4. `scripts/build_bioproject_table.py` → per-BioProject WGS-run counts +
   name composition (`metadata/bioproject_table.tsv`), plus db=bioproject
   esummary titles for all projects ≥ 100 frozen WGS runs (≤ 1 req/s,
   sha256-logged → `metadata/bioproject_titles.tsv`).

## 5. Tier E1 rule (mechanical, frozen before any download)

E1 candidates are BioProjects (≥ 100 frozen WGS runs, title known) matching
exactly one archetype by **case-insensitive title/attribute substring**:

| archetype | match strings (any) |
|---|---|
| ST131 | `st131` |
| global-diversity | `global`, `worldwide`, `international`, `diversity` |
| longitudinal | `longitudinal`, `timecourse`, `time course`, `time-series`, `persistence`, `follow-up` |

Selection: within each archetype, rank projects by frozen WGS run count
(desc, accession asc); greedily take whole projects in archetype order
(ST131 → global-diversity → longitudinal) while the cumulative union stays
≤ **12,000** runs (a project that would exceed the cap is skipped, not
truncated). Accept union in the **8,000–12,000** window; if the matching
union is < 8,000, take it whole and document. Tier E1 = this exact run set;
manifest frozen (query strings, per-project justification, run list,
sha256) in `manifests/tierE1_manifest.json` **before any download**.
Sampling seed (for any residual stochasticity): **20260822**.

## 6. Tier E2 rule (stratified random sample, frozen before any download)

- Frame: all frozen WGS runs (§4) **minus** Tier E1 runs (project-level
  exclusion of E1 projects' runs is NOT applied outside E1 — only the E1
  *runs* are removed; remaining runs of an E1 project stay eligible? **No**:
  to keep tiers disjoint by design intent, E2 excludes runs from E1
  BioProjects entirely).
- Stratify by BioProject; per-project cap **C = 60** runs.
- Sampling: global seed **20260822**; per-project RNG =
  `random.Random(f"{seed}|{bioproject}")`, sample
  `min(n_avail, 60)` runs per project (sorted run list, `rng.sample`);
  projects processed in `random.Random(seed).shuffle` order of the
  project list sorted by (run count desc, accession asc); accumulate until
  exactly **12,000** E2 runs (the crossing project contributes only the
  remaining quota, taken as the first q of its per-project sample).
- Exclusion rules: library_strategy ≠ WGS; scientific_name filter (§4);
  Tier E1 runs/projects; duplicate run accessions. **No mid-run scope
  growth without a logged amendment.**
- Manifest frozen (seed, strata sizes, caps, run list, sha256) in
  `manifests/tierE2_manifest.json` before any download.

## 7. Control gates (evaluated FIRST; any failure → STOP + document)

1. **Positive controls:** `ECBAIT_POSCON_LAMBDA_01` must reach
   `kmer_coverage ≥ 0.7` against full NC_001416.1 **and** minimap2 identity
   ≥ 0.7 on the matched sequence; same for `ECBAIT_POSCON_HK97_02` vs full
   NC_002167.1. (References: panel `public_refs/`, sha256 in provenance.)
2. **Negative controls:** `ECBAIT_SHUF_CONTROL_{01,02}` must produce **0
   screen hits ≥ 0.7** against the two control references **and every
   screened accession** (checked continuously; any violation ⇒ STOP).
3. **Host controls (diagnostic, panel-expected `host_like`):**
   `ECBAIT_HOST_CONTROL_{01,02}` are E. coli host-derived; in an E. coli
   screen they are **expected to hit many accessions**. They are excluded
   from bait-hit classification, gates, and confirmation volume; their
   accessions-hit rate is reported as the host-mask sanity diagnostic. A
   zero-everywhere host-control rate is a documented anomaly (would suggest
   the screened set is not E. coli), not a pass/fail gate.
4. No threshold adjustment is permitted between gate evaluation and
   reporting.

## 8. Execution (identical protocol to local-ntm-logan)

- **Screen**: k = 31 canonical k-mers; per-(bait, accession)
  `kmer_coverage = |bait 31-mers present| / |distinct bait 31-mers|`;
  hit ⇔ ≥ 0.7 (0.5/0.9 reported for context only). Target:
  `https://s3.amazonaws.com/logan-pub/c/<acc>/<acc>.contigs.fa.zst`.
- **Confirmation** (per non-control hit): b2s localization
  (`--out-sequences --output-mapping-positions`) → minimap2 `-c --eqx`
  ladder `sr → asm5 → asm20` → query coverage + identity + single-contig.
  Junction/synteny claims require **same-contig** support with a
  full-length bait alignment (a full-length junction-bait alignment spans
  the junction it encodes — NTM §8 logic, unchanged).
- **Volume guard (preregistered):** if non-control (bait, accession) pairs
  ≥ 0.7 exceed **40,000**, confirmation proceeds in deterministic seeded
  order (seed 20260822) and truncation is logged as an amendment (no
  threshold change).
- **Classification ladder (frozen):** `no hit` → `module-only homology` →
  `junction/synteny supported` → `near-complete analogue`
  (**≥ 2 non-control baits ≥ 0.7 in one accession**, each
  alignment-confirmed — task-definition; the NTM run's stricter
  junction+2-module+20 kb rule is reported alongside as a secondary,
  non-promoting statistic). "Not detected" ≠ absence: subset + assembly-only
  caveats carried; effective-independence audit collapses hits by
  BioProject and by identical zst sha256.
- **Budgets (hard stops, checked before every tier and continuously):**

| Budget | Value |
|---|---|
| Downloaded bytes (zst, cumulative) | 400 GB |
| Distinct accessions downloaded | 30,000 |
| Wall clock from RUN_START | 48 h |
| S3 concurrent connections | ≤ 4 |
| Retries | ≤ 3 per file, backoff 5/20/60 s |
| EUtils request rate | ≤ 1 req/s (no API key) |
| Logan Search spend | 0 |

- **Disk preflight before each tier**: S3 HEAD stratified sample (≥ 500,
  seed 20260822) → projected volume vs NVMe free space.
- Every downloaded file records URL, HTTP status, size, ETag,
  Last-Modified, sha256, timestamp in `downloads/download_ledger.jsonl`
  (append-only). 404s recorded as `not_in_logan` (never retried).
- Bulky files live only under
  `/mnt/nvme3n1/erikg/phind-genome-work/ecoli_screen/`; the repo mirrors
  code, manifests, ledgers, hit tables, summaries, and the report
  (`NVMe_MANIFEST.tsv` receipt).

## 9. Toolchain diffs vs `ntm/v2/external_validation/local_screen` (path/generalization only)

| script | diff |
|---|---|
| `fetch_manifest_ecoli.py` | E. coli query table (§3) + `usehistory`/`retstart` pagination (universe > 100k); writes universe manifest |
| `fetch_runinfo.py` | verbatim copy |
| `freeze_manifest.py` | filter constants → E. coli (§4); tier names E1/E2; no A2bare clause |
| `build_bioproject_table.py` | NEW — §4 step 4 + §5 archetype matching + §6 E2 sampler (stratification did not exist in NTM) |
| `preflight_disk.py` | verbatim copy (defaults: seed 20260822) |
| `controls_gate.py` | bait ids parameterized (λ/HK97 poscons, SHUF negs); HOST controls screened as diagnostics |
| `screen_accessions.py` | verbatim copy (budget defaults 400 GB / 30k) |
| `confirm_tier_hits.py` | bait-class tokens (`ECBAIT`, `SHUF`/`HOST` exclusions); near-complete rule per §8 task definition (NTM rule kept as secondary statistic) |
| `summarize_screen.py` | ROLES → E. coli controls; pilot merge → per-genome functional report (`research/phage_annotation/per_genome_annotation_qc.tsv`) |
| `known_phage_calibration.py` | NEW — deliverable 4 (§10) |
| `emit_nvme_manifest.py`, `setup_env.sh` | verbatim copies, NVMe root → ecoli_screen |

Screen/gate/confirmation math, thresholds, ladder order, budgets, and
provenance records are byte-comparable in intent to the validated NTM run.

## 10. Deliverables

1. `manifests/tierE1_manifest.json`, `tierE2_manifest.json` (queries, seed,
   strata, sha256) + frozen TSVs — preregistered before downloads.
2. Controls-gate receipts; download ledger + availability; disk preflights.
3. Per-bait hit tables, alignment summaries, full-provenance TSV
   (accession → S3 URL → zst sha256 → bait → contig + span + CIGAR).
4. Per-clade external-evidence update merged with the per-genome functional
   report (`research/phage_annotation/per_genome_annotation_qc.tsv`).
5. Calibration vs NTM: junction-support rate (E. coli 59 junction baits vs
   NTM 5), near-complete analogue rate (NTM 0), clone-collapse audit
   (BioProject + sha256 collapse).
6. Known-phage calibration: hits matching λ/HK97/P2-class references vs
   none. P2 NC_001895.1 fetched from NCBI EFetch (sha256-logged); a
   confirmed hit is "known-phage-matched" iff the bait aligns to a known
   reference with qcov ≥ 0.5 AND identity ≥ 0.7; else "novel (no
   known-phage match)". Poscon-bait hits are λ/HK97-positive accessions by
   construction.
7. `REPORT.md` with comparison-vs-NTM section + GO/NO-GO (full coverage or
   102-bait full panel) with budget math; `NVMe_MANIFEST.tsv` receipt.

## 11. Environment (reproducible)

`scripts/setup_env.sh` recreates the micromamba env at
`<nvme_root>/env` with pinned `back_to_sequences=0.8.4`, `minimap2=2.31`
(conda-forge/bioconda) — identical versions to the validated NTM env
(`env/VERSIONS.txt` receipt).
