# NTM v3 full-cohort run — stage playbook (no-WG run)

**Purpose.** Execute the NTM v3 all-inclusive campaign end-to-end on the
collaborator delivery (`/home/erikg/phind/NTM_Data`, 2026-09-22), using **only
`pi process` + `pi subagents`**. The WG graph is deliberately frozen for this
run and must not be touched.

**Orchestration split (do not blur):**

- **Parent chat = orchestrator / "guiding light".** Owns the DAG, launches and
  monitors the heavy detached `process` jobs, holds decisions and provenance,
  validates at stage barriers. Parent runs in `/home/erikg/phind`.
- **Subagents = bounded stage work.** Validation, reconciliation, report
  writing, shard preparation, failure diagnosis, independent review. Run async,
  keyed, dependency-ordered through one workflow script.
- **`pi process` = the long compute.** Nothing >5 min runs in a foreground
  call, in the parent *or* a child, or it re-creates the `ntm-v3-ml` provider
  stream-timeout death.

**Model routing.** Grinder/validator lanes → `lunaroute/deepseek-4.1-flash`.
Independent verification + release reconciliation → keep a strong tier
(`lunaroute/deepseek-4.1-flash` is fine for grinders; use `:high` thinking or
`glm-5.3` for the release/verify lanes). Verify exact routes with
`subagent({action:"models"})` before launch.

---

## Data layout

| path | meaning |
|---|---|
| `/home/erikg/phind/NTM_Data/*.zip` | collaborator delivery, 42 GB (frozen source) |
| `/home/erikg/phind/NTM_Data/NTM_master_table_complete_20260922.tsv` | 26,499 genomes × 28 cols |
| NVMe `.../ntm/v3/genomes/delivery_raw/{sra_assembled,ncbi_bvbrc}/<genome_id>.fasta` | Stage 0a output |
| NVMe `.../ntm/v3/genomes/canonical_objects/<canonical_acc>/<canonical_acc>.pansn.fa.gz{,.fai,.gzi}` | Stage 0b output (pipeline input) |
| NVMe `.../ntm/v3/{host_clades,mash_clades,clades,ml,full_prophages.fa}` | downstream outputs |
| repo `ntm/v3/inputs/*` | frozen manifests + source CSVs |

`NVME=/mnt/nvme3n1/erikg/phind-genome-work`

**Provenance decision (open, parent).** Completing v3 in place supersedes the
committed public-portion `clade_summary.tsv` / `partition_report.md` /
`partition_validation.txt` and the smoke ML output. Snapshot the current
public-only generation before Stage 2 starts.

---

## Reconciliation invariants (must hold after ingest)

1. delivery FASTAs == 26,499 (17,920 sra_assembled + 8,579 ncbi_bvbrc).
2. every `genome_id` in `ntm_qc_passed_phigaro_coordinates_20260909.csv`
   (17,117 unique) has a delivered FASTA — **0 orphans** (verified pre-ingest).
3. scaffold names resolve verbatim per namespace: `NODE_*` (ASSEMBLY),
   `accn|*` (BV-BRC), versioned accession (NCBI) — **verified on 3 samples**.
4. `canonical_objects/` has one object per unique `canonical_acc` from
   `v3_acquisition_manifest.tsv.gz`; no duplicate/missing object.
5. `blocked_run_rows` 27,494 → 0; extraction errors 0.

---

## Scope decision (2026-10-06) — unified Mycobacteriaceae catalog

**DECIDED: one catalog spanning the whole family. Tag, do not filter.**

- Analysis universe = the union (34,163 canonical objects: 26,499 export
  + 8,212 v2-only). No population is excluded.
- `host_scope` ∈ {NTM, MTC} is a **label column** on every ML/ancestral genome
  and every prophage clade — never a filter.
- NTM/MTC split: 26,890 objects are NTM, 7,273 are MTC (all MTC are v2-only;
  0 MTC in the export). 391 prophages / 12 wholly-MTC clades are entries with
  `host_scope=MTC`, not casualties.
- Safety: tight clades are MASH-defined, so a dissimilar MTC phage cannot drag
  unrelated NTM phages into a clade — family-wide clustering is safe by
  construction.
- Required release artifact: a **cross-boundary report** — prophage clades
  whose hosts span NTM↔MTC, plus shared-vs-distinct prophage structure by
  scope. This is the mutual-information payoff and is impossible in split runs.
- Release notes must carry an explicit scoping statement (counts by scope) so
  an NTM-only subset remains defensible while MTC data stays available.
- Supersedes the earlier "drop MTC at release" proposal.

---

## Stages

### S0a — raw ingest  *(PARENT process — RUNNING)*
- cmd: `bash ntm/v3/scripts/ingest_delivery_raw.sh`
- out: `$NVME/ntm/v3/genomes/delivery_raw/…`; expect `INGEST_DONE … total=26499`
- process name: `ntm-v3-ingest-raw`

### S0b — PanSN canonical objects  *(PARENT process, resumable)*
- build `canonical_objects/<canonical_acc>/<canonical_acc>.pansn.fa.gz{,.fai,.gzi}`
- reuse `rename_fasta_to_pansn()` + `write_bgzip_faidx()` from
  `ntm/v3/scripts/acquire_v3_genomes.py`
- map `genome_id → canonical_acc` from `v3_acquisition_manifest.tsv.gz`
  (dedup: one object per canonical_acc; GCA-preferred on twins)
- skip objects that already exist with a valid `.fai` (resume)
- hard rule: keep existing v1/v2 objects untouched; write only into v3 tree
- **interface risk to re-check at scale:** `accn|`/`NZ_`/unversioned normalization
  in the extractor's `load_fai` resolution (already proven on 3 samples)

### S0c — manifest rebuild (delivery)  *(SUBAGENT worker — DONE)*
- new files: `ntm/v3/inputs/v3_acquisition_manifest_delivery.tsv.gz`,
  `v3_prophage_manifest_delivery.tsv.gz`; NVMe working copies staged
- flipped 17,920 `blocked_run → delivered_run`; 135 V2_RUN stay blocked
  (genuinely undelivered)
- unified manifest: extractable 9,446 → **36,857**, blocked_run 27,494 → **83**,
  total 36,940 unchanged; report `ntm/v3/s0c_manifest_rebuild_report.md`
- new CLI overrides (defaults unchanged): `--acquisition`/`--out-manifest`
  (builder), `--manifest` (extractor)

### S1 — ingest verification  *(SUBAGENT, read-only — DONE, PASS)*
- assert all 5 reconciliation invariants; sha256 the zips; spot-extract 10
  prophages across namespaces; emit `ntm/v3/ingest_validation.md` PASS/FAIL

### S2 — prophage extraction (full)  *(PARENT process)*
- cmd: `python3 ntm/v3/scripts/extract_v3_prophages.py --manifest ntm/v3/inputs/v3_prophage_manifest_delivery.tsv.gz --work-dir /mnt/nvme3n1/erikg/phind-genome-work`
- out: `$NVME/ntm/v3/full_prophages.fa` + `.manifest.tsv`; expect ~36,857 records
  (9,446 → full unified manifest), 0 errors; 83 V2-only rows stay blocked
- dep: S0b, S0c

### S3 — prophage MASH clades  *(PARENT process)*
- cmd: `python3 ntm/v3/scripts/full_ntm_mash_clades_v3.py --phase all --threads 64`
- then `python3 ntm/v3/scripts/validate_mash_clades_v3.py`
- out: `clades/` (tight_clades.json.gz, clade_summary.tsv, SHA256SUMS)
- dep: S2

### S4 — per-clade allwave + impg partition  *(SHARDED SUBAGENT lanes + PARENT processes)*
- per clade: `scripts/per_clade_alignment_pipeline.py` (allwave `tree:k:0:0`,
  window 500, no stranger-joining)
- resume helper: `ntm/v3/scripts/resume_per_clade_partition.py`
- shard unit = contiguous clade-id range; one lane per shard, **disjoint clade
  dirs** (no shared writes). Merge only `partition_summary.tsv` shards in parent.
- validate: `ntm/v3/scripts/validate_v3_partitions.py --summary-tsv …`
- this is the long pole (v3-public: 65,970 s CPU for 456 clades) → chunked,
  detached, resume-first
- dep: S3

### S5 — ML + ancestral traversal  *(PARENT process + SUBAGENT verify)*
- cmd: `python3 ntm/v3/scripts/build_v3_ml_genomes.py --jobs 16 --n-samples 5 --seed 42`
- out: `ml/all_ntm_v3_{ml,ancestral}_phage_genomes.fa`, `release_manifest.tsv`,
  `per_clade_stats.tsv`, `warnings.tsv`
- validate: `ntm/v3/scripts/validate_v3_ml.py`; invariants ML==clades,
  ancestral==alignable, determinism re-run 10 clades × 2 modes byte-identical
- dep: S4

### S6 — Pharokka annotation  *(SHARDED SUBAGENT lanes)*
- per ML genome GFF3 + committed index; prior art
  `ntm/v2/scripts/run_annotation_v2.py` (stages/resume)
- dep: S5

### S7 — release build  *(SUBAGENT, strong tier)*
- join host clades (`host_clades_mash_v3.py` output) × ML genomes; adapt
  `ntm/v2/scripts/build_v2_release.py`; emit RELEASE.md + manifest + labels
- dep: S5, S6, S1(host clades)
- **host clades (S-H) run in parallel:** `ntm/v3/scripts/host_clades_mash_v3.py`
  on the 26,499 cohort (parent process), independent of S2–S6

### S8 — end-to-end independent verification  *(SUBAGENT, fresh context)*
- recompute counts from artifacts, compare vs every stage report, no number
  drift, sha receipts present; emit PASS/FAIL

---

## Work-graph edges

```
S0a → S0b → S1 → {S2, S-H}
S2 → S3 → S4 → S5 → S6 → S7 → S8
S-H ─────────────────────→ S7
```

`S-H` (host clades) is the only true parallel branch; S2–S6 are a chain.

## Failure / resume rules

- Every stage writes a machine-checkable summary + sha receipts on NVMe.
- Drivers are resume-first (`--no-resume` only for deliberate redo).
- A killed lane never invalidates completed per-unit dirs.
- On a provider timeout, the parent re-attaches to the detached process; the
  lane is restarted from its shard boundary, not from zero.
