# S0c — NTM v3 delivery manifest rebuild report

Lane: S0c (manifest rebuild only; no genome ingest, no extraction run).
Date: 2026-10-06. cwd `/home/erikg/phind`.
`$NVME = /mnt/nvme3n1/erikg/phind-genome-work`.

## 1. Deliverable paths

| artifact | path | sha256 |
|---|---|---|
| CLI override — acquisition + out-manifest | `ntm/v3/scripts/build_v3_prophage_manifest.py` | (uncommitted edit) |
| CLI override — manifest | `ntm/v3/scripts/extract_v3_prophages.py` | (uncommitted edit) |
| delivery acquisition-manifest builder | `ntm/v3/scripts/build_v3_acquisition_manifest_delivery.py` | (new) |
| **new acquisition manifest** | `ntm/v3/inputs/v3_acquisition_manifest_delivery.tsv.gz` | `cfa51b1b7c0fa9544f52650992309898bc731fa1c4a7c980b4102824b00428b4` |
| **new prophage manifest (repo)** | `ntm/v3/inputs/v3_prophage_manifest_delivery.tsv.gz` | `1615e0c7fff2cf31f9e0c32e30bd1e8df3aa42ca43ab6b5d85349540ba4e7c6c` |
| NVMe working copy (regenerated) | `$NVME/ntm/v3/inputs/v3_prophage_manifest.tsv` | contents sha256 `b8770bf76328176096825dd1fbcb99d38e6c37957761d882ae8267faf54625f7` (== decompressed delivery manifest) |

Frozen repo inputs were **not** modified: `git status --porcelain` on
`v3_acquisition_manifest.tsv.gz`, `v3_prophage_manifest.tsv.gz` and
`coverage.tsv` is clean. Only the NVMe *working* copy
(`v3_prophage_manifest.tsv`) was overwritten, as instructed.

## 2. Exact commands

```bash
cd /home/erikg/phind

# 2a. build the delivery-updated acquisition manifest (gzip -n, deterministic)
python3 ntm/v3/scripts/build_v3_acquisition_manifest_delivery.py
#   -> delivered run FASTAs in .../delivery_raw/sra_assembled: 17920
#   -> flipped blocked_run -> delivered_run: 17920 ({'ASSEMBLY': 17920})
#   -> still blocked_run: 135 ({'V2_RUN': 135})
#   -> wrote ntm/v3/inputs/v3_acquisition_manifest_delivery.tsv.gz (34846 rows)

# 2b. regenerate the unified prophage manifest from the new acquisition manifest
python3 ntm/v3/scripts/build_v3_prophage_manifest.py \
  --acquisition ntm/v3/inputs/v3_acquisition_manifest_delivery.tsv.gz \
  --out-manifest ntm/v3/inputs/v3_prophage_manifest_delivery.tsv.gz
#   -> unified_manifest_rows 36940
#   -> unified_extractable_rows 36857 / unified_blocked_run_rows 83
#   -> wrote repo delivery gz + plain copy $NVME/ntm/v3/inputs/v3_prophage_manifest.tsv
```

CLI additions (defaults = previous hardcoded paths, so a no-arg run is
unchanged):
- `build_v3_prophage_manifest.py`: `--acquisition` (default `ACQ_TSVGZ` =
  `ntm/v3/inputs/v3_acquisition_manifest.tsv.gz`), `--out-manifest`
  (default `ntm/v3/inputs/v3_prophage_manifest.tsv.gz`). The NVMe plain copy
  keeps its fixed path `$NVME/ntm/v3/inputs/v3_prophage_manifest.tsv`.
- `extract_v3_prophages.py`: `--manifest` (default `MANIFEST_GZ` =
  `ntm/v3/inputs/v3_prophage_manifest.tsv.gz`).

## 3. Acquisition-manifest flip (deliverable 2)

Baseline (frozen `v3_acquisition_manifest.tsv.gz`, 34,846 rows):
`download_bvbrc 70`, `link_v1v2 16688`, `download_ncbi 33`,
`blocked_run 18055` (17,920 `ASSEMBLY` + 135 `V2_RUN`); every blocked row had
an empty `canonical_acc`.

Delivered FASTAs found: **17,920** in
`$NVME/ntm/v3/genomes/delivery_raw/sra_assembled/`, exactly the 17,920
`blocked_run` `ASSEMBLY` accessions.

| result | count |
|---|---:|
| blocked_run rows flipped → `delivered_run` | **17,920** (all `source_ns=ASSEMBLY`) |
| blocked_run `ASSEMBLY` rows with FASTA **not** found | **0** |
| blocked_run `V2_RUN` rows (no delivered FASTA) left untouched | **135** |

Flipped rows set `plan=delivered_run`, `resolution_method=collaborator_delivery`,
`canonical_acc=accession` (non-empty). Field-level diff old vs new: exactly
17,920 rows changed, and **only** the fields
`{canonical_acc, plan, resolution_method}` — verified programmatically. All
other rows/columns are byte-identical (the CSV round-trip reproduces the
frozen file byte-for-byte when no flip is applied). Re-running the builder
reproduced identical bytes (sha256 stable, gzip `mtime=0`).

**V2_RUN decision:** the 135 `V2_RUN` blocked rows have no delivered FASTA in
`sra_assembled/`, so they are left at `plan=blocked_run` with empty
`canonical_acc` (the collaborator delivery covers only export `ASSEMBLY` runs).
They are genuinely blocked, not silently dropped.

Note: flipped rows keep their original `notes` text ("collaborator run assembly
not delivered locally … PENDING"), which is now stale. Per the task contract
("leave every other column untouched") it was not rewritten — flagged as a
residual risk below.

## 4. Prophage-manifest before/after (deliverable 3)

The builder asserted, with no failures: `(canonical_acc, prophage_id)`
uniqueness, one caller source per canonical genome, and
`length == end-begin+1` for all rows.

| status | before (`v3_prophage_manifest.tsv.gz`) | after (`..._delivery.tsv.gz`) | delta |
|---|---:|---:|---:|
| `extractable` | 9,446 | **36,857** | +27,411 |
| `blocked_run` | 27,494 | **83** | −27,411 |
| **total rows** | **36,940** | **36,940** | **0** |

The 27,411 flipped prophage rows are all `source=BV-BRC` phigaro rows whose
`canonical_acc` is a run accession (ERR 14,323 / SRR 11,358 / DRR 1,730
distinct runs = 11,323 distinct objects; the remaining 17,920−11,323 = 6,597
delivered runs carry `prophage_count=0` and produce no manifest rows). The 83
remaining `blocked_run` rows are all `source=V2` (36 v2-only prophage-bearing
runs — no delivered FASTA).

Reconciliation: `9,446 + 27,411 = 36,857` extractable, total unchanged at
36,940.

## 5. Verification results

1. **Flipped rows have non-empty `canonical_acc == accession` and a delivered
   FASTA.** PASS — of 27,411 rows that flipped `blocked_run → extractable`,
   0 had an empty/mismatched `canonical_acc`; 0 lacked a delivered
   `$NVME/.../sra_assembled/<acc>.fasta`.
2. **`extractable == 9,446 + flipped` and total rows unchanged at 36,940.**
   PASS — expected 36,857, observed 36,857; total 36,940 → 36,940.
3. **`(canonical_acc, prophage_id)` uniqueness.** PASS — builder assertion
   passed; independently re-checked on the output: 36,940 keys, all unique.
   No collision between the 11,323 run canonical objects and the pre-existing
   extractable canonical objects (intersection = 0); no run accession equals a
   GCA_/GCF_ object.
4. **Defaults unchanged.** PASS (dry inspection) — module constants
   `ACQ_TSVGZ = <repo>/ntm/v3/inputs/v3_acquisition_manifest.tsv.gz` and
   `MANIFEST_GZ = <repo>/ntm/v3/inputs/v3_prophage_manifest.tsv.gz` are the
   argparse defaults; the out-manifest default literal is the frozen repo gz
   path. No old output regenerated.

Note on concurrency: `$NVME/ntm/v3/genomes/canonical_objects/` already holds
4,311 run-like dirs (4,288 with `*.pansn.fa.gz`) — a parallel ingest lane is
materialising the delivered runs. Of the 11,323 newly-extractable run objects,
**5,766** currently have a local pansn object; **5,557** do not yet.

## 6. Proposed `coverage.tsv` values (do NOT overwrite yet)

Manifest-derived metrics that change (values below are from the new delivery
manifest):

| metric | old | proposed | note |
|---|---:|---:|---|
| `extractable_rows` | 9,446 | **36,857** | genome local / delivered |
| `blocked_run_rows` | 27,494 | **83** | 36 v2-only runs, no delivery |
| `distinct_objects_read` | 3,994 | **15,317** | distinct `canonical_acc` among extractable rows (extraction re-run confirms) |
| `v2_run_rows_in_union` | 14,741 | 14,741 | unchanged (14,658 + 83) |
| `fasta_single_source_per_genome` | 36,940 | 36,940 | unchanged (asserted at build) |
| `length_equals_end_minus_begin_plus_1` | 36,940 | 36,940 | unchanged |

Extraction-dependent metrics (`extracted_records`, `extraction_errors`,
`total_bp`, `len_min/median/mean/max`, `begin_zero_rows_extracted`,
`out_of_range_rows`, `end_beyond_contig_rows`, `fasta_sha256`,
`fasta_records`) are **pending the extraction re-run** on the new manifest —
the projected values are `extracted_records = 36,857`, `extraction_errors = 0`,
`fasta_records = 36,857`; length/byte totals and the FASTA sha256 must come
from the actual extract run. `coverage.tsv` remains untouched.

## 7. Residual risks

- **V2_RUN handling** (135 rows, 83 prophage rows, 36 objects): left blocked;
  no delivered FASTA exists for them under `sra_assembled/`. They remain a
  genuine gap unless a separate V2-run delivery arrives.
- **Stale `notes`** on the 17,920 flipped acquisition rows still say the
  assembly was not delivered (task contract required leaving the column
  untouched). If downstream consumers read `notes`, the delivery run needs a
  follow-up rewrite.
- **canonical_acc collisions**: none found — run accessions (ERR/SRR/DRR)
  cannot collide with GCA_/GCF_ objects, and the 11,323 run objects do not
  intersect the 9,446 old extractable canonical objects. Uniqueness is
  enforced by the builder assertion regardless.
- **Extraction not yet possible for 5,557 of 11,323 newly-extractable run
  objects** (no local pansn object yet). The extraction re-run must wait for
  the ingest lane to finish, or will report scaffold/faidx errors for those.
- **NVMe working copy overwritten**: `$NVME/ntm/v3/inputs/v3_prophage_manifest.tsv`
  now holds the delivery manifest. If any consumer expected the old copy, it
  must use the frozen repo gz instead.
- Pre-existing `SyntaxWarning` (`\|` escape) in `extract_v3_prophages.py` line
  377 is unrelated to this change and unchanged.

## 8. Recommended next step

Once the ingest lane finishes materialising all 11,323 run pansn objects, run
the extraction re-run:

```bash
python3 ntm/v3/scripts/extract_v3_prophages.py \
  --manifest ntm/v3/inputs/v3_prophage_manifest_delivery.tsv.gz
```

and emit a delivery `coverage.tsv` from the resulting stats (do not overwrite
the frozen `coverage.tsv`).
