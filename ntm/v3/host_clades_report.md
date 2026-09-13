# NTM v3 — host clades (MASH) on the all-inclusive cohort

Generated: 2026-09-13 (task `ntm-v3-host`). Artifacts on NVMe:
`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/host_clades/` (sha256 receipts
below and in `SHA256SUMS` on NVMe). Driver:
`ntm/v3/scripts/host_clades_mash_v3.py`.

## Cohort

- **16,243 canonical genome objects** from `acquire-ntm-v3`
  (`ntm/v3/genomes/canonical_objects/<acc>/<acc>.pansn.fa.gz`) — every genome
  physically on this host. All 16,243 sketched, all 16,243 clustered, all
  16,243 labelled.
- **18,055 run assemblies (ERR/SRR/DRR) excluded** — no FASTA exists anywhere
  on this host (collaborator upload pending; see `ntm/v3/download_report.md`
  "blocked" section). They cannot be sketched; when they land, re-running the
  driver regenerates the full table unchanged.

## Method (v2 parity)

Identical MASH parameters and clustering rule to v1/v2:

1. `mash sketch -k 21 -s 10000` per genome (48 parallel workers)
2. `mash paste` → combined `host.msh` (1,304 MB, 16,243 sketches)
3. `mash triangle -k 21 -s 10000 -p 64 host.msh` → `host.dist`
   (complete lower triangle, **131,909,403 pairs = n(n-1)/2** for n=16,243,
   0 malformed rows; 5.0 min)
4. Union-find connected components (single-linkage) at **dist ≤ 0.05**
   (~95% ANI, species level — the v1/v2 convention); components sorted
   size-desc; fresh ids `host_clade_0001..host_clade_0411`.

One real bug was found and fixed during validation: `mash triangle` prints
near-identical distances (< ~1e-4) in scientific notation (`6.91981e-05`);
a string fast-path initially treated those as ≥0.1 and skipped 490,759
union operations. The v2-cohort parity check (below) surfaced it (one
2-member v2 clade split); after the fix, clustering uses `float()` for all
non-plain-format values. Final: 37,208,619 unions, 411 clades (pre-fix
458 — 47 spurious splits).

### Nohup/resume protection (lesson from operations.jsonl 2026-08-16)

The v2 host task died because a provider stream timeout killed the attempt
AFTER compute completed. v3 runs the heavy phase (sketch + paste + triangle)
detached (`nohup setsid … --phase heavy`), with per-genome sketches as
resume units on disk; the parse/cluster/emit phase (`--phase light`) reads
only artifacts already on disk and is a pure function of them.

**Resume proof (executed, not just designed):** the driver was `kill -9`-ed
(process group) mid-sketch at 14,550/16,243 sketches; re-invocation logged

```
sketch: 16243 genomes total, 14550 already on disk (resume), 1693 to sketch
```

and completed the remaining 1,693 before pasting — no sketch recomputed.
`state.json` records `sketch.resumed = 14550`.

### Species labels

- 8,031 genomes: species from the v3 acquisition manifest (BV-BRC/NCBI
  export label; BV-BRC type-strain strings like
  "Mycobacteroides abscessus ATCC 19977" are kept verbatim)
- 8,212 genomes: v2 QC-passed accession list
  (`ntm/v2/inputs/NTM_QC_passed_accession_list.tsv`) — V2_NUMERIC-only
  objects whose manifest row carries no species (0 conflicts measured)
- genus = first token of species; organism = species (v2 `host_clades.tsv`
  schema, unchanged)

## Results

- **411 clades**, all 16,243 genomes assigned, exactly one clade per genome
  (singletons: 99)

Clade size distribution (top; 317 further clades of size 1–6):

| n | # clades | dominant species |
|---:|---:|---|
| 7465 | 1 | Mycobacterium tuberculosis (TB complex — 46.0% of cohort, v2-verbatim policy) |
| 4080 | 1 | Mycobacteroides abscessus ATCC 19977 (abscessus complex, 25.1% of cohort) |
| 1189 | 1 | Mycobacterium ulcerans (+ M. marinum / M. pseudoshottsii / M. liflandii) |
| 675 | 1 | M. avium subsp. paratuberculosis ATCC 19698 (+ MAC) |
| 269 | 1 | Mycobacterium paraintracellulare (+ M. intracellulare) |
| 227 | 1 | Mycobacterium kansasii ATCC 12478 (+ M. persicum) |
| 136 | 1 | Mycobacteroides chelonae subsp. bovistauri QIA-37 (+ M. chelonae) |
| 129 | 1 | Mycolicibacterium smegmatis |
| 121 | 1 | Mycolicibacterium fortuitum subsp. fortuitum |
| 80 | 1 | Mycolicibacterium septicum DSM 44393 (+ M. nivoides) |
| 41 | 1 | Mycolicibacterium senegalense |

Cohort composition: TB-complex 7,239 (44.6%), abscessus complex 4,081
(25.1%), distinct species strings 190+ (full table: `host_clade_summary.tsv`).

Top 15 species (label strings verbatim): M. tuberculosis 7,233 |
M. abscessus ATCC 19977 2,591 | M. abscessus subsp. massiliense CIP 108297
1,252 | M. ulcerans 1,047 | Mycobacterium sp. 358 | M. abscessus subsp.
bolletii BD 219 | Mycolicibacterium sp. 217 | M. kansasii ATCC 12478 210 |
M. avium subsp. paratuberculosis ATCC 19698 206 | M. avium subsp. avium 151 |
M. avium subsp. hominissuis ATCC 700898 139 | M. smegmatis 129 |
M. paraintracellulare 123 | M. chelonae subsp. bovistauri QIA-37 119 |
M. timonense 114.

## Validation

- [x] **Row completeness:** `host_clades.tsv` has 16,243 data rows ==
  16,243 cohort canonical objects; 0 duplicates; set equality with the
  filelist (`cluster.log`: "validation: PASS")
- [x] **Triangle complete:** 131,909,403 pairs == n(n-1)/2 for n=16,243;
  every row has exactly i fields; 0 malformed
- [x] **Spot-check:** 20/20 sampled pairs re-measured with direct
  `mash dist` on the individual sketches match the triangle values with
  delta = 0.00e+00 (`spotcheck.tsv`, seed 42)
- [x] **No genome in >1 clade:** union-find guarantees it; re-asserted from
  the emitted rows (0 violations)
- [x] **Clade size distribution + per-clade dominant species reported**
  (`host_clade_summary.tsv` — full distribution; top table above)
- [x] **v2↔v3 parity:** of the 764 v2 clades containing genomes shared with
  v3 (11,772 shared), **764/764 nest into exactly one v3 clade** — v3 only
  merges (bridges via new genomes), never splits. v2's TB-complex clade
  (7,465) reappears at exactly the same size; abscessus 2,336→4,080 and
  avium 507→675 grew by export additions without splitting
- [x] **Resume pattern proven:** kill + re-invoke resumed from 14,550
  on-disk sketches (log excerpt above)

## Files (NVMe, sha256 receipts)

| file | bytes | sha256 |
|---|---:|---|
| filelist.txt | 1,850,894 | `4bf541fd…8ca62de` |
| host.msh | 1,303,818,104 | `39fecf03…0d40260` |
| host.dist | 1,266,954,497 | `1c267b90…bcd0e15` |
| host_clades.tsv | 1,797,474 | `adff5fd4…c3de32d` |
| host_clade_summary.tsv | 31,107 | `dfd0f852…cd41f906` |
| cluster.log | 2,377 | `25d1077b…da560` |
| spotcheck.tsv | 1,330 | `f1127e23…7a7d8f4d` |

Full hashes in `SHA256SUMS` next to the artifacts; per-genome sketches in
`host_clades/sketches/` (16,243 × ~80 KB, resume units).

`host_clades.tsv` schema (v2, unchanged): `accession | host_clade_id |
species | genus | organism` — these clade labels go on the final release ML
genomes (downstream: `ntm-v3-release`).

## Extending to run-assemblies

When the 18,055 blocked run assemblies are delivered and linked into
`ntm/v3/genomes/`, re-run the driver end-to-end (`--phase all` under nohup);
sketches, paste, triangle, clustering and labels regenerate with fresh
clade ids over the extended cohort, resuming all existing per-genome
sketches from disk.
