# NTM v3 — all-inclusive cohort acquisition report

Generated: 2026-09-12T22:27:30.506374Z (task `acquire-ntm-v3`)

Cohort: union of the BV-BRC phigaro QC-passed export (26,499 genomes,
2026-09-09) and the v2 local holdings (13,122 NCBI assembly entries +
9,543 run assemblies). Genomes live under
`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/genomes/canonical_objects/`
(PanSN bgzip + faidx, linked objects are symlinks into frozen v1/v2).

## Counts

| metric | count |
|---|---:|
| cohort total (manifest rows) | 34846 |
| — export rows | 26499 |
| — v2-only assembly rows added | 8347 (8,212 numerics + 135 runs) |
| linked from v1/v2 | 16688 rows / 16148 objects |
| — of which BV-BRC ids linked by contig identity | 178 |
| downloaded | 103 rows / 95 objects |
| — from BV-BRC API | 70 |
| — from NCBI (datasets ZIP) | 33 rows / 25 assemblies |
| blocked (run assemblies, collaborator pending) | 18055 |
| failed (acquisition errors) | 0 |
| resolved (linked + downloaded + blocked) | 34846 |
| resolved == cohort total | YES (34846 == 34846) |
| unique namespace-union assemblies | 31,931 (13,628 numerics + 18,055 runs + 248 BV-BRC) |
| unique identity-deduplicated assemblies | 31753 (13,628 GC numerics + 70 BV-BRC-only + 18,055 runs) |
| canonical objects on NVMe | 16243 |
| total bp | 80,772,491,382 |
| total contigs | 1,623,753 |

## Validation

- [x] resolved == cohort total: 34846 == 34846 (every manifest row has a terminal state; zero silent skips)
- [x] `samtools faidx` region query ok on all 16243/16243 canonical objects (existing index used; frozen v1/v2 `.fai` never rewritten)
- [x] `gzip -t` full-stream: 300/300 sampled objects pass (seed=42 sample)
- [x] contig names vs coordinates scaffolds: 10540 matched / 0 mismatched over 5794 non-run prophage-bearing genomes (accn| stripped, unversioned, NZ_-strip fallback — the v2 extraction convention)
- [x] all 13,086 v2 NCBI numerics represented (13086; 0 missing) and all 9,543 v2 run ids present (0 missing); none dropped
- [x] coverage.tsv reconciles exactly with this report (same generated counts)

## Per-source breakdown

| source | rows | linked | downloaded | blocked |
|---|---:|---:|---:|---:|
| NCBI | 8331 | 8298 | 33 | 0 |
| BV-BRC | 248 | 178 | 70 | 0 |
| ASSEMBLY | 17920 | 0 | 0 | 17920 |
| V2_NUMERIC | 8212 | 8212 | 0 | 0 |
| V2_RUN | 135 | 0 | 0 | 135 |

### BV-BRC genome_id resolution (248 taxid.version ids)

All 248 resolved via the BV-BRC API (`/api/genome_sequence/`): 178 were
already held locally and are **linked** by full contig-identity (every
BV-BRC sequence accession present in one local v1/v2 object — see
`inputs/bvbrc_resolution.tsv`); 70 were **downloaded** from the BV-BRC
API (contig accessions equal the coordinates-CSV scaffolds minus `accn|`).
Unresolvable BV-BRC ids: **none**.

### NCBI delta (25 assemblies downloaded)

Only 12 export numerics had no local object at all (20 rows), plus
5 export rows whose exact assembly version is not held locally
(version-skew; the local object under the same numeric is an older
version), plus 8 export rows whose same-version local twin
is the RefSeq copy with renamed contigs (scaffolds unaddressable there).
These 25 assemblies were downloaded at the **exact export accession**
via the NCBI Datasets v2 API (v2 source chain: FTP → datasets ZIP → ENA;
datasets ZIP served all 25). The vast majority of export NCBI numerics
(5,411/5,423) were already on NVMe from v1/v2 and are linked, not
re-downloaded.

## Failures / blocked genomes (explicit, no silent skips)

**18055 run assemblies (ERR/SRR/DRR) are blocked**, not failed
for a fixable reason: the phigaro coordinates were computed on the
collaborator's SPAdes assemblies (scaffold tokens `NODE_*`), which exist
only on the collaborator's cluster. The collaborator upload to this host
is still **PENDING** (`ntm/v2/run_assemblies/REQUEST.md`, 0 of 9,543 v2
run-assembly FASTAs delivered; verified absent from NVMe). ENA submitted
assemblies exist for ~79% of runs but are **content-mismatched** to the
coordinate-bearing assemblies (v2 `ena_backfill` study: 79/79 mismatches,
contig tokens differ), so they are deliberately NOT downloaded as
substitutes. Full accession list with reason:
`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3/genomes/blocked_run_accessions.txt`.
Of the 18,055 blocked rows, 9,408 are shared with the v2 cohort and 8,512
are export-only; 135 further v2 runs are outside the export (also blocked).

Acquisition errors: **0** (all 95 download targets and all 16148 link
targets succeeded).

## Deviations from the task's expected shape (documented)

1. **The task framing assumed the 9,543 v2 run assemblies were already
   on local NVMe** ('v2 downloads (11,972) … including the 9,543
   ERR/SRR/DRR collaborator run assemblies'). They are not — the v2
   collaborator upload is still PENDING (REQUEST.md) and no run-assembly
   FASTA exists anywhere on this host (verified). They are therefore
   *accounted for* as blocked rows rather than linked: all 9,543 v2 runs
   appear as manifest rows (9,408 via export rows + 135 v2-only rows),
   none dropped. When the collaborator delivery lands, linking them into
   the v3 layout is a mechanical rerun of `acquire_v3_genomes.py`.
2. **The expected download delta (~16k) did not materialize** because
   (a) the export is not a pure BV-BRC snapshot — 17,920 of 26,499 rows
   are run accessions with no public sequence source, and (b) v1
   holdings already cover almost all export NCBI numerics. The true
   acquirable delta was 87 objects (70 BV-BRC + 17 NCBI).
3. **563 export NCBI rows link a GCA/GCF twin object** and 178
   BV-BRC rows link a local object by contig identity (RefSeq or GenBank
   twin): extraction must use the documented NZ_-strip fallback
   (`ntm/v2/scripts/extract_full_ntm_prophages.py` convention); all
   10,540 non-run prophage scaffolds verified addressable this way.

## Notes for the prophage-extraction stage (contig-name addressing)

All 10540 non-run prophage scaffolds are addressable, but with
three name-resolution cases the extraction code must handle:

1. **exact name** — 9793 scaffolds are literal contig names;
2. **version suffix** — 168 scaffolds are unversioned while the
   canonical object's contigs are versioned (mostly BV-BRC `accn|` ids
   linked to GCA objects whose GenBank contigs carry `.1`); resolve by
   unversioned-prefix match against the object `.fai`;
3. **NZ_ prefix** — 579 scaffolds need the v2 NZ_-strip fallback
   (RefSeq-style scaffold vs GenBank twin object);
4. **version differs** — 0 scaffolds carry a different
   version suffix than the local contig of the same accession.

The run-assembly (ASSEMBLY/`NODE_*`) scaffolds are only addressable once the
collaborator run assemblies are delivered and linked (see blocked section).

## Artifacts

- `ntm/v3/inputs/v3_acquisition_manifest.tsv` — one row per cohort row with terminal state
- `ntm/v3/inputs/bvbrc_resolution.tsv` — the 248 BV-BRC genome_id resolutions
- `ntm/v3/inputs/coverage.tsv` — machine-readable counts (v2 convention)
- `ntm/v3/scripts/` — build/`acquire`/validate pipeline
- NVMe: `ntm/v3/genomes/` (canonical_objects, acquisition_log.jsonl, progress.json, blocked_run_accessions.txt, validation.json)
