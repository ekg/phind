# NTM v3 — all-inclusive NTM expansion

**v3 = the all-inclusive NTM cohort**: the union of the BV-BRC phigaro
QC-passed export (26,499 genomes / 47,994 prophages, dated 2026-09-09) and the
existing v2 local holdings. See `inputs/v2_overlap_report.md` for the measured
overlap and union math (v3 union ≈ 31,975 unique assemblies: 13,672 NCBI
numeric + 18,055 runs + 248 BV-BRC).

## Status

| step | status |
|---|---|
| Import source inputs (summary + coordinates CSVs, provenance, overlap report) | ✅ done — `inputs/` |
| Acquire cohort (link existing v2 holdings, download delta) | ✅ done — `download_report.md`, `inputs/v3_acquisition_manifest.tsv.gz`, `inputs/coverage.tsv` (34,846 rows resolved: 16,148 objects linked + 95 downloaded; 18,055 run assemblies blocked pending collaborator delivery, none dropped) |
| Host clades (MASH, host structure) | ✅ done — `host_clades_report.md` (16,243/16,243 genomes labelled, 411 clades at dist ≤ 0.05; artifacts + sha256 receipts on NVMe `ntm/v3/host_clades/`) |
| Prophage extraction / downstream pipeline | ⏳ to run |

## Inputs

- `inputs/INPUTS.md` — provenance, schemas, sha256 receipts, identifier caveats
  (mixed `ASSEMBLY`/`NCBI`/`BV-BRC` namespaces; BV-BRC type-strain species labels)
- `inputs/v2_overlap_report.md` — v2↔v3 contig-bridge and union analysis
- `inputs/v3_acquisition_manifest.tsv.gz` — acquisition manifest: one row per
  cohort row (export ∪ v2 holdings) with terminal state and resolution method
  (gzip -n, deterministic; plain copy on NVMe at `ntm/v3/genomes/`)
- `inputs/bvbrc_resolution.tsv` — the 248 BV-BRC `taxid.version` genome_id
  resolutions (178 linked by contig identity, 70 downloaded from the BV-BRC API)
- `inputs/coverage.tsv` — acquisition coverage metrics (v2 convention)
- `download_report.md` — acquisition report with validation results and the
  blocked run-assembly accounting

## Working rules

- Genome data will be acquired under `ntm/v3` on NVMe:
  `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3`. The repo holds only code,
  manifests and small reports (same rule as v1/v2).
- v1/v2 outputs are frozen inputs to v3 — never modify them in place.
