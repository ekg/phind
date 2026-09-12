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
| Acquire cohort (link existing v2 holdings, download delta) | ⏳ next — task `acquire-ntm-v3` |
| Prophage extraction / downstream pipeline | ⏳ to run |

## Inputs

- `inputs/INPUTS.md` — provenance, schemas, sha256 receipts, identifier caveats
  (mixed `ASSEMBLY`/`NCBI`/`BV-BRC` namespaces; BV-BRC type-strain species labels)
- `inputs/v2_overlap_report.md` — v2↔v3 contig-bridge and union analysis

## Working rules

- Genome data will be acquired under `ntm/v3` on NVMe:
  `/mnt/nvme3n1/erikg/phind-genome-work/ntm/v3`. The repo holds only code,
  manifests and small reports (same rule as v1/v2).
- v1/v2 outputs are frozen inputs to v3 — never modify them in place.
