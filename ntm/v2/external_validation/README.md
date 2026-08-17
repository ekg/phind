# Public mycobacteriophage reference panel (external validation)

Reproducible public reference panel of actinobacteriophage genomes for
comparison with the PHIND NTM v2 (ML / ancestral) phage genomes.

**Panel (build of 2026-08-17):** 6,266 retained records — 3,065
mycobacteriophages (hosts in *Mycobacteriaceae*) + 3,200 other
actinobacteriophages (PhagesDB leg: *Gordonia*, *Arthrobacter*,
*Microbacterium*, *Streptomyces*, …) + 1 host-unresolved record.

## Sources (all responses cached; see PROVENANCE.md)

| source | artifact | scope |
|---|---|---|
| PhagesDB | `/data/?set=seq&type=full` (metadata TSV), `/media/Actinobacteriophages-All.fasta` (bulk FASTA) | **all** sequenced actinobacteriophages (broad host taxonomy) |
| INPHARED (Millard Lab, release 7Apr2026) | `7Apr2026_millardlab_website_table.txt.gz`, `7Apr2026_genomes.fa.gz` (s3.climb.ac.uk) | mycobacteriophages only (Mycobacteriaceae hosts / mycobacteriophage naming) |
| NCBI GenBank/RefSeq (E-utilities) | esearch `mycobacteriophage[All Fields] OR "Mycobacterium phage"[Title]` (nuccore), batched efetch FASTA (versioned accessions) + esummary JSON + backfill-by-id + GenBank remainder for `/host=` qualifiers | mycobacteriophages + text-search artifacts screened by explicit rules |

Dedup: accession-base aliasing across sources + exact-sequence sha256
merging (RefSeq↔GenBank twins). Every alias/duplicate/variant stays in
`crosswalk.tsv`; every drop has an explicit reason in `exclusions.tsv`.

## Layout

```
scripts/panel_lib.py      pure parse/merge/emit library (no network, no clock)
scripts/fetch_sources.py  idempotent cache fetcher (bulk | ncbi | taxonomy |
                          backfill | gb stages), writes cache/fetch_manifest.json
scripts/build_panel.py    offline deterministic builder -> external panel/
scripts/verify_panel.py   acceptance checks V1–V6 (incl. byte-identical rerun)
scripts/test_panel.py     unit tests on synthetic fixtures (pytest)
scripts/sync_repo.py      copies summaries + provenance into this directory

summaries/REPORT.md       coverage, host distribution, evidence classes, caveats
summaries/PROVENANCE.md   source URLs, retrieval timestamps, sha256 per artifact
summaries/manifest_panel.tsv.gz   full manifest (one row per retained record)
summaries/crosswalk.tsv.gz        alias / duplicate / variant crosswalk
summaries/exclusions.tsv          explicit exclusion reasons
summaries/source_summary.tsv      count reconciliation
summaries/host_distribution.tsv   host-label distribution
summaries/CHECKSUMS.sha256        sha256 of every external panel artifact
summaries/panel_files.txt         external panel file sizes
```

Bulky FASTA lives externally only, at
`/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation/panel/`:

- `panel.full.fasta` — all records (unique `PHPUB-######` IDs)
- `panel.mycobacteriophage.fasta`, `panel.other_actinobacteriophage.fasta`,
  `panel.phage_host_unresolved.fasta`

## Reproduce

```bash
python3 scripts/fetch_sources.py --stage all   # populates ../external cache (idempotent)
python3 scripts/build_panel.py                 # offline, deterministic
python3 scripts/verify_panel.py --determinism  # V1–V6 acceptance checks
python3 scripts/sync_repo.py                   # refresh repo summaries
```

`fetch_sources.py --stage all` needs network access to phagesdb.org,
s3.climb.ac.uk and eutils.ncbi.nlm.nih.gov. Everything downstream of the
cache is pure: a rerun from the same cached responses reproduces manifests,
reports and checksums byte-for-byte (verify V6).

## Manifest columns (headline)

`panel_id`, `primary_accession` (NCBI **versioned** where available),
`primary_source`, `scope`, `length_bp`, `seq_sha256`,
`fasta_bytes_sha256`, `host_reported` (+ `host_evidence` — always a
**reported label**, never verified host range), `host_genus`, `host_taxid`,
`evidence_class` (`isolated_sequenced` / `deposited_isolate_sequence` /
`predicted_prophage` / `metagenome_assembled` / `synthetic_construct`),
`topology`, `completeness`, phagesdb cluster/isolation metadata,
inphared accession/family, ncbi title/taxid, `gb_host_qualifier`,
`sources`, `seq_source`, `retrieval_timestamp`, `source_url`,
`literature_reference`.
